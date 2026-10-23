#!/usr/bin/env python3
"""
Off-target selectivity screen for bridge_48_s01.

Tests the bridge binder against a panel of HMG-box and non-HMG nuclear proteins
to verify it does not cross-react with non-target proteins.

On-target controls (2-chain: HMG + bridge):
  SOX2  → expect high ipTM (designed target)
  SOX9  → expect high ipTM (designed target)

Off-target panel (2-chain: off-target HMG + bridge):
  SOX3  → most critical (SOX B1 group like SOX2; expressed in neural progenitors)
  SOX4  → SOX C group; expressed in brain
  SOX11 → SOX C group; expressed in neurons
  SOX17 → SOX F group; vasculature (brain vessel risk)
  SRY   → SOX A group; most divergent SOX HMG
  HMGB1 → non-SOX HMG box; ubiquitous (structural negative control)

Pass criterion (non-cross-reactive): ipTM < 0.40
(same threshold used for Track 1 individual binder selectivity in Stage 3)

Output: analysis/offtarget_selectivity_results.tsv
"""

import os, json, subprocess
import pandas as pd
from pathlib import Path

BASE      = Path("/mnt/d/gbm-protein-switch")
OUT_DIR   = BASE / "analysis" / "offtarget_selectivity_outputs"
BOLTZ_BIN = Path("/home/ayushregmibagale/miniconda3/envs/boltz2/bin/boltz")

BRIDGE_SEQ = "GAAAAAAAEAAAAAAAAAAAAAAAAAAAAEAARAAAAAAAAAAAAAAAAAAQAAAAQAYLAQAQAAAAA"

# On-target sequences (same as used for MD design)
SOX2_HMG = "DRVKRPMNAFMVWSRGQRRKMAQENPKMHNSEISKRLGAEWKLLSETEKRPFIDEAKRLRALHMKEHPDYKYRPRRKTKTLMKKDKYTLPGGLLAPG"
SOX9_HMG = "PHVKRPMNAFMVWAQAARRKLADQYPHLHNAELSKTLGKLWRLLNESEKRPFVEEAERLRVQHKKDHPDYKY"

# Off-target HMG domain sequences (from UniProt, annotated HMG box regions)
OFF_TARGETS = {
    "SOX3":  "VKRPMNAFMVWSRGQRRKMALENPKMHNSEISKRLGADWKLLTDAEKRPFIDEAKRLRAVHMKEYPDYK",   # P41225 aa139-207
    "SOX4":  "IKRPMNAFMVWSQIERRKIMEQSPDMHNAEISKRLGKRWKLLKDSDKIPFIREAERLRLKHMADYPDYK",   # Q06945 aa59-127
    "SOX11": "IKRPMNAFMVWSKIERRKIMEQSPDMHNAEISKRLGKRWKMLKDSEKIPFIREAERLRLKHMADYPDYK",   # P35716 aa49-117
    "SOX17": "IRRPMNAFMVWAKDERKRLAQQNPDLHNAELSKMLGKSWKALTLAEKRPFVEEAERLRVQHMQDHPNYK",   # Q9H6I2 aa68-136
    "SRY":   "VKRPMNAFIVWSRDQRRKMALENPRMRNSEISKQLGYQWKMLTEAEKWPFFQEAQKLQAMHREKYPNYK",   # Q05066 aa60-128
    "HMGB1": "PRGKMSSYAFFVQTCREEHKKKHPDASVNFSEFSKKCSERWKTMSAKEKGKFEDMAKADKARYEREMKTYI", # P09429 aa9-79
}

SEEDS = [0, 1, 2]


def run_boltz_2chain(target_seq: str, bridge_seq: str, run_dir: Path, seed: int) -> float | None:
    run_dir.mkdir(parents=True, exist_ok=True)
    yaml = f"""sequences:
  - protein:
      id: A
      sequence: {target_seq}
      msa: empty
  - protein:
      id: B
      sequence: {bridge_seq}
      msa: empty
"""
    yaml_path = run_dir / "input.yaml"
    yaml_path.write_text(yaml)

    cmd = [str(BOLTZ_BIN), "predict", str(yaml_path),
           "--out_dir", str(run_dir),
           "--recycling_steps", "3", "--sampling_steps", "200",
           "--output_format", "pdb", "--no_kernels", "--seed", str(seed)]
    env = os.environ.copy()
    env["LD_LIBRARY_PATH"] = "/usr/lib/wsl/lib:" + env.get("LD_LIBRARY_PATH", "")
    env["XLA_PYTHON_CLIENT_PREALLOCATE"] = "false"
    r = subprocess.run(cmd, capture_output=True, text=True, env=env)
    if r.returncode != 0:
        print(f"    ERROR: {r.stderr[-200:]}")
        return None
    conf_files = list(run_dir.rglob("confidence_*.json"))
    if not conf_files:
        return None
    with open(conf_files[0]) as f:
        c = json.load(f)
    return c.get("iptm", c.get("protein_iptm"))


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    panel = {"SOX2": SOX2_HMG, "SOX9": SOX9_HMG}
    panel.update(OFF_TARGETS)

    rows = []
    total = len(panel) * len(SEEDS)
    done  = 0

    print(f"=== Bridge_48_s01 off-target selectivity screen ===")
    print(f"Bridge: {len(BRIDGE_SEQ)} aa")
    print(f"Panel: {list(panel.keys())}")
    print(f"Seeds: {SEEDS}\n")

    for target_name, target_seq in panel.items():
        iptms = []
        for seed in SEEDS:
            done += 1
            run_dir = OUT_DIR / f"{target_name}_s{seed}"
            print(f"[{done}/{total}] {target_name} ({len(target_seq)} aa) seed={seed}", flush=True)

            # Use cached result if available
            cached = list(run_dir.rglob("confidence_*.json"))
            if cached:
                with open(cached[0]) as f:
                    c = json.load(f)
                iptm = c.get("iptm", c.get("protein_iptm"))
                print(f"  (cached) ipTM={iptm:.3f}")
            else:
                iptm = run_boltz_2chain(target_seq, BRIDGE_SEQ, run_dir, seed)
                print(f"  ipTM={iptm:.3f}" if iptm is not None else "  FAILED")

            iptms.append(iptm)
            rows.append({
                "target": target_name,
                "target_len": len(target_seq),
                "seed": seed,
                "iptm": iptm,
                "is_on_target": target_name in ("SOX2", "SOX9"),
            })

        valid = [x for x in iptms if x is not None]
        if valid:
            mean_iptm = sum(valid) / len(valid)
            passes = mean_iptm < 0.40
            flag = "ON-TARGET" if target_name in ("SOX2","SOX9") else ("PASS (non-reactive)" if passes else "FAIL (cross-reactive!)")
            print(f"  → {target_name}: mean ipTM={mean_iptm:.3f}  [{flag}]\n")

    df = pd.DataFrame(rows)
    out_path = BASE / "analysis" / "offtarget_selectivity_results.tsv"
    df.to_csv(out_path, sep="\t", index=False)

    # Summary
    summary = df.groupby("target").agg(
        mean_iptm=("iptm","mean"), min_iptm=("iptm","min"), max_iptm=("iptm","max"),
        is_on_target=("is_on_target","first")
    ).reset_index().sort_values("mean_iptm", ascending=False)

    print(f"\n{'='*65}")
    print(f"  {'Target':<8} {'Mean ipTM':>10} {'Min':>7} {'Max':>7}  Verdict")
    print(f"  {'-'*63}")
    for _, row in summary.iterrows():
        verdict = ("ON-TARGET" if row["is_on_target"]
                   else ("PASS" if row["mean_iptm"] < 0.40 else "FAIL — CROSS-REACTIVE"))
        print(f"  {row['target']:<8} {row['mean_iptm']:>10.3f} {row['min_iptm']:>7.3f} {row['max_iptm']:>7.3f}  {verdict}")
    print(f"{'='*65}")
    print(f"\nResults saved to {out_path}")


if __name__ == "__main__":
    main()
