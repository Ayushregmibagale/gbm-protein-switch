#!/usr/bin/env python3
"""
3-chain selectivity screen for bridge_48_s01.

Proper assay for a bridge binder: replace ONE target at a time with an off-target
and measure whether the complex ipTM drops vs the on-target co-complex.

On-target baseline (cached, boltz_01):
  SOX2 + SOX9 + bridge  →  ipTM = 0.471

Tests (replace SOX2 with off-target):
  SOX11  + SOX9  + bridge   (SOX11 expressed in neurons — main concern)
  HMGB1  + SOX9  + bridge   (ubiquitous nuclear HMG)
  SRY    + SOX9  + bridge   (most divergent SOX HMG)

Tests (replace SOX9 with off-target):
  SOX2   + SOX11 + bridge
  SOX2   + HMGB1 + bridge
  SOX2   + SRY   + bridge

Historical heuristic: 3-chain ipTM < 0.35.
This is structural confidence, not a validated specificity or affinity threshold.

Seeds: 3 per condition to capture Boltz-2 stochasticity.

Output: analysis/bridge_3chain_selectivity_results.tsv
"""

import os, json, subprocess
import argparse
import pandas as pd
from pathlib import Path

BASE      = Path(__file__).resolve().parents[1]
OUT_DIR   = BASE / "analysis" / "bridge_3chain_selectivity_outputs"
BOLTZ_BIN = Path("/home/ayushregmibagale/miniconda3/envs/boltz2/bin/boltz")

BRIDGE_SEQ = "GAAAAAAAEAAAAAAAAAAAAAAAAAAAAEAARAAAAAAAAAAAAAAAAAAQAAAAQAYLAQAQAAAAA"
SOX2_HMG   = "DRVKRPMNAFMVWSRGQRRKMAQENPKMHNSEISKRLGAEWKLLSETEKRPFIDEAKRLRALHMKEHPDYKYRPRRKTKTLMKKDKYTLPGGLLAPG"
SOX9_HMG   = "PHVKRPMNAFMVWAQAARRKLADQYPHLHNAELSKTLGKLWRLLNESEKRPFVEEAERLRVQHKKDHPDYKY"

# Off-targets that flagged in 2-chain screen (HMG domain sequences)
SOX11_HMG  = "IKRPMNAFMVWSKIERRKIMEQSPDMHNAEISKRLGKRWKMLKDSEKIPFIREAERLRLKHMADYPDYK"   # P35716 aa49-117
HMGB1_HMG  = "PRGKMSSYAFFVQTCREEHKKKHPDASVNFSEFSKKCSERWKTMSAKEKGKFEDMAKADKARYEREMKTYI" # P09429 aa9-79
SRY_HMG    = "VKRPMNAFIVWSRDQRRKMALENPRMRNSEISKQLGYQWKMLTEAEKWPFFQEAQKLQAMHREKYPNYK"   # Q05066 aa60-128



# Each condition: (label, chain_A_seq, chain_B_seq)
# Chain A = position-1 target, Chain B = position-2 target, Chain C = bridge
CONDITIONS = [
    # On-target (will be re-run to confirm baseline with same seeds)
    ("SOX2+SOX9",   SOX2_HMG,  SOX9_HMG),
    # Replace chain A (SOX2 position) with off-target
    ("SOX11+SOX9",  SOX11_HMG, SOX9_HMG),
    ("HMGB1+SOX9",  HMGB1_HMG, SOX9_HMG),
    ("SRY+SOX9",    SRY_HMG,   SOX9_HMG),
    # Replace chain B (SOX9 position) with off-target
    ("SOX2+SOX11",  SOX2_HMG,  SOX11_HMG),
    ("SOX2+HMGB1",  SOX2_HMG,  HMGB1_HMG),
    ("SOX2+SRY",    SOX2_HMG,  SRY_HMG),
]

SEEDS = [0, 1, 2]
PASS_THRESHOLD = 0.35


def run_boltz_3chain(seqA: str, seqB: str, seqC: str, run_dir: Path, seed: int) -> float | None:
    run_dir.mkdir(parents=True, exist_ok=True)
    yaml = f"""sequences:
  - protein:
      id: A
      sequence: {seqA}
      msa: empty
  - protein:
      id: B
      sequence: {seqB}
      msa: empty
  - protein:
      id: C
      sequence: {seqC}
      msa: empty
"""
    (run_dir / "input.yaml").write_text(yaml)
    cmd = [str(BOLTZ_BIN), "predict", str(run_dir / "input.yaml"),
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
    conf = list(run_dir.rglob("confidence_*.json"))
    if not conf:
        return None
    with open(conf[0]) as f:
        c = json.load(f)
    return c.get("iptm", c.get("protein_iptm"))


def main():
    global BASE, OUT_DIR, BOLTZ_BIN
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--base', type=Path, default=BASE)
    parser.add_argument('--boltz', type=Path, default=Path(os.environ.get('BOLTZ_BIN', 'boltz')))
    args = parser.parse_args()
    BASE = args.base
    OUT_DIR = BASE / 'analysis' / 'bridge_3chain_selectivity_outputs'
    BOLTZ_BIN = args.boltz
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    print(f"=== Bridge_48_s01 — 3-chain selectivity screen ===")
    print("On-target baseline will be estimated from the same seeds as off-targets.")
    print(f"Descriptive confidence threshold: {PASS_THRESHOLD}; not a specificity test.\n")

    rows = []
    total = len(CONDITIONS) * len(SEEDS)
    done  = 0

    for label, seqA, seqB in CONDITIONS:
        iptms = []
        for seed in SEEDS:
            done += 1
            run_dir = OUT_DIR / f"{label.replace('+','_vs_')}_s{seed}"
            print(f"[{done}/{total}] {label} | seed {seed}", flush=True)

            cached = list(run_dir.rglob("confidence_*.json"))
            if cached:
                with open(cached[0]) as f:
                    c = json.load(f)
                iptm = c.get("iptm", c.get("protein_iptm"))
                print(f"  (cached) ipTM={iptm:.3f}")
            else:
                iptm = run_boltz_3chain(seqA, seqB, BRIDGE_SEQ, run_dir, seed)
                print(f"  ipTM={iptm:.3f}" if iptm is not None else "  FAILED")

            iptms.append(iptm)
            rows.append({"condition": label, "seed": seed, "iptm": iptm})

        valid = [x for x in iptms if x is not None]
        print(f"  {label}: completed {len(valid)}/{len(SEEDS)} seeds")

    df = pd.DataFrame(rows)
    out_path = BASE / "analysis" / "bridge_3chain_selectivity_results.tsv"
    df.to_csv(out_path, sep="\t", index=False)

    if df['iptm'].isna().any():
        raise RuntimeError(f"Incomplete screen; recorded partial results at {out_path}. No comparison issued.")
    baseline = float(df.loc[df['condition'] == 'SOX2+SOX9', 'iptm'].mean())
    print(f"Measured matched-seed on-target mean: {baseline:.3f}")
    print("ipTM is structural confidence; these comparisons do not establish specificity.")

    summary = df.groupby("condition").agg(
        mean_iptm=("iptm", "mean"), min_iptm=("iptm", "min"), max_iptm=("iptm", "max")
    ).reset_index().sort_values("mean_iptm", ascending=False)

    print(f"\n{'='*72}")
    print(f"  {'Condition':<15} {'Mean ipTM':>10} {'Min':>7} {'Max':>7}  {'Drop vs on-target':>18}  Verdict")
    print(f"  {'-'*70}")
    for _, row in summary.iterrows():
        is_on  = row["condition"] == "SOX2+SOX9"
        drop   = baseline - row["mean_iptm"]
        pct    = drop / baseline * 100 if baseline else float("nan")
        if is_on:
            verdict = "ON-TARGET"
            drop_str = "—"
        else:
            verdict  = "LOWER CONFIDENCE (heuristic)" if row["mean_iptm"] < PASS_THRESHOLD else "REVIEW STRUCTURE"
            drop_str = f"-{drop:.3f} ({pct:.0f}%)"
        print(f"  {row['condition']:<15} {row['mean_iptm']:>10.3f} {row['min_iptm']:>7.3f} "
              f"{row['max_iptm']:>7.3f}  {drop_str:>18}  {verdict}")
    print(f"{'='*72}")
    print(f"\nResults saved to {out_path}")


if __name__ == "__main__":
    main()
