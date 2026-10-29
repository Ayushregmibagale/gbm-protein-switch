#!/usr/bin/env python3
"""
Stage 3: Selectivity filter.
For every unique binder in Stage 2 passing pairs, predict against the WRONG target.
A selective binder should have: correct_iptm >> wrong_iptm (wrong < 0.40).

SOX2 binders tested vs SOX9 HMG.
SOX9 binders tested vs SOX2 HMG.

Output: stage3_results.tsv — one row per pair, with selectivity annotations.
"""

import os, json, subprocess
import pandas as pd
from pathlib import Path

SOX2_TARGET = "DRVKRPMNAFMVWSRGQRRKMAQENPKMHNSEISKRLGAEWKLLSETEKRPFIDEAKRLRALHMKEHPDYKYRPRRKTKTLMKKDKYTLPGGLLAPG"
SOX9_TARGET = "PHVKRPMNAFMVWAQAARRKLADQYPHLHNAELSKTLGKLWRLLNESEKRPFVEEAERLRVQHKKDHPDYKY"

BASE      = Path("/mnt/d/gbm-protein-switch")
STAGE1    = BASE / "selection" / "stage1_results.tsv"
STAGE2    = BASE / "selection" / "stage2_results.tsv"
OUT_DIR   = BASE / "selection" / "stage3_boltz_outputs"
BOLTZ_BIN = Path("/home/ayushregmibagale/miniconda3/envs/boltz2/bin/boltz")

WRONG_IPTM_MAX = 0.40   # binder vs wrong target must be below this


def run_boltz(binder_seq: str, target_seq: str, name: str, out_dir: Path) -> float | None:
    yaml_path = out_dir / "input.yaml"
    yaml_path.write_text(f"""sequences:
  - protein:
      id: A
      sequence: {binder_seq}
      msa: empty
  - protein:
      id: B
      sequence: {target_seq}
      msa: empty
""")
    cmd = [
        str(BOLTZ_BIN), "predict", str(yaml_path),
        "--out_dir", str(out_dir),
        "--recycling_steps", "3",
        "--sampling_steps", "200",
        "--output_format", "pdb",
        "--no_kernels",
    ]
    env = os.environ.copy()
    env["LD_LIBRARY_PATH"] = "/usr/lib/wsl/lib:" + env.get("LD_LIBRARY_PATH", "")
    env["XLA_PYTHON_CLIENT_PREALLOCATE"] = "false"
    result = subprocess.run(cmd, capture_output=True, text=True, env=env)
    if result.returncode != 0:
        print(f"  ERROR boltz failed for {name}:\n{result.stderr[-300:]}")
        return None
    conf_files = list(out_dir.rglob("confidence_*.json"))
    if not conf_files:
        return None
    with open(conf_files[0]) as f:
        conf = json.load(f)
    return conf.get("iptm", conf.get("protein_iptm", None))


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    df1 = pd.read_csv(STAGE1, sep="\t")
    df2 = pd.read_csv(STAGE2, sep="\t")
    passed = df2[df2["pass_stage2"] == True].copy()

    # Build name → binder_seq lookup from stage1
    seq_lookup = dict(zip(df1["name"], df1["binder_seq"]))
    # Also build name → correct_iptm
    iptm_lookup = dict(zip(df1["name"], df1["iptm"]))

    sox2_binders = sorted(passed["sox2_name"].unique())
    sox9_binders = sorted(passed["sox9_name"].unique())
    print(f"=== Stage 3: Selectivity ===")
    print(f"  SOX2 binders to test vs SOX9: {len(sox2_binders)}")
    print(f"  SOX9 binders to test vs SOX2: {len(sox9_binders)}")

    wrong_iptm = {}   # name → wrong-target ipTM

    # SOX2 binders vs SOX9 (wrong target)
    for i, name in enumerate(sox2_binders):
        seq = seq_lookup.get(name)
        if not seq:
            print(f"  WARN: no seq for {name}")
            wrong_iptm[name] = None
            continue
        d = OUT_DIR / f"{name}__vs__SOX9_wrong"
        d.mkdir(exist_ok=True)
        # Skip if already done
        if list(d.rglob("confidence_*.json")):
            with open(list(d.rglob("confidence_*.json"))[0]) as f:
                c = json.load(f)
            iptm = c.get("iptm", c.get("protein_iptm", None))
            print(f"[{i+1}/{len(sox2_binders)}] {name} vs SOX9 (cached): ipTM={iptm:.3f}")
            wrong_iptm[name] = iptm
            continue
        print(f"[{i+1}/{len(sox2_binders)}] {name} vs SOX9 (wrong target)...")
        iptm = run_boltz(seq, SOX9_TARGET, name, d)
        print(f"  wrong-target ipTM = {iptm:.3f}" if iptm else "  FAILED")
        wrong_iptm[name] = iptm

    # SOX9 binders vs SOX2 (wrong target)
    for i, name in enumerate(sox9_binders):
        seq = seq_lookup.get(name)
        if not seq:
            print(f"  WARN: no seq for {name}")
            wrong_iptm[name] = None
            continue
        d = OUT_DIR / f"{name}__vs__SOX2_wrong"
        d.mkdir(exist_ok=True)
        if list(d.rglob("confidence_*.json")):
            with open(list(d.rglob("confidence_*.json"))[0]) as f:
                c = json.load(f)
            iptm = c.get("iptm", c.get("protein_iptm", None))
            print(f"[{i+1}/{len(sox9_binders)}] {name} vs SOX2 (cached): ipTM={iptm:.3f}")
            wrong_iptm[name] = iptm
            continue
        print(f"[{i+1}/{len(sox9_binders)}] {name} vs SOX2 (wrong target)...")
        iptm = run_boltz(seq, SOX2_TARGET, name, d)
        print(f"  wrong-target ipTM = {iptm:.3f}" if iptm else "  FAILED")
        wrong_iptm[name] = iptm

    # Build output table
    rows = []
    for _, row in passed.iterrows():
        s2_name = row["sox2_name"]
        s9_name = row["sox9_name"]
        s2_correct = iptm_lookup.get(s2_name)
        s9_correct = iptm_lookup.get(s9_name)
        s2_wrong = wrong_iptm.get(s2_name)
        s9_wrong = wrong_iptm.get(s9_name)

        s2_selective = s2_wrong is not None and s2_wrong < WRONG_IPTM_MAX
        s9_selective = s9_wrong is not None and s9_wrong < WRONG_IPTM_MAX
        both_selective = s2_selective and s9_selective

        rows.append({
            "pair": row["pair"],
            "sox2_name": s2_name,
            "sox9_name": s9_name,
            "sox2_source": row["sox2_source"],
            "sox9_source": row["sox9_source"],
            "pair_iptm": row["pair_iptm"],
            "terminus_dist_A": row["terminus_dist_A"],
            "sox2_correct_iptm": s2_correct,
            "sox2_wrong_iptm": s2_wrong,
            "sox2_selective": s2_selective,
            "sox9_correct_iptm": s9_correct,
            "sox9_wrong_iptm": s9_wrong,
            "sox9_selective": s9_selective,
            "both_selective": both_selective,
        })

    df3 = pd.DataFrame(rows).sort_values("pair_iptm", ascending=False)
    out_path = BASE / "selection" / "stage3_results.tsv"
    df3.to_csv(out_path, sep="\t", index=False)

    sel = df3[df3["both_selective"] == True]
    print(f"\n=== Stage 3 complete ===")
    print(f"  Both selective: {len(sel)}/{len(df3)} pairs")
    print(f"  Results: {out_path}")
    if len(sel):
        print("\nSelective pairs:")
        cols = ["pair","sox2_source","sox9_source","pair_iptm","terminus_dist_A",
                "sox2_wrong_iptm","sox9_wrong_iptm"]
        print(sel[cols].to_string(index=False))
    else:
        print("\nNo pairs fully selective — showing selectivity breakdown:")
        cols = ["sox2_name","sox2_correct_iptm","sox2_wrong_iptm","sox2_selective",
                "sox9_name","sox9_correct_iptm","sox9_wrong_iptm","sox9_selective"]
        print(df3[cols].drop_duplicates(subset=["sox2_name","sox9_name"]).head(15).to_string(index=False))


if __name__ == "__main__":
    main()
