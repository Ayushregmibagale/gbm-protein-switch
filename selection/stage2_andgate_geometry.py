#!/usr/bin/env python3
"""
Stage 2: AND-gate geometry filter.
For all Stage 1 passing designs, test pairs (SOX2 binder + SOX9 binder) in a
4-chain Boltz-2 prediction: SOX2:SOX9:BinderA:BinderB.
Filter: terminus(BinderA) ↔ terminus(BinderB) < 40 Å AND pair ipTM > 0.6.

Usage:
  python stage2_andgate_geometry.py
"""

from selection_support import require_execution, reference_sequence, validate_recorded_candidates, campaign_output, campaign_table
require_execution(__name__)

import os, json, subprocess, itertools
import numpy as np
import pandas as pd
from pathlib import Path
from Bio.PDB import PDBParser

SOX2_TARGET = "DRVKRPMNAFMVWSRGQRRKMAQENPKMHNSEISKRLGAEWKLLSETEKRPFIDEAKRLRALHMKEHPDYKYRPRRKTKTLMKKDKYTLPGGLLAPG"
SOX9_TARGET = "PHVKRPMNAFMVWAQAARRKLADQYPHLHNAELSKTLGKLWRLLNESEKRPFVEEAERLRVQHKKDHPDYKY"

BASE      = Path(os.environ.get("GBM_ROOT", Path(__file__).resolve().parents[1]))
STAGE1    = Path(os.environ.get("GBM_STAGE1_TABLE", BASE / "selection" / "stage1_results.tsv"))
OUT_DIR   = campaign_output(BASE / "selection" / "stage2_boltz_outputs")
BOLTZ_BIN = Path(os.environ.get("BOLTZ_BIN", "boltz"))

TERM_DIST_THRESH = 40.0   # Å — caspase reconstitution geometry
PAIR_IPTM_THRESH = 0.60   # pair-level ipTM

pdb_parser = PDBParser(QUIET=True)


def write_quartet_yaml(name: str, sox2_binder: str, sox9_binder: str, yaml_path: Path):
    """4-chain input: SOX2 HMG (A) + SOX9 HMG (B) + BinderSOX2 (C) + BinderSOX9 (D)."""
    content = f"""sequences:
  - protein:
      id: A
      sequence: {SOX2_TARGET}
      msa: empty
  - protein:
      id: B
      sequence: {SOX9_TARGET}
      msa: empty
  - protein:
      id: C
      sequence: {sox2_binder}
      msa: empty
  - protein:
      id: D
      sequence: {sox9_binder}
      msa: empty
"""
    yaml_path.write_text(content)


def terminus_distance(pdb_path: Path, chain_c: str = "C", chain_d: str = "D") -> float | None:
    """Compute C-alpha distance between C-terminus of chain C and N-terminus of chain D."""
    try:
        struct = pdb_parser.get_structure("x", pdb_path)
        model = struct[0]
        c_residues = list(model[chain_c].get_residues())
        d_residues = list(model[chain_d].get_residues())
        if not c_residues or not d_residues:
            return None
        ca_c = c_residues[-1]["CA"].get_vector().get_array()  # C-terminus of binder C
        ca_d = d_residues[0]["CA"].get_vector().get_array()   # N-terminus of binder D
        return float(np.linalg.norm(ca_c - ca_d))
    except Exception as e:
        print(f"  WARN terminus_distance: {e}")
        return None


def run_boltz_quartet(yaml_path: Path, out_dir: Path) -> tuple[dict | None, Path | None]:
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
        print(f"  ERROR boltz failed:\n{result.stderr[-500:]}")
        return None, None

    conf_files = list(out_dir.rglob("confidence_*.json"))
    pdb_files  = list(out_dir.rglob("*.pdb"))
    if not conf_files:
        return None, None

    with open(conf_files[0]) as f:
        conf = json.load(f)
    pdb_file = pdb_files[0] if pdb_files else None
    return conf, pdb_file


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    results_path = campaign_table(BASE / "selection" / "stage2_results.tsv")

    df1 = pd.read_csv(STAGE1, sep="\t")
    validate_recorded_candidates(df1, BASE)
    sox2 = df1[(df1["target"] == "SOX2") & (df1["pass_stage1"] == True)]
    sox9 = df1[(df1["target"] == "SOX9") & (df1["pass_stage1"] == True)]

    print(f"=== Stage 2: AND-gate geometry filter ===")
    print(f"  SOX2 designs passing Stage 1: {len(sox2)}")
    print(f"  SOX9 designs passing Stage 1: {len(sox9)}")
    print(f"  Total pairs to test: {len(sox2) * len(sox9)}")

    rows = []
    for i, (_, s2) in enumerate(sox2.iterrows()):
        for j, (_, s9) in enumerate(sox9.iterrows()):
            pair_name = f"{s2['name']}__x__{s9['name']}"
            print(f"\n[{i*len(sox9)+j+1}/{len(sox2)*len(sox9)}] {pair_name}")

            pair_dir = OUT_DIR / pair_name
            pair_dir.mkdir(exist_ok=True)
            yaml_path = pair_dir / "input.yaml"

            write_quartet_yaml(pair_name, s2["binder_seq"], s9["binder_seq"], yaml_path)
            conf, pdb_file = run_boltz_quartet(yaml_path, pair_dir)

            if conf is None:
                rows.append({"pair": pair_name, "sox2_name": s2["name"], "sox9_name": s9["name"],
                            "sox2_source": s2["source"], "sox9_source": s9["source"],
                            "pair_iptm": None, "terminus_dist_A": None, "pass_stage2": False})
                continue

            pair_iptm = conf.get("iptm", conf.get("protein_iptm", None))
            term_dist = terminus_distance(pdb_file) if pdb_file else None

            passed = (pair_iptm is not None and pair_iptm >= PAIR_IPTM_THRESH and
                      term_dist is not None and term_dist <= TERM_DIST_THRESH)

            print(f"  pair_ipTM={pair_iptm:.3f}  term_dist={term_dist:.1f}Å  pass={passed}")
            rows.append({
                "pair": pair_name,
                "sox2_name": s2["name"], "sox9_name": s9["name"],
                "sox2_source": s2["source"], "sox9_source": s9["source"],
                "sox2_iptm": s2["iptm"], "sox9_iptm": s9["iptm"],
                "pair_iptm": pair_iptm,
                "terminus_dist_A": term_dist,
                "pass_stage2": passed,
                "pdb": str(pdb_file) if pdb_file else None,
            })

    df2 = pd.DataFrame(rows).sort_values("pair_iptm", ascending=False)
    df2.to_csv(results_path, sep="\t", index=False)

    passed = df2[df2["pass_stage2"] == True]
    print(f"\n=== Stage 2 complete ===")
    print(f"  Passed: {len(passed)}/{len(df2)} pairs")
    print(f"  Results: {results_path}")
    if len(passed) > 0:
        print("\nTop pairs:")
        print(passed[["pair","sox2_source","sox9_source","pair_iptm","terminus_dist_A"]].head(10).to_string(index=False))


if __name__ == "__main__":
    main()
