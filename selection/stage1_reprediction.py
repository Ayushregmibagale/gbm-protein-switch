#!/usr/bin/env python3
"""
Stage 1: Boltz-2 re-prediction of all designs from BindCraft, PH, BoltzGen.
Normalises scores across all three methods on a common scorer.
Output: stage1_results.tsv, ranked by ipTM.

Usage:
  python stage1_reprediction.py [--dry_run]
"""

from selection_support import require_execution, reference_sequence
require_execution(__name__)

import os, sys, json, shutil, subprocess, argparse
import pandas as pd
from pathlib import Path
from Bio.PDB import PDBParser, PPBuilder

# ── Target sequences (HMG domains from crystal structures) ──────────────────
SOX2_TARGET = "DRVKRPMNAFMVWSRGQRRKMAQENPKMHNSEISKRLGAEWKLLSETEKRPFIDEAKRLRALHMKEHPDYKYRPRRKTKTLMKKDKYTLPGGLLAPG"
SOX9_TARGET = "PHVKRPMNAFMVWAQAARRKLADQYPHLHNAELSKTLGKLWRLLNESEKRPFVEEAERLRVQHKKDHPDYKY"

# ── Thresholds ───────────────────────────────────────────────────────────────
IPTM_THRESH  = 0.70
PLDDT_THRESH = 70.0

# ── Paths ────────────────────────────────────────────────────────────────────
BASE      = Path(os.environ.get("GBM_ROOT", Path(__file__).resolve().parents[1]))
OUT_DIR   = BASE / "selection" / "stage1_boltz_outputs"
BOLTZ_BIN = Path(os.environ.get("BOLTZ_BIN", "boltz"))

pdb_parser = PDBParser(QUIET=True)
ppb        = PPBuilder()


def extract_binder_seq_pdb(pdb_path: Path) -> str:
    """Read explicit BindCraft binder chain B; target length must not decide identity."""
    from Bio.Data.PDBData import protein_letters_3to1
    struct = pdb_parser.get_structure("x", pdb_path)
    if "B" not in struct[0]:
        raise ValueError(f"{pdb_path}: expected binder chain B is absent")
    residues = [r for r in struct[0]["B"] if r.id[0] == " "]
    if not residues or any("CA" not in r or r.resname not in protein_letters_3to1 for r in residues):
        raise ValueError(f"{pdb_path}: empty, incomplete or noncanonical binder chain")
    return "".join(protein_letters_3to1[r.resname] for r in residues)


def collect_designs() -> list[dict]:
    designs = []

    # ── BindCraft SOX2 ──────────────────────────────────────────────────────
    bc_sox2 = BASE / "bindcraft/designs/sox2_binder_xtal/Accepted"
    for pdb in sorted(bc_sox2.glob("*.pdb")):
        try:
            seq = extract_binder_seq_pdb(pdb)
            designs.append({"name": pdb.stem, "source": "bindcraft", "target": "SOX2",
                            "binder_seq": seq, "pdb": str(pdb)})
        except Exception as e:
            print(f"  WARN {pdb.name}: {e}")

    # ── BindCraft SOX9 ──────────────────────────────────────────────────────
    bc_sox9 = BASE / "bindcraft/designs/sox9_binder_xtal_v3/Accepted"
    for pdb in sorted(bc_sox9.glob("*.pdb")):
        try:
            seq = extract_binder_seq_pdb(pdb)
            designs.append({"name": pdb.stem, "source": "bindcraft", "target": "SOX9",
                            "binder_seq": seq, "pdb": str(pdb)})
        except Exception as e:
            print(f"  WARN {pdb.name}: {e}")

    # ── Protein Hunter SOX2 — use best_seq from summary CSV ────────────────
    ph_sox2_csv = BASE / "protein_hunter/sox2_binder/summary_all_runs.csv"
    if ph_sox2_csv.exists():
        df = pd.read_csv(ph_sox2_csv)
        for _, row in df.iterrows():
            seq = str(row.get("best_seq", "")).strip()
            if seq and seq != "nan" and len(seq) > 10:
                designs.append({"name": f"ph_sox2_run{int(row['run_id'])}", "source": "ph",
                                "target": "SOX2", "binder_seq": seq, "pdb": None,
                                "ph_iptm": row.get("best_iptm", None)})

    # ── Protein Hunter SOX9 HMG — use best_seq from summary CSV ────────────
    ph_sox9_csv = BASE / "protein_hunter/sox9_binder_hmg/summary_all_runs.csv"
    if ph_sox9_csv.exists():
        df = pd.read_csv(ph_sox9_csv)
        for _, row in df.iterrows():
            seq = str(row.get("best_seq", "")).strip()
            if seq and seq != "nan" and len(seq) > 10:
                designs.append({"name": f"ph_sox9_run{int(row['run_id'])}", "source": "ph",
                                "target": "SOX9", "binder_seq": seq, "pdb": None,
                                "ph_iptm": row.get("best_iptm", None)})

    # ── BoltzGen SOX2 — read sequence from metrics CSV ─────────────────────
    bg_sox2_csv = BASE / "boltzgen/output/sox2/final_ranked_designs/final_designs_metrics_20.csv"
    if bg_sox2_csv.exists():
        df = pd.read_csv(bg_sox2_csv)
        for _, row in df.iterrows():
            seq = str(row.get("sequence", "")).strip()
            if seq and seq != "nan" and len(seq) > 10:
                designs.append({"name": f"boltzgen_sox2_{row['id']}", "source": "boltzgen",
                                "target": "SOX2", "binder_seq": seq, "pdb": None,
                                "bg_iptm": row.get("design_to_target_iptm", None)})

    # ── BoltzGen SOX9 — read sequence from metrics CSV ─────────────────────
    bg_sox9_csv = BASE / "boltzgen/output/sox9/final_ranked_designs/final_designs_metrics_20.csv"
    if bg_sox9_csv.exists():
        df = pd.read_csv(bg_sox9_csv)
        for _, row in df.iterrows():
            seq = str(row.get("sequence", "")).strip()
            if seq and seq != "nan" and len(seq) > 10:
                designs.append({"name": f"boltzgen_sox9_{row['id']}", "source": "boltzgen",
                                "target": "SOX9", "binder_seq": seq, "pdb": None,
                                "bg_iptm": row.get("design_to_target_iptm", None)})

    return designs


def write_boltz_yaml(name: str, binder_seq: str, target: str, yaml_path: Path):
    target_seq = SOX2_TARGET if target == "SOX2" else SOX9_TARGET
    content = f"""sequences:
  - protein:
      id: A
      sequence: {binder_seq}
      msa: empty
  - protein:
      id: B
      sequence: {target_seq}
      msa: empty
"""
    yaml_path.write_text(content)


def run_boltz(yaml_path: Path, out_dir: Path) -> dict | None:
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
        return None

    # Parse confidence JSON
    conf_files = list(out_dir.rglob("confidence_*.json"))
    if not conf_files:
        print(f"  ERROR no confidence JSON found in {out_dir}")
        return None
    with open(conf_files[0]) as f:
        conf = json.load(f)
    return conf


def main():
    global OUT_DIR
    parser = argparse.ArgumentParser()
    parser.add_argument("--dry_run", action="store_true", help="Collect designs only, don't run Boltz")
    parser.add_argument("--results", type=Path, help="New output table; historical results are preserved")
    parser.add_argument("--run-dir", type=Path, help="New inference output directory; prevents stale-cache reuse")
    args = parser.parse_args()

    if not args.dry_run:
        if args.results is None or args.run_dir is None:
            parser.error('--results and --run-dir are required for an inference campaign')
        if args.results.exists() or args.run_dir.exists():
            parser.error('Use new result and run directories to preserve history and avoid stale predictions')
        OUT_DIR = args.run_dir
    results_path = args.results

    print("=== Stage 1: Collecting designs ===")
    designs = collect_designs()
    print(f"  Total: {len(designs)} designs")
    for src in ["bindcraft", "ph", "boltzgen"]:
        n = sum(1 for d in designs if d["source"] == src)
        sox2 = sum(1 for d in designs if d["source"] == src and d["target"] == "SOX2")
        sox9 = sum(1 for d in designs if d["source"] == src and d["target"] == "SOX9")
        print(f"  {src}: {n} total (SOX2={sox2}, SOX9={sox9})")

    if args.dry_run:
        print("\nDry run — exiting before Boltz prediction.")
        return

    OUT_DIR.mkdir(parents=True)
    results_path.parent.mkdir(parents=True, exist_ok=True)
    print("\n=== Stage 1: Boltz-2 re-prediction ===")
    rows = []
    for i, d in enumerate(designs):
        name = d["name"]
        print(f"[{i+1}/{len(designs)}] {name} ({d['source']}, {d['target']})")

        yaml_dir = OUT_DIR / name
        yaml_dir.mkdir(exist_ok=True)
        yaml_path = yaml_dir / "input.yaml"
        write_boltz_yaml(name, d["binder_seq"], d["target"], yaml_path)

        conf = run_boltz(yaml_path, yaml_dir)
        if conf is None:
            rows.append({**d, "iptm": None, "plddt": None, "ptm": None, "pass_stage1": False})
            continue

        iptm  = conf.get("iptm",  conf.get("protein_iptm", None))
        plddt = conf.get("complex_plddt", conf.get("plddt", None))
        ptm   = conf.get("ptm", None)

        passed = (iptm  is not None and iptm  >= IPTM_THRESH and
                  plddt is not None and 0 <= plddt <= 1 and plddt * 100 >= PLDDT_THRESH)

        print(f"  ipTM={iptm}  pLDDT(fraction)={plddt}  pass={passed}")
        rows.append({**d, "iptm": iptm, "plddt": plddt, "plddt_percent": plddt * 100 if plddt is not None else None, "ptm": ptm, "pass_stage1": passed})

    df = pd.DataFrame(rows).sort_values("iptm", ascending=False)
    df.to_csv(results_path, sep="\t", index=False)

    passed = df[df["pass_stage1"] == True]
    print(f"\n=== Stage 1 complete ===")
    print(f"  Passed: {len(passed)}/{len(df)} (ipTM≥{IPTM_THRESH}, pLDDT≥{PLDDT_THRESH})")
    print(f"  Results: {results_path}")
    print(f"\nTop 10:")
    print(df[["name","source","target","iptm","plddt","pass_stage1"]].head(10).to_string(index=False))


if __name__ == "__main__":
    main()
