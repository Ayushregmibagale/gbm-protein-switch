#!/usr/bin/env python3
"""
Stage 5: Caspase-9 geometry modelling.
Full 4-chain Boltz-2 assembly:
  A = SOX2 HMG, B = SOX9 HMG
  C = BinderSOX2 + GGS(n) + Caspase9_N-half (P55211 aa 140-315)
  D = BinderSOX9 + GGS(n) + Caspase9_C-half (P55211 aa 316-416)

Linker sweep: n = 1,2,3,4 (3,6,9,12 aa linkers).
Checks: inter-half distance < 15 Å, pair_iptm of full assembly.
"""

from selection_support import require_execution, reference_sequence, validate_recorded_candidates, campaign_output, campaign_table
require_execution(__name__)

import os, json, subprocess, urllib.request
import numpy as np, pandas as pd
from pathlib import Path
from Bio.PDB import PDBParser

SOX2_TARGET = "DRVKRPMNAFMVWSRGQRRKMAQENPKMHNSEISKRLGAEWKLLSETEKRPFIDEAKRLRALHMKEHPDYKYRPRRKTKTLMKKDKYTLPGGLLAPG"
SOX9_TARGET = "PHVKRPMNAFMVWAQAARRKLADQYPHLHNAELSKTLGKLWRLLNESEKRPFVEEAERLRVQHKKDHPDYKY"

BASE      = Path(os.environ.get("GBM_ROOT", Path(__file__).resolve().parents[1]))
STAGE1    = Path(os.environ.get("GBM_STAGE1_TABLE", BASE / "selection" / "stage1_results.tsv"))
STAGE4    = BASE / "selection" / "stage4_results.tsv"
OUT_DIR   = campaign_output(BASE / "selection" / "stage5_caspase_outputs")
BOLTZ_BIN = Path(os.environ.get("BOLTZ_BIN", "boltz"))

INTER_HALF_THRESH = 15.0   # Å — caspase half proximity for reconstitution
LINKER_REPEATS    = [1, 2, 3, 4]   # GGS × n

pdb_parser = PDBParser(QUIET=True)


def fetch_casp9_halves() -> tuple[str, str]:
    """Fetch P55211 from UniProt and extract N-half (140-315) and C-half (316-416)."""
    seq = reference_sequence()
    print(f"  P55211 length: {len(seq)} aa")
    # aa 140-315 (1-based → 0-based: 139:315)
    n_half = seq[139:315]
    # aa 316-416 (1-based → 0-based: 315:416)
    c_half = seq[315:416]
    print(f"  N-half (140-315): {len(n_half)} aa")
    print(f"  C-half (316-416): {len(c_half)} aa")
    return n_half, c_half


def make_linker(n: int) -> str:
    return "GGS" * n


def write_assembly_yaml(name: str, sox2_binder: str, sox9_binder: str,
                        casp_n: str, casp_c: str, linker: str,
                        yaml_path: Path):
    chain_c = sox2_binder + linker + casp_n
    chain_d = sox9_binder + linker + casp_c
    yaml_path.write_text(f"""sequences:
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
      sequence: {chain_c}
      msa: empty
  - protein:
      id: D
      sequence: {chain_d}
      msa: empty
""")


def run_boltz(yaml_path: Path, out_dir: Path) -> tuple[dict | None, Path | None]:
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
        print(f"  ERROR boltz:\n{result.stderr[-400:]}")
        return None, None
    conf_files = list(out_dir.rglob("confidence_*.json"))
    pdb_files  = list(out_dir.rglob("*.pdb"))
    if not conf_files:
        return None, None
    with open(conf_files[0]) as f:
        conf = json.load(f)
    return conf, pdb_files[0] if pdb_files else None


def inter_half_distance(pdb_path: Path, chain_c_id: str = "C", chain_d_id: str = "D",
                        casp_n_offset: int = 0, casp_c_offset: int = 0) -> float | None:
    """
    Distance between C-terminus of CaspaseN (end of chain C) and
    N-terminus of CaspaseC (start of chain D after linker).
    We use last residue of chain C and first residue of chain D as proxy
    (the linker-caspase junction geometry).
    """
    try:
        struct = pdb_parser.get_structure("x", pdb_path)
        model = struct[0]
        c_res = list(model[chain_c_id].get_residues())
        d_res = list(model[chain_d_id].get_residues())
        if not c_res or not d_res:
            return None
        ca_c = c_res[-1]["CA"].get_vector().get_array()   # C-terminus of chain C (end of CaspN)
        ca_d = d_res[0]["CA"].get_vector().get_array()    # N-terminus of chain D (start of CaspC)
        return float(np.linalg.norm(ca_c - ca_d))
    except Exception as e:
        print(f"  WARN inter_half_distance: {e}")
        return None


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    print("=== Stage 5: Caspase geometry modelling ===")
    casp_n, casp_c = fetch_casp9_halves()

    df4 = pd.read_csv(STAGE4, sep="\t")
    df1 = pd.read_csv(STAGE1, sep="\t")
    validate_recorded_candidates(df1, BASE)
    seq_lookup = dict(zip(df1["name"], df1["binder_seq"]))

    top5 = df4.head(5)
    rows = []

    total = len(top5) * len(LINKER_REPEATS)
    done = 0

    for _, pair in top5.iterrows():
        s2_name = pair["sox2_name"]
        s9_name = pair["sox9_name"]
        s2_seq  = seq_lookup.get(s2_name, "")
        s9_seq  = seq_lookup.get(s9_name, "")

        if not s2_seq or not s9_seq:
            print(f"  WARN: missing sequence for {s2_name} or {s9_name}")
            continue

        print(f"\n{'='*60}")
        print(f"Pair: {s2_name} × {s9_name}")
        print(f"  Stage4 score={pair['composite_score']:.4f}  pair_ipTM={pair['pair_iptm']:.3f}  term={pair['terminus_dist_A']:.1f}Å")

        for n in LINKER_REPEATS:
            linker = make_linker(n)
            name = f"casp_{s2_name}__{s9_name}__GGS{n}"
            done += 1
            print(f"\n[{done}/{total}] {name}")
            print(f"  Linker: GGS×{n} ({len(linker)} aa)")
            print(f"  Chain C length: {len(s2_seq)+len(linker)+len(casp_n)} aa")
            print(f"  Chain D length: {len(s9_seq)+len(linker)+len(casp_c)} aa")

            run_dir = OUT_DIR / name
            # Skip if already done
            if list(run_dir.rglob("confidence_*.json")):
                print("  (cached)")
                conf_files = list(run_dir.rglob("confidence_*.json"))
                pdb_files  = list(run_dir.rglob("*.pdb"))
                with open(conf_files[0]) as f:
                    conf = json.load(f)
                pdb_file = pdb_files[0] if pdb_files else None
            else:
                run_dir.mkdir(exist_ok=True)
                yaml_path = run_dir / "input.yaml"
                write_assembly_yaml(name, s2_seq, s9_seq, casp_n, casp_c, linker, yaml_path)
                conf, pdb_file = run_boltz(yaml_path, run_dir)

            if conf is None:
                rows.append({"name": name, "sox2_name": s2_name, "sox9_name": s9_name,
                             "linker_n": n, "pair_iptm": None, "inter_half_dist_A": None,
                             "pass_geometry": False})
                continue

            pair_iptm = conf.get("iptm", conf.get("protein_iptm", None))
            inter_dist = inter_half_distance(pdb_file) if pdb_file else None

            passed = (pair_iptm is not None and
                      inter_dist is not None and inter_dist <= INTER_HALF_THRESH)

            print(f"  pair_ipTM={pair_iptm:.3f}  inter_half={inter_dist:.1f}Å  pass={passed}" if inter_dist else
                  f"  pair_ipTM={pair_iptm:.3f}  inter_half=None  pass=False")

            rows.append({
                "name": name,
                "sox2_name": s2_name, "sox9_name": s9_name,
                "linker_n": n,
                "pair_iptm": pair_iptm,
                "inter_half_dist_A": inter_dist,
                "pass_geometry": passed,
                "pdb": str(pdb_file) if pdb_file else None,
            })

    df5 = pd.DataFrame(rows).sort_values(["pass_geometry", "inter_half_dist_A"],
                                          ascending=[False, True])
    out_path = campaign_table(BASE / "selection" / "stage5_results.tsv")
    df5.to_csv(out_path, sep="\t", index=False)

    passed = df5[df5["pass_geometry"] == True]
    print(f"\n=== Stage 5 complete ===")
    print(f"  Passed (inter-half < {INTER_HALF_THRESH}Å): {len(passed)}/{len(df5)}")
    print(f"  Results: {out_path}")
    if len(passed):
        print("\nPassing assemblies:")
        print(passed[["name","linker_n","pair_iptm","inter_half_dist_A"]].to_string(index=False))
    else:
        print("\nAll assemblies — inter-half distances:")
        print(df5[["name","linker_n","pair_iptm","inter_half_dist_A","pass_geometry"]].to_string(index=False))


if __name__ == "__main__":
    main()
