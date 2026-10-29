#!/usr/bin/env python3
"""
Stage 5 flipped — ph_sox2_run18 × ph_sox9_run18 (Stage 4 rank-1 pair).

This pair scored highest in composite selectivity (SOX2 3.21×, SOX9 6.69×)
but was never run through Stage 5 caspase geometry. Running it now.

Flipped orientation:
  Chain C: ph_sox2_run18 + GGS(n) + CaspN (aa140-315)
  Chain D: CaspC (aa316-416) + GGS(n) + ph_sox9_run18

Results appended to: selection/stage5_ph_sox2_run18_results.tsv
"""

from selection_support import require_execution, reference_sequence, validate_recorded_candidates
require_execution(__name__)

import os, json, subprocess
import numpy as np, pandas as pd
from pathlib import Path
from Bio.PDB import PDBParser

SOX2_TARGET = "DRVKRPMNAFMVWSRGQRRKMAQENPKMHNSEISKRLGAEWKLLSETEKRPFIDEAKRLRALHMKEHPDYKYRPRRKTKTLMKKDKYTLPGGLLAPG"
SOX9_TARGET = "PHVKRPMNAFMVWAQAARRKLADQYPHLHNAELSKTLGKLWRLLNESEKRPFVEEAERLRVQHKKDHPDYKY"

BASE      = Path(os.environ.get("GBM_ROOT", Path(__file__).resolve().parents[1]))
OUT_DIR   = BASE / "selection" / "stage5_ph_sox2_run18_outputs"
BOLTZ_BIN = Path(os.environ.get("BOLTZ_BIN", "boltz"))

SOX2_NAME = "ph_sox2_run18"
SOX2_SEQ  = "SRLEELVDHIVDSMEEANQGGLTDEAREMIRNNLLEDLLDVNPNVEAELDAASDAEFAELLAGVDRIIREGTTSDLLGFLEALAEGGVAGALAELAK"
SOX9_NAME = "ph_sox9_run18"
SOX9_SEQ  = "MEEVRELLERARELEEKVRELEDDPSKREERERLLEEAERLREEALDLLVEQIREPLLQWMREEGLL"

LINKER_NS = [1, 2, 3]
SEEDS     = [0, 1, 2, 3, 4]   # 5 seeds — more than original to get stable hit rate

pdb_parser = PDBParser(QUIET=True)


def get_casp_halves() -> tuple[str, str]:
    cache = Path(__file__).resolve().parents[1] / "inputs/reference/P55211.fasta"
    seq = "".join(cache.read_text().strip().split("\n")[1:])
    return seq[139:315], seq[315:416]


def recon_dist_flipped(pdb_path: Path) -> float | None:
    """C-term(chain C) <-> N-term(chain D) distance."""
    try:
        s = pdb_parser.get_structure("x", pdb_path)
        m = s[0]
        c_res = list(m["C"].get_residues())
        d_res = list(m["D"].get_residues())
        ca_c = c_res[-1]["CA"].get_vector().get_array()
        ca_d = d_res[0]["CA"].get_vector().get_array()
        return float(np.linalg.norm(ca_c - ca_d))
    except Exception as e:
        print(f"  WARN: {e}")
        return None


def run_boltz(yaml_path: Path, out_dir: Path, seed: int) -> tuple[float | None, float | None]:
    cmd = [str(BOLTZ_BIN), "predict", str(yaml_path),
           "--out_dir", str(out_dir),
           "--recycling_steps", "3", "--sampling_steps", "200",
           "--output_format", "pdb", "--no_kernels", "--seed", str(seed)]
    env = os.environ.copy()
    env["LD_LIBRARY_PATH"] = "/usr/lib/wsl/lib:" + env.get("LD_LIBRARY_PATH", "")
    env["XLA_PYTHON_CLIENT_PREALLOCATE"] = "false"
    r = subprocess.run(cmd, capture_output=True, text=True, env=env)
    if r.returncode != 0:
        print(f"  ERROR: {r.stderr[-300:]}")
        return None, None
    conf_files = list(out_dir.rglob("confidence_*.json"))
    pdb_files  = list(out_dir.rglob("*.pdb"))
    if not conf_files: return None, None
    with open(conf_files[0]) as f: c = json.load(f)
    iptm = c.get("iptm", c.get("protein_iptm"))
    dist = recon_dist_flipped(pdb_files[0]) if pdb_files else None
    return iptm, dist


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    casp_n, casp_c = get_casp_halves()

    print(f"=== Stage 5 flipped: {SOX2_NAME} × {SOX9_NAME} ===")
    print(f"Stage 4 rank-1 pair (composite score 5.996, highest selectivity)")
    print(f"Chain C: {SOX2_NAME}({len(SOX2_SEQ)}aa) + GGS(n) + CaspN({len(casp_n)}aa)")
    print(f"Chain D: CaspC({len(casp_c)}aa) + GGS(n) + {SOX9_NAME}({len(SOX9_SEQ)}aa)")
    print(f"Linkers: GGS×{LINKER_NS}  Seeds: {SEEDS}")

    rows = []
    total = len(LINKER_NS) * len(SEEDS)
    done  = 0

    for n in LINKER_NS:
        linker  = "GGS" * n
        chain_c = SOX2_SEQ + linker + casp_n
        chain_d = casp_c + linker + SOX9_SEQ

        for seed in SEEDS:
            done += 1
            name    = f"ph_sox2_run18_ph_sox9_run18_GGS{n}_s{seed}"
            run_dir = OUT_DIR / name
            run_dir.mkdir(exist_ok=True)

            print(f"\n[{done}/{total}] GGS×{n} | seed {seed}")

            cached_conf = list(run_dir.rglob("confidence_*.json"))
            if cached_conf:
                cached_pdb = list(run_dir.rglob("*.pdb"))
                with open(cached_conf[0]) as f: c = json.load(f)
                iptm = c.get("iptm", c.get("protein_iptm"))
                dist = recon_dist_flipped(cached_pdb[0]) if cached_pdb else None
                status = f"recon={dist:.1f}Å  iptm={iptm:.3f}" if dist else "no PDB"
                print(f"  (cached) {status}")
                rows.append({"sox2": SOX2_NAME, "sox9": SOX9_NAME, "linker_n": n, "seed": seed,
                             "recon_dist_A": dist, "iptm": iptm,
                             "pass": dist is not None and dist < 15.0})
                continue

            yaml = f"""sequences:
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
"""
            (run_dir / "input.yaml").write_text(yaml)
            iptm, dist = run_boltz(run_dir / "input.yaml", run_dir, seed)
            passed = dist is not None and dist < 15.0
            print(f"  recon={dist:.1f}Å  iptm={iptm:.3f}  pass={passed}" if dist else "  FAILED")
            rows.append({"sox2": SOX2_NAME, "sox9": SOX9_NAME, "linker_n": n, "seed": seed,
                         "recon_dist_A": dist, "iptm": iptm, "pass": passed})

    df = pd.DataFrame(rows)
    out_path = BASE / "selection" / "stage5_ph_sox2_run18_results.tsv"
    df.to_csv(out_path, sep="\t", index=False)

    print(f"\n=== Results: {SOX2_NAME} × {SOX9_NAME} ===")
    summary = df.groupby("linker_n").agg(
        mean_dist=("recon_dist_A", "mean"), min_dist=("recon_dist_A", "min"),
        n_pass=("pass", "sum"), n_seeds=("pass", "count")
    ).reset_index()
    print(summary.to_string(index=False))
    print(f"\nTotal passing: {df['pass'].sum()}/{len(df)}")

    passing = df[df["pass"]]
    if len(passing):
        print(f"\nPassing seeds:")
        print(passing[["linker_n", "seed", "recon_dist_A", "iptm"]].to_string(index=False))
    else:
        print("\nNo seeds passed <15 Å.")


if __name__ == "__main__":
    main()
