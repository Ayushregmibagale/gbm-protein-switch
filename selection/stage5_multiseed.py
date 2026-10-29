#!/usr/bin/env python3
"""
Stage 5 multi-seed confirmation for lead construct:
  BC mpnn18 (SOX2_mini_l71_s878170_mpnn18_model2) × ph_sox9_run18, GGS×1

Runs seeds 1–4 (seed 0 already done in stage5_caspase_geometry.py).
Reports distribution of reconstitution distances.
Pass criterion: ≥3/5 seeds give recon dist < 15 Å.
"""

from selection_support import require_execution, reference_sequence, validate_recorded_candidates, campaign_output, campaign_table
require_execution(__name__)

import os, json, subprocess, urllib.request
import numpy as np, pandas as pd
from pathlib import Path
from Bio.PDB import PDBParser

BASE      = Path(os.environ.get("GBM_ROOT", Path(__file__).resolve().parents[1]))
OUT_DIR   = campaign_output(BASE / "selection" / "stage5_multiseed_outputs")
BOLTZ_BIN = Path(os.environ.get("BOLTZ_BIN", "boltz"))

SOX2_TARGET = "DRVKRPMNAFMVWSRGQRRKMAQENPKMHNSEISKRLGAEWKLLSETEKRPFIDEAKRLRALHMKEHPDYKYRPRRKTKTLMKKDKYTLPGGLLAPG"
SOX9_TARGET = "PHVKRPMNAFMVWAQAARRKLADQYPHLHNAELSKTLGKLWRLLNESEKRPFVEEAERLRVQHKKDHPDYKY"

# Lead construct sequences
SOX2_BINDER = "PYEIPEEDREWDEYSFRMVEVHERLMEARRELRLLAPEETEKKRELEEEIARAEKELEEVHKEGEAIERKM"
SOX9_BINDER = "MEEVRELLERARELEEKVRELEDDPSKREERERLLEEAERLREEALDLLVEQIREPLLQWMREEGLL"
LINKER      = "GGS"   # GGS×1

SEEDS = [1, 2, 3, 4]   # seed 0 already exists in stage5_caspase_outputs

pdb_parser = PDBParser(QUIET=True)


def fetch_casp9_halves() -> tuple[str, str]:
    seq = reference_sequence()
    return seq[139:315], seq[315:416]   # N-half (140-315), C-half (316-416)


def recon_distance(pdb_path: Path, binder_d_len: int, linker_len: int) -> float | None:
    casp_c_idx = binder_d_len + linker_len
    try:
        s = pdb_parser.get_structure("x", pdb_path)
        m = s[0]
        c_res = list(m["C"].get_residues())
        d_res = list(m["D"].get_residues())
        ca_c = c_res[-1]["CA"].get_vector().get_array()          # end of CaspN
        ca_d = d_res[casp_c_idx]["CA"].get_vector().get_array()  # start of CaspC
        return float(np.linalg.norm(ca_c - ca_d))
    except Exception as e:
        print(f"  WARN recon_distance: {e}")
        return None


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    print("=== Stage 5 multi-seed confirmation ===")
    print("  Lead: SOX2_mini_l71_s878170_mpnn18_model2 × ph_sox9_run18, GGS×1")

    casp_n, casp_c = fetch_casp9_halves()
    print(f"  CaspN: {len(casp_n)} aa | CaspC: {len(casp_c)} aa")

    chain_c_seq = SOX2_BINDER + LINKER + casp_n
    chain_d_seq = SOX9_BINDER + LINKER + casp_c
    print(f"  Chain C: {len(chain_c_seq)} aa | Chain D: {len(chain_d_seq)} aa")
    print(f"  CaspC starts at D[{len(SOX9_BINDER)+len(LINKER)}]")

    yaml_template = f"""sequences:
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
      sequence: {chain_c_seq}
      msa: empty
  - protein:
      id: D
      sequence: {chain_d_seq}
      msa: empty
"""

    rows = []

    # Include seed 0 result from stage5
    seed0_pdb = list((BASE / "selection" / "stage5_caspase_outputs" /
                      "casp_SOX2_mini_l71_s878170_mpnn18_model2__ph_sox9_run18__GGS1").rglob("*.pdb"))
    if seed0_pdb:
        d0 = recon_distance(seed0_pdb[0], len(SOX9_BINDER), len(LINKER))
        conf0 = list((BASE / "selection" / "stage5_caspase_outputs" /
                      "casp_SOX2_mini_l71_s878170_mpnn18_model2__ph_sox9_run18__GGS1").rglob("confidence_*.json"))
        iptm0 = None
        if conf0:
            with open(conf0[0]) as f: c = json.load(f)
            iptm0 = c.get("iptm", c.get("protein_iptm"))
        print(f"\n[seed 0 / cached] recon={d0:.1f}Å  iptm={iptm0:.3f}" if d0 else "  seed 0: no distance")
        rows.append({"seed": 0, "recon_dist_A": d0, "iptm": iptm0, "pass": d0 is not None and d0 < 15.0})

    for seed in SEEDS:
        name = f"lead_mpnn18_ph_sox9run18_GGS1_seed{seed}"
        run_dir = OUT_DIR / name
        run_dir.mkdir(exist_ok=True)
        yaml_path = run_dir / "input.yaml"
        yaml_path.write_text(yaml_template)

        print(f"\n[seed {seed}/4] {name}")

        # Skip if cached
        cached_conf = list(run_dir.rglob("confidence_*.json"))
        if cached_conf:
            cached_pdb = list(run_dir.rglob("*.pdb"))
            with open(cached_conf[0]) as f: c = json.load(f)
            iptm = c.get("iptm", c.get("protein_iptm"))
            dist = recon_distance(cached_pdb[0], len(SOX9_BINDER), len(LINKER)) if cached_pdb else None
            print(f"  (cached) recon={dist:.1f}Å  iptm={iptm:.3f}" if dist else "  cached, no dist")
            rows.append({"seed": seed, "recon_dist_A": dist, "iptm": iptm, "pass": dist is not None and dist < 15.0})
            continue

        cmd = [
            str(BOLTZ_BIN), "predict", str(yaml_path),
            "--out_dir", str(run_dir),
            "--recycling_steps", "3",
            "--sampling_steps", "200",
            "--output_format", "pdb",
            "--no_kernels",
            "--seed", str(seed),
        ]
        env = os.environ.copy()
        env["LD_LIBRARY_PATH"] = "/usr/lib/wsl/lib:" + env.get("LD_LIBRARY_PATH", "")
        env["XLA_PYTHON_CLIENT_PREALLOCATE"] = "false"
        result = subprocess.run(cmd, capture_output=True, text=True, env=env)
        if result.returncode != 0:
            print(f"  ERROR: {result.stderr[-300:]}")
            rows.append({"seed": seed, "recon_dist_A": None, "iptm": None, "pass": False})
            continue

        conf_files = list(run_dir.rglob("confidence_*.json"))
        pdb_files  = list(run_dir.rglob("*.pdb"))
        if not conf_files:
            rows.append({"seed": seed, "recon_dist_A": None, "iptm": None, "pass": False})
            continue

        with open(conf_files[0]) as f: c = json.load(f)
        iptm = c.get("iptm", c.get("protein_iptm"))
        dist = recon_distance(pdb_files[0], len(SOX9_BINDER), len(LINKER)) if pdb_files else None

        passed = dist is not None and dist < 15.0
        print(f"  recon={dist:.1f}Å  iptm={iptm:.3f}  pass={passed}" if dist else f"  iptm={iptm:.3f}  dist=None")
        rows.append({"seed": seed, "recon_dist_A": dist, "iptm": iptm, "pass": passed})

    df = pd.DataFrame(rows)
    out_path = campaign_table(BASE / "selection" / "stage5_multiseed_results.tsv")
    df.to_csv(out_path, sep="\t", index=False)

    n_pass = df["pass"].sum()
    print(f"\n=== Multi-seed confirmation ===")
    print(f"  Seeds passing <15 Å: {n_pass}/5")
    print(f"  Recon distances: {df['recon_dist_A'].dropna().tolist()}")
    print(f"  Mean ± SD: {df['recon_dist_A'].mean():.1f} ± {df['recon_dist_A'].std():.1f} Å")
    print(f"\n  {'CONFIRMED — proceed to MD' if n_pass >= 3 else 'NOT CONFIRMED — geometry unreliable'}")
    print(f"\n{df.to_string(index=False)}")


if __name__ == "__main__":
    main()
