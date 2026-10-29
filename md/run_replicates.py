#!/usr/bin/env python3
"""
Run MD replicates 2 and 3 for the GBM lead construct.

Each replicate starts from the same equilibrated npt.gro (rep1) but
re-assigns Maxwell-Boltzmann velocities at 300 K with a different random seed.
This is the standard approach for GROMACS replicates from a shared equilibrated start.

Usage:
  python run_replicates.py          # run/resume all replicates
  python run_replicates.py --status # show progress
"""

import os, sys, json, shutil, subprocess, argparse, shlex
import numpy as np
from pathlib import Path
from workflow_support import read_state, write_state, run_gmx, production, periodic_analysis

BASE     = Path(os.environ.get("GBM_ROOT", Path(__file__).resolve().parents[1]))
REP1_DIR = BASE / "md" / "lead_construct"
MD_DIR   = BASE / "md"
MDP_DIR  = BASE / "md" / "mdp"
GMX      = os.environ.get("GMX_BIN", "gmx")
PYTHON   = "/home/ayushregmibagale/.conda/envs/gromacs/bin/python3"

REPLICATES = [
    {"id": "rep2", "seed": 12345},
    {"id": "rep3", "seed": 98765},
]

CHAIN_C_LAST_RESNUM  = 419   # CaspN C-term (chain C, after renumbering)
CHAIN_D_FIRST_RESNUM = 420   # CaspC N-term (chain D, after renumbering)


def get_state(rep_dir):
    return read_state(rep_dir)



def set_state(rep_dir, s):
    write_state(rep_dir, s)



def gmx_run(args, log, cwd, stdin=None):
    run_gmx(GMX, args, log, cwd, stdin)



def setup_rep(rep_dir: Path, seed: int):
    """Set up a replicate directory from rep1 equilibration files."""
    rep_dir.mkdir(parents=True, exist_ok=True)

    # Copy topology and equilibration files from rep1
    for f in ["topol.top", "npt.gro", "start.pdb"]:
        shutil.copy2(REP1_DIR / f, rep_dir / f)
    for f in REP1_DIR.glob("topol_*.itp"):
        shutil.copy2(f, rep_dir / f.name)
    for f in REP1_DIR.glob("posre_*.itp"):
        shutil.copy2(f, rep_dir / f.name)

    # Write production MDP with this replicate's seed
    mdp_template = (MDP_DIR / "prod_20ns_newvel.mdp").read_text()
    mdp_text = mdp_template.replace("SEED_PLACEHOLDER", str(seed))
    (rep_dir / "prod_newvel.mdp").write_text(mdp_text)

    set_state(rep_dir, "SETUP_DONE")


def run_prod(rep_dir):
    prepare = ['grompp', '-f', 'prod_newvel.mdp', '-c', 'npt.gro', '-p', 'topol.top', '-o', 'prod.tpr']
    production(GMX, rep_dir, prepare, 'production_checked')
    set_state(rep_dir, 'PROD_DONE')



def run_analysis(rep_dir):
    metrics = periodic_analysis(rep_dir, 'resnum 419 and name CA', 'resnum 420 and name CA', rep_dir.name, 'near')
    temporary = rep_dir / 'metrics.json.tmp'
    temporary.write_text(json.dumps(metrics, indent=2) + '\n')
    temporary.replace(rep_dir / 'metrics.json')
    set_state(rep_dir, 'ANALYSIS_DONE')



def print_summary():
    for directory in [REP1_DIR] + [MD_DIR / row['id'] for row in REPLICATES]:
        path = directory / 'metrics.json'
        if path.exists():
            metrics = json.loads(path.read_text())
            method = 'validated periodic' if 'topology_sha256' in metrics else 'LEGACY: requires periodic reanalysis'
            print(directory.name, method, json.dumps(metrics))



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--status", action="store_true")
    parser.add_argument('--execute', action='store_true', help='Run or resume stages; otherwise only report status')
    args = parser.parse_args()
    if not args.execute:
        args.status = True

    if args.status:
        for r in REPLICATES:
            rep_dir = MD_DIR / r["id"]
            print(f"{r['id']}: {get_state(rep_dir)}")
        print_summary()
        return

    # Run replicates sequentially (single GPU)
    for r in REPLICATES:
        rep_id  = r["id"]
        seed    = r["seed"]
        rep_dir = MD_DIR / rep_id
        state   = get_state(rep_dir)

        print(f"\n{'='*40}")
        print(f"  Replicate {rep_id} (seed={seed})  state={state}")

        if state == "PENDING":
            print("  Setting up from rep1 equilibration...")
            setup_rep(rep_dir, seed)
            state = "SETUP_DONE"

        if state == "SETUP_DONE":
            print("  Running production MD (20 ns)...")
            run_prod(rep_dir)
            state = "PROD_DONE"

        if state == "PROD_DONE":
            print("  Running analysis...")
            run_analysis(rep_dir)

        if get_state(rep_dir) == "ANALYSIS_DONE":
            print(f"  {rep_id} complete.")

    print_summary()


if __name__ == "__main__":
    main()
