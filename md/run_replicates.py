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

import os, sys, json, shutil, subprocess, argparse
import numpy as np
from pathlib import Path

BASE     = Path("/mnt/d/gbm-protein-switch")
REP1_DIR = BASE / "md" / "lead_construct"
MD_DIR   = BASE / "md"
MDP_DIR  = BASE / "md" / "mdp"
GMX      = "/home/ayushregmibagale/.conda/envs/gromacs/bin/gmx"
PYTHON   = "/home/ayushregmibagale/.conda/envs/gromacs/bin/python3"

REPLICATES = [
    {"id": "rep2", "seed": 12345},
    {"id": "rep3", "seed": 98765},
]

CHAIN_C_LAST_RESNUM  = 419   # CaspN C-term (chain C, after renumbering)
CHAIN_D_FIRST_RESNUM = 420   # CaspC N-term (chain D, after renumbering)


def get_state(rep_dir: Path) -> str:
    s = rep_dir / "state.txt"
    return s.read_text().strip() if s.exists() else "PENDING"


def set_state(rep_dir: Path, s: str):
    (rep_dir / "state.txt").write_text(s + "\n")
    print(f"  → {s}")


def gmx_run(args: str, log: Path, cwd: Path, stdin: str = None):
    env = os.environ.copy()
    env["GMX_MAXBACKUP"] = "-1"
    with open(log, "a") as f:
        r = subprocess.run(f"{GMX} {args}", shell=True, cwd=cwd, env=env,
                           input=stdin, text=True, stdout=f, stderr=subprocess.STDOUT)
    if r.returncode != 0:
        print(f"  FAILED (rc={r.returncode}) — see {log.name}")
        sys.exit(1)


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


def run_prod(rep_dir: Path):
    """Run grompp + mdrun for production."""
    print("  grompp...")
    gmx_run(f"grompp -f prod_newvel.mdp -c npt.gro -p topol.top -o prod.tpr -maxwarn 2",
            log=rep_dir / "grompp_prod.log", cwd=rep_dir)
    print("  mdrun (20 ns)...")
    gmx_run(f"mdrun -v -deffnm prod -ntmpi 1 -ntomp 8 -gpu_id 0 -nb gpu -pme gpu -bonded gpu",
            log=rep_dir / "mdrun_prod.log", cwd=rep_dir)
    set_state(rep_dir, "PROD_DONE")


def run_analysis(rep_dir: Path):
    """MDAnalysis: inter-half distance across production frames."""
    try:
        import MDAnalysis as mda
    except ImportError:
        print("  MDAnalysis not found")
        sys.exit(1)

    u = mda.Universe(str(rep_dir / "prod.tpr"), str(rep_dir / "prod.xtc"))

    c_term = u.select_atoms(f"resnum {CHAIN_C_LAST_RESNUM} and name CA")
    n_term = u.select_atoms(f"resnum {CHAIN_D_FIRST_RESNUM} and name CA")

    if len(c_term) != 1 or len(n_term) != 1:
        print(f"  ERROR: c_term={len(c_term)}, n_term={len(n_term)}")
        print("  Segments:", sorted(set(u.atoms.segids)))
        sys.exit(1)

    dists = [np.linalg.norm(c_term.positions[0] - n_term.positions[0])
             for _ in u.trajectory]

    arr = np.array(dists)
    pct = float((arr < 15.0).mean() * 100)
    metrics = {
        "n_frames": int(len(arr)),
        "time_ns": float(u.trajectory.n_frames * u.trajectory.dt / 1000),
        "mean_dist_A": float(arr.mean()),
        "min_dist_A":  float(arr.min()),
        "max_dist_A":  float(arr.max()),
        "std_dist_A":  float(arr.std()),
        "pct_below_15A": pct,
        "pass": bool(pct > 50.0),
    }
    (rep_dir / "metrics.json").write_text(json.dumps(metrics, indent=2))
    np.savetxt(rep_dir / "recon_dist.txt",
               np.column_stack([np.arange(len(arr)) * u.trajectory.dt / 1000, arr]),
               header="time_ns  dist_A", fmt="%.3f  %.3f")

    print(f"  Frames: {metrics['n_frames']}  ({metrics['time_ns']:.1f} ns)")
    print(f"  Inter-half: {metrics['mean_dist_A']:.1f} ± {metrics['std_dist_A']:.1f} Å")
    print(f"  Min: {metrics['min_dist_A']:.1f} Å   Max: {metrics['max_dist_A']:.1f} Å")
    print(f"  % < 15 Å: {pct:.1f}%  PASS: {metrics['pass']}")
    set_state(rep_dir, "ANALYSIS_DONE")


def print_summary():
    """Print combined results across all replicates including rep1."""
    rows = []
    for rep_id, rep_dir in [("rep1", REP1_DIR)] + [(r["id"], MD_DIR / r["id"]) for r in REPLICATES]:
        m_path = rep_dir / "metrics.json"
        if m_path.exists():
            m = json.loads(m_path.read_text())
            rows.append((rep_id, m["mean_dist_A"], m["std_dist_A"],
                         m["min_dist_A"], m["pct_below_15A"], m["pass"]))
    if rows:
        print(f"\n{'='*65}")
        print(f"  {'Rep':<6} {'Mean±SD':>12}  {'Min':>6}  {'%<15Å':>7}  Pass")
        print(f"  {'-'*55}")
        for rep_id, mean, sd, mn, pct, passed in rows:
            print(f"  {rep_id:<6} {mean:6.1f}±{sd:<5.1f}  {mn:6.1f}  {pct:7.1f}%  {'✓' if passed else '✗'}")
        all_pct = [r[4] for r in rows]
        print(f"  {'-'*55}")
        print(f"  {'Mean':<6} {sum(r[1] for r in rows)/len(rows):>12.1f}         {sum(all_pct)/len(all_pct):7.1f}%")
        print(f"{'='*65}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--status", action="store_true")
    args = parser.parse_args()

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
