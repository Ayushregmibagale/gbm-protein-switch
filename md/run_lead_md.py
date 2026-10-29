#!/usr/bin/env python3
"""
GROMACS MD for GBM protein switch lead construct.

Lead: BC mpnn17 × ph_sox9_run18, GGS×1, flipped orientation
Chain A: SOX2 target (97 aa)
Chain B: SOX9 target (72 aa)
Chain C: mpnn17(71aa) + GGS + CaspN(aa140-315, 176aa) = 250 aa
Chain D: CaspC(aa316-416, 101aa) + GGS + ph_sox9_run18(67aa) = 171 aa

Pass criterion: inter-half distance (C-term chain C ↔ N-term chain D) < 15 Å
in > 50% of production frames.

Pipeline (fully resumable):
  PENDING         → initial state
  STRUCTURE_DONE  → PDB copied + TER fixed
  TOPOLOGY_DONE   → pdb2gmx + box + solvate + ions
  EM_DONE         → energy minimisation
  NVT_DONE        → 200 ps NVT (restrained)
  NPT_DONE        → 500 ps NPT (restrained)
  PROD_DONE       → 20 ns production MD
  ANALYSIS_DONE   → inter-half distance → metrics.json

Usage:
  python run_lead_md.py          # run / resume
  python run_lead_md.py --status # print current stage
"""

import os, sys, json, subprocess, argparse
import numpy as np
from pathlib import Path

BASE    = Path("/mnt/d/gbm-protein-switch")
MD_DIR  = BASE / "md" / "lead_construct"
MDP_DIR = BASE / "md" / "mdp"
SRC_PDB = (BASE / "selection" / "stage5_flipped_outputs"
           / "flip_n17_model2_ph_sox9_run18_GGS1_s2"
           / "boltz_results_input" / "predictions" / "input" / "input_model_0.pdb")
GMX     = "/home/ayushregmibagale/.conda/envs/gromacs/bin/gmx"
STATE   = MD_DIR / "state.txt"

CHAIN_A_RES = 97   # SOX2
CHAIN_B_RES = 72   # SOX9
CHAIN_C_RES = 250  # mpnn17 + GGS + CaspN
CHAIN_D_RES = 171  # CaspC + GGS + ph_sox9_run18


def get_state() -> str:
    return STATE.read_text().strip() if STATE.exists() else "PENDING"


def set_state(s: str):
    STATE.write_text(s + "\n")
    print(f"  → state: {s}")


def gmx(args: str, log: Path, stdin: str = None):
    env = os.environ.copy()
    env["GMX_MAXBACKUP"] = "-1"
    cmd = f"{GMX} {args}"
    with open(log, "a") as f:
        r = subprocess.run(cmd, shell=True, cwd=MD_DIR, env=env,
                           input=stdin, text=True,
                           stdout=f, stderr=subprocess.STDOUT)
    if r.returncode != 0:
        print(f"  FAILED (rc={r.returncode}) — see {log.name}")
        sys.exit(1)


# ── Stage 1: Fix PDB ──────────────────────────────────────────────────────────
def stage_structure():
    """Copy PDB, add missing TER for chain D, renumber residues to be unique."""
    MD_DIR.mkdir(parents=True, exist_ok=True)
    print("  Fixing PDB: adding TER for chain D, making residues unique...")

    lines = SRC_PDB.read_text().splitlines()
    offsets = {"A": 0, "B": CHAIN_A_RES, "C": CHAIN_A_RES + CHAIN_B_RES,
               "D": CHAIN_A_RES + CHAIN_B_RES + CHAIN_C_RES}

    fixed = []
    for line in lines:
        if line.startswith("ATOM") or line.startswith("HETATM"):
            chain = line[21]
            if chain in offsets:
                orig_resnum = int(line[22:26])
                new_resnum  = orig_resnum + offsets[chain]
                line = line[:22] + f"{new_resnum:4d}" + line[26:]
        if line.startswith("END"):
            fixed.append("TER")   # add missing chain D TER
        fixed.append(line)

    (MD_DIR / "start.pdb").write_text("\n".join(fixed) + "\n")

    # Verify key atoms
    for ln in (MD_DIR / "start.pdb").read_text().splitlines():
        if " CA " in ln and " C 419" in ln:
            print(f"  CaspN C-term (C419 CA): {ln[30:54].strip()}")
        if " CA " in ln and " D 420" in ln:
            print(f"  CaspC N-term (D420 CA): {ln[30:54].strip()}")

    set_state("STRUCTURE_DONE")


# ── Stage 2: Topology ─────────────────────────────────────────────────────────
def stage_topology():
    print("  pdb2gmx (AMBER99SB-ILDN, TIP3P)...")
    gmx("pdb2gmx -f start.pdb -o protein.gro -p topol.top -i posre.itp "
        "-ff amber99sb-ildn -water tip3p -ignh -missing",
        log=MD_DIR / "pdb2gmx.log")

    print("  editconf (dodecahedron, 1.2 nm)...")
    gmx("editconf -f protein.gro -o boxed.gro -bt dodecahedron -d 1.2",
        log=MD_DIR / "editconf.log")

    print("  solvate...")
    gmx("solvate -cp boxed.gro -cs spc216.gro -o solvated.gro -p topol.top",
        log=MD_DIR / "solvate.log")

    print("  genion (0.15 M NaCl)...")
    gmx(f"grompp -f {MDP_DIR}/em.mdp -c solvated.gro -p topol.top "
        "-o ions.tpr -maxwarn 2",
        log=MD_DIR / "grompp_ions.log")

    # Find SOL group index
    r = subprocess.run(f"{GMX} make_ndx -f solvated.gro -o /dev/null",
                       shell=True, cwd=MD_DIR, input="q\n", text=True,
                       capture_output=True, env={**os.environ, "GMX_MAXBACKUP": "-1"})
    sol_group = "15"
    for line in (r.stdout + r.stderr).split("\n"):
        if "SOL" in line and line.strip() and line.strip()[0].isdigit():
            sol_group = line.strip().split()[0]
            break

    gmx(f"genion -s ions.tpr -o ions.gro -p topol.top "
        "-pname NA -nname CL -neutral -conc 0.15",
        log=MD_DIR / "genion.log", stdin=f"{sol_group}\n")

    set_state("TOPOLOGY_DONE")


# ── Stage 3: EM ───────────────────────────────────────────────────────────────
def stage_em():
    print("  Energy minimisation...")
    gmx(f"grompp -f {MDP_DIR}/em.mdp -c ions.gro -p topol.top -o em.tpr -maxwarn 2",
        log=MD_DIR / "grompp_em.log")
    gmx("mdrun -v -deffnm em -ntmpi 1 -ntomp 8 -gpu_id 0 -nb gpu",
        log=MD_DIR / "mdrun_em.log")
    set_state("EM_DONE")


# ── Stage 4: NVT ──────────────────────────────────────────────────────────────
def stage_nvt():
    print("  NVT equilibration (200 ps)...")
    gmx(f"grompp -f {MDP_DIR}/nvt.mdp -c em.gro -r em.gro -p topol.top "
        "-o nvt.tpr -maxwarn 2",
        log=MD_DIR / "grompp_nvt.log")
    gmx("mdrun -v -deffnm nvt -ntmpi 1 -ntomp 8 -gpu_id 0 -nb gpu -pme gpu -bonded gpu",
        log=MD_DIR / "mdrun_nvt.log")
    set_state("NVT_DONE")


# ── Stage 5: NPT ──────────────────────────────────────────────────────────────
def stage_npt():
    print("  NPT equilibration (500 ps)...")
    gmx(f"grompp -f {MDP_DIR}/npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt "
        "-p topol.top -o npt.tpr -maxwarn 2",
        log=MD_DIR / "grompp_npt.log")
    gmx("mdrun -v -deffnm npt -ntmpi 1 -ntomp 8 -gpu_id 0 -nb gpu -pme gpu -bonded gpu",
        log=MD_DIR / "mdrun_npt.log")
    set_state("NPT_DONE")


# ── Stage 6: Production ───────────────────────────────────────────────────────
def stage_prod():
    print("  Production MD (20 ns)...")
    gmx(f"grompp -f {MDP_DIR}/prod_20ns.mdp -c npt.gro -t npt.cpt "
        "-p topol.top -o prod.tpr -maxwarn 2",
        log=MD_DIR / "grompp_prod.log")
    gmx("mdrun -v -deffnm prod -ntmpi 1 -ntomp 8 -gpu_id 0 -nb gpu -pme gpu -bonded gpu",
        log=MD_DIR / "mdrun_prod.log")
    set_state("PROD_DONE")


# ── Stage 7: Analysis ─────────────────────────────────────────────────────────
def stage_analysis():
    """MDAnalysis: track inter-half distance across all production frames."""
    try:
        import MDAnalysis as mda
    except ImportError:
        print("  MDAnalysis not found — install with: pip install MDAnalysis")
        sys.exit(1)

    print("  Loading trajectory with MDAnalysis...")
    u = mda.Universe(str(MD_DIR / "prod.tpr"), str(MD_DIR / "prod.xtc"))

    # After pdb2gmx with renumbered start.pdb:
    #   Chain C (PROC): residues 170-419 (original per-chain 1-250)
    #   Chain D (PROD): residues 420-590 (original per-chain 1-171)
    # But pdb2gmx preserves the *input* residue numbers, which after renumbering are:
    #   Chain C last residue = 170+250-1 = 419
    #   Chain D first residue = 420
    #
    # MDAnalysis resnum matches PDB residue numbers preserved by pdb2gmx.
    # segid in GROMACS TPR = PROA/PROB/PROC/PROD in chain order.

    # Try segid selection first
    c_term = u.select_atoms("segid PROC and resnum 419 and name CA")
    n_term = u.select_atoms("segid PROD and resnum 420 and name CA")

    # Fallback: use bynum based on cumulative atom count (estimate after H addition)
    if len(c_term) == 0 or len(n_term) == 0:
        print("  segid selection failed, trying resnum only...")
        c_term = u.select_atoms("resnum 419 and name CA")
        n_term = u.select_atoms("resnum 420 and name CA")

    if len(c_term) == 0 or len(n_term) == 0:
        print(f"  ERROR: Could not find CA atoms. c_term={len(c_term)}, n_term={len(n_term)}")
        print("  Segments:", sorted(set(u.atoms.segids)))
        print("  Residues near junction:")
        for seg in u.segments:
            res_ids = [r.resnum for r in seg.residues]
            if 419 in res_ids or 420 in res_ids:
                print(f"    segid={seg.segid} resids={min(res_ids)}-{max(res_ids)}")
        sys.exit(1)

    print(f"  C-term (chain C res 419): {c_term}")
    print(f"  N-term (chain D res 420): {n_term}")

    dists = []
    for ts in u.trajectory:
        d = np.linalg.norm(c_term.positions[0] - n_term.positions[0])
        dists.append(d)

    arr = np.array(dists)
    pct_pass = float((arr < 15.0).mean() * 100)
    metrics = {
        "n_frames": int(len(arr)),
        "time_ns": float(u.trajectory.n_frames * u.trajectory.dt / 1000),
        "mean_dist_A": float(arr.mean()),
        "min_dist_A":  float(arr.min()),
        "max_dist_A":  float(arr.max()),
        "std_dist_A":  float(arr.std()),
        "pct_below_15A": pct_pass,
        "pass": bool(pct_pass > 50.0),
    }
    (MD_DIR / "metrics.json").write_text(json.dumps(metrics, indent=2))

    # Save raw distances for plotting
    np.savetxt(MD_DIR / "recon_dist.txt",
               np.column_stack([np.arange(len(arr)) * u.trajectory.dt / 1000, arr]),
               header="time_ns  dist_A", fmt="%.3f  %.3f")

    print(f"\n{'='*60}")
    print(f"  LEAD CONSTRUCT MD RESULTS")
    print(f"  Frames: {metrics['n_frames']}  ({metrics['time_ns']:.1f} ns)")
    print(f"  Inter-half dist: {metrics['mean_dist_A']:.1f} ± {metrics['std_dist_A']:.1f} Å")
    print(f"  Min: {metrics['min_dist_A']:.1f} Å   Max: {metrics['max_dist_A']:.1f} Å")
    print(f"  % frames < 15 Å: {pct_pass:.1f}%")
    print(f"  PASS (>50%): {metrics['pass']}")
    print(f"{'='*60}")

    set_state("ANALYSIS_DONE")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--status", action="store_true")
    args = parser.parse_args()

    state = get_state()
    if args.status:
        print(f"State: {state}")
        if (MD_DIR / "metrics.json").exists():
            print(json.dumps(json.loads((MD_DIR / "metrics.json").read_text()), indent=2))
        return

    print(f"GBM lead construct MD  (state: {state})")

    steps = [
        ("PENDING",        "STRUCTURE_DONE",  "[1/7] Fixing PDB...",          stage_structure),
        ("STRUCTURE_DONE", "TOPOLOGY_DONE",   "[2/7] Building topology...",    stage_topology),
        ("TOPOLOGY_DONE",  "EM_DONE",         "[3/7] Energy minimisation...",  stage_em),
        ("EM_DONE",        "NVT_DONE",        "[4/7] NVT equilibration...",    stage_nvt),
        ("NVT_DONE",       "NPT_DONE",        "[5/7] NPT equilibration...",    stage_npt),
        ("NPT_DONE",       "PROD_DONE",       "[6/7] Production MD (20 ns)...",stage_prod),
        ("PROD_DONE",      "ANALYSIS_DONE",   "[7/7] Analysis...",             stage_analysis),
    ]

    for from_state, to_state, label, fn in steps:
        state = get_state()
        if state == from_state:
            print(f"\n{label}")
            fn()

    state = get_state()
    if state == "ANALYSIS_DONE":
        print("\nAll stages complete.")


if __name__ == "__main__":
    main()
