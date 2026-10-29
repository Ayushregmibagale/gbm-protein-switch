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

import os, sys, json, subprocess, argparse, shlex
import numpy as np
from pathlib import Path
from workflow_support import read_state, write_state, run_gmx, production, periodic_analysis

BASE    = Path(os.environ.get("GBM_ROOT", Path(__file__).resolve().parents[1]))
MD_DIR  = BASE / "md" / "lead_construct"
MDP_DIR = BASE / "md" / "mdp"
SRC_PDB = (BASE / "selection" / "stage5_flipped_outputs"
           / "flip_n17_model2_ph_sox9_run18_GGS1_s2"
           / "boltz_results_input" / "predictions" / "input" / "input_model_0.pdb")
GMX     = os.environ.get("GMX_BIN", "gmx")
STATE   = MD_DIR / "state.txt"

CHAIN_A_RES = 97   # SOX2
CHAIN_B_RES = 72   # SOX9
CHAIN_C_RES = 250  # mpnn17 + GGS + CaspN
CHAIN_D_RES = 171  # CaspC + GGS + ph_sox9_run18


def get_state():
    return read_state(MD_DIR)



def set_state(s):
    write_state(MD_DIR, s)



def gmx(args, log, stdin=None):
    run_gmx(GMX, args, log, MD_DIR, stdin)



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
        "-ff amber99sb-ildn -water tip3p -ignh",
        log=MD_DIR / "pdb2gmx.log")

    print("  editconf (dodecahedron, 1.2 nm)...")
    gmx("editconf -f protein.gro -o boxed.gro -bt dodecahedron -d 1.2",
        log=MD_DIR / "editconf.log")

    print("  solvate...")
    gmx("solvate -cp boxed.gro -cs spc216.gro -o solvated.gro -p topol.top",
        log=MD_DIR / "solvate.log")

    print("  genion (0.15 M NaCl)...")
    gmx(f"grompp -f {shlex.quote(str(MDP_DIR / 'em.mdp'))} -c solvated.gro -p topol.top "
        "-o ions.tpr ",
        log=MD_DIR / "grompp_ions.log")

    gmx("genion -s ions.tpr -o ions.gro -p topol.top -pname NA -nname CL -neutral -conc 0.15",
        log=MD_DIR / "genion.log", stdin="SOL\n")

    set_state("TOPOLOGY_DONE")


# ── Stage 3: EM ───────────────────────────────────────────────────────────────
def stage_em():
    print("  Energy minimisation...")
    gmx(f"grompp -f {shlex.quote(str(MDP_DIR / 'em.mdp'))} -c ions.gro -p topol.top -o em.tpr ",
        log=MD_DIR / "grompp_em.log")
    gmx("mdrun -v -deffnm em -ntmpi 1 -ntomp 8 -gpu_id 0 -nb gpu",
        log=MD_DIR / "mdrun_em.log")
    set_state("EM_DONE")


# ── Stage 4: NVT ──────────────────────────────────────────────────────────────
def stage_nvt():
    print("  NVT equilibration (200 ps)...")
    gmx(f"grompp -f {shlex.quote(str(MDP_DIR / 'nvt.mdp'))} -c em.gro -r em.gro -p topol.top "
        "-o nvt.tpr ",
        log=MD_DIR / "grompp_nvt.log")
    gmx("mdrun -v -deffnm nvt -ntmpi 1 -ntomp 8 -gpu_id 0 -nb gpu -pme gpu -bonded gpu",
        log=MD_DIR / "mdrun_nvt.log")
    set_state("NVT_DONE")


# ── Stage 5: NPT ──────────────────────────────────────────────────────────────
def stage_npt():
    print("  NPT equilibration (500 ps)...")
    gmx(f"grompp -f {shlex.quote(str(MDP_DIR / 'npt.mdp'))} -c nvt.gro -r nvt.gro -t nvt.cpt "
        "-p topol.top -o npt.tpr ",
        log=MD_DIR / "grompp_npt.log")
    gmx("mdrun -v -deffnm npt -ntmpi 1 -ntomp 8 -gpu_id 0 -nb gpu -pme gpu -bonded gpu",
        log=MD_DIR / "mdrun_npt.log")
    set_state("NPT_DONE")


# ── Stage 6: Production ───────────────────────────────────────────────────────
def stage_prod():
    prepare = ['grompp', '-f', str(MDP_DIR / 'prod_20ns.mdp'), '-c', 'npt.gro', '-t', 'npt.cpt', '-p', 'topol.top', '-o', 'prod.tpr']
    production(GMX, MD_DIR, prepare, 'production_checked')
    set_state('PROD_DONE')



# ── Stage 7: Analysis ─────────────────────────────────────────────────────────
def stage_analysis():
    metrics = periodic_analysis(MD_DIR, 'resnum 419 and name CA', 'resnum 420 and name CA', 'lead_bound', 'near')
    temporary = MD_DIR / 'metrics.json.tmp'
    temporary.write_text(json.dumps(metrics, indent=2) + '\n')
    temporary.replace(MD_DIR / 'metrics.json')
    set_state('ANALYSIS_DONE')



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--status", action="store_true")
    parser.add_argument('--execute', action='store_true', help='Run or resume stages; otherwise only report status')
    args = parser.parse_args()
    if not args.execute:
        args.status = True

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
