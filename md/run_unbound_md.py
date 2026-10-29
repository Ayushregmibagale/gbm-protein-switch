#!/usr/bin/env python3
"""
Unbound-state MD for both leads — leakiness check.

Tests whether caspase halves spontaneously reconstitute when binders are NOT
bound to their SOX2/SOX9 targets.

Descriptive distance criterion: inter-half CA distance >25 Å in ≥95% of frames.
This does not establish biological leakiness or safety.

Track 1 (two-binder split-caspase):
  Chain A: mpnn17(71aa) + GGS + CaspN(176aa)  = 250 aa
  Chain B: CaspC(101aa) + GGS + ph_sox9_run18(67aa) = 171 aa
  → 2-chain simulation; measure last-CA(A) ↔ first-CA(B)

Track 2 (bridge binder):
  Single chain: CaspN(176aa) + GGS + bridge_48_s01(69aa) + GGS + CaspC(101aa) = 352 aa
  → 1-chain simulation; measure CA(res176) ↔ CA(res252)

Pipeline (resumable):
  PENDING → STRUCTURE_DONE → TOPOLOGY_DONE → EM_DONE → NVT_DONE → NPT_DONE
  → PROD_DONE → ANALYSIS_DONE

Usage:
  python run_unbound_md.py track1 --execute  # explicitly run/resume
  python run_unbound_md.py track2 --execute  # explicitly run/resume
  python run_unbound_md.py --status     # print both states
"""

import os, sys, json, subprocess, argparse, shlex
import numpy as np
from pathlib import Path
from workflow_support import read_state, write_state, run_gmx, production, periodic_analysis

BASE    = Path(os.environ.get("GBM_ROOT", Path(__file__).resolve().parents[1]))
MD_DIR  = BASE / "md" / "unbound"
MDP_DIR = BASE / "md" / "mdp"
BOLTZ   = Path(os.environ.get("BOLTZ_BIN", "boltz"))
GMX     = os.environ.get("GMX_BIN", "gmx")
PYTHON  = "/home/ayushregmibagale/.conda/envs/gromacs/bin/python3"

# ── Sequences ─────────────────────────────────────────────────────────────────

MPNN17      = "PYEIPEEDKEWDEYSFAMVEVHERLMEAERQLRLLSPEEEEKKRELEEEIERARKELEEVHKRGEAIERKM"
PH_SOX9_R18 = "MEEVRELLERARELEEKVRELEDDPSKREERERLLEEAERLREEALDLLVEQIREPLLQWMREEGLL"
BRIDGE_48   = "GAAAAAAAEAAAAAAAAAAAAAAAAAAAAEAARAAAAAAAAAAAAAAAAAAQAAAAQAYLAQAQAAAAA"

_casp = "".join((Path(__file__).resolve().parents[1] / "inputs/reference/P55211.fasta").read_text().splitlines()[1:])
if len(_casp) != 416 or _casp[314:316] != "DA":
    raise ValueError("Unexpected reference sequence length or split identity")
CASP_N = _casp[139:315]   # aa140-315, 176 aa
CASP_C = _casp[315:416]   # aa316-416, 101 aa

TRACK1_A = MPNN17 + "GGS" + CASP_N          # 250 aa: mpnn17-GGS-CaspN
TRACK1_B = CASP_C + "GGS" + PH_SOX9_R18    # 171 aa: CaspC-GGS-ph_sox9_run18
TRACK2   = CASP_N + "GGS" + BRIDGE_48 + "GGS" + CASP_C  # 352 aa single chain

# Residue numbers for analysis (1-indexed, post-renumbering in TPR)
# Track 1: chain A ends at res 250; chain B starts at res 251
T1_CTERM = 250   # last residue of chain A (CaspN C-term)
T1_NTERM = 251   # first residue of chain B (CaspC N-term)
# Track 2: CaspN ends at res 176; CaspC starts at res 252 (176+3+69+3+1)
T2_CTERM = 176   # CaspN C-term
T2_NTERM = 252   # CaspC N-term

TRACKS = {
    "track1": {
        "label":   "Track1_unbound",
        "chains":  [("A", TRACK1_A), ("B", TRACK1_B)],
        "cterm":   T1_CTERM,
        "nterm":   T1_NTERM,
        "nchains": 2,
    },
    "track2": {
        "label":   "Track2_unbound",
        "chains":  [("A", TRACK2)],
        "cterm":   T2_CTERM,
        "nterm":   T2_NTERM,
        "nchains": 1,
    },
}


# ── Utilities ─────────────────────────────────────────────────────────────────

def get_state(d):
    return read_state(d)


def set_state(d, s):
    write_state(d, s)


def gmx(cmd, log, cwd):
    run_gmx(GMX, cmd, log, cwd)



# ── Stage: predict structure with Boltz-2 ─────────────────────────────────────

def stage_structure(md_dir: Path, cfg: dict):
    print("  Predicting unbound structure with Boltz-2...")
    pred_dir = md_dir / "boltz_pred"
    pred_dir.mkdir(exist_ok=True)

    yaml_lines = ["sequences:"]
    for chain_id, seq in cfg["chains"]:
        yaml_lines += [
            f"  - protein:",
            f"      id: {chain_id}",
            f"      sequence: {seq}",
            f"      msa: empty",
        ]
    (pred_dir / "input.yaml").write_text("\n".join(yaml_lines) + "\n")

    best_iptm, best_pdb = -1.0, None
    for seed in [0, 1, 2]:
        run_dir = pred_dir / f"seed{seed}"
        run_dir.mkdir(exist_ok=True)
        env = os.environ.copy()
        env["LD_LIBRARY_PATH"] = "/usr/lib/wsl/lib:" + env.get("LD_LIBRARY_PATH","")
        env["XLA_PYTHON_CLIENT_PREALLOCATE"] = "false"
        r = subprocess.run(
            [str(BOLTZ), "predict", str(pred_dir/"input.yaml"),
             "--out_dir", str(run_dir), "--recycling_steps", "3",
             "--sampling_steps", "200", "--output_format", "pdb",
             "--no_kernels", "--seed", str(seed)],
            capture_output=True, text=True, env=env)
        if r.returncode != 0:
            print(f"    seed {seed}: inference failed; ignoring any stale outputs")
            continue
        conf = sorted(run_dir.rglob("confidence_*model_0.json"))
        if len(conf) != 1:
            print(f"    seed {seed}: FAILED")
            continue
        with open(conf[0]) as f: c = json.load(f)
        iptm = c.get("iptm", c.get("protein_iptm"))
        if not isinstance(iptm, (int, float)) or not np.isfinite(iptm) or not 0 <= iptm <= 1:
            raise ValueError(f"Invalid ipTM for seed {seed}")
        pdb  = sorted(conf[0].parent.glob("*_model_0.pdb"))
        if len(pdb) != 1:
            raise ValueError(f"Expected exactly one model_0 structure for seed {seed}")
        print(f"    seed {seed}: ipTM={iptm:.3f}")
        if iptm > best_iptm and pdb:
            best_iptm, best_pdb = iptm, pdb[0]

    if best_pdb is None:
        print("  ERROR: all Boltz-2 seeds failed"); sys.exit(1)

    import shutil
    shutil.copy(best_pdb, md_dir / "start.pdb")
    print(f"  Best structure: seed ipTM={best_iptm:.3f} → start.pdb")
    set_state(md_dir, "STRUCTURE_DONE")


def stage_topology(md_dir: Path):
    print("  Building topology...")
    gmx(f"pdb2gmx -f start.pdb -o protein.gro -p topol.top "
        f"-ff amber99sb-ildn -water tip3p -ignh",
        md_dir / "pdb2gmx.log", md_dir)
    gmx(f"editconf -f protein.gro -o boxed.gro -c -d 1.2 -bt dodecahedron",
        md_dir / "editconf.log", md_dir)
    gmx(f"solvate -cp boxed.gro -cs spc216.gro -o solvated.gro -p topol.top",
        md_dir / "solvate.log", md_dir)
    gmx(f"grompp -f {shlex.quote(str(MDP_DIR / 'em.mdp'))} -c solvated.gro -p topol.top "
        f"-o ions.tpr ",
        md_dir / "grompp_ions.log", md_dir)
    run_gmx(GMX, ['genion', '-s', 'ions.tpr', '-o', 'ions.gro', '-p', 'topol.top', '-pname', 'NA', '-nname', 'CL', '-neutral'], md_dir / 'genion.log', md_dir, 'SOL\n')
    set_state(md_dir, "TOPOLOGY_DONE")


def stage_em(md_dir: Path):
    print("  Energy minimisation...")
    gmx(f"grompp -f {shlex.quote(str(MDP_DIR / 'em.mdp'))} -c ions.gro -p topol.top -o em.tpr ",
        md_dir / "grompp_em.log", md_dir)
    gmx("mdrun -v -deffnm em -ntmpi 1 -ntomp 8 -gpu_id 0 -nb gpu",
        md_dir / "mdrun_em.log", md_dir)
    set_state(md_dir, "EM_DONE")


def stage_nvt(md_dir: Path):
    print("  NVT equilibration (200 ps)...")
    gmx(f"grompp -f {shlex.quote(str(MDP_DIR / 'nvt.mdp'))} -c em.gro -r em.gro -p topol.top "
        f"-o nvt.tpr ",
        md_dir / "grompp_nvt.log", md_dir)
    gmx("mdrun -v -deffnm nvt -ntmpi 1 -ntomp 8 -gpu_id 0 -nb gpu -pme gpu -bonded gpu",
        md_dir / "mdrun_nvt.log", md_dir)
    set_state(md_dir, "NVT_DONE")


def stage_npt(md_dir: Path):
    print("  NPT equilibration (500 ps)...")
    gmx(f"grompp -f {shlex.quote(str(MDP_DIR / 'npt.mdp'))} -c nvt.gro -r nvt.gro -t nvt.cpt "
        f"-p topol.top -o npt.tpr ",
        md_dir / "grompp_npt.log", md_dir)
    gmx("mdrun -v -deffnm npt -ntmpi 1 -ntomp 8 -gpu_id 0 -nb gpu -pme gpu -bonded gpu",
        md_dir / "mdrun_npt.log", md_dir)
    set_state(md_dir, "NPT_DONE")


def stage_prod(md_dir):
    prepare = ['grompp', '-f', str(MDP_DIR / 'prod_20ns.mdp'), '-c', 'npt.gro', '-t', 'npt.cpt', '-p', 'topol.top', '-o', 'prod.tpr']
    production(GMX, md_dir, prepare, 'production_checked')
    set_state(md_dir, 'PROD_DONE')



def stage_analysis(md_dir, cterm_res, nterm_res, label):
    metrics = periodic_analysis(md_dir, f'resid {cterm_res} and name CA', f'resid {nterm_res} and name CA', label, 'far')
    temporary = md_dir / 'metrics.json.tmp'
    temporary.write_text(json.dumps(metrics, indent=2) + '\n')
    temporary.replace(md_dir / 'metrics.json')
    set_state(md_dir, 'ANALYSIS_DONE')
    return metrics



# ── Main ──────────────────────────────────────────────────────────────────────

def run_track(track_key: str):
    cfg     = TRACKS[track_key]
    md_dir  = MD_DIR / cfg["label"]
    md_dir.mkdir(parents=True, exist_ok=True)

    stages = {
        "PENDING":        (stage_structure,  (md_dir, cfg),                         "STRUCTURE_DONE"),
        "STRUCTURE_DONE": (stage_topology,   (md_dir,),                             "TOPOLOGY_DONE"),
        "TOPOLOGY_DONE":  (stage_em,         (md_dir,),                             "EM_DONE"),
        "EM_DONE":        (stage_nvt,        (md_dir,),                             "NVT_DONE"),
        "NVT_DONE":       (stage_npt,        (md_dir,),                             "NPT_DONE"),
        "NPT_DONE":       (stage_prod,       (md_dir,),                             "PROD_DONE"),
        "PROD_DONE":      (stage_analysis,   (md_dir, cfg["cterm"], cfg["nterm"],
                                              cfg["label"]),                         "ANALYSIS_DONE"),
    }

    print(f"\n{'='*60}")
    print(f"  {cfg['label']}")
    print(f"  Chains: {[c for c,_ in cfg['chains']]}")
    for cid, seq in cfg["chains"]:
        print(f"    Chain {cid}: {len(seq)} aa")
    print(f"  Split site: res {cfg['cterm']} ↔ res {cfg['nterm']}")
    print(f"{'='*60}")

    while True:
        state = get_state(md_dir)
        print(f"\n[{cfg['label']}] State: {state}")
        if state == "ANALYSIS_DONE":
            print("Done.")
            m = json.loads((md_dir / "metrics.json").read_text())
            print(f"  Distance criterion: {m.get('distance_criterion_met', 'LEGACY')} | "
                  f"mean={m['mean_dist_A']:.1f}Å | "
                  f"%>25Å={m.get('pct_above_25A', float('nan')):.1f}%")
            break
        if state not in stages:
            raise ValueError(f"Unknown state: {state}")
        fn, fn_args, _ = stages[state]
        fn(*fn_args)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("track", nargs="?", choices=["track1","track2","both"],
                        default="both", help="Which track to run")
    parser.add_argument("--status", action="store_true")
    parser.add_argument('--execute', action='store_true', help='Run or resume stages; otherwise only report status')
    args = parser.parse_args()
    if not args.execute:
        args.status = True

    if args.execute:
        MD_DIR.mkdir(parents=True, exist_ok=True)

    if args.status:
        for key, cfg in TRACKS.items():
            d = MD_DIR / cfg["label"]
            state = get_state(d)
            print(f"{cfg['label']}: {state}")
            if state == "ANALYSIS_DONE" and (d/"metrics.json").exists():
                m = json.loads((d/"metrics.json").read_text())
                print(f"  mean={m['mean_dist_A']:.1f}Å  %>25Å={m.get('pct_above_25A', float('nan')):.1f}%  pass={m.get('distance_criterion_met', 'LEGACY')}")
        return

    keys = ["track1","track2"] if args.track == "both" else [args.track]
    for key in keys:
        run_track(key)


if __name__ == "__main__":
    main()
