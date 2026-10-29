#!/usr/bin/env python3
"""
Unbound-state MD for both leads — leakiness check.

Tests whether caspase halves spontaneously reconstitute when binders are NOT
bound to their SOX2/SOX9 targets.

Pass criterion (non-leaky): inter-half CA distance > 25 Å in ≥ 95% of frames.

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
  python run_unbound_md.py track1       # run Track 1 unbound
  python run_unbound_md.py track2       # run Track 2 unbound
  python run_unbound_md.py --status     # print both states
"""

import os, sys, json, subprocess, argparse
import numpy as np
from pathlib import Path

BASE    = Path("/mnt/d/gbm-protein-switch")
MD_DIR  = BASE / "md" / "unbound"
MDP_DIR = BASE / "md" / "mdp"
BOLTZ   = Path("/home/ayushregmibagale/miniconda3/envs/boltz2/bin/boltz")
GMX     = "/home/ayushregmibagale/.conda/envs/gromacs/bin/gmx"
PYTHON  = "/home/ayushregmibagale/.conda/envs/gromacs/bin/python3"

# ── Sequences ─────────────────────────────────────────────────────────────────

MPNN17      = "PYEIPEEDKEWDEYSFAMVEVHERLMEAERQLRLLSPEEEEKKRELEEEIERARKELEEVHKRGEAIERKM"
PH_SOX9_R18 = "MEEVRELLERARELEEKVRELEDDPSKREERERLLEEAERLREEALDLLVEQIREPLLQWMREEGLL"
BRIDGE_48   = "GAAAAAAAEAAAAAAAAAAAAAAAAAAAAEAARAAAAAAAAAAAAAAAAAAQAAAAQAYLAQAQAAAAA"

_casp = "".join(open(BASE/"selection"/"stage5_caspase_outputs"/"P55211.fasta").read().strip().split("\n")[1:])
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

def get_state(d: Path) -> str:
    f = d / "state.txt"
    return f.read_text().strip() if f.exists() else "PENDING"

def set_state(d: Path, s: str):
    (d / "state.txt").write_text(s)

def gmx(cmd: str, log: Path, cwd: Path):
    full = f"{GMX} {cmd}"
    with open(log, "w") as fh:
        r = subprocess.run(full, shell=True, cwd=str(cwd),
                           stdout=fh, stderr=subprocess.STDOUT)
    if r.returncode != 0:
        print(f"  ERROR running: {full}")
        print(open(log).read()[-500:])
        sys.exit(1)


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
        conf = list(run_dir.rglob("confidence_*.json"))
        if not conf:
            print(f"    seed {seed}: FAILED")
            continue
        with open(conf[0]) as f: c = json.load(f)
        iptm = c.get("iptm", c.get("protein_iptm", 0))
        pdb  = list(run_dir.rglob("*.pdb"))
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
    gmx(f"grompp -f {MDP_DIR}/em.mdp -c solvated.gro -p topol.top "
        f"-o ions.tpr -maxwarn 2",
        md_dir / "grompp_ions.log", md_dir)
    r = subprocess.run(f"echo 'SOL' | {GMX} genion -s ions.tpr -o ions.gro "
                       f"-p topol.top -pname NA -nname CL -neutral",
                       shell=True, cwd=str(md_dir),
                       capture_output=True, text=True)
    (md_dir / "genion.log").write_text(r.stdout + r.stderr)
    set_state(md_dir, "TOPOLOGY_DONE")


def stage_em(md_dir: Path):
    print("  Energy minimisation...")
    gmx(f"grompp -f {MDP_DIR}/em.mdp -c ions.gro -p topol.top -o em.tpr -maxwarn 2",
        md_dir / "grompp_em.log", md_dir)
    gmx("mdrun -v -deffnm em -ntmpi 1 -ntomp 8 -gpu_id 0 -nb gpu",
        md_dir / "mdrun_em.log", md_dir)
    set_state(md_dir, "EM_DONE")


def stage_nvt(md_dir: Path):
    print("  NVT equilibration (200 ps)...")
    gmx(f"grompp -f {MDP_DIR}/nvt.mdp -c em.gro -r em.gro -p topol.top "
        f"-o nvt.tpr -maxwarn 2",
        md_dir / "grompp_nvt.log", md_dir)
    gmx("mdrun -v -deffnm nvt -ntmpi 1 -ntomp 8 -gpu_id 0 -nb gpu -pme gpu -bonded gpu",
        md_dir / "mdrun_nvt.log", md_dir)
    set_state(md_dir, "NVT_DONE")


def stage_npt(md_dir: Path):
    print("  NPT equilibration (500 ps)...")
    gmx(f"grompp -f {MDP_DIR}/npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt "
        f"-p topol.top -o npt.tpr -maxwarn 2",
        md_dir / "grompp_npt.log", md_dir)
    gmx("mdrun -v -deffnm npt -ntmpi 1 -ntomp 8 -gpu_id 0 -nb gpu -pme gpu -bonded gpu",
        md_dir / "mdrun_npt.log", md_dir)
    set_state(md_dir, "NPT_DONE")


def stage_prod(md_dir: Path):
    print("  Production MD (20 ns)...")
    gmx(f"grompp -f {MDP_DIR}/prod_20ns.mdp -c npt.gro -t npt.cpt "
        f"-p topol.top -o prod.tpr -maxwarn 2",
        md_dir / "grompp_prod.log", md_dir)
    gmx("mdrun -v -deffnm prod -ntmpi 1 -ntomp 8 -gpu_id 0 -nb gpu -pme gpu -bonded gpu",
        md_dir / "mdrun_prod.log", md_dir)
    set_state(md_dir, "PROD_DONE")


def stage_analysis(md_dir: Path, cterm_res: int, nterm_res: int, label: str):
    print("  Analysing inter-half distance...")
    try:
        import MDAnalysis as mda
    except ImportError:
        print("  MDAnalysis not found"); sys.exit(1)

    u = mda.Universe(str(md_dir/"prod.tpr"), str(md_dir/"prod.xtc"))
    c_sel = u.select_atoms(f"resid {cterm_res} and name CA")
    n_sel = u.select_atoms(f"resid {nterm_res} and name CA")

    if len(c_sel) != 1 or len(n_sel) != 1:
        print(f"  ERROR: cterm={len(c_sel)}, nterm={len(n_sel)} atoms found")
        print(f"  Available resids: {sorted(set(u.atoms.resids))[:20]} ...")
        sys.exit(1)

    dists = []
    for _ in u.trajectory:
        dists.append(float(np.linalg.norm(c_sel.positions[0] - n_sel.positions[0])))

    arr = np.array(dists)
    pct_gt25 = float(np.mean(arr > 25.0) * 100)
    pct_gt15 = float(np.mean(arr > 15.0) * 100)
    passed   = pct_gt25 >= 95.0

    metrics = {
        "label":        label,
        "n_frames":     int(len(arr)),
        "mean_dist_A":  round(float(arr.mean()), 2),
        "std_dist_A":   round(float(arr.std()),  2),
        "min_dist_A":   round(float(arr.min()),  2),
        "max_dist_A":   round(float(arr.max()),  2),
        "pct_gt_25A":   round(pct_gt25, 1),
        "pct_gt_15A":   round(pct_gt15, 1),
        "pass_nonleaky": passed,
    }
    (md_dir / "metrics.json").write_text(json.dumps(metrics, indent=2))
    np.savetxt(md_dir / "interhalfdist.txt",
               np.column_stack([np.arange(len(arr)) * u.trajectory.dt / 1000, arr]),
               header="time_ns  dist_A", fmt="%.3f  %.3f")

    print(f"\n  {'='*50}")
    print(f"  {label} — UNBOUND STATE")
    print(f"  Inter-half: {metrics['mean_dist_A']:.1f} ± {metrics['std_dist_A']:.1f} Å")
    print(f"  Min: {metrics['min_dist_A']:.1f} Å  Max: {metrics['max_dist_A']:.1f} Å")
    print(f"  % frames > 25 Å: {pct_gt25:.1f}%  (pass ≥ 95%)")
    print(f"  % frames > 15 Å: {pct_gt15:.1f}%")
    verdict = "NON-LEAKY ✓" if passed else "LEAKY — REDESIGN LINKER"
    print(f"  Verdict: {verdict}")
    print(f"  {'='*50}\n")

    set_state(md_dir, "ANALYSIS_DONE")
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
            print(f"  Non-leaky: {m['pass_nonleaky']} | "
                  f"mean={m['mean_dist_A']:.1f}Å | "
                  f"%>25Å={m['pct_gt_25A']:.1f}%")
            break
        if state not in stages:
            print(f"Unknown state: {state}"); break
        fn, fn_args, _ = stages[state]
        fn(*fn_args)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("track", nargs="?", choices=["track1","track2","both"],
                        default="both", help="Which track to run")
    parser.add_argument("--status", action="store_true")
    args = parser.parse_args()

    MD_DIR.mkdir(parents=True, exist_ok=True)

    if args.status:
        for key, cfg in TRACKS.items():
            d = MD_DIR / cfg["label"]
            state = get_state(d)
            print(f"{cfg['label']}: {state}")
            if state == "ANALYSIS_DONE" and (d/"metrics.json").exists():
                m = json.loads((d/"metrics.json").read_text())
                print(f"  mean={m['mean_dist_A']:.1f}Å  %>25Å={m['pct_gt_25A']:.1f}%  pass={m['pass_nonleaky']}")
        return

    keys = ["track1","track2"] if args.track == "both" else [args.track]
    for key in keys:
        run_track(key)


if __name__ == "__main__":
    main()
