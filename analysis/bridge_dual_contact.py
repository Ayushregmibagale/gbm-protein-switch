#!/usr/bin/env python3
"""
Dual-contact AND-gate analysis for bridge binder MD trajectories.

For each passing bridge design, measures per-frame simultaneous contact
with both SOX2 (chain A) and SOX9 (chain B).

Metric: bridge Cα within 8 Å of SOX2 AND bridge Cα within 8 Å of SOX9
in the same MD frame.

Contact thresholds tested: 0.05, 0.10, 0.15 (fraction of bridge Cα contacting target).

Output: analysis/dual_contact_results.tsv  +  per-design time series .txt files
"""

import sys
import argparse
import numpy as np
import json
from pathlib import Path

BASE    = Path(__file__).resolve().parents[1]
MD_ROOT = BASE / "md" / "bridge_binder_md"
OUT_DIR = BASE / "analysis"

CHAIN_A_RES = 97   # SOX2
CHAIN_B_RES = 72   # SOX9
CONTACT_CUTOFF = 8.0  # Å

THRESHOLDS = [0.05, 0.10, 0.15]

# Passing designs from bridge MD campaign
PASSING_DESIGNS = [
    "bridge_48_s01",
    "bridge_73_s09",
    "bridge_136_s01",
    "bridge_101_s06",
    "bridge_90_s01",
]


def infer_bridge_len(md_dir: Path) -> int:
    """Count residues in chain C of start.pdb."""
    residues = set()
    for line in (md_dir / "start.pdb").read_text().splitlines():
        if line.startswith("ATOM  ") and line[21:22] == "C" and line[12:16].strip() == "CA":
            residues.add((line[22:26], line[26:27]))
    if not residues:
        raise ValueError(f"{md_dir}: chain C has no CA residues")
    return len(residues)


def analyse_design(design: str) -> dict | None:
    md_dir = MD_ROOT / design
    tpr = md_dir / "prod.tpr"
    xtc = md_dir / "prod.xtc"
    if not (tpr.exists() and xtc.exists()):
        raise FileNotFoundError(f"{design}: missing tpr/xtc; no complete report issued")

    bridge_len = infer_bridge_len(md_dir)
    c_start = CHAIN_A_RES + CHAIN_B_RES + 1
    c_end   = CHAIN_A_RES + CHAIN_B_RES + bridge_len

    print(f"  {design}: bridge_len={bridge_len}, chain C resid {c_start}-{c_end}")

    try:
        import MDAnalysis as mda
    except ImportError:
        print("MDAnalysis not found"); sys.exit(1)

    u = mda.Universe(str(tpr), str(xtc))

    sox2_sel   = u.select_atoms(f"resid 1:{CHAIN_A_RES} and name CA")
    sox9_sel   = u.select_atoms(f"resid {CHAIN_A_RES+1}:{CHAIN_A_RES+CHAIN_B_RES} and name CA")
    bridge_sel = u.select_atoms(f"resid {c_start}:{c_end} and name CA")

    print(f"    SOX2={len(sox2_sel)} CA, SOX9={len(sox9_sel)} CA, bridge={len(bridge_sel)} CA")

    if (len(sox2_sel), len(sox9_sel), len(bridge_sel)) != (CHAIN_A_RES, CHAIN_B_RES, bridge_len):
        raise ValueError(f"{design}: CA selections do not match expected residue counts")
    from MDAnalysis.lib.distances import distance_array
    s2_fracs, s9_fracs, times = [], [], []

    for ts in u.trajectory:
        b_pos  = bridge_sel.positions
        s2_pos = sox2_sel.positions
        s9_pos = sox9_sel.positions

        box = ts.dimensions
        if box is None or not np.isfinite(box).all() or np.any(box[:3] <= 0) or np.any(box[3:] <= 0) or np.any(box[3:] >= 180):
            raise ValueError(f"{design}: missing or invalid periodic box")
        if not all(np.isfinite(x).all() for x in (b_pos, s2_pos, s9_pos)):
            raise ValueError(f"{design}: nonfinite coordinates")
        time = float(ts.time / 1000)
        if times and time <= times[-1]:
            raise ValueError(f"{design}: non-increasing trajectory times")
        times.append(time)
        s2_f = float(np.mean(distance_array(b_pos, s2_pos, box=box).min(axis=1) < CONTACT_CUTOFF))
        s9_f = float(np.mean(distance_array(b_pos, s9_pos, box=box).min(axis=1) < CONTACT_CUTOFF))
        s2_fracs.append(s2_f)
        s9_fracs.append(s9_f)

    s2 = np.array(s2_fracs)
    s9 = np.array(s9_fracs)
    n  = len(s2)

    if not n:
        raise ValueError(f"{design}: empty trajectory")

    # Save per-frame time series
    ts_file = OUT_DIR / f"{design}_dual_contact_timeseries.txt"
    np.savetxt(ts_file,
               np.column_stack([times, s2, s9]),
               header="time_ns  sox2_frac  sox9_frac", fmt="%.4f  %.4f  %.4f")

    result = {
        "design": design,
        "n_frames": n,
        "bridge_len": bridge_len,
        "sox2_mean": round(float(s2.mean()), 3),
        "sox9_mean": round(float(s9.mean()), 3),
    }

    for thr in THRESHOLDS:
        dual   = float(np.mean((s2 > thr) & (s9 > thr)) * 100)
        sox2_only = float(np.mean((s2 > thr) & (s9 <= thr)) * 100)
        sox9_only = float(np.mean((s2 <= thr) & (s9 > thr)) * 100)
        neither   = float(np.mean((s2 <= thr) & (s9 <= thr)) * 100)
        key = f"thr{int(thr*100):02d}"
        result[f"{key}_dual_pct"]      = round(dual, 1)
        result[f"{key}_sox2only_pct"]  = round(sox2_only, 1)
        result[f"{key}_sox9only_pct"]  = round(sox9_only, 1)
        result[f"{key}_neither_pct"]   = round(neither, 1)

    return result


def main():
    global MD_ROOT, OUT_DIR
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--md-root', type=Path, default=MD_ROOT)
    parser.add_argument('--out-dir', type=Path, required=True,
                        help='Use a new directory to preserve historical results')
    parser.add_argument('--designs', nargs='+', default=PASSING_DESIGNS)
    args = parser.parse_args()
    MD_ROOT, OUT_DIR = args.md_root, args.out_dir
    if len(set(args.designs)) != len(args.designs) or any(Path(x).name != x for x in args.designs):
        parser.error('Design names must be unique directory basenames')
    try:
        import MDAnalysis  # noqa
    except ImportError:
        print("ERROR: MDAnalysis not installed"); sys.exit(1)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    results = []

    for design in args.designs:
        print(f"\n[{design}]")
        r = analyse_design(design)
        if r:
            results.append(r)
            print(f"  SOX2 mean={r['sox2_mean']:.3f}  SOX9 mean={r['sox9_mean']:.3f}")
            for thr in THRESHOLDS:
                key = f"thr{int(thr*100):02d}"
                print(f"  thr={thr:.2f}: dual={r[f'{key}_dual_pct']:.1f}%  "
                      f"sox2only={r[f'{key}_sox2only_pct']:.1f}%  "
                      f"sox9only={r[f'{key}_sox9only_pct']:.1f}%  "
                      f"neither={r[f'{key}_neither_pct']:.1f}%")

    if not results:
        print("No results."); return

    # Write TSV
    tsv_path = OUT_DIR / "dual_contact_results.tsv"
    import csv
    fieldnames = list(results[0].keys())
    with open(tsv_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        w.writerows(results)
    print(f"\nResults saved to {tsv_path}")

    # Summary table
    print(f"\n{'='*75}")
    print(f"  {'Design':<20} {'SOX2':>6} {'SOX9':>6}  "
          f"{'Dual(0.10)':>10} {'Dual(0.05)':>10} {'Dual(0.15)':>10}")
    print(f"  {'-'*73}")
    for r in results:
        print(f"  {r['design']:<20} {r['sox2_mean']:>6.3f} {r['sox9_mean']:>6.3f}  "
              f"{r['thr10_dual_pct']:>10.1f}% {r['thr05_dual_pct']:>10.1f}% "
              f"{r['thr15_dual_pct']:>10.1f}%")
    print(f"{'='*75}")


if __name__ == "__main__":
    main()
