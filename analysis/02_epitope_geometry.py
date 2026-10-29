#!/usr/bin/env python3
"""
Check spatial relationship between SOX2 and SOX9 epitope patches.
Key question: are the epitopes close enough that binders fused to split-caspase-9
halves will reconstitute when both TFs are co-present?
Epitope centroid distances alone cannot establish split-caspase reconstitution.
"""
import numpy as np
from pathlib import Path
from Bio import PDB
import argparse


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--pdb", type=Path, required=True)
    cli.add_argument("--out-dir", type=Path, default=Path(__file__).resolve().parent)
    args = cli.parse_args()
    PDB_PATH = args.pdb
    OUT_DIR = args.out_dir
    if not PDB_PATH.is_file():
        cli.error(f"Structure not found: {PDB_PATH}")
    OUT_DIR.mkdir(parents=True, exist_ok=True)



    parser = PDB.PDBParser(QUIET=True)
    struct = parser.get_structure("sox", PDB_PATH)
    model  = struct[0]
    chain_a = model["A"]  # SOX2
    chain_b = model["B"]  # SOX9

    # Epitope residues from analysis
    SOX2_EPITOPE = [44, 59, 60, 62, 63, 66, 70, 74, 77, 81, 83]
    SOX9_EPITOPE = [103, 108, 141, 145, 147, 148, 149, 150, 151, 152,
                    155, 156, 159, 160, 163, 166, 167, 168, 169, 170, 171, 172]

    # SOX9 high-confidence patch only
    SOX9_PATCH = [r for r in SOX9_EPITOPE if r >= 141]

    def ca_coords(chain, resnums):
        coords = []
        for resnum in resnums:
            key = (" ", resnum, " ")
            if key not in chain or "CA" not in chain[key]:
                raise ValueError(f"Chain {chain.id}: missing CA at residue {resnum}")
            coords.append(chain[key]["CA"].get_vector().get_array())
        if not coords:
            raise ValueError("Epitope residue list is empty")
        return np.array(coords)


    def centroid(coords):
        return coords.mean(axis=0)

    sox2_ca  = ca_coords(chain_a, SOX2_EPITOPE)
    sox9_ca  = ca_coords(chain_b, SOX9_PATCH)

    sox2_cen = centroid(sox2_ca)
    sox9_cen = centroid(sox9_ca)

    dist = np.linalg.norm(sox2_cen - sox9_cen)

    print(f"SOX2 epitope centroid (A:{SOX2_EPITOPE}): {sox2_cen.round(1)}")
    print(f"SOX9 epitope centroid (B:{SOX9_PATCH}):  {sox9_cen.round(1)}")
    print(f"Inter-epitope centroid distance: {dist:.1f} Å")

    # Min distance between any pair of epitope Cα
    min_d, min_pair = 9999, None
    for i, (r2, c2) in enumerate(zip(SOX2_EPITOPE, sox2_ca)):
        for j, (r9, c9) in enumerate(zip(SOX9_PATCH, sox9_ca)):
            d = np.linalg.norm(c2 - c9)
            if d < min_d:
                min_d, min_pair = d, (r2, r9)

    print(f"Minimum Cα-Cα across epitopes: {min_d:.1f} Å  (SOX2:{min_pair[0]} — SOX9:{min_pair[1]})")

    print()
    print("Split-caspase-9 reconstitution geometry:")
    print(f"  Centroid separation: {dist:.0f} Å")
    if dist < 60:
        print("  Centroids are within 60 Å; activity and fusion geometry remain untested.")
        print("    Binder termini need to be oriented outward from epitopes")
    else:
        print("  Centroids exceed 60 Å; this descriptive distance alone cannot assess activity.")

    # All pairwise Cα distances for inter-epitope proximity map
    print()
    print("Closest 10 SOX2-SOX9 epitope Cα pairs:")
    pairs = []
    for r2, c2 in zip(SOX2_EPITOPE, sox2_ca):
        for r9, c9 in zip(SOX9_PATCH, sox9_ca):
            pairs.append((np.linalg.norm(c2-c9), r2, r9))
    pairs.sort()
    for d, r2, r9 in pairs[:10]:
        print(f"  A:{r2:3d} — B:{r9:3d}: {d:.1f} Å")

    # Check N/C terminus positions (for understanding binder orientation)
    def terminus_coord(chain, term="N"):
        residues = [r for r in chain.get_residues() if r.id[0]==" " and "CA" in r]
        res = residues[0] if term == "N" else residues[-1]
        return res["CA"].get_vector().get_array(), res.id[1]

    sox2_n, sox2_n_rnum = terminus_coord(chain_a, "N")
    sox2_c, sox2_c_rnum = terminus_coord(chain_a, "C")
    sox9_n, sox9_n_rnum = terminus_coord(chain_b, "N")
    sox9_c, sox9_c_rnum = terminus_coord(chain_b, "C")

    print()
    print("Termini positions (for binder fusion orientation):")
    for label, coord, rnum, epitope_cen in [
        ("SOX2 N-term", sox2_n, sox2_n_rnum, sox2_cen),
        ("SOX2 C-term", sox2_c, sox2_c_rnum, sox2_cen),
        ("SOX9 N-term", sox9_n, sox9_n_rnum, sox9_cen),
        ("SOX9 C-term", sox9_c, sox9_c_rnum, sox9_cen),
    ]:
        d_to_sox2_ep = np.linalg.norm(coord - sox2_cen)
        d_to_sox9_ep = np.linalg.norm(coord - sox9_cen)
        print(f"  {label} (res {rnum:4d}): d_to_SOX2_ep={d_to_sox2_ep:.0f}Å  d_to_SOX9_ep={d_to_sox9_ep:.0f}Å")

    # Distance between termini relevant to linker design
    print()
    print("Distances relevant to split-caspase-9 linker design:")
    combos = [
        ("SOX2_Cterm — SOX9_Nterm", sox2_c, sox9_n),
        ("SOX2_Nterm — SOX9_Cterm", sox2_n, sox9_c),
        ("SOX2_Cterm — SOX9_Cterm", sox2_c, sox9_c),
        ("SOX2_Nterm — SOX9_Nterm", sox2_n, sox9_n),
    ]
    for label, c1, c2 in combos:
        print(f"  {label}: {np.linalg.norm(c1-c2):.0f} Å")


if __name__ == "__main__":
    main()
