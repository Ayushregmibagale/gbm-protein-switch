#!/usr/bin/env python3
"""
Stage 4: Composite ranking of selective pairs.

Score = pair_iptm
      × selectivity_ratio_sox2  (correct_iptm / wrong_iptm)
      × selectivity_ratio_sox9  (correct_iptm / wrong_iptm)
      × geometry_bonus           exp(-(terminus_dist / 40))

Higher score = better. Output ranked TSV + human-readable summary.
"""

import pandas as pd
import numpy as np
from pathlib import Path

BASE   = Path("/mnt/d/gbm-protein-switch")
STAGE1 = BASE / "selection" / "stage1_results.tsv"
STAGE3 = BASE / "selection" / "stage3_results.tsv"

GEOM_SCALE = 40.0   # Å — terminus distance scale for geometry bonus

def main():
    df1 = pd.read_csv(STAGE1, sep="\t")
    df3 = pd.read_csv(STAGE3, sep="\t")

    sel = df3[df3["both_selective"] == True].copy()
    print(f"=== Stage 4: Composite ranking ===")
    print(f"  Selective pairs to rank: {len(sel)}")

    # Selectivity ratios (correct / wrong — higher = more specific)
    sel["sox2_sel_ratio"] = sel["sox2_correct_iptm"] / sel["sox2_wrong_iptm"]
    sel["sox9_sel_ratio"] = sel["sox9_correct_iptm"] / sel["sox9_wrong_iptm"]

    # Geometry bonus: closer termini = better for caspase reconstitution
    sel["geometry_bonus"] = np.exp(-sel["terminus_dist_A"] / GEOM_SCALE)

    # Composite score
    sel["composite_score"] = (
        sel["pair_iptm"]
        * sel["sox2_sel_ratio"]
        * sel["sox9_sel_ratio"]
        * sel["geometry_bonus"]
    )

    sel = sel.sort_values("composite_score", ascending=False).reset_index(drop=True)
    sel.index += 1  # 1-based rank

    out_path = BASE / "selection" / "stage4_results.tsv"
    sel.to_csv(out_path, sep="\t", index_label="rank")

    print(f"\n{'Rank':<5} {'SOX2 binder':<42} {'SOX9 binder':<18} {'pair_ipTM':<11} {'term(Å)':<9} {'S2_ratio':<10} {'S9_ratio':<10} {'geom':<7} {'SCORE'}")
    print("-" * 130)
    for rank, row in sel.iterrows():
        print(f"{rank:<5} {row['sox2_name']:<42} {row['sox9_name']:<18} "
              f"{row['pair_iptm']:<11.3f} {row['terminus_dist_A']:<9.1f} "
              f"{row['sox2_sel_ratio']:<10.2f} {row['sox9_sel_ratio']:<10.2f} "
              f"{row['geometry_bonus']:<7.3f} {row['composite_score']:.4f}")

    print(f"\nResults saved: {out_path}")
    print(f"\n=== Top 5 for caspase geometry modelling ===")
    for rank, row in sel.head(5).iterrows():
        print(f"  {rank}. {row['sox2_name']} × {row['sox9_name']}")
        print(f"     pair_ipTM={row['pair_iptm']:.3f}  term={row['terminus_dist_A']:.1f}Å  score={row['composite_score']:.4f}")

if __name__ == "__main__":
    main()
