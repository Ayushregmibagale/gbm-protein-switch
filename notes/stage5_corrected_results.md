# Stage 5 — Corrected Results & Diagnosis (Sep 2 2026)

## Bug found and fixed
The original Stage 5 inter-half distance was measured WRONG:
- Was measuring: C-term(CaspN, chain C) ↔ N-term(BinderSOX9, chain D)  ← meaningless
- Should be:      C-term(CaspN, chain C) ↔ N-term(CaspC, chain D)       ← reconstitution point

The N-term of chain D is the START of BinderSOX9, not the start of CaspC. CaspC starts at
index len(BinderSOX9) + len(linker) into chain D.

The "21 Å best result" (BC mpnn17 × ph_sox9_run30, GGS×3) was actually 53.6 Å when measured correctly.

---

## Corrected reconstitution distances (all 20 runs)

| Pair | GGS | Correct recon dist | Pass <15Å? |
|------|-----|--------------------|------------|
| **BC mpnn18 × ph_sox9_run18** | **×1** | **8.3 Å** | **✓ YES** |
| BC mpnn18 × ph_sox9_run18 | ×4 | 22.4 Å | ✗ |
| BC mpnn17 × ph_sox9_run18 | ×3 | 28.6 Å | ✗ |
| BC mpnn17 × ph_sox9_run30 | ×1 | 29.8 Å | ✗ |
| BC mpnn17 × ph_sox9_run1  | ×3 | 32.1 Å | ✗ |
| ph_sox2_run18 × ph_sox9_run18 | ×1 | 33.3 Å | ✗ |
| BC mpnn17 × ph_sox9_run1  | ×4 | 34.4 Å | ✗ |
| BC mpnn17 × ph_sox9_run18 | ×1 | 35.0 Å | ✗ |
| BC mpnn17 × ph_sox9_run30 | ×2 | 35.3 Å | ✗ |
| BC mpnn17 × ph_sox9_run18 | ×2 | 35.6 Å | ✗ |
| ph_sox2_run18 × ph_sox9_run18 | ×3 | 33.7 Å | ✗ |
| BC mpnn18 × ph_sox9_run18 | ×2 | 46.4 Å | ✗ |
| BC mpnn18 × ph_sox9_run18 | ×3 | 46.1 Å | ✗ |
| BC mpnn17 × ph_sox9_run18 | ×4 | 49.3 Å | ✗ |
| BC mpnn17 × ph_sox9_run1  | ×2 | 50.9 Å | ✗ |
| BC mpnn17 × ph_sox9_run1  | ×1 | 52.7 Å | ✗ |
| BC mpnn17 × ph_sox9_run30 | ×3 | 53.6 Å | ✗ |
| ph_sox2_run18 × ph_sox9_run18 | ×2 | 58.4 Å | ✗ |
| ph_sox2_run18 × ph_sox9_run18 | ×4 | 44.9 Å | ✗ |
| BC mpnn17 × ph_sox9_run30 | ×4 | 65.6 Å | ✗ |

---

## Lead construct identified
**BC mpnn18 × ph_sox9_run18, GGS×1 linker (3 aa)**
- Reconstitution distance: **8.3 Å** (well within 15 Å threshold)
- pair_iptm: 0.279
- Full names:
  - SOX2 binder: `SOX2_mini_l71_s878170_mpnn18_model2`
  - SOX9 binder: `ph_sox9_run18`

## Diagnosis — why GGS×1 wins
- In Stage 2, the C-term(SOX2-binder) ↔ C-term(SOX9-binder) distance is 50.8 Å (far apart)
- But GGS×1 + CaspN (176 aa) folds back and positions the split site at 8.3 Å from CaspC start
- SHORTER linkers outperform longer ones here because more slack → more non-productive conformations
- GGS×2,3,4 for this same pair gave 46, 46, 22 Å — GGS×1 is the clear winner

## Stage 2 filter was measuring the wrong distance for this design
Stage 2 measured C-term(BinderSOX2) ↔ N-term(BinderSOX9) — correct for an N+C mixed fusion design.
Our design uses C-terminal fusion for BOTH binders. The correct Stage 2 metric for this
design would have been C-term(C) ↔ C-term(D). For mpnn18 × ph_sox9_run18, that distance
is 50.8 Å in the binder model, yet GGS×1 works because CaspN folds back.

## Caveat
8.3 Å is from a SINGLE Boltz-2 seed. The model is stochastic — need multi-seed
confirmation to verify this isn't a lucky outlier.

## Next step
Multi-seed confirmation: run seeds 1–4 for BC mpnn18 × ph_sox9_run18, GGS×1.
If ≥3/4 seeds give <15 Å → proceed to GROMACS MD.
Script: `selection/stage5_multiseed.py`
Output: `selection/stage5_multiseed_results.tsv`
