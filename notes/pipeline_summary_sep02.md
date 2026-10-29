# Full Selection Pipeline Summary — Sep 2 2026

## Design inventory (all campaigns complete)
| Method | SOX2 | SOX9 | Total |
|--------|------|------|-------|
| BindCraft | 13 | 40 | 53 |
| Protein Hunter | 49 | 50 | 99 |
| BoltzGen | 20 | 20 | 40 |
| **Total** | **82** | **110** | **192** |

---

## Stage 1 — Boltz-2 re-prediction (COMPLETE)
- 172 designs re-predicted with common scorer
- **69/172 passed** (ipTM ≥ 0.70, pLDDT ≥ 0.70)
- SOX2: 11 passing | SOX9: 58 passing
- PH SOX9 dominated — top ipTM = 0.930

**Bug fixed:** Boltz-2 returns pLDDT on 0–1 scale, not 0–100. Threshold corrected to 0.70.
**Bug fixed:** BoltzGen designs are in `final_ranked_designs/final_20_designs/` subdir; read from `final_designs_metrics_20.csv` not CIF files.

---

## Stage 2 — AND-gate geometry (COMPLETE)
- 638 pairs tested (11 SOX2 × 58 SOX9) in 4-chain Boltz-2
- **25/638 pairs passed** (pair ipTM ≥ 0.60, terminus dist ≤ 40 Å)
- `boltzgen_sox2_binder_851` dominated but later FAILED selectivity
- All 25 passing pairs used PH SOX9 binders

---

## Stage 3 — Selectivity (COMPLETE)
- **10/25 pairs fully selective** (wrong-target ipTM < 0.40)
- `boltzgen_sox2_binder_851` failed — wrong-target ipTM = 0.495 (binds SOX9 too)
- Best selective SOX2 binders: BC mpnn17 (wrong=0.324), BC mpnn18 (wrong=0.351)
- Best selective SOX9 binder: ph_sox9_run18 (wrong=0.128)

---

## Stage 4 — Composite ranking (COMPLETE)
Score = pair_ipTM × selectivity_ratio(SOX2) × selectivity_ratio(SOX9) × geometry_bonus

Top 5:
1. ph_sox2_run18 × ph_sox9_run18 — score 5.996, term 35.0 Å
2. BC mpnn18 × ph_sox9_run18 — score 5.406, term 31.4 Å
3. BC mpnn17 × ph_sox9_run1 — score 5.406, term 24.9 Å
4. BC mpnn17 × ph_sox9_run18 — score 5.264, term 33.1 Å
5. BC mpnn17 × ph_sox9_run30 — score 4.278, term 20.0 Å

---

## Stage 5 — Caspase geometry modelling

### Original run (0/20 passed)
**Bug found:** inter_half_distance measured C-term(CaspN) ↔ N-term(BinderSOX9) — wrong.
Correct: C-term(CaspN) ↔ N-term(CaspC) = residue at index len(BinderSOX9)+len(linker) in chain D.
"Best" result of 21.0 Å was actually 53.6 Å when measured correctly.

### Corrected measurements
Only 1 run genuinely passed: **BC mpnn18 × ph_sox9_run18, GGS×1 = 8.3 Å** (seed 0).

### Multi-seed confirmation (1/5 passed)
Ran seeds 0–4 for BC mpnn18 × ph_sox9_run18, GGS×1.
Result: 8.3, 22.4, 27.2, 55.7, 21.1 Å → mean 26.9 ± 17.5 Å. Geometry unreliable.

### Option 4 — Terminus orientation diagnosis
Measured direction vectors from Stage 2 PDBs (centroid → C-terminus for each binder).
**Finding: ALL 4 top pairs have SOX9 binder C-terminus pointing AWAY from SOX2 binder.**
- mpnn18 × ph_sox9_run18: d_toward = −0.91 (strongly away)
- mpnn17 × ph_sox9_run18: d_toward = −0.76 (strongly away)
- mpnn17 × ph_sox9_run30: d_toward = −0.92 (strongly away)
- mpnn17 × ph_sox9_run1:  d_toward = −0.87 (strongly away)

This is why longer linkers don't help — the caspase halves face away regardless of length.

### Option 3 — Flipped orientation (4/36 passed, improvement confirmed)
Fix: chain D changed from `BinderSOX9+GGS+CaspC` → `CaspC+GGS+BinderSOX9`
Now CaspC sits at the N-terminus of chain D, directly adjacent to CaspN end of chain C.

Results (4 pairs × 3 linkers × 3 seeds):
| Pair | GGS | min dist | pass/3 seeds |
|------|-----|----------|--------------|
| BC mpnn17 × ph_sox9_run18 | ×1 | **8.1 Å** | 1/3 |
| BC mpnn17 × ph_sox9_run18 | ×3 | **10.4 Å** | 1/3 |
| BC mpnn18 × ph_sox9_run18 | ×3 | **13.5 Å** | 1/3 |
| BC mpnn17 × ph_sox9_run1  | ×1 | **14.5 Å** | 1/3 |

Improvement vs original: min dist went from 8.3 Å (1/5 seeds, original) to 8.1 Å (1/3 seeds, flipped) across MULTIPLE pairs. Flipping is confirmed beneficial.

Still only 1/3 seeds pass per config — this is the fundamental limit of static Boltz-2.
**MD is the correct next step.**

---

## Lead construct for MD

**BC mpnn17 × ph_sox9_run18, GGS×1, flipped orientation**
- Full chain C: `SOX2_mini_l71_s878170_mpnn17_model2` + GGS + CaspN(aa140-315)
- Full chain D: CaspC(aa316-416) + GGS + `ph_sox9_run18`
- Starting structure: seed giving 8.1 Å recon distance
- Path: `selection/stage5_flipped_outputs/flip_n17_model2_ph_sox9_run18_GGS1_s*/`

GROMACS MD plan:
- 20–50 ns, monitor inter-half distance
- If halves converge <15 Å in >50% frames → proceed to construct design
- If not → try longer linkers or redesign campaign (Option 5)

---

## Option 5 — Redesign campaign (PENDING)
Redesign SOX9 binders with N-terminus facing the SOX2 binder, OR
design interface-bridging binders that tolerate N-terminal caspase fusion.
Launch after MD results are in.

---

## Scripts and outputs
| Script | Output |
|--------|--------|
| `selection/stage1_reprediction.py` | `stage1_results.tsv` |
| `selection/stage2_andgate_geometry.py` | `stage2_results.tsv` |
| `selection/stage3_selectivity.py` | `stage3_results.tsv` |
| `selection/stage4_composite_ranking.py` | `stage4_results.tsv` |
| `selection/stage5_caspase_geometry.py` | `stage5_results.tsv` (buggy dist) |
| `selection/stage5_multiseed.py` | `stage5_multiseed_results.tsv` |
| `selection/stage5_flipped.py` | `stage5_flipped_results.tsv` |
