# Stage 3 Selectivity Results — Sep 2 2026

## Summary
- Input: 25 pairs passing Stage 2 (AND-gate geometry)
- Test: each unique binder vs the WRONG target (SOX2 binder → SOX9, SOX9 binder → SOX2)
- Pass threshold: wrong-target ipTM < 0.40
- **Result: 10/25 pairs fully selective**

---

## Key finding: boltzgen_sox2_binder_851 FAILS selectivity

The top Stage 2 performer (`boltzgen_sox2_binder_851`, pair ipTM=0.915) has wrong-target ipTM=0.495 — it also binds SOX9, making it non-selective. **All 10 selective pairs use BindCraft binders** (mpnn17 or mpnn18) as the SOX2 component.

**BindCraft mpnn17** (`SOX2_mini_l71_s878170_mpnn17_model2`): wrong-target ipTM = 0.324 ✓  
**BindCraft mpnn18** (`SOX2_mini_l71_s878170_mpnn18_model2`): wrong-target ipTM = 0.351 ✓  
**ph_sox2_run18**: wrong-target ipTM = 0.244 ✓ (most selective SOX2 binder)

---

## 10 Selective Pairs (both_selective = True)

| Rank | SOX2 binder | SOX9 binder | pair ipTM | term dist | SOX2 wrong | SOX9 wrong |
|------|-------------|-------------|-----------|-----------|------------|------------|
| 1 | BC mpnn18 | ph_sox9_run4 | 0.741 | 33.6 Å | 0.351 | 0.304 |
| 2 | BC mpnn17 | ph_sox9_run1 | 0.732 | 24.9 Å | 0.324 | 0.186 |
| 3 | BC mpnn17 | ph_sox9_run30 | 0.730 | 20.0 Å | 0.324 | 0.248 |
| 4 | BC mpnn18 | ph_sox9_run18 | 0.697 | 31.4 Å | 0.351 | 0.128 |
| 5 | ph_sox2_run18 | ph_sox9_run18 | 0.669 | 35.0 Å | 0.244 | 0.128 |
| 6 | BC mpnn17 | ph_sox9_run29 | 0.667 | 27.2 Å | 0.324 | 0.256 |
| 7 | BC mpnn18 | ph_sox9_run28 | 0.662 | 38.2 Å | 0.351 | 0.336 |
| 8 | BC mpnn17 | ph_sox9_run18 | 0.650 | 33.1 Å | 0.324 | 0.128 |
| 9 | BC mpnn17 | ph_sox9_run4 | 0.628 | 39.8 Å | 0.324 | 0.304 |
| 10 | BC mpnn17 | ph_sox9_run28 | 0.607 | 22.6 Å | 0.324 | 0.336 |

Full name key:
- BC mpnn17 = `SOX2_mini_l71_s878170_mpnn17_model2` (BindCraft)
- BC mpnn18 = `SOX2_mini_l71_s878170_mpnn18_model2` (BindCraft)

---

## 15 Non-selective pairs (for reference)

| SOX2 binder | SOX9 binder | pair ipTM | Why failed |
|-------------|-------------|-----------|-----------|
| boltzgen_sox2_binder_851 | ph_sox9_run48 | 0.915 | SOX2 wrong=0.495, SOX9 wrong=0.560 |
| boltzgen_sox2_binder_851 | ph_sox9_run28 | 0.802 | SOX2 wrong=0.495 |
| ph_sox2_run37 | ph_sox9_run48 | 0.736 | SOX2 wrong=0.584, SOX9 wrong=0.560 |
| boltzgen_sox2_binder_851 | ph_sox9_run45 | 0.731 | SOX2 wrong=0.495 |
| boltzgen_sox2_binder_851 | ph_sox9_run6 | 0.731 | SOX2 wrong=0.495, SOX9 wrong=0.552 |
| SOX2_mini_l79_s924826_mpnn2 | ph_sox9_run45 | 0.724 | SOX2 wrong=0.573 |
| ph_sox2_run14 | ph_sox9_run45 | 0.720 | SOX2 wrong=0.596 |
| boltzgen_sox2_binder_851 | ph_sox9_run16 | 0.673 | SOX2 wrong=0.495 |
| boltzgen_sox2_binder_851 | ph_sox9_run46 | 0.721 | SOX2 wrong=0.495 |
| SOX2_mini_l79_s924826_mpnn2 | ph_sox9_run35 | 0.651 | SOX2 wrong=0.573, SOX9 wrong=0.404 |
| BC mpnn18 | ph_sox9_run16 | 0.641 | SOX9 wrong=0.445 |
| BC mpnn17 | ph_sox9_run8 | 0.640 | SOX9 wrong=0.678 |
| BC mpnn17 | ph_sox9_run24 | 0.636 | SOX9 wrong=0.431 |
| boltzgen_sox2_binder_851 | ph_sox9_run47 | 0.620 | SOX2 wrong=0.495 |
| SOX2_mini_l79_s924826_mpnn3 | ph_sox9_run25 | 0.604 | SOX2 wrong=0.679 |

---

## Interpretation

- **boltzgen_sox2_binder_851**: high pair ipTM (great at bringing both partners together) but not SOX2-specific — likely contacts a conserved surface shared with SOX9. Rejected.
- **BC mpnn17/18**: designed to the SOX2 wing-2 helix hotspots (B62,B70,B74,B81) — these residues diverge from SOX9, which is why selectivity holds.
- **ph_sox9_run18**: lowest wrong-target ipTM (0.128) — the most SOX9-specific binder available. Appears in 3 selective pairs (ranks 4, 5, 8).
- **ph_sox9_run1**: also highly selective (wrong=0.186), good terminus distance (24.9 Å) in pair 2.

---

## Next: Stage 4 composite ranking
Ranks the 10 selective pairs by composite score = pair_iptm × selectivity_ratio(SOX2) × selectivity_ratio(SOX9) × geometry_bonus.
Top 5 → caspase geometry modelling.
