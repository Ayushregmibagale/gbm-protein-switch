# GBM Protein Switch — Campaign Progress (Sep 2 2026)

## Design inventory (all campaigns COMPLETE)

| Method | SOX2 | SOX9 |
|--------|------|------|
| BindCraft | 13 | 40 |
| Protein Hunter | 49 | 50 |
| BoltzGen | 20 | 20 |
| **Total** | **82** | **110** |

Grand total: **172 designs**

---

## Selection pipeline

### Stage 1 — Boltz-2 re-prediction (COMPLETE, 00:28 BST Sep 2)
- **69/172 passed** (ipTM ≥ 0.70, pLDDT ≥ 0.70)
- SOX2: 11 passing | SOX9: 58 passing
- Results: `selection/stage1_results.tsv`

Bugs fixed:
- BoltzGen designs: must read from `final_ranked_designs/final_designs_metrics_20.csv`, not glob of CIF files
- pLDDT scale: Boltz-2 returns 0–1, not 0–100 — threshold corrected to 0.70

### Stage 2 — AND-gate geometry filter (COMPLETE, 16:05 BST Sep 2)
- Tested **638 pairs** (11 SOX2 × 58 SOX9) in 4-chain Boltz-2
- **25/638 pairs passed** (pair ipTM ≥ 0.60 AND terminus dist ≤ 40 Å)
- Results: `selection/stage2_results.tsv`

Top 10 passing pairs:
| Rank | SOX2 binder | SOX9 binder | pair ipTM | term dist |
|------|-------------|-------------|-----------|-----------|
| 1 | boltzgen_sox2_binder_851 | ph_sox9_run48 | 0.915 | 31.7 Å |
| 2 | boltzgen_sox2_binder_851 | ph_sox9_run28 | 0.802 | 32.5 Å |
| 3 | BC mpnn18_model2 | ph_sox9_run4 | 0.741 | 33.6 Å |
| 4 | ph_sox2_run37 | ph_sox9_run48 | 0.736 | 37.9 Å |
| 5 | BC mpnn17_model2 | ph_sox9_run1 | 0.732 | 24.9 Å |
| 6 | boltzgen_sox2_binder_851 | ph_sox9_run45 | 0.731 | 33.9 Å |
| 7 | boltzgen_sox2_binder_851 | ph_sox9_run6 | 0.731 | 26.0 Å |
| 8 | BC mpnn17_model2 | ph_sox9_run30 | 0.730 | 20.0 Å |
| 9 | BC mpnn2_model2 | ph_sox9_run45 | 0.724 | 33.0 Å |
| 10 | boltzgen_sox2_binder_851 | ph_sox9_run46 | 0.721 | 37.7 Å |

Key findings:
- `boltzgen_sox2_binder_851` is the best SOX2 binder (6/10 top pairs)
- ALL 25 passing pairs use PH SOX9 binders — PH dominated SOX9 design
- Terminus distances are well under 40 Å threshold for top pairs → good caspase geometry

### Stage 3 — Selectivity filter (RUNNING, launched 16:18 BST Sep 2)
- Script: `selection/stage3_selectivity.py`
- Tests each unique binder (8 SOX2 + 16 SOX9 = 24 total) against wrong target
- Threshold: wrong-target ipTM < 0.40
- ETA: ~16:55 BST Sep 2

### Stage 4 — Composite ranking (TODO)
- score = ipTM × ShapeComplementarity × exp(-ddG/kT)
- Take top 5 pairs per source combination

### Stage 5 — Caspase geometry modelling (TODO)
- Boltz-2 full assembly: SOX2 + SOX9 + Binder-GGS-CaspHalf-A + Binder-GGS-CaspHalf-B
- GGS linker sweep: ×1/2/3/4
- Check: caspase catalytic residues not buried, inter-half distance < 15 Å
- GROMACS MD 20–50 ns on best pose
- Plan: `notes/caspase_geometry_modelling_plan.md`

---

## Unique binders in top pairs

**SOX2 binders (8 unique):**
- boltzgen_sox2_sox2_binder_851 ← star performer
- SOX2_mini_l71_s878170_mpnn17_model2 (BindCraft)
- SOX2_mini_l71_s878170_mpnn18_model2 (BindCraft)
- SOX2_mini_l79_s924826_mpnn2_model2 (BindCraft)
- SOX2_mini_l79_s924826_mpnn3_model2 (BindCraft)
- ph_sox2_run14, ph_sox2_run18, ph_sox2_run37

**SOX9 binders (16 unique, all PH):**
- ph_sox9_run48, ph_sox9_run28, ph_sox9_run4, ph_sox9_run1 ← top 4

---

## Logs
- Stage 1: `selection/logs/stage1.log`
- Stage 2: `selection/logs/stage2.log`
- Stage 3: `selection/logs/stage3.log`
