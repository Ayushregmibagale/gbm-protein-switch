# Design Selection Pipeline

## Why three methods need a shared funnel
BindCraft, Protein Hunter, and BoltzGen each use different scoring functions. You cannot compare their raw scores directly. All designs must pass through the same scorer before ranking.

## Pipeline overview

```
BindCraft (13 SOX2 + 40 SOX9)  ──┐
PH SOX2 + SOX9                  ──┼──→ Stage 1 → Stage 2 → Stage 3 → Stage 4 → MD
BoltzGen SOX2 + SOX9            ──┘
```

---

## Stage 1 — Boltz-2 re-prediction (normalise)
**Script:** `selection/stage1_reprediction.py`
**Output:** `selection/stage1_results.tsv`

For every design from every method, run a fresh Boltz-2 prediction (binder + target HMG).
Pass threshold: ipTM ≥ 0.70 AND pLDDT ≥ 70

Why Boltz-2: it's already installed and is the same model used by PH, making it the fairest common scorer.
Expected pass rate: ~30–50% (BindCraft tends to do best here since its acceptance filter already requires good AF2 scores).

---

## Stage 2 — AND-gate geometry filter
**Script:** `selection/stage2_andgate_geometry.py`
**Output:** `selection/stage2_results.tsv`

Tests pairs of designs (one SOX2 binder + one SOX9 binder) in a 4-chain Boltz-2 prediction:
Chain A = SOX2 HMG, Chain B = SOX9 HMG, Chain C = BinderSOX2, Chain D = BinderSOX9

Pass threshold:
- pair ipTM ≥ 0.60
- C-terminus of BinderSOX2 ↔ N-terminus of BinderSOX9 < 40 Å (caspase reconstitution geometry)

Expected pass rate: ~20–30% of pairs.

---

## Stage 3 — Selectivity check (TODO)
For surviving pairs, verify binders are selective:
- BinderSOX2 vs SOX2 alone: should score well (ipTM > 0.7)
- BinderSOX2 vs SOX9 alone: should score poorly (ipTM < 0.4)
- BinderSOX9 vs SOX9 alone: should score well
- BinderSOX9 vs SOX2 alone: should score poorly

Discard any binder with cross-reactivity.

---

## Stage 4 — Composite ranking (TODO)
For all Stage 3 passing pairs, compute composite score:
```
score = ipTM × ShapeComplementarity × exp(-ddG / kT)
```
Take top 5 pairs per source combination (BC/BC, BC/PH, BC/BG, PH/PH, etc.)

---

## Stage 5 — Caspase geometry modelling (see caspase_geometry_modelling_plan.md)
Top pairs go to full-assembly Boltz-2/AF3 prediction with GGS linkers and caspase halves.
Linker length sweep: GGS×1, ×2, ×3, ×4.
Best poses → GROMACS MD 20–50 ns.

---

## Running the pipeline
```bash
cd ~/gbm-protein-switch/selection

# Dry run to check all design files are found
conda activate boltz2
python stage1_reprediction.py --dry_run

# Full Stage 1 (run after BoltzGen completes)
nohup python stage1_reprediction.py > logs/stage1.log 2>&1 &

# Stage 2 (run after Stage 1)
nohup python stage2_andgate_geometry.py > logs/stage2.log 2>&1 &
```

## Target sequences used
- SOX2: 6WX8 chain B HMG domain (auth 39-120, extracted from crystal structure)
- SOX9: 4EUW chain A HMG domain (auth 103-174, 72 aa):
  PHVKRPMNAFMVWAQAARRKLADQYPHLHNAELSKTLGKLWRLLNESEKRPFVEEAERLRVQHKKDHPDYKY
