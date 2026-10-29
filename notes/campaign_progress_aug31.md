# GBM Protein Switch — Campaign Progress (Aug 31 2026)

## Project summary
Split-caspase-9 AND gate targeting SOX2/SOX9 co-complex in glioblastoma stem cells.
Inspired by Koeber et al. Nature 2026. Both binders must be co-present for caspase reconstitution.

Construct architecture:
```
NLS-[BinderSOX2]-GGS5-CaspHalf1-T2A-NLS-[BinderSOX9]-GGS5-CaspHalf2
```

---

## Structures used
| Target | PDB | Chain | Resolution |
|--------|-----|-------|------------|
| SOX2 | 6WX8 | B | 2.3 Å |
| SOX9 | 4EUW | A | 2.77 Å |

Co-complex geometry verified with AF3 HMG-only model (seed-3_sample-3, iptm=0.63).

---

## Hotspot selection (v3, final)

### SOX2 hotspots: B62, B70, B74, B81
- HMG wing-2 face (H2-H3 helix)
- Not at SOX2-SOX9 PPI interface (interface = S69/K117-119)
- Safe throughout all hotspot iterations

### SOX9 hotspots: A168, A169, A171, A172, A173
- C-terminal cluster of HMG domain
- Angles to SOX2 epitope centroid: 62–89° → same face → caspase halves can meet
- Clear of PPI interface (E159/R162/K166/Y174 are the interface residues)
- BoltzGen label_seq_id (auth−102): 66, 67, 69, 70, 71

### Inter-epitope geometry
- Centroid-centroid distance: 33.5 Å → feasible for GGS5 bridging (~18 Å reach per arm)

### Iteration history (SOX9)
- v1: A148,A152,A159,A166,A171 — FAILED: A159/A166 at PPI interface
- v2: A148,A152,A156,A157,A171 — FAILED: angles 130–165° (wrong face, caspase halves point away)
- v3: A168,A169,A171,A172,A173 — CORRECT: angles 62–89°

---

## Design campaigns

### BindCraft (AF2-based hallucination)

| Campaign | Hotspots | Filter | Accepted | Status |
|----------|----------|--------|----------|--------|
| SOX2 xtal | B62,B70,B74,B81 | sox2_relaxed | **13** | DONE |
| SOX9 xtal v3 | A168,A169,A171,A172,A173 | default | **40** | DONE |
| SOX2 groove | B59,B61,B62,B63,B70 | sox2_relaxed | 0 | NOT LAUNCHED |

SOX2 relaxed filters (`sox2_relaxed_filters.json`):
- ShapeComplementarity: avg ≥0.45 (default 0.60) — flat helix face scores 0.4–0.5 inherently
- n_InterfaceHbonds: avg ≥2 (default 3)
- i_pAE: avg ≤0.42 (default 0.35)

SOX9 accepted 40 designs in 88 trajectories (45% rate) — concave surface scores well.
SOX2 accepted 13 in ~2700 trajectories (0.5% rate) — flat surface, hard.

### Protein Hunter (Boltz-2 + LigandMPNN cycles)
- 50 designs, 7 cycles, ipTM/pLDDT threshold 0.7, single-sequence MSA
- Bug fixed: `export MPLBACKEND=Agg` added to both run scripts (matplotlib headless crash)

| Campaign | Status |
|----------|--------|
| PH SOX2 | **RUNNING** (Run ~24/49, ~13:30 BST ETA) |
| PH SOX9 | Queued (~15:30 BST ETA) |

### BoltzGen
- YAML-based binder design using Boltz-2
- Failed in first pipeline pass — errors TBD when queue reaches it
- Specs: `boltzgen/specs/sox2_binder.yaml` and `sox9_binder.yaml`

---

## Key bugs fixed this campaign

| Bug | Fix |
|-----|-----|
| SOX2 BindCraft DSSP segfault (WSL2 vsyscall) | Watchdog wrapper auto-restarts (`run_sox2_watchdog.sh`) |
| Rogue processes using wrong filter (default vs relaxed) | Killed PIDs 1531796, 1711865 |
| PH matplotlib GUI crash in headless env | `export MPLBACKEND=Agg` in both PH scripts |
| PH SOX2 target → 20 (too slow, flat surface) | Lowered target to 13; unblocked pipeline |

---

## Pipeline sequence
```
SOX2 BindCraft (DONE, 13) ──┐
                             ├──→ PH SOX2 → PH SOX9 → BoltzGen SOX2 → BoltzGen SOX9
SOX9 BindCraft (DONE, 40) ──┘
```
Pipeline log: `~/gbm-protein-switch/bindcraft/logs/pipeline_launcher_v6.log`

---

## Post-design validation (TODO — not started)
1. Quaternary ColabFold: SOX2:SOX9:BinderA:BinderB
2. Filter pairs where BinderA terminus ↔ BinderB terminus < 40 Å
3. AF3 or pRosettaC quaternary docking
4. GROMACS MD on top candidates

## SOX2 groove campaign (TODO — not started)
- Script ready: `run_sox2_groove_watchdog.sh`
- Settings: `SOX2_groove_miniprotein.json`, hotspots B59,B61,B62,B63,B70
- Launch when GPU is free after PH/BoltzGen completes
