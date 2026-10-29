# Hotspot Selection Rationale — SOX2 and SOX9 Binders

## Design constraint
Both binders must engage their targets **simultaneously while SOX2 and SOX9 are in complex**.
Two requirements follow:
1. Epitopes must be **exposed in the co-complex** (not buried at the SOX2-SOX9 PPI)
2. Epitopes must be on the **same face** of the co-complex so binder termini point toward each other → caspase halves can reconstitute

## AF3 co-complex model used for analysis
- Input: SOX2 HMG (82 aa, auth 39-120) + SOX9 HMG (72 aa, auth 103-174) + DNA (12+12 nt)
- Best model: seed-3_sample-3, ranking=0.703, iptm=0.63, chain_pair_iptm A-B=0.48
- Output: `~/gbm-protein-switch/structures/af3_cocomplex/hmg_output/SOX2_HMG_SOX9_HMG_DNA_20260828_174417/`
- AF3 residue numbering is 1-indexed from input; offsets: SOX2 +38, SOX9 +102

## SOX2-SOX9 PPI interface (5 Å cutoff, AF3 model)
- SOX2 side (auth): S69, K117, T118, L119
- SOX9 side (auth): E159, R162, K166, Y174

## Iteration history

### v1 (initial, from Boltz-2 model)
SOX9 hotspots: A148, A152, A159, A166, A171
Problem: A159 and A166 are at the SOX2-SOX9 PPI interface → blocked in co-complex

### v2 (PPI fix)
SOX9 hotspots: A148, A152, A156, A157, A171
Problem: 148, 152, 156, 157 are at angles 130–165° from SOX2 epitope → wrong face of co-complex → caspase halves point away from each other

### v3 (geometry fix — CURRENT)
SOX9 hotspots: **A168, A169, A171, A172, A173**
Rationale: all at angle 62–89° from SOX2 epitope centroid → same face of co-complex

| Auth | AA  | Angle to SOX2 epi | DNA dist | PPI safe |
|------|-----|-------------------|----------|----------|
| 168  | ASP |  89°              | 12.8 Å   | yes      |
| 169  | HIS |  80°              | 11.2 Å   | yes      |
| 171  | ASP |  68°              | 8.7 Å    | yes      |
| 172  | TYR |  69°              | 7.5 Å    | yes      |
| 173  | LYS |  62°              | 6.6 Å    | yes      |

BoltzGen label_seq_id (auth − 102): 66, 67, 69, 70, 71

## SOX2 hotspots — confirmed correct throughout
A62, A70, A74, A81 — HMG wing-2 face
- Not at SOX2-SOX9 PPI interface (interface is S69/K117-119)
- Inter-epitope centroid distance to SOX9 v3 cluster: ~33.5 Å → feasible for GGS5 bridging

## Figures
- `figures/cocomplex_hotspot_geometry.png` — AF3 co-complex CA trace (PC1/PC2 and PC2/PC3 views), hotspots marked, inter-epitope vector and distance shown
- `figures/andgate_mechanism_schematic.png` — AND-gate logic flow + construct architecture + cell-type specificity table

## Post-design filtering requirement
Hotspot geometry is necessary but not sufficient. Must also filter designs by:
- **Terminus proximity**: in quaternary ColabFold model (SOX2:SOX9:BinderA:BinderB),
  the terminus of Binder-SOX2 and terminus of Binder-SOX9 must be < ~40 Å apart
- **ipTM > 0.7** and **pLDDT > 70** for each binder in the co-complex context
- Test pairs of designs, not each binder individually
