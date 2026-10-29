# Bridge Binder Rescue + Caspase Redesign Plan — 2026-09-06

## Diagnosis

| Campaign | Trajectories | Accepted | Root cause |
|----------|-------------|---------|------------|
| BindCraft bridge (frag) | 233 | 0 | i_pTM ~0.12–0.16; AF2 optimizer can't satisfy 5 hard hotspots across 2 chains simultaneously |
| BindCraft bridge (relaxed rerun) | ~92 | 0 | Same — all 20 MPNN seqs fail base AF2 filter for every trajectory |
| Protein Hunter bridge | 50 | 0 | best_ipTM ~0.46, pLDDT ~0.35 — structure disordered, not folded |
| Caspase Stage 5 flipped | 36 seeds | 1 at 8.1 Å | SOX9 binder C-term points away from SOX2 (d_toward = −0.76 to −0.92); only N-terminal caspase fusion geometry works |

MD replicates for lead construct (BC mpnn17 × ph_sox9_run18, GGS×1, flipped):
- rep1: 76.3% < 15 Å ✓  rep2: 80.1% ✓  rep3: 87.8% ✓  Mean: 81.4% — ROBUST PASS

---

## 5 Approaches (ordered by GPU priority)

### Phase 1 — RFdiffusion dual-chain bridge (30–60 min, first on GPU)
Script: `rfdiffusion/run_rfdiffusion_bridge_both.sh`
Output: `rfdiffusion/bridge_both_chains/`

Design a SINGLE binder targeting both chain A (SOX2, res 175-215) and chain B (SOX9, res 450-509)
simultaneously. RFdiffusion handles multi-chain targets far better than BindCraft's AF2 optimizer.

- Contig: `A175-215/0 B450-509/0 60-100`
- Hotspots: `A193,A199,B466,B505,B508`
- 200 designs → filter by Boltz-2 ipTM ≥ 0.5 + contacts to both chains

### Phase 2 — Boltz-2 two-chain linker screen (1–2 hrs, CPU-light, parallel with Phase 1)
Script: `selection/bridge_linker_screen.py`
Output: `selection/bridge_linker_outputs/`

Fuse the top stage-3-selective SOX2+SOX9 binder pair into a single chain with GGS linkers.
Test both orientations and 3 linker lengths:
- Pairs tested: mpnn17×ph_sox9_run1, mpnn18×ph_sox9_run4, ph_sox2_run18×ph_sox9_run18
- Linkers: GGS, GGSGGS, GGSGGSGGS (×1/×2/×3)
- Orientations: SOX2binder→SOX9binder and reversed
- Boltz-2 3-chain prediction (SOX2_HMG + SOX9_HMG + fused_binder)
- Measure: ipTM, individual-chain contacts, pLDDT of fused chain

Top sequences:
- ph_sox2_run18 (97 AA): SRLEELVDHIVDSMEEANQGGLTDEAREMIRNNLLEDLLDVNPNVEAELDAASDAEFAELLAGVDRIIREGTTSDLLGFLEALAEGGVAGALAELAK
- SOX2_mpnn17 (71 AA): PYEIPEEDKEWDEYSFAMVEVHERLMEAERQLRLLSPEEEEKKRELEEEIERARKELEEVHKRGEAIERKM
- SOX2_mpnn18 (71 AA): PYEIPEEDREWDEYSFRMVEVHERLMEARRELRLLAPEETEKKRELEEEIARAEKELEEVHKEGEAIERKM
- ph_sox9_run1 (77 AA): PPLSPEEQELEEIEEGYRNLTEKSGEALLKALGITSKNEEEREKKMRELALELAREVTRKLSEKTGVRYAEDGTRLP
- ph_sox9_run4 (86 AA): SEGIELMKELGLPVKDGKPYIEDLENLPELTNEEKLELGLLYGELNGGEKKAREAVEKAKKEGTEEAIERAFEEAVEAGDEGLKKK
- ph_sox9_run18 (67 AA): MEEVRELLERARELEEKVRELEDDPSKREERERLLEEAERLREEALDLLVEQIREPLLQWMREEGLL

### Phase 3 — BindCraft no hotspot constraint (overnight, after Phase 1 frees GPU)
Script: `bindcraft/run_bridge_nohot.sh`
Settings: `SOX2_SOX9_bridge_nohot.json` (no target_hotspot_residues)
Filter: `very_relaxed_bridge_filters.json` (i_pTM ≥ 0.20, pLDDT ≥ 0.70)
Output: `bindcraft/designs/sox2_sox9_bridge_nohot/`

Without hotspot pinning, BindCraft's optimizer is free to find natural contact geometry.
The very-relaxed filter (i_pTM 0.20 vs 0.35 before) allows marginal interfaces to survive to MPNN.

### Phase 4 — SOX9 N-terminal redesign (after Phase 1 structures available)
Scripts: `rfdiffusion/run_rfdiffusion_sox9_nterm.sh` + `bindcraft/run_sox9_nterm.sh`
Output: `rfdiffusion/sox9_nterm/` + `bindcraft/designs/sox9_nterm/`

Goal: produce SOX9 binders where the N-terminus is on the SOX2-PROXIMAL face.
This is required for CaspC+GGS+BinderSOX9 fusion to work geometrically.

Strategy: bias hotspots to the N-terminal half of SOX9 HMG fragment (B450–B470).
Hotspots: A193,A199 (SOX2 anchor for co-design) OR B453,B457,B463,B468 (SOX9 N-term face only).

After design: run terminus-direction diagnosis (Stage 5 script) to select designs
with d_toward > 0 (N-term pointing toward SOX2 binder).

### Phase 5 — Caspase geometry + MD (after Phase 4)
Workflow:
1. Take new SOX9 N-terminal binders
2. Run Stage 5 flipped caspase geometry screen (selection/stage5_flipped.py)
3. MD replicates on best hits (target: min_dist < 15 Å in >50% frames)

---

## GPU sequencing
RFdiffusion bridge (~30 min) → BindCraft no-hotspot (overnight) → RFdiffusion SOX9 N-term (~30 min) → BindCraft SOX9 N-term

## Key insight
The fundamental problem with BindCraft bridge is NOT the filter thresholds — it's that
the AF2-based trajectory optimizer achieves i_pTM ~0.12 during diffusion, meaning it never
generates a good enough backbone to pass even the relaxed base filter. RFdiffusion (backbone
diffusion only, no sequence) is much better suited to generate viable multi-chain bridging geometries.
