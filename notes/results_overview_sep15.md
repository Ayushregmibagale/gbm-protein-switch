# GBM Protein Switch — Results Overview (2026-09-15)

## Project goal
AND-gate protein kill switch for glioblastoma stem cells (GSCs).
Logic: SOX2-binder + SOX9-binder, each fused to a split-caspase-9 half.
When both SOX2 AND SOX9 are co-expressed (as in GSCs), binders engage simultaneously → caspase reconstitution → apoptosis.

---

## Design approaches (parallel)

### Track 1 — Split-caspase construct (lead)
Two separate binders fused by a linker to split-caspase-9 halves:
- Chain C: mpnn17 (SOX2-binder, 71 aa) + GGS + CaspN (aa140–315)
- Chain D: CaspC (aa316–416) + GGS + ph_sox9_run18 (SOX9-binder, 67 aa)
- SOX2 and SOX9 targets held in solution; binding brings caspase halves together.

### Track 2 — RFdiffusion bridge binder
Single de-novo protein designed to simultaneously contact both SOX2 and SOX9.
N-term fused to CaspN, C-term fused to CaspC → binding bridges caspase halves.

---

## Computational results

### Individual binder quality (Stage 1 re-prediction)
| Binder | Target | iptm | Selective? |
|--------|--------|------|-----------|
| mpnn17 (SOX2_mini_l71) | SOX2 | 0.897 | Yes (2.77× vs SOX9) |
| ph_sox9_run18 | SOX9 | 0.856 | Yes (6.69× vs SOX2) |
| ph_sox2_run18 | SOX2 | 0.784 | Yes (3.21× vs SOX9) |

### AND-gate geometry — Stage 2 (binder pair terminus distance)
Top pairs passing (both_selective=True), from Stage 4 composite ranking:
| Rank | Pair | Terminus dist (Å) | Composite score |
|------|------|-------------------|-----------------|
| 1 | ph_sox2_run18 × ph_sox9_run18 | 35.0 | 5.996 |
| 2 | mpnn18 × ph_sox9_run18 | 31.4 | 5.406 |
| 3 | mpnn17 × ph_sox9_run1 | 24.9 | 5.406 |
| 4 | **mpnn17 × ph_sox9_run18 (selected lead)** | 33.1 | 5.264 |

### Split-caspase geometry — Stage 5 flipped (Boltz-2, CaspN/C fused)
| Pair | Seeds passing <15 Å | Best dist | Notes |
|------|---------------------|-----------|-------|
| mpnn17 × ph_sox9_run18, GGS×1 | 1/3 (8.05 Å) | 8.05 Å | lead |
| mpnn17 × ph_sox9_run18, GGS×3 | 1/3 (10.35 Å) | 10.35 Å | |
| mpnn17 × ph_sox9_run1, GGS×1 | 1/3 (14.54 Å) | 14.54 Å | |
| ph_sox2_run18 × ph_sox9_run18 | NOT RUN | — | rank 1 pair, gap |

### GROMACS MD — split-caspase lead construct (3 replicates × 20 ns)
| Rep | % frames <15 Å (inter-half) | Pass |
|-----|----------------------------|------|
| rep1 | 76.3% | ✓ |
| rep2 | 80.1% | ✓ |
| rep3 | 87.8% | ✓ |
| **Mean** | **81.4%** | **✓** |

The lead construct maintains reconstitution-compatible geometry in all 3 MD replicates; functional caspase-9 activity is not yet demonstrated.

### GROMACS MD — bridge binder campaign (32 designs × 20 ns, 3 replicates for 5 best)
5/32 designs pass (15.6%). Top performers:

| Design | Rep1 | Rep2 | Rep3 | SOX2 contact | SOX9 contact |
|--------|------|------|------|-------------|-------------|
| bridge_48_s01 | 98.5% | 98.6% | 100.0% | 0.376 | 0.166 |
| bridge_73_s09 | 86.7% | 63.4% | 90.2% | 0.218 | 0.229 |
| bridge_136_s01 | 86.0% | 73.1% | 74.6% | 0.143 | 0.172 |
| bridge_101_s06 | 79.3% | 100.0% | 79.0% | 0.274 | 0.255 |
| bridge_90_s01 | 55.5% | 87.1% | 60.1% | 0.234 | 0.162 |

bridge_48_s01 has outstanding N-C geometry but lower SOX9 contact (0.166).
**Open question**: does it contact both targets simultaneously in each frame?

---

## Gaps before lab

### Computational (active)
1. **Dual-contact analysis for bridge_48_s01** (running: `analysis/bridge_dual_contact.py`)
   - Per-frame AND-gate: bridge contacts SOX2 AND SOX9 simultaneously?
   - Critical for validating Track 2 AND-gate mechanism
2. **Stage5 Boltz-2 for ph_sox2_run18 × ph_sox9_run18** (running: `selection/run_stage5_ph_sox2_run18.py`)
   - Rank-1 composite pair was never put through caspase geometry
   - May reveal a better lead than current mpnn17 pair

### For lab readiness
3. **Which track to pursue first** — decision needed once dual-contact analysis done
   - Track 1 (split-caspase): stronger end-to-end computational evidence
   - Track 2 (bridge binder): better N-C geometry but AND-gate unverified
4. **Individual binder affinities** — no Kd estimates; lab needs SPR/MST
5. **Full caspase construct synthesis** — sequence ready; needs gene synthesis + cloning
6. **DEVDase assay** — validate caspase reconstitution in vitro before cell assays

---

## Lead construct sequence (Track 1)
- SOX2-binder (mpnn17): `PYEIPEEDKEWDEYSFAMVEVHERLMEAERQLRLLSPEEEEKKRELEEEIERARKELEEVHKRGEAIERKM`
- SOX9-binder (ph_sox9_run18): `MEEVRELLERARELEEKVRELEDDPSKREERERLLEEAERLREEALDLLVEQIREPLLQWMREEGLL`
- Linker: GGS (3 aa)
- Split site: caspase-9 D315/G316 (P55211); N-half aa140–315, C-half aa316–416
- Full construct: [mpnn17]–GGS–[CaspN aa140–315] + [CaspC aa316–416]–GGS–[ph_sox9_run18]
