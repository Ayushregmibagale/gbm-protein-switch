# GBM Project — Results & Lab Readiness — 2026-09-06

## Summary of all results to date

### Stage 1 — Boltz-2 re-prediction
- 172 designs re-predicted; **69/172 passed** (ipTM ≥ 0.70, pLDDT ≥ 0.70)
- SOX2: 11 passing | SOX9: 58 passing

### Stage 2 — AND-gate geometry
- 638 pairs tested (11 SOX2 × 58 SOX9); **25/638 passed** (pair ipTM ≥ 0.60, terminus dist ≤ 40 Å)

### Stage 3 — Selectivity
- **10/25 pairs fully selective** (wrong-target ipTM < 0.40 for both binders)
- Top pair: mpnn18 × ph_sox9_run4 (ipTM=0.741, sox2_wrong=0.351, sox9_wrong=0.304)

### Stage 4 — Composite ranking
- Top: ph_sox2_run18 × ph_sox9_run18 (score=5.996), mpnn18 × ph_sox9_run18 (5.406), mpnn17 × ph_sox9_run1 (5.406)

### Stage 5 — Caspase geometry (flipped orientation)
- Lead: **BC mpnn17 × ph_sox9_run18, GGS×1, flipped** — 8.1 Å (1/3 seeds pass)
- All top-4 pairs: SOX9 C-term points AWAY from SOX2 (d_toward = −0.76 to −0.92)
- Fix: flipped chain D = CaspC + GGS + SOX9binder (N-terminal caspase fusion)

### MD replicates (lead construct)
| Rep | Mean dist | %<15Å | Pass |
|-----|-----------|-------|------|
| rep1 | 31.0±40.4 Å | 76.3% | ✓ |
| rep2 | 25.7±39.8 Å | 80.1% | ✓ |
| rep3 | 20.9±31.6 Å | 87.8% | ✓ |
| **Mean** | **25.9 Å** | **81.4%** | **ROBUST PASS** |

### Interface bridge binder — all attempts failed
| Method | Runs | Result |
|--------|------|--------|
| BindCraft (hotspot constrained) | 233+92 traj | 0 accepted; i_pTM ~0.12–0.16 |
| Protein Hunter | 50 runs | 0 accepted; best_ipTM 0.46, pLDDT 0.35 (disordered) |
| Boltz-2 linker screen | 18 combos | dual contact ~0.20 each — too weak |

Best linker screen hit: mpnn17×sox9run1 fwd GGS (SOX2=0.21, SOX9=0.21, ipTM=0.46)
Verdict: fused single-chain doesn't co-fold reliably; inter-domain distances 22–53 Å.

---

## Lab Readiness Assessment

### ✓ READY — Individual binders (binding validation)
10 fully selective SOX2/SOX9 binder pairs validated in silico through 3 stages.
Express as His-SUMO fusions in *E. coli*; validate binding to purified HMG domains by SPR/MST.

**Top sequences to order/express:**
| Binder | Target | Length | Sequence |
|--------|--------|--------|----------|
| SOX2_mpnn17 | SOX2 | 71 AA | PYEIPEEDKEWDEYSFAMVEVHERLMEAERQLRLLSPEEEEKKRELEEEIERARKELEEVHKRGEAIERKM |
| SOX2_mpnn18 | SOX2 | 71 AA | PYEIPEEDREWDEYSFRMVEVHERLMEARRELRLLAPEETEKKRELEEEIARAEKELEEVHKEGEAIERKM |
| ph_sox2_run18 | SOX2 | 97 AA | SRLEELVDHIVDSMEEANQGGLTDEAREMIRNNLLEDLLDVNPNVEAELDAASDAEFAELLAGVDRIIREGTTSDLLGFLEALAEGGVAGALAELAK |
| ph_sox9_run18 | SOX9 | 67 AA | MEEVRELLERARELEEKVRELEDDPSKREERERLLEEAERLREEALDLLVEQIREPLLQWMREEGLL |
| ph_sox9_run1 | SOX9 | 77 AA | PPLSPEEQELEEIEEGYRNLTEKSGEALLKALGITSKNEEEREKKMRELALELAREVTRKLSEKTGVRYAEDGTRLP |
| ph_sox9_run4 | SOX9 | 86 AA | SEGIELMKELGLPVKDGKPYIEDLENLPELTNEEKLELGLLYGELNGGEKKAREAVEKAKKEGTEEAIERAFEEAVEAGDEGLKKK |

### ✓ READY — Two-chain caspase construct (in vitro validation)
Lead construct: BC mpnn17 × ph_sox9_run18, GGS×1, flipped orientation
- Chain C: SOX2_mpnn17 (71 AA) + GGS + CaspN (aa 140–315, 176 AA) = 250 AA total
- Chain D: CaspC (aa 316–416, 101 AA) + GGS + ph_sox9_run18 (67 AA) = 171 AA total

MD: 81.4% of frames have inter-half distance < 15 Å across 3 independent replicates → robust.

**In vitro validation plan:**
1. Clone both chains separately into His-SUMO-pET28 vector
2. Express in BL21(DE3), purify by Ni-NTA + size exclusion
3. Mix both chains (equimolar) + purified SOX2 HMG + SOX9 HMG in vitro
4. Measure caspase-3 activity: DEVDase assay (DEVD-AMC substrate, fluorescence)
5. AND-gate controls: chains alone, with only SOX2, with only SOX9, with both → activity only with both

### ✗ NOT READY — Interface bridge binder
0 accepted designs from all methods. RFdiffusion (200 backbones, ETA Sep 6 evening) is the most
promising approach. Results needed before any lab decision.

### ✗ NOT ADDRESSED — Cellular delivery
Proteins cannot enter nuclei without a delivery strategy:
- Options: CPP fusion (e.g. TAT, penetratin), LNP-mRNA encoding, lentiviral/AAV expression
- Nuclear localisation signal (NLS) needed for SOX2/SOX9 targeting
- This is a full campaign in itself — defer until binding validated

---

## Recommended lab timeline

| Week | Action |
|------|--------|
| Now | Order synthetic peptides: top 3 SOX2 + top 3 SOX9 binders |
| Week 1–2 | CD spectroscopy (confirm folding), MST/SPR (confirm binding to HMG) |
| Week 1–2 | Clone caspase chains (gene synthesis → pET28) |
| Week 3–4 | Express + purify caspase chains |
| Week 4–5 | In vitro DEVDase assay (AND-gate validation) |
| Week 5+ | Design delivery strategy based on binding results |
| Ongoing | Wait for RFdiffusion/N-term campaign results → update construct if better binders emerge |

---

---

## RFdiffusion bridge results (Sep 6, 200 designs complete)

| Metric | Value |
|--------|-------|
| Total designs | 200 |
| pLDDT >0.80 | 191/200 (96%) — all backbones well-folded |
| SOX2 hotspot contact <8 Å | 109/200 (54%) |
| SOX9 hotspot contact <8 Å | 36/200 (18%) |
| **Both chains contacted** | **5/200 (2.5%) strictly; 20/200 within 10 Å** |

**Top 5 dual-contact backbones (both hotspots <8 Å):**
| Design | Length | pLDDT | SOX2 dist | SOX9 dist |
|--------|--------|-------|-----------|-----------|
| bridge_147 | 85 AA | 0.833 | 5.72 Å | 6.05 Å |
| bridge_167 | 84 AA | 0.846 | 6.86 Å | 7.01 Å |
| bridge_60 | 60 AA | 0.888 | 7.09 Å | 7.11 Å |
| bridge_197 | 61 AA | 0.844 | 6.81 Å | 7.72 Å |
| bridge_45 | 97 AA | 0.819 | 7.31 Å | 7.74 Å |

**bridge_147 is the lead backbone** — both hotspots under 6.1 Å, 85 AA, pLDDT=0.833.

**ProteinMPNN + Boltz-2 pipeline queued (PID 24658, waits for BindCraft PID 12669):**
- Script: `rfdiffusion/run_mpnn_boltz_bridge.py`
- Top 20 backbones × 10 MPNN seqs = 200 Boltz-2 predictions
- Filter: ipTM ≥ 0.40, SOX2 contact ≥ 0.12, SOX9 contact ≥ 0.12
- Results: `rfdiffusion/bridge_mpnn_boltz/mpnn_boltz_results.tsv`

---

## ProteinMPNN + Boltz-2 bridge results (Sep 7, COMPLETE)

220 sequences tested (top 20 backbones × 11 MPNN seqs). **5/220 passed** (ipTM ≥ 0.40, SOX2_contact ≥ 0.12, SOX9_contact ≥ 0.12).

| Design | ipTM | SOX2_c | SOX9_c | Len | d_toward | Orientation | Sequence |
|--------|------|--------|--------|-----|----------|-------------|----------|
| bridge_145_s06 | **0.572** | 0.224 | 0.153 | 85 | −13.04 | N-term→SOX2 | AAAAAAAAAAAAAAAAAAALAAALAARAAALAAARAREERVRAEAAAEIARRRAAIEAADPERAPALAAAAEQEVEAAALQRLLG |
| bridge_90_s01 | 0.474 | 0.222 | **0.178** | 90 | −4.47 | ambiguous | GSRLDAEIARARARAEAREAEIAALRAAAAAEAAALEAKAAEDAEHGEELLAEAAAVRERAELEAEALEQEVAAERALADARAAAAAAAR |
| bridge_55_s06 | 0.453 | 0.197 | 0.167 | 66 | −25.82 | N-term→SOX2 | AEAAAREAAAAAAAALAAEAARRAAQRAAEEAAAAAAEAAEAAARAAAEAAAKAREEEEKRLRELE |
| bridge_85_s07 | 0.428 | **0.305** | 0.122 | 82 | +6.07 | C-term→SOX2 | TLVLTKSILPSPEILPELPKRVEQVKKKAEEELARLQQQWQQVLAWIQAQDPSNQQAPQYWQQHLQQRYQTKSSEAHPVIKA |
| **bridge_123_s02** | 0.416 | 0.169 | 0.182 | 77 | −25.74 | **N-term→SOX2** | ELAALAARQEALAAQAEAARAELAALLAAGAPQAEIDAAQARLDELEAARDAAAEALARLEAERAAAAAAAAAAAAA |

**Terminus direction summary (all 5 designs):**
- 4/5: N-term→SOX2 / C-term→SOX9 (d_toward < 0)
- 1/5: C-term→SOX2 / N-term→SOX9 (bridge_85_s07)

**Lead for caspase fusion: bridge_123_s02**
Cleanest geometry: N-term only 19.3 Å from SOX2 centroid, C-term only 16.7 Å from SOX9 centroid.
Proposed single-chain construct: **CaspN + GGS + bridge_123_s02 + GGS + CaspC**
→ When bridge contacts both SOX2+SOX9, CaspN (near SOX2) and CaspC (near SOX9) reconstitute.

**Bugs fixed in MPNN pipeline (Sep 7):**
1. `--chain_list=A B` → `--chain_list=C` (was designing targets, not binder)
2. Sequence parser: `if "/" in seq_full` → now handles single-chain MPNN output

## Bridge binder MD results (Sep 8, COMPLETE)

| Design | N-to-C mean | % <15 Å | SOX2_c | SOX9_c | Verdict |
|--------|------------|---------|--------|--------|---------|
| bridge_123_s02 | 34.1±15.3 Å | 0.0% | 0.148 | 0.171 | FAIL |
| **bridge_90_s01** | 50.2±52.7 Å | **55.5%** | **0.234** | 0.162 | **PASS** |

**bridge_90_s01 is the lead for single-chain caspase construct.**
The large std (52.7 Å) reflects bimodal dynamics: bridge alternates between compact (<5 Å min) and extended states, spending 55.5% of time <15 Å. Both targets remain engaged (SOX2=0.234, SOX9=0.162).

**Lead single-chain construct:**
`CaspN (aa140-315) + GGS + bridge_90_s01 (90 AA) + GGS + CaspC (aa316-416)`
Sequence: GSRLDAEIARARARAEAREAEIAALRAAAAAEAAALEAKAAEDAEHGEELLAEAAAVRERAELEAEALEQEVAAERALADARAAAAAAAR

MD replicates (rep2 + rep3) running to confirm (ETA Sep 8 evening).

## Bridge binder MD (Sep 7, RUNNING — PID 250072)

Running MD on 2 top bridge designs to validate caspase reconstitution geometry.

**Setup:** 3-chain complex (A=SOX2, B=SOX9, C=bridge_binder), 20 ns production, AMBER99SB-ILDN/TIP3P.
**Key metric:** bridge N-term Cα ↔ C-term Cα distance (where CaspN and CaspC would attach).
Pass = >50% frames <15 Å. Promising = >50% frames <30 Å (bridgeable with GGS×5 linkers).

**Designs running:**
1. **bridge_123_s02** (77 AA, N-to-C=29.7 Å static) — proper traversal: N-term→SOX2, C-term→SOX9
2. **bridge_90_s01** (90 AA, N-to-C=9.5 Å static) — both termini near SOX9 (compact fold)

Script: `md/bridge_binder_md/run_bridge_md.py`
Results: `md/bridge_binder_md/<design>/metrics.json` — ETA Sep 8 morning.

**N-to-C distances (static Boltz-2 structures):**
| Design | ipTM | N-to-C (static) | N-term near | C-term near |
|--------|------|-----------------|------------|------------|
| bridge_123_s02 | 0.416 | 29.7 Å | SOX2 (19 Å) | SOX9 (17 Å) |
| bridge_90_s01 | 0.474 | **9.5 Å** | SOX9 (23 Å) | SOX9 (21 Å) |
| bridge_145_s06 | 0.572 | 54.3 Å | ambiguous | ambiguous |
| bridge_85_s07 | 0.428 | 12.2 Å | SOX2 (14 Å) | SOX2 (17 Å) |
| bridge_55_s06 | 0.453 | 95.9 Å | SOX2 (18 Å) | dangling |

## Ongoing computational campaigns (as of Sep 7 23:10 BST)
- Ph1: RFdiffusion bridge both chains — **DONE**
- Ph3: BindCraft no-hotspot — **KILLED** (Sep 7; 2 SOX9-only designs, not bridges)
- ProteinMPNN+Boltz-2 bridge — **DONE** (5 accepted, Sep 7 19:27)
- Bridge binder MD — **RUNNING** (PID 250072, ETA Sep 8 AM)
- Ph4: RFdiffusion SOX9 N-term — **still queued** (PID 7084 master)
- Ph A/B/C: PH + BoltzGen + BindCraft SOX9 N-term — queued after Ph4 (PID 9949)
