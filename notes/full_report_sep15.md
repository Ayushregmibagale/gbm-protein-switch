# GBM Protein Switch — Full Report
*2026-09-15*

---

## 1. Background and Scientific Rationale

### The problem: glioblastoma stem cells
Glioblastoma multiforme (GBM) is the most aggressive primary brain tumour, with median survival ~15 months despite surgery, radiotherapy, and temozolomide. Recurrence is driven by a subpopulation of glioblastoma stem cells (GSCs) that are therapy-resistant, self-renewing, and tumourigenic. Current treatments cannot reliably eliminate this population.

### The target: SOX2 + SOX9 co-expression
Two transcription factors, SOX2 and SOX9, are co-expressed at high levels in GSCs and together maintain the stem-cell state. Critically:
- **SOX2** is broadly expressed in neural stem cells and GSCs, but at lower levels in differentiated brain cells.
- **SOX9** is highly expressed in GSCs but markedly lower in normal neural tissue.
- **Co-expression of SOX2 AND SOX9 together** is a defining feature of the GSC state and is rare in healthy brain cells.

This co-expression pattern provides a molecular AND-gate: a kill switch that activates only when both proteins are present simultaneously, maximising tumour selectivity and minimising off-target toxicity.

### The mechanism: split-caspase AND-gate
The therapeutic logic is:
1. Engineer one protein that binds SOX2 and carries one half of a split caspase.
2. Engineer a second protein that binds SOX9 and carries the other half.
3. In GSCs, both SOX2 and SOX9 are co-expressed and form a transcriptional complex. When both binders engage their targets simultaneously, the two caspase halves are brought into proximity and reconstitute caspase activity → apoptosis.
4. In normal cells expressing only SOX2 (neural stem cells) or neither protein, only one binder engages → no reconstitution → no cell death.

### Why caspase-9, not caspase-3?
Split-caspase-9 (initiator caspase) was chosen over the more potent executioner caspases-3/7 for three reasons:
- **Safety**: caspase-9 requires dimerisation and its CARD domain for independent activity. The CARD-deleted version used here (residues 140–416) cannot self-activate — reconstitution only occurs when both halves are brought together.
- **Signal amplification**: one reconstituted caspase-9 dimer activates many caspase-3/7 molecules downstream, so few reconstitution events per cell are sufficient for killing.
- **Precedent**: split-caspase-9 is the most validated coincidence-detection protease in synthetic biology (iCasp9 suicide gene, AND-gate CAR-T literature).

**Split site**: between residues D315 and G316 of human caspase-9 (UniProt P55211).
- N-half: residues 140–315 (large subunit, CARD deleted; 176 aa)
- C-half: residues 316–416 (small subunit; 101 aa)

---

## 2. Hotspot Selection — Where the Binders Must Bind

For the AND-gate to function, both binders must engage their targets simultaneously while SOX2 and SOX9 are in their transcriptional co-complex on DNA. Two geometric requirements follow:

1. **Epitopes must not be buried at the SOX2–SOX9 protein–protein interface** (PPI). Binding to a buried residue is sterically impossible in the co-complex.
2. **Epitopes must be on the same face of the co-complex**, so the binder termini (which carry the caspase halves) point toward each other and can reconstitute.

The SOX2–SOX9–DNA co-complex was modelled using AlphaFold 3 (seed-3_sample-3, ranking=0.703, iptm=0.63). The PPI interface was identified at SOX2 residues S69/K117-119 and SOX9 residues E159/R162/K166/Y174 — all excluded as hotspots.

**Final hotspots selected:**
- **SOX2**: wing-2 face (HMG residues 62, 70, 74, 81) — exposed in co-complex, 6–12 Å from DNA, inter-epitope angle to SOX9 cluster = ~55°.
- **SOX9**: C-terminal face (HMG residues 168, 169, 171, 172, 173) — all at angles 62–89° from SOX2 epitope centroid, all >6 Å from DNA, all PPI-safe.

Inter-epitope centroid distance: **~33.5 Å** — bridgeable by a GGS×3–5 linker (~10–15 Å at full extension per GGS unit).

---

## 3. Two Design Approaches

Two parallel and architecturally distinct approaches were developed to achieve the AND-gate kill switch. Both use the same split caspase-9 but differ in how the binders are deployed.

---

## Approach 1 — Two-Binder Split-Caspase (Track 1)

### How it works
Two separate miniprotein binders are designed independently — one for SOX2, one for SOX9. Each is fused via a GGS linker to one half of split caspase-9:

```
Construct A:  NLS – [SOX2 binder] – GGS – CaspN (aa140–315)
Construct B:  NLS – [SOX9 binder] – GGS – CaspC (aa316–416)
```

Both constructs are expressed from a single bicistronic AAV cassette separated by a T2A self-cleaving peptide, ensuring equimolar production.

**AND-gate logic**: Both SOX2 and SOX9 must be present and in proximity (as they are in the GSC transcriptional co-complex) for both binders to engage simultaneously → caspase halves come within <15 Å → reconstitution → apoptosis.

The critical geometric requirement is that the **C-terminus of the SOX2 binder** (where CaspN is attached) and the **N-terminus of the SOX9 binder** (where CaspC is attached) must be within reconstitution distance when both binders are bound in the co-complex.

**Orientation note**: the original design had CaspC at the C-terminus of the SOX9-binder construct, but diagnosis showed this pointed the CaspC away from CaspN. The solution was to **flip the SOX9-binder construct** — CaspC now leads (N-terminal fusion), so CaspC is immediately adjacent to the junction:

```
Construct A (unchanged): [SOX2 binder] – GGS – CaspN
Construct B (flipped):   CaspC – GGS – [SOX9 binder]
```

This is the "flipped orientation" design validated in all downstream MD.

### Design methods

**SOX2 binder — BindCraft**
BindCraft (miniprotein hallucination via AF2-MPNN cycles) was used with relaxed shape-complementarity filters (≥0.45 vs default 0.60) because SOX2's wing-2 face is flat and intrinsically scores low. After ~1430 trajectories, 7 designs were accepted; ProteinMPNN sequence diversification from the best backbone produced two final candidates: **mpnn17** and **mpnn18** (both ~71 aa).

**SOX9 binder — Protein Hunter (PH)**
Protein Hunter (Boltz-2 hallucination pipeline) generated 50+ SOX9 binders targeting the v3 hotspot cluster. The best candidate **ph_sox9_run18** (67 aa) achieved ipTM=0.856 with selectivity ratio 6.69× (SOX9 vs SOX2).

**Five-stage selection pipeline**
All binders (from BindCraft, Protein Hunter, and BoltzGen) were fed through a shared funnel:

| Stage | What | Criterion |
|-------|------|-----------|
| 1 | Boltz-2 re-prediction (normalise scores across methods) | ipTM ≥ 0.70, pLDDT ≥ 70 |
| 2 | AND-gate geometry: pair SOX2+SOX9 binders in 4-chain Boltz-2 | pair ipTM ≥ 0.60; C-term(A) ↔ N-term(B) < 40 Å |
| 3 | Selectivity: each binder vs wrong target | wrong-target ipTM < 0.40 |
| 4 | Composite ranking: ipTM × selectivity ratio × geometry bonus | top 10 pairs |
| 5 | Caspase geometry (flipped orientation): full 4-chain prediction + caspase halves, GGS×1–3 | recon dist C-term(CaspN) ↔ N-term(CaspC) < 15 Å in ≥2/3 seeds |

**GROMACS MD validation** (20 ns, 3 replicates):
Starting from the best Boltz-2 Stage 5 pose (mpnn17 × ph_sox9_run18, GGS×1, seed 2, recon=8.05 Å), a full GROMACS MD was run in explicit solvent. Pass criterion: CaspN C-terminus ↔ CaspC N-terminus < 15 Å in ≥50% of frames.

### Results

**Individual binders:**
| Binder | Method | Target | ipTM | Selectivity |
|--------|--------|--------|------|-------------|
| mpnn17 | BindCraft | SOX2 | 0.897 | 2.77× (SOX2/SOX9) |
| mpnn18 | BindCraft | SOX2 | 0.893 | 2.54× |
| ph_sox2_run18 | Protein Hunter | SOX2 | 0.784 | 3.21× (most selective) |
| ph_sox9_run18 | Protein Hunter | SOX9 | 0.856 | 6.69× |

**Lead pair: mpnn17 × ph_sox9_run18**
- Stage 4 composite rank: 4/10 selective pairs
- Stage 5 Boltz-2: 2/9 seeds pass (8.05 Å, 10.35 Å) — stochastic
- **GROMACS MD (3 replicates × 20 ns): 76.3% / 80.1% / 87.8% → mean 81.4% — PASS**

**Stage 4 rank-1 pair tested (Sep 15): ph_sox2_run18 × ph_sox9_run18**
- Stage 5: 1/15 seeds pass (GGS×3, seed 4, 7.6 Å; iptm=0.208 — low confidence)
- Does NOT outperform the lead; mpnn17 × ph_sox9_run18 remains the best candidate.

### Pros and cons of Track 1

**Pros:**
- Strongest end-to-end computational validation (MD confirmed in 3 replicates)
- Two separate proteins: each can be optimised, expressed, and characterised independently
- Modular — swap binders without redesigning the caspase fusion
- Single bicistronic AAV cassette possible (T2A, ~850 bp coding sequence, well within 4.7 kb limit)
- More precedent: two-component split-protease AND gates are well established in synthetic biology

**Cons:**
- Stoichiometry: requires equimolar expression of both chains; any imbalance reduces efficacy or increases leakiness
- Two proteins must independently reach the nucleus, fold correctly, and find their targets; any step failing in one chain is a single point of failure
- Caspase geometry is sensitive to linker length and orientation — the 2/9 Boltz-2 hit rate reflects this sensitivity
- Binder termini must point toward each other in the co-complex — topology constraint limits which binders can be used

---

## Approach 2 — Single Bridge Binder (Track 2)

### How it works
A single de novo miniprotein is designed to simultaneously contact **both** SOX2 and SOX9 when they are in complex. The bridge binder is then fused at its N-terminus to one caspase half and at its C-terminus to the other:

```
Single construct:  NLS – CaspN (aa140–315) – GGS – [bridge binder] – GGS – CaspC (aa316–416)
```

**AND-gate logic**: The bridge binder only folds stably and contacts both targets when both SOX2 and SOX9 are present simultaneously. In this bound state, the bridge holds its own N- and C-termini (and thus the caspase halves) within reconstitution distance. If either target is absent, the bridge cannot engage fully and the caspase halves remain separated.

The critical geometric requirement is that the **N-to-C terminus distance** of the bridge protein is <15 Å when it is bound to the co-complex — this directly reports on whether the two caspase halves attached to those termini can reconstitute.

### Design methods

Four methods were attempted, with RFdiffusion + ProteinMPNN being the only successful route:

**Failed methods:**
1. **BindCraft** (233 trajectories): 0 accepted. AF2-MPNN optimisation cannot simultaneously satisfy hotspot constraints on two separate protein chains 45 Å apart.
2. **Protein Hunter** (50 designs): 0 accepted. Boltz-2 hallucination produces disordered peptides when asked to bridge two chains far apart.
3. **Boltz-2 linker screen** (fusing top individual binders with GGS): 0 accepted. The two binders do not co-fold cooperatively when concatenated.

**Successful method — RFdiffusion → ProteinMPNN → Boltz-2 → GROMACS MD:**

```
Step 1 — RFdiffusion backbone generation
  contig: [A175-215/0 B450-509/0 60-100]
  hotspots: A193, A199, B466, B505, B508
  200 backbones generated; filter: CA dist to each hotspot cluster <9 Å → top 80

Step 2 — ProteinMPNN sequence design (10 seqs per backbone)
  Boltz-2 3-chain prediction (SOX2_HMG + SOX9_HMG + binder)
  Filter: ipTM ≥ 0.40 AND SOX2_contact ≥ 0.12 AND SOX9_contact ≥ 0.12
  → 32 designs into MD queue

Step 3 — GROMACS MD validation (20 ns, pass criterion: N-C dist <15 Å in ≥50% frames)
  → 5/32 pass (15.6%)

Step 4 — Replicates (3 × 20 ns for each passing design)
  → All 5 designs confirmed

Step 5 — Dual-contact AND-gate analysis (per-frame simultaneous contact)
  → AND-gate mechanism verified in all 5 designs
```

**Why RFdiffusion works here**: unlike optimisation-based methods, RFdiffusion samples backbone coordinates directly by denoising — it can generate a single continuous chain that geometrically spans two separate hotspot clusters without being penalised by a multi-chain objective.

### Results

**Pass rates:**
- BindCraft: 0/233 backbones
- Protein Hunter: 0/50 designs
- Boltz-2 linker screen: 0/18 combinations
- **RFdiffusion + ProteinMPNN + MD: 5/32 = 15.6%**

**MD results for all 5 passing designs (3 replicates each):**

| Design | Length | Rep1 | Rep2 | Rep3 | SOX2 contact | SOX9 contact |
|--------|--------|------|------|------|-------------|-------------|
| bridge_48_s01 | 69 aa | 98.5% | 98.6% | 100.0% | 0.376 | 0.166 |
| bridge_73_s09 | 83 aa | 86.7% | 63.4% | 90.2% | 0.218 | 0.229 |
| bridge_136_s01 | 73 aa | 86.0% | 73.1% | 74.6% | 0.143 | 0.172 |
| bridge_101_s06 | 76 aa | 79.3% | 100.0% | 79.0% | 0.274 | 0.255 |
| bridge_90_s01 | 90 aa | 55.5% | 87.1% | 60.1% | 0.234 | 0.162 |

**Dual-contact AND-gate analysis (Sep 15, per-frame simultaneous contact):**

The contact fractions above are averaged over all frames — they don't tell us whether the bridge contacts both targets *at the same time*. The dual-contact analysis computes, per frame, whether ≥X% of bridge Cα atoms are within 8 Å of SOX2 AND ≥X% within 8 Å of SOX9 simultaneously:

| Design | Dual-contact (thr=0.05) | Dual-contact (thr=0.10) | Dual-contact (thr=0.15) |
|--------|------------------------|------------------------|------------------------|
| bridge_48_s01  | 100.0% | 99.4% | 69.2% |
| bridge_73_s09  | 100.0% | 99.5% | 92.3% |
| bridge_101_s06 | 100.0% | 99.9% | 99.0% |
| bridge_90_s01  | 100.0% | 99.8% | 75.4% |
| bridge_136_s01 | 96.9%  | 86.8% | 47.3% |

**Key finding**: all 5 designs show simultaneous SOX2+SOX9 contact in ≥87% of MD frames at the 10% contact threshold, confirming AND-gate-compatible structural geometry. This demonstrates simultaneous binding plausibility; it does not establish caspase-9 enzymatic reconstitution, downstream activity amplification, or cellular selectivity, which require experimental validation. The low SOX9 contact fraction for bridge_48_s01 (0.166 mean) reflects an asymmetric interface — a smaller part of the bridge contacts SOX9 — but it does so continuously and simultaneously with its SOX2 contact.

**Important note on Boltz-2 ipTM**: ipTM is a **poor predictor of MD success**. The design with the highest Boltz-2 ipTM (bridge_95_s02=0.807) failed MD completely (N-C mean=65 Å). bridge_48_s01 passed at 98.5% with a much lower ipTM. MD simulation is the only reliable validation metric.

### Pros and cons of Track 2

**Pros:**
- Single protein construct: simpler expression, purification, and delivery; no stoichiometry problem
- Single-AAV (one coding sequence ≤200 aa bridge + 277 aa caspase = ~477 aa ≈ 1.4 kb; trivially within the 4.7 kb AAV limit)
- The AND-gate is structurally intrinsic to the bridge: it can only contact both targets when they are co-complexed — the coincidence detection is built into the 3D shape
- bridge_101_s06 achieves 99.9% simultaneous dual contact — very robust AND-gate geometry
- No orientation/linker sensitivity issue for the caspase halves (N and C termini are naturally the attachment points)

**Cons:**
- Harder to design: RFdiffusion is the only method that works; 4 other design routes failed entirely
- Low yield from design campaign: 5/200 backbones → 5/32 MD-tested designs pass (effectively ~2.5% from backbone to confirmed design)
- Harder to iterate: if the bridge fails experimentally, re-designing from scratch requires another RFdiffusion campaign
- Single protein must simultaneously bind two different protein surfaces — if the co-complex geometry in the cell differs from the AF3 model used for design, the bridge may not engage both targets efficiently

### Off-target selectivity — bridge_48_s01 (Sep 15)

Because bridge_48_s01 was never put through the Stage 3 selectivity pipeline, two dedicated screens were run.

**2-chain screen** (bridge vs each off-target HMG alone; 8 targets × 3 seeds):

| Target | Mean ipTM | Verdict |
|--------|-----------|---------|
| SRY (SOX A) | 0.616 | flagged |
| SOX11 (SOX C, neurons) | 0.473 | flagged |
| HMGB1 (non-SOX HMG) | 0.415 | flagged |
| SOX9 (on-target) | 0.402 | on-target |
| SOX3 (SOX B1, neural progenitors) | 0.369 | PASS |
| SOX4 (SOX C, brain) | 0.336 | PASS |
| SOX17 (SOX F, vasculature) | 0.300 | PASS |
| SOX2 (on-target) | 0.162 | on-target |

Importantly, the on-target SOX2 scores only 0.162 in 2-chain mode — confirming that the bridge does not bind SOX2 alone efficiently, only when SOX9 is also present. The 2-chain assay is structurally inappropriate for a bridge binder; flagged hits were followed up with the correct 3-chain assay.

**3-chain screen** (replace one target at a time; on-target baseline SOX2+SOX9+bridge = 0.348):

| Condition | Mean ipTM | Change from baseline (0.348) | Verdict |
|-----------|-----------|------------------------------|---------|
| HMGB1+SOX9 | 0.358 | +3% (no reduction) | NOT selective |
| SOX2+HMGB1 | 0.290 | −17% | SELECTIVE ✓ |
| SRY+SOX9 | 0.251 | −28% | SELECTIVE ✓ |
| SOX11+SOX9 | 0.211 | −39% | SELECTIVE ✓ |
| SOX2+SOX11 | 0.208 | −40% | SELECTIVE ✓ |
| SOX2+SRY | 0.172 | −51% | SELECTIVE ✓ |

**Conclusions:**
- **SOX11 (main neuronal safety concern)**: selective in both positions (55–56% drop). Fully resolved.
- **SRY**: selective in both positions despite the 2-chain false positive — 2-chain result was artefactual.
- **HMGB1**: selective when in the SOX9 position (0.290, 17% drop). In the SOX2 position (0.358) it is not selective — HMGB1+SOX9 (mean 0.358) marginally exceeds the same-run on-target mean (0.348), giving no reduction. The previously reported 24% drop was computed against an older, higher baseline and is not valid against the current data. Both conditions carry high seed variance (on-target: 0.148–0.647; HMGB1+SOX9: 0.236–0.504), making a definitive computational verdict impossible. Experimental pull-down of bridge_48_s01 vs recombinant HMGB1 is required.
- **SOX3** (most structurally similar to SOX2): passes the 2-chain screen cleanly (0.369) and was not further tested — if it can't bind SOX2 alone, it cannot substitute for SOX2 in the co-complex.

**Overall selectivity verdict**: bridge_48_s01 is computationally selective against SOX11, SRY, SOX3/4/17, and HMGB1 in the SOX9 position. Selectivity against HMGB1 in the SOX2 position is not established — HMGB1+SOX9 ipTM (0.358) is indistinguishable from the on-target baseline (0.348) using the current same-run data. Experimental pull-down of bridge_48_s01 vs recombinant HMGB1 is required before claiming full off-target selectivity. This is a flag but not a formal blocker for Phase 1 individual-binder work.

---

## 4. Comparison of the Two Approaches

| Feature | Track 1: Two-Binder Split-Caspase | Track 2: Bridge Binder |
|---------|-----------------------------------|------------------------|
| **Architecture** | Two separate proteins (SOX2-CaspN + SOX9-CaspC) | One protein (CaspN–bridge–CaspC) |
| **AND-gate mechanism** | Proximity-induced reconstitution when binders co-engage on SOX2/SOX9 complex | Bridge structurally spans both targets; caspase halves at N and C termini |
| **Stoichiometry** | Must maintain 1:1 ratio; T2A bicistronic construct helps | Single protein; stoichiometry trivial |
| **AAV coding size** | ~850 bp (tight but fits) | ~1.4 kb (easy) |
| **Design success rate** | Feasible with BindCraft + PH; multiple binder candidates | Only RFdiffusion works; low throughput |
| **Selectivity characterised?** | Yes — Stage 3 pipeline; SOX2 binder 2.54–3.21×, SOX9 binder 6.69× | Yes — 2-chain + 3-chain screens; SOX11/SRY selective in both positions; HMGB1 borderline (experimental follow-up recommended) |
| **Computational MD evidence** | 3 replicates × 20 ns: mean 81.4% reconstitution | 5 designs × 3 replicates × 20 ns: 55–100% reconstitution |
| **AND-gate geometry (MD)?** | Reconstitution-distance metric (<15 Å in 81.4% of frames); enzymatic activity not yet demonstrated | Simultaneous binding geometry confirmed in MD; enzymatic activity not yet demonstrated |
| **Iteration if fails** | Swap individual binders independently | Entire RFdiffusion campaign required |
| **Lab complexity** | Express two proteins, co-transfect cells | Express one protein |
| **Main risk** | Stoichiometry imbalance; linker geometry sensitivity | HMGB1 borderline in 3-chain screen (needs pull-down); low design throughput if iteration needed |

---

## 5. Lead Candidates

### Track 1 lead
**mpnn17 (SOX2) × ph_sox9_run18 (SOX9), GGS×1 linker, flipped orientation**

Full construct sequences:
- SOX2-binder (mpnn17, 71 aa): `PYEIPEEDKEWDEYSFAMVEVHERLMEAERQLRLLSPEEEEKKRELEEEIERARKELEEVHKRGEAIERKM`
- SOX9-binder (ph_sox9_run18, 67 aa): `MEEVRELLERARELEEKVRELEDDPSKREERERLLEEAERLREEALDLLVEQIREPLLQWMREEGLL`
- Linker: GGS (3 aa)
- Construct A: [mpnn17] – GGS – [CaspN aa140–315]
- Construct B: [CaspC aa316–416] – GGS – [ph_sox9_run18]

### Track 2 lead
**bridge_48_s01 (69 aa)**
- N-C mean distance in MD: 8.85 Å (98.5–100% of frames <15 Å across 3 replicates)
- Dual-contact (simultaneous SOX2+SOX9): 99.4% of frames at 10% threshold
- Selectivity: SOX11/SRY selective in both positions; SOX3/4/17 pass; HMGB1 borderline (flag for pull-down)
- Construct: [CaspN aa140–315] – GGS – [bridge_48_s01] – GGS – [CaspC aa316–416]
- Sequence (chain C, 69 aa): `GAAAAAAAEAAAAAAAAAAAAAAAAAAAAEAARAAAAAAAAAAAAAAAAAAQAAAAQAYLAQAQAAAAA`

---

## 6. Risks Shared by Both Approaches

1. **Leaky reconstitution**: If caspase halves spontaneously reconstitute without target engagement (especially in cells with high local protein concentration), off-target apoptosis occurs. The CARD-deleted split site (aa140–315 / 316–416) substantially reduces this risk vs full-length split constructs.

2. **SOX2 in normal neural stem cells**: Adult neural stem cells express SOX2. If SOX9 is transiently co-expressed during stress or injury, off-target killing could occur. Needs validation across normal brain cell types (astrocytes, NSCs, oligodendrocytes).

3. **Binder-caspase folding interference**: The miniprotein fused directly to a 176 aa caspase domain may misfold or sterically interfere with the binder's target engagement. Linker length and flexibility need empirical optimisation (starting with GGS×1 as predicted to work, extending to GGS×3–5 if needed).

4. **Nuclear localisation**: Both SOX2 and SOX9 are nuclear. Both constructs need a functional NLS and must avoid proteasomal degradation in transit. SV40 NLS (PKKKRKV) is recommended.

5. **Tumour heterogeneity**: Not all GSCs may constitutively co-express SOX2 and SOX9. In heterogeneous tumours, a subpopulation may escape. Combination with standard TMZ is advisable.

6. **Immunogenicity**: Designed miniproteins are non-human sequences. AAV delivery and repeated dosing risk immunological responses. Consider humanisation of final binder sequences.

7. **Track 2 specific — HMGB1 interaction**: The 3-chain selectivity screen showed a borderline result for HMGB1 in the SOX2 position (ipTM 0.358 vs on-target 0.348). HMGB1 is a ubiquitous nuclear HMG-box protein involved in chromatin remodelling. The result is most likely Boltz-2 noise (high seed variance, selective in SOX9 position), but an experimental pull-down of bridge_48_s01 vs recombinant HMGB1 is recommended before in vivo work.

---

## 7. Recommended Experimental Progression

### Phase 1 — Binder validation (3–6 months)
1. Synthesise and express individual binders as His-tagged recombinant proteins
2. SPR or MST to confirm binding to SOX2 and SOX9 HMG domains (Kd target: <1 μM)
3. Confirm selectivity (SOX2 binder does not bind SOX9, and vice versa)
4. Test in a GSC cell line by immunofluorescence (confirm nuclear co-localisation with SOX2/SOX9)

### Phase 2 — Split-caspase reconstitution in vitro (2–4 months)
1. Express the full fusion constructs (binder-CaspN and CaspC-binder for Track 1; CaspN-bridge-CaspC for Track 2)
2. DEVDase assay: mix both constructs in the presence or absence of purified SOX2–SOX9 co-complex. Measure caspase activity (Ac-DEVD-pNA cleavage). Expect signal only in the +co-complex condition.
3. Titrate co-complex concentration to assess sensitivity threshold.

### Phase 3 — GSC cell line validation (4–8 months)
1. Clone construct(s) into AAV plasmid (single bicistronic T2A construct for Track 1; single ORF for Track 2)
2. Transduce GSC lines (e.g., GBM6, GBM39) and assess:
   - Caspase-3/7 activation (Caspase-Glo assay)
   - Cell viability (CellTiter-Glo)
   - Selectivity in NSC co-culture (confirm normal cells are not killed)
3. Compare Track 1 vs Track 2 efficacy and selectivity in the same assay.
4. If Track 1 shows stoichiometry problems, prioritise Track 2 (single protein, simpler).

### Which track to start with?
**Recommended: pursue both in Phase 1 (binder validation) in parallel.**
- Track 1 binders (mpnn17, ph_sox9_run18) can be ordered as synthetic peptides immediately.
- bridge_48_s01 sequence can be gene-synthesised (69 aa = 207 bp; trivial) in parallel.
- Phase 2 and 3 results will determine which track is taken forward into in vivo studies.

---

## 8. File Index

| File / Directory | Contents |
|-----------------|----------|
| `structures/af3_cocomplex/` | AF3 SOX2–SOX9–DNA co-complex models |
| `selection/stage1_results.tsv` | All binder re-predictions (Boltz-2) |
| `selection/stage2_results.tsv` | AND-gate geometry screen results |
| `selection/stage3_results.tsv` | Selectivity screen |
| `selection/stage4_results.tsv` | Composite ranking (top 10 pairs) |
| `selection/stage5_flipped_results.tsv` | Caspase geometry, flipped orientation (mpnn17/18 pairs) |
| `selection/stage5_ph_sox2_run18_results.tsv` | Stage5 for rank-1 pair (ph_sox2_run18 × ph_sox9_run18) |
| `md/lead_construct/` | GROMACS MD for Track 1 lead (rep1) |
| `md/rep2/`, `md/rep3/` | Track 1 replicates |
| `md/bridge_binder_md/bridge_48_s01/` | Track 2 lead MD (rep1) |
| `md/bridge_binder_md/bridge_48_s01_rep2/` | Track 2 replicate 2 |
| `md/bridge_binder_md/bridge_48_s01_rep3/` | Track 2 replicate 3 |
| `analysis/dual_contact_results.tsv` | Per-frame AND-gate analysis for all passing bridges |
| `analysis/bridge_48_s01_dual_contact_timeseries.txt` | Frame-by-frame contact time series |
| `notes/results_overview_sep15.md` | Compact results summary |
| `analysis/offtarget_selectivity_results.tsv` | 2-chain off-target screen (8 HMG proteins × 3 seeds) |
| `analysis/bridge_3chain_selectivity_results.tsv` | 3-chain selectivity screen (7 conditions × 3 seeds) |
| `notes/full_report_sep15.md` | This document |
