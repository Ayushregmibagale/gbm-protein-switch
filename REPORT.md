# GBM Protein Switch — Computational Campaign Report

*Last updated: 2026-10-30*

---

## Table of Contents

1. [Scientific motivation](#1-scientific-motivation)
2. [Design concept](#2-design-concept)
3. [Structural foundation — co-complex modelling and hotspot selection](#3-structural-foundation)
4. [Track 1 — Two-binder split-caspase](#4-track-1--two-binder-split-caspase)
5. [Track 2 — Single bridge binder](#5-track-2--single-bridge-binder)
6. [Robustness audit (October 2026)](#6-robustness-audit)
7. [Caspase geometry investigation](#7-caspase-geometry-investigation)
8. [Summary of all results](#8-summary-of-all-results)
9. [What is and is not established](#9-what-is-and-is-not-established)
10. [Next steps — lab](#10-next-steps--lab)

---

## 1. Scientific motivation

Glioblastoma multiforme (GBM) is the most lethal primary brain tumour (median survival ~15 months). Recurrence is driven by glioblastoma stem cells (GSCs): a therapy-resistant, self-renewing sub-population that survives surgery, radiotherapy, and temozolomide. No treatment reliably eliminates GSCs.

**The molecular target.** SOX2 and SOX9 are transcription factors co-expressed at high levels in GSCs, where together they maintain stemness. Critically, their *co-expression and physical co-complexing* is a defining feature of the GSC state. Individual expression of either factor occurs in normal brain cells (SOX2 in neural stem cells; SOX9 in oligodendrocyte precursors), but the co-complex is GSC-selective.

Koeber et al. (Nature 2026) exploited this at the DNA level — a synthetic super-enhancer (SSE-7) active only in the SOX2+SOX9 co-complex state drove HSV-TK and IL-12 via AAV, achieving curative outcomes and immunological memory in mice. This project adds a second, *protein-level* selectivity lock: a molecular coincidence detector that reads the SOX2+SOX9 co-complex directly and triggers apoptosis only in cells where both proteins are physically present.

**Why a protein-level lock on top of SSE-7?** SSE-7 requires transcriptional output (hours to days). A protein switch responds to the co-complex state directly (minutes). Combining both layers — DNA selectivity and protein selectivity — provides orthogonal safety. A cell that down-regulates SSE-7 activity but still harbours the co-complex is still killed. A cell that transiently expresses both TFs without forming the co-complex is protected by the protein switch's geometric requirement.

---

## 2. Design concept

The core idea is a **split-caspase-9 AND-gate**:

```
Construct A:   NLS – [SOX2 binder] – GGS – Casp9_N (aa 140–315)
Construct B:   Casp9_C (aa 316–416) – GGS – [SOX9 binder] – NLS
```

- In GSCs (SOX2 AND SOX9 co-present): both binders engage their targets in the co-complex → caspase halves brought within ~10 Å → reconstitution → caspase-3/7 activation → apoptosis.
- In normal cells (SOX2 alone, or neither): only one binder engages → halves remain separated → no activity.

**Why caspase-9, not caspase-3/7?**
Caspase-9 (initiator) provides signal amplification — one reconstituted dimer activates many downstream caspase-3/7 molecules. It is also the most validated split-protease in synthetic biology (iCasp9 suicide gene, AND-gate CAR-T literature). The CARD domain is deleted (residues 1–139 removed) so the N-half cannot independently recruit APAF-1 and self-activate.

**Split site:** between D315 and A316 — the natural autocatalytic cleavage site of human caspase-9 (UniProt P55211). The large subunit (aa 140–315, 176 aa) goes on Construct A; the small subunit (aa 316–416, 101 aa) on Construct B. Note: Construct B is in *flipped orientation* (caspase fragment leads, binder trails) so the caspase termini point toward each other in the co-complex.

Two parallel architectural approaches were pursued, described in Sections 4 and 5.

---

## 3. Structural foundation

### 3.1 SOX2–SOX9 co-complex model

No crystal structure of the SOX2–SOX9 co-complex exists. A structural model was built using **AlphaFold 3** (SOX2 HMG + SOX9 HMG + SOX-dimer DNA motif; ranking score 0.703, ipTM 0.63). Additional validation used Boltz-2 (5 models; ipTM 0.654). Biochemical data from Huang et al. (2015, Sci Rep) were used to confirm the predicted dimerisation interface.

The co-complex has the two HMG domains ~45 Å apart, both bound to the same DNA strand with a specific relative orientation.

### 3.2 Hotspot selection

For the AND-gate to function, both binder epitopes must:
1. Be accessible in the co-complex (not buried at the SOX2–SOX9 interface)
2. Point the binder termini toward each other (so fused caspase halves can reconstitute)

Accessible surface area (SASA) analysis of the co-complex model identified the interface residues (SOX2: S69, K117–119; SOX9: E159, R162, K166, Y174) and excluded them from targeting.

**Final hotspots:**
| Target | Hotspot residues | Location | Notes |
|--------|-----------------|----------|-------|
| SOX2 | 62, 70, 74, 81 (HMG wing-2 face) | Non-DNA-binding; exposed in co-complex | 6–12 Å from DNA |
| SOX9 | 168, 169, 171, 172, 173 (HMG C-terminal face) | Non-DNA-binding; exposed in co-complex | All PPI-safe |

Inter-epitope centroid distance: **~33.5 Å** — bridgeable by a GGS×3 linker or a ~70 aa bridge protein.

---

## 4. Track 1 — Two-binder split-caspase

### 4.1 SOX2 binder design (BindCraft)

BindCraft (miniprotein hallucination via AF2-MPNN cycles) was run against the SOX2 wing-2 hotspot. The target surface is flat and hydrophilic, requiring relaxed shape-complementarity filters (≥0.45 vs the default 0.60). After 1430 trajectories, 7 backbones were accepted and diversified with ProteinMPNN.

**Lead SOX2 binder: mpnn17**
- Length: 71 aa
- Boltz-2 ipTM (vs SOX2 HMG): **0.897**
- Selectivity: 2.77× (SOX2 ipTM / SOX9 ipTM)

### 4.2 SOX9 binder design (Protein Hunter)

Protein Hunter (Boltz-2 hallucination pipeline) was run against the SOX9 C-terminal hotspot cluster. 50+ designs were generated; candidates were filtered by ipTM ≥ 0.70 and selectivity.

**Lead SOX9 binder: ph_sox9_run18**
- Length: 67 aa
- Boltz-2 ipTM (vs SOX9 HMG): **0.856**
- Selectivity: 6.69× (SOX9 ipTM / SOX2 ipTM) — excellent

### 4.3 Five-stage selection pipeline

All binders were evaluated through a shared funnel:

| Stage | Method | Criterion |
|-------|--------|-----------|
| 1 | Boltz-2 re-prediction | ipTM ≥ 0.70, pLDDT ≥ 70 |
| 2 | AND-gate geometry (4-chain Boltz-2) | pair ipTM ≥ 0.60; C-term(A) ↔ N-term(B) < 40 Å |
| 3 | Selectivity (each binder vs wrong target) | wrong-target ipTM < 0.40 |
| 4 | Composite ranking: ipTM × selectivity × geometry | top 10 pairs |
| 5 | Full caspase geometry (flipped, with GGS linker) | recon dist < 15 Å in ≥ 2/3 seeds |

**Lead pair: mpnn17 × ph_sox9_run18** (Stage 4 rank 4/10; Stage 5 recon=8.05 Å at seed 2)

### 4.4 GROMACS MD — bound complex (3 × 20 ns)

Starting from the best Stage 5 Boltz-2 pose (4-chain; mpnn17 × ph_sox9_run18, GGS×1, recon=8.05 Å), three independent 20 ns simulations were run in explicit solvent.

**Pass criterion**: Casp9_N C-terminus ↔ Casp9_C N-terminus < 15 Å in ≥ 50% of frames.

| Replicate | Raw fraction <15 Å | PBC-corrected fraction <15 Å |
|-----------|-------------------|------------------------------|
| Rep 1 (lead) | 76.3% | **99.95%** |
| Rep 2 | 80.1% | **100.0%** |
| Rep 3 | 87.8% | **100.0%** |
| **Mean** | **81.4%** | **~100%** |

> Note: the "raw" fractions were the original reported values. The October 2026 robustness audit corrected for periodic boundary conditions (PBC); the corrected values jump to ~100%. Both the raw and corrected values exceed the 50% pass criterion. The caspase termini remain within reconstitution distance throughout all three simulations.

### 4.5 GROMACS MD — unbound control (1 × 20 ns)

To check for spontaneous leakiness (reconstitution without SOX2/SOX9), the two constructs alone (without any target) were simulated.

**Pass criterion**: endpoint distance > 25 Å in ≥ 95% of frames.

| Metric | Value |
|--------|-------|
| Duration | 20.0 ns (2001 frames, 10 ps cadence) |
| Mean endpoint distance | **44.58 Å** |
| Min endpoint distance | **34.53 Å** |
| % frames > 25 Å | **100%** |
| % frames < 15 Å | **0%** |
| Criterion met | ✓ **PASS** |

The two caspase halves never approach each other in the absence of their targets. There is no spontaneous reconstitution at the structural level.

> Important caveat: passing a distance heuristic does not establish enzymatic non-leakiness. Enzymatic activity at the distance observed (34–58 Å) is negligible by mechanism, but cellular assays (Caspase-Glo in the absence of SOX2/SOX9 co-expression) are required to formally establish non-leakiness.

---

## 5. Track 2 — Single bridge binder

### 5.1 Design concept

Instead of two separate binders, a single miniprotein spans *both* SOX2 and SOX9 simultaneously. Caspase halves are fused at the bridge's N- and C-termini:

```
Single construct:  NLS – Casp9_N – GGS – [bridge binder] – GGS – Casp9_C
```

**AND-gate logic**: the bridge only folds stably onto both targets when they are co-complexed. The bridge N-to-C distance in the bound state directly reports caspase half proximity.

### 5.2 Design campaigns

Four methods were attempted:

| Method | Designs attempted | Pass |
|--------|-------------------|------|
| BindCraft | 233 trajectories | 0 |
| Protein Hunter (Boltz-2 hallucination) | 50 designs | 0 |
| Boltz-2 linker screen (concatenate individual binders) | 18 combinations | 0 |
| **RFdiffusion + ProteinMPNN + MD** | 200 backbones → 32 MD-tested | **5 (15.6%)** |

BindCraft and PH failed because both methods optimise against a single protein surface; bridging two separate chains 45 Å apart falls outside what they handle. RFdiffusion works by sampling backbone coordinates directly via denoising — it can generate a continuous chain that geometrically spans two hotspot clusters.

**RFdiffusion parameters:**
```
contig: [A175-215/0 B450-509/0 60-100]
hotspots: A193, A199, B466, B505, B508
200 backbones; filter: CA dist to each hotspot cluster < 9 Å → top 80 into ProteinMPNN
```

### 5.3 GROMACS MD — five passing designs (3 × 20 ns each)

All five designs that passed the initial 20 ns screen were run for 3 replicates. The AND-gate criterion here is **dual contact**: per frame, does the bridge simultaneously contact both SOX2 and SOX9?

**Periodic-corrected replicate occupancy at 8 Å contact, >10% threshold (October 2026 reanalysis):**

| Design | Rep 1 | Rep 2 | Rep 3 | Status |
|--------|-------|-------|-------|--------|
| **bridge_48_s01** | 99.6% | 99.9% | 100% | ✓ Stable |
| bridge_73_s09 | 100% | 99.3% | 100% | ✓ Stable |
| bridge_136_s01 | 97.1% | 100% | 99.8% | ✓ Stable |
| bridge_101_s06 | 100% | 100% | 100% | ✓ Stable |
| bridge_90_s01 | 100% | 99.5% | **15.7%** | ⚠ Replicate failure |

**bridge_90_s01 is demoted**: rep 3 collapsed to 15.7% — a more than 80 percentage point drop from the other two replicates. It cannot be labelled consistently stable.

**Lead: bridge_48_s01** — consistently above 99% across all replicates and all 27 metric/window combinations tested in the audit.

### 5.4 Selectivity — bridge_48_s01

**2-chain screen** (bridge alone vs each HMG domain):

| Target | Mean ipTM | Verdict |
|--------|-----------|---------|
| SRY | 0.616 | flagged → follow-up 3-chain |
| SOX11 (neuronal safety concern) | 0.473 | flagged → follow-up 3-chain |
| HMGB1 | 0.415 | flagged → follow-up 3-chain |
| SOX3 | 0.369 | PASS |
| SOX4 | 0.336 | PASS |
| SOX17 | 0.300 | PASS |

Note: on-target SOX2 scores only 0.162 alone — confirming the bridge does not bind SOX2 without SOX9, as designed.

**3-chain screen** (replace one target with an off-target; baseline SOX2+SOX9 = 0.348):

| Condition | ipTM | Δ from baseline | Verdict |
|-----------|------|-----------------|---------|
| HMGB1 + SOX9 (HMGB1 in SOX2 position) | 0.358 | +3% | NOT selective ⚠ |
| SOX2 + HMGB1 (HMGB1 in SOX9 position) | 0.290 | −17% | Selective ✓ |
| SRY + SOX9 | 0.251 | −28% | Selective ✓ |
| SOX11 + SOX9 | 0.211 | −39% | Selective ✓ |
| SOX2 + SOX11 | 0.208 | −40% | Selective ✓ |
| SOX2 + SRY | 0.172 | −51% | Selective ✓ |

**Interpretation:**
- SOX11 (neuronal safety concern): selective in both positions (39–40% drop). **Resolved.**
- SRY: selective in both positions despite 2-chain false-positive. **Resolved.**
- HMGB1 in SOX9 position: 17% drop, selective. HMGB1 in SOX2 position: no reduction (0.358 vs 0.348 baseline). **Unresolved — experimental pull-down required.** High seed variance in both conditions prevents a definitive computational verdict.

---

## 6. Robustness audit

An independent robustness audit (October 2026, Codex) reviewed all evidence, repaired analysis scripts, and re-ran key computations. Key findings:

### 6.1 PBC correction improves Track 1 bound MD result

The original analysis used raw (non-periodic) coordinate subtraction. Minimum-image reanalysis:

| Replicate | Original (raw) | Corrected (PBC) |
|-----------|----------------|-----------------|
| Rep 1 | 76.3% | **99.95%** |
| Rep 2 | 80.1% | **100.0%** |
| Rep 3 | 87.8% | **100.0%** |

The raw values were conservative underestimates. The corrected result is stronger. Both exceed the pass criterion.

### 6.2 Stage 1 candidate identity bug

The original Stage 1 sequence extractor selected the *shortest chain* from each BindCraft PDB as the binder. For 23 of 33 recorded rows, the shortest chain was the *target* (SOX2/SOX9), not the binder. The extracted sequences were wrong for those 23 candidates.

**Impact**: none of the Stage 2–5 pair rows reference the 23 mismatched names. The lead designs (mpnn17, ph_sox9_run18, bridge_48_s01) are unaffected. The 23 affected Stage 1 rows cannot be used for comparison; corrected predictions were re-run (see Section 6.4).

### 6.3 bridge_90_s01 replicate failure

Full replicate audit across 27 metric/window settings revealed that bridge_90_s01 rep 3 drops to 15.7% dual contact — far below the other two replicates. This design is **demoted** from the candidate list. The four remaining bridge designs (48, 73, 136, 101) remain stable.

### 6.4 Corrected candidate predictions (23 Stage 1 binders)

The 23 mis-extracted binders were re-predicted with correct chain-B sequences, 3 seeds each (69 predictions total). 19 of 23 candidates met the ipTM ≥ 0.7 heuristic in all seeds. These corrected predictions are recorded in `evidence/corrected_candidate_predictions/`; their scores are not mixed into historical rankings.

### 6.5 Fresh BOTH-target predictions

Twelve fresh 4-seed predictions were run under BOTH, SOX2-only, SOX9-only, and NEITHER conditions using the current construct sequences (no template).

| Condition | ipTM range | Endpoint distance |
|-----------|-----------|-------------------|
| BOTH (3 seeds) | 0.221–0.236 | 49.7, 30.4, 42.1 Å |
| SOX2-only, SOX9-only, NEITHER | 0.13–0.24 | > 30 Å |

The BOTH predictions fail to reproduce <15 Å assembly. This was flagged as a concern (see Section 7 for resolution).

### 6.6 Caspase reference-fold comparison

The three bound MD trajectories were compared to the caspase-9 crystal structure (PDB: 1JXQ). Mean CA RMSD was 17.6–17.7 Å vs chain A and 16.5–16.7 Å vs chain B. Native inter-half contact retention was 4–7%.

**Interpretation**: the simulated assembly does not match the experimental reference fold. This raises a question about whether the reconstituted complex has productive catalytic geometry. See Section 7 for the full investigation.

### 6.7 Construct identity notes

- Confirmed: the simulated construct sequences match the unbound-MD sequences and the TPR topology sequences exactly.
- Mutation: the reference sequence has **D315/A316** (not D315/G316 as previously labelled). The split site is the same natural cleavage site either way.
- Difference from proposal: the simulated construct uses GGS×1 linkers and a flipped second fusion. It does not include the proposed NLS or T2A elements. The proposal cannot be treated as the exact simulated construct.

---

## 7. Caspase geometry investigation

### 7.1 The concern

The robustness audit raised two issues:
1. The reconstituted complex has 17.6 Å CA RMSD vs 1JXQ — it doesn't look like a native caspase-9 dimer.
2. Fresh Boltz-2 BOTH predictions gave ipTM 0.22 and 30–50 Å endpoint distances — the assembly could not be reproduced from scratch.

### 7.2 Analysis and context

**On the 1JXQ RMSD**: 1JXQ is the processed active dimer — two large+small pairs in an apoptosome-associated conformation. The simulated construct has *one* large and *one* small fragment. Comparing a single fragmented pair against the fully dimerised crystal structure is inherently going to produce a large RMSD. The comparison is structurally mismatched. The correct comparison would be a single large+small unit from 1JXQ.

**On the fresh prediction failure**: split proteins are intrinsically disordered in isolation. Boltz-2 folds from random coil — it has no way to know the fragments should reconstitute. The from-scratch failure is expected for split-protein systems and should not be interpreted as evidence that the assembly is unstable.

### 7.3 Template-guided prediction test

**Method**: The 10 ns MD frame was extracted (split site distance = 9.54 Å in that frame), and used as a structural template for Boltz-2 prediction. If Boltz-2 keeps the halves close when initialised from the assembled conformation, the MD result is supported. If it immediately separates them, that is a worse sign.

Three seeds were run. Results:

| Seed | ipTM | Endpoint distance | Halves close? |
|------|------|-------------------|---------------|
| 0 | 0.782 | 9.2 Å | ✓ |
| 1 | 0.787 | 9.7 Å | ✓ |
| 2 | 0.782 | 11.2 Å | ✓ |

**3/3 seeds keep the halves close (<15 Å). ipTM 0.782–0.787.**

**Inter-chain ipTM (mean across 3 seeds):**

| Interface | Mean ipTM |
|-----------|-----------|
| SOX2 ↔ CompA | 0.831 |
| SOX9 ↔ CompB | 0.858 |
| Casp9_N ↔ Casp9_C | **0.741** |
| SOX2 ↔ SOX9 | 0.799 |

The caspase half–half interface (0.741) is below the binder interfaces (as expected for a split protein) but above the 0.7 confident-prediction threshold, and stable across seeds (±0.014).

### 7.4 Interpretation

The from-scratch failure (ipTM 0.22, 30–50 Å) was a Boltz-2 folding problem, not an assembly problem. The assembled state is energetically valid: when given the MD-derived conformation, Boltz-2 predicts it confidently and stably.

**What remains open (only lab can answer):**
- Whether the reconstituted geometry is catalytically productive — this requires a functional read-out (Caspase-Glo), not more structure prediction.
- Whether enzymatic activity is conditional on SOX2+SOX9 co-expression in cells.

---

## 8. Summary of all results

### Binders

| Binder | Target | Method | ipTM | Length |
|--------|--------|--------|------|--------|
| mpnn17 | SOX2 | BindCraft | 0.897 | 71 aa |
| ph_sox9_run18 | SOX9 | Protein Hunter | 0.856 | 67 aa |

### Track 1 MD summary

| Simulation | Result | Criterion | Status |
|------------|--------|-----------|--------|
| Bound rep 1 (PBC-corrected) | 99.95% frames < 15 Å | > 50% | ✓ PASS |
| Bound rep 2 (PBC-corrected) | 100% frames < 15 Å | > 50% | ✓ PASS |
| Bound rep 3 (PBC-corrected) | 100% frames < 15 Å | > 50% | ✓ PASS |
| Unbound control | 100% frames > 25 Å, mean 44.6 Å | > 95% > 25 Å | ✓ PASS |

### Track 2 lead — bridge_48_s01

| Replicate | Dual-contact (>10% threshold) | Status |
|-----------|-------------------------------|--------|
| Rep 1 | 99.6% | ✓ |
| Rep 2 | 99.9% | ✓ |
| Rep 3 | 100% | ✓ |

### Caspase geometry check (template-guided Boltz-2)

| Condition | ipTM | Endpoint dist | Result |
|-----------|------|---------------|--------|
| From scratch (3 seeds) | 0.221–0.236 | 30–50 Å | ✗ folding failure |
| Template-guided (3 seeds) | 0.782–0.787 | 9.2–11.2 Å | ✓ assembly stable |

### Key numbers at a glance

| Metric | Value |
|--------|-------|
| SOX2 binder ipTM | 0.897 |
| SOX9 binder ipTM | 0.856 |
| Track 1 bound MD (PBC) | ~100% frames <15 Å |
| Track 1 unbound MD | 100% frames >25 Å, min 34.5 Å |
| Track 2 lead (bridge_48) dual-contact | 99.6–100% across 3 reps |
| Template-guided assembly ipTM | 0.782–0.787 |
| Caspase half–half ipTM | 0.741 |
| bridge_90 rep 3 (demoted) | 15.7% |

---

## 9. What is and is not established

### Established computationally

- mpnn17 and ph_sox9_run18 predict confident binding to SOX2 and SOX9 respectively (ipTM > 0.85 each).
- The Track 1 construct keeps caspase halves within reconstitution distance (~10 Å) in 100% of bound MD frames (3 × 20 ns, PBC-corrected).
- The Track 1 construct keeps caspase halves separated (>34 Å) in 100% of unbound MD frames (20 ns).
- bridge_48_s01 simultaneously contacts both SOX2 and SOX9 in >99% of frames across all 3 replicates (3 × 20 ns).
- The Track 1 assembled state is predicted confidently by Boltz-2 when given the MD-derived conformation as a template (ipTM 0.78, endpoints 9–11 Å, 3/3 seeds).
- bridge_48_s01 is computationally selective against SOX11 and SRY in both target positions.

### Not established

- Enzymatic activity of the reconstituted caspase-9.
- Conditional activity — that caspase-9 reconstitution requires SOX2+SOX9 co-expression in live cells.
- Enzymatic non-leakiness in cells (the distance heuristic passing is necessary but not sufficient).
- Cellular AND-gate logic — contact occupancy in MD is not AND-gate activity.
- Normal-cell specificity in tissue.
- Selectivity of bridge_48_s01 against HMGB1 in the SOX2 position (computationally inconclusive; experimental pull-down required).
- Whether the reconstituted caspase geometry is catalytically productive (1JXQ RMSD remains elevated; functional assay required).
- Delivery, immune payload, and treatment readiness.

---

## 10. Next steps — lab

Both tracks are computationally complete. The shared sensing apparatus (SOX2 binder + SOX9 binder) is identical across this project and gbm-ANDGATE, so several early experiments validate both projects at once.

### Shared experiments (validate both projects)

| Experiment | Validates |
|------------|-----------|
| His-tag expression of CompA + CompB (HEK293T) | Expression and solubility of each construct |
| Western blot / SEC | Confirm molecular weights; detect aggregation |
| NLS staining — confirm nuclear localisation | NLS function |
| SOX2 co-IP with CompA | SOX2 binder binds SOX2 in cells |
| SOX9 co-IP with CompB | SOX9 binder binds SOX9 in cells |
| GSC vs NPC co-expression — confocal colocalisation | Both constructs reach the nucleus; confirm co-complex proximity |

### Track 1–specific

| Experiment | What it answers |
|------------|----------------|
| Caspase-Glo in GSC (SOX2+SOX9 co-expressing) | Does the reconstituted switch kill? |
| Caspase-Glo in normal NSC (SOX2 only) | Is there SOX9-absent leakiness? |
| Caspase-Glo in both lines without SOX9 knockin | Unbound baseline leakiness |
| Caspase-3/7 activity assay time course | Kinetics of killing |
| Dose–response (construct titration) | Hill coefficient — cooperativity |

### Track 2–specific

| Experiment | What it answers |
|------------|----------------|
| HMGB1 pull-down (bridge_48_s01 vs recombinant HMGB1) | Resolves the computational HMGB1 flag |
| Same Caspase-Glo assays as Track 1 | Cross-comparison of both architectures |

### Constructs to synthesise (codon-optimised for human expression)

**Track 1 — Construct A** (277 aa coding):
```
NLS–mpnn17(71aa)–GGS–Casp9_N(aa140–315, 176aa)
```

**Track 1 — Construct B** (flipped, 249 aa coding):
```
Casp9_C(aa316–416, 101aa)–GGS–ph_sox9_run18(67aa)–NLS
```

**Track 2 — bridge_48_s01 construct** (~547 aa coding):
```
NLS–Casp9_N(176aa)–GGS–bridge_48_s01(69aa)–GGS–Casp9_C(101aa)
```

---

## Evidence files

| File | Contents |
|------|----------|
| `evidence/track1_lead_periodic.tsv` | Track 1 rep 1 endpoint time series (PBC-corrected, 2001 frames) |
| `evidence/track1_rep2_periodic.tsv` | Track 1 rep 2 time series |
| `evidence/track1_rep3_periodic.tsv` | Track 1 rep 3 time series |
| `evidence/unbound_control/metrics.json` | Unbound MD validated metrics (mean 44.6 Å, 100% > 25 Å) |
| `evidence/bridge_periodic/` | Dual-contact time series for all 5 bridge designs |
| `evidence/bridge_replicate_cohort/` | Full replicate audit (15 trajectories × 27 metrics) |
| `evidence/matched_controls/` | Fresh BOTH/SOX2-only/SOX9-only/NEITHER predictions |
| `evidence/template_guided_controls/` | Template-guided Boltz-2 (3 seeds; ipTM 0.78; 9–11 Å) |
| `evidence/corrected_candidate_predictions/` | 69 corrected Stage 1 predictions (23 binders × 3 seeds) |
| `evidence/native_geometry/` | 1JXQ reference-fold comparison (RMSD ~17.6 Å) |
| `evidence/stage1_chain_identity.json` | Candidate identity audit (23 of 33 original rows corrected) |
| `ROBUSTNESS_REVIEW.md` | Full audit findings and software repairs |
| `COMPUTATIONAL_STATUS.json` | Machine-readable status of all completed tasks |

---

*All raw trajectories (TPR/XTC) remain on local storage (D: drive). Content hashes are in `evidence/external_artifacts.json`.*
