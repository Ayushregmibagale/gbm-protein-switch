# How Caspase-9 is Added to the Binder — Fusion Workflow

## Overview

The caspase-9 halves are NOT designed by BindCraft. They are fused to accepted binder
sequences as a post-design computational step.

## Step-by-step

### Step 1 — Get binder sequences
BindCraft / Protein Hunter produces accepted miniprotein sequences (60–100 aa) with:
- ipTM > 0.7 (good binding to SOX2 or SOX9 epitope)
- pLDDT > 70 (well-folded binder)
- No steric clash with DNA or co-complex interface

### Step 2 — Build fusion constructs
Append caspase half to binder C-terminus via a flexible linker:

  Construct 1: NLS – SOX2_binder – [GGS]₅ – Casp9[140–315]
  Construct 2: NLS – SOX9_binder – [GGS]₅ – Casp9[316–416]

Caspase-9 sequence source: UniProt P55211 (human CASP9)
- N-half (large subunit): residues 140–315 (176 aa, CARD deleted)
- C-half (small subunit): residues 316–416 (101 aa)

### Step 3 — Predict folding of each fusion
Use AlphaFold3 server (alphafoldserver.com) to predict each full fusion construct alone.
Criterion: binder pLDDT should NOT drop significantly (<5 units) vs binder-only prediction.
Reject candidates where caspase half destabilises the binder or vice versa.

### Step 4 — Quaternary complex modelling
Dock the full system: binder1-Casp9N : SOX2 : SOX9 : binder2-Casp9C
Tools: pRosettaC or RosettaDock (same pipeline as PROTAC ternary complexes)
Screen criterion: D315/G316 split interface must be within ~10 Å in top-scoring poses.

### Step 5 — MD validation (Aim 4)
Run GROMACS MD on top quaternary complex poses:
- Bound state (SOX2/SOX9 co-complex present):
  - Criterion: split-half termini <10 Å in ≥20% of frames = reconstitution-capable
- Unbound state (binder fusion alone, no SOX2/SOX9):
  - Criterion: split-half termini >25 Å in >95% of frames = not leaky

### Step 6 — Linker optimisation
Start with [GGS]₅ (flexible). If unbound leakiness > 5%, test rigid helical linker.
If bound reconstitution < 20%, try longer [GGS]₇ to allow halves to find each other.

## Key parameters

| Parameter | Value |
|-----------|-------|
| N-half (Casp9) | P55211 residues 140–315 |
| C-half (Casp9) | P55211 residues 316–416 |
| Split site | D315/G316 (most validated) |
| CARD deletion | N-half starts at 140 (no CARD domain 1-97) |
| Linker (starting) | [GGS]₅ (~15 aa, ~45 Å fully extended) |
| Bound criterion | Termini <10 Å in ≥20% MD frames |
| Unbound criterion | Termini >25 Å in >95% MD frames |
| NLS sequence | SV40 NLS (PKKKRKV) prepended to each construct |

## Full bicistronic construct (single AAV)

  NLS – SOX2_binder – [GGS]₅ – Casp9[140–315] – T2A – NLS – SOX9_binder – [GGS]₅ – Casp9[316–416]

T2A self-cleaving peptide ensures equimolar expression of both halves from one ORF.
Approximate size: ~850 bp coding (well within 4.7 kb AAV limit).
