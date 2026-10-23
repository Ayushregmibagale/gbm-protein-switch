# GBM Protein Switch — Project Aims

## One-line concept
Design a dual-selectivity AAV therapy that kills glioblastoma stem cells (GSCs) via a protein coincidence detector (split-caspase-9) and activates durable immune memory (IL-12), both locked to the SOX2/SOX9 co-complex that defines the GSC state.

## Scientific rationale
Koeber et al. (Nature, 2026 — s41586-026-10329-6) showed that SOX2 and SOX9 form higher-order multimeric complexes specifically in GSCs, driving SSE (synthetic super-enhancer) activity. They exploited this at the DNA level (SSE-7 drives HSV-TK + IL-12 via AAV) achieving curative outcomes + immunological memory in mice. This project adds a second, protein-level selectivity lock — making the system safer and more mechanistically novel.

## Delivery system — two-component AAV

### Component 1 (novel) — Protein coincidence detector
- Split-caspase-9 **CARD-deleted** design (see notes/risks_and_split_sites.md):
  - N-half: residues 140–315 (large subunit only, no CARD) fused to SOX2 binder C-terminus
  - C-half: residues 316–416 (small subunit) fused to SOX9 binder C-terminus
  - CARD deletion prevents N-half from interacting with APAF-1 independently
  - Split site: D315/G316 — most validated in literature
- In GSCs: SOX2 and SOX9 are co-complexed (epitopes 45 Å apart) → binder arms forced into proximity → caspase-9 reconstitutes → activates caspase-3/7 → apoptosis + immunogenic antigen release
- In normal cells: TFs not co-complexed → no reconstitution → no toxicity
- Construct architecture (single bicistronic AAV):
  `NLS-[SOX2_binder]-[GGS linker]-[Casp9 140-315]-T2A-NLS-[SOX9_binder]-[GGS linker]-[Casp9 316-416]`
- T2A self-cleaving peptide ensures equimolar expression of both halves
- Encoded as AAV transgene (no systemic drug required — replaces HSV-TK/GCV)

**Why caspase-9 not caspase-3/7:**
Caspase-9 is the initiator; one reconstituted dimer activates many caspase-3/7 molecules (signal amplification). Caspase-3/7 as the split protease risks leaky killing from partial reconstitution. See notes/caspase_choice_rationale.md.

### Component 2 (from Koeber et al.) — Immune activation payload
- IL-12 expression cassette driven by SSE-7 (validated GSC-selective enhancer from the paper)
- Optionally: add IL-15 under SSE-7 to sustain CD8+ memory T cells beyond what IL-12 alone achieves
- Both IL-12 and the protein switch transgene delivered in the same AAV vector
- IL-12 → activates NK cells + CD8+ T cells → immunological memory → prevents tumour recurrence on re-challenge

### Why two locks are better than one
| | Koeber et al. | This project |
|--|--|--|
| Selectivity mechanism | SSE-7 (DNA level only) | SSE-7 (DNA) + split-caspase switch (protein level) |
| Cytotoxicity | HSV-TK + systemic GCV | Direct caspase-9 (no systemic drug needed) |
| Immune activation | IL-12 via SSE-7 | IL-12 (±IL-15) via SSE-7 |
| Risk of normal brain toxicity | CMV-IL-12 lethal; SSE-7 safe | Double-locked — both DNA and protein selectivity required |

## Available structures

### SOX2 (UniProt P48431)
| PDB | Resolution | What it shows | Use in this project |
|-----|-----------|--------------|-------------------|
| 1O4X | 2.9 Å | OCT1+SOX2 HMG+DNA ternary | HMG domain in native DNA-binding pose |
| 6T7B | cryo-EM | SOX2 HMG + nucleosome | Chromatin-bound SOX2 surface |
| 6WX7 | crystal | SOX2 + Importin-α2 | **Best for binder design** — shows accessible non-DNA-binding faces |
| 6WX8 | crystal | SOX2 + Importin-α3 | Alternative accessible surface view |
| 8BX1/8BX2 | crystal | OCT4/SOX2/DNA | Most recent; shows co-complex interface |
| AF-P48431 | AlphaFold | Full-length SOX2 | Covers disordered transactivation domain |

### SOX9 (UniProt P48436)
| PDB | Resolution | What it shows | Use in this project |
|-----|-----------|--------------|-------------------|
| 4EUW | 2.77 Å | SOX9 HMG + DNA | **Best resolution HMG structure** — primary input for binder design |
| 4S2Q | crystal | SOX9 HMG + ChIP-seq DNA element | Biologically relevant bound conformation |
| AF-P48436 | AlphaFold | Full-length SOX9 incl. dimerization domain | Critical — DIM domain drives higher-order complex formation |

### SOX2-SOX9 co-complex (no crystal structure exists)
- Huang et al. 2015 (Sci Rep, PMC4445065) characterised the SOXE dimerisation interface biochemically
- Model using AlphaFold3 multimer (SOX2:SOX9:DNA) as starting point
- This is a key deliverable — the co-complex model defines which surfaces to target

## Aims

### Aim 1 — Structural modelling of the SOX2-SOX9 co-complex
- Download 6WX7 (SOX2) and 4EUW (SOX9) as primary structural inputs
- Run AlphaFold3 multimer prediction: SOX2 HMG + SOX9 HMG + SOX dimer DNA motif
- Identify accessible surfaces on each TF that are (a) not at the DNA-binding interface and (b) not buried at the SOX2-SOX9 dimer interface
- These surfaces become the BindCraft target epitopes

### Aim 2 — Design and validate SOX2 and SOX9 binders
- Run BindCraft/RFdiffusion targeting the epitopes identified in Aim 1
- Prioritise binders with predicted Kd < 100 nM and no steric clash with DNA or co-complex interface
- MD simulation of binder:TF binary complexes (20 ns minimum) to assess stability
- Validate: contact frequency, ΔΔG, RMSD of binder on TF surface

### Aim 3 — Model and validate the quaternary complex
- Dock binder-A:SOX2:SOX9:binder-B using RosettaDock / pRosettaC
- Run MD on top poses — key criterion: split-caspase-9 C-termini (from half-A and half-B) must come within ~10 Å for reconstitution
- Confirm geometric compatibility across ≥20% of MD frames (equivalent to P_stable criterion from PROTAC work)

### Aim 4 — Design the full switch construct
- Split site: D315/G316, CARD-deleted (N-half = residues 140-315; C-half = 316-416)
- Architecture: NLS-[binder]-[GGS linker]-[casp9 half]-T2A-NLS-[binder]-[GGS linker]-[casp9 half]
- Design and test GGS vs rigid linkers between binder and caspase half
- Add SV40 NLS to both halves for nuclear localisation
- Run MD on full construct in absence of SOX2/SOX9 to check for spontaneous reconstitution (leakiness)
- Leakiness criterion: split-half termini must stay >25 Å apart in >95% of frames when unbound
- Run MD in presence of SOX2/SOX9 co-complex — reconstitution criterion: termini <10 Å in >20% of frames

### Aim 5 — ML-guided optimisation
- Adapt XGBoost/GNN pipeline from PROTAC project to rank binder candidates
- Features: binder:TF binding energy, co-complex geometry, split-half proximity distance distribution from MD, binder surface area, linker geometry
- Use LOPO-CV for validation (as per PROTAC ML approach)

### Aim 6 — Experimental validation (in collaboration)
- Encode full system (protein switch + SSE-7-IL-12 ± IL-15) in AAV1 vector
- Test in GSC panel (G7, E28, E55, E31) vs normal controls (ioGlutamatergic neurons, ioMicroglia, fibroblasts, normal NSCs)
- Readouts: caspase-3/7 activation, cell viability, IL-12 secretion (ELISA), selectivity index
- Re-challenge experiment to confirm immunological memory
- Compare double-lock system vs SSE-7-only (Koeber et al.) for selectivity in normal NSCs

## Key risks and mitigations

| Risk | Mitigation |
|------|-----------|
| SOX2/SOX9 also expressed in normal NSCs | Target co-complex interface surface (only forms in GSCs); double-lock with SSE-7 provides backup selectivity |
| CARD domain triggers APAF-1 independently | Use CARD-deleted N-half (residues 140-315 only) — eliminates independent apoptosome activation |
| Split-caspase reconstitutes spontaneously | CARD-deleted design + MD leakiness check in Aim 4 (split-half termini must stay >25 Å in >95% frames when unbound) |
| SOX9 transiently expressed in injured NSCs | Validate co-complex specificity in stress conditions; SSE-7 double-lock provides backup |
| Tumour heterogeneity / TF downregulation | Combine with standard-of-care TMZ + immune component (IL-12) for multi-modal attack |
| No crystal structure of SOX2-SOX9 co-complex | Boltz-2 model (ipTM=0.654, 5 models); AlphaFold Server for AF3 validation; Huang et al. 2015 biochemical data |
| Binder access blocked by DNA binding | Epitopes A:59-83 and B:141-172 confirmed non-DNA-contacting by SASA analysis of Boltz-2 model |
| Binder-caspase folding interference | Test GGS vs rigid linkers; validate each half independently before combining |
| AAV cargo size limit (~4.7 kb) | Single bicistronic construct with T2A; miniproteins ~60-100 aa each; split-caspase halves ~175+100 aa = ~850 bp total coding |
| Stoichiometry of two halves | T2A ensures equimolar expression from single ORF; avoids dual-AAV stoichiometry problem |
| Immunogenicity of designed miniproteins | Humanise binder sequences post-design; screen against human proteome for off-target homology |
| IL-12 toxicity in normal CNS | SSE-7 selectivity validated in Koeber et al.; protein-level double-lock adds additional safety layer |

## Key literature
- Koeber et al., Nature 2026 — SSE-7, SOX2/SOX9 co-complex in GSCs, IL-12 immune memory (primary inspiration)
- Huang et al., Sci Rep 2015 (PMC4445065) — SOXE dimerisation interface
- Scott et al., Nat Neurosci 2010 — SOX9 in neural stem cells
- Bhattacharyya et al. — split-protein coincidence detectors
- Diao et al. — split-caspase-9 design and split site characterisation

## Computational tools
- BindCraft — curated miniprotein binder design with AF2 scoring (Aims 1-2)
- Protein Hunter (Boltz edition) — high-throughput binder design via Boltz-2 hallucination (Aims 1-2)
- AlphaFold Server (alphafoldserver.com) — SOX2:SOX9:DNA co-complex (AF3, 10 free jobs/day)
- pRosettaC / RosettaDock — quaternary complex docking (Aim 3)
- GROMACS + AMBER ff — MD simulations (Aims 2-4)
- XGBoost + GNN — ML-guided candidate ranking (Aim 5)

## Notes files
- notes/caspase_choice_rationale.md — why caspase-9 not caspase-3/7
- notes/risks_and_split_sites.md — full risk table + validated split sites + recommended construct architecture

## Deliverables (2-year PhD horizon)
1. AlphaFold3 model of SOX2:SOX9 co-complex with validated epitope surfaces
2. ≥3 validated miniprotein binders each for SOX2 and SOX9
3. MD-verified quaternary complex geometry (binder-A:SOX2:SOX9:binder-B)
4. Full switch construct design with leakiness assessment
5. ML model for binder optimisation
6. Experimental collaboration data: GSC selectivity + immune activation + re-challenge
7. Paper: "A dual-selectivity protein switch kills glioblastoma stem cells and drives immune memory via SOX2/SOX9 co-complex sensing"

## Notes
- Inspired by: Koeber et al. Nature 2026 (s41586-026-10329-6)
- Project started: 2026-08-26
- Primary investigator: Ayush Regmi Bagale (PhD, University of Edinburgh)
- Project folder: ~/gbm-protein-switch/
