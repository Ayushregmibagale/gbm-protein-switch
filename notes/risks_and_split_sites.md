# Risks and Split Sites — Split-Caspase-9 AND Gate

## Risks

**1. Leaky reconstitution (false positives)**
The N-half contains the CARD domain, which can interact with endogenous APAF-1
and trigger the apoptosome independently of the AND gate. Even without SOX9
binding, the N-half alone could cause background apoptosis in normal cells.
Mitigation: use a CARD-deleted N-half (residues 140-315 only).

**2. SOX2 expression in normal NSCs**
SOX2 is expressed in adult neural stem cells. If SOX9 is even transiently
co-expressed in NSCs (e.g. during injury/stress responses), off-target killing
could occur. Needs careful validation across normal brain cell types before
any in vivo work.

**3. Tumour heterogeneity / escape**
Not all GSCs may have the SOX2/SOX9 co-complex at all times. Tumours can
evolve to downregulate either TF and escape. Should be combined with other
modalities (e.g. standard-of-care TMZ, immune component).

**4. AAV cargo and delivery**
Two constructs (binder-casp9-N + binder-casp9-C) must fit in AAV (~4.7 kb
limit). Options:
- Dual-AAV system (each half in a separate vector) — stoichiometry problem
- Single bicistronic construct with T2A self-cleaving peptide (preferred)
BBB penetration and uniform GSC transduction are also not guaranteed.
Best serotypes for CNS: AAV9, AAVrh10, AAV-PHP.eB.

**5. Binder-caspase folding interference**
The miniprotein fused to a caspase half may misfold or sterically block
reconstitution. Linker length and rigidity between binder and caspase half
is critical and will need empirical optimisation (GGS vs rigid linkers).

**6. Nuclear trafficking**
Both TFs are nuclear. Fusion constructs need a functional NLS and must avoid
proteasomal degradation before reaching the nucleus.

**7. Immunogenicity**
Designed miniproteins are non-human sequences. Repeated AAV dosing could
trigger neutralising antibodies against both the vector and binder sequences.
Consider humanisation of binder sequences post-design.


## Caspase-9 Domain Structure

Human caspase-9 (416 aa):

  |--CARD--|--linker--|----Large (p35)----|--linker--|--Small (p10)--|
  1       97         140                315  330                   416


## Validated Split Sites

| Split         | N-half          | C-half       | Notes                                      |
|---------------|-----------------|--------------|--------------------------------------------|
| D315          | 1-315           | 316-416      | Most common; N-half can bind APAF-1 via CARD |
| D330          | 1-330           | 331-416      | Alternative processing site                |
| CARD-deleted  | 140-315 (large) | 316-416      | Safer — no APAF-1 interaction; RECOMMENDED |


## Recommended Design: CARD-deleted split at D315

- N-half (residues 140-315, large subunit only) fused to SOX2 binder C-terminus
- C-half (residues 316-416, small subunit) fused to SOX9 binder C-terminus
- Neither half can trigger the apoptosome independently
- Reconstitution only occurs when both binders co-localise on SOX2/SOX9 complex
- Reduces false-positive apoptosis risk significantly vs. full-length split

## Construct architecture (single AAV, bicistronic)

  NLS-[SOX2_binder]-[GGS linker]-[Casp9 140-315]-T2A-NLS-[SOX9_binder]-[GGS linker]-[Casp9 316-416]

The T2A self-cleaving peptide ensures equimolar expression of both halves
from a single ORF, avoiding stoichiometry problems of the dual-AAV approach.
