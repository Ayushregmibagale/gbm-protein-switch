# GBM protein-switch robustness review — 2026-10-29

Scope: this project only. Recovered sources were read from
`/mnt/d/gbm-protein-switch`; the WSL home checkout is the intended source-code
location. Large simulation data remain on D:. Identity-corrected and matched-control predictions were run in fresh seeded
batches. No new MD production was launched; the externally resumed unbound run
was allowed to finish. See [WORKFLOW.md](WORKFLOW.md) for operating instructions.

## Current findings

### Candidate identities need quarantine

The original Stage 1 extractor chose the shortest PDB chain as the binder. Among
33 recorded BindCraft rows, **23 instead contain target chain A's sequence**, not
binder chain B. This includes equal-length target/binder cases; chain length
cannot identify the binder. The repaired extractor explicitly reads chain B and
rejects missing or incomplete chains.

`analysis/audit_candidate_identity.py` reproduces the identity comparison from the
accepted PDBs. `evidence/stage1_chain_identity.json` records all 33 comparisons and
source hashes. `selection/validated_stage1_results.tsv` excludes the 23 mismatches;
no scores have been attached to replacement sequences or invented.

None of the recorded Stage 2–5 pair rows references the 23 mismatched names, so this
finding does not by itself invalidate those recorded pairs or the selected lead.
It does invalidate interpretation of the affected Stage 1 rows and comparisons
that include them. The repaired downstream launchers refuse mismatched candidate
identities. An explicit `GBM_STAGE1_TABLE` can select the quarantined table.

The Stage 1 source also compared fractional Boltz pLDDT with a threshold of 70.
It now converts the fractional score to percent for that comparison and records
both scales. Inference now requires new run/result paths to avoid reusing old
predictions after changing candidate identities. The 23 corrected binders were re-predicted using three seeds each; see
`evidence/corrected_candidate_predictions/` for the inputs, models, confidence
files, commands, software versions and checkpoint identity. Historical scores
remain unchanged; new results are not silently mixed into historical rankings.

### Periodic boundaries change the bound distance result

The original scripts use direct coordinate subtraction. Minimum-image reanalysis
of the same three trajectories gives:

| Replicate | Historical raw fraction <15 Å | Periodic fraction <15 Å | Elapsed time |
|---|---:|---:|---:|
| lead_construct | 76.262% | 99.950% | 20.000 ns |
| rep2 | 80.110% | 100.000% | 20.000 ns |
| rep3 | 87.756% | 100.000% | 20.000 ns |

Equal-weight replicate averages are 81.376% raw and 99.983% periodic. The original
raw values reproduce. Each trajectory has 2001 frames spanning 20 ns; the historical
20.01 ns duration counted frames rather than intervals.

Saved endpoint time series and provenance are in `evidence/track1_*_periodic.*`.
They establish proximity of the selected CA atoms. Minimum-image distances refer
to the closest periodic images; productive enzyme geometry still requires an
appropriate assessment of the reconstructed complex. Proximity does not establish
conditional caspase activity or safety.

### All five recorded passing bridge designs have been reanalysed

Historical rounded time series reproduce their original tables. Separately,
periodic raw-trajectory reanalysis gives the following dual-contact occupancies
at the historical fraction threshold of strictly greater than 0.10:

| Design | Historical dual contact | Periodic dual contact |
|---|---:|---:|
| bridge_48_s01 | 99.4% | 99.6% |
| bridge_73_s09 | 99.5% | 100.0% |
| bridge_136_s01 | 86.8% | 97.1% |
| bridge_101_s06 | 99.9% | 100.0% |
| bridge_90_s01 | 99.8% | 100.0% |

`evidence/bridge_periodic/` contains new time series, summaries, and SHA-256 hashes
of the source PDB/TPR/XTC files. Historical tables are preserved. All 15 available trajectories (three per design) have additionally been audited
under 27 cutoff/fraction/window combinations each; see
`evidence/bridge_replicate_cohort/`. At 8 Å / fraction >0.10, replicate occupancies
are 99.6/99.9/100% for bridge_48, 100/99.3/100% for bridge_73,
97.1/100/99.8% for bridge_136, 100/100/100% for bridge_101, and
100/99.5/**15.7%** for bridge_90. The last design has substantial replicate
instability and must not be labelled consistently stable. Frames are correlated; 2001 frames
are not 2001 independent experiments. Contact occupancy does not establish AND logic.

### The complete unbound control passes its distance heuristic

The resumed production reached 20 ns. Validation retained old frames through
4.57 ns and resumed frames from 4.58 ns, excluding the corrupt and overlapping old
tail. The audit verified 2001 frames spanning exactly 0–20 ns at 10 ps cadence,
with source hashes and no gap or duplicated frame times.

Periodic endpoint distances have mean 44.58 Å, minimum 34.53 Å and maximum
58.45 Å. Every frame exceeds 25 Å and none is below 15 Å. The specified distance
criterion is met. This single trajectory does not establish enzymatic non-leakiness
or normal-cell safety. `evidence/unbound_control/` contains the validated metrics,
complete endpoint series and the exact splice/source provenance.

The legacy state file remains `NPT_DONE`; the repaired status command recognises
the validated report after verifying topology identity, exported-series hash, and
source sizes/mtimes. It refuses an unnecessary restart of this audited control. The old log watcher searched the console log for a completion marker
that is actually in `prod.part0002.log`. Repaired analysis uses trajectory coverage
instead of that log-marker assumption. Historical corrupted files and logs were
preserved.

### Caspase-core reference geometry raises a mechanistic concern

The recorded bound actuator is compared with both sequence-aligned caspase-9
conformations in experimental PDB 1JXQ. Missing reference residues are excluded;
each simulated half is reconstructed through consecutive CA minimum-image vectors.
The reference was chosen as an experimental structural comparison, not an activity
classifier. [Experimental entry and primary citation](https://www.rcsb.org/structure/1JXQ)

Across the three bound trajectories, mean assembly CA RMSDs are about 17.6–17.7 Å
against chain A and 16.5–16.7 Å against chain B. Mean retention of reference
inter-half CA contacts is only about 4–7%. Individually fitted halves also have
large deviations. These observations raise a substantial concern about the
predicted fold and assembly, even though endpoint proximity is high. They do not
support a native-like catalytic-core claim. See `evidence/native_geometry/`.

The experimental study concerns dimer-associated activation. Inference: the one
large/small fragment pair represented here is a single catalytic core, so its
proximity metric does not establish the functional oligomeric context. More MD
of the same questionable starting fold would not validate that missing mechanism.
[1JXQ primary-study record](https://www.rcsb.org/structure/1JXQ)

### Construct identities are reconciled against simulation inputs

The bound actuator chains (C/D) exactly match the unbound chains (A/B), and both
match the protein sequences in the respective TPR topologies. The supplied 416-aa
reference has **D315/A316**, correcting the historical D315/G316 label. This is a
check against the supplied reference, not independent verification of a validated
split site in the literature.

The simulated architecture has GGS1 linkers and a flipped second fusion. It does
not include the proposal's added NLS/T2A elements. The original proposal cannot be
treated as the exact simulated construct. See `inputs/README.md` and
`evidence/construct_identity.json` for hashes and full identity details.

### Structural-confidence selectivity remains uncertain

The three recorded on-target seeds have mean ipTM **0.3480**, sample SD **0.2636**,
rather than the historical hard-coded baseline 0.471. For each replacement condition,
baseline-minus-replacement changes sign across seeds. Those seeds do not establish
a consistent direction of confidence change. The launcher now compares with the
measured on-target mean and refuses incomplete comparisons. The audit reports
matched-seed differences. ipTM is structural confidence, not measured affinity,
specificity, or cellular selectivity.

## Software repairs and verification

- Recovered 66 missing source/configuration/result files with original-byte hashes
  in `evidence/recovery_manifest.json`. Their subsequent changes are tracked in Git.
- MD entry points default to status and require explicit `--execute`. Environment
  variables select data roots and tools. Shell command execution was replaced by
  checked argument-list calls; ionisation failures cannot silently advance state.
- State transitions require nonempty expected outputs, reject unknown states, and
  write atomically. Existing production checkpoints resume with the existing TPR;
  ambiguous partial outputs are refused. Fragments require explicit validation.
- Production preprocessing no longer automatically suppresses two GROMACS warnings
  or asks pdb2gmx to tolerate missing atoms. Failures now require investigation.
- Fixed duplicate first epitope patches, residue/coordinate pairing, and misleading
  target-terminus labels. Structure scripts take explicit paths and do not run on import.
- Bridge contacts use periodic distances, actual times and checked atom counts.
- Stage 4 ranking reproduces all ten historical ranks to floating-point precision;
  zero/nonfinite denominators fail instead of creating infinite ranking scores.
- Regression tests cover coverage failures, real XTC fragments, periodic boundaries,
  command failures, checkpoint preservation, artifact hashes and binder identity.
  GitHub Actions checks both dependency-free audit code and the MD analysis environment.

## Cache integrity and matched controls

The historical-cache audit covered 972 runs. All inspected model_0 structures
match their submitted protein input sequences and have valid confidence files.
This includes internally consistent caches whose submitted candidate identity was
wrong; cache consistency does not repair the 23 Stage 1 mistakes. Historical
caches lack complete verified software/checkpoint/seed provenance and are not
promoted to trusted reusable predictions. New legacy inference requires a fresh
`GBM_RUN_ROOT`, so historical caches and result tables cannot be reused or replaced
implicitly. `evidence/historical_cache_audit.json` records every inspected run.

Twelve fresh predictions examine BOTH, SOX2_ONLY, SOX9_ONLY and NEITHER using the
same recorded actuator sequences and seed set 0/1/2. These are structural-confidence
and endpoint-distance controls, not cellular activity or leakage assays. Inputs,
models and raw confidence files are recorded under `evidence/matched_controls/`.
In the three BOTH predictions, endpoint distances are 49.7, 30.4 and 42.1 Å;
none meets the historical <15 Å proximity cutoff. Global ipTM is 0.221–0.236,
and C-to-D actor-pair ipTM is about 0.141–0.144. These fresh predictions do not
independently reproduce the selected MD starting assembly. Across conditions,
global ipTM is not directly an activity metric; directional actor-pair confidence
is reported separately. A seed label does not imply identical
random noise for inputs with different shapes.

`evidence/external_artifacts.json` identifies all external production trajectories
used here by content hashes. Their data remain on D:, rather than in Git. The
manifest is not an off-machine data backup.

## Remaining scientific limits

- The current caspase-core models have substantial reference-fold disagreement.
  A structurally supported assembly and its relevant oligomeric context are needed
  before further production MD could strengthen a mechanism claim. No parameter
  choice in this analysis makes a failed structural observation a positive result.
- Presence/absence predictions and one unbound MD trajectory do not establish
  enzymatic leakiness, cellular AND logic or normal-cell specificity. These remain
  unsupported rather than passed by computational substitution.
- Delivery, immune activation and treatment readiness are not established by the
  available structure/trajectory data. The cited SSE study does not validate this
  protein actuator; see `notes/literature_scope_review.md`.
- Suitable independent functional labels are absent. Training an ML activity
  predictor on these heuristic scores would not independently validate activity.
- The 31 historical shell launchers are archived and inactive by default. Repaired
  entry points and explicit fresh-output protocols are the supported routes.

The branch records reproducible computations, repaired workflow behavior, failed
observations and evidence limits. It does not assert biological function or safety.
