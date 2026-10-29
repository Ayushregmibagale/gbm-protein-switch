# GBM protein-switch robustness review — 2026-10-29

Scope: this project only. Recovered sources were read from
`/mnt/d/gbm-protein-switch`; the WSL home checkout is the intended source-code
location. Large simulation data remain on D:. No new inference campaign or MD
production was launched during this review. The externally resumed unbound run
was left undisturbed. See [WORKFLOW.md](WORKFLOW.md) for operating instructions.

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
predictions after changing candidate identities. Corrected predictions for the
excluded binders have not been run.

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
of the source PDB/TPR/XTC files. Historical tables are preserved. Additional bridge
replicate trajectories have not been reanalysed. Frames are correlated; 2001 frames
are not 2001 independent experiments. Contact occupancy does not establish AND logic.

### The unbound control remains unfinished

The state file reads `NPT_DONE`, while the resumed production log and trajectory
are actively advancing. No completed control metrics were verified. The old
trajectory has valid frames through 5.47 ns and an incomplete frame at 5.48 ns.
`gmx check` reports this warning but returns zero. The resumed run begins from a
checkpoint at 4.5726 ns and overlaps the surviving old trajectory.

The legacy completion launcher ignores concatenation failure and can substitute
only the resumed portion. Its output must not be accepted without coverage and
integrity checks. The new `analysis/audit_control.py` checks exact frame times,
gaps, overlaps, corruption in the retained portion, atom counts, periodic boxes,
input stability and full expected coverage. It writes no completed metrics for
partial data. `WORKFLOW.md` documents an explicit, provenance-recorded splice at
4.57/4.58 ns after the writer finishes. A snapshot is recorded in
`evidence/unbound_integrity_status.json`.

There is currently no verified “non-leaky” result. Even meeting the >25 Å/95%
distance heuristic would not establish absence of enzymatic leakiness.

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

## Remaining evidence requirements

1. Let the existing unbound writer finish, then validate full coverage and compute
   corrected control metrics. No partial trajectory should be substituted.
2. Review individual-input and target-absent controls, productive catalytic geometry,
   starting-state bias, replicate dependence and trajectory-window sensitivity before
   claiming AND-gate activity. The current observations cannot answer those questions.
3. Re-run excluded Stage 1 candidates only if needed, using the corrected identities,
   normalized confidence scale, new output directories and recorded provenance.
4. Archive complete external data durably. The new bridge results have full raw-input
   content hashes; older Track 1 periodic exports identify trajectories by size/mtime
   plus topology hashes and still need complete trajectory archival manifests.
5. Review primary literature and experimental evidence for mechanism, normal-cell
   specificity, delivery and immune-payload claims. These were not independently
   established by this code/data review. Historical shell launchers and inference
   caches are preserved records, not fully validated campaign infrastructure.

The repository now has recovered sources, stronger execution checks and auditable
analysis. Biological function, non-leakiness and treatment readiness remain unproven.
