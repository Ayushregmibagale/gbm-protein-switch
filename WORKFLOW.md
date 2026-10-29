# Operating the recovered workflow

The WSL checkout is the source-code location. Large existing simulation data remain
under `/mnt/d/gbm-protein-switch`. Set `GBM_ROOT` to that data root when inspecting
or deliberately resuming the old campaigns. Set `GMX_BIN` and `BOLTZ_BIN` to the
executables in your own environments. No environment paths are inferred from a
particular home directory in the repaired MD entry points.

## Status and analysis

These commands only report status; no `--execute` means no MD launch:

```bash
GBM_ROOT=/mnt/d/gbm-protein-switch python md/run_lead_md.py --status
GBM_ROOT=/mnt/d/gbm-protein-switch python md/run_replicates.py --status
GBM_ROOT=/mnt/d/gbm-protein-switch python md/run_unbound_md.py --status
```

State files are checked against required nonempty stage outputs. Historical
`ANALYSIS_DONE` and legacy `pass` fields still need the periodic/coverage audit;
file presence alone is not scientific validation.

`analysis/audit_control.py` checks topology atom counts, unique endpoint selections,
periodic boxes, exact frame cadence, start/end coverage, ordered fragments, and
input immutability during the read. It hashes complete source files and writes
to a new output directory only after validation. A failed/partial input produces
no completed metrics. Files alongside the original trajectories are not created.

## Interrupted Track 1 unbound control

On 2026-10-29, the existing resumed run was actively writing `prod.part0002.xtc`.
The old `prod_partial_corrupted.xtc` has valid data through 5.47 ns and an incomplete
frame at 5.48 ns. `gmx check` reports that truncation but returns exit code zero:
return-code checking alone cannot establish trajectory integrity.

The resumed log records continuation from step 2286300 (4.5726 ns). Therefore a
deliberate splice can retain old frames through 4.57 ns and resumed frames from
4.58 ns onward. The audit must confirm those exact frame times and final coverage.
Do not concatenate the overlapping files blindly or substitute only the resumed
part. Do not run this against a file that is still being written.

The completed fragments were validated with the following command (choose a new
output directory for another audit):

```bash
python analysis/audit_control.py \
  --topology /mnt/d/gbm-protein-switch/md/unbound/Track1_unbound/prod.tpr \
  --fragments /mnt/d/gbm-protein-switch/md/unbound/Track1_unbound/prod_partial_corrupted.xtc \
              /mnt/d/gbm-protein-switch/md/unbound/Track1_unbound/prod.part0002.xtc \
  --first-fragment-end-ns 4.57 \
  --first 'resid 250 and name CA' --second 'resid 251 and name CA' \
  --label Track1_unbound --criterion far \
  --output /tmp/gbm-unbound-validated
```

This explicitly discards the old post-checkpoint tail. The original files remain
unchanged. The command fails for gaps, duplicate times, corruption in the retained
portion, ongoing writes, or missing 0–20 ns coverage. Review its topology hashes,
endpoint identities and fragment provenance before interpretation. Even a met
distance criterion is not evidence of absent enzymatic leakiness.

## Deliberate execution

The three recovered MD launchers require `--execute` to run stages. GROMACS is
invoked with argument lists and checked return codes. Stage transitions require
their outputs and are written atomically. Production resumes from `prod.cpt` and
the existing `prod.tpr`, with `-noappend`; it refuses partial outputs without a
checkpoint. Fragmented output requires the explicit validation above before a
complete production state is recorded. Run only one writer per campaign. The
existing external unbound writer was left undisturbed during this work.

`md/watch_and_analyse.sh` now forwards explicit inputs to the coverage audit; it
never infers completion from a PID or writes a production state.

Recovered selection campaigns require `--execute` as well, except Stage 1's
`--dry_run`. Their historical scoring, candidate selection, and inference caches
still require scientific review. The legacy shell campaign launchers remain
workstation records, not validated entry points. Do not run them as a substitute
for the repaired MD commands.

The repaired Stage 4 ranking is read-only unless given a new output path:

```bash
python selection/stage4_composite_ranking.py \
  --input selection/stage3_results.tsv --output /tmp/stage4-reproduced.tsv
```

It reproduces all ten historical ranks, rejects zero or nonfinite denominators,
and treats the score as a model-confidence heuristic rather than specificity.


## Candidate identity quarantine

The recorded Stage 1 table has 23 BindCraft rows containing the target sequence
instead of the binder sequence. The historical table is preserved; new inference
stages reject those identities. A filtered table is supplied without inventing
replacement scores:

```bash
export GBM_STAGE1_TABLE=/path/to/checkout/selection/validated_stage1_results.tsv
```

Regenerate the identity report and filtered table, if needed, with:

```bash
python analysis/audit_candidate_identity.py --output /tmp/candidate-identities.json \
  --quarantine-output /tmp/stage1-without-mismatches.tsv
```

New Stage 1 inference requires `--execute --run-dir NEW_DIRECTORY --results NEW_TABLE`.
Its `--dry_run` only collects designs and creates no output directories. Changing
candidate sequences requires new predictions; renaming existing scores is invalid.


## Fresh prediction campaigns and archived launchers

Supported prediction commands require fresh output directories. Legacy Python
campaigns additionally require a new `GBM_RUN_ROOT`; their results and predictions
are redirected there, preserving historical tables and caches. Existing directories
are rejected. The 31 historical shell launchers are inactive unless
`GBM_ALLOW_LEGACY_CAMPAIGNS=1` is deliberately set. They remain archived workflows,
not supported alternatives to the repaired commands.

The corrected-candidate and matched-control prediction protocols are recorded in
`analysis/repair_candidate_predictions.py` and `analysis/predict_existing_controls.py`.
Their manifests record exact commands, model checkpoint identity, seeds and
protocol settings. Predictions use existing sequences, empty MSAs, three recycling
steps and 200 sampling steps. Fresh outputs are checked against requested sequence
identities before scores are accepted. The evidence exports include inputs,
confidence files, PDB models and result tables. Inference runtimes are separate
from the analysis environment; see `evidence/prediction_environment.json` and the
batch manifests.

`selection/repaired_stage1_candidates.tsv` contains only the 23 re-predicted
candidates, with mean confidence, between-seed SD and an all-seed confidence
heuristic. It is not silently merged with historical scores or treated as an
affinity/activity ranking.
