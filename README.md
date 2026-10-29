# GBM protein switch

Exploratory computational research on SOX2/SOX9 sensing and a split-caspase actuator.
The repository contains design outputs, repaired workflows and versioned evidence. Conditional activity,
normal-cell specificity, and safety remain unestablished.

Start with [the robustness review](ROBUSTNESS_REVIEW.md) and the [operating instructions](WORKFLOW.md). [aims.md](aims.md) describes
the original proposal; its numerical criteria are hypotheses and historical choices.

## Reproduce the evidence audit

Python 3.11 or later; the recorded-table audit needs only the standard library:

```bash
python analysis/audit_evidence.py --output /tmp/recorded_evidence.json
diff evidence/recorded_evidence.json /tmp/recorded_evidence.json
python -m unittest discover -s tests -v
```

The audit checks table/time-series agreement, complete seed coverage, valid scores,
and input hashes. Successful validation means the recorded data are internally
consistent; it does not validate the biological mechanism.

## Reanalyse existing trajectories

Install `requirements-audit.txt` in an isolated environment. Raw TPR/XTC files remain external to this checkout. Missing small workflow files
have been recovered with hashes in `evidence/recovery_manifest.json`. The following example reads them without launching MD:

```bash
python analysis/audit_periodic_distance.py /path/to/md/lead_construct \
  --first 'resnum 419 and name CA' --second 'resnum 420 and name CA' \
  --output /tmp/lead_periodic.tsv
```

Selections are specific to the audited topology: verify atom identities before
using them for another construct. Output includes raw and minimum-image distances,
actual trajectory times, topology hashes, and software versions. Trajectory size
and modification time are recorded, but are not a trajectory content hash.

Legacy structure analyses now require `--pdb /path/to/model.pdb` and accept
`--out-dir /path/to/new/results`. `bridge_dual_contact.py` requires `--out-dir`
and accepts `--md-root`; it applies periodic boundaries and checks CA counts.
It retains the historical 97/72 target-residue assumptions: verify these against
the topology. Historical bridge tables are preserved. Periodic reanalysis of all five recorded
passing bridge designs is saved under `evidence/bridge_periodic/`, with raw input
content hashes. Additional bridge replicate trajectories have not been reanalysed.

The three-chain prediction launcher accepts `--base` and `--boltz` (or `BOLTZ_BIN`).
It can launch expensive inference jobs; it is separate from the read-only audit.
Recovered MD runners use `GBM_ROOT`, `GMX_BIN`, and `BOLTZ_BIN`, default to status,
and require `--execute` to run. See `WORKFLOW.md` for checkpoint handling and
strict control validation. Historical shell launchers remain workstation records.


## Recorded computational closure

The review branch includes the validated 20 ns unbound control, all 15 available
bridge trajectory reanalyses, reference-fold diagnostics, a 972-run cache audit,
69 corrected-binder predictions and 12 matched presence/absence predictions.
Consult `COMPUTATIONAL_STATUS.json` for outcomes and evidence paths. Passing code
checks do not make the structural observations positive: bridge_90 has substantial
replicate instability, and the current caspase-core models disagree with the
experimental reference fold.
