# GBM protein switch

Exploratory computational research on SOX2/SOX9 sensing and a split-caspase actuator.
The repository contains design outputs and analysis scripts. Conditional activity,
normal-cell specificity, and safety remain unestablished.

Start with [the robustness review](ROBUSTNESS_REVIEW.md). [aims.md](aims.md) describes
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

Install `requirements-audit.txt` in an isolated environment. Raw TPR/XTC files are
external to this checkout. The following example reads them without launching MD:

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
the topology. Historical bridge tables are preserved. A separate periodic reanalysis of
bridge_48_s01 is saved under `evidence/bridge48_periodic/`; remaining bridge
designs have not been reanalysed.

The three-chain prediction launcher accepts `--base` and `--boltz` (or `BOLTZ_BIN`).
It can launch expensive inference jobs; it is separate from the read-only audit.
Remaining campaign launchers depend on workstation paths and external software,
and are not a portable workflow. Do not interpret them as a documented entry point.
