# Evidence and robustness review — 2026-10-29

Scope: gbm-protein-switch only. GitHub checkout plus read-only inspection of the
working data at `/mnt/d/gbm-protein-switch`. Neither RIPTAC nor ANDGate has been
reviewed. No new designs, inference campaigns, or MD production were run.

## Findings that change interpretation

1. **The reported Track 1 proximity occupancy contains a periodic-boundary artifact.**
   The original code uses direct coordinate subtraction. Reanalysis of the same
   three trajectories with minimum-image distances gives:

   | Replicate | Historical raw fraction <15 Å | Periodic fraction <15 Å | Elapsed time |
   |---|---:|---:|---:|
   | lead_construct | 76.262% | 99.950% | 20.000 ns |
   | rep2 | 80.110% | 100.000% | 20.000 ns |
   | rep3 | 87.756% | 100.000% | 20.000 ns |

   Equal-weight replicate averages are 81.376% raw and 99.983% periodic. The
   original raw fractions reproduce. The corrected observations establish only
   proximity of the selected endpoint CA atoms in these trajectories. They do not
   establish productive enzyme assembly or conditional apoptosis. Minimum-image
   distances describe the closest periodic images; catalytic geometry still
   requires whole-complex reconstruction and an appropriate structural assessment.
   Each trajectory has 2001 frames spanning 20 ns; the old 20.01 ns duration counted
   frames rather than intervals. See `evidence/track1_*_periodic.{tsv,json}`.

2. **The unbound-control completion claim is unsupported by the files inspected.**
   `md/unbound/Track1_unbound/state.txt` reads `NPT_DONE`; no completed
   `metrics.json` was found there. A resume/analysis launcher refers to a partial
   trajectory named `prod_partial_corrupted.xtc`, tolerates concatenation failure,
   and can fall back to only the resumed segment. A completion flag alone would
   therefore be insufficient: continuity, total analysed time, and checkpoint
   provenance must be checked before interpreting this control. No running jobs
   or original files were changed. Absence of completed evidence is not a measured
   leakiness result. The >25 Å / 95% criterion itself is a distance heuristic.

3. **The three-chain screen uses the wrong baseline in its historical verdicts.**
   Recorded on-target seeds are 0.2491, 0.6469, and 0.1482: mean **0.3480**, sample
   SD approximately **0.2636**, rather than the hard-coded 0.471. For every
   replacement condition, baseline-minus-replacement changes sign across seeds.
   These three seeds do not support a consistent direction of confidence change.
   The new report lists matched seed differences. The launcher now uses the
   measured on-target mean and refuses a comparison when a seed fails.
   ipTM thresholds are confidence heuristics, not measured specificity or affinity.

4. **Track 1 cannot be reconstructed from the current GitHub checkout.**
   The D: working copy has selection scripts, MD scripts, MDP files, structures,
   and trajectories missing from GitHub. The Linux working copy contains empty
   directories for several of these components. Essential code and small inputs
   need deliberate recovery and versioning; raw trajectories need an external
   manifest and durable storage. Existing launchers contain missing script paths.
   Around 24,000 files are tracked, largely generated structure outputs, while
   core workflow sources are missing. The ignore patterns also miss actual
   `boltzgen/output/` and BindCraft `Trajectory/` output names.

5. **The checked-in bridge contact summaries reproduce from rounded time series.**
   All five designs match the means and four contact states at thresholds 0.05,
   0.10, and 0.15 within their published rounding. This is an internal-consistency
   check, not a raw-trajectory validation. Historical contact code also omitted
   periodic boundaries. The updated code fixes this; a single-design check of bridge_48_s01 gives 99.6% dual occupancy at 0.10,
   versus the historical 99.4% (see `evidence/bridge48_periodic/`). Other designs
   and bridge replicates have not been reanalysed. Frames are temporally correlated;
   2001 frames are not 2001 independent experiments.

## Repairs in this branch

- Deterministic standard-library evidence audit with input SHA-256 hashes, strict
  seed completeness, duplicate/nonfinite rejection, and table/time-series checks.
- Periodic endpoint reanalysis with exact atom identities, actual times, and saved
  results for all three existing Track 1 bound trajectories.
- Bridge contact analysis now uses periodic distances, fails on missing inputs or
  invalid selections, counts chain-C residues independently of numbering, and
  requires an explicit output directory to avoid accidentally replacing history.
- Fixed duplicate first epitope patch in `01_analyse_cocomplex.py`.
- Fixed epitope residue/coordinate pairing in `02_epitope_geometry.py`; missing
  residues now fail clearly. Its terminal distances are labelled as target termini,
  and centroid separation is no longer presented as proof of reconstitution.
- Explicit structure paths replace workstation paths in both structure scripts;
  imports no longer create outputs or run analyses.
- Regression checks and GitHub Actions for the recorded evidence audit. Optional
  trajectory checks cover orthogonal and triclinic periodic boxes and missing boxes.

## Outstanding work before stronger claims

1. Recover and review the missing Track 1 sources and small inputs, with provenance.
   Archive complete trajectories and their content hashes separately from Git.
2. Establish the actual status and integrity of the unbound control. Review matched
   target-present, target-absent, and individual-input controls before claiming AND
   logic. Structural confidence and endpoint distances alone cannot establish this.
3. Reconcile construct identities: proposal uses GGS3/NLS language while the local
   unbound script assembles GGS linkers and different fusion orientation. The audited
   split endpoint atoms are ASP and ALA; verify the proposal's D315/G316 wording
   against the exact reference sequence and construct. Treat these as unresolved
   identity mismatches, not as grounds to alter sequences automatically.
4. Reanalyse bridge trajectories with periodic boundaries; assess independent
   replicates and sensitivity to time window and contact definitions. Preserve
   candidate-selection history, including failures, rather than only passing designs.
5. Review primary literature and experimental evidence for mechanistic, normal-cell,
   delivery, and immune-payload claims. Those claims were not independently verified
   in this code/data review. Nothing here establishes treatment readiness.

The current evidence supports a computational proximity observation and a set of
testable hypotheses. It does not yet support “computationally complete,” “non-leaky,”
or “ready for lab” as verified conclusions.
