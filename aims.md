# GBM protein-switch — current computational aims

The project tests whether recorded SOX2/SOX9 binders and caspase-fragment models
support a conditional protein-sensing hypothesis. It is exploratory research;
there is no established therapeutic function or safety result.

The original proposal is preserved in [docs/historical_aims.md](docs/historical_aims.md).
Current conclusions come from [ROBUSTNESS_REVIEW.md](ROBUSTNESS_REVIEW.md) and the
versioned evidence files, rather than the proposal's anticipated outcomes.

## Computational questions

1. Are candidate identities, input sequences, structural models and trajectory
   topologies consistent? Repair mistaken binder extraction and use fresh,
   seeded predictions rather than relabelling historical scores.
2. Are contact and endpoint-distance observations reproducible across independent
   trajectory runs, periodic images, time windows and metric definitions?
3. Do assembled caspase cores resemble experimental reference conformations,
   rather than merely bringing selected termini close together?
4. Does the complete unbound control meet the stated distance heuristic after
   its interrupted trajectory is validated? This observation cannot by itself
   establish absence of enzymatic leakiness.
5. Do matched controls support a reproducible difference in model confidence?
   Confidence is not affinity, target specificity or cellular selectivity.
6. Are conclusions restricted to what the simulated constructs actually contain?
   The recorded simulations omit the original proposal's added NLS/T2A elements
   and use GGS1 linkers with a flipped second fusion.

## Decision rules

Preserve failed candidates, missing evidence and replicate variation. Treat
reference-fold disagreement as a mechanistic concern. Do not equate contact
occupancy, ipTM, or a short endpoint distance with functional AND-gate activity.
The reconstituted catalytic core and its functional oligomeric context require
separate support. Normal-cell specificity, delivery, immune payloads and biological
function remain outside what these structure/trajectory calculations establish.

ML ranking cannot be validated as an activity predictor without suitable outcome
labels and independent evaluation data. An ML score trained on these structural
heuristics would reproduce those assumptions rather than validate the mechanism.
