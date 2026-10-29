# Recorded input identity

These files preserve inputs of simulations already performed. They are not a new
construct design or an experimentally validated therapeutic sequence.

`track1_bound.yaml` records the prediction used for the bound lead. Its actuator
chains are C (250 aa) and D (171 aa). `track1_unbound.yaml` records the corresponding
two-chain control, with identical actuator sequences as A and B.

The observed layout uses one GGS repeat with a flipped second fusion:

- Bound C / unbound A: recorded SOX2 binder, GGS, reference residues 140–315.
- Bound D / unbound B: reference residues 316–416, GGS, recorded SOX9 binder.

The simulated actuator chains do not contain the proposal's added NLS, GGS3 linker,
or T2A cassette. Claims about that intended full construct cannot be transferred
automatically from these simulations.

The supplied reference `reference/P55211.fasta` is 416 residues long and has **D315
and A316**, correcting the historical D315/G316 label. This is a sequence-identity
check against the recovered file, not independent literature verification of a
validated split site. `evidence/construct_identity.json` compares both topologies
with both recorded prediction inputs and records their hashes.
