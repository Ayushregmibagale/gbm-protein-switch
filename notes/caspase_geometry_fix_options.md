# Caspase Geometry — Fix Options (Sep 2 2026)

## Problem
Stage 5 Boltz-2 full assembly predictions: 0/20 passed 15 Å inter-half threshold.
Closest: BC mpnn17 × ph_sox9_run30, GGS×3 at 21.0 Å.

## Why this is expected (not a failure)
Boltz-2 is a static snapshot — it has no driving force to pull the caspase halves together.
The 21 Å gap is the resting/starting geometry. MD simulation then lets the linkers flex
and sample conformations, driving the halves toward each other over 20–50 ns.
The 15 Å threshold is the activation criterion for reconstitution, checked AFTER MD.

## Fix options (in order of priority)

### Option 4 — Diagnose terminus orientations (DO FIRST, 5 min)
Look at the Stage 2 PDB structures for the top 5 pairs. Measure:
- C-terminus of chain C (SOX2 binder) position
- C-terminus of chain D (SOX9 binder) position
- Direction the termini point relative to each other and the SOX surface

If C-termini point TOWARD each other → longer linkers will close the gap (use Option 1).
If C-termini point AWAY from each other → need orientation flip (use Option 3).
This is a 5-minute diagnosis that determines which fix to apply.

### Option 1 — Longer linkers (30 min, conditional on Option 4)
Only useful if termini point toward each other.
Try GGS×5, ×6, ×7, ×8 (15–24 aa linkers) for the best pair (BC mpnn17 × ph_sox9_run30).
GGS×3 gave 21 Å — a few more residues of reach may close the gap.

### Option 2 — Multiple Boltz-2 seeds (1 hr, can run in parallel)
Run 3–5 seeds per linker for the top 2 candidates.
The 21 Å result is one sample from a stochastic model — another seed may find
a lower-energy conformation where the halves are already at 14 Å.
Command: add --num_seeds 3 to boltz predict.

### Option 3 — Flip fusion orientation (conditional on Option 4)
If termini point away from each other, the caspase halves will always face away.
Fix: change one or both fusions from C-terminal to N-terminal:
  Current:  BinderSOX2-Cterm → linker → CaspN
  Flipped:  CaspN → linker → BinderSOX2-Nterm (N-terminal fusion)
This requires rewriting the chain C/D sequences in the YAML (swap binder+caspase order).
Feasibility depends on whether the N-terminus of each binder is accessible (not buried).

### Option 5 — Redesign with terminus proximity constraint (days, last resort)
Go back to BindCraft/PH and add a constraint that the binder's C-terminus is near
the other binder's terminus during design. This is the clean solution but requires
launching a new design campaign.

## What MD adds (after these fixes)
MD is validation, not a fix. It tells you whether the static geometry is
compatible with reconstitution — i.e., whether the halves can find each other
dynamically given the binder-imposed geometry. Even 21 Å can reconstitute if
the linkers are flexible enough and the halves are in the right hemisphere.

## Action plan
1. Option 4 → diagnosis (running now)
2. Option 1 → longer linkers if termini face each other
3. Option 2 → multiple seeds in parallel
4. Option 3 → flip orientation if termini face away
5. Best geometry → GROMACS MD 20–50 ns
