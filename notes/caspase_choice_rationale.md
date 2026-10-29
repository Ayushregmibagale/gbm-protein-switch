# Why Split-Caspase-9 (not Caspase-3/7)

## The caspase cascade

Caspase-9 (initiator) → activates → Caspase-3/7 (executioners) → cell death

## Rationale for caspase-9 as the AND gate

**1. Constitutive activity risk**
Caspase-3 and -7 are executioner caspases. Even partial/leaky reconstitution
of a split executioner caspase could kill the cell. Caspase-9 requires
dimerisation AND its CARD domain to function — more layers of control.

**2. Signal amplification**
One reconstituted caspase-9 dimer activates many caspase-3/7 molecules.
Signal amplification means fewer reconstitution events are needed to commit
the cell to death — lower threshold for tumour killing.

**3. Precedent and validated split site**
Split-caspase-9 is the most validated split-protease for coincidence detection
in synthetic biology (iCasp9 suicide gene work, AND-gate CAR-T papers).
The split site between the CARD domain and catalytic domain (~residue 316)
is well-characterised and gives tight off/on control.

**4. Reversibility window**
Caspase-9 activation can still be aborted early (e.g. by IAPs).
Caspase-3/7 activation is essentially the point of no return.
Having caspase-9 as the AND gate gives a small safety buffer against
transient spurious co-localisation in normal cells.

## When caspase-3 might be preferred

If guaranteed lethality with minimal reconstitution events is the priority,
split caspase-3 is more potent — but the leakiness risk in normal cells
makes it less suitable for a therapeutic context.

## Bottom line

Caspase-9: upstream, amplified, precedented, safer.
Caspase-3/7: downstream, potent, leaky risk — better for research tools
than clinical-grade AND gates.
