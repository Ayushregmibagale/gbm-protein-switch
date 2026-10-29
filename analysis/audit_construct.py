#!/usr/bin/env python3
"""Compare recorded construct inputs with the protein sequences in MD topologies."""
import argparse
import hashlib
import json
from pathlib import Path


def topology_sequences(path):
    import MDAnalysis as mda
    from Bio.Data.PDBData import protein_letters_3to1
    universe = mda.Universe(str(path))
    sequences = []
    for segment in universe.segments:
        residues = segment.atoms.select_atoms('protein').residues
        if not len(residues):
            continue
        try:
            sequence = ''.join(protein_letters_3to1[r.resname] for r in residues)
        except KeyError as exc:
            raise ValueError(f'Noncanonical residue requires explicit mapping: {exc}') from exc
        sequences.append({'segment': segment.segid, 'length': len(sequence),
                          'first_resid': int(residues[0].resid), 'last_resid': int(residues[-1].resid),
                          'sequence_sha256': hashlib.sha256(sequence.encode()).hexdigest(),
                          'sequence': sequence})
    return sequences


def compare(bound_input, unbound_input, bound_topology, unbound_topology, reference):
    import yaml
    def read(path):
        doc = yaml.safe_load(path.read_text())
        rows = doc['sequences']
        pairs = [(row['protein']['id'], row['protein']['sequence']) for row in rows]
        if len({key for key, _ in pairs}) != len(pairs):
            raise ValueError('Duplicate chain ID in recorded input')
        return dict(pairs)
    bound, unbound = read(bound_input), read(unbound_input)
    observed_bound, observed_unbound = topology_sequences(bound_topology), topology_sequences(unbound_topology)
    if [row['sequence'] for row in observed_bound] != list(bound.values()):
        raise ValueError('Bound topology protein sequences differ from recorded prediction input order')
    if [row['sequence'] for row in observed_unbound] != list(unbound.values()):
        raise ValueError('Unbound topology protein sequences differ from recorded prediction input order')
    if bound.get('C') != unbound.get('A') or bound.get('D') != unbound.get('B'):
        raise ValueError('Bound and unbound actuator chains differ')
    ref = ''.join(reference.read_text().splitlines()[1:])
    if len(ref) != 416 or ref[314:316] != 'DA':
        raise ValueError('Unexpected reference length or split-site identity')
    if not bound['C'].endswith(ref[139:315]) or not bound['D'].startswith(ref[315:416]):
        raise ValueError('Observed construct does not match reference caspase halves')
    return {
        'bound_actuator_matches_unbound': True,
        'topologies_match_recorded_inputs': True,
        'reference_length': len(ref), 'reference_315_316': ref[314:316],
        'observed_bound_segments': observed_bound, 'observed_unbound_segments': observed_unbound,
        'scope': 'Identity of existing simulation inputs only; no activity or safety inference',
        'inputs': [{'path': str(path.resolve()), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}
                   for path in [bound_input, unbound_input, bound_topology, unbound_topology, reference]],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    root = Path(__file__).resolve().parents[1]
    parser.add_argument('--bound-input', type=Path, default=root/'inputs/track1_bound.yaml')
    parser.add_argument('--unbound-input', type=Path, default=root/'inputs/track1_unbound.yaml')
    parser.add_argument('--bound-topology', type=Path, required=True)
    parser.add_argument('--unbound-topology', type=Path, required=True)
    parser.add_argument('--reference', type=Path, default=root/'inputs/reference/P55211.fasta')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    try:
        report = compare(args.bound_input, args.unbound_input, args.bound_topology, args.unbound_topology, args.reference)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        with args.output.open('x') as handle:
            json.dump(report, handle, indent=2)
            handle.write('\n')
    except (OSError, ValueError, KeyError) as exc:
        parser.exit(1, f'Construct identity audit failed: {exc}\n')
    print('Recorded inputs and bound/unbound topology sequences agree; reference split is D315/A316.')


if __name__ == '__main__':
    main()
