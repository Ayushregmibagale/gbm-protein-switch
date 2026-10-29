#!/usr/bin/env python3
"""Audit recorded BindCraft identities and downstream table exposure, without inference."""
import argparse
import csv
import hashlib
import json
from pathlib import Path


def audit(root):
    from Bio.PDB import PDBParser
    from Bio.Data.PDBData import protein_letters_3to1
    parser = PDBParser(QUIET=True)
    def table(name):
        with (root/'selection'/name).open() as handle:
            return list(csv.DictReader(handle, delimiter='\t'))
    rows = table('stage1_results.tsv')
    checked = []
    for row in rows:
        if row['source'] != 'bindcraft':
            continue
        directory = {'SOX2':'sox2_binder_xtal','SOX9':'sox9_binder_xtal_v3'}[row['target']]
        path = Path('bindcraft/designs')/directory/'Accepted'/(row['name']+'.pdb')
        model = parser.get_structure(row['name'], root/path)[0]
        def sequence(chain):
            return ''.join(protein_letters_3to1[r.resname] for r in model[chain] if r.id[0] == ' ')
        binder, target = sequence('B'), sequence('A')
        checked.append({'name':row['name'],'target':row['target'],
                        'recorded_length':len(row['binder_seq']),'chain_B_length':len(binder),
                        'recorded_sequence_matches_chain_B':row['binder_seq']==binder,
                        'recorded_sequence_matches_chain_A':row['binder_seq']==target,
                        'source_path':str(path), 'source_sha256':hashlib.sha256((root/path).read_bytes()).hexdigest()})
    invalid = {row['name'] for row in checked if not row['recorded_sequence_matches_chain_B']}
    exposure = {}
    for name in ['stage2_results.tsv','stage3_results.tsv','stage4_results.tsv','stage5_results.tsv','stage5_flipped_results.tsv']:
        data = table(name)
        affected = [row for row in data if row.get('sox2_name',row.get('sox2')) in invalid or row.get('sox9_name',row.get('sox9')) in invalid]
        exposure[name] = {'rows':len(data),'rows_referencing_mismatched_candidates':len(affected)}
    return {'scope':'Historical table identity audit; affected prediction scores cannot be repaired by relabelling',
            'checked_bindcraft_rows':len(checked),'mismatched_rows':len(invalid),
            'checked':checked,'downstream_exposure':exposure,
            'stage1_sha256':hashlib.sha256((root/'selection/stage1_results.tsv').read_bytes()).hexdigest()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[1])
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--quarantine-output', type=Path, help='New table excluding mismatched candidates; no prediction scores are invented')
    args=parser.parse_args()
    report=audit(args.root)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open('w') as handle:
        json.dump(report,handle,indent=2);handle.write('\n')
    if args.quarantine_output:
        invalid = {row['name'] for row in report['checked'] if not row['recorded_sequence_matches_chain_B']}
        with (args.root/'selection/stage1_results.tsv').open() as handle:
            reader = csv.DictReader(handle, delimiter='\t')
            fields = reader.fieldnames
            retained = [row for row in reader if row['name'] not in invalid]
        args.quarantine_output.parent.mkdir(parents=True, exist_ok=True)
        with args.quarantine_output.open('x', newline='') as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t')
            writer.writeheader();writer.writerows(retained)
    print(f"{report['mismatched_rows']}/{report['checked_bindcraft_rows']} recorded BindCraft identities differ from chain B")


if __name__=='__main__':
    main()
