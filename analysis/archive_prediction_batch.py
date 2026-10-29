#!/usr/bin/env python3
"""Validate and export a completed prediction batch without weight/cache dumps."""
import argparse
import csv
import json
from pathlib import Path
import shutil
import statistics
import tempfile
from audit_control import digest


def archive(source,destination,kind):
    import numpy as np
    from Bio.PDB import PDBParser
    manifest=json.loads((source/'manifest.json').read_text())
    if manifest.get('status')!='complete':raise ValueError('Only completed batches can be archived')
    if destination.exists():raise FileExistsError('Choose a new evidence directory')
    with (source/'results.tsv').open() as handle:rows=list(csv.DictReader(handle,delimiter='\t'))
    expected=69 if kind=='candidates' else 12
    group_key='name' if kind=='candidates' else 'condition'
    if len(rows)!=expected:raise ValueError('Incomplete result table')
    for group in {row[group_key] for row in rows}:
        seeds=[row['seed'] for row in rows if row[group_key]==group]
        if sorted(seeds)!=['0','1','2']:raise ValueError('Incomplete or duplicated seed coverage')
    destination.parent.mkdir(parents=True,exist_ok=True)
    with tempfile.TemporaryDirectory(dir=destination.parent,prefix='.prediction-export-') as temporary:
        out=Path(temporary)/'export';out.mkdir()
        shutil.copytree(source/'inputs',out/'inputs')
        parser=PDBParser(QUIET=True)
        for row in rows:
            for path_key,hash_key in [('confidence_path','confidence_sha256'),('structure_path','structure_sha256')]:
                original=Path(row[path_key]);target=out/f"seed{row['seed']}"/original.name
                target.parent.mkdir(exist_ok=True);shutil.copy2(original,target)
                if digest(target)!=row[hash_key]:raise ValueError('Prediction changed or its hash disagrees')
                row[path_key]=str(target.relative_to(out))
            model=parser.get_structure('model',out/row['structure_path'])
            if not np.isfinite(np.array([atom.coord for atom in model.get_atoms()])).all():
                raise ValueError('Nonfinite predicted atom coordinates')
        with (out/'results.tsv').open('w',newline='') as handle:
            writer=csv.DictWriter(handle,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n');writer.writeheader();writer.writerows(rows)
        summary=[]
        for group in sorted({r[group_key] for r in rows}):
            selected=[r for r in rows if r[group_key]==group]
            fields=['iptm','plddt_fraction'] if kind=='candidates' else ['global_iptm','C_to_D_iptm','D_to_C_iptm','endpoint_CA_distance_A']
            entry={group_key:group,'n_seeds':len(selected)}
            for field in fields:
                values=[float(row[field]) for row in selected]
                if not all(np.isfinite(values)):raise ValueError('Nonfinite summary data')
                entry[field]={'values':values,'mean':statistics.mean(values),'sd':statistics.stdev(values),'min':min(values),'max':max(values)}
            if kind=='candidates':
                entry['target']=selected[0]['target']
                entry['n_seeds_meeting_historical_confidence_cutoffs']=sum(float(r['iptm'])>=.7 and float(r['plddt_fraction'])>=.7 for r in selected)
            else:entry['n_predicted_endpoints_below_15A']=sum(float(r['endpoint_CA_distance_A'])<15 for r in selected)
            summary.append(entry)
        (out/'summary.json').write_text(json.dumps({'scope':'Descriptive model confidence/geometry, not affinity, activity or safety','groups':summary},indent=2)+'\n')
        for path in source.glob('seed*.log'):shutil.copy2(path,out/path.name)
        shutil.copy2(source/'manifest.json',out/'manifest.json')
        (out/'export_manifest.json').write_text(json.dumps({'kind':kind,'n_predictions':len(rows),
            'files_sha256':{str(p.relative_to(out)):digest(p) for p in sorted(out.rglob('*')) if p.is_file()}},indent=2)+'\n')
        out.rename(destination)
    return summary


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source',type=Path,required=True);parser.add_argument('--destination',type=Path,required=True)
    parser.add_argument('--kind',choices=['candidates','controls'],required=True)
    args=parser.parse_args();archive(args.source,args.destination,args.kind)
    print(f'Validated and archived {args.kind}: {args.destination}')


if __name__=='__main__':main()
