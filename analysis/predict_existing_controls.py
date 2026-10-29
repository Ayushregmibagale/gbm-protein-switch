#!/usr/bin/env python3
"""Matched-seed presence/absence predictions of the existing recorded actuator."""
import argparse
import sys
import csv
import json
import os
from pathlib import Path
import subprocess
import time
from audit_control import digest


def main():
    import yaml
    from Bio.PDB import PDBParser
    from Bio.Data.PDBData import protein_letters_3to1
    import numpy as np
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--boltz',type=Path,required=True)
    parser.add_argument('--checkpoint',type=Path,required=True)
    parser.add_argument('--cache',type=Path,required=True)
    parser.add_argument('--wait-manifest',type=Path)
    args=parser.parse_args();root=Path(__file__).resolve().parents[1]
    if args.output.exists():parser.error('Use a new directory')
    args.output.mkdir(parents=True);inputs=args.output/'inputs';inputs.mkdir()
    source=root/'inputs/track1_bound.yaml'
    proteins=yaml.safe_load(source.read_text())['sequences']
    sets={'BOTH':{'A','B','C','D'},'SOX2_ONLY':{'A','C','D'},'SOX9_ONLY':{'B','C','D'},'NEITHER':{'C','D'}}
    requested={}
    for name,ids in sets.items():
        rows=[row for row in proteins if row['protein']['id'] in ids]
        requested[name]={row['protein']['id']:row['protein']['sequence'] for row in rows}
        (inputs/f'{name}.yaml').write_text(yaml.safe_dump({'sequences':rows},sort_keys=False))
    manifest={'status':'prepared','source_sha256':digest(source),'checkpoint_sha256':digest(args.checkpoint),
              'seeds':[0,1,2],'recycling_steps':3,'sampling_steps':200,'model':'boltz2','commands':[],
              'scope':'Existing sequences only; presence/absence structure-confidence observations, not activity or cellular selectivity'}
    path=args.output/'manifest.json'
    def save():path.write_text(json.dumps(manifest,indent=2)+'\n')
    save()
    if args.wait_manifest:
        manifest['status']='waiting_for_candidate_batch';save();deadline=time.monotonic()+3600
        while json.loads(args.wait_manifest.read_text()).get('status')!='complete':
            if time.monotonic()>deadline:raise TimeoutError('Candidate batch did not complete')
            time.sleep(10)
    manifest['status']='running';save();records=[];pdbparser=PDBParser(QUIET=True)
    env=os.environ.copy();env['NUMBA_CACHE_DIR']=str(args.output/'numba_cache');env['MPLCONFIGDIR']=str(args.output/'matplotlib')
    for seed in (0,1,2):
        output=args.output/f'seed{seed}'
        cmd=[str(args.boltz),'predict',str(inputs),'--out_dir',str(output),'--cache',str(args.cache),
             '--checkpoint',str(args.checkpoint),'--model','boltz2','--seed',str(seed),'--recycling_steps','3',
             '--sampling_steps','200','--output_format','pdb','--no_kernels','--num_workers','0','--preprocessing-threads','1']
        manifest['commands'].append(cmd);save();print(f'Control predictions seed {seed}',flush=True)
        with (args.output/f'seed{seed}.log').open('w') as handle:
            subprocess.run(cmd,check=True,stdout=handle,stderr=subprocess.STDOUT,env=env)
        for name,expected in requested.items():
            files=list(output.rglob(f'confidence_{name}_model_0.json'))
            if len(files)!=1:raise ValueError('Missing/ambiguous control confidence')
            structure=files[0].parent/f'{name}_model_0.pdb';model=pdbparser.get_structure(name,structure)[0]
            observed={chain.id:''.join(protein_letters_3to1[r.resname] for r in chain if r.id[0]==' ') for chain in model}
            if observed!=expected:raise ValueError('Control structure differs from requested sequences')
            conf=json.loads(files[0].read_text());iptm=conf['iptm']
            if not np.isfinite(iptm) or not 0<=iptm<=1:raise ValueError('Invalid confidence')
            c=[r['CA'].coord for r in model['C'] if r.id[0]==' ' and 'CA' in r]
            d=[r['CA'].coord for r in model['D'] if r.id[0]==' ' and 'CA' in r]
            ids=list(expected);pair=conf['pair_chains_iptm']
            records.append({'condition':name,'seed':seed,'global_iptm':iptm,
                            'C_to_D_iptm':pair[str(ids.index('C'))][str(ids.index('D'))],
                            'D_to_C_iptm':pair[str(ids.index('D'))][str(ids.index('C'))],
                            'endpoint_CA_distance_A':float(np.linalg.norm(c[-1]-d[0])),
                            'confidence_path':str(files[0]),'structure_path':str(structure),
                            'confidence_sha256':digest(files[0]),'structure_sha256':digest(structure)})
        with (args.output/'results.tsv').open('w',newline='') as handle:
            writer=csv.DictWriter(handle,fieldnames=list(records[0]),delimiter='\t',lineterminator='\n');writer.writeheader();writer.writerows(records)
        manifest['completed_seeds']=seed+1;save()
    manifest['status']='complete';manifest['n_predictions']=len(records);save()
    print('Completed 12 matched control predictions',flush=True)


if __name__=='__main__':
    try:
        main()
    except Exception as exc:
        if '--output' in sys.argv:
            output=Path(sys.argv[sys.argv.index('--output')+1])/'manifest.json'
            if output.exists():
                state=json.loads(output.read_text())
                if state.get('status')!='complete':
                    state['status']='failed';state['failure']=str(exc)
                    output.write_text(json.dumps(state,indent=2)+'\n')
        raise
