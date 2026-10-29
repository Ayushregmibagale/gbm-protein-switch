#!/usr/bin/env python3
"""Re-predict recorded misidentified candidates in isolated, seeded Boltz batches."""
import argparse
import csv
import json
import os
from pathlib import Path
import subprocess
import sys
import time
from audit_control import digest


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[1])
    parser.add_argument('--output',required=True,type=Path)
    parser.add_argument('--boltz',required=True,type=Path)
    parser.add_argument('--checkpoint',required=True,type=Path)
    parser.add_argument('--cache',required=True,type=Path)
    parser.add_argument('--wait-for-log',type=Path)
    parser.add_argument('--execute',action='store_true')
    args=parser.parse_args()
    import yaml
    from Bio.PDB import PDBParser
    from Bio.Data.PDBData import protein_letters_3to1
    report=json.loads((args.root/'evidence/stage1_chain_identity.json').read_text())
    selected=[row for row in report['checked'] if not row['recorded_sequence_matches_chain_B']]
    if args.output.exists():parser.error('Use a new output directory')
    args.output.mkdir(parents=True)
    inputs=args.output/'inputs';inputs.mkdir()
    # Use the recorded target definitions, without altering candidate design sequences.
    sys.path.insert(0,str(args.root/'selection'))
    from stage1_reprediction import SOX2_TARGET,SOX9_TARGET
    pdbparser=PDBParser(QUIET=True)
    records=[]
    for row in selected:
        source=args.root/row['source_path']
        if digest(source)!=row['source_sha256']:raise ValueError('Accepted source structure changed')
        model=pdbparser.get_structure(row['name'],source)[0]
        sequence=''.join(protein_letters_3to1[r.resname] for r in model['B'] if r.id[0]==' ')
        target=SOX2_TARGET if row['target']=='SOX2' else SOX9_TARGET
        path=inputs/(row['name']+'.yaml')
        path.write_text(yaml.safe_dump({'sequences':[{'protein':{'id':'A','sequence':sequence,'msa':'empty'}},
                                                     {'protein':{'id':'B','sequence':target,'msa':'empty'}}]},sort_keys=False))
        records.append({'name':row['name'],'target':row['target'],'binder_length':len(sequence),
                        'input_sha256':digest(path),'binder_sequence':sequence,'target_sequence':target,'accepted_pdb_sha256':digest(source)})
    manifest={'status':'prepared','checkpoint':str(args.checkpoint),'checkpoint_sha256':digest(args.checkpoint),
              'model':'boltz2','seeds':[0,1,2],'recycling_steps':3,'sampling_steps':200,'candidates':records,
              'scope':'Corrections to recorded binder identity, not a new design campaign'}
    manifest_path=args.output/'manifest.json'
    def save():manifest_path.write_text(json.dumps(manifest,indent=2)+'\n')
    save()
    if not args.execute:return
    if args.wait_for_log:
        manifest['status']='waiting_for_existing_md';save()
        deadline=time.monotonic()+7200
        while True:
            if args.wait_for_log.exists():
                with args.wait_for_log.open('rb') as handle:
                    handle.seek(max(0,args.wait_for_log.stat().st_size-65536));tail=handle.read()
                if b'Finished mdrun' in tail:break
                if b'Fatal error' in tail:raise RuntimeError('Existing MD failed; inspect it before prediction')
            if time.monotonic()>deadline:raise TimeoutError('MD wait exceeded two hours')
            time.sleep(15)
    python=args.boltz.parent/'python'
    versions=subprocess.check_output([str(python),'-c',"import json,importlib.metadata as m; print(json.dumps({p:m.version(p) for p in ['boltz','torch','numpy','pytorch-lightning']}))"],text=True)
    manifest['software']=json.loads(versions)
    manifest['status']='running';manifest['commands']=[];save()
    rows=[]
    env=os.environ.copy();env['NUMBA_CACHE_DIR']=str(args.output/'numba_cache');env['MPLCONFIGDIR']=str(args.output/'matplotlib')
    for seed in (0,1,2):
        output=args.output/f'seed{seed}'
        cmd=[str(args.boltz),'predict',str(inputs),'--out_dir',str(output),'--cache',str(args.cache),
             '--checkpoint',str(args.checkpoint),'--model','boltz2','--recycling_steps','3','--sampling_steps','200',
             '--seed',str(seed),'--output_format','pdb','--no_kernels','--num_workers','0','--preprocessing-threads','1']
        manifest['commands'].append(cmd);save()
        print(f'Running seed {seed}, {len(records)} corrected candidates',flush=True)
        with (args.output/f'seed{seed}.log').open('w') as handle:
            subprocess.run(cmd,check=True,stdout=handle,stderr=subprocess.STDOUT,env=env)
        for record in records:
            files=list(output.rglob(f"confidence_{record['name']}_model_0.json"))
            if len(files)!=1:raise RuntimeError(f"Missing/ambiguous prediction for {record['name']} seed {seed}")
            conf=json.loads(files[0].read_text());iptm=conf.get('iptm',conf.get('protein_iptm'))
            plddt=conf.get('complex_plddt',conf.get('plddt'))
            import math
            if not all(isinstance(x,(int,float)) and math.isfinite(x) and 0<=x<=1 for x in [iptm,plddt]):
                raise ValueError('Invalid confidence scores')
            structures=list(files[0].parent.glob(f"{record['name']}_model_0.pdb"))
            if len(structures)!=1:raise RuntimeError('Missing/ambiguous matched structure')
            model=pdbparser.get_structure(record['name'],structures[0])[0]
            observed={chain.id:''.join(protein_letters_3to1[r.resname] for r in chain if r.id[0]==' ') for chain in model}
            if observed!={'A':record['binder_sequence'],'B':record['target_sequence']}:
                raise ValueError('Prediction structure sequences differ from requested corrected input')
            rows.append({'name':record['name'],'target':record['target'],'seed':seed,'iptm':iptm,
                         'plddt_fraction':plddt,'confidence_sha256':digest(files[0]),'confidence_path':str(files[0]),'structure_path':str(structures[0]),'structure_sha256':digest(structures[0])})
        with (args.output/'results.tsv').open('w',newline='') as handle:
            writer=csv.DictWriter(handle,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n');writer.writeheader();writer.writerows(rows)
        manifest['completed_seeds']=seed+1;save()
    manifest['status']='complete';manifest['n_predictions']=len(rows);save()
    print('Completed 69 identity-corrected predictions',flush=True)


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
