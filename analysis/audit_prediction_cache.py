#!/usr/bin/env python3
"""Inventory historical cache identity/integrity; never promote it to trusted provenance."""
import argparse
import json
import math
from pathlib import Path
from audit_control import digest


def audit(root):
    import yaml
    from Bio.PDB import PDBParser
    from Bio.Data.PDBData import protein_letters_3to1
    parser=PDBParser(QUIET=True)
    parents=[root/'selection'/name for name in ['stage1_boltz_outputs','stage2_boltz_outputs','stage3_boltz_outputs',
             'stage5_caspase_outputs','stage5_flipped_outputs','stage5_multiseed_outputs','stage5_ph_sox2_run18_outputs','bridge_linker_outputs']]
    parents += [root/'analysis'/name for name in ['bridge_3chain_selectivity_outputs','offtarget_selectivity_outputs']]
    rows=[]
    for parent in parents:
        if not parent.exists():continue
        for directory in sorted(p for p in parent.iterdir() if p.is_dir()):
            input_path=directory/'input.yaml'
            if not input_path.exists():continue
            row={'run':str(directory.relative_to(root)),'issues':[],'input_sha256':digest(input_path)}
            try:
                inputs=yaml.safe_load(input_path.read_text())['sequences']
                expected={item['protein']['id']:item['protein']['sequence'] for item in inputs}
                files=sorted(directory.rglob('confidence_*model_0.json'))
                if len(files)!=1:raise ValueError(f'Expected one confidence model_0, found {len(files)}')
                conf=json.loads(files[0].read_text());score=conf.get('iptm',conf.get('protein_iptm'))
                if not isinstance(score,(int,float)) or not math.isfinite(score) or not 0<=score<=1:
                    raise ValueError('Invalid/missing ipTM')
                structures=sorted(files[0].parent.glob('*_model_0.pdb'))
                if len(structures)!=1:raise ValueError(f'Expected one matched structure, found {len(structures)}')
                model=parser.get_structure(directory.name,structures[0])[0]
                observed={chain.id:''.join(protein_letters_3to1[r.resname] for r in chain if r.id[0]==' ') for chain in model}
                if expected!=observed:row['issues'].append('PDB protein chain sequences differ from input YAML')
                row.update(iptm=score,confidence_sha256=digest(files[0]),structure_sha256=digest(structures[0]))
            except (ValueError,KeyError,OSError,TypeError) as exc:
                row['issues'].append(str(exc))
            rows.append(row)
    return {'scope':'Historical identity/integrity audit only. Missing model, software, seed and checkpoint fingerprints prevent verified cache reuse.',
            'n_runs':len(rows),'n_runs_with_issues':sum(bool(row['issues']) for row in rows),'runs':rows}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-root',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();report=audit(args.data_root)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(report,indent=2)+'\n')
    print(f"Audited {report['n_runs']} historical runs; {report['n_runs_with_issues']} have identity/integrity issues")


if __name__=='__main__':main()
