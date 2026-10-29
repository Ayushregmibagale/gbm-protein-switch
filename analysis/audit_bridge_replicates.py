#!/usr/bin/env python3
"""Audit every available bridge replicate with periodic contacts and window sensitivity."""
import argparse
import json
from pathlib import Path
import statistics
import numpy as np
from audit_control import check_times, digest, identity

DESIGNS = ['bridge_48_s01', 'bridge_73_s09', 'bridge_136_s01', 'bridge_101_s06', 'bridge_90_s01']


def analyse(directory, output, expected_sequences=None):
    import MDAnalysis as mda
    from MDAnalysis.lib.formats.libmdaxdr import XTCFile
    from MDAnalysis.lib.distances import distance_array
    from MDAnalysis.lib.mdamath import triclinic_box
    from Bio.Data.PDBData import protein_letters_3to1
    top, trajectory = directory/'prod.tpr', directory/'prod.xtc'
    universe = mda.Universe(str(top))
    groups = [segment.atoms.select_atoms('protein and name CA') for segment in universe.segments]
    groups = [group for group in groups if len(group)]
    if len(groups) != 3 or [len(g) for g in groups[:2]] != [97,72]:
        raise ValueError(f'{directory}: expected three protein segments, target lengths 97/72')
    sequences = [''.join(protein_letters_3to1[r.resname] for r in g.residues) for g in groups]
    if expected_sequences is not None and sequences != expected_sequences:
        raise ValueError(f'{directory}: replicate protein identities differ')
    before = identity(trajectory)
    rows = []
    with XTCFile(str(trajectory),'r') as reader:
        for frame in reader:
            if len(frame.x) != universe.atoms.n_atoms:
                raise ValueError('Trajectory atom count differs from topology')
            box = triclinic_box(*(np.asarray(frame.box)*10))
            if not np.isfinite(box).all() or np.any(box[:3]<=0):
                raise ValueError('Invalid periodic box')
            positions = [np.asarray(frame.x)[g.indices]*10 for g in groups]
            if not all(np.isfinite(p).all() for p in positions):
                raise ValueError('Nonfinite coordinates')
            minima = [distance_array(positions[2],positions[index],box=box).min(axis=1) for index in (0,1)]
            rows.append([frame.time/1000,*[float(np.mean(d<cutoff)) for cutoff in (6,8,10) for d in minima]])
    check_times([row[0] for row in rows])
    checksum = digest(trajectory)
    if identity(trajectory) != before:
        raise RuntimeError('Trajectory changed during audit')
    values = np.array(rows)
    sensitivity = []
    for discard in (0,5,10):
        window = values[values[:,0]>=discard]
        for index,cutoff in enumerate((6,8,10)):
            a,b=window[:,1+2*index],window[:,2+2*index]
            for threshold in (.05,.10,.15):
                sensitivity.append({'discard_ns':discard,'cutoff_A':cutoff,'fraction_threshold':threshold,
                                    'n_frames':len(window),'dual_pct':float(np.mean((a>threshold)&(b>threshold))*100)})
    output.mkdir(parents=True)
    series=output/'contacts.tsv'
    np.savetxt(series,values,delimiter='\t',header='time_ns s2_6A s9_6A s2_8A s9_8A s2_10A s9_10A',fmt='%.8f')
    report={'directory':str(directory),'n_frames':len(rows),'elapsed_ns':rows[-1][0]-rows[0][0],
            'sequences':sequences,'sensitivity':sensitivity,'topology_sha256':digest(top),
            'trajectory_sha256':checksum,'trajectory_bytes':before[0],
            'timeseries_sha256':digest(series),'software':{'MDAnalysis':mda.__version__,'numpy':np.__version__}}
    (output/'metrics.json').write_text(json.dumps(report,indent=2)+'\n')
    return report


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--md-root',required=True,type=Path)
    parser.add_argument('--output',required=True,type=Path)
    args=parser.parse_args()
    if args.output.exists():parser.error('Choose a new output directory')
    cohort=[]
    for design in DESIGNS:
        names=[design,*[f"{design if design!='bridge_90_s01' else 'bridge_90'}_rep{i}" for i in (2,3)]]
        reports=[]
        for name in names:
            print(f'Analysing {name}',flush=True)
            reports.append(analyse(args.md_root/name,args.output/name,reports[0]['sequences'] if reports else None))
        for index,setting in enumerate(reports[0]['sensitivity']):
            occupancies=[report['sensitivity'][index]['dual_pct'] for report in reports]
            cohort.append({'design':design,**{k:v for k,v in setting.items() if k!='dual_pct'},
                           'replicate_values_pct':occupancies,'n_replicates':len(reports),
                           'replicate_mean_pct':statistics.mean(occupancies),
                           'replicate_sd_pct':statistics.stdev(occupancies),
                           'replicate_min_pct':min(occupancies),'replicate_max_pct':max(occupancies)})
    (args.output/'cohort.json').write_text(json.dumps({'scope':'Replicate-level summaries; frames are correlated; shared equilibrated starts limit independence',
                                                     'sensitivity':cohort},indent=2)+'\n')
    print('Completed 15 trajectories',flush=True)


if __name__=='__main__':main()
