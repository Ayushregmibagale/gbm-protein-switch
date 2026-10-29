#!/usr/bin/env python3
"""Describe simulated caspase-core similarity to experimental 1JXQ conformations.

Sequence-aligned CA comparisons and native inter-half contacts are descriptive,
not a validated enzyme-activity classifier. Reference missing residues are excluded.
"""
import argparse
import json
from pathlib import Path
import numpy as np
from audit_control import check_times,digest,identity


def rmsd_fit(moving,fixed):
    a,b=moving-moving.mean(axis=0),fixed-fixed.mean(axis=0)
    u,_,vt=np.linalg.svd(a.T@b)
    if np.linalg.det(u@vt)<0:u[:,-1]*=-1
    return float(np.sqrt(np.mean(np.sum((a@u@vt-b)**2,axis=1))))


def references(path,sequence):
    from Bio.PDB import PDBParser
    from Bio.Data.PDBData import protein_letters_3to1
    from Bio.Align import PairwiseAligner
    model=PDBParser(QUIET=True).get_structure('1JXQ',path)[0]
    result={}
    for chain in ['A','B']:
        residues=[r for r in model[chain] if r.id[0]==' ' and 'CA' in r]
        native=''.join(protein_letters_3to1[r.resname] for r in residues)
        aligner=PairwiseAligner(mode='global',match_score=2,mismatch_score=-4,open_gap_score=-6,extend_gap_score=-.5)
        alignment=aligner.align(sequence,native)[0]
        mapped=[(int(a),int(b)) for a,b in zip(*alignment.indices) if a>=139 and a<416 and b>=0 and sequence[a]==native[b]]
        if len(mapped)<200:raise ValueError('Insufficient reference sequence alignment')
        result[chain]=[(a+1,residues[b]['CA'].coord.tolist()) for a,b in mapped]
    return result


def analyse(directory,reference,sequence,output):
    import MDAnalysis as mda
    from MDAnalysis.lib.formats.libmdaxdr import XTCFile
    from MDAnalysis.lib.distances import distance_array,minimize_vectors
    from MDAnalysis.lib.mdamath import triclinic_box
    u=mda.Universe(str(directory/'prod.tpr'))
    groups=[seg.atoms.select_atoms('protein and name CA') for seg in u.segments]
    groups=[g for g in groups if len(g)]
    if len(groups)!=4 or [len(g) for g in groups]!=[97,72,250,171]:raise ValueError('Unexpected bound construct topology')
    core=np.concatenate([groups[2].indices[-176:],groups[3].indices[:101]])
    refs=references(reference,sequence)
    prepared={}
    for chain,pairs in refs.items():
        nums=np.array([p[0] for p in pairs]);coords=np.array([p[1] for p in pairs]);n=nums<=315;c=nums>=316
        contacts=np.linalg.norm(coords[n,None,:]-coords[None,c,:],axis=2)<8
        prepared[chain]=(nums,coords,n,c,contacts)
    rows=[];times=[];before=identity(directory/'prod.xtc')
    with XTCFile(str(directory/'prod.xtc'),'r') as reader:
        for frame_index,frame in enumerate(reader):
            times.append(frame.time/1000)
            if frame_index%10:continue
            box=triclinic_box(*(np.asarray(frame.box)*10))
            if not np.isfinite(box).all() or np.any(box[:3]<=0):raise ValueError('Invalid box')
            raw=np.asarray(frame.x)[core]*10
            halves=[]
            for part in (raw[:176],raw[176:]):
                differences=minimize_vectors(part[1:]-part[:-1],box)
                halves.append(np.concatenate([part[:1],part[:1]+np.cumsum(differences,axis=0)]))
            shift=halves[0][-1]+minimize_vectors((halves[1][0]-halves[0][-1])[None,:],box)[0]-halves[1][0]
            halves[1]+=shift
            positions=np.concatenate(halves)
            for chain,(nums,coords,n,c,contacts) in prepared.items():
                current=positions[nums-140]
                contact_dist=distance_array(current[n],current[c],box=box)
                rows.append({'time_ns':frame.time/1000,'reference_chain':chain,
                             'assembly_CA_rmsd_A':rmsd_fit(current,coords),
                             'N_half_CA_rmsd_A':rmsd_fit(current[n],coords[n]),
                             'C_half_CA_rmsd_A':rmsd_fit(current[c],coords[c]),
                             'native_interhalf_contact_fraction':float(np.mean(contact_dist[contacts]<8))})
    check_times(times)
    checksum=digest(directory/'prod.xtc')
    if identity(directory/'prod.xtc')!=before:raise RuntimeError('Trajectory changed')
    summary={}
    for chain in refs:
        selected=[row for row in rows if row['reference_chain']==chain]
        summary[chain]={key:{'mean':float(np.mean([row[key] for row in selected])),
                              'min':float(np.min([row[key] for row in selected])),
                              'max':float(np.max([row[key] for row in selected]))}
                        for key in selected[0] if key not in ['reference_chain','time_ns']}
    output.mkdir(parents=True)
    report={'scope':__doc__,'directory':str(directory),'reference_url':'https://www.rcsb.org/structure/1JXQ',
            'reference_sha256':digest(reference),'trajectory_sha256':checksum,'topology_sha256':digest(directory/'prod.tpr'),
            'sample_stride_frames':10,'reference_mapping':refs,'summary':summary,'observations':rows,
            'reconstruction':'Each half reconstructed by consecutive CA minimum-image displacements; C-half placed at closest image relative to N-half split endpoint. Native contacts use PBC independently.'}
    (output/'geometry.json').write_text(json.dumps(report,indent=2)+'\n')
    print(directory.name,summary,flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--md-root',type=Path,required=True);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();root=Path(__file__).resolve().parents[1]
    if args.output.exists():parser.error('Use a new output directory')
    sequence=''.join((root/'inputs/reference/P55211.fasta').read_text().splitlines()[1:])
    for name in ['lead_construct','rep2','rep3']:
        analyse(args.md_root/name,root/'inputs/reference/1JXQ.pdb',sequence,args.output/name)


if __name__=='__main__':main()
