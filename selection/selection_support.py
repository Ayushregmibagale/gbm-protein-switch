"""Configuration and explicit execution gate for recovered prediction campaigns."""
import argparse
import os
import json
import sys
from pathlib import Path


def require_execution(module_name):
    if module_name != '__main__':
        return
    parser = argparse.ArgumentParser(
        description='Recovered prediction campaign. Set GBM_ROOT and BOLTZ_BIN for your environment. '
                    'Historical confidence thresholds are heuristics, not measured selectivity.',
        add_help=True)
    parser.add_argument('--execute', action='store_true', help='Explicitly run the legacy campaign')
    args, remaining = parser.parse_known_args()
    dry_run = Path(sys.argv[0]).name == 'stage1_reprediction.py' and '--dry_run' in remaining
    if not args.execute and not dry_run:
        parser.error('Use --execute to launch inference, or stage1 --dry_run to inspect inputs')
    if args.execute and Path(sys.argv[0]).name != 'stage1_reprediction.py':
        root = os.environ.get('GBM_RUN_ROOT')
        if not root:
            parser.error('Set GBM_RUN_ROOT to a new campaign directory; historical caches cannot be reused implicitly')
        output = Path(root).resolve()
        if output.exists():
            parser.error(f'Campaign directory already exists: {output}. Choose a new directory.')
        output.mkdir(parents=True)
        (output/'campaign.json').write_text(json.dumps({'script':str(Path(sys.argv[0]).resolve()),
            'arguments':sys.argv[1:],'data_root':os.environ.get('GBM_ROOT'),
            'stage1_table':os.environ.get('GBM_STAGE1_TABLE'),
            'cache_policy':'fresh outputs; no implicit historical cache reuse'},indent=2)+'\n')
        os.environ['GBM_ACTIVE_RUN_DIR']=str(output)
    sys.argv = [sys.argv[0], *remaining]


def reference_sequence():
    path = Path(__file__).resolve().parents[1] / 'inputs/reference/P55211.fasta'
    sequence = ''.join(path.read_text().splitlines()[1:]).strip()
    if len(sequence) != 416 or sequence[314:316] != 'DA':
        raise ValueError('Reference sequence identity differs from the audited 416-aa input')
    return sequence


def validate_recorded_candidates(frame, root):
    """Refuse stale candidate identities before starting another prediction stage."""
    from Bio.PDB import PDBParser
    from Bio.Data.PDBData import protein_letters_3to1
    parser = PDBParser(QUIET=True)
    invalid = []
    for _, row in frame.iterrows():
        if row['source'] != 'bindcraft':
            continue
        directory = {'SOX2': 'sox2_binder_xtal', 'SOX9': 'sox9_binder_xtal_v3'}[row['target']]
        path = root/'bindcraft/designs'/directory/'Accepted'/(row['name']+'.pdb')
        model = parser.get_structure(row['name'], path)[0]
        if 'B' not in model:
            raise ValueError(f'Missing explicit binder chain B: {path}')
        expected = ''.join(protein_letters_3to1[r.resname] for r in model['B'] if r.id[0] == ' ')
        if not expected or expected != row['binder_seq']:
            invalid.append(row['name'])
    if invalid:
        raise ValueError(f'{len(invalid)} candidate identities disagree with binder chain B. '
                         'Regenerate Stage 1 from explicit binder chains before downstream inference: ' + ', '.join(invalid))


def campaign_output(default):
    root = os.environ.get('GBM_ACTIVE_RUN_DIR')
    return Path(root)/'predictions' if root else default


def campaign_table(default):
    root = os.environ.get('GBM_ACTIVE_RUN_DIR')
    return Path(root)/Path(default).name if root else default
