"""Shared checks for the recovered workstation workflows."""
import json
import os
from pathlib import Path
import shlex
import subprocess
import hashlib


def recorded_control_audit(directory):
    """Recognise the validated control without trusting a stale legacy state flag."""
    path=Path(__file__).resolve().parents[1]/'evidence/unbound_control/metrics.json'
    if not path.exists():return None
    report=json.loads(path.read_text())
    topology=directory/'prod.tpr'
    if not topology.exists() or topology.resolve()!=Path(report['topology_path']).resolve():return None
    if hashlib.sha256(topology.read_bytes()).hexdigest()!=report['topology_sha256']:return None
    for source in report['sources']:
        candidate=Path(source['path'])
        if not candidate.exists():return None
        stat=candidate.stat()
        if stat.st_size!=source['bytes'] or stat.st_mtime_ns!=source['mtime_ns']:return None
    series=path.parent/'distances.tsv'
    if hashlib.sha256(series.read_bytes()).hexdigest()!=report['timeseries_sha256']:return None
    return report

STAGE_OUTPUTS = {
    'STRUCTURE_DONE': ['start.pdb'],
    'TOPOLOGY_DONE': ['topol.top', 'ions.gro'],
    'EM_DONE': ['em.gro'],
    'NVT_DONE': ['nvt.gro', 'nvt.cpt'],
    'NPT_DONE': ['npt.gro', 'npt.cpt'],
    'SETUP_DONE': ['topol.top', 'npt.gro', 'prod_newvel.mdp'],
    'PROD_DONE': ['prod.tpr', 'prod.xtc'],
    'ANALYSIS_DONE': ['metrics.json'],
}


def require_files(directory, names):
    missing = [str(directory / name) for name in names
               if not (directory / name).is_file() or (directory / name).stat().st_size == 0]
    if missing:
        raise RuntimeError('Missing or empty required files: ' + ', '.join(missing))


def read_state(directory):
    path = directory / 'state.txt'
    state = path.read_text().strip() if path.exists() else 'PENDING'
    if state not in {'PENDING', *STAGE_OUTPUTS}:
        raise ValueError(f'{path}: unknown state {state!r}')
    if state != 'PENDING':
        require_files(directory, STAGE_OUTPUTS[state])
    return state


def write_state(directory, state):
    if state not in STAGE_OUTPUTS:
        raise ValueError(f'Unknown output state: {state}')
    require_files(directory, STAGE_OUTPUTS[state])
    directory.mkdir(parents=True, exist_ok=True)
    temporary = directory / 'state.txt.tmp'
    temporary.write_text(state + '\n')
    temporary.replace(directory / 'state.txt')


def run_gmx(executable, arguments, log, directory, stdin=None):
    """No shell, checked return code, append diagnostics instead of erasing history."""
    argv = [str(executable), *shlex.split(arguments)] if isinstance(arguments, str) else [str(executable), *arguments]
    with Path(log).open('a') as handle:
        handle.write('\nCOMMAND ' + shlex.join(argv) + '\n')
        handle.flush()
        subprocess.run(argv, cwd=directory, input=stdin, text=True,
                       stdout=handle, stderr=subprocess.STDOUT, check=True)


def production(executable, directory, prepare, log_prefix, resources=None):
    """Resume a checkpoint using its existing TPR; refuse ambiguous partial output."""
    resources = resources or ['-ntmpi', '1', '-ntomp', '8', '-gpu_id', '0', '-nb', 'gpu', '-pme', 'gpu', '-bonded', 'gpu']
    checkpoint = directory / 'prod.cpt'
    tpr = directory / 'prod.tpr'
    existing = [p for p in directory.glob('prod*') if p.suffix in {'.xtc', '.trr', '.edr', '.log'}]
    if checkpoint.exists():
        require_files(directory, ['prod.cpt', 'prod.tpr'])
        argv = ['mdrun', '-v', '-s', 'prod.tpr', '-deffnm', 'prod', '-cpi', 'prod.cpt', '-noappend', *resources]
    elif tpr.exists() or existing:
        raise RuntimeError('Partial production output without a checkpoint. Inspect it; no files were overwritten.')
    else:
        run_gmx(executable, prepare, directory / f'{log_prefix}_grompp.log', directory)
        require_files(directory, ['prod.tpr'])
        argv = ['mdrun', '-v', '-deffnm', 'prod', *resources]
    run_gmx(executable, argv, directory / f'{log_prefix}_mdrun.log', directory)
    # -noappend can produce fragments: do not pretend these are a complete prod.xtc.
    if list(directory.glob('prod.part*.xtc')):
        raise RuntimeError('Production wrote trajectory fragments. Validate and assemble them before marking PROD_DONE.')
    require_files(directory, ['prod.tpr', 'prod.xtc'])


def periodic_analysis(directory, first, second, label, criterion):
    """Validate full 0–20 ns coverage before issuing descriptive distance metrics."""
    import sys
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'analysis'))
    from audit_control import analyse_control
    output = directory / 'validated_analysis'
    return analyse_control(directory / 'prod.tpr', [directory / 'prod.xtc'],
                           first, second, output, label=label, criterion=criterion)
