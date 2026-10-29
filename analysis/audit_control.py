#!/usr/bin/env python3
"""Validate trajectory coverage and compute periodic endpoint-distance observations.

All fragments must be immutable and explicitly ordered. A deliberately trimmed
first fragment is recorded; gaps, overlap, corruption and partial coverage fail.
No activity, specificity, or safety verdict is inferred.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import statistics


def check_times(times, start=0.0, end=20.0, cadence=0.01):
    if not times or cadence <= 0 or end <= start:
        raise ValueError('Empty trajectory or invalid expected interval/cadence')
    tolerance = max(1e-6, cadence * 0.001)
    if not all(math.isfinite(t) for t in times):
        raise ValueError('Nonfinite trajectory time')
    if abs(times[0] - start) > tolerance or abs(times[-1] - end) > tolerance:
        raise ValueError(f'Incomplete coverage: {times[0]:.6f}–{times[-1]:.6f} ns; expected {start}–{end} ns')
    for a, b in zip(times, times[1:]):
        if abs((b - a) - cadence) > tolerance:
            raise ValueError(f'Gap, overlap, or unexpected cadence at {a:.6f} -> {b:.6f} ns')
    expected = round((end-start)/cadence) + 1
    if len(times) != expected:
        raise ValueError(f'Expected {expected} frames, found {len(times)}')


def digest(path):
    result = hashlib.sha256()
    with path.open('rb') as handle:
        for chunk in iter(lambda: handle.read(4 * 1024 * 1024), b''):
            result.update(chunk)
    return result.hexdigest()


def identity(path):
    stat = path.stat()
    return stat.st_size, stat.st_mtime_ns


def analyse_control(topology, fragments, first, second, output, label='control',
                    criterion='far', start=0.0, end=20.0, cadence=0.01,
                    first_fragment_end=None):
    import MDAnalysis as mda
    import numpy as np
    from MDAnalysis.lib.formats.libmdaxdr import XTCFile
    from MDAnalysis.lib.distances import distance_array
    from MDAnalysis.lib.mdamath import triclinic_box

    if not fragments or len({p.resolve() for p in fragments}) != len(fragments):
        raise ValueError('Provide unique, ordered fragments')
    if criterion not in {'near', 'far'}:
        raise ValueError('Criterion must be near or far')
    if first_fragment_end is not None and (len(fragments) < 2 or not start < first_fragment_end < end):
        raise ValueError('A first-fragment trim requires multiple fragments and an interior boundary')
    if output.exists():
        raise FileExistsError(f'Output already exists; use a new audit directory: {output}')
    universe = mda.Universe(str(topology))
    selections = [universe.select_atoms(s) for s in (first, second)]
    if any(len(selection) != 1 for selection in selections):
        raise ValueError('Each endpoint selection must identify exactly one atom')
    indices = [int(selection.indices[0]) for selection in selections]
    rows, sources = [], []
    for index, fragment in enumerate(fragments):
        before = identity(fragment)
        accepted = 0
        with XTCFile(str(fragment), 'r') as reader:
            # No offset/index files are written alongside the original trajectory.
            for frame in reader:
                time = float(frame.time / 1000)
                if index == 0 and first_fragment_end is not None and time > first_fragment_end + 1e-6:
                    break
                if len(frame.x) != universe.atoms.n_atoms:
                    raise ValueError(f'{fragment}: atom count differs from topology')
                box = triclinic_box(*(np.asarray(frame.box) * 10))
                if not np.isfinite(box).all() or np.any(box[:3] <= 0) or np.any(box[3:] <= 0) or np.any(box[3:] >= 180):
                    raise ValueError(f'{fragment}: invalid periodic box')
                positions = np.asarray(frame.x)[indices] * 10
                if not np.isfinite(positions).all():
                    raise ValueError(f'{fragment}: nonfinite endpoint coordinates')
                raw = float(np.linalg.norm(positions[0] - positions[1]))
                periodic = float(distance_array(positions[:1], positions[1:], box=box)[0, 0])
                rows.append([time, raw, periodic])
                accepted += 1
        checksum = digest(fragment)
        if identity(fragment) != before:
            raise RuntimeError(f'{fragment}: changed during audit; wait for the writer to finish')
        sources.append({'path': str(fragment.resolve()), 'bytes': before[0],
                        'mtime_ns': before[1], 'sha256': checksum, 'accepted_frames': accepted})
    check_times([row[0] for row in rows], start, end, cadence)
    distances = [row[2] for row in rows]
    near = sum(d < 15 for d in distances) / len(rows) * 100
    far = sum(d > 25 for d in distances) / len(rows) * 100
    report = {
        'label': label, 'scope': 'Descriptive endpoint distances; biological activity and safety are untested',
        'n_frames': len(rows), 'start_ns': rows[0][0], 'end_ns': rows[-1][0],
        'elapsed_ns': rows[-1][0] - rows[0][0], 'cadence_ns': cadence,
        'mean_dist_A': statistics.mean(distances), 'std_dist_A': statistics.pstdev(distances),
        'min_dist_A': min(distances), 'max_dist_A': max(distances),
        'pct_below_15A': near, 'pct_above_25A': far,
        'distance_criterion': '<15 A in >50% frames' if criterion == 'near' else '>25 A in >=95% frames',
        'distance_criterion_met': near > 50 if criterion == 'near' else far >= 95,
        'first_fragment_end_ns': first_fragment_end,
        'topology_path': str(topology.resolve()), 'topology_sha256': digest(topology),
        'selections': [first, second], 'atoms': [str(s[0]) for s in selections],
        'sources': sources, 'software': {'MDAnalysis': mda.__version__, 'numpy': np.__version__},
    }
    output.mkdir(parents=True)
    series = output / 'distances.tsv'
    np.savetxt(series, rows, delimiter='\t', header='time_ns\traw_distance_A\tminimum_image_distance_A', fmt='%.8f')
    report['timeseries_sha256'] = digest(series)
    (output / 'metrics.json').write_text(json.dumps(report, indent=2, allow_nan=False) + '\n')
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--topology', required=True, type=Path)
    parser.add_argument('--fragments', required=True, nargs='+', type=Path)
    parser.add_argument('--first', required=True)
    parser.add_argument('--second', required=True)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--criterion', choices=['near', 'far'], default='far')
    parser.add_argument('--label', default='control')
    parser.add_argument('--expected-start-ns', type=float, default=0)
    parser.add_argument('--expected-end-ns', type=float, default=20)
    parser.add_argument('--cadence-ns', type=float, default=0.01)
    parser.add_argument('--first-fragment-end-ns', type=float)
    args = parser.parse_args()
    try:
        result = analyse_control(args.topology, args.fragments, args.first, args.second,
                                 args.output, label=args.label, criterion=args.criterion,
                                 start=args.expected_start_ns, end=args.expected_end_ns,
                                 cadence=args.cadence_ns, first_fragment_end=args.first_fragment_end_ns)
    except (OSError, ValueError, RuntimeError) as exc:
        parser.exit(1, f'Control audit failed: {exc}\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
