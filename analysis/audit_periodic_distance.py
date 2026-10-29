#!/usr/bin/env python3
"""Reanalyse existing endpoint distances with minimum-image periodic boundaries.

Reports descriptive distances only; does not infer enzyme activity or leakiness.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path


def measure_frames(universe, a, b):
    import numpy as np
    from MDAnalysis.lib.distances import distance_array
    rows = []
    for frame in universe.trajectory:
        box = frame.dimensions
        if box is None or not np.isfinite(box).all() or np.any(box[:3] <= 0) or np.any(box[3:] <= 0) or np.any(box[3:] >= 180):
            raise ValueError(f'Missing or invalid box at frame {frame.frame}')
        raw = float(np.linalg.norm(a.positions[0] - b.positions[0]))
        periodic = float(distance_array(a.positions, b.positions, box=box)[0, 0])
        if not np.isfinite([frame.time, raw, periodic]).all():
            raise ValueError('Nonfinite time or distance')
        if rows and frame.time / 1000 <= rows[-1][0]:
            raise ValueError('Trajectory times must strictly increase')
        rows.append((frame.time / 1000, raw, periodic))
    if not rows:
        raise ValueError('Empty trajectory')
    return rows


def analyse(directory, first, second, output):
    import MDAnalysis as mda
    import numpy as np

    topology, trajectory = directory / 'prod.tpr', directory / 'prod.xtc'
    universe = mda.Universe(str(topology), str(trajectory))
    a, b = universe.select_atoms(first), universe.select_atoms(second)
    if len(a) != 1 or len(b) != 1:
        raise ValueError(f'Expected one atom per selection, got {len(a)}, {len(b)}')
    rows = measure_frames(universe, a, b)
    data = np.array(rows)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open('w') as handle:
        writer = csv.writer(handle, delimiter='\t')
        writer.writerow(['time_ns', 'raw_distance_A', 'minimum_image_distance_A'])
        writer.writerows(rows)
    report = {
        'source_directory': str(directory.resolve()),
        'selections': [first, second],
        'selected_atoms': [str(a[0]), str(b[0])],
        'topology_sha256': hashlib.sha256(topology.read_bytes()).hexdigest(),
        'trajectory_bytes': trajectory.stat().st_size,
        'trajectory_mtime_ns': trajectory.stat().st_mtime_ns,
        'trajectory_identity_note': 'Size and mtime only; not a content hash.',
        'timeseries_sha256': hashlib.sha256(output.read_bytes()).hexdigest(),
        'MDAnalysis_version': mda.__version__, 'numpy_version': np.__version__,
        'n_frames': len(rows), 'elapsed_ns': float(data[-1, 0] - data[0, 0]),
        'raw_pct_below_15A': float(np.mean(data[:, 1] < 15) * 100),
        'minimum_image_pct_below_15A': float(np.mean(data[:, 2] < 15) * 100),
        'raw_mean_A': float(data[:, 1].mean()),
        'minimum_image_mean_A': float(data[:, 2].mean()),
        'interpretation': 'Endpoint proximity only. Does not establish catalytic geometry, activity, or safety.',
    }
    output.with_suffix('.json').write_text(json.dumps(report, indent=2, allow_nan=False) + '\n')
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory', type=Path)
    parser.add_argument('--first', required=True)
    parser.add_argument('--second', required=True)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    print(json.dumps(analyse(args.directory, args.first, args.second, args.output), indent=2))


if __name__ == '__main__':
    main()
