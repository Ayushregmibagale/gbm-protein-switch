#!/usr/bin/env python3
"""Audit recorded evidence without running design jobs or interpreting scores as activity."""
import argparse
import csv
import hashlib
import json
import math
import statistics
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SEEDS = {0, 1, 2}
CONDITIONS = {'SOX2+SOX9', 'SOX11+SOX9', 'HMGB1+SOX9', 'SRY+SOX9',
              'SOX2+SOX11', 'SOX2+HMGB1', 'SOX2+SRY'}


def number(value, low=0, high=1):
    result = float(value)
    if not math.isfinite(result) or not low <= result <= high:
        raise ValueError(f'Invalid numeric value: {value}')
    return result


def read_series(path, columns):
    rows = []
    for line in path.read_text().splitlines():
        if not line.strip() or line.lstrip().startswith('#'):
            continue
        row = [float(x) for x in line.split()]
        if len(row) != columns or not all(math.isfinite(x) for x in row):
            raise ValueError(f'{path}: malformed or nonfinite row')
        if row[0] < 0 or (rows and row[0] <= rows[-1][0]):
            raise ValueError(f'{path}: times must be nonnegative and strictly increasing')
        rows.append(row)
    if not rows:
        raise ValueError(f'{path}: empty time series')
    return rows


def control_summary(rows):
    groups = {}
    for row in rows:
        condition, seed = row['condition'], int(row['seed'])
        if seed not in SEEDS or condition not in CONDITIONS:
            raise ValueError('Unexpected condition or seed')
        group = groups.setdefault(condition, {})
        if seed in group:
            raise ValueError('Duplicate condition/seed')
        group[seed] = number(row['iptm'])
    if set(groups) != CONDITIONS or any(set(g) != SEEDS for g in groups.values()):
        raise ValueError('Missing condition or seed; comparison is incomplete')
    baseline = groups['SOX2+SOX9']
    result = {}
    for condition, values in sorted(groups.items()):
        scores = list(values.values())
        deltas = [baseline[s] - values[s] for s in sorted(SEEDS)]
        result[condition] = {
            'n_seeds': len(scores), 'mean_iptm': statistics.mean(scores),
            'sd_iptm': statistics.stdev(scores), 'min_iptm': min(scores),
            'max_iptm': max(scores), 'paired_baseline_minus_condition': deltas,
            'mean_paired_difference': statistics.mean(deltas),
            'all_paired_differences_positive': all(d > 0 for d in deltas),
        }
    return result


def audit(root):
    inputs = []
    def source(relative):
        path = root / relative
        inputs.append({'path': relative, 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()})
        return path
    with source('analysis/bridge_3chain_selectivity_results.tsv').open() as handle:
        controls = control_summary(csv.DictReader(handle, delimiter='\t'))
    with source('analysis/dual_contact_results.tsv').open() as handle:
        recorded = list(csv.DictReader(handle, delimiter='\t'))
    contacts, seen = {}, set()
    for row in recorded:
        design = row['design']
        if design in seen or Path(design).name != design:
            raise ValueError('Duplicate or invalid design name')
        seen.add(design)
        series = read_series(source(f'analysis/{design}_dual_contact_timeseries.txt'), 3)
        if len(series) != int(row['n_frames']):
            raise ValueError(f'{design}: frame count mismatch')
        for _, a, b in series:
            number(a); number(b)
        checks = {}
        for label, col in [('sox2', 1), ('sox9', 2)]:
            mean = statistics.mean(r[col] for r in series)
            if abs(mean - number(row[f'{label}_mean'])) > 0.0006:
                raise ValueError(f'{design}: {label} mean mismatch')
            checks[f'{label}_mean'] = mean
        for threshold in (0.05, 0.10, 0.15):
            key = f'thr{round(threshold*100):02d}'
            states = {'dual': 0, 'sox2only': 0, 'sox9only': 0, 'neither': 0}
            for _, a, b in series:
                state = ('dual' if a > threshold and b > threshold else
                         'sox2only' if a > threshold else
                         'sox9only' if b > threshold else 'neither')
                states[state] += 1
            checks[key] = {k: v / len(series) * 100 for k, v in states.items()}
            for state, pct in checks[key].items():
                if abs(pct - number(row[f'{key}_{state}_pct'], 0, 100)) > 0.051:
                    raise ValueError(f'{design}: {key}/{state} percentage mismatch')
        contacts[design] = {'n_frames': len(series),
                            'elapsed_ns': series[-1][0] - series[0][0], **checks}
    if not contacts:
        raise ValueError('No recorded designs')
    return {'scope': 'Recorded tables and rounded time series; no activity or safety verdict',
            'three_chain_controls': controls, 'dual_contact': contacts,
            'limitations': ['Frames are correlated, not independent replicates.',
                            'Contact occupancy is not evidence of conditional caspase activity.',
                            'ipTM is structural confidence, not affinity or measured specificity.',
                            'Raw trajectories and Track 1 sources are absent from this checkout.'],
            'inputs': inputs}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=ROOT)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    try:
        report = audit(args.root)
    except (OSError, ValueError, KeyError) as exc:
        parser.exit(1, f'Evidence audit failed: {exc}\n')
    text = json.dumps(report, indent=2, allow_nan=False) + '\n'
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text)
    else:
        print(text, end='')


if __name__ == '__main__':
    main()
