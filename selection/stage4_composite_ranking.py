#!/usr/bin/env python3
"""Reproduce the historical ranking heuristic from explicitly supplied recorded data.

Confidence ratios are descriptive model-score ratios, not measured specificity.
"""
import argparse
import csv
import math
from pathlib import Path


def score_rows(rows):
    result, seen = [], set()
    for row in rows:
        pair = row['pair']
        if pair in seen:
            raise ValueError(f'Duplicate pair: {pair}')
        seen.add(pair)
        if row['both_selective'] not in ('True', 'False'):
            raise ValueError('both_selective must be an explicit True/False historical flag')
        if row['both_selective'] != 'True':
            continue
        values = {}
        for key in ['pair_iptm', 'sox2_correct_iptm', 'sox2_wrong_iptm', 'sox9_correct_iptm', 'sox9_wrong_iptm', 'terminus_dist_A']:
            value = float(row[key])
            if not math.isfinite(value) or value < 0 or (key != 'terminus_dist_A' and value > 1):
                raise ValueError(f'{pair}: invalid {key}')
            values[key] = value
        for key in ['sox2_wrong_iptm', 'sox9_wrong_iptm']:
            if values[key] == 0:
                raise ValueError(f'{pair}: zero confidence denominator; ratio ranking undefined')
        s2 = values['sox2_correct_iptm'] / values['sox2_wrong_iptm']
        s9 = values['sox9_correct_iptm'] / values['sox9_wrong_iptm']
        geometry = math.exp(-values['terminus_dist_A'] / 40)
        score = values['pair_iptm'] * s2 * s9 * geometry
        if not math.isfinite(score):
            raise ValueError(f'{pair}: nonfinite composite score')
        result.append({**row, 'sox2_sel_ratio': s2, 'sox9_sel_ratio': s9,
                       'geometry_bonus': geometry, 'composite_score': score})
    if not result:
        raise ValueError('No rows passed the recorded historical filter')
    return sorted(result, key=lambda row: (-row['composite_score'], row['pair']))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    try:
        with args.input.open() as handle:
            rows = score_rows(list(csv.DictReader(handle, delimiter='\t')))
        args.output.parent.mkdir(parents=True, exist_ok=True)
        with args.output.open('x', newline='') as handle:
            writer = csv.DictWriter(handle, fieldnames=['rank', *rows[0]], delimiter='\t')
            writer.writeheader()
            writer.writerows({'rank': rank, **row} for rank, row in enumerate(rows, 1))
    except (OSError, ValueError, KeyError) as exc:
        parser.exit(1, f'Ranking failed: {exc}\n')
    print(f'Recorded heuristic ranking: {len(rows)} pairs -> {args.output}')


if __name__ == '__main__':
    main()
