import csv
import importlib.util
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('ranking', ROOT/'selection/stage4_composite_ranking.py')
ranking = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ranking)


class RankingTests(unittest.TestCase):
    def rows(self):
        with (ROOT/'selection/stage3_results.tsv').open() as handle:
            return list(csv.DictReader(handle, delimiter='\t'))

    def test_historical_ranking_reproduces(self):
        result = ranking.score_rows(self.rows())
        with (ROOT/'selection/stage4_results.tsv').open() as handle:
            old = list(csv.DictReader(handle, delimiter='\t'))
        self.assertEqual([row['pair'] for row in result], [row['pair'] for row in old])
        for actual, expected in zip(result, old):
            self.assertAlmostEqual(actual['composite_score'], float(expected['composite_score']))

    def test_bad_scores_fail_instead_of_becoming_top_ranked(self):
        for value in ['0', 'nan', '-1', 'inf', '']:
            row = next(row for row in self.rows() if row['both_selective'] == 'True')
            row['sox2_wrong_iptm'] = value
            with self.subTest(value=value), self.assertRaises(ValueError):
                ranking.score_rows([row])
