import copy
import importlib.util
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('audit', ROOT / 'analysis/audit_evidence.py')
audit = importlib.util.module_from_spec(spec)
spec.loader.exec_module(audit)


class EvidenceTests(unittest.TestCase):
    def rows(self):
        return [dict(condition=c, seed=s, iptm=0.5 if c == 'SOX2+SOX9' else 0.3)
                for c in sorted(audit.CONDITIONS) for s in sorted(audit.SEEDS)]

    def test_measured_baseline_and_paired_seeds(self):
        result = audit.control_summary(self.rows())
        self.assertAlmostEqual(result['SOX11+SOX9']['mean_paired_difference'], 0.2)
        self.assertTrue(result['SOX11+SOX9']['all_paired_differences_positive'])

    def test_missing_duplicate_and_unexpected_seeds(self):
        rows = self.rows()
        bad = copy.deepcopy(rows); bad[0]['seed'] = 10
        for case in (rows[:-1], rows + [rows[0]], bad):
            with self.subTest(case=case[-1]):
                with self.assertRaises(ValueError):
                    audit.control_summary(case)

    def test_invalid_scores(self):
        for value in ('nan', 'inf', '-0.1', '1.1', ''):
            rows = self.rows(); rows[0]['iptm'] = value
            with self.subTest(value=value), self.assertRaises(ValueError):
                audit.control_summary(rows)

    def test_bad_timeseries(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'series.txt'
            for text in ('', '0 nan 1\n', '0 1\n', '0 1 1\n0 1 1\n', '-1 1 1\n'):
                path.write_text(text)
                with self.subTest(text=text), self.assertRaises(ValueError):
                    audit.read_series(path, 3)

    def test_committed_evidence_reproduces(self):
        report = audit.audit(ROOT)
        self.assertEqual(len(report['dual_contact']), 5)
        baseline = report['three_chain_controls']['SOX2+SOX9']['mean_iptm']
        self.assertAlmostEqual(baseline, 0.34804320335388184)
        self.assertFalse(report['three_chain_controls']['SOX11+SOX9']['all_paired_differences_positive'])


if __name__ == '__main__':
    unittest.main()
