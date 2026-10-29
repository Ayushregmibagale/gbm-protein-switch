import csv
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]


class CandidateTests(unittest.TestCase):
    def test_quarantine_retains_no_mismatched_names(self):
        report = json.loads((ROOT/'evidence/stage1_chain_identity.json').read_text())
        invalid = {row['name'] for row in report['checked'] if not row['recorded_sequence_matches_chain_B']}
        self.assertEqual(len(invalid), 23)
        self.assertTrue(all(row['recorded_sequence_matches_chain_A'] for row in report['checked'] if row['name'] in invalid))
        with (ROOT/'selection/validated_stage1_results.tsv').open() as handle:
            retained = list(csv.DictReader(handle,delimiter='\t'))
        self.assertTrue(retained)
        self.assertFalse(invalid & {row['name'] for row in retained})

    def test_candidate_report_reproduces_from_accepted_structures(self):
        try:
            import Bio
        except ImportError:
            self.skipTest('Optional Biopython dependency is absent')
        spec = importlib.util.spec_from_file_location('candidate',ROOT/'analysis/audit_candidate_identity.py')
        module = importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
        self.assertEqual(module.audit(ROOT),json.loads((ROOT/'evidence/stage1_chain_identity.json').read_text()))

    def test_campaign_help_cannot_launch_inference(self):
        result = subprocess.run([sys.executable,str(ROOT/'selection/stage5_flipped.py'),'--help'],capture_output=True,text=True)
        self.assertEqual(result.returncode,0)
        self.assertIn('--execute', result.stdout)
        result = subprocess.run([sys.executable,str(ROOT/'selection/stage5_flipped.py')],capture_output=True,text=True)
        self.assertNotEqual(result.returncode,0)
        self.assertIn('Use --execute',result.stderr)

    def test_legacy_inference_requires_fresh_campaign_root(self):
        import os
        env=os.environ.copy()
        env.pop('GBM_RUN_ROOT',None)
        result=subprocess.run([sys.executable,str(ROOT/'selection/stage5_flipped.py'),'--execute'],
                              capture_output=True,text=True,env=env)
        self.assertNotEqual(result.returncode,0)
        self.assertIn('Set GBM_RUN_ROOT',result.stderr)
