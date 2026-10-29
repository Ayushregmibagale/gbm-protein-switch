import importlib.util
import json
from pathlib import Path
import sys
import unittest

ROOT=Path(__file__).resolve().parents[1]


class CohortTests(unittest.TestCase):
    def test_matched_controls_have_complete_verified_exports(self):
        import csv,hashlib
        root=ROOT/'evidence/matched_controls'
        with (root/'results.tsv').open() as handle:rows=list(csv.DictReader(handle,delimiter='\t'))
        self.assertEqual(len(rows),12)
        self.assertEqual({row['condition'] for row in rows},{'BOTH','SOX2_ONLY','SOX9_ONLY','NEITHER'})
        for name in {row['condition'] for row in rows}:
            self.assertEqual({row['seed'] for row in rows if row['condition']==name},{'0','1','2'})
        report=json.loads((root/'export_manifest.json').read_text())
        for name,checksum in report['files_sha256'].items():
            self.assertEqual(hashlib.sha256((root/name).read_bytes()).hexdigest(),checksum)

    def test_prediction_model_sequences_match_archived_inputs(self):
        try:
            import yaml
            from Bio.PDB import PDBParser
            from Bio.Data.PDBData import protein_letters_3to1
        except ImportError:self.skipTest('Optional structural dependencies are absent')
        import csv
        parser=PDBParser(QUIET=True)
        for directory,key in [('corrected_candidate_predictions','name'),('matched_controls','condition')]:
            root=ROOT/'evidence'/directory
            with (root/'results.tsv').open() as handle:rows=list(csv.DictReader(handle,delimiter='\t'))
            for row in rows:
                proteins=yaml.safe_load((root/'inputs'/f"{row[key]}.yaml").read_text())['sequences']
                expected={p['protein']['id']:p['protein']['sequence'] for p in proteins}
                model=parser.get_structure(row[key],root/row['structure_path'])[0]
                actual={c.id:''.join(protein_letters_3to1[r.resname] for r in c if r.id[0]==' ') for c in model}
                self.assertEqual(actual,expected)

    def test_corrected_predictions_have_complete_seeds_and_valid_hashes(self):
        import csv,hashlib
        root=ROOT/'evidence/corrected_candidate_predictions'
        with (root/'results.tsv').open() as handle:rows=list(csv.DictReader(handle,delimiter='\t'))
        self.assertEqual(len(rows),69)
        self.assertEqual(len({row['name'] for row in rows}),23)
        for name in {row['name'] for row in rows}:
            self.assertEqual({row['seed'] for row in rows if row['name']==name},{'0','1','2'})
        for row in rows:
            for path_key,hash_key in [('confidence_path','confidence_sha256'),('structure_path','structure_sha256')]:
                self.assertEqual(hashlib.sha256((root/row[path_key]).read_bytes()).hexdigest(),row[hash_key])

    def test_legacy_shell_launchers_are_inactive_by_default(self):
        import os,subprocess
        env=os.environ.copy();env.pop('GBM_ALLOW_LEGACY_CAMPAIGNS',None)
        for pattern in ['run_*.sh','bindcraft/*.sh','boltzgen/*.sh','protein_hunter/*.sh','rfdiffusion/*.sh','mpnn/*.sh']:
            for path in ROOT.glob(pattern):
                result=subprocess.run(['bash',str(path)],capture_output=True,text=True,env=env)
                self.assertEqual(result.returncode,2,str(path))
                self.assertIn('Archived campaign launcher',result.stderr)

    def test_unbound_control_has_complete_spliced_coverage(self):
        import hashlib
        root=ROOT/'evidence/unbound_control'
        report=json.loads((root/'metrics.json').read_text())
        self.assertEqual(report['n_frames'],2001)
        self.assertEqual(report['elapsed_ns'],20)
        self.assertEqual(sum(source['accepted_frames'] for source in report['sources']),2001)
        self.assertEqual(report['pct_above_25A'],100)
        self.assertEqual(hashlib.sha256((root/'distances.tsv').read_bytes()).hexdigest(),report['timeseries_sha256'])

    def test_replicates_are_summarised_as_replicates(self):
        rows=json.loads((ROOT/'evidence/bridge_replicate_cohort/cohort.json').read_text())['sensitivity']
        self.assertEqual(len(rows),5*3*3*3)
        selected=[r for r in rows if r['discard_ns']==0 and r['cutoff_A']==8 and r['fraction_threshold']==.10]
        self.assertEqual(len(selected),5)
        for row in selected:
            self.assertEqual(row['n_replicates'],3)
            self.assertEqual(len(row['replicate_values_pct']),3)
            self.assertAlmostEqual(row['replicate_mean_pct'],sum(row['replicate_values_pct'])/3)
        unstable=next(row for row in selected if row['design']=='bridge_90_s01')
        self.assertLess(unstable['replicate_min_pct'],20)
        self.assertGreater(unstable['replicate_max_pct'],99)

    def test_geometry_fit_is_rigid_transform_invariant(self):
        try:
            import numpy as np
        except ImportError:self.skipTest('Optional numpy dependency is absent')
        sys.path.insert(0,str(ROOT/'analysis'))
        spec=importlib.util.spec_from_file_location('geometry',ROOT/'analysis/audit_native_geometry.py')
        module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
        coords=np.array([[0,0,0],[1,0,0],[0,1,0],[0,0,1]],dtype=float)
        rotation=np.array([[0,-1,0],[1,0,0],[0,0,1]])
        self.assertAlmostEqual(module.rmsd_fit(coords@rotation+10,coords),0,places=10)
        stretched=coords.copy();stretched[3]*=10
        self.assertGreater(module.rmsd_fit(stretched,coords),1)

    def test_native_mapping_and_outputs_are_finite(self):
        for name in ['lead_construct','rep2','rep3']:
            report=json.loads((ROOT/'evidence/native_geometry'/name/'geometry.json').read_text())
            for chain,mapping in report['reference_mapping'].items():
                self.assertGreater(len(mapping),200)
                numbers=[row[0] for row in mapping]
                self.assertEqual(len(numbers),len(set(numbers)))
                self.assertTrue(all(140<=n<=416 for n in numbers))
            self.assertEqual(len(report['observations']),402)

    def test_cache_inventory_does_not_claim_verified_provenance(self):
        report=json.loads((ROOT/'evidence/historical_cache_audit.json').read_text())
        self.assertEqual(report['n_runs'],len(report['runs']))
        self.assertIn('prevent verified cache reuse',report['scope'])
