import hashlib
import json
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[1]


class ProvenanceTests(unittest.TestCase):
    def test_saved_periodic_series_match_hashes(self):
        for path in sorted((ROOT/'evidence').glob('track1_*_periodic.json')):
            report = json.loads(path.read_text())
            self.assertEqual(hashlib.sha256(path.with_suffix('.tsv').read_bytes()).hexdigest(),
                             report['timeseries_sha256'])
            self.assertEqual(report['n_frames'], 2001)
            self.assertEqual(report['elapsed_ns'], 20)

    def test_bridge_exports_match_manifest(self):
        folder = ROOT/'evidence/bridge_periodic'
        report = json.loads((folder/'provenance.json').read_text())
        self.assertEqual(len(report['source_files']), 5)
        for name, digest in report['outputs'].items():
            self.assertEqual(hashlib.sha256((folder/name).read_bytes()).hexdigest(), digest)

    def test_recorded_construct_inputs_match_provenance(self):
        sources = json.loads((ROOT/'evidence/construct_input_sources.json').read_text())
        for source in sources:
            self.assertEqual(hashlib.sha256((ROOT/source['input']).read_bytes()).hexdigest(), source['sha256'])
