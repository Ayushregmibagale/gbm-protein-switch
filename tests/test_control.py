import importlib.util
from pathlib import Path
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('control', ROOT/'analysis/audit_control.py')
control = importlib.util.module_from_spec(spec)
spec.loader.exec_module(control)


class CoverageTests(unittest.TestCase):
    def test_complete_interval(self):
        control.check_times([0, .01, .02], end=.02)

    def test_partial_overlap_gap_and_nonfinite_rejected(self):
        for times in [[.01,.02], [0,.01], [0,0,.02], [0,.02], [0,float('nan'),.02], [0,.005,.02]]:
            with self.subTest(times=times), self.assertRaises(ValueError):
                control.check_times(times, end=.02)

    def test_real_fragment_reader_and_periodic_distance(self):
        try:
            import MDAnalysis as mda
        except ImportError:
            self.skipTest('Optional MDAnalysis dependency is absent')
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            top = root/'top.pdb'
            top.write_text('ATOM      1  CA  ALA A   1       1.000   0.000   0.000  1.00  0.00           C\n'
                           'ATOM      2  CA  ALA B   2      99.000   0.000   0.000  1.00  0.00           C\nEND\n')
            universe = mda.Universe(str(top))
            universe.dimensions = [100,100,100,90,90,90]
            def fragment(name, times):
                path = root/name
                with mda.Writer(str(path), n_atoms=2) as writer:
                    for time in times:
                        universe.trajectory.ts.time = time
                        writer.write(universe.atoms)
                return path
            a = fragment('a.xtc', [0,10,20])
            b = fragment('b.xtc', [20,30])
            with self.assertRaises(ValueError):
                control.analyse_control(top, [a,b], 'index 0','index 1',root/'bad', end=.03)
            self.assertFalse((root/'bad').exists())
            result = control.analyse_control(top,[a,b],'index 0','index 1',root/'good',end=.03,first_fragment_end=.01)
            self.assertEqual(result['n_frames'], 4)
            self.assertAlmostEqual(result['mean_dist_A'], 2, places=4)
            self.assertFalse(result['distance_criterion_met'])
            self.assertEqual(result['sources'][0]['accepted_frames'], 2)
            self.assertFalse(list(root.glob('.*offset*')))
