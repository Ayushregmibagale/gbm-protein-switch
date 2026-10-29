import importlib.util
import unittest
from pathlib import Path

try:
    import MDAnalysis as mda
    import numpy as np
except ImportError:
    mda = None

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('periodic', ROOT / 'analysis/audit_periodic_distance.py')
periodic = importlib.util.module_from_spec(spec)
spec.loader.exec_module(periodic)


@unittest.skipIf(mda is None, 'Optional MDAnalysis dependency is not installed')
class PeriodicTests(unittest.TestCase):
    def universe(self, box):
        universe = mda.Universe.empty(2, trajectory=True)
        universe.atoms.positions = [[1, 0, 0], [99, 0, 0]]
        universe.dimensions = box
        return universe

    def test_wrapped_pair_is_close(self):
        u = self.universe([100, 100, 100, 90, 90, 90])
        row = periodic.measure_frames(u, u.atoms[:1], u.atoms[1:])[0]
        self.assertAlmostEqual(row[1], 98)
        self.assertAlmostEqual(row[2], 2, places=4)

    def test_triclinic_box(self):
        u = self.universe([100, 100, 100, 60, 60, 60])
        row = periodic.measure_frames(u, u.atoms[:1], u.atoms[1:])[0]
        self.assertAlmostEqual(row[2], 2, places=4)

    def test_missing_box_is_error(self):
        u = self.universe(None)
        with self.assertRaises(ValueError):
            periodic.measure_frames(u, u.atoms[:1], u.atoms[1:])


if __name__ == '__main__':
    unittest.main()
