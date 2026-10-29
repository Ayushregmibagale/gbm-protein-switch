"""Exercise original failure cases without running top-level structure analyses."""
import ast
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def function(path, name, namespace=None):
    tree = ast.parse(path.read_text())
    node = next(n for n in ast.walk(tree) if isinstance(n, ast.FunctionDef) and n.name == name)
    scope = {} if namespace is None else namespace
    exec(compile(ast.Module(body=[node], type_ignores=[]), str(path), 'exec'), scope)
    return scope[name]


class RegressionTests(unittest.TestCase):
    def test_first_patch_is_not_duplicated(self):
        patches = function(ROOT / 'analysis/01_analyse_cocomplex.py', 'find_patches')
        self.assertEqual(patches([{'resnum': 1}, {'resnum': 2}, {'resnum': 20}]), [[1, 2], [20]])
        self.assertEqual(patches([]), [])
        self.assertEqual(patches([{'resnum': 1}]), [[1]])

    def test_chain_length_with_restarted_numbering(self):
        import tempfile
        infer = function(ROOT / 'analysis/bridge_dual_contact.py', 'infer_bridge_len', {'Path': Path})
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)
            lines = [f'ATOM  {i:5d}  CA  ALA C{i:4d}    0.000   0.000   0.000\n' for i in (1, 2, 3)]
            (path / 'start.pdb').write_text(''.join(lines))
            self.assertEqual(infer(path), 3)
            (path / 'start.pdb').write_text('')
            with self.assertRaises(ValueError):
                infer(path)

    def test_residue_coordinates_keep_requested_order_and_reject_missing(self):
        try:
            import numpy as np
            from Bio.PDB import Atom, Chain, Residue
        except ImportError:
            self.skipTest('Optional structural dependencies are absent')
        coords = function(ROOT / 'analysis/02_epitope_geometry.py', 'ca_coords', {'np': np})
        chain = Chain.Chain('A')
        for i in (2, 1):
            res = Residue.Residue((' ', i, ' '), 'ALA', '')
            res.add(Atom.Atom('CA', np.array([i, 0, 0]), 0, 1, ' ', ' CA ', i, element='C'))
            chain.add(res)
        np.testing.assert_equal(coords(chain, [1, 2])[:, 0], [1, 2])
        with self.assertRaises(ValueError):
            coords(chain, [1, 3])
