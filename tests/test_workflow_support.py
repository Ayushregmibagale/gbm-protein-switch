import importlib.util
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('md_support', ROOT/'md/workflow_support.py')
support = importlib.util.module_from_spec(spec)
spec.loader.exec_module(support)


class WorkflowTests(unittest.TestCase):
    def test_state_requires_nonempty_artifacts(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root/'state.txt').write_text('NPT_DONE')
            with self.assertRaises(RuntimeError):
                support.read_state(root)
            for name in ['npt.gro', 'npt.cpt']:
                (root/name).write_text('valid')
            self.assertEqual(support.read_state(root), 'NPT_DONE')
            (root/'state.txt').write_text('nonsense')
            with self.assertRaises(ValueError):
                support.read_state(root)

    def test_failed_stage_does_not_advance_state(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            with self.assertRaises(RuntimeError):
                support.write_state(root, 'PROD_DONE')
            self.assertFalse((root/'state.txt').exists())

    def test_command_uses_arguments_and_checks_failure(self):
        with tempfile.TemporaryDirectory(prefix='gbm space ') as directory:
            root = Path(directory)
            with patch.object(support.subprocess, 'run', side_effect=subprocess.CalledProcessError(1, 'gmx')) as runner:
                with self.assertRaises(subprocess.CalledProcessError):
                    support.run_gmx('/path with spaces/gmx', ['genion', '-f', 'file with spaces'], root/'log', root, 'SOL\n')
                args, kwargs = runner.call_args
                self.assertEqual(args[0], ['/path with spaces/gmx', 'genion', '-f', 'file with spaces'])
                self.assertTrue(kwargs['check'])
                self.assertNotIn('shell', kwargs)

    def test_checkpoint_resume_does_not_rebuild_tpr(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            for name in ['prod.tpr', 'prod.cpt', 'prod.xtc']:
                (root/name).write_text('existing')
            with patch.object(support, 'run_gmx') as run:
                support.production('gmx', root, ['grompp', 'UNEXPECTED'], 'test')
                self.assertEqual(run.call_count, 1)
                argv = run.call_args.args[1]
                self.assertIn('-cpi', argv)
                self.assertIn('-noappend', argv)
                self.assertNotIn('grompp', argv)
            self.assertEqual((root/'prod.tpr').read_text(), 'existing')

    def test_ambiguous_partial_run_is_not_overwritten(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root/'prod.tpr').write_text('old')
            with patch.object(support, 'run_gmx') as run:
                with self.assertRaises(RuntimeError):
                    support.production('gmx', root, ['grompp'], 'test')
                run.assert_not_called()

    def test_fragments_do_not_become_complete_production(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            for name in ['prod.tpr', 'prod.cpt', 'prod.part0002.xtc']:
                (root/name).write_text('old')
            with patch.object(support, 'run_gmx'):
                with self.assertRaisesRegex(RuntimeError, 'fragments'):
                    support.production('gmx', root, ['grompp'], 'test')
