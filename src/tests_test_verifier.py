"""Regression checks for rejection of corrupt physical data and field arithmetic."""
import ast
import copy
import json
import os
from pathlib import Path
import subprocess
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from rlf import Ring, binary_vector, reconstruct, support, verify


class VerificationTests(unittest.TestCase):
    def setUp(self):
        self.path = ROOT / 'counterexamples/lp-binary/lp-f2-20-4-3'
        self.code = json.loads((self.path / 'code.json').read_text())
        self.checks = json.loads((self.path / 'checks.json').read_text())
        self.frame = json.loads((self.path / 'frames/paper.json').read_text())

    def test_small_lp_analytical_distance(self):
        hx, hz, _ = reconstruct(self.code)
        # Nonzero, distinct syndrome columns exclude every weight-one/two word.
        for rows in (hx, hz):
            cols = [sum(((row >> j) & 1) << i for i, row in enumerate(rows)) for j in range(20)]
            self.assertNotIn(0, cols)
            self.assertEqual(len(set(cols)), 20)
        self.assertEqual(verify(self.code, self.checks, self.frame)['status'], 'PASS')

    def test_corrupt_seed_rejected(self):
        changed = copy.deepcopy(self.frame)
        changed['seeds'][0]['Z'] = []
        with self.assertRaisesRegex(ValueError, 'construction recipe'):
            verify(self.code, self.checks, changed)

    def test_invalid_witness_rejected(self):
        changed = copy.deepcopy(self.code)
        changed['witnesses'][0]['support'] = [0]
        with self.assertRaisesRegex(ValueError, 'distance witness'):
            verify(changed, self.checks, self.frame)

    def test_corrupt_check_rejected(self):
        changed = copy.deepcopy(self.checks)
        changed['X'][0] = []
        with self.assertRaisesRegex(ValueError, 'construction recipe'):
            verify(self.code, changed, self.frame)

    def test_trace_basis_and_square_zero_final_lift(self):
        ring = Ring(64)
        self.assertEqual(binary_vector([{0: 1}], 64, 4), (1 << 64) + 1)
        self.assertEqual(binary_vector([{0: 2}], 64, 4), 1)
        self.assertEqual(binary_vector([{0: 3}], 64, 4), 1 << 64)
        self.assertEqual(ring.mul({18: 1, 50: 1}, {18: 1, 50: 1}), {})
        self.assertEqual(ring.mul({18: 1, 36: 1, 50: 1}, {18: 1, 36: 1, 50: 1}), {8: 1})

    def test_nonmonomial_pivot_inverse(self):
        ring = Ring(16)
        determinant = {e: 1 for e in support(3023)}
        inverse = {e: 1 for e in support(60408)}
        self.assertEqual(ring.inverse(determinant), inverse)

    def test_optimized_python_rejected(self):
        # Cover every distinct assertion-dependent module, including direct
        # checker calls that bypass the release's top-level entry point.
        targets = [ROOT / 'verify_all.py']
        seen = set()
        for path in sorted(ROOT.rglob('*.py')):
            source = path.read_text()
            if source in seen:
                continue
            seen.add(source)
            if any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source))):
                targets.append(path)
        self.assertGreater(len(targets), 1)
        for flags, optimization in [(['-O'], None), (['-OO'], None), ([], '1'), ([], '2')]:
            env = os.environ.copy()
            env.pop('PYTHONOPTIMIZE', None)
            env['PYTHONDONTWRITEBYTECODE'] = '1'
            if optimization is not None:
                env['PYTHONOPTIMIZE'] = optimization
            for path in targets:
                with self.subTest(path=str(path.relative_to(ROOT)), flags=flags,
                                  optimization=optimization):
                    result = subprocess.run([sys.executable, '-B', *flags, str(path)],
                                            env=env, capture_output=True, text=True,
                                            timeout=10)
                    self.assertNotEqual(result.returncode, 0)
                    self.assertIn('Verification requires Python assertions', result.stderr)


if __name__ == '__main__':
    unittest.main()
