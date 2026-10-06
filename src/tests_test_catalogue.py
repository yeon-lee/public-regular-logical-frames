"""Scientific regression checks for odd periods, partial frames, and shears."""
import copy
import json
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from rlf import Ring, reconstruct, verify
from comparisons.check import verify_comparison
from codes.reproduce_distance import completed


def load(identifier, comparison=False):
    cat = json.loads((ROOT/'catalog.json').read_text())
    entry = next(e for e in cat['comparisons' if comparison else 'entries'] if e['id'] == identifier)
    path = ROOT/entry['path']
    code = json.loads(path.read_text())
    keys = ['checks', 'frame', 'detectors'] if comparison else ['checks', 'frame']
    return (code, *(json.loads((path.parent/code['files'][k]).read_text()) for k in keys))


class CatalogueTests(unittest.TestCase):
    def test_incomplete_distance_runs_are_rejected(self):
        terminal = 'COMPLETE W=9 roots=[0,2) children=[0,1073741824) mod(D=-1,M=1,r=0) nodes=50 t=0.1\n'
        roots = 'root block 0 done: nodes 20\nroot block 1 done: nodes 30\n'
        self.assertTrue(completed(terminal, roots, 9, 2, 0))
        self.assertFalse(completed('TIMEOUT nodes=50 t=0.1\n', roots, 9, 2, 0))
        self.assertFalse(completed(terminal, roots.splitlines()[0], 9, 2, 0))
        self.assertFalse(completed(terminal, roots, 9, 2, -9))
        self.assertFalse(completed(terminal.replace('M=1', 'M=2'), roots, 9, 2, 0))

    def test_monomial_inverse_odd_period(self):
        ring = Ring(7)
        for coefficient in (1, 2, 3):
            polynomial = {3: coefficient}
            self.assertEqual(ring.mul(polynomial, ring.inverse(polynomial)), {0: 1})
        with self.assertRaises(ValueError):
            ring.inverse({0: 1, 1: 1})

    def test_odd_period_complete_quaternary_frames(self):
        for d in (12, 13):
            data = load(f'lp-f4-476-56-{d}')
            self.assertEqual(verify(*data)['status'], 'PASS')

    def test_pp_paired_shear_improves_and_is_involutory(self):
        code, checks, frame = load('pp-f4-1024-256-18-w9')
        self.assertEqual(verify(code, checks, frame)['frame_width'], 26)
        raw = copy.deepcopy(code)
        operation = raw['construction'].pop('paired_binary_shears')
        original = reconstruct(raw)
        self.assertEqual(max(v.bit_count() for pair in original[2] for v in pair), 28)
        twice = copy.deepcopy(code)
        twice['construction']['paired_binary_shears'] = operation*2
        self.assertEqual(reconstruct(twice), original)
        with self.assertRaisesRegex(ValueError, 'construction recipe'):
            verify(raw, checks, frame)

    def test_partial_frame_residuals_are_checked(self):
        code, checks, frame, detectors = load('lp-f2-1122-148-20-partial', True)
        report = verify_comparison(code, checks, frame, detectors)
        self.assertEqual(report['partial_frame_width'], 44)
        self.assertEqual(report['residual_pairs'], 16)
        damaged = copy.deepcopy(detectors)
        damaged['X'] = damaged['X'][:-16]
        with self.assertRaisesRegex(ValueError, 'full detector'):
            verify_comparison(code, checks, frame, damaged)
        damaged = copy.deepcopy(code)
        damaged['partial_frame']['width'] = 12
        with self.assertRaisesRegex(ValueError, 'partial width'):
            verify_comparison(damaged, checks, frame, detectors)

    def test_public_catalogue_coverage(self):
        cat = json.loads((ROOT/'catalog.json').read_text())
        route = json.loads((ROOT/'docs/data-map.json').read_text())
        entries = cat['entries'] + cat['comparisons']
        self.assertEqual(len(entries), len({e['id'] for e in entries}))
        self.assertEqual({e['id'] for e in entries}, {e['id'] for e in route['codes']})
        self.assertEqual(len(entries), route['entry_count'])
        for entry in route['codes']:
            self.assertIn(entry['paper_status'], ('current', 'current-comparison', 'supplementary'))
            code = json.loads((ROOT/entry['code_record']).read_text())
            self.assertEqual(entry['id'], code['id'])
            for path in entry['files'].values():
                self.assertTrue((ROOT/path).is_file(), path)
            for path in entry['evidence']['referenced_paths']:
                self.assertTrue((ROOT/path).is_file(), path)


if __name__ == '__main__':
    unittest.main()
