#!/usr/bin/env python3
"""Verify odd-period LP comparisons, including all coordinate-pivot solutions.

Physical algebra and upper witnesses are checked afresh; this does not rerun
the recorded exhaustive distance exclusions. Standard-library Python only.
"""
import argparse
from itertools import combinations
import json
from pathlib import Path
import sys
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from rlf import binary_checks, echelon, from_support, in_span, reconstruct, support, translate


def affine_solutions(rows, rhs, n):
    """Enumerate the small affine solution space of rows*x=rhs over F2."""
    piv = {}
    for row, b in zip(rows, rhs):
        while row:
            p = row.bit_length()-1
            if p not in piv:
                piv[p] = (row, b)
                break
            row ^= piv[p][0]
            b ^= piv[p][1]
        if not row and b:
            return []
    free = [j for j in range(n) if j not in piv]
    if len(free) > 12:
        raise ValueError('Unexpected affine nullity')
    ans = []
    for assignment in range(1 << len(free)):
        x = sum(((assignment >> j) & 1) << f for j, f in enumerate(free))
        for p, (row, b) in sorted(piv.items()):
            if ((row & x).bit_count() & 1) != b:
                x ^= 1 << p
        if not all(((r & x).bit_count() & 1) == b for r, b in zip(rows, rhs)):
            raise ValueError('Internal affine solution check failed')
        ans.append(x)
    return ans


def optimized_pivot_width(code):
    spec = code['construction']
    P = code['group']['order']
    A = [[dict(entry) for entry in row] for row in spec['constituent']]
    r, s = len(A), len(A[0])
    outcomes = []
    for piv in combinations(range(s), r):
        rows = binary_checks([[row[j] for j in piv] for row in A], P, 2)
        widths, counts = [], []
        for j in range(s):
            if j in piv:
                continue
            rhs = [A[i][j].get(g, 0) for i in range(r) for g in range(P)]
            sol = affine_solutions(rows, rhs, r*P)
            if not sol:
                break
            widths.append(1 + min(x.bit_count() for x in sol))
            counts.append(len(sol))
        if len(widths) == s-r:
            outcomes.append(dict(pivot=list(piv), width=max(widths), affine_solution_counts=counts))
    if not outcomes:
        raise ValueError('No admissible coordinate pivot')
    return min(x['width'] for x in outcomes), outcomes


def verify_comparison(code, checks, frame, detectors):
    def require(ok, message):
        if not ok:
            raise ValueError(code['id'] + ': ' + message)
    n, k, P = code['n'], code['k'], code['group']['order']
    hx = [from_support(v, n) for v in checks['X']]
    hz = [from_support(v, n) for v in checks['Z']]
    seeds = [(from_support(v['X'], n), from_support(v['Z'], n)) for v in frame['seeds']]
    require((hx, hz, seeds) == reconstruct(code), 'physical data differ from recipe')
    require(all(not ((a & b).bit_count() & 1) for a in hx for b in hz), 'CSS commutation')
    bx, bz = echelon(hx), echelon(hz)
    require(n-len(bx)-len(bz) == k, 'dimension')
    xx = [translate(x, g, P) for x, z in seeds for g in range(P)]
    zz = [translate(z, g, P) for x, z in seeds for g in range(P)]
    require(len(xx) == len(zz) == code['partial_frame']['pairs'] == 4*P, 'partial pair count')
    require(all(not ((x & h).bit_count() & 1) for x in xx for h in hz), 'X syndrome')
    require(all(not ((z & h).bit_count() & 1) for z in zz for h in hx), 'Z syndrome')
    require(all(((x & z).bit_count() & 1) == (i == j) for i, x in enumerate(xx) for j, z in enumerate(zz)), 'partial canonical Gram')
    width = max(v.bit_count() for pair in seeds for v in pair)
    require(width == frame['width'] == code['partial_frame']['width'], 'partial width')
    require(k-len(xx) == code['partial_frame']['residual_pairs'] == 16, 'residual dimension')
    require(k % P != 0 and code['full_regular_frame_exists'] is False, 'full regular orbit divisibility obstruction')
    full = [from_support(v, n) for v in detectors['X']]
    require(len(full) == k and all(not ((x & h).bit_count() & 1) for x in full for h in hz), 'full detector syndrome/count')
    require(len(echelon(hx+full)) == len(bx)+k, 'full detector quotient rank including residual classes')

    def trans(v):
        out = 0
        for i in support(v):
            block, g = divmod(i, P)
            if block < 25:
                a, b = divmod(block, 5)
                dst = b*5+a
            else:
                a, b = divmod(block-25, 3)
                dst = 25+b*3+a
            out |= 1 << (dst*P+g)
        return out

    require(all(in_span(trans(h), bz) for h in hx) and all(in_span(trans(h), bx) for h in hz), 'sector transpose isometry')
    for witness in code['witnesses']:
        v = from_support(witness['support'], n)
        H, B = (hx, bz) if witness['sector'] == 'Z' else (hz, bx)
        require(v.bit_count() == witness['weight'] == code['distance']['upper'], 'witness weight')
        require(all(not ((v & h).bit_count() & 1) for h in H) and not in_span(v, B), 'logical witness')
    stats = {}
    for sector, H in [('X', hx), ('Z', hz)]:
        col = [sum((h >> i) & 1 for h in H) for i in range(n)]
        stats[sector] = dict(rank=len(echelon(H)), maximum_check_weight=max(map(int.bit_count, H)), maximum_column_weight=max(col))
        require(stats[sector]['maximum_check_weight'] == 8 and stats[sector]['maximum_column_weight'] == 5, 'check weights')
    best, outcomes = optimized_pivot_width(code)
    require(best == width, 'optimality within coordinate-pivot partial construction')
    return dict(id=code['id'], status='PASS', n=n, k=k, partial_pairs=len(xx), residual_pairs=16, partial_frame_width=width, full_regular_frame_exists=False, full_opposite_detector_rank=k, check_statistics=stats, coordinate_pivot_options=outcomes, distance_exhaustive_search_rerun=False)


def run(entry):
    folder = (ROOT/entry['path']).parent
    code = json.loads((ROOT/entry['path']).read_text())
    args = [json.loads((folder/code['files'][key]).read_text()) for key in ['checks', 'frame', 'detectors']]
    return verify_comparison(code, *args)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--code')
    parser.add_argument('--report', type=Path)
    a = parser.parse_args()
    entries = json.loads((ROOT/'catalog.json').read_text())['comparisons']
    if a.code:
        entries = [e for e in entries if e['id'] == a.code]
    if not entries:
        parser.error('Unknown comparison ID')
    reports = [run(e) for e in entries]
    for r in reports:
        print(r['id'] + ': PASS; partial width ' + str(r['partial_frame_width']) + '; 16 residual pairs')
    if a.report:
        a.report.write_text(json.dumps(dict(status='PASS', comparisons=reports, distance_searches_rerun=False), indent=2)+'\n')


if __name__ == '__main__':
    main()
