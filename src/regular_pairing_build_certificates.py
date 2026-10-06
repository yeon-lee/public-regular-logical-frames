#!/usr/bin/env python3
"""Regenerate small finite certificates from the preserved physical orbit pools.

This script is not the enumerator and makes no claim to prove completeness.
It uses translated dot products; check.py independently verifies polynomials
using within-block position differences.
"""
import json
from pathlib import Path
from check import HERE, coordinates, exchange, mask, parity, shift


def build(name):
    data = json.loads((HERE / 'data' / (name + '_physical.json')).read_text())
    pool = json.loads((HERE / 'data' / (name + '_pool.json')).read_text())
    p = data['P']
    zpool = list(map(mask, pool['representatives']))
    xpool = [exchange(z, p) for z in zpool]
    cx = [shift(mask(s['X_support']), p, h) for s in data['frame_seeds'] for h in range(p)]
    cz = [shift(mask(s['Z_support']), p, h) for s in data['frame_seeds'] for h in range(p)]
    gram = [[sum(parity(x, shift(z, p, h)) << h for h in range(p))
             for z in zpool] for x in xpool]
    vertices = [[i, j, entry.bit_length() - 1] for i, row in enumerate(gram)
                for j, entry in enumerate(row) if entry.bit_count() == 1]
    adjacency = [set() for _ in vertices]
    for a, (i, j, _) in enumerate(vertices):
        for b, (u, v, _) in enumerate(vertices[:a]):
            if gram[i][v] == gram[u][j] == 0:
                adjacency[a].add(b)
                adjacency[b].add(a)
    colors = {}
    while len(colors) < len(vertices):
        v = max((j for j in range(len(vertices)) if j not in colors),
                key=lambda j: (len({colors[a] for a in adjacency[j] if a in colors}),
                               len(adjacency[j]), -j))
        used = {colors[a] for a in adjacency[v] if a in colors}
        colors[v] = next(c for c in range(len(vertices)) if c not in used)
    unions = [[min((x | shift(z, p, h)).bit_count() for h in range(p))
               for z in zpool] for x in xpool]
    result = dict(code=name, n=data['n'], k=data['k'], P=p, cap=pool['cap'],
                  Z_physical_masks_hex=list(map(hex, zpool)),
                  X_physical_masks_hex=list(map(hex, xpool)),
                  logical_coordinate_convention='Bit i is pairing with canonical opposite-sector frame operator i. Index i = 16*seed + translate.',
                  Z_logical_coordinates_hex=[[hex(coordinates(shift(z, p, h), cx)) for h in range(p)] for z in zpool],
                  X_logical_coordinates_hex=[[hex(coordinates(shift(x, p, h), cz)) for h in range(p)] for x in xpool],
                  polynomial_convention='Bit h is <x_i,t^h z_j>; coefficients in F2 and h=0,...,15.',
                  pairing_polynomial_masks=gram,
                  union_convention='Entry (i,j) is min over ALL h=0,...,15 of |supp(x_i) union supp(t^h z_j)|. A common translation preserves weight.',
                  minimum_union_by_orbit_pair=unions, minimum_union=min(map(min, unions)),
                  vertex_convention='Lexicographic triples (i,j,h) such that G_ij=t^h; normalize the canonical pair by a relative translation.',
                  vertices=vertices, adjacency=[sorted(row) for row in adjacency],
                  colors=[colors[i] for i in range(len(vertices))],
                  coloring_method='DSATUR: maximum distinct neighbor colors, then larger degree, then smaller vertex index; use smallest available nonnegative color.',
                  canonical_frame_upper_width=max(v.bit_count() for v in cx + cz))
    (HERE / 'data' / (name + '_certificate.json')).write_text(json.dumps(result, indent=2) + '\n')
    print(name, len(zpool), 'orbits;', len(vertices), 'vertices;', max(colors.values()) + 1, 'colors')


if __name__ == '__main__':
    for name in ['LPQ6', 'LPQ5']:
        build(name)
