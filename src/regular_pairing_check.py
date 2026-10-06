#!/usr/bin/env python3
"""Small finite checker; this does not replace the exhaustive pool enumeration.

No search-engine code is imported. Python >=3.10, standard library only.
Checks raise exceptions explicitly, so Python optimization cannot disable them.
"""
import collections
import hashlib
import json
from pathlib import Path
import re

HERE = Path(__file__).resolve().parent


def require(condition, message):
    if not condition:
        raise ValueError(message)


def mask(support):
    require(len(support) == len(set(support)), 'duplicate coordinate')
    return sum(1 << q for q in support)


def support(value):
    result = []
    while value:
        bit = value & -value
        result.append(bit.bit_length() - 1)
        value ^= bit
    return result


def rank(rows):
    pivots = {}
    for row in rows:
        while row:
            p = row.bit_length() - 1
            if p not in pivots:
                pivots[p] = row
                break
            row ^= pivots[p]
    return len(pivots)


def parity(x, y):
    return (x & y).bit_count() % 2


def shift(value, period, h):
    return sum(1 << (q // period * period + (q % period + h) % period)
               for q in support(value))


def exchange(value, period, s=5, r=3):
    """Transpose each quaternary tensor block; retain the trace component."""
    result = 0
    for q in support(value):
        block, tail = divmod(q, 2 * period)
        if block < s * s:
            a, b = divmod(block, s)
            block = b * s + a
        else:
            a, b = divmod(block - s * s, r)
            block = s * s + b * r + a
        result |= 1 << (block * 2 * period + tail)
    return result


def read_input(path):
    lines = path.read_text().splitlines()
    n, p, blocks = map(int, lines[0].split())
    pos, groups = 1, []
    for _ in range(3):
        count = int(lines[pos])
        pos += 1
        group = [list(map(int, row.split())) for row in lines[pos:pos + count]]
        pos += count
        require(len(group) == count, 'truncated input')
        require(all(0 <= q < n for row in group for q in row), 'coordinate range')
        groups.append(list(map(mask, group)))
    require(pos == len(lines), 'unexpected input suffix')
    return n, p, blocks, groups


def coordinates(value, dual_basis):
    return sum(parity(value, dual) << i for i, dual in enumerate(dual_basis))


def polynomial_by_differences(x, z, period):
    """Coefficient h is <x,t^h z>; compute by position differences."""
    result = 0
    for a in support(x):
        for b in support(z):
            if a // period == b // period:
                result ^= 1 << ((a - b) % period)
    return result


def reconstruct_checks(data):
    """Build LP checks directly from the supplied F4 coefficient/exponent arrays.

    rho(c) represents multiplication in the self-dual trace basis (omega,omega^2).
    The constituent is A_ij = C_ij*t^E_ij. No search or exported check matrix is
    used in this reconstruction.
    """
    rho = [((0, 0), (0, 0)), ((1, 0), (0, 1)),
           ((0, 1), (1, 1)), ((1, 1), (1, 0))]
    c, e, p = data['C'], data['E'], data['P']
    r, s = len(c), len(c[0])
    hx, hz = [], []
    for i in range(r):
        for b in range(s):
            hx.append([(a*s+b, c[i][a], e[i][a]) for a in range(s)]
                      + [(s*s+i*r+j, c[j][b], -e[j][b]) for j in range(r)])
    for a in range(s):
        for j in range(r):
            hz.append([(a*s+b, c[j][b], e[j][b]) for b in range(s)]
                      + [(s*s+i*r+j, c[i][a], -e[i][a]) for i in range(r)])

    def expand(rows):
        raw, balanced = [], []
        for monomials in rows:
            first = [0, 0]
            for block, coeff, exponent in monomials:
                require(0 <= coeff <= 3, 'F4 coefficient range')
                for a in range(2):
                    for b in range(2):
                        if rho[coeff][a][b]:
                            first[a] ^= 1 << ((2*block+b)*p + (-exponent) % p)
            raw.extend(shift(v, p, h) for v in first for h in range(p))
            three = [first[0], first[1], first[0] ^ first[1]]
            chosen = sorted(range(3), key=lambda i: (three[i].bit_count(), i))[:2]
            balanced.extend(shift(three[i], p, h) for i in chosen for h in range(p))
        return raw, balanced
    return expand(hx), expand(hz)


def load_code(name):
    data = json.loads((HERE / 'data' / (name + '_physical.json')).read_text())
    n, k, p = data['n'], data['k'], data['P']
    require((n, k, p) == (1088, 128, 16), 'unexpected code parameters')
    for key in ['HX_balanced', 'HZ_balanced']:
        require(all(0 <= q < n for row in data[key] for q in row), 'coordinate range')
    hx = list(map(mask, data['HX_balanced']))
    hz = list(map(mask, data['HZ_balanced']))
    for key, (raw, balanced) in zip(['HX', 'HZ'], reconstruct_checks(data)):
        require(raw == list(map(mask, data[key + '_raw'])), key + ' constituent reconstruction')
        require(balanced == list(map(mask, data[key + '_balanced'])), key + ' balanced presentation')
    require(rank(hx) == rank(hz) == 480, 'stabilizer ranks')
    require(all(not parity(x, z) for x in hx for z in hz), 'CSS commutation')
    require(n - rank(hx) - rank(hz) == k, 'encoded dimension')
    require(rank(hx + [shift(x, p, 1) for x in hx]) == 480, 'X translation symmetry')
    require(rank(hz + [shift(z, p, 1) for z in hz]) == 480, 'Z translation symmetry')
    require(rank(hz + [exchange(x, p) for x in hx]) == 480, 'sector exchange')
    require(rank(hx + [exchange(z, p) for z in hz]) == 480, 'inverse sector exchange')
    require(all(0 <= q < n for row in data['frame_seeds']
                for key in ['X_support', 'Z_support'] for q in row[key]), 'canonical coordinate range')
    seedx = [mask(row['X_support']) for row in data['frame_seeds']]
    seedz = [mask(row['Z_support']) for row in data['frame_seeds']]
    xbasis = [shift(x, p, h) for x in seedx for h in range(p)]
    zbasis = [shift(z, p, h) for z in seedz for h in range(p)]
    require(len(xbasis) == len(zbasis) == k, 'canonical basis size')
    require(all(not parity(x, z) for x in hx for z in zbasis), 'canonical Z syndrome')
    require(all(not parity(x, z) for x in xbasis for z in hz), 'canonical X syndrome')
    require(all(parity(x, z) == (i == j)
                for i, x in enumerate(xbasis) for j, z in enumerate(zbasis)), 'canonical Gram')
    return data, hx, hz, xbasis, zbasis


def check_execution(name, cap, normalized):
    """Bind released pool to complete transcripts; do not certify the traversal."""
    inp = HERE / 'inputs' / (name + '.txt')
    engine = HERE / 'src' / ('kernel_minima_reverse.cpp' if name == 'LPQ6' else 'kernel_minima_native.cpp')
    sha = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()
    summary = []
    for era in ['records', 'replay']:
        directory = HERE / 'records' / era
        record = json.loads((directory / (name + '.json')).read_text())
        require(record['input_sha256'] == sha(inp), era + ' input hash')
        require(record.get('source_sha256', record.get('engine_sha256')) == sha(engine), era + ' engine hash')
        raw = directory / (name + '.stdout')
        lines = raw.read_text().splitlines()
        require(lines[-1] == record['terminal'], era + ' terminal record')
        pattern = (rf'ENUMERATION_COMPLETE W={cap} roots=\[0,68\) '
                   r'children=\[0,1073741824\) mod\(D=-1,M=1,r=0\) nodes=(\d+) t=([\d.]+)')
        match = re.fullmatch(pattern, lines[-1])
        require(match is not None, era + ' complete, unpartitioned root coverage')
        words = set()
        for line in lines[:-1]:
            require(line.startswith('WORD '), era + ' unexpected output')
            values = list(map(int, line.split()[1:]))
            require(values[0] == len(values[1:]) <= cap, era + ' word length')
            value = mask(values[1:])
            words.add(min(shift(value, 16, h) for h in range(16)))
        require(words == set(normalized), era + ' output pool comparison')
        if era == 'replay':
            require(record['complete'] and record['returncode'] == 0, 'replay return state')
            require(record['stdout_sha256'] == sha(raw), 'replay stdout hash')
            require(record['stderr_sha256'] == sha(directory / (name + '.stderr')), 'replay stderr hash')
        summary.append(dict(record=era, nodes=int(match.group(1)), reported_seconds=float(match.group(2))))
    return summary


def check(name):
    data, hx, hz, canonical_x, canonical_z = load_code(name)
    p, k, n = data['P'], data['k'], data['n']
    pool = json.loads((HERE / 'data' / (name + '_pool.json')).read_text())
    cert = json.loads((HERE / 'data' / (name + '_certificate.json')).read_text())
    cap = pool['cap']
    n0, p0, blocks, (checks, stabilizers, detectors) = read_input(HERE / 'inputs' / (name + '.txt'))
    require((n0, p0, blocks) == (n, p, n // p), 'search dimensions')
    require(rank(hx + checks) == 480 == rank(checks), 'search check row space')
    require(rank(hz + stabilizers) == 480, 'search stabilizer row space')
    require(all(not parity(x, z) for x in detectors for z in hz), 'detector syndrome')
    require(rank(hx + detectors) == 480 + k, 'complete logical detectors')
    require(all(any(row >> q & 1 for row in checks) for q in range(n)), 'no zero columns')
    if name == 'LPQ6':
        require(not stabilizers, 'Q6 must disable stabilizer-overlap pruning: rZ = 0')
    require(all(0 <= q < n for row in pool['representatives'] for q in row), 'pool coordinate range')
    zpool = list(map(mask, pool['representatives']))
    require(zpool == sorted(set(zpool)), 'normalized mask order and uniqueness')
    require(all(z == min(shift(z, p, h) for h in range(p)) for z in zpool), 'orbit normalization')
    xpool = [exchange(z, p) for z in zpool]
    for sector, values, checks_, dual in [('Z', zpool, hx, canonical_x), ('X', xpool, hz, canonical_z)]:
        all_words = [shift(value, p, h) for value in values for h in range(p)]
        require(len(set(all_words)) == p * len(values), sector + ' free physical orbits')
        require(all(0 < value.bit_count() <= cap for value in all_words), sector + ' weight cap')
        require(all(not parity(value, row) for value in all_words for row in checks_), sector + ' syndrome')
        labels = [[hex(coordinates(shift(value, p, h), dual)) for h in range(p)] for value in values]
        require(all(int(label, 16) for row in labels for label in row), sector + ' nontriviality')
        require(cert[sector + '_logical_coordinates_hex'] == labels, sector + ' logical coordinates')
    require(cert['Z_physical_masks_hex'] == list(map(hex, zpool)), 'physical Z masks')
    require(cert['X_physical_masks_hex'] == list(map(hex, xpool)), 'physical X masks')
    gram = [[polynomial_by_differences(x, z, p) for z in zpool] for x in xpool]
    require(cert['pairing_polynomial_masks'] == gram, 'pairing polynomials')
    # Translation of both supports preserves the union size. Thus one relative
    # shift is enough, but EVERY i, j and h must be tested.
    union_minima = [[min((x | shift(z, p, h)).bit_count() for h in range(p))
                    for z in zpool] for x in xpool]
    require(cert['minimum_union_by_orbit_pair'] == union_minima, 'all relative-shift support unions')
    union_min = min(map(min, union_minima))
    require(union_min == cert['minimum_union'] > cap, 'mixed-component exclusion')
    vertices = [[i, j, entry.bit_length() - 1] for i, row in enumerate(gram)
                for j, entry in enumerate(row) if entry.bit_count() == 1]
    require(cert['vertices'] == vertices, 'complete monomial-pair vertex set')
    adjacency = [[] for _ in vertices]
    for a, (i, j, _) in enumerate(vertices):
        for b in range(a):
            u, v, _ = vertices[b]
            if gram[i][v] == gram[u][j] == 0:
                adjacency[a].append(b)
                adjacency[b].append(a)
    adjacency = [sorted(row) for row in adjacency]
    require(cert['adjacency'] == adjacency, 'compatibility graph')
    colors = cert['colors']
    require(len(colors) == len(vertices), 'color count')
    require(all(isinstance(c, int) and c >= 0 for c in colors), 'nonnegative integer colors')
    require(all(colors[a] != colors[b] for a, neighbors in enumerate(adjacency) for b in neighbors), 'proper coloring')
    color_count = max(colors) + 1
    require(color_count < k // p, 'coloring excludes eight-clique')
    upper = max(value.bit_count() for value in canonical_x + canonical_z)
    require(upper == cert['canonical_frame_upper_width'], 'canonical upper bound')
    expected = {'LPQ6': (45, 26, 177, 5, 18), 'LPQ5': (75, 32, 88, 3, 22)}[name]
    require((len(zpool), union_min, len(vertices), color_count, upper) == expected, 'reported quantities')
    separate = json.loads((HERE / 'data' / (name + '_sector_bases.json')).read_text())
    require(all(0 <= q < n for row in separate['X_seeds'] + separate['Z_seeds'] for q in row),
            'separate sector-basis coordinate range')
    minimum_weight = {'LPQ6': 16, 'LPQ5': 20}[name]
    for sector, checks_, stabilizers in [('X', hz, hx), ('Z', hx, hz)]:
        basis_ = [shift(mask(seed), p, h) for seed in separate[sector + '_seeds'] for h in range(p)]
        require(len(basis_) == k and all(v.bit_count() == minimum_weight for v in basis_),
                sector + ' minimum-weight regular basis')
        require(all(not parity(v, c) for v in basis_ for c in checks_), sector + ' basis syndrome')
        require(rank(stabilizers + basis_) == rank(stabilizers) + k, sector + ' complete logical span')
    executions = check_execution(name, cap, zpool)
    return dict(code=name, checked=True, physical_orbits=len(zpool),
                relative_shift_unions_checked=len(zpool) ** 2 * p,
                minimum_union=union_min, vertices=len(vertices), colors=color_count,
                regular_width_lower=cap + 1, regular_width_upper=upper, executions=executions,
                scope='Finite consequences checked independently. Pool completeness still relies on the released exhaustive enumeration algorithm and completed runs. No unrestricted-width lower bound above distance is inferred.')


if __name__ == '__main__':
    print(json.dumps([check(name) for name in ['LPQ6', 'LPQ5']], indent=2))
