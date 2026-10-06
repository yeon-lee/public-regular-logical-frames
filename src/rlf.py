"""Reconstruct and verify the manuscript's cyclic CSS codes (Python 3.10+).

Only the Python standard library is required. Polynomials are sparse maps from
exponent to coefficient; coefficients 0,1,2,3 mean 0,1,omega,omega^2.
"""
from itertools import permutations

MUL = ((0, 0, 0, 0), (0, 1, 2, 3), (0, 2, 3, 1), (0, 3, 1, 2))
COORD = (0, 3, 1, 2)  # coordinates in the self-dual basis (omega, omega^2)


class Ring:
    def __init__(self, period):
        self.P = period

    def add(self, *args):
        out = {}
        for a in args:
            for e, c in a.items():
                e %= self.P
                out[e] = out.get(e, 0) ^ c
        return {e: c for e, c in out.items() if c}

    def mul(self, a, b):
        out = {}
        for e, c in a.items():
            for f, d in b.items():
                g = (e + f) % self.P
                out[g] = out.get(g, 0) ^ MUL[c][d]
        return {e: c for e, c in out.items() if c}

    def bar(self, a):
        return {(-e) % self.P: c for e, c in a.items()}

    def det(self, a):
        out = {}
        for perm in permutations(range(len(a))):
            term = {0: 1}
            for i, j in enumerate(perm):
                term = self.mul(term, a[i][j])
            out = self.add(out, term)
        return out

    def inverse_monomial(self, a):
        if len(a) != 1:
            raise ValueError("This catalogue recipe requires a monomial inverse")
        (e, c), = a.items()
        return {(-e) % self.P: (0, 1, 3, 2)[c]}

    def inverse(self, a):
        """Invert a cyclic-ring unit, retaining the power-of-two fast path."""
        # A nonzero monomial is a unit for every cyclic period, including odd P.
        if len(a) == 1:
            return self.inverse_monomial(a)
        augmentation = 0
        for coefficient in a.values():
            augmentation ^= coefficient
        if not augmentation or self.P < 1:
            raise ValueError('Polynomial is not a cyclic-ring unit')
        if self.P & (self.P - 1):
            # Solve a*b=1 in the regular F4 representation. Unlike the local
            # power-of-two case, nonzero augmentation alone does not establish
            # invertibility at odd or general composite periods.
            matrix = [[a.get((i-j) % self.P, 0) for j in range(self.P)]
                      + [int(i == 0)] for i in range(self.P)]
            for j in range(self.P):
                pivot = next((i for i in range(j, self.P) if matrix[i][j]), None)
                if pivot is None:
                    raise ValueError('Polynomial is not a cyclic-ring unit')
                matrix[j], matrix[pivot] = matrix[pivot], matrix[j]
                scalar = (0, 1, 3, 2)[matrix[j][j]]
                matrix[j] = [MUL[scalar][c] for c in matrix[j]]
                for i in range(self.P):
                    if i != j and matrix[i][j]:
                        scalar = matrix[i][j]
                        matrix[i] = [x ^ MUL[scalar][y] for x, y in zip(matrix[i], matrix[j])]
            result = {i: row[-1] for i, row in enumerate(matrix) if row[-1]}
            if self.mul(a, result) != {0: 1}:
                raise ValueError('Polynomial inverse check failed')
            return result
        scalar_inverse = {0: (0, 1, 3, 2)[augmentation]}
        term = self.add(self.mul(a, scalar_inverse), {0: 1})
        value, power = {0: 1}, 1
        while power < self.P:
            value = self.mul(value, self.add({0: 1}, term))
            term = self.mul(term, term)
            power *= 2
        result = self.mul(value, scalar_inverse)
        if self.mul(a, result) != {0: 1}:
            raise ValueError('Polynomial inverse check failed')
        return result

    def adj(self, a):
        n = len(a)
        return [[self.det([[a[u][v] for v in range(n) if v != i]
                          for u in range(n) if u != j])
                 for j in range(n)] for i in range(n)]

    def matvec(self, a, v):
        return [self.add(*(self.mul(x, y) for x, y in zip(row, v))) for row in a]


def support(v):
    out = []
    while v:
        low = v & -v
        out.append(low.bit_length() - 1)
        v ^= low
    return out


def from_support(indices, n):
    if indices != sorted(set(indices)) or any(type(i) is not int or not 0 <= i < n for i in indices):
        raise ValueError("Support must contain sorted, unique, in-range integers")
    return sum(1 << i for i in indices)


def echelon(rows):
    pivots = {}
    for value in rows:
        while value:
            p = value.bit_length() - 1
            if p not in pivots:
                pivots[p] = value
                break
            value ^= pivots[p]
    return pivots


def in_span(value, pivots):
    while value:
        p = value.bit_length() - 1
        if p not in pivots:
            return False
        value ^= pivots[p]
    return True


def translate(v, shift, P):
    return sum(1 << ((i // P) * P + (i + shift) % P) for i in support(v))


def binary_vector(vector, P, field, scalar=1):
    out = 0
    degree = 1 if field == 2 else 2
    for block, poly in enumerate(vector):
        for e, c in poly.items():
            bits = c if field == 2 else COORD[MUL[c][scalar]]
            for b in range(degree):
                if bits >> b & 1:
                    out ^= 1 << ((degree * block + b) * P + e)
    return out


def binary_checks(matrix, P, field):
    degree = 1 if field == 2 else 2
    beta = (1,) if field == 2 else (2, 3)
    out = []
    for row in matrix:
        for a in range(degree):
            seed = 0
            for block, poly in enumerate(row):
                for e, c in poly.items():
                    for b, scalar in enumerate(beta):
                        bits = c if field == 2 else COORD[MUL[c][scalar]]
                        if bits >> a & 1:
                            seed ^= 1 << ((degree * block + b) * P + (-e) % P)
            out.extend(translate(seed, g, P) for g in range(P))
    return out


def balance(rows, P, choices=None):
    out = []
    for h in range(len(rows) // (2 * P)):
        a, b = rows[2*h*P], rows[(2*h+1)*P]
        three = (a, b, a ^ b)
        chosen = choices[h] if choices is not None else sorted(range(3), key=lambda i: (three[i].bit_count(), i))[:2]
        for i in chosen:
            out.extend(translate(three[i], g, P) for g in range(P))
    return out


def reconstruct(code):
    spec = code['construction']
    P, field = code['group']['order'], code['field']
    ring = Ring(P)
    if spec['family'] in ('lp', 'lp_partial'):
        A = [[ring.add(*({e: c} for e, c in entry)) for entry in row]
             for row in spec['constituent']]
        r, s = len(A), len(A[0])
        L = s*s + r*r
        HX, HZ = [], []
        for i in range(r):
            for b in range(s):
                row = [{} for _ in range(L)]
                for a in range(s):
                    row[a*s+b] = A[i][a]
                for c in range(r):
                    row[s*s+i*r+c] = ring.bar(A[c][b])
                HX.append(row)
        for a in range(s):
            for j in range(r):
                row = [{} for _ in range(L)]
                for b in range(s):
                    row[a*s+b] = A[j][b]
                for i in range(r):
                    row[s*s+i*r+j] = ring.bar(A[i][a])
                HZ.append(row)
        pivots = spec['pivot_columns']
        if len(pivots) != r or len(set(pivots)) != r or any(not 0 <= i < s for i in pivots):
            raise ValueError('Invalid pivot columns')
        free = [i for i in range(s) if i not in pivots]
        if spec['family'] == 'lp_partial':
            U = [[ring.add(*({e: c} for e, c in entry)) for entry in column]
                 for column in spec['kernel_columns']]
            if len(U) != s-r or any(len(column) != s for column in U):
                raise ValueError('Invalid partial kernel columns')
            for j, column in enumerate(U):
                if any(ring.matvec(A, column)) or any(column[a] != ({0: 1} if i == j else {}) for i, a in enumerate(free)):
                    raise ValueError('Partial kernel columns violate syndrome/free-coordinate conditions')
        else:
            D = [[row[i] for i in pivots] for row in A]
            inverse = ring.inverse(ring.det(D))
            adj = ring.adj(D)
            U = []
            for j in range(s-r):
                top = ring.matvec(adj, [row[free[j]] for row in A])
                column = [{} for _ in range(s)]
                for i, value in zip(pivots, top):
                    column[i] = ring.mul(value, inverse)
                column[free[j]] = {0: 1}
                U.append(column)
        seeds = []
        for i in range(s-r):
            for j in range(s-r):
                x, z = [{} for _ in range(L)], [{} for _ in range(L)]
                for a in range(s):
                    x[free[i]*s+a] = U[j][a]
                    z[a*s+free[j]] = U[i][a]
                for scalar in ((1,) if field == 2 else (2, 3)):
                    seeds.append((binary_vector(x, P, field, scalar), binary_vector(z, P, field, scalar)))
    else:
        E, C = spec['exponents'], spec['coefficients']
        L = len(E[0])
        HX = [[{} if c == 0 else {e: c} for e, c in zip(er, cr)] for er, cr in zip(E, C)]
        HZ = [[ring.bar(row[(j+L//2) % L]) for j in range(L)] for row in HX]
        if field == 2:
            z = [[{e: 1 for e in support(int(mask, 16))} for mask in word] for word in spec['Z_hex']]
            order = spec.get('partner_seed_order', list(reversed(range(len(z)))))
            x = [[ring.bar(z[i][(j+L//2) % L]) for j in range(L)] for i in order]
            seeds = [(binary_vector(a, P, field), binary_vector(b, P, field)) for a, b in zip(x, z)]
        else:
            # Disjoint correction pivots for a reversing r-by-2M PP array.
            # The 3-by-8 construction is the default special case.
            r, M = len(HX), L // 2
            p = tuple(spec.get('pivot_columns', range(r)))
            if L % 2 or len(p) != r or len(set(p)) != r or any(not 0 <= j < L for j in p):
                raise ValueError('Invalid reversing-PP pivot columns')
            pp = tuple((j + M) % L for j in p)
            if set(p) & set(pp):
                raise ValueError('Reversing-PP correction pivots overlap')
            free = tuple(j for j in range(L) if j not in p and j not in pp)
            SX, SZ = [[row[j] for j in p] for row in HX], [[row[j] for j in pp] for row in HZ]
            dx, dz = ring.det(SX), ring.det(SZ)
            inverse = ring.inverse_monomial(ring.mul(ring.bar(dz), dx))
            seeds = []
            scalars = spec.get('reciprocal_seed_scalars', [1] * len(free))
            if len(scalars) != len(free) or any(c not in (1, 2, 3) for c in scalars):
                raise ValueError('Invalid reciprocal seed scalars')
            for a, mu in zip(free, scalars):
                x, z = [{} for _ in range(L)], [{} for _ in range(L)]
                for j, value in zip(p, ring.matvec(ring.adj(SX), [row[a] for row in HX])):
                    z[j] = value
                for j, value in zip(pp, ring.matvec(ring.adj(SZ), [row[a] for row in HZ])):
                    x[j] = value
                x[a], z[a] = dz, dx
                z = [ring.mul(value, inverse) for value in z]
                x = [ring.mul(value, {0: (0, 1, 3, 2)[mu]}) for value in x]
                z = [ring.mul(value, {0: mu}) for value in z]
                for scalar in (2, 3):
                    seeds.append((binary_vector(x, P, field, scalar), binary_vector(z, P, field, scalar)))
    for operation in spec.get('paired_binary_shears', []):
        # x_i <- x_i + t^h x_j, z_j <- z_j + t^-h z_i preserves
        # every translated canonical pairing, and leaves the checks unchanged.
        i, j, shift = operation['target'], operation['source'], operation.get('shift', 0)
        if i == j or not 0 <= i < len(seeds) or not 0 <= j < len(seeds):
            raise ValueError('Invalid paired binary shear')
        xi, zi = seeds[i]
        xj, zj = seeds[j]
        seeds[i] = (xi ^ translate(xj, shift, P), zi)
        seeds[j] = (xj, zj ^ translate(zi, -shift, P))
    hx, hz = binary_checks(HX, P, field), binary_checks(HZ, P, field)
    if spec.get('check_basis') == 'balanced':
        hx, hz = balance(hx, P), balance(hz, P)
    elif spec.get('check_basis') == 'specified':
        hx, hz = balance(hx, P, spec['row_choices']), balance(hz, P, spec['row_choices'])
    return hx, hz, seeds


def verify(code, checks, frame):
    def require(ok, message):
        if not ok:
            raise ValueError(f"{code['id']}: {message}")
    n, k, P = code['n'], code['k'], code['group']['order']
    hx = [from_support(row, n) for row in checks['X']]
    hz = [from_support(row, n) for row in checks['Z']]
    seeds = [(from_support(s['X'], n), from_support(s['Z'], n)) for s in frame['seeds']]
    require((hx, hz, seeds) == reconstruct(code), 'stored checks/frame differ from construction recipe')
    require(all((x & z).bit_count() % 2 == 0 for x in hx for z in hz), 'CSS commutation')
    xx = [translate(x, g, P) for x, _ in seeds for g in range(P)]
    zz = [translate(z, g, P) for _, z in seeds for g in range(P)]
    require(len(xx) == len(zz) == k, 'seed orbit count')
    require(all(not (x & h).bit_count() % 2 for x in xx for h in hz), 'X logical syndrome')
    require(all(not (z & h).bit_count() % 2 for z in zz for h in hx), 'Z logical syndrome')
    require(all((x & z).bit_count() % 2 == (i == j) for i, x in enumerate(xx) for j, z in enumerate(zz)), 'complete translated canonical pairing')
    bx, bz = echelon(hx), echelon(hz)
    require(n - len(bx) - len(bz) == k, 'code dimension')
    require(len(echelon(hx+xx)) == len(bx)+k and len(echelon(hz+zz)) == len(bz)+k, 'logical completeness')
    width = max(v.bit_count() for pair in seeds for v in pair)
    require(width == frame['width'] == code['frame_width']['upper'], 'frame width')
    witnesses = []
    for record in code['witnesses']:
        v = from_support(record['support'], n)
        h, b = (hx, bz) if record['sector'] == 'Z' else (hz, bx)
        require(not any((v & row).bit_count() % 2 for row in h) and not in_span(v, b), 'distance witness')
        require(v.bit_count() == record['weight'] == code['distance']['upper'], 'witness weight')
        witnesses.append(record['weight'])
    require(bool(witnesses), 'missing physical upper-bound witness')
    if code['construction']['family'] == 'pp':
        blocks = n // P
        perm = lambda v: sum(1 << (((i//P+blocks//2) % blocks)*P+(-i)%P) for i in support(v))
    else:
        r = len(code['construction']['constituent'])
        s = len(code['construction']['constituent'][0])
        degree = 1 if code['field'] == 2 else 2
        def perm(v):
            result = 0
            for i in support(v):
                block, rem = divmod(i, degree*P)
                if block < s*s:
                    a, b = divmod(block, s)
                    dst = b*s+a
                else:
                    a, b = divmod(block-s*s, r)
                    dst = s*s+b*r+a
                result |= 1 << (dst*degree*P+rem)
            return result
    require(all(in_span(perm(h), bz) for h in hx) and all(in_span(perm(h), bx) for h in hz), 'sector isometry')
    stats = {}
    for sector, rows in [('X', hx), ('Z', hz)]:
        columns = [0]*n
        for row in rows:
            for i in support(row):
                columns[i] += 1
        stats[sector] = {'rank': len(echelon(rows)), 'maximum_check_weight': max(map(int.bit_count, rows)), 'maximum_column_weight': max(columns)}
    if 'maximum_check_weight' in code:
        require(all(v['maximum_check_weight'] == code['maximum_check_weight'] and v['maximum_column_weight'] == code['maximum_column_weight'] for v in stats.values()), 'check/column statistics')
    return {'id': code['id'], 'status': 'PASS', 'n': n, 'k': k, 'frame_width': width, 'translated_pairings_checked': k*k, 'checks': stats, 'witness_weights': witnesses, 'distance_exhaustive_search_rerun': False}
