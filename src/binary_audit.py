"""Independent integer-bitset CSS, cyclic-frame, and witness verification."""

if not __debug__:
    raise RuntimeError(
        'Verification requires Python assertions; run without -O/-OO and unset PYTHONOPTIMIZE.'
    )

import hashlib
import json
from collections import Counter
from pathlib import Path



# Public releases omit binaries; their identities are declarations only.
_public_root = next((q for q in Path(__file__).resolve().parents
                     if (q / 'public_integrity.py').is_file()), None)
if _public_root is not None:
    import sys as _public_sys
    _public_sys.path.insert(0, str(_public_root))
    from public_integrity import evidence_digest
else:
    def evidence_digest(path):
        return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def sha(path):
    """File digest or explicitly declared omitted executable identity."""
    return evidence_digest(path)


def basis(rows):
    out = {}
    for v in rows:
        while v:
            p = v.bit_length() - 1
            if p in out:
                v ^= out[p]
            else:
                out[p] = v
                break
    return out


def remainder(v, b):
    while v:
        p = v.bit_length() - 1
        if p not in b:
            return v
        v ^= b[p]
    return 0


def bits(support, n):
    assert len(support) == len(set(support))
    assert all(isinstance(j, int) and 0 <= j < n for j in support)
    return sum(1 << j for j in support)


def support(v):
    out = []
    while v:
        low = v & -v
        out.append(low.bit_length() - 1)
        v ^= low
    return out


def read_matrix(path):
    lines = Path(path).read_text().splitlines()
    n, m, g = map(int, lines[0].split())
    assert n > 0 and m >= 0 and g >= 0
    rows = []
    for line in lines[1:m + g + 1]:
        a = list(map(int, line.split()))
        assert a[0] == len(a) - 1
        rows.append(bits(a[1:], n))
    assert len(rows) == m + g
    # Samplers may have a trailing witness list; the CSS itself ends above.
    return n, rows[:m], rows[m:]


def write_matrix(path, n, H, G, witnesses=None):
    lines = [f'{n} {len(H)} {len(G)}']
    for v in H + G:
        a = support(v)
        lines.append(' '.join(map(str, [len(a)] + a)))
    if witnesses is not None:
        lines.append(str(len(witnesses)))
        for v in witnesses:
            a = support(v)
            lines.append(' '.join(map(str, [len(a)] + a)))
    Path(path).write_text('\n'.join(lines) + '\n')


def dot(a, b):
    return (a & b).bit_count() & 1


def rotate(v, n, P, t=1):
    assert n % P == 0
    t %= P
    mask = (1 << P) - 1
    out = 0
    for b in range(n // P):
        a = (v >> (b * P)) & mask
        a = ((a << t) | (a >> (P - t))) & mask
        out |= a << (b * P)
    return out


def check_witness(v, n, H, G, claimed_weight=None):
    assert 0 < v < 1 << n
    assert all(not dot(h, v) for h in H), 'nonzero syndrome'
    assert remainder(v, basis(G)), 'stabilizer witness'
    w = v.bit_count()
    if claimed_weight is not None:
        assert w == claimed_weight
    return {'weight': w, 'support': support(v), 'syndrome_zero': True,
            'outside_stabilizer_span': True}


def audit_css(n, HX, HZ, P, Xseeds=None, Zseeds=None):
    assert n % P == 0
    assert all(0 <= v < 1 << n for v in HX + HZ)
    bx, bz = basis(HX), basis(HZ)
    assert all(not dot(x, z) for x in HX for z in HZ), 'noncommuting CSS'
    for rows, b in ((HX, bx), (HZ, bz)):
        assert all(not remainder(rotate(v, n, P), b) for v in rows)
    k = n - len(bx) - len(bz)
    out = dict(n=n, k=k, P=P, ranks=[len(bx), len(bz)],
               commuting=True, translation_invariant=True,
               check_weight_X=dict(sorted(Counter(v.bit_count() for v in HX).items())),
               check_weight_Z=dict(sorted(Counter(v.bit_count() for v in HZ).items())))
    if Xseeds is not None and Zseeds is not None:
        assert len(Xseeds) == len(Zseeds) == k // P and k % P == 0
        assert all(not dot(x, z) for x in Xseeds for z in HZ)
        assert all(not dot(x, z) for x in HX for z in Zseeds)
        # Translation invariance turns these k^2/P correlations into the full
        # k-by-k translated cross-Gram test without assuming any sector map.
        for i, x in enumerate(Xseeds):
            for j, z in enumerate(Zseeds):
                for t in range(P):
                    assert dot(x, rotate(z, n, P, t)) == int(i == j and t == 0), (i, j, t)
        out.update(full_canonical_frame=True,
                   seed_weights_X=[x.bit_count() for x in Xseeds],
                   seed_weights_Z=[z.bit_count() for z in Zseeds],
                   displayed_width=max(v.bit_count() for v in Xseeds + Zseeds))
    return out


if __name__ == '__main__':
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument('matrix', type=Path)
    p.add_argument('--period', type=int, required=True)
    a = p.parse_args()
    n, H, G = read_matrix(a.matrix)
    print(json.dumps(audit_css(n, H, G, a.period), indent=2))
