"""Independent binary/GF(4) reconstruction of the supplied PP weight-nine code.

The arithmetic helpers are frozen from an earlier independent audit.  The
supplied PP builder is used only for a separate cross-check after reconstruction.
This script verifies algebra and upper-bound witnesses; run_exclusion.py records
the separate exhaustive lower-bound computation.
"""

if not __debug__:
    raise RuntimeError(
        'Verification requires Python assertions; run without -O/-OO and unset PYTHONOPTIMIZE.'
    )

from pathlib import Path
import hashlib
import json
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE / 'source'))
import independent_baseline_audit as a

C = [[3, 1, 2, 1, 1, 1, 1, 1],
     [1, 3, 1, 1, 1, 1, 1, 2],
     [1, 1, 1, 1, 2, 1, 3, 1]]
E = [[60, 25, 58, 14, 44, 42, 46, 31],
     [32, 57, 62, 13, 39, 51, 22, 2],
     [43, 8, 37, 56, 26, 61, 56, 5]]
WITNESS = [0, 41, 64, 68, 140, 163, 227, 231, 258,
           435, 629, 633, 658, 662, 722, 726, 883, 930]


def read_exclusion_input(path):
    lines = iter(path.read_text().splitlines())
    dimensions = tuple(map(int, next(lines).split()))
    matrices = []
    for _ in range(3):
        size = int(next(lines))
        rows = [a.from_support(list(map(int, next(lines).split())))
                for _ in range(size)]
        matrices.append(rows)
    assert not list(lines)
    return dimensions, matrices


def main():
    P, n = 64, 1024
    HX, HZ = a.expand(C, E), a.expand(C, E, True)
    BX, BZ = a.basis(HX), a.basis(HZ)
    assert len(BX) == len(BZ) == 384
    assert not any((x & z).bit_count() % 2 for x in HX for z in HZ)
    assert {a.reflection(x) for x in HX} == set(HZ)
    assert {a.reflection(z) for z in HZ} == set(HX)
    assert {a.translate(x, 1) for x in HX} == set(HX)
    assert {a.translate(z, 1) for z in HZ} == set(HZ)

    delta = a.determinant(C, E, [0, 1, 2])
    assert delta == {2: 2, 26: 2, 34: 2}
    g = a.pmul(delta, delta)
    assert g == {52: 3}
    ginv = {12: 2}
    assert a.pmul(g, ginv) == {0: 1}
    rawz = []
    for free in [3, 7]:
        polynomials = [{} for _ in range(8)]
        polynomials[free] = delta
        for j in range(3):
            columns = [0, 1, 2]
            columns[j] = free
            polynomials[j] = a.determinant(C, E, columns)
        rawz.append(polynomials)
    X, Z = [], []
    for i in range(2):
        other = rawz[1 - i]
        X.extend(a.binary_seed([
            {(-e) % P: c for e, c in other[(j + 4) % 8].items()}
            for j in range(8)]))
        Z.extend(a.binary_seed([a.pmul(poly, ginv) for poly in rawz[i]]))
    assert [v.bit_count() for v in X] == [26, 26, 22, 24]
    assert [v.bit_count() for v in Z] == [24, 26, 26, 28]
    TX = [a.translate(x, t) for x in X for t in range(P)]
    TZ = [a.translate(z, t) for z in Z for t in range(P)]
    assert not any((v & h).bit_count() % 2 for v in TX for h in HZ)
    assert not any((v & h).bit_count() % 2 for v in TZ for h in HX)
    assert all((x & z).bit_count() % 2 == (i == j)
               for i, x in enumerate(TX) for j, z in enumerate(TZ))
    assert len(a.basis(HX + TX)) == len(a.basis(HZ + TZ)) == 640

    source = HERE / 'input/final_PPw9.txt'
    dims, supplied = read_exclusion_input(source)
    assert dims == (1024, 64, 16)
    assert supplied == [HX, HZ, TX]
    assert hashlib.md5(source.read_bytes()).hexdigest() == 'b5459cf86ad06d5a9772810a6ac4a2e6'
    row_hist = {s: a.hist(h.bit_count() for h in H)
                for s, H in [('X', HX), ('Z', HZ)]}
    col_hist = {s: a.hist(sum(h >> j & 1 for h in H) for j in range(n))
                for s, H in [('X', HX), ('Z', HZ)]}
    assert row_hist == {'X': {'9': 384}, 'Z': {'9': 384}}
    assert col_hist == {'X': {'3': 640, '4': 384},
                        'Z': {'3': 640, '4': 384}}
    wz = a.from_support(WITNESS)
    wx = a.reflection(wz)
    for v, H, B in [(wz, HX, BZ), (wx, HZ, BX)]:
        assert v.bit_count() == 18
        assert not any((v & h).bit_count() % 2 for h in H)
        assert a.remainder(v, B)
    assert any((wz & x).bit_count() % 2 for x in TX)
    assert any((wx & z).bit_count() % 2 for z in TZ)

    # This check uses the vendor builder only after the independent assertions.
    from pp import PP
    code = PP(C, E, P)
    vendor_hx, vendor_hz = code.check_rows()
    vendor_z, vendor_x = code.seeds()
    assert list(map(a.support, HX)) == vendor_hx
    assert list(map(a.support, HZ)) == vendor_hz
    assert list(map(a.support, X)) == vendor_x
    assert list(map(a.support, Z)) == vendor_z

    report = dict(
        status='PASS', n=n, k=n-len(BX)-len(BZ), verified_distance_upper_bound=18,
        distance_lower_bound='Separate exhaustive exclusion required',
        P=P, C=C, E=E, field_encoding='0,1,2,3 = 0,1,omega,omega^2',
        coordinate='(2*quaternary_column+trace_bit)*64+cyclic_position',
        ranks={'HX':384, 'HZ':384, 'HX_plus_X_frame':640, 'HZ_plus_Z_frame':640},
        CSS_commutation=True, translation_invariance=True, reflection_sector_exchange=True,
        reflection_coordinate='(binary_block,j) -> ((binary_block+8)%16,-j mod 64)',
        row_weight_histograms=row_hist, column_weight_histograms=col_hist,
        cofactor_frame=dict(delta=delta, pairing=g, normalization=ginv,
            X_weights=[v.bit_count() for v in X], Z_weights=[v.bit_count() for v in Z],
            X_seed_supports=list(map(a.support,X)), Z_seed_supports=list(map(a.support,Z)),
            width=28, unrestricted_optimality_established=False,
            regular_optimality_established=False, translated_Gram='I_256',
            pairings_checked=65536, seed_order='free column 3,7; trace basis omega,omega^2'),
        witnesses={'Z':a.support(wz), 'X':a.support(wx), 'weight':18,
                   'zero_syndrome_and_stabilizer_nonmembership_independently_checked':True},
        supplied_builder_matches_independent_reconstruction=True,
        exhaustive_input_exactly_matches_checks_and_complete_X_detectors=True,
        input_sha256=a.sha(source), input_md5=hashlib.md5(source.read_bytes()).hexdigest(),
        source_sha256={str(p.relative_to(HERE)):a.sha(p) for p in
            [Path(__file__), HERE/'source/independent_baseline_audit.py',
             HERE/'source/pp.py', HERE/'source/lpq.py']})
    (HERE/'results/algebraic_validation.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:report[k] for k in ['status','n','k','verified_distance_upper_bound',
                                          'row_weight_histograms','column_weight_histograms']}))
    print('Full translated Gram: I_256; cofactor frame width 28; input matches exactly.')


if __name__ == '__main__':
    main()
