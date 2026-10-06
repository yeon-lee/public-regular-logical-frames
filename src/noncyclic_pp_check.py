#!/usr/bin/env python3
"""Independent, standard-library audit of one noncyclic PP code and its receipts.

This verifies physical algebra and completed-execution metadata. It does not
replay the exhaustive search or turn its logs into a formal proof object.
"""

if not __debug__:
    raise RuntimeError(
        'Verification requires Python assertions; run without -O/-OO and unset PYTHONOPTIMIZE.'
    )

import argparse
import hashlib
import itertools
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
IDENTIFIER = "pp-f4-c16xc2-512-128-19-w10"



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


def support(v):
    out = []
    while v:
        q = v & -v
        out.append(q.bit_length() - 1)
        v ^= q
    return out


def bits(items):
    assert len(items) == len(set(items)) and all(0 <= j < 512 for j in items)
    return sum(1 << j for j in items)


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


def remainder(v, base):
    while v:
        p = v.bit_length() - 1
        if p not in base:
            return v
        v ^= base[p]
    return 0


def dot(a, b):
    return (a & b).bit_count() % 2


def matrix(path):
    lines = path.read_text().splitlines()
    n, m, g = map(int, lines[0].split())
    assert n == 512 and len(lines) == 1 + m + g
    rows = []
    for line in lines[1:]:
        v = list(map(int, line.split()))
        assert v[0] == len(v) - 1
        rows.append(bits(v[1:]))
    return rows[:m], rows[m:]


def gf(a, b):
    """F4 = F2[omega]/(omega^2+omega+1); 2=omega, 3=omega^2."""
    v = 0
    for _ in range(2):
        if b & 1:
            v ^= a
        b >>= 1
        a <<= 1
        if a & 4:
            a ^= 7
    return v


def coords(c):
    """Coordinates in the trace-self-dual basis (omega, omega^2)."""
    return (((c >> 1) ^ c) & 1, c & 1)


def dec(e):
    return e % 16, e // 16


def enc(g):
    return g[0] + 16 * g[1]


def addg(a, b):
    return (a[0] + b[0]) % 16, (a[1] + b[1]) % 2


def neg(g):
    return (-g[0]) % 16, g[1]


GROUP = [(a, b) for b in range(2) for a in range(16)]


def add(*polys):
    result = {}
    for p in polys:
        for g, c in p.items():
            result[g] = result.get(g, 0) ^ c
    return {g: c for g, c in result.items() if c}


def mul(p, q):
    result = {}
    for g, a in p.items():
        for h, b in q.items():
            k = addg(g, h)
            result[k] = result.get(k, 0) ^ gf(a, b)
    return {g: c for g, c in result.items() if c}


def perm(v, translate=(0, 0), reflect=False):
    out = 0
    for j in support(v):
        block, e = divmod(j, 32)
        g = dec(e)
        if reflect:
            block, g = (block + 8) % 16, neg(g)
        out |= 1 << (block * 32 + enc(addg(g, translate)))
    return out


def check():
    data, completion = HERE / "data", HERE / "completion"
    read = lambda p: json.loads(p.read_text())
    raw = read(data / "code.json")
    cert = read(data / "certificate.json")
    receipt = read(completion / "receipt.json")
    cov = read(completion / "coverage.json")
    provenance = read(HERE / "provenance.json")
    assert raw["group"] == "C16xC2" and raw["P"] == 32
    assert raw["cyclic_subgroup_period"] == 16
    C, E = raw["C"], raw["E"]
    assert len(C) == len(E) == 3
    assert all(len(row) == 8 for row in C + E)
    assert all(c in (1, 2, 3) for row in C for c in row)
    assert all(0 <= e < 32 for row in E for e in row)

    # Reconstruct the physical matrices directly from tuple group arithmetic.
    def expand(reverse):
        result = []
        for i in range(3):
            for a in range(2):
                for h in GROUP:
                    row = 0
                    for ell in range(8):
                        j = (ell + 4) % 8 if reverse else ell
                        e = neg(dec(E[i][j])) if reverse else dec(E[i][j])
                        for b, beta in enumerate((2, 3)):
                            if coords(gf(C[i][j], beta))[a]:
                                row ^= 1 << ((2 * ell + b) * 32 + enc(addg(h, neg(e))))
                    result.append(row)
        return result

    HX, HZ = expand(False), expand(True)
    assert (HX, HZ) == matrix(data / "code.mat")
    BX, BZ = basis(HX), basis(HZ)
    assert len(BX) == len(BZ) == 192
    assert all(dot(x, z) == 0 for x in HX for z in HZ)
    assert set(perm(h, reflect=True) for h in HX) == set(HZ)
    for rows, base in [(HX, BX), (HZ, BZ)]:
        for t in [(1, 0), (0, 1)]:
            assert all(remainder(perm(h, t), base) == 0 for h in rows)
        assert set(h.bit_count() for h in rows) == {10}
        assert max(sum((h >> j) & 1 for h in rows) for j in range(512)) == 4

    A = [[{dec(E[i][j]): C[i][j]} for j in range(8)] for i in range(3)]

    def minor(cols):
        result = {}
        for pp in itertools.permutations(cols):
            term = {(0, 0): 1}
            for i, j in enumerate(pp):
                term = mul(term, A[i][j])
            result = add(result, term)
        return result

    assert raw["pivot_columns"] == [0, 1, 2] and raw["free_columns"] == [3, 7]
    delta = minor((0, 1, 2))
    pairing = mul(delta, delta)
    assert len(pairing) == 1
    gp, cp = next(iter(pairing.items()))
    inverse = {neg(gp): gf(cp, cp)}
    assert mul(pairing, inverse) == {(0, 0): 1}
    for p, key in [(delta, "delta"), (pairing, "pairing"), (inverse, "pairing_inverse")]:
        assert {str(enc(g)): c for g, c in p.items()} == raw[key]
    polys = []
    for free in (3, 7):
        p = [{} for _ in range(8)]
        p[free] = delta
        for j in range(3):
            cols = [0, 1, 2]
            cols[j] = free
            p[j] = minor(cols)
        # Cofactor syzygy, checked in the entire product-group ring.
        assert all(not add(*(mul(A[i][j], p[j]) for j in range(8))) for i in range(3))
        polys.append(p)

    def binary(p, beta):
        row = 0
        for ell, q in enumerate(p):
            for g, c in q.items():
                for b, val in enumerate(coords(gf(beta, c))):
                    if val:
                        row ^= 1 << ((2 * ell + b) * 32 + enc(g))
        return row

    X, Z = [], []
    for i in range(2):
        for beta in (2, 3):
            X.append(perm(binary(polys[1 - i], beta), reflect=True))
            Z.append(binary([mul(p, inverse) for p in polys[i]], beta))
    assert list(map(support, X)) == raw["X_seed_supports"]
    assert list(map(support, Z)) == raw["Z_seed_supports"]
    XX = [perm(v, t) for v in X for t in GROUP]
    ZZ = [perm(v, t) for v in Z for t in GROUP]
    assert len(XX) == len(set(XX)) == len(ZZ) == len(set(ZZ)) == 128
    assert all(dot(x, z) == int(i == j) for i, x in enumerate(XX) for j, z in enumerate(ZZ))
    assert all(dot(v, h) == 0 for v in XX for h in HZ)
    assert all(dot(v, h) == 0 for v in ZZ for h in HX)
    assert len(basis(HX + XX)) == len(basis(HZ + ZZ)) == 320
    assert max(v.bit_count() for v in XX + ZZ) == 25
    witness = bits(raw["witness_Z"])
    assert witness.bit_count() == 19 and all(dot(witness, h) == 0 for h in HX)
    assert remainder(witness, BZ) != 0
    assert cert["best_witness"]["support"] == support(witness)
    assert all(dot(perm(witness, reflect=True), h) == 0 for h in HZ)
    assert remainder(perm(witness, reflect=True), BX) != 0

    # Bind the redundant search matrices to the reconstructed physical code.
    RH, RZ = matrix(completion / "input.mat")
    for original, red in [(HX, RH), (HZ, RZ)]:
        base = basis(original)
        assert len(basis(red)) == len(base)
        assert all(remainder(v, base) == 0 for v in red)
    assert max(v.bit_count() for v in RH) <= 64
    expected_pi = [j ^ 16 for j in range(512)]
    assert cov["C2_coordinate_permutation"] == expected_pi
    assert all(perm(1 << j, (0, 1)) == 1 << expected_pi[j] for j in range(512))
    # Every position can be moved to zero inside its original 32-site block.
    for j in range(512):
        block, e = divmod(j, 32)
        assert perm(1 << j, neg(dec(e))) == 1 << (32 * block)

    assert receipt["mode"] == "exact_verified_product_group_representatives"
    assert receipt["sector"] == "Z" and receipt["lower_bound_this_sector"] == 19
    for obj in (receipt, cov):
        assert obj["original_matrix_sha256"] == sha(data / "code.mat")
        assert obj["source_sha256"] == sha(HERE / "src/exact.cpp")
        assert obj["engine_sha256"] == provenance["historical_engine_binary_sha256"]
    assert cov["input_sha256"] == receipt["execution_input_sha256"] == sha(completion / "input.mat")
    assert receipt["coverage_sha256"] == sha(completion / "coverage.json")
    roots = list(range(0, 32, 2))
    assert cov["required_roots"] == cov["covered_roots"] == roots and cov["complete"]
    assert cov["native_period"] == receipt["final"]["native_period"] == 16
    assert receipt["final"]["required_group_representatives"] == roots
    assert receipt["final"]["completed_group_representatives"] == roots
    nodes, seconds = 0, 0.0
    assert len(cov["jobs"]) == 16
    for root, job in zip(roots, cov["jobs"]):
        stem = f"root{root:02d}"
        log = completion / (stem + ".stdout.jsonl")
        rec = read(completion / (stem + ".json"))
        events = [json.loads(line) for line in log.read_text().splitlines()]
        assert rec == job and len(events) == 1 and rec["final"] == events[0]
        result = events[0]
        assert rec["returncode"] == 0 and rec["stdout_sha256"] == sha(log)
        assert (completion / (stem + ".stderr")).read_text() == ""
        assert result["status"] == "excluded" and result["W"] == 18
        assert result["root"] == root and result["roots_completed"] == 1
        assert result["partition"] == "root_branch_v1" and result["coset_pruning"]
        assert result["branch_begin"] == 0 and result["branch_end"] == 64
        assert result["shards"] == 1 and result["shard"] == 0
        cmd = rec["command"]
        # Normalized paths are metadata only; resolve all inputs locally.
        assert Path(cmd[0]).name == "exact" and Path(cmd[1]).name == "input.mat"
        assert list(map(int, cmd[2:4])) == [16, 18]
        assert list(map(int, cmd[5:8])) == [root, 0, 64]
        assert int(cmd[8]) == 2**64 - 1 and len(cmd) == 9
        nodes += result["nodes"]
        seconds += result["seconds"]

    assert cert["code_json_sha256"] == sha(data / "code.json")
    assert cert["audit_sha256"] == sha(data / "audit.json")
    assert cert["matrix_sha256"] == sha(data / "code.mat")
    assert raw["source_sha256"] == sha(data / "source.json")
    assert len(cert["receipts"]) == 1
    assert cert["receipts"][0]["sha256"] == sha(completion / "receipt.json")
    assert cert["distance_lower_bound"] == cert["distance_upper_bound"] == cert["distance_exact"] == 19
    assert cert["displayed_frame_width"] == cert["canonical_width_upper_bound"] == 25
    assert not cert["optimal_width_claimed"]
    for name, digest in provenance["source_files"].items():
        assert sha(HERE / "src" / name) == digest
    if (HERE / "SHA256SUMS").exists():
        for line in (HERE / "SHA256SUMS").read_text().splitlines():
            digest, name = line.split("  ", 1)
            assert sha(HERE / name) == digest, name
    return {
        "status": "PASS", "id": IDENTIFIER, "group": "C16xC2", "group_order": 32,
        "not_a_folded_code": True, "n": 512, "k": 128, "ranks": [192, 192],
        "d_X": 19, "d_Z": 19, "distance_status": "witness plus completed exhaustive exclusion through 18",
        "check_weight_each_sector": 10, "maximum_column_weight_each_sector": 4,
        "displayed_frame_width": 25, "optimized_frame_width_interval": [19, 25],
        "binary_seed_pairs": 4, "translates_per_seed": 32, "canonical_pairs": 128,
        "X_seed_weights": list(map(int.bit_count, X)), "Z_seed_weights": list(map(int.bit_count, Z)),
        "all_16384_cross_pairings_checked": True, "augmented_ranks": [320, 320],
        "native_search_period": 16, "completed_even_roots": roots, "node_visits": nodes,
        "sum_recorded_search_seconds": seconds,
        "scope": "Independent physical/algebra/receipt audit; no exhaustive replay and no formal node-by-node certificate.",
        "proof_of_root_coverage": "Use the first occupied original 32-site block; C2 swaps its halves, then C16 anchors an occupied site at zero. Earlier original blocks remain empty. Native root 2*b fixes exactly that 32*b coordinate. All 16 even roots and all first-branch choices [0,64) are complete.",
    }


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--report", type=Path)
    args = ap.parse_args()
    result = check()
    output = json.dumps(result, indent=2) + "\n"
    if args.report:
        args.report.write_text(output)
    print(output, end="")
