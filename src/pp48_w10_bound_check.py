#!/usr/bin/env python3
"""Audit the P48, check-weight-ten code and its completed distance exclusion."""

if not __debug__:
    raise RuntimeError(
        'Verification requires Python assertions; run without -O/-OO and unset PYTHONOPTIMIZE.'
    )

from pathlib import Path
import hashlib
import json
import sys

PACKAGE = Path(__file__).resolve().parent
REPO = PACKAGE.parents[1]
PROVENANCE = PACKAGE / 'completed'
sys.path.insert(0, str(REPO))
import rlf
sys.path.insert(0, str(PROVENANCE / 'src'))
from binary_audit import read_matrix


def read(path):
    return json.loads(path.read_text())



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


def completed_weight21():
    root = PACKAGE / 'completed'
    reused = set(read(root / 'reused_exclusions.json')['shards'])
    for origin in (root, root / 'reused'):
        for rel, digest in read(origin / 'source_manifest.json').items():
            assert sha(origin / rel) == digest, rel
        preflight = read(origin / 'runs/preflight/receipt.json')
        assert preflight['status'] == 'PASS'
        assert preflight['reflection_isometry_verified']
        assert preflight['matrix_sha256'] == sha(root / 'inputs/exact_Z.mat')
        assert preflight['manifest_sha256'] == sha(origin / 'source_manifest.json')
        assert preflight['binary_sha256']['exact'] == sha(origin / 'bin/exact')
        for filename in ('exact.cpp', 'common.hpp'):
            assert sha(origin / 'src' / filename) == sha(PROVENANCE / 'src' / filename)
    assert sha(root / 'inputs/original.mat') == sha(PROVENANCE / 'inputs/original.mat')
    assert sha(root / 'inputs/exact_Z.mat') == sha(PROVENANCE / 'inputs/exact_Z.mat')
    results = []
    for shard in range(64):
        origin = root / 'reused' if shard in reused else root
        prefix = origin / f'runs/exact/shard_{shard:03d}'
        receipt = read(prefix.with_suffix('.receipt.json'))
        result = receipt['result']
        assert receipt['returncode'] == 0 and receipt['status'] == 'excluded'
        assert not receipt['interrupted_by_controller']
        assert result['status'] == 'excluded' and result['W'] == 21
        assert result['roots_completed'] == 16
        assert result['partition'] == 'hash_prefix_v1'
        assert (result['shard'], result['shards'], result['split_weight']) == (shard, 64, 5)
        assert receipt['command'][2:4] == ['48', '21']
        assert receipt['command'][5:9] == ['--shard', str(shard), '64', '5']
        bindings = {
            'input_sha256': root / 'inputs/exact_Z.mat',
            'manifest_sha256': origin / 'source_manifest.json',
            'executable_sha256': origin / 'bin/exact',
            'stdout_sha256': prefix.with_suffix('.jsonl'),
            'stderr_sha256': prefix.with_suffix('.stderr'),
        }
        for key, path in bindings.items():
            assert receipt[key] == sha(path), (shard, key)
        assert json.loads(prefix.with_suffix('.jsonl').read_text().splitlines()[-1]) == result
        results.append(result)
    assert len({r['frontier_digest'] for r in results}) == 1
    assert len({r['frontier_seen'] for r in results}) == 1
    assert sum(r['frontier_owned'] for r in results) == results[0]['frontier_seen']
    return sum(r['nodes'] for r in results)


def main():
    for rel, value in read(PACKAGE / 'SHA256SUMS.json')['files'].items():
        assert sha(PACKAGE / rel) == value, rel
    for rel, value in read(PROVENANCE / 'source_manifest.json').items():
        assert sha(PROVENANCE / rel) == value, rel
    location = REPO / 'codes/pp-quaternary/pp-f4-768-192-w10'
    code = read(location / 'code.json')
    physical = rlf.verify(code, read(location / 'checks.json'),
                          read(location / 'frames/paper.json'))
    n, hx, hz = read_matrix(PROVENANCE / 'inputs/original.mat')
    rx, rz, _ = rlf.reconstruct(code)
    assert (n, hx, hz) == (768, rx, rz)
    nn, search_h, search_g = read_matrix(PROVENANCE / 'inputs/exact_Z.mat')
    assert nn == n and search_h == hx
    bz, bg = rlf.echelon(hz), rlf.echelon(search_g)
    assert len(bz) == len(bg)
    assert all(rlf.in_span(v, bz) for v in search_g)
    assert all(rlf.in_span(v, bg) for v in hz)
    preflight = read(PROVENANCE / 'runs/preflight/receipt.json')
    assert preflight['status'] == 'PASS' and preflight['reflection_isometry_verified']
    assert preflight['binary_sha256']['exact'] == sha(PROVENANCE / 'bin/exact')
    bounds = read(PACKAGE / 'bounds.json')
    nodes = completed_weight21()
    assert (bounds['distance_lower'], bounds['distance_upper']) == (22, 22)
    assert bounds['distance_exact'] and not bounds['frame_width_optimal']
    assert (code['distance']['lower'], code['distance']['upper']) == (22, 22)
    assert code['distance']['status'] == 'exact'
    assert physical['frame_width'] == 29 and physical['witness_weights'] == [22]
    print(json.dumps(dict(status='PASS', code_id=code['id'],
                          distance_interval=[22, 22], distance_exact=True,
                          reported_frame_width=29,
                          translated_pairings_checked=192**2,
                          completed_exclusion_weight=21,
                          completed_shards=64,
                          total_node_visits_including_shared_prefixes=nodes,
                          full_search_replayed=False), indent=2))


if __name__ == '__main__':
    main()
