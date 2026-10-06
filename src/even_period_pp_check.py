"""Independently reconstruct the code and audit complete distance-search coverage.

The quick audit checks completed execution records, not every search node.
--replay compiles the frozen source and runs an unsplit full-code exclusion.
"""

if not __debug__:
    raise RuntimeError(
        'Verification requires Python assertions; run without -O/-OO and unset PYTHONOPTIMIZE.'
    )

from pathlib import Path
import argparse
import hashlib
import json
import subprocess
import sys
import tempfile

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE / 'vendor'))
from audit_pp import reconstruct
from binary_audit import read_matrix, bits, check_witness, basis, remainder, rotate, support



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

def sha(p):
    """File digest or explicitly declared omitted executable identity."""
    return evidence_digest(p)


def read(p):
    return json.loads(p.read_text())


def validate_record(sid, receipt, events, expected):
    f = receipt['final']
    assert len(events) == 1 and events[0] == f
    assert receipt['returncode'] == 0 and not receipt['witnesses']
    assert f['event'] == 'result' and f['status'] == 'excluded'
    assert f['W'] == 19 and f['coset_pruning'] is True
    assert f['partition'] == 'hash_prefix_v1'
    assert (f['shard'], f['shards'], f['split_weight']) == (sid, 64, 4)
    assert (f['root'], f['branch_begin'], f['branch_end']) == (0, 0, 64)
    assert f['roots_completed'] == 16
    assert f['nodes'] >= f['common_nodes'] > 0
    assert 0 <= f['frontier_owned'] <= f['frontier_seen']
    for key, value in expected.items():
        assert receipt[key] == value, key
    c = receipt['command']
    assert len(c) == 10 and Path(c[0]).name == 'exact'
    assert Path(c[1]).name == 'execution.mat'
    assert c[2:4] == ['48', '19'] and float(c[4]) > f['seconds']
    assert c[5:] == ['--shard', str(sid), '64', '4', '18446744073709551615']
    return f


def check():
    meta = read(HERE / 'INPUT.json')
    assert meta['code_tag'] == 'b61735013df0f220'
    assert (meta['n'], meta['k'], meta['P'], meta['W']) == (768, 192, 48, 19)
    for name, value in meta['inputs'].items():
        assert sha(HERE / name) == value, name
    expected = dict(
        execution_input_sha256=sha(HERE / 'data/execution.mat'),
        engine_sha256=sha(HERE / 'vendor/exact'),
        source_sha256=sha(HERE / 'vendor/exact.cpp'),
        header_sha256=sha(HERE / 'vendor/common.hpp'))
    # Rebuild checks and normalized cofactor seeds from C,E, using independent
    # F4 polynomial arithmetic. No local-ring test is assumed at period 48.
    with tempfile.TemporaryDirectory(prefix='pp-p48-audit-') as tmp:
        out = Path(tmp) / 'reconstructed'
        rec, algebra = reconstruct(HERE / 'data/code.json', out)
        assert (out / 'code.mat').read_bytes() == (HERE / 'data/code.mat').read_bytes()
        assert (out / 'redundant3.mat').read_bytes() == (HERE / 'data/execution.mat').read_bytes()
    assert algebra['ranks'] == [288, 288] and algebra['k'] == 192
    assert algebra['physical_components'] == [768]
    assert algebra['displayed_width'] == 30
    assert algebra['X_check_max'] == algebra['Z_check_max'] == 9
    assert algebra['X_column_max'] == algebra['Z_column_max'] == 4
    assert algebra['full_translated_pairing_entries'] == 192**2
    n, HX, HZ = read_matrix(HERE / 'data/code.mat')
    nn, RX, RZ = read_matrix(HERE / 'data/execution.mat')
    assert nn == n == 768
    for original, redundant in [(HX, RX), (HZ, RZ)]:
        b = basis(original)
        assert len(b) == len(basis(redundant)) == 288
        assert all(not remainder(v, b) for v in redundant)

    def reflect(v):
        return sum(1 << (((j // 48 + 8) % 16) * 48 + (-j) % 48) for j in support(v))

    assert {reflect(v) for v in HX} == set(HZ)
    assert {reflect(v) for v in HZ} == set(HX)
    pilot = read(HERE / 'provenance/pilot_results.json')
    w = bits(pilot['witness']['support'], n)
    upper = check_witness(w, n, HX, HZ, 20)
    upper_X = check_witness(reflect(w), n, HZ, HX, 20)
    detector = pilot['explicit_nontriviality_detector']
    x = rotate(bits(rec['X_seed_supports'][detector['X_seed_index']], n), n, 48, detector['translation'])
    assert support(x) == detector['support']
    assert all((x & h).bit_count() % 2 == 0 for h in HZ)
    assert (x & w).bit_count() % 2 == 1
    wd = HERE / 'provenance/weight20_witness'
    wr = read(wd / 'receipt.json')
    assert sha(wd / 'receipt.json') == pilot['witness_receipt_sha256']
    assert sha(wd / 'stdout.jsonl') == wr['stdout_sha256']
    we = [json.loads(l) for l in (wd / 'stdout.jsonl').read_text().splitlines()]
    assert wr['returncode'] == 10 and wr['final']['status'] == 'witness'
    assert we[-1] == wr['final'] and any(e.get('support') == upper['support'] for e in we)
    for key, value in expected.items():
        assert wr[key] == value

    records = []
    frontier = set()
    for sid in range(64):
        d = HERE / 'shards' / f'{sid:02d}'
        r = read(d / 'receipt.json')
        assert sha(d / 'stdout.jsonl') == r['stdout_sha256']
        assert not (d / 'stderr.txt').read_text().strip()
        e = [json.loads(line) for line in (d / 'stdout.jsonl').read_text().splitlines()]
        f = validate_record(sid, r, e, expected)
        frontier.add((f['frontier_seen'], f['frontier_digest'], f['common_nodes']))
        records.append(dict(shard=sid, receipt=f'shards/{sid:02d}/receipt.json',
                            receipt_sha256=sha(d / 'receipt.json'), stdout_sha256=r['stdout_sha256'],
                            nodes=f['nodes'], seconds=f['seconds'], frontier_owned=f['frontier_owned']))
    assert len(frontier) == 1
    seen, digest, common = next(iter(frontier))
    assert sum(r['frontier_owned'] for r in records) == seen
    assert seen == 1780 and digest == '1378211316050488293'
    assert sorted(r['shard'] for r in records) == list(range(64))
    for reused in meta['reused_shards']:
        assert records[reused['shard']]['receipt_sha256'] == reused['receipt_sha256']
    return dict(status='CERTIFIED_EXACT', code_tag=meta['code_tag'], n=n, k=192, P=48,
                distance_exact=20, distance_lower=20, distance_upper=20,
                maximum_check_weight=9, maximum_column_weight=4, reported_frame_width=30,
                frame_width_optimality_proved=False,
                independent_algebra=algebra,
                pivot_determinant=rec['delta'], pairing=rec['pairing'], pairing_inverse=rec['pairing_inverse'],
                witness_Z=upper, witness_X=upper_X, nontriviality_detector=detector,
                reflection_exchanges_sectors=True, excluded_through=19,
                completed_shards=64, translation_roots_per_shard=16,
                frontier_seen=seen, frontier_digest=digest, common_nodes_per_shard=common,
                total_nodes=sum(r['nodes'] for r in records),
                sum_worker_seconds=sum(r['seconds'] for r in records),
                maximum_shard_seconds=max(r['seconds'] for r in records),
                input_hashes=meta['inputs'], records=records,
                evidence_kind='Completed exhaustive executions with independently audited physical algebra and partition coverage; not a formal node-by-node proof trace.',
                deep_search_replayed=False)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--report', type=Path)
    p.add_argument('--replay', action='store_true')
    p.add_argument('--seconds', type=float, default=1800)
    args = p.parse_args()
    result = check()
    if args.replay:
        with tempfile.TemporaryDirectory(prefix='pp-p48-replay-') as tmp:
            exe = Path(tmp) / 'exact'
            subprocess.run(['c++', '-O3', '-std=c++17', str(HERE / 'vendor/exact.cpp'), '-o', str(exe)], check=True)
            p = subprocess.run([str(exe), str(HERE / 'data/execution.mat'), '48', '19', str(args.seconds),
                                '--shard', '0', '1', '4', '18446744073709551615'],
                               capture_output=True, text=True)
            assert p.returncode in [0, 75], p.stdout + p.stderr
            final = json.loads(p.stdout.splitlines()[-1])
            result['replay_result'] = final
            result['deep_search_replayed'] = (p.returncode == 0 and final['status'] == 'excluded'
                                             and final['roots_completed'] == 16)
    if args.report:
        args.report.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: result[k] for k in ['status', 'n', 'k', 'distance_exact', 'maximum_check_weight',
                     'maximum_column_weight', 'reported_frame_width', 'completed_shards', 'total_nodes',
                     'sum_worker_seconds', 'deep_search_replayed']}, indent=2))
    if args.replay and not result['deep_search_replayed']:
        raise SystemExit(75)


if __name__ == '__main__':
    main()
