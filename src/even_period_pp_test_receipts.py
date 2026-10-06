"""Reject incomplete, mismatched, or overlapping execution records."""
from copy import deepcopy
import json
import check


def main():
    d = check.HERE / 'shards/00'
    good = check.read(d / 'receipt.json')
    expected = {k: good[k] for k in ['execution_input_sha256', 'engine_sha256', 'source_sha256', 'header_sha256']}
    check.validate_record(0, good, [good['final']], expected)
    mutations = [
        ('status', lambda r: r['final'].update(status='incomplete')),
        ('cutoff', lambda r: r['final'].update(W=18)),
        ('roots', lambda r: r['final'].update(roots_completed=15)),
        ('duplicate-shard', lambda r: r['final'].update(shard=1)),
        ('partition-count', lambda r: r['final'].update(shards=63)),
        ('partition-depth', lambda r: r['final'].update(split_weight=5)),
        ('returncode', lambda r: r.update(returncode=75)),
        ('input-hash', lambda r: r.update(execution_input_sha256='0'*64)),
        ('source-hash', lambda r: r.update(source_sha256='0'*64)),
        ('witness', lambda r: r.update(witnesses=[{'weight': 19}])),
        ('command-cutoff', lambda r: r['command'].__setitem__(3, '18')),
    ]
    passed = []
    for name, mutate in mutations:
        r = deepcopy(good)
        mutate(r)
        try:
            check.validate_record(0, r, [r['final']], expected)
        except AssertionError:
            passed.append(name)
        else:
            raise AssertionError('Accepted invalid record: ' + name)
    try:
        check.validate_record(0, good, [dict(event='witness', weight=19), good['final']], expected)
    except AssertionError:
        passed.append('contradictory-transcript')
    else:
        raise AssertionError('Accepted contradictory transcript')
    print(json.dumps(dict(status='PASS', rejected_invalid_records=passed)))


if __name__ == '__main__':
    main()
