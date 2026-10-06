#!/usr/bin/env python3
"""Verify catalogue integrity and every stored code/frame; no distance search."""
import argparse
import hashlib
import json
from pathlib import Path
import sys

if sys.version_info < (3, 10):
    raise SystemExit('Python 3.10 or newer is required.')
from rlf import verify
from comparisons.check import run as verify_comparison

ROOT = Path(__file__).resolve().parent


def read(path):
    return json.loads(path.read_text())


def check_integrity():
    manifest = read(ROOT / 'SHA256SUMS.json')
    for name, expected in manifest['files'].items():
        actual = hashlib.sha256((ROOT / name).read_bytes()).hexdigest()
        if actual != expected:
            raise ValueError(f'Integrity mismatch: {name}')
    return len(manifest['files'])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--code', help='Verify one catalogue ID; default is all entries')
    parser.add_argument('--report', type=Path, help='Optional JSON report destination')
    args = parser.parse_args()
    count = check_integrity()
    catalogue = read(ROOT / 'catalog.json')
    entries = catalogue['entries'] + catalogue.get('comparisons', [])
    if not entries or len({e['id'] for e in entries}) != len(entries):
        raise ValueError('Expected nonempty catalogue with unique IDs')
    if args.code:
        entries = [e for e in entries if e['id'] == args.code]
        if not entries:
            parser.error('Unknown catalogue ID')
    reports = []
    for entry in entries:
        path = ROOT / entry['path']
        code = read(path)
        for key in ('id', 'n', 'k', 'distance', 'frame_width', 'paper', 'collection'):
            if entry[key] != code[key]:
                raise ValueError(f'Catalogue metadata mismatch: {entry["id"]}/{key}')
        for evidence in code['evidence']:
            if not (ROOT / evidence).is_file():
                raise ValueError(f'Missing evidence: {evidence}')
        if entry['collection'] == 'comparisons':
            result = verify_comparison(entry)
        else:
            result = verify(code, read(path.parent / code['files']['checks']),
                            read(path.parent / code['files']['frame']))
        reports.append(result)
        width_label = 'partial_width' if entry['collection'] == 'comparisons' else 'width'
        width = result.get('frame_width', result.get('partial_frame_width'))
        print(f'{entry["id"]}: PASS (n={code["n"]}, k={code["k"]}, {width_label}={width})', flush=True)
    report = {'status': 'PASS', 'integrity_files_checked': count, 'codes': reports,
              'distance_searches_rerun': False,
              'distance_scope': 'Upper witnesses checked afresh; exact lower bounds are the preserved manuscript certificates or analytical proof.'}
    if args.report:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(json.dumps(report, indent=2) + '\n')
    print(f'PASS: {len(reports)} codes and {count} file hashes. No exhaustive distance searches rerun.')


if __name__ == '__main__':
    try:
        main()
    except (ValueError, KeyError, OSError) as exc:
        raise SystemExit(str(exc)) from exc
