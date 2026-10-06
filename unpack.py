#!/usr/bin/env python3
"""Expand paper proof data into a new directory outside this release.

Source files and proof records are stored once. Expanded physical check matrices
are reconstructed from the catalogue and must match their recorded SHA-256.
No discovery data or additional constructions are generated.
"""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path, PurePosixPath
import sys

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parent

def digest(data):
    return hashlib.sha256(data).hexdigest()

def read(path):
    return json.loads(path.read_text())

def safe_path(root, name):
    if not isinstance(name, str):
        raise ValueError('Package paths must be strings')
    parts = PurePosixPath(name)
    if (not name or not parts.parts or parts.is_absolute()
            or '..' in parts.parts or '\\' in name or str(parts) != name):
        raise ValueError(f'Invalid package path: {name!r}')
    path = root.joinpath(*parts.parts)
    if not path.resolve().is_relative_to(root.resolve()):
        raise ValueError(f'Path escapes destination: {name!r}')
    return path

def check_integrity(root=ROOT):
    expected = read(root / 'SHA256SUMS.json')['files']
    actual = set()
    for path in root.rglob('*'):
        rel = path.relative_to(root)
        if '.git' in rel.parts or '__pycache__' in rel.parts or path.name == '.DS_Store':
            continue
        if path.is_symlink():
            raise ValueError(f'Symlink is outside release format: {path}')
        if path.is_file() and path != root / 'SHA256SUMS.json':
            actual.add(str(path.relative_to(root)))
    if actual != set(expected):
        raise ValueError('Release inventory differs from SHA256SUMS.json')
    for name, want in expected.items():
        if digest(safe_path(root, name).read_bytes()) != want:
            raise ValueError(f'Integrity mismatch: {name}')
    return len(expected)

def json_bytes(value):
    return (json.dumps(value, indent=2) + '\n').encode()

def expand(destination, root=ROOT):
    if sys.version_info < (3, 10):
        raise ValueError('Python 3.10 or newer is required')
    check_integrity(root)
    destination = Path(destination).resolve()
    if destination.is_relative_to(root.resolve()):
        raise ValueError('Choose a new directory outside this public release')
    if destination.exists():
        raise ValueError('Destination already exists; choose a new directory')
    layout = read(root / 'layout.json')
    if layout['schema'] != 'paper-proof-layout-1':
        raise ValueError('Unrecognized package layout')
    destination.mkdir(parents=True)
    written = set()

    def put(name, data, expected=None):
        path = safe_path(destination, name)
        if name in written:
            raise ValueError(f'Duplicate expanded path: {name}')
        if expected is not None and digest(data) != expected:
            raise ValueError(f'Record integrity mismatch: {name}')
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(data)
        written.add(name)

    for source in layout['sources']:
        data = safe_path(root, source['source']).read_bytes()
        for name in source['paths']:
            put(name, data, source['sha256'])
    for bundle_name in layout['proof_bundles']:
        bundle = read(safe_path(root, bundle_name))
        if bundle['schema'] != 'paper-proof-records-1':
            raise ValueError(f'Unrecognized proof bundle: {bundle_name}')
        for record in bundle['records']:
            if record['format'] == 'json':
                data = json_bytes(record['content'])
            elif record['format'] == 'text':
                data = record['content'].encode()
            else:
                raise ValueError('Unrecognized record encoding')
            for name in record['paths']:
                put(name, data, record['sha256'])
    catalogue = read(root / 'catalogue.json')
    if catalogue['schema'] != 'paper-constructions-1':
        raise ValueError('Unrecognized construction format')
    put('catalog.json', json_bytes(catalogue['catalogue']))
    spec = importlib.util.spec_from_file_location('paper_rlf', destination / 'rlf.py')
    algebra = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(algebra)
    for record in catalogue['constructions']:
        code = record['code']
        put(record['path'], json_bytes(code))
        put(record['frame_path'], json_bytes(record['frame']))
        hx, hz, _ = algebra.reconstruct(code)
        checks = dict(record['checks_template'])
        checks.update(X=[algebra.support(row) for row in hx],
                      Z=[algebra.support(row) for row in hz])
        put(record['checks_path'], json_bytes(checks), record['checks_sha256'])
    if digest((destination / 'SHA256SUMS.json').read_bytes()) != layout['expanded_manifest_sha256']:
        raise ValueError('Expanded inventory identity differs from layout')
    check_integrity(destination)
    return destination

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('destination', type=Path, help='New directory, outside the public release')
    args = parser.parse_args()
    result = expand(args.destination)
    print(f'Expanded and checked the paper proof workspace: {result}')
    print('See docs/VERIFICATION.md for per-claim replay commands.')

if __name__ == '__main__':
    try:
        main()
    except (ValueError, OSError, KeyError) as exc:
        raise SystemExit(str(exc)) from exc
