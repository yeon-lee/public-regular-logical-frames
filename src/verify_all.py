#!/usr/bin/env python3
"""Audit packaged routes and run all finite checks; no deep search replay.

Requires Python >=3.10 and a C++17 compiler. Commands run from this checkout;
compiler tests use temporary directories. --report optionally writes a JSON
record outside the checkout. Scientific inputs and evidence are never changed.
"""

if not __debug__:
    raise RuntimeError(
        'Verification requires Python assertions; run without -O/-OO and unset PYTHONOPTIMIZE.'
    )

import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time
import zipfile

ROOT = Path(__file__).resolve().parent
COMMANDS = [
    ['verify.py'],
    ['-m', 'unittest', 'discover', '-s', 'tests', '-v'],
    ['verification/regular_pairing/check.py'],
    ['verification/pp_width/check.py', '--replay-prefix'],
    ['verification/pp_distance/check.py'],
    ['verification/small_quaternary_pp/check.py'],
    ['verification/noncyclic_pp/check.py'],
    ['verification/lpq_extended/independent544/check.py'],
    ['verification/check8_pp/check.py'],
    ['verification/even_period_pp/check.py'],
    ['verification/even_period_pp/test_receipts.py'],
    ['verification/check8_pp_period48/check.py'],
    ['verification/pp48_w10_bound/check.py'],
    *[[str(p.relative_to(ROOT))] for p in sorted(ROOT.glob('codes/*/*/proof/replay.py'))],
    ['verification/pp_width/test_small.py'],
]


def read(path):
    return json.loads(path.read_text())


def packaged_path(name, manifest):
    path = ROOT / name
    if Path(name).is_absolute() or not path.resolve().is_relative_to(ROOT):
        raise ValueError(f'Nonportable data route: {name}')
    if not path.exists() or path.is_symlink():
        raise ValueError(f'Missing file or symlink: {name}')
    if path.is_file() and name not in manifest:
        raise ValueError(f'Unmanifested data route: {name}')
    return path


def audit_package():
    manifest = read(ROOT/'SHA256SUMS.json')['files']
    for name, expected in manifest.items():
        path = packaged_path(name, manifest)
        if hashlib.sha256(path.read_bytes()).hexdigest() != expected:
            raise ValueError(f'Integrity mismatch: {name}')
        if path.suffix.lower() in ('.pdf', '.bib', '.tex') and name != 'verification/pp_width/records/coverage.tex':
            raise ValueError(f'Manuscript/figure material is outside public package scope: {name}')
        if path.suffix == '.zip':
            with zipfile.ZipFile(path) as archive:
                for member in archive.namelist():
                    if Path(member).suffix.lower() in ('.pdf', '.bib', '.tex'):
                        raise ValueError(f'Manuscript material inside archive: {name}/{member}')
    guide = read(ROOT/'docs/data-map.json')
    catalogue = read(ROOT/'catalog.json')
    if {c['id'] for c in guide['codes']} != {c['id'] for c in catalogue['entries'] + catalogue['comparisons']}:
        raise ValueError('Data map does not cover the catalogue')
    for code in guide['codes']:
        for name in [code['code_record'], *code['files'].values(), *code['evidence']['referenced_paths']]:
            packaged_path(name, manifest)
        for key in ('construction', 'witnesses', 'evidence'):
            location = code[key]
            value = read(packaged_path(location['path'], manifest))
            for component in location['json_pointer'].strip('/').split('/'):
                component = component.replace('~1', '/').replace('~0', '~')
                value = value[int(component)] if isinstance(value, list) else value[component]
    for claim in guide['claims']:
        for name in [claim['readme'], *claim['data_routes']]:
            packaged_path(name, manifest)
    return {'manifest_files': len(manifest), 'catalogue_ids': len(guide['codes']),
            'claim_packages': len(guide['claims'])}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--report', type=Path)
    args = parser.parse_args()
    if sys.version_info < (3, 10):
        parser.error('Python 3.10 or newer is required')
    if args.report and args.report.resolve().is_relative_to(ROOT):
        parser.error('--report must be outside the checkout')
    report = {'status': 'PASS', 'package': audit_package(), 'commands': [],
              'manifest_sha256': hashlib.sha256((ROOT/'SHA256SUMS.json').read_bytes()).hexdigest(),
              'deep_searches_rerun': False, 'binary_contents_verified': False,
              'record_metadata_normalized': True}
    env = os.environ.copy()
    env.pop('PYTHONPATH', None)
    env['PYTHONDONTWRITEBYTECODE'] = '1'
    for command in COMMANDS:
        started = time.monotonic()
        result = subprocess.run([sys.executable, '-B', *command], cwd=ROOT, env=env,
                                capture_output=True, text=True)
        report['commands'].append({'command': ['python3', '-B', *command],
                                   'returncode': result.returncode,
                                   'elapsed_seconds': time.monotonic() - started,
                                   'stdout': result.stdout, 'stderr': result.stderr})
        print(f"{'PASS' if result.returncode == 0 else 'FAIL'}: {' '.join(command)}", flush=True)
        if result.returncode:
            report['status'] = 'FAIL'
            print(result.stdout + result.stderr, file=sys.stderr)
            break
    if args.report:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(json.dumps(report, indent=2) + '\n')
    print(f"{report['status']}: {len(report['commands'])} checks; no deep searches rerun.")
    return 0 if report['status'] == 'PASS' else 1


if __name__ == '__main__':
    sys.exit(main())
