#!/usr/bin/env python3
"""Compile and replay a complete enumeration; never infer completion from timeout.

Python 3.10+, a C++17 compiler, and the standard library suffice. Output goes to
a caller-selected directory; no compiled binary is retained in the package.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import resource
import shutil
import subprocess
import tempfile
import time

HERE = Path(__file__).resolve().parent



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


def normalize(coords, period):
    return min(sum(1 << (q // period * period + (q % period + h) % period)
                   for q in coords) for h in range(period))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--code', choices=['LPQ6', 'LPQ5'], required=True)
    ap.add_argument('--output', type=Path, required=True)
    ap.add_argument('--compiler', default=os.environ.get('CXX', 'c++'))
    ap.add_argument('--seconds', type=float, default=3600)
    args = ap.parse_args()
    name = args.code
    cap = {'LPQ6': 17, 'LPQ5': 20}[name]
    source = HERE / 'src' / ('kernel_minima_reverse.cpp' if name == 'LPQ6' else 'kernel_minima_native.cpp')
    inp = HERE / 'inputs' / (name + '.txt')
    args.output.mkdir(parents=True, exist_ok=True)
    output = args.output.resolve()
    for suffix in ['stdout', 'stderr', 'json']:
        if (output / (name + '.' + suffix)).exists():
            raise SystemExit('Refusing to overwrite an existing execution record; choose a new --output directory.')
    compiler = shutil.which(args.compiler)
    if compiler is None:
        raise SystemExit('C++17 compiler not found: ' + args.compiler)
    metadata = dict(code=name, cutoff=cap, input_sha256=sha(inp), source_sha256=sha(source),
                    compiler_flags=['-O3', '-std=c++17'], single_thread=True,
                    stabilizer_overlap_pruning=(name != 'LPQ6'),
                    root_range=[0, 68], first_child_partition=None, modular_partition=None,
                    evidence_kind='Completed exhaustive execution record, not a formal proof object.')
    with tempfile.TemporaryDirectory(prefix='regular-pairing-') as temp:
        binary = Path(temp) / 'enumerate'
        subprocess.run([compiler, '-O3', '-std=c++17', str(source), '-o', str(binary)], check=True)
        metadata['binary_sha256'] = sha(binary)
        command = [str(binary), str(inp), str(cap), '0', '68', str(args.seconds)]
        metadata['portable_command'] = ['ENUMERATE', f'inputs/{name}.txt', str(cap), '0', '68', str(args.seconds)]
        before = resource.getrusage(resource.RUSAGE_CHILDREN)
        start = time.monotonic()
        try:
            with (output / (name + '.stdout')).open('w') as out, (output / (name + '.stderr')).open('w') as err:
                result = subprocess.run(command, stdout=out, stderr=err, timeout=args.seconds + 60)
            metadata['returncode'] = result.returncode
        except subprocess.TimeoutExpired:
            metadata['returncode'] = 'external timeout'
        metadata['wall_seconds'] = time.monotonic() - start
        after = resource.getrusage(resource.RUSAGE_CHILDREN)
        metadata['user_seconds'] = after.ru_utime - before.ru_utime
        metadata['system_seconds'] = after.ru_stime - before.ru_stime
    lines = (output / (name + '.stdout')).read_text().splitlines()
    terminal = lines[-1] if lines else ''
    metadata['terminal'] = terminal
    pattern = (rf'ENUMERATION_COMPLETE W={cap} roots=\[0,68\) '
               r'children=\[0,1073741824\) mod\(D=-1,M=1,r=0\) nodes=(\d+) t=([\d.]+)')
    match = re.fullmatch(pattern, terminal)
    metadata['complete'] = bool(match and metadata['returncode'] == 0)
    metadata['nodes'] = int(match.group(1)) if match else None
    masks = set()
    for line in lines:
        if line.startswith('WORD '):
            values = list(map(int, line.split()[1:]))
            if not (values[0] == len(values[1:]) == len(set(values[1:])) and values[0] <= cap):
                raise ValueError('Malformed word in search output')
            masks.add(normalize(values[1:], 16))
    expected = json.loads((HERE / 'data' / (name + '_pool.json')).read_text())
    metadata['normalized_orbits'] = len(masks)
    metadata['matches_released_pool'] = masks == {sum(1 << q for q in support) for support in expected['representatives']}
    metadata['stdout_sha256'] = sha(output / (name + '.stdout'))
    metadata['stderr_sha256'] = sha(output / (name + '.stderr'))
    (output / (name + '.json')).write_text(json.dumps(metadata, indent=2) + '\n')
    print(json.dumps(metadata, indent=2))
    if not metadata['complete'] or not metadata['matches_released_pool']:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
