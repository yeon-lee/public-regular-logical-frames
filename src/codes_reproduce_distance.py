#!/usr/bin/env python3
"""Export a verified full-code exclusion input and optionally rerun its search.

This conservative binary search is portable across the catalogue. Specialized
PP searches can use stronger pruning and may be substantially faster. A timeout
never establishes a lower bound. Python 3.10+ and, for --run, a C++17 compiler.
"""
import argparse
import hashlib
import json
from pathlib import Path
import re
import subprocess
import sys
import time
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from rlf import reconstruct, support, translate, verify
from comparisons.check import verify_comparison

SOURCE = ROOT/'evidence/records/construction_records/claude_lp_frames/binary_audit/excl.cpp'


def completed(stdout, stderr, cutoff, blocks, returncode):
    expected = rf'COMPLETE W={cutoff} roots=\[0,{blocks}\) children=\[0,1073741824\) mod\(D=-1,M=1,r=0\) nodes=\d+ t=[0-9.]+'
    terminal = stdout.strip().splitlines()
    roots = [int(x) for x in re.findall(r'^root block (\d+) done:', stderr, re.M)]
    return returncode == 0 and len(terminal) == 1 and re.fullmatch(expected, terminal[0]) is not None and roots == list(range(blocks))


def export_input(entry, output):
    path = ROOT/entry['path']
    code = json.loads(path.read_text())
    read = lambda key: json.loads((path.parent/code['files'][key]).read_text())
    checks, frame = read('checks'), read('frame')
    hx, hz, seeds = reconstruct(code)
    P = code['group']['order']
    if code['collection'] == 'comparisons':
        detectors = read('detectors')
        verify_comparison(code, checks, frame, detectors)
        logical_rows = detectors['X']
    else:
        verify(code, checks, frame)
        logical_rows = [support(translate(x, g, P)) for x, z in seeds for g in range(P)]
    lines = [f'{code["n"]} {P} {code["n"]//P}']
    for rows in [list(map(support, hx)), list(map(support, hz)), logical_rows]:
        lines.append(str(len(rows)))
        lines.extend(' '.join(map(str, row)) for row in rows)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text('\n'.join(lines)+'\n')
    return code


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('code', help='Stable catalogue ID')
    parser.add_argument('--output', required=True, type=Path, help='Output input file; run logs use its stem')
    parser.add_argument('--run', action='store_true')
    parser.add_argument('--seconds', type=float, default=3600)
    parser.add_argument('--compiler', default='c++')
    a = parser.parse_args()
    if a.seconds <= 0:
        parser.error('--seconds must be positive')
    cat = json.loads((ROOT/'catalog.json').read_text())
    entry = next((e for e in cat['entries']+cat.get('comparisons', []) if e['id'] == a.code), None)
    if entry is None:
        parser.error('Unknown code')
    code = export_input(entry, a.output)
    print('Verified physical input:', a.output)
    if not a.run:
        return
    if code['distance']['lower'] is None:
        parser.error('This entry has no asserted exhaustive lower bound to reproduce')
    cutoff = code['distance']['lower']-1
    binary = a.output.with_suffix('.excl')
    build = [a.compiler, '-O3', '-std=c++17', str(SOURCE), '-o', str(binary)]
    subprocess.run(build, check=True)
    command = [str(binary.resolve()), str(a.output.resolve()), str(cutoff), '0', str(code['n']//code['group']['order']), str(a.seconds)]
    start = time.time()
    proc = subprocess.run(command, capture_output=True, text=True)
    outpath, errpath = a.output.with_suffix('.stdout'), a.output.with_suffix('.stderr')
    outpath.write_text(proc.stdout)
    errpath.write_text(proc.stderr)
    done = completed(proc.stdout, proc.stderr, cutoff, code['n']//code['group']['order'], proc.returncode)
    digest = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
    record = dict(code=a.code, status='COMPLETE' if done else 'INCOMPLETE_OR_FOUND', cutoff=cutoff,
                  distance_lower_bound=cutoff+1 if done else None, elapsed_seconds=time.time()-start,
                  sector='Z; independently verified sector isometry gives the X result',
                  returncode=proc.returncode, command=command, build_command=build,
                  input_sha256=digest(a.output), engine_sha256=digest(SOURCE),
                  stdout_sha256=digest(outpath), stderr_sha256=digest(errpath),
                  evidence_type='Exhaustive execution record, not a standalone formal proof object')
    a.output.with_suffix('.receipt.json').write_text(json.dumps(record, indent=2)+'\n')
    print(record['status'], proc.stdout.strip())
    if not done:
        raise SystemExit(2)


if __name__ == '__main__':
    main()
