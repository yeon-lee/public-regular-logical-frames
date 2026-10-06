#!/usr/bin/env python3
"""Check the paper-only release using a temporary expanded workspace.

Default: all finite proof audits and small-instance tests. Large exhaustive
searches are not repeated. Use unpack.py for explicit deep replay commands.
"""

if not __debug__:
    raise RuntimeError(
        'Verification requires Python assertions; run without -O/-OO and unset PYTHONOPTIMIZE.'
    )

import argparse
import os
from pathlib import Path
import subprocess
import sys
import tempfile

sys.dont_write_bytecode = True
from unpack import expand, ROOT

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--quick', action='store_true', help='Only physical catalogue reconstruction and witnesses')
    parser.add_argument('--code', help='Only one catalogue ID')
    parser.add_argument('--report', type=Path, help='Optional output report outside the public folder')
    args = parser.parse_args()
    if sys.version_info < (3, 10):
        parser.error('Python 3.10 or newer is required')
    report = args.report.resolve() if args.report else None
    if report and report.is_relative_to(ROOT):
        parser.error('--report must be outside the public folder')
    env = os.environ.copy()
    env.pop('PYTHONPATH', None)
    env['PYTHONDONTWRITEBYTECODE'] = '1'
    with tempfile.TemporaryDirectory(prefix='regular-frames-') as tmp:
        work = expand(Path(tmp) / 'proofs')
        command = [sys.executable, '-B', 'verify.py' if (args.quick or args.code) else 'verify_all.py']
        if args.code:
            command += ['--code', args.code]
        if report:
            command += ['--report', str(report)]
        result = subprocess.run(command, cwd=work, env=env)
    return result.returncode

if __name__ == '__main__':
    raise SystemExit(main())
