#!/usr/bin/env python3
"""Run the preserved detector tests on a fresh temporary build, without altering evidence."""
from pathlib import Path
import argparse, hashlib, json, shutil, subprocess, sys, tempfile
HERE=Path(__file__).resolve().parent
ap=argparse.ArgumentParser();ap.add_argument('--compiler',default='c++');ap.add_argument('--output');args=ap.parse_args()
with tempfile.TemporaryDirectory() as tmp:
    tmp=Path(tmp)
    for name in ['exact_target.cpp','common.hpp','tests.py']:shutil.copyfile(HERE/'records'/name,tmp/name)
    subprocess.run([args.compiler,'-O3','-std=c++17',str(tmp/'exact_target.cpp'),'-o',str(tmp/'exact_target')],check=True)
    subprocess.run([sys.executable,str(tmp/'tests.py')],check=True,capture_output=True,text=True)
    report=json.loads((tmp/'validation.json').read_text())
    if report['status']!='passed_detector_target_completeness_tests':raise ValueError('Small tests failed')
    report.update(scope='Independent small-code enumeration; not replay of the deep 1024-qubit exclusion',compiler_flags=['-O3','-std=c++17'])
    text=json.dumps(report,indent=2)+'\n'
    if args.output:Path(args.output).write_text(text)
    print(text,end='')
