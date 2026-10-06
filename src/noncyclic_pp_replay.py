#!/usr/bin/env python3
"""Explicitly rerun the noncyclic code's 16 covered native roots.

Default cutoff18 reproduces the distance lower bound and can be expensive.
--cutoff2 is only a quick compiled-engine smoke check, not a distance19 proof.
All compilation and replay files go outside the package.
"""

if not __debug__:
    raise RuntimeError(
        'Verification requires Python assertions; run without -O/-OO and unset PYTHONOPTIMIZE.'
    )

import argparse
import json
import os
from pathlib import Path
import subprocess
import tempfile

HERE = Path(__file__).resolve().parent
ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument("--cutoff", type=int, default=18)
ap.add_argument("--seconds", type=float, default=3600, help="time limit per root")
ap.add_argument("--output", type=Path, help="new output directory; default is a temporary directory")
ap.add_argument("--cxx", default=os.environ.get("CXX", "c++"))
args = ap.parse_args()
assert 1 <= args.cutoff <= 18 and args.seconds > 0
if args.output:
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=False)
else:
    out = Path(tempfile.mkdtemp(prefix="noncyclic-pp-replay-"))
binary = out / "exact"
build = [args.cxx, "-O3", "-std=c++17", str(HERE / "src/exact.cpp"), "-o", str(binary)]
subprocess.run(build, check=True)
results = []
for root in range(0, 32, 2):
    command = [str(binary), str(HERE / "completion/input.mat"), "16", str(args.cutoff),
               str(args.seconds), str(root), "0", "64", str(2**64 - 1)]
    run = subprocess.run(command, capture_output=True, text=True)
    stem = out / f"root{root:02d}"
    stem.with_suffix(".stdout.jsonl").write_text(run.stdout)
    stem.with_suffix(".stderr").write_text(run.stderr)
    events = [json.loads(s) for s in run.stdout.splitlines()]
    result = dict(root=root, command=command, returncode=run.returncode, events=events)
    stem.with_suffix(".json").write_text(json.dumps(result, indent=2) + "\n")
    results.append(result)
    if run.returncode or len(events) != 1 or events[0].get("status") != "excluded":
        raise SystemExit(f"Root {root} did not complete exclusion; see {out}")
summary = dict(status="COMPLETE", cutoff=args.cutoff, roots=list(range(0, 32, 2)),
               nodes=sum(x["events"][0]["nodes"] for x in results),
               seconds=sum(x["events"][0]["seconds"] for x in results),
               build_command=build,
               output_directory=str(out),
               scope="Distance >=19 replay" if args.cutoff == 18 else "Smoke check only; not the distance19 lower-bound replay")
(out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
print(json.dumps(summary, indent=2))
