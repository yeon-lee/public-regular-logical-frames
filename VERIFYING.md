# Format and evidence

`verify.py` first checks the public file inventory. It then uses `unpack.py` to
create a temporary proof workspace and run the retained finite verification
suite. No original research directory, cluster account, or third-party Python
package is required. `--report /path/outside/the/release.json` saves the result.

Run Python without `-O` or `-OO`, and leave `PYTHONOPTIMIZE` unset. Tools that
depend on assertions reject optimized Python explicitly, including direct
calls to individual checkers.

`catalogue.json` contains construction records and explicit frame seeds. The
binary check matrices are reproducible expansions: each regenerated file must
match the included checksum before any checker runs. `layout.json` maps shared
source files to the paths expected by the checkers. Each `proofs/*.json` file
contains exact records with their paths, SHA-256, and content. JSON records are
restored with Python `json.dumps(value, indent=2) + "\n"`; text records are
restored verbatim. Every expanded file must match the expanded inventory too.
The bundles are plain JSON, not opaque compressed archives.

The sources in `src/` are inspected directly. Run their commands from an
unpacked workspace so that their data paths resolve. For example:

```sh
python3 unpack.py /tmp/regular-frames-proof
cd /tmp/regular-frames-proof
python3 verification/regular_pairing/check.py
python3 verification/regular_pairing/replay.py --code LPQ6 --output /tmp/lp-pool-replay --seconds 1800
```

The last command freshly enumerates a complete pool; it is not part of the
routine finite suite. Other package READMEs give the corresponding distance
and detector-exclusion commands. An incomplete run supplies no new lower bound.

Physical witnesses, canonical frames, detector functionals, and proper graph
colorings are finite objects checked independently. Pool completeness and
exhaustive lower bounds also rely on the supplied algorithms and completed
execution records. The records bind inputs, source identities, terminal
results, and full root/partition coverage; they are not formal proofs checking
every visited search node. Some interrupted-prefix records remain because only
their completed subtrees enter a checked coverage union.

The regular LP pairing bounds apply to the specified cyclic action. The PP
frame optimum uses a detector-targeted exclusion distinct from its ordinary
distance exclusion. Other frame widths are achieved upper bounds unless an
optimum is explicitly stated. The published odd-period comparison has a partial
regular frame and complete logical detectors; its partial width does not bound
full-code frame width.

Public records omit calendar dates, scheduler and workstation metadata, and
candidate-discovery history. Elapsed computational costs and required partition
indices are retained. Public hashes bind these normalized records, not untouched
private logs. Omitted executable fingerprints are declarations, not checks of
binary contents or proof of a source-to-binary build. Fresh replays compile the
supplied C++17 source. Packaging makes no additional mathematical claim.
