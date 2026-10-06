#!/usr/bin/env python3
"""Run one fresh complete certificate attempt; never overwrite shard evidence."""

if not __debug__:
    raise RuntimeError(
        'Verification requires Python assertions; run without -O/-OO and unset PYTHONOPTIMIZE.'
    )

from concurrent.futures import ThreadPoolExecutor,as_completed
import hashlib,json,subprocess,time
from pathlib import Path

HERE=Path(__file__).resolve().parent

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

def sha(p):
    """File digest or explicitly declared omitted executable identity."""
    return evidence_digest(p)
SOURCE={p.name:sha(p) for p in (HERE/'engine/exact.cpp',HERE/'engine/common.hpp')}
EXPECTED={'exact.cpp':'9f733fcdc935f7011b3128216cecc58ed8cb3949e7d422431a85ab5cb5f152b6',
          'common.hpp':'e9233ec63460a789be4bcaa9ecc37075bd59470773c09a898ab968daf1bad667'}
assert SOURCE==EXPECTED
EXE=HERE/'engine/exact';MATRIX=HERE/'inputs/exact_Z.mat'
INPUT_HASH=sha(MATRIX);EXE_HASH=sha(EXE)
assert INPUT_HASH=='8775a4f953b045f77f8942734bf78466cc048e4fd2d5167e3bbab8ae72b9966d'

def run(shard):
    prefix=HERE/'receipts'/f'shard_{shard:02d}'
    stdout=prefix.with_suffix('.stdout');stderr=prefix.with_suffix('.stderr');receipt=prefix.with_suffix('.receipt.json')
    assert not any(p.exists() for p in (stdout,stderr,receipt)),f'Existing shard {shard}'
    cmd=[str(EXE),str(MATRIX),'16','21','300','--shard',str(shard),'64','5','1000000000000']
    start=time.monotonic()
    with stdout.open('x') as out,stderr.open('x') as err:
        process=subprocess.run(cmd,stdout=out,stderr=err,check=False)
    elapsed=time.monotonic()-start
    result={'command':cmd,'returncode':process.returncode,'elapsed_seconds':elapsed,
            'input_sha256':INPUT_HASH,'source_sha256':SOURCE,'executable_sha256':EXE_HASH,
            'stdout_sha256':sha(stdout),'stderr_sha256':sha(stderr),'driver_sha256':sha(__file__)}
    with receipt.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
    events=[json.loads(line) for line in stdout.read_text().splitlines() if line.strip()]
    final=events[-1] if events else {}
    return {'shard':shard,'returncode':process.returncode,'status':final.get('status'),
            'seconds':elapsed,'roots_completed':final.get('roots_completed'),'nodes':final.get('nodes')}

if __name__=='__main__':
    completed=[];start=time.monotonic()
    with ThreadPoolExecutor(max_workers=4) as pool:
        futures=[pool.submit(run,shard) for shard in range(64)]
        for future in as_completed(futures):
            result=future.result();completed.append(result)
            print(json.dumps(dict(event='shard_complete',completed=len(completed),**result)),flush=True)
    summary={'workers':4,'wall_seconds':time.monotonic()-start,'runs':sorted(completed,key=lambda r:r['shard'])}
    (HERE/'run_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps({'event':'all_complete','wall_seconds':summary['wall_seconds'],'excluded':sum(r['status']=='excluded' for r in completed)}),flush=True)
