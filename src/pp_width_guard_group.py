if not __debug__:
    raise RuntimeError(
        'Verification requires Python assertions; run without -O/-OO and unset PYTHONOPTIMIZE.'
    )

from pathlib import Path
import hashlib,json,subprocess,sys,time
H=Path(__file__).resolve().parent;group=int(sys.argv[1]);assert 0<=group<4
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
b=json.loads((H/'build.json').read_text());assert sha(H/'exact_target')==b['executable_sha256']
assert sha(H/'bundle_sha256.json')==b['bundle_sha256']
out=H/'results';out.mkdir(exist_ok=True)
assert not list(out.glob(f'group_{group:02d}.*'))
ids=range(16*group,16*group+16)
assert all(not(out/f'target_W24_shard_{i:05d}.jsonl').exists() for i in ids)
start=time.time()
args=[sys.executable,str(H/'grouped.py'),'--start',str(16*group),'--count','16','--shards','64','--seconds','6000','--split-weight','4','--output','results']
p=subprocess.run(args,cwd=H)
assert sha(H/'exact_target')==b['executable_sha256']
r=dict(group=group,shards=list(ids),returncode=p.returncode,wall_seconds=time.time()-start,build=b,argv=args)
(out/f'group_{group:02d}.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r))
raise SystemExit(p.returncode)
