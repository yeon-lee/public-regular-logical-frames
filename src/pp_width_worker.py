"""One auditable augmentation-functional exclusion shard."""

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
import time
HERE=Path(__file__).resolve().parent
ap=argparse.ArgumentParser();ap.add_argument('--shard',type=int,required=True);ap.add_argument('--shards',type=int,required=True)
ap.add_argument('--seconds',type=float,default=7100);ap.add_argument('--split-weight',type=int,default=4)
ap.add_argument('--node-cap',type=int,default=10**15);ap.add_argument('--output',default='results');args=ap.parse_args()
assert 0<=args.shard<args.shards
manifest=json.loads((HERE/'bundle_sha256.json').read_text())
for name,digest in manifest['files'].items():assert hashlib.sha256((HERE/name).read_bytes()).hexdigest()==digest,name
data=json.loads((HERE/'data/detector.json').read_text());folder=Path(args.output)
if not folder.is_absolute():folder=HERE/folder
folder.mkdir(parents=True,exist_ok=True);stem=f'target_W24_shard_{args.shard:05d}';log=folder/(stem+'.jsonl')
cmd=[str(HERE/'exact_target'),str(HERE/'data/search.mat'),'64','24',str(args.seconds),'--shard',str(args.shard),str(args.shards),str(args.split_weight),str(args.node_cap),'--detector',hex(data['detector']['mask'])]
start=time.time()
with log.open('w') as f:proc=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT)
events=[]
for line in log.read_text().splitlines():
    try:events.append(json.loads(line))
    except json.JSONDecodeError:pass
result=next((r for r in reversed(events) if r.get('event')=='result'),None)
receipt=dict(claim_type='augmentation_functional_exclusion_not_distance',argv=cmd,returncode=proc.returncode,
             wall_seconds=time.time()-start,result=result,source_and_input_sha256=manifest['files'],
             log_sha256=hashlib.sha256(log.read_bytes()).hexdigest(),detector_mask=data['detector']['mask'])
path=folder/(stem+'.json');temp=path.with_suffix('.tmp');temp.write_text(json.dumps(receipt,indent=2)+'\n');os.replace(temp,path)
print(json.dumps(dict(shard=args.shard,returncode=proc.returncode,status=None if result is None else result['status'],seconds=receipt['wall_seconds'])),flush=True)
raise SystemExit(0 if proc.returncode in(0,10,75) else proc.returncode)
