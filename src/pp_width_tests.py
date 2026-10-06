"""Independent exhaustive validation of detector-root normalization and pruning."""

if not __debug__:
    raise RuntimeError(
        'Verification requires Python assertions; run without -O/-OO and unset PYTHONOPTIMIZE.'
    )

import hashlib
import json
from pathlib import Path
import random
import subprocess
import tempfile
import time
HERE=Path(__file__).resolve().parent
rng=random.Random(261031);start=time.time();runs=0

def basis(rows):
    out={}
    for v in rows:
        while v:
            p=v.bit_length()-1
            if p not in out:out[p]=v;break
            v^=out[p]
    return out
def remainder(v,rows):
    for p,r in sorted(rows.items(),reverse=True):
        if(v>>p)&1:v^=r
    return v
def expand(mask,n,P):return sum(1<<j for j in range(n) if(mask>>(j//P))&1)

cases=[]
for _ in range(40):
    n=rng.randrange(5,12);H=[rng.randrange(1,1<<n) for _ in range(rng.randrange(1,n))]
    ker=[v for v in range(1<<n) if all(not((v&h).bit_count()%2) for h in H)]
    G=rng.sample(ker,min(len(ker),rng.randrange(4)));hb=basis(H)
    masks=[v for v in range(1,1<<n) if all(not((v&g).bit_count()%2) for g in G) and remainder(v,hb)]
    if not masks:continue
    detector=rng.choice(masks);target=[v for v in ker if(v&detector).bit_count()%2]
    cases.append((n,1,H,G,detector,min(v.bit_count() for v in target)))
for P in (1,2):
    n=7*P;H=[sum(1<<(j*P+t) for j in range(7) if((j+1)>>i)&1) for i in range(3) for t in range(P)];hb=basis(H)
    for mask in (7,25):
        d=expand(mask,n,P)
        if any((d&h).bit_count()%2 for h in H) or not remainder(d,hb):continue
        target=[v for v in range(1<<n) if all(not((v&h).bit_count()%2) for h in H) and(v&d).bit_count()%2]
        cases.append((n,P,H,H,mask,min(v.bit_count() for v in target)))
with tempfile.TemporaryDirectory() as td:
    path=Path(td)/'code.mat'
    for n,P,H,G,mask,d in cases:
        rows=[f'{n} {len(H)} {len(G)}']
        for v in H+G:
            S=[j for j in range(n) if(v>>j)&1];rows.append(' '.join(map(str,[len(S)]+S)))
        path.write_text('\n'.join(rows)+'\n');detector=expand(mask,n,P)
        for W in sorted(set((max(1,d-1),d))):
            for pruning in (True,False):
                for depth in (1,2,4):
                    results=[];hits=0
                    for shard in range(3):
                        cmd=[str(HERE/'exact_target'),str(path),str(P),str(W),'30','--shard',str(shard),'3',str(depth),str(10**10),'--detector',hex(mask)]
                        if not pruning:cmd.append('no-coset')
                        proc=subprocess.run(cmd,capture_output=True,text=True);assert proc.returncode in(0,10),(cmd,proc.stderr);runs+=1
                        events=[json.loads(s) for s in proc.stdout.splitlines()];r=events[-1];results.append(r)
                        for e in events[:-1]:
                            assert e['event']=='witness';v=sum(1<<j for j in e['support'])
                            assert v.bit_count()<=W and(v&detector).bit_count()%2 and all(not((v&h).bit_count()%2) for h in H);hits+=1
                    assert bool(hits)==(d<=W)
                    if not hits:
                        assert all(r['roots_completed']==mask.bit_count() for r in results)
                        assert len({r['frontier_seen'] for r in results})==1 and len({r['frontier_digest'] for r in results})==1
                        assert sum(r['frontier_owned'] for r in results)==results[0]['frontier_seen']
result=dict(status='passed_detector_target_completeness_tests',fixtures=len(cases),exact_runs=runs,seconds=time.time()-start,
            tests=['literal minimum detector-odd kernel weights','pruning on/off','detector-only root normalization','cyclic periods1/2','complete 3-shard prefix ownership/digests','prefix depths1/2/4'],
            hashes={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (HERE/'exact_target.cpp',HERE/'common.hpp',Path(__file__))})
(HERE/'validation.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
