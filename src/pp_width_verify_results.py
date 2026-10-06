"""Aggregate detector-target receipts; never convert them to a distance bound."""

if not __debug__:
    raise RuntimeError(
        'Verification requires Python assertions; run without -O/-OO and unset PYTHONOPTIMIZE.'
    )

import argparse
import hashlib
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent
ap=argparse.ArgumentParser();ap.add_argument('--results',default='results');ap.add_argument('--shards',type=int,default=128)
ap.add_argument('--output',default='target_certificate.json');args=ap.parse_args()
folder=Path(args.results)
if not folder.is_absolute():folder=HERE/folder
manifest=json.loads((HERE/'bundle_sha256.json').read_text())
for name,digest in manifest['files'].items():assert hashlib.sha256((HERE/name).read_bytes()).hexdigest()==digest,name
data=json.loads((HERE/'data/detector.json').read_text());mask=data['detector']['mask'];P=data['P'];W=data['W']
lines=(HERE/'data/code.mat').read_text().splitlines();n,m,g=map(int,lines[0].split())
def word(S):
    assert len(set(S))==len(S) and all(0<=j<n for j in S);return sum(1<<j for j in S)
def read(s):
    x=list(map(int,s.split()));assert x[0]==len(x)-1;return word(x[1:])
rows=list(map(read,lines[1:]));H=rows[:m];G=rows[m:]
def basis(rows):
    b={}
    for v in rows:
        while v:
            p=v.bit_length()-1
            if p not in b:b[p]=v;break
            v^=b[p]
    return b
def remainder(v,b):
    for p,r in sorted(b.items(),reverse=True):
        if(v>>p)&1:v^=r
    return v
def translate(v,t):
    out=0;pmask=(1<<P)-1;t%=P
    for b in range(n//P):
        x=(v>>(P*b))&pmask;out|=(((x<<t)|(x>>(P-t)))&pmask)<<(P*b)
    return out
gb=basis(G);hb=basis(H);functional=sum(1<<j for j in range(n) if(mask>>(j//P))&1)
assert all(not((functional&v).bit_count()%2) for v in G) and remainder(functional,hb)
search_lines=(HERE/'data/search.mat').read_text().splitlines();sn,sm,sg=map(int,search_lines[0].split());assert sn==n
search_rows=list(map(read,search_lines[1:]));SH,SG=search_rows[:sm],search_rows[sm:];assert len(SG)==sg
shb,sgb=basis(SH),basis(SG);assert len(shb)==len(hb) and len(sgb)==len(gb)
assert all(not remainder(v,hb) for v in SH) and all(not remainder(v,shb) for v in H)
assert all(not remainder(v,gb) for v in SG) and all(not remainder(v,sgb) for v in G)
assert all(not((h&g).bit_count()%2) for h in H for g in G)
for supp in data['known_short_supports']:
    v=word(supp);assert v.bit_count()<=W and all(not((h&v).bit_count()%2) for h in H) and remainder(v,gb)
    assert not((v&functional).bit_count()%2)
attain=word(data['attaining_target_support']);assert attain.bit_count()==25 and(attain&functional).bit_count()%2
assert all(not((attain&h).bit_count()%2) for h in H)
upper=json.loads((HERE/'frame_upper_bound.json').read_text());X=[word(s) for s in upper['X_seed_supports']];Z=[word(s) for s in upper['Z_seed_supports']]
assert max(v.bit_count() for v in X+Z)==25
assert all(not((x&h).bit_count()%2) for x in X for h in G) and all(not((z&h).bit_count()%2) for z in Z for h in H)
XX=[translate(v,t) for v in X for t in range(P)];ZZ=[translate(v,t) for v in Z for t in range(P)]
assert all((x&z).bit_count()%2==(i==j) for i,x in enumerate(XX) for j,z in enumerate(ZZ))
assert len(basis(H+XX))-len(hb)==256 and len(basis(G+ZZ))-len(gb)==256
complete=[];missing=[];incomplete=[];witnesses=[];receipts=[]
for shard in range(args.shards):
    path=folder/f'target_W24_shard_{shard:05d}.json';log=path.with_suffix('.jsonl')
    if not path.exists():missing.append(shard);continue
    rec=json.loads(path.read_text());assert rec['claim_type']=='augmentation_functional_exclusion_not_distance'
    assert rec['source_and_input_sha256']==manifest['files'] and rec['detector_mask']==mask
    assert hashlib.sha256(log.read_bytes()).hexdigest()==rec['log_sha256']
    events=[json.loads(s) for s in log.read_text().splitlines() if s.startswith('{')];result=events[-1]
    assert result==rec['result'] and result['event']=='result'
    assert result['shard']==shard and result['shards']==args.shards and result['W']==W and result['detector_mask']==mask
    assert result['partition']=='hash_prefix_v1' and result['coset_pruning']
    for e in events[:-1]:
        if e.get('event')!='witness':continue
        v=word(e['support']);assert v.bit_count()==e['weight']<=W
        assert all(not((v&h).bit_count()%2) for h in H) and(v&functional).bit_count()%2 and remainder(v,gb)
        witnesses.append(dict(shard=shard,weight=v.bit_count(),support=e['support']))
    if result['status']=='excluded':
        assert rec['returncode']==0 and result['roots_completed']==data['expected_roots_completed'];complete.append(result)
    elif result['status']=='witness':assert rec['returncode']==10
    else:assert result['status']=='incomplete';incomplete.append(shard)
    receipts.append(dict(shard=shard,path=str(path),sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
certified=len(complete)==args.shards and not witnesses
if certified:
    for key in('frontier_seen','frontier_digest','common_nodes','split_weight'):assert len({r[key] for r in complete})==1,key
    assert sum(r['frontier_owned'] for r in complete)==complete[0]['frontier_seen']
report=dict(claim_type='augmentation_functional_exclusion_not_distance',status='complete_detector_exclusion' if certified else 'target_counterexample' if witnesses else 'incomplete',
            detector_mask=mask,weight_excluded_through=W,completed_shards=len(complete),total_shards=args.shards,
            missing_shards=missing,incomplete_shards=incomplete,witnesses=witnesses,receipts=receipts,
            completed_nodes=sum(r['nodes'] for r in complete),complete_width25_frame_independently_verified=True,
            bundle_sha256=hashlib.sha256((HERE/'bundle_sha256.json').read_bytes()).hexdigest(),verifier_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
if certified:report.update(exact_detector_distance=25,exact_complete_canonical_width=25,ordinary_distance_lower_bound_claimed=False,
                           frontier_seen=complete[0]['frontier_seen'],frontier_digest=complete[0]['frontier_digest'],frontier_owned=sum(r['frontier_owned'] for r in complete))
out=Path(args.output)
if not out.is_absolute():out=HERE/out
out.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items() if k not in('receipts','witnesses')},indent=2))
