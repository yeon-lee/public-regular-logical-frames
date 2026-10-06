#!/usr/bin/env python3
"""Audit the flagship parent CSS distance and all completed W=23 shards.

This checker independently validates finite algebra and coverage/identity of
execution records. It does not re-execute the large exhaustive search trees.
"""
from pathlib import Path
import argparse, hashlib, json, zipfile
HERE=Path(__file__).resolve().parent

def need(x,msg):
    if not x:raise ValueError(msg)
def digest(data):return hashlib.sha256(data).hexdigest()
def sha(p):return digest(Path(p).read_bytes())
def load(p):return json.loads(Path(p).read_text())
def bits(support,n):
    need(len(support)==len(set(support)) and all(0<=q<n for q in support),'Invalid support')
    return sum(1<<q for q in support)
def basis(rows):
    out={}
    for v in rows:
        while v:
            q=v.bit_length()-1
            if q not in out:out[q]=v;break
            v^=out[q]
    return out
def rem(v,b):
    while v:
        q=v.bit_length()-1
        if q not in b:break
        v^=b[q]
    return v
def same_span(a,b):
    aa,bb=basis(a),basis(b)
    return len(aa)==len(bb) and all(rem(v,bb)==0 for v in a)
def binary_matrix(path):
    lines=Path(path).read_text().splitlines();n,m,g=map(int,lines[0].split());rows=[]
    for line in lines[1:]:
        a=list(map(int,line.split()));need(a[0]==len(a)-1,'Invalid binary row');rows.append(bits(a[1:],n))
    need(len(rows)==m+g,'Wrong row count')
    return n,rows[:m],rows[m:]
def phi(q):
    block,t=divmod(q,64)
    return ((block+8)%16)*64+(-t)%64

def parent_audit():
    p=HERE/'records';manifest=load(p/'bundle.json')
    for rel,want in manifest['files'].items():need(sha(p/rel)==want,'Parent original hash '+rel)
    build=load(p/'results/build.json')
    need(build['bundle_sha256']==sha(p/'bundle.json') and build['source_sha256']==sha(p/'exact.cpp') and build['matrix_sha256']==sha(p/'search.mat'),'Parent build binding')
    n,HX,HZ=binary_matrix(p/'data/code.mat');sn,SX,SZ=binary_matrix(p/'search.mat')
    need(n==sn==1024 and same_span(HX,SX) and same_span(HZ,SZ),'Parent search checks changed')
    need(len(basis(HX))==len(basis(HZ))==384,'Parent rank')
    need(all((x&z).bit_count()%2==0 for x in HX for z in HZ),'Parent commutation')
    witness=load(p/'data/witness24.json');v=bits(witness['support'],n)
    need(v.bit_count()==24 and all((v&h).bit_count()%2==0 for h in HX) and rem(v,basis(HZ))!=0,'Parent weight24 witness')
    reflected=[bits([phi(q) for q in range(n) if (v>>q)&1],n) for v in HX]
    need(same_span(reflected,HZ),'Parent reflection does not exchange sectors')
    results=[]
    for i in range(64):
        file=p/f'results/W23/{i:03d}.json';r=load(file)
        need(r['bundle_sha256']==sha(p/'bundle.json') and r['matrix_sha256']==sha(p/'search.mat'),'Parent receipt input binding')
        need(r['returncode']==0 and r['status']=='excluded' and r['witnesses']==[],'Parent incomplete')
        need(r['executable_sha256']==build['executable_sha256'],'Parent worker binary differs from build')
        need(sha(file.with_suffix('.jsonl'))==r['log_sha256'],'Parent log digest')
        events=[json.loads(line) for line in file.with_suffix('.jsonl').read_text().splitlines()]
        need(events==[r['result']],'Parent log/result disagreement')
        z=r['result'];expect={'event':'result','status':'excluded','W':23,'partition':'hash_prefix_v1','shard':i,'shards':64,'split_weight':6,'roots_completed':16,'coset_pruning':True,'root':0,'branch_begin':0,'branch_end':64}
        for k,v in expect.items():need(z[k]==v,'Parent result setting '+k)
        cmd=r['command'];need(cmd[2:4]==['64','23'] and cmd[5:9]==['--shard',str(i),'64','6'],'Parent command')
        results.append(z)
    for k in ['frontier_seen','frontier_digest','common_nodes']:
        need(len({r[k] for r in results})==1,'Parent inconsistent frontiers')
    need(sum(r['frontier_owned'] for r in results)==results[0]['frontier_seen'],'Parent incomplete frontier')
    return HX,HZ

if __name__ == "__main__":
    HX,HZ=parent_audit()
    print(json.dumps({"status":"PASS","n":1024,"k":256,"distance":24,"completed_shards":64,"physical_witness_checked":True,"deep_search_replayed":False,"formal_proof_object":False},indent=2))
