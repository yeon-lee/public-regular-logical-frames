#!/usr/bin/env python3
"""Strict, independent union validator for the preserved hash_prefix_v1 engine.

No partial, timed-out, missing, or merely sampled run earns lower-bound credit.
Recomputes physical algebra from supplied matrices and constituent coefficients.
Checks retained source, receipts, stdout and stderr. Public executable identities
are declarations; fresh replay hashes its locally compiled executable.
"""
import argparse
import hashlib
import json
import math
import sys
from pathlib import Path

TRUSTED_SOURCE={
 'exact.cpp':'9f733fcdc935f7011b3128216cecc58ed8cb3949e7d422431a85ab5cb5f152b6',
 'common.hpp':'e9233ec63460a789be4bcaa9ecc37075bd59470773c09a898ab968daf1bad667'}
RHO={0:((0,0),(0,0)),1:((1,0),(0,1)),2:((0,1),(1,1)),3:((1,1),(1,0))}

def need(condition,message):
    if not condition:raise ValueError(message)

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
def load(p):return json.loads(Path(p).read_text())
def vec(supp,n):
    need(isinstance(supp,list),'Support must be a list')
    need(all(type(j) is int and 0<=j<n for j in supp),'Invalid coordinate')
    need(len(set(supp))==len(supp),'Repeated coordinate')
    return sum(1<<j for j in supp)
def basis(rows):
    out={}
    for v in rows:
        while v:
            j=v.bit_length()-1
            if j not in out:out[j]=v;break
            v^=out[j]
    return out
def outside(v,b):
    while v:
        j=v.bit_length()-1
        if j not in b:return True
        v^=b[j]
    return False
def equalspan(a,b):
    aa,bb=basis(a),basis(b)
    return len(aa)==len(bb) and all(not outside(v,bb) for v in aa.values())
def orthogonal(a,b):return all((x&y).bit_count()%2==0 for x in a for y in b)
def shifted(v,P,n,t=1):
    mask=(1<<P)-1;t%=P;out=0
    for b in range(n//P):
        a=(v>>(b*P))&mask;out|=(((a<<t)|(a>>(P-t)))&mask)<<(b*P)
    return out
def transpose(v,P,r,s):
    out=0
    while v:
        bit=v&-v;j=bit.bit_length()-1;v^=bit
        block,t=divmod(j,P);block,beta=divmod(block,2)
        if block<s*s:a,b=divmod(block,s);dst=b*s+a
        else:a,b=divmod(block-s*s,r);dst=s*s+b*r+a
        out|=1<<((2*dst+beta)*P+t)
    return out

def weighted(path,initials=False):
    vals=list(map(int,Path(path).read_text().split()));at=0
    def take():
        nonlocal at
        need(at<len(vals),'Truncated matrix');x=vals[at];at+=1;return x
    n,m,g=take(),take(),take();need(n>0 and m>=0 and g>=0,'Invalid matrix header')
    def row():
        w=take();need(0<=w<=n,'Invalid row weight');return vec([take() for _ in range(w)],n)
    h=[row() for _ in range(m)];gg=[row() for _ in range(g)];ii=[]
    if initials:
        count=take();need(count>=0,'Invalid witness count');ii=[row() for _ in range(count)]
    need(at==len(vals),'Unexpected matrix suffix')
    return n,h,gg,ii

def constituent(data,P):
    if isinstance(data,dict):
        C,E=data['C'],data['E'];need(len(C)==len(E),'Constituent row mismatch')
        need(all(len(c)==len(e)==len(C[0]) for c,e in zip(C,E)),'Constituent column mismatch')
        out=[]
        for cr,er in zip(C,E):
            row=[]
            for c,e in zip(cr,er):
                need(type(c) is int and 0<=c<=3 and type(e) is int,'Invalid F4 monomial')
                row.append({e%P:c} if c else {})
            out.append(row)
        return out
    need(isinstance(data,list) and len(data)>0,'Invalid polynomial constituent')
    out=[]
    for rr in data:
        row=[]
        for p in rr:
            need(len(p)==P and all(type(c) is int and 0<=c<=3 for c in p),'Invalid F4 polynomial')
            row.append({e:c for e,c in enumerate(p) if c})
        out.append(row)
    need(all(len(row)==len(out[0]) for row in out),'Ragged constituent')
    return out

def reconstruct(data,P,r,s):
    A=constituent(data.get('A',data),P);B=constituent(data.get('B',data.get('A',data)),P)
    need(len(A)==len(B)==r and all(len(a)==s for a in A+B),'Constituent shape mismatch')
    def bar(p):return {(-e)%P:c for e,c in p.items()}
    def expand(entries):
        out=[]
        for alpha in (0,1):
            for t in range(P):
                row=0
                for block,p in entries:
                    for e,c in p.items():
                        for beta in (0,1):
                            if RHO[c][alpha][beta]:row^=1<<((2*block+beta)*P+(t-e)%P)
                out.append(row)
        return out
    hx=[];hz=[]
    for i in range(r):
        for b in range(s):hx+=expand([(a*s+b,A[i][a]) for a in range(s)]+[(s*s+i*r+c,bar(B[c][b])) for c in range(r)])
    for a in range(s):
        for j in range(r):hz+=expand([(a*s+b,B[j][b]) for b in range(s)]+[(s*s+i*r+j,bar(A[i][a])) for i in range(r)])
    return hx,hz

def audit_physical(cdir,sector):
    phy=load(cdir/'physical.json');data=load(cdir/'input.json')
    n,k,P,r,s=[phy[a] for a in ('n','k','P','r','s')]
    need(P>0 and n==2*P*(s*s+r*r),'Invalid dimensions')
    hx,hz,x,z,xs,zs=([vec(v,n) for v in phy[key]] for key in ('HX','HZ','X_basis','Z_basis','X_seeds','Z_seeds'))
    bx,bz=basis(hx),basis(hz)
    need(n-len(bx)-len(bz)==k==2*P*(s-r)**2,'Dimension/rank mismatch')
    need(orthogonal(hx,hz),'Noncommuting CSS matrices')
    ax,az=reconstruct(data,P,r,s)
    need(equalspan(hx,ax) and equalspan(hz,az),'Physical code differs from constituent input')
    need(len(x)==len(z)==k,'Incomplete frame')
    need(orthogonal(x,hz) and orthogonal(z,hx),'Frame has nonzero syndrome')
    need(all((a&b).bit_count()%2==(i==j) for i,a in enumerate(x) for j,b in enumerate(z)),'Frame Gram is not identity')
    need(len(basis(hx+x))==len(bx)+k and len(basis(hz+z))==len(bz)+k,'Logical basis rank mismatch')
    need(x==[shifted(v,P,n,t) for v in xs for t in range(P)],'X basis is not the translated seed frame')
    need(z==[shifted(v,P,n,t) for v in zs for t in range(P)],'Z basis is not the translated seed frame')
    need(equalspan(hx,[shifted(v,P,n) for v in hx]) and equalspan(hz,[shifted(v,P,n) for v in hz]),'Missing common cyclic symmetry')
    iso=equalspan([transpose(v,P,r,s) for v in hx],hz) and equalspan([transpose(v,P,r,s) for v in hz],hx)
    exact=cdir/f'exact_{sector}.mat';nn,h,g,_=weighted(exact)
    h_expected,g_expected=(hx,hz) if sector=='Z' else (hz,hx)
    need(nn==n and equalspan(h,h_expected) and equalspan(g,g_expected),'Exact matrix represents a different CSS sector')
    upper={};supports={}
    for sec,hh,gg in [('Z',hx,hz),('X',hz,hx)]:
        nn,hs,gs,initial=weighted(cdir/f'sample_{sec}.mat',True)
        need(nn==n and equalspan(hs,hh) and equalspan(gs,gg),'Sampling/witness matrix mismatch')
        need(len(initial)>0,'Missing explicit upper-bound witness')
        bb=basis(gg)
        for v in initial:need(v and orthogonal([v],hh) and outside(v,bb),'Invalid physical logical witness')
        best=min(initial,key=int.bit_count);upper[sec]=best.bit_count();supports[sec]=[j for j in range(n) if(best>>j)&1]
    return {'n':n,'k':k,'P':P,'r':r,'s':s,'rank_HX':len(bx),'rank_HZ':len(bz),
            'full_Gram_identity':True,'common_cyclic_action':True,'constituent_reconstruction_matches':True,
            'sector_transpose_isometry':iso,'frame_width':max(v.bit_count() for v in xs+zs),
            'balanced_check_max':max(v.bit_count() for v in hx+hz),'verified_upper_bounds':upper,'witnesses':supports}

def log_paths(receipt_path,receipt):
    # Standard root protocol: PREFIX.receipt.json, PREFIX.stdout, PREFIX.stderr.
    name=receipt_path.name
    need(name.endswith('.receipt.json'),'Use the .receipt.json protocol')
    prefix=name[:-len('.receipt.json')]
    return receipt_path.with_name(prefix+'.stdout'),receipt_path.with_name(prefix+'.stderr')

def validate(args):
    cdir=Path(args.candidate_dir).resolve();rdir=Path(args.receipts_dir).resolve();edir=Path(args.engine_dir).resolve()
    source={name:sha(edir/name) for name in TRUSTED_SOURCE}
    need(source==TRUSTED_SOURCE,'Source is not the independently reviewed exact engine')
    exe=sha(args.executable);matrix=cdir/f'exact_{args.sector}.mat';matrix_hash=sha(matrix)
    audit=audit_physical(cdir,args.sector);P=audit['P'];roots=audit['n']//P
    paths=sorted(rdir.glob(args.pattern));need(len(paths)==args.shards,f'Expected {args.shards} receipts; found {len(paths)}')
    results={};receipts=[]
    for path in paths:
        rec=load(path);cmd=rec['command']
        need(isinstance(cmd,list) and len(cmd) in (10,11) and all(isinstance(v,str) for v in cmd),'Invalid command record')
        need(cmd[5]=='--shard' and int(cmd[2])==P and int(cmd[3])==args.weight,'Wrong command period/weight/partition')
        need(Path(cmd[1]).name==matrix.name,'Wrong command input sector')
        need(math.isfinite(float(cmd[4])) and float(cmd[4])>0 and int(cmd[9])>0,'Invalid command resource bounds')
        shard,count,split=map(int,cmd[6:9]);need(0<=shard<args.shards and count==args.shards and split==args.split_weight,'Wrong shard partition')
        need(len(cmd)==10 or cmd[10]=='no-coset','Unrecognized extra engine argument')
        need(shard not in results,'Duplicate shard receipt')
        need(rec['returncode']==0,'Nonzero exit status earns no exclusion credit')
        need(rec['input_sha256']==matrix_hash,'Receipt input hash mismatch')
        need(rec['executable_sha256']==exe,'Retained executable hash mismatch')
        src=rec['source_sha256']
        if isinstance(src,dict):src={Path(k).name:v for k,v in src.items()}
        else:src={'exact.cpp':src,'common.hpp':rec.get('header_sha256')}
        need(all(src.get(name)==value for name,value in source.items()),'Receipt source/header hash mismatch')
        stdout,stderr=log_paths(path,rec)
        need(sha(stdout)==rec['stdout_sha256'] and sha(stderr)==rec['stderr_sha256'],'Log hash mismatch')
        events=[json.loads(line) for line in stdout.read_text().splitlines() if line.strip()]
        need(len(events)==1 and events[0].get('event')=='result','Unexpected witness or log events')
        e=events[0]
        need(e['status']=='excluded','Incomplete/witness result earns no exclusion credit')
        need(e['W']==args.weight and e['partition']=='hash_prefix_v1','Result weight/partition mismatch')
        need(e['shard']==shard and e['shards']==count and e['split_weight']==split,'Result shard mismatch')
        need(e['root']==0 and e['branch_begin']==0 and e['branch_end']==64,'Unexpected root/branch restrictions')
        need(e['roots_completed']==roots,'Not all normalized roots completed')
        need(e['coset_pruning']==(len(cmd)==10),'Pruning flag differs from command')
        need(1<=e['nodes']<=int(cmd[9]) and e['seconds']>=0,'Invalid node/time statistics')
        need(0<=e['frontier_owned']<=e['frontier_seen']<=e['common_nodes'],'Invalid frontier statistics')
        need(0<=int(e['frontier_digest'])<2**64,'Invalid frontier digest')
        results[shard]=e;receipts.append({'path':str(path),'receipt_sha256':sha(path),'stdout_sha256':sha(stdout),'stderr_sha256':sha(stderr)})
    need(set(results)==set(range(args.shards)),'Shard union has missing IDs')
    first=results[0]
    for e in results.values():
        need(all(e[key]==first[key] for key in ('frontier_seen','frontier_digest','common_nodes')),'Common frontier differs across shards')
    need(sum(e['frontier_owned'] for e in results.values())==first['frontier_seen'],'Owned frontier counts do not cover the union exactly')
    lower=args.weight+1;upper=audit['verified_upper_bounds'][args.sector]
    need(upper>=lower,'Verified witness contradicts exclusion result')
    full=audit['sector_transpose_isometry']
    if full:need(min(audit['verified_upper_bounds'].values())>=lower,'Opposite-sector witness contradicts isometric lower bound')
    return {'status':'CERTIFIED_EXACT_DISTANCE' if full and upper==lower else 'CERTIFIED_INTERVAL' if full else 'CERTIFIED_ONE_SECTOR_LOWER_BOUND',
            'distance_lower_bound':lower if full else None,'sector_lower_bound':{args.sector:lower},
            'distance_upper_bound':min(audit['verified_upper_bounds'].values()),
            'exact_distance':lower if full and upper==lower else None,'algebra':audit,
            'partition':{'name':'hash_prefix_v1','shards':args.shards,'split_weight':args.split_weight,'roots':roots,
                         'frontier_seen':first['frontier_seen'],'frontier_digest':first['frontier_digest'],'owned_sum':sum(e['frontier_owned'] for e in results.values())},
            'total_nodes':sum(e['nodes'] for e in results.values()),'summed_seconds':sum(e['seconds'] for e in results.values()),
            'source_sha256':source,'executable_sha256':exe,'input_sha256':matrix_hash,'receipts':receipts,
            'validator_sha256':sha(__file__),'physical_artifact_sha256':{f:sha(cdir/f) for f in ['input.json','physical.json',f'exact_{args.sector}.mat','sample_X.mat','sample_Z.mat']}}

if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    for field in ('candidate-dir','receipts-dir','engine-dir','executable','output'):ap.add_argument('--'+field,required=True)
    ap.add_argument('--weight',type=int,required=True);ap.add_argument('--shards',type=int,default=64)
    ap.add_argument('--split-weight',type=int,default=5);ap.add_argument('--sector',choices=['X','Z'],default='Z')
    ap.add_argument('--pattern',default='**/*.receipt.json');args=ap.parse_args()
    try:
        need(args.weight>0 and args.shards>0 and 1<=args.split_weight<=args.weight,'Invalid target partition')
        result=validate(args);rc=0
    except Exception as exc:
        result={'status':'VALIDATION_FAILED_NO_LOWER_BOUND_CREDIT','error':str(exc),'validator_sha256':sha(__file__)};rc=2
    Path(args.output).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k in ('status','error','exact_distance','distance_lower_bound','distance_upper_bound','total_nodes','summed_seconds')}))
    sys.exit(rc)
