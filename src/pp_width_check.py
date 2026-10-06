#!/usr/bin/env python3
"""Independent algebra/receipt audit; optionally replay only the split frontier.

This does not replay the 39-billion-node exhaustive exclusion. Completion
receipts are evidence of execution, not a formal proof object for every node.
"""
from pathlib import Path
import argparse, hashlib, json, subprocess, tempfile

HERE = Path(__file__).resolve().parent
H = HERE / 'records'
def need(x, message):
    if not x: raise ValueError(message)

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

def sha(path):
    """File digest or explicitly declared omitted executable identity."""
    return evidence_digest(path)
def load(path): return json.loads(Path(path).read_text())
def bits(support, n):
    need(len(set(support)) == len(support) and all(0 <= q < n for q in support), 'Invalid support')
    return sum(1 << q for q in support)
def basis(rows):
    out = {}
    for v in rows:
        while v:
            q = v.bit_length()-1
            if q not in out: out[q] = v; break
            v ^= out[q]
    return out
def rem(v, b):
    while v:
        q = v.bit_length()-1
        if q not in b: break
        v ^= b[q]
    return v
def matrix(path):
    lines = Path(path).read_text().splitlines(); n, m, g = map(int, lines[0].split()); rows=[]
    for line in lines[1:]:
        s = list(map(int, line.split())); need(s[0] == len(s)-1, 'Invalid matrix row'); rows.append(bits(s[1:], n))
    need(len(rows) == m+g, 'Wrong matrix row count')
    return n, rows[:m], rows[m:]
def same_span(a,b):
    aa,bb=basis(a),basis(b)
    return len(aa)==len(bb) and all(rem(v,bb)==0 for v in a)
def shift(v,t,P,n):
    t%=P; mask=(1<<P)-1; out=0
    for b in range(n//P):
        x=(v>>(P*b))&mask; out|=(((x<<t)|(x>>(P-t)))&mask)<<(P*b)
    return out

def audit():
    manifest=load(H/'bundle_sha256.json')['files']
    for rel,want in manifest.items(): need(sha(H/rel)==want, 'Public input/source hash: '+rel)
    build=load(H/'build.json'); need(build['source_sha256']==sha(H/'exact_target.cpp'), 'Build source mismatch')
    need(build['bundle_sha256']==sha(H/'bundle_sha256.json'), 'Build manifest mismatch')
    for rel,want in load(H/'orchestration_manifest.json').items(): need(sha(H/rel)==want, 'Orchestration mismatch')
    for i in range(4):
        g=load(H/f'results/group_{i:02d}.json')
        need(g['build']==build and g['returncode']==0 and g['group']==i, 'Group failed')
        need(g['shards']==list(range(16*i,16*(i+1))), 'Group partition mismatch')
        need(g['argv'][2:8]==['--start',str(16*i),'--count','16','--shards','64'], 'Group command mismatch')
    n,HX,HZ=matrix(H/'data/code.mat'); sn,SX,SZ=matrix(H/'data/search.mat')
    need(n==sn==1024 and len(basis(HX))==len(basis(HZ))==384, 'Code dimensions/ranks')
    code=load(H/'data/code.json');coeff=code['coefficients'];exponents=code['E']
    rho=[((0,0),(0,0)),((1,0),(0,1)),((0,1),(1,1)),((1,1),(1,0))]
    expanded=[]
    for i in range(3):
        for component in range(2):
            for t in range(64):
                expanded.append(bits([(2*ell+b)*64+(t-exponents[8*i+ell])%64
                    for ell in range(8) for b in range(2) if rho[coeff[i][ell]][component][b]],n))
    need(sorted(HX)==sorted(expanded), 'Coefficient/exponent expansion mismatch')
    reflected=[bits([((q//64+8)%16)*64+(-q%64) for q in range(n) if (v>>q)&1],n) for v in HX]
    need(sorted(reflected)==sorted(HZ), 'Reflection does not produce Z checks')
    need(all((x&z).bit_count()%2==0 for x in HX for z in HZ), 'Noncommuting checks')
    need(same_span(HX,SX) and same_span(HZ,SZ), 'Search matrix changes code')
    # The original frozen matrix uses both trace components. The manuscript's
    # weight-eleven presentation selects the two lighter independent rows from
    # {a,b,a+b} at each quaternary row and translate, without changing the span.
    for raw in [HX,HZ]:
        balanced=[]
        for i in range(3):
            for t in range(64):
                a,b=raw[128*i+t],raw[128*i+64+t]
                balanced.extend(sorted([a,b,a^b],key=int.bit_count)[:2])
        need(same_span(raw,balanced) and max(v.bit_count() for v in balanced)==11, 'Balanced check presentation')
    need(same_span(HX,[shift(v,1,64,n) for v in HX]) and same_span(HZ,[shift(v,1,64,n) for v in HZ]), 'Translation not symmetry')
    d=load(H/'data/detector.json'); need(d['P']==64 and d['W']==24 and d['normalization_roots']==[6,15], 'Detector roots')
    mask=d['detector']['mask']; need(mask==(1<<6)|(1<<15), 'Detector mask')
    functional=sum(1<<j for j in range(n) if (mask>>(j//64))&1)
    need(all((g&functional).bit_count()%2==0 for g in HZ), 'Detector nonzero on stabilizer')
    need(rem(functional,basis(HX))!=0, 'Detector vanishes on entire kernel')
    frame=load(H/'frame_upper_bound.json'); X=[bits(s,n) for s in frame['X_seed_supports']];Z=[bits(s,n) for s in frame['Z_seed_supports']]
    XX=[shift(v,t,64,n) for v in X for t in range(64)]; ZZ=[shift(v,t,64,n) for v in Z for t in range(64)]
    need(len(XX)==len(ZZ)==256 and max(v.bit_count() for v in XX+ZZ)==25, 'Frame size/width')
    need(all((x&h).bit_count()%2==0 for x in XX for h in HZ), 'X frame syndrome')
    need(all((z&h).bit_count()%2==0 for z in ZZ for h in HX), 'Z frame syndrome')
    need(all((x&z).bit_count()%2==(i==j) for i,x in enumerate(XX) for j,z in enumerate(ZZ)), 'Canonical Gram mismatch')
    v=bits(d['attaining_target_support'],n)
    need(v.bit_count()==25 and (v&functional).bit_count()%2==1 and all((v&h).bit_count()%2==0 for h in HX), 'Attaining target witness')
    need(any((z&functional).bit_count()%2 for z in Z), 'Detector trivial on frame')
    for support in d['known_short_supports']:
        v=bits(support,n);need(v.bit_count()==24 and all((v&h).bit_count()%2==0 for h in HX) and rem(v,basis(HZ))!=0, 'Distance upper witness')
    results=[]
    for i in range(64):
        p=H/f'results/target_W24_shard_{i:05d}.json'; r=load(p)
        need(r['source_and_input_sha256']==manifest and r['returncode']==0 and r['detector_mask']==mask, 'Shard identity')
        need(sha(p.with_suffix('.jsonl'))==r['log_sha256'], 'Shard log hash')
        events=[json.loads(s) for s in p.with_suffix('.jsonl').read_text().splitlines() if s.startswith('{')]
        need(len(events)==1 and events[0]==r['result'], 'Unexpected events')
        z=r['result'];expect={'status':'excluded','W':24,'roots_completed':2,'detector_mask':mask,'partition':'hash_prefix_v1','shard':i,'shards':64,'split_weight':4,'root':0,'branch_begin':0,'branch_end':64,'coset_pruning':True}
        for k,v in expect.items():need(z[k]==v,'Shard setting '+k)
        a=r['argv'];need(len(a)==12 and Path(a[0]).name=='exact_target' and Path(a[1]).name=='search.mat' and a[2:4]==['64','24'] and a[5:9]==['--shard',str(i),'64','4'] and a[-2:]==['--detector','0x8040'], 'Shard command')
        results.append(z)
    for key in ['frontier_seen','frontier_digest','common_nodes','split_weight']:
        need(len({r[key] for r in results})==1,'Different shard frontiers')
    need(sum(r['frontier_owned'] for r in results)==results[0]['frontier_seen']==856,'Incomplete frontier partition')
    return results

def replay_prefix(results, compiler):
    """Replay all prefix decisions only; cannot certify omitted deep subtrees."""
    s=(H/'exact_target.cpp').read_text(); marker='        if(bad){++coset_prunes;return;}'
    need(s.count(marker)==1,'Unrecognized source for prefix replay')
    s=s.replace(marker,'        if(hashed&&wt==split_weight)return; // PREFIX ONLY: no exclusion claim\n'+marker)
    old='incomplete?"incomplete":"excluded"';need(s.count(old)==1,'Unrecognized result tag')
    s=s.replace(old,'incomplete?"incomplete":"prefix_complete"')
    with tempfile.TemporaryDirectory() as tmp:
        tmp=Path(tmp);(tmp/'prefix.cpp').write_text(s);(tmp/'common.hpp').write_bytes((H/'common.hpp').read_bytes())
        subprocess.run([compiler,'-O3','-std=c++17',str(tmp/'prefix.cpp'),'-o',str(tmp/'prefix')],check=True)
        for i,want in enumerate(results):
            cmd=[str(tmp/'prefix'),str(H/'data/search.mat'),'64','24','60','--shard',str(i),'64','4',str(2**64-1),'--detector','0x8040']
            got=json.loads(subprocess.run(cmd,check=True,capture_output=True,text=True).stdout)
            need(got['status']=='prefix_complete' and got['roots_completed']==2,'Prefix did not complete')
            for k in ['frontier_seen','frontier_owned','frontier_digest','common_nodes']:need(got[k]==want[k],'Prefix mismatch '+k)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--replay-prefix',action='store_true');ap.add_argument('--compiler',default='c++');ap.add_argument('--output');args=ap.parse_args()
    results=audit()
    if args.replay_prefix:replay_prefix(results,args.compiler)
    report={'status':'passed','independent_algebraic_checks':True,'completed_execution_receipts':64,'required_detector_roots':[6,15],'weight_excluded_through':24,'attaining_detector_weight':25,'full_canonical_frame_width':25,'reported_search_nodes':sum(r['nodes'] for r in results),'reported_sum_process_seconds':sum(r['seconds'] for r in results),'reported_shard_seconds_range':[min(r['seconds'] for r in results),max(r['seconds'] for r in results)],'prefix_replayed':args.replay_prefix,'deep_exclusion_replayed':False,'formal_proof_object':False,'execution_compiler':'g++ -O3 -std=c++17','execution_compiler_version':None,'execution_processor_model':None}
    text=json.dumps(report,indent=2)+'\n'
    if args.output:Path(args.output).write_text(text)
    print(text,end='')
