"""Portable physical validation and completed-exclusion record audit.
This quick check does not replay the deep search or verify a node-by-node proof.
"""

if not __debug__:
    raise RuntimeError(
        'Verification requires Python assertions; run without -O/-OO and unset PYTHONOPTIMIZE.'
    )

from pathlib import Path
import argparse,json,sys,subprocess,tempfile
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1];sys.path.insert(0,str(ROOT))
import rlf
from binary_audit import basis,remainder,read_matrix,sha,audit_css,check_witness,bits

def read(p):return json.loads(p.read_text())
def same_spaces(a,b):
 assert len(basis(a))==len(basis(b)) and all(not remainder(x,basis(a)) for x in b)
def transcript(p):return [json.loads(line) for line in p.read_text().splitlines()]
def check_receipt(candidate,path,n,P,H,G):
 """Read preserved native leaves recursively; never credit an interrupted root."""
 candidate=candidate.resolve();path=path.resolve();source_hash=sha(HERE/'src/exact.cpp');memo={};visiting=set()
 def locate(ref):
  # Original absolute paths are metadata; only preserved files in this case are read.
  p=Path(ref);q=candidate/p.parent.name/p.name
  assert q.is_file() and q.resolve().is_relative_to(candidate)
  return q.resolve()
 def common(p):
  rec=read(p);d=p.parent;f=rec['final'];assert f['W']==W and rec['sector']=='Z'
  assert rec['source_sha256']==source_hash and rec['engine_sha256']==engine_hash
  if 'header_sha256' in rec:assert rec['header_sha256']==sha(HERE/'src/common.hpp')
  assert rec['original_matrix_sha256']==sha(candidate/'code.mat')
  assert rec['execution_input_sha256']==input_hash==sha(d/'input.mat')
  nn,A,B=read_matrix(d/'input.mat');assert nn==n;same_spaces(H,A);same_spaces(G,B);audit_css(n,A,B,P)
  return rec,f
 def leaf(jp):
  jr=read(jp);log=jp.with_name(jp.stem+'.stdout.jsonl');events=transcript(log)
  assert jr['stdout_sha256']==sha(log) and events[-1]==jr['final']
  f=jr['final'];r,b,e=f['root'],f['branch_begin'],f['branch_end']
  assert jp.stem==f'root{r:02d}_b{b:02d}_{e:02d}' and f['event']=='result' and f['W']==W and f['partition']=='root_branch_v1'
  assert 0<=r<n//P and 0<=b<e<=64 and f['coset_pruning']
  cmd=jr['command'];assert Path(cmd[1]).name=='input.mat' and list(map(int,cmd[2:4]))==[P,W] and list(map(int,cmd[5:8]))==[r,b,e] and len(cmd)==9
  assert sha(jp.parent/'input.mat')==input_hash
  if f['status']=='incomplete':
   assert jr['returncode']==75 and not any(x['event']=='witness' for x in events);return None
  assert f['status']=='excluded' and jr['returncode']==0 and f['roots_completed']==1 and len(events)==1
  return (r,b,e,f['nodes'])
 def collect(p):
  if p in memo:return memo[p]
  assert p not in visiting;visiting.add(p);rec,f=common(p);d=p.parent;leaves={}
  if rec['mode'] not in ('exact_partition_completion','exact_native_interval_union'):
   assert rec['mode'] in ('exact','external_exact_independently_audited')
   log=d/'stdout.jsonl';assert rec['stdout_sha256']==sha(log);events=transcript(log)
   assert len(events)==1 and events[0]==f and f['event']=='result' and f['partition']=='hash_prefix_v1' and f['shards']==1 and f['shard']==0 and f['frontier_owned']==f['frontier_seen']
   cmd=rec['command'];assert list(map(int,cmd[2:4]))==[P,W] and cmd[5:9]==['--shard','0','1','4']
   count=f['roots_completed'];assert 0<=count<=n//P
   if f['status']=='excluded':assert rec['returncode']==0 and count==n//P
   else:assert f['status']=='incomplete' and rec['returncode']==75 and count<n//P
   # A partial all-root run supplies completed-prefix coverage, but its aggregate
   # node count includes the interrupted root and is not credited as completed work.
   for r in range(count):leaves[(str(p),r)]=(r,0,64,f['nodes'] if f['status']=='excluded' and r==0 else 0)
  else:
   cov=read(d/'coverage.json');assert rec['coverage_sha256']==sha(d/'coverage.json') and cov['input_sha256']==input_hash and cov['source_sha256']==source_hash and cov['engine_sha256']==engine_hash
   if rec['mode']=='exact_partition_completion':
    prefix=locate(cov['prefix_receipt']);assert sha(prefix)==cov['prefix_receipt_sha256'];pr=read(prefix);assert pr['final']['status']=='incomplete'
    leaves.update(collect(prefix))
    for j in cov['jobs']:
     jp=d/f"root{j['root']:02d}_b{j['branch_begin']:02d}_{j['branch_end']:02d}.json";assert read(jp)==j
     v=leaf(jp)
     if v is not None:assert str(jp) not in leaves;leaves[str(jp)]=v
   else:
    for j in cov['records']:
     jp=locate(j['native_receipt']);log=locate(j['stdout']);assert sha(jp)==j['receipt_sha256'] and sha(log)==j['stdout_sha256'] and log==jp.with_name(jp.stem+'.stdout.jsonl')
     v=leaf(jp);assert v is not None and v[:3]==(j['root'],j['branch_begin'],j['branch_end']) and v[3]==j['nodes']
     if jp.parent!=d:
      previous=collect(jp.parent/'receipt.json');assert previous.get(str(jp))==v
     assert str(jp) not in leaves;leaves[str(jp)]=v
   assert cov['complete']==(f['status']=='excluded')
  coverage={r:[] for r in range(n//P)}
  for r,b,e,nodes in leaves.values():coverage[r].append((b,e))
  full=[]
  for r,intervals in coverage.items():
   last=0;contiguous=True
   for b,e in sorted(intervals):
    assert b>=last and b<e<=64
    if b!=last:contiguous=False
    last=e
   if contiguous and last==64:full.append(r)
  assert f['roots_completed']==len(full)
  if f['status']=='excluded':assert full==list(range(n//P)) and rec['lower_bound_this_sector']==W+1
  else:assert rec.get('lower_bound_this_sector') is None
  if 'covered_roots' in f:assert f['covered_roots']==full
  visiting.remove(p);memo[p]=leaves;return leaves
 rec=read(path);W=rec['final']['W'];engine_hash=rec['engine_sha256'];input_hash=rec['execution_input_sha256'];assert rec['final']['status']=='excluded' and rec['lower_bound_this_sector']==W+1
 return sum(v[3] for v in collect(path).values())

def audit(entry):
 code=read(ROOT/entry['path']);directory=(ROOT/entry['path']).parent;checks=read(directory/'checks.json');frame=read(directory/'frames/paper.json');report=rlf.verify(code,checks,frame)
 cp=next(ROOT/p for p in code['evidence'] if p.endswith('/certificate.json'));c=cp.parent;cert=read(cp);assert cert['status']=='CERTIFIED_EXACT' and cert['code_json_sha256']==sha(c/'code.json') and cert['audit_sha256']==sha(c/'audit.json') and cert['matrix_sha256']==sha(c/'code.mat')
 n,H,G=read_matrix(c/'code.mat');P=code['group']['order'];audit_css(n,H,G,P);same_spaces(H,[bits(x,n) for x in checks['X']]);same_spaces(G,[bits(x,n) for x in checks['Z']]);assert code['distance']['lower']==code['distance']['upper']==cert['distance_exact'];check_witness(bits(cert['best_witness']['support'],n),n,H,G,cert['distance_exact'])
 proven=[]
 for item in cert['receipts']:
  if item['lower'] is None:continue
  p=c/item['file'];assert item['sha256']==sha(p);rec=read(p);assert item['sector']=='Z' and item['lower']==rec['lower_bound_this_sector'];proven.append((item['lower'],check_receipt(c,p,n,P,H,G)))
 assert proven and max(x[0] for x in proven)==code['distance']['lower'];report.update(completed_exclusion_records_checked=True,nodes_in_successful_branch_records=sum(x[1] for x in proven),deep_search_replayed=False);return report

if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--code');p.add_argument('--report',type=Path);p.add_argument('--replay',action='store_true',help='Compile and rerun the full exclusion for the selected code.');p.add_argument('--seconds',type=float,default=60);a=p.parse_args()
 entries=[]
 for e in read(ROOT/'catalog.json')['entries']:
  c=read(ROOT/e['path'])
  if any('small_pp_records/small_quaternary_pp/' in s for s in c.get('evidence',[])) and (a.code is None or e['id']==a.code):entries.append(e)
 assert entries,'No matching campaign code';reports=[audit(e) for e in entries]
 if a.replay:
  assert a.code and len(entries)==1,'Select one code for replay';code=read(ROOT/entries[0]['path']);c=(ROOT/next(s for s in code['evidence'] if s.endswith('/certificate.json'))).parent
  with tempfile.TemporaryDirectory(prefix='pp-small-replay-') as tmp:
   exe=Path(tmp)/'exact';subprocess.run(['c++','-O3','-std=c++17',str(HERE/'src/exact.cpp'),'-o',str(exe)],check=True)
   result=subprocess.run([str(exe),str(c/'redundant3.mat'),str(code['group']['order']),str(code['distance']['lower']-1),str(a.seconds),'--shard','0','1','4','18446744073709551615'],capture_output=True,text=True);events=[json.loads(x) for x in result.stdout.splitlines()];f=events[-1];reports[0]['replay_result']=f;reports[0]['deep_search_replayed']=result.returncode==0 and f['status']=='excluded' and f['roots_completed']==code['n']//code['group']['order'];assert result.returncode in [0,75], 'Replay failed or found a counterexample: '+result.stdout+result.stderr
 if a.report:a.report.write_text(json.dumps(reports,indent=2)+'\n')
 print(json.dumps(reports,indent=2))
 if a.replay and not reports[0]['deep_search_replayed']:raise SystemExit(75)
