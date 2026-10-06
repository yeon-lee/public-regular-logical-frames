#!/usr/bin/env python3
"""Reconstruct both check8 PP codes and optionally replay exhaustive exclusions."""

if not __debug__:
    raise RuntimeError(
        'Verification requires Python assertions; run without -O/-OO and unset PYTHONOPTIMIZE.'
    )

from pathlib import Path
import argparse,json,sys,tempfile,subprocess,hashlib
P=Path(__file__).resolve().parent;sys.path.insert(0,str(P/'src'))
from audit_pp import reconstruct,gf
from binary_audit import read_matrix,check_witness,bits,sha,basis,remainder,rotate,support
p=argparse.ArgumentParser();p.add_argument('--replay',action='store_true',help='Compile the supplied engine and freshly replay both full exclusions');a=p.parse_args()
manifest=json.loads((P/'SHA256SUMS.json').read_text())
for name,h in manifest['files'].items():assert sha(P/name)==h,('integrity mismatch',name)
results=[]
with tempfile.TemporaryDirectory(prefix='check8_pp_')as td:
 T=Path(td)
 if a.replay:subprocess.run(['c++','-O3','-std=c++17',str(P/'src/exact.cpp'),'-o',str(T/'exact')],check=True)
 for name in ['pp512_128_11_check8','pp1024_256_12_check8']:
  o=P/name;c=json.loads((o/'code.json').read_text());cert=json.loads((o/'certificate.json').read_text());P0=c['P'];n=c['n'];d=cert['d']
  for f,h in cert['hashes'].items():assert sha(P/f)==h,('certificate binding',f)
  rec,au=reconstruct(o/'code.json',T/name,with_cofactors=False)
  assert (T/name/'code.mat').read_bytes()==(o/'code.mat').read_bytes()
  assert (T/name/'redundant3.mat').read_bytes()==(o/'redundant3.mat').read_bytes()
  n,HX,HZ=read_matrix(o/'code.mat');_,HX3,HZ3=read_matrix(o/'redundant3.mat')
  for H,H3 in [(HX,HX3),(HZ,HZ3)]:
   assert len(basis(H))==len(basis(H3)) and all(not remainder(v,basis(H))for v in H3)
   assert all(not remainder(rotate(v,n,P0,1),basis(H))for v in H)
  def reflect(v):
   out=0
   for bit in support(v):
    b,t=divmod(bit,P0);col,trace=divmod(b,2);out^=1<<((2*((col+4)%8)+trace)*P0+(-t)%P0)
   return out
  assert set(map(reflect,HX))==set(HZ) and set(map(reflect,HZ))==set(HX)
  w=check_witness(bits(c['witness_Z'],n),n,HX,HZ,d)
  check_witness(reflect(bits(c['witness_Z'],n)),n,HZ,HX,d)
  assert max(map(int.bit_count,HX+HZ))==8
  assert au['physical_components']==[n] and au['k']==c['k']
  cy=c['nontrivial_coefficient_cycle'];i,j=cy['rows'];u,v=cy['columns'];C=c['C'];assert gf(C[i][u],C[j][v])!=gf(C[i][v],C[j][u])
  events=[json.loads(l)for l in(o/'exact.stdout').read_text().splitlines()if l.startswith('{')];e=events[-1]
  assert e==cert['exact_result'];assert e['status']=='excluded'and e['W']==d-1 and e['roots_completed']==n//P0==16
  assert e['partition']=='hash_prefix_v1'and e['shards']==1 and e['shard']==0 and e['split_weight']==1
  assert e['frontier_seen']==e['frontier_owned']==16
  r=dict(name=name,n=n,k=c['k'],d=d,frame_width=au['displayed_width'],check_max=8,column_max=max(au['X_column_max'],au['Z_column_max']),all_checks_passed=True)
  if a.replay:
   cmd=[str(T/'exact'),str(o/'redundant3.mat'),str(P0),str(d-1),'60','--shard','0','1','1','10000000000'];q=subprocess.run(cmd,capture_output=True,text=True)
   assert q.returncode==0,(q.returncode,q.stdout,q.stderr)
   fresh=[json.loads(l)for l in q.stdout.splitlines()if l.startswith('{')][-1]
   assert fresh['status']=='excluded'and fresh['W']==d-1 and fresh['roots_completed']==16
   r['fresh_exclusion']=fresh
  results.append(r)
print(json.dumps(dict(status='PASS',replay=a.replay,codes=results),indent=2))
