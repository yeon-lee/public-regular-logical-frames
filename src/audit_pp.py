"""Independent F4 cyclic reversing-PP reconstruction and physical export.
No imports from construction/search code. Works at any period when the chosen
pivot is explicitly invertible; augmentation alone is never a unit proof.
"""

if not __debug__:
    raise RuntimeError(
        'Verification requires Python assertions; run without -O/-OO and unset PYTHONOPTIMIZE.'
    )

from pathlib import Path
import argparse,json,itertools,collections,functools,sys,time
HERE=Path(__file__).resolve().parent;sys.path.insert(0,str(HERE/'vendor'))
from binary_audit import audit_css,basis,bits,check_witness,dot,read_matrix,remainder,rotate,sha,support,write_matrix

def gf(a,b):
 c=0
 for _ in range(2):
  if b&1:c^=a
  b>>=1;a<<=1
  if a&4:a^=7
 return c

def coords(c):return [((c>>1)^c)&1,c&1]

class Ring:
 def __init__(self,P):self.P=P
 def add(self,*ps):
  v={}
  for p in ps:
   for e,c in p.items():v[e]=v.get(e,0)^c
  return {e:c for e,c in v.items() if c}
 def mul(self,a,b):
  v={}
  for e,c in a.items():
   for f,d in b.items():
    j=(e+f)%self.P;v[j]=v.get(j,0)^gf(c,d)
  return {e:c for e,c in v.items() if c}
 def star(self,p):return {(-e)%self.P:c for e,c in p.items()}
 def inverse(self,p):
  if len(p)==1:
   e,c=next(iter(p.items()));v={(-e)%self.P:gf(c,c)};assert self.mul(p,v)=={0:1};return v
  # Exact multiplication matrix over F4, valid at odd and composite periods.
  P=self.P;A=[[p.get((i-j)%P,0) for j in range(P)]+[int(i==0)] for i in range(P)]
  for j in range(P):
   k=next((k for k in range(j,P) if A[k][j]),None)
   if k is None:raise ValueError('Nonunit cyclic pivot pairing')
   A[j],A[k]=A[k],A[j];iv=gf(A[j][j],A[j][j]);A[j]=[gf(iv,c) for c in A[j]]
   for i in range(P):
    if i!=j and A[i][j]:
     c=A[i][j];A[i]=[x^gf(c,y) for x,y in zip(A[i],A[j])]
  v={i:A[i][-1] for i in range(P) if A[i][-1]};assert self.mul(p,v)=={0:1};return v

def reconstruct(source,out,with_cofactors=True,provided_frame=False,balanced=False):
 start=time.monotonic();source=Path(source);out=Path(out);out.mkdir(parents=True,exist_ok=False);source_bytes=source.read_bytes();(out/'source.json').write_bytes(source_bytes);raw=json.loads(source_bytes);P=raw['P'];C=raw.get('C',raw.get('coefficients'));assert C is not None
 flat_width=raw.get('L',8)
 if isinstance(C[0],int):
  assert len(C)%flat_width==0;C=[C[i*flat_width:(i+1)*flat_width] for i in range(len(C)//flat_width)]
 r,L=len(C),len(C[0]);assert L%2==0 and all(len(row)==L for row in C);M=L//2;assert r<M
 E=raw['E'];E=[E[i*L:(i+1)*L] for i in range(r)] if isinstance(E[0],int) else E;assert len(E)==r and all(len(row)==L for row in E)
 assert all(c in range(4) for row in C for c in row);E=[[e%P for e in row] for row in E];R=Ring(P);n=2*L*P;A=[[{E[i][l]:C[i][l]} if C[i][l] else {} for l in range(L)] for i in range(r)]
 def expand(reverse):
  H=[]
  for i in range(r):
   for a in range(2):
    for t in range(P):
     v=0
     for l in range(L):
      j=(l+M)%L if reverse else l;e=(-E[i][j] if reverse else E[i][j])%P
      for b in range(2):
       if coords(gf(C[i][j],[2,3][b]))[a]:v^=1<<((2*l+b)*P+(t-e)%P)
     H.append(v)
  return H
 HX,HZ=expand(False),expand(True)
 rawHX,rawHZ=HX[:],HZ[:]
 if balanced:
  balanced_rows=[]
  for H in [HX,HZ]:
   HH=[]
   for i in range(r):
    rows=[H[2*i*P:((2*i)+1)*P],H[((2*i)+1)*P:(2*i+2)*P]];rows.append([a^b for a,b in zip(*rows)])
    picks=sorted(range(3),key=lambda j:(rows[j][0].bit_count(),j))[:2]
    for j in picks:HH.extend(rows[j])
   assert len(basis(H))==len(basis(HH)) and all(not remainder(v,basis(H)) for v in HH);balanced_rows.append(HH)
  HX,HZ=balanced_rows
  write_matrix(out/'raw.mat',n,rawHX,rawHZ)
 def reflect(v):
  z=0
  for bit in support(v):
   l,t=divmod(bit,P);col,b=divmod(l,2);z^=1<<((2*((col+M)%L)+b)*P+(-t)%P)
  return z
 assert set(map(reflect,HX))==set(HZ)
 @functools.lru_cache(None)
 def minor(cols):
  v={}
  for perm in itertools.permutations(cols):
   term={0:1}
   for i,j in enumerate(perm):term=R.mul(term,A[i][j])
   v=R.add(v,term)
  return v
 pv=raw.get('pivot_columns',raw.get('pivotA',list(range(r))));assert len(pv)==r and len(set(pv))==r and all(0<=x<L for x in pv) and not set(pv)&{(x+M)%L for x in pv}
 free=[l for l in range(L) if l not in pv and (l+M)%L not in pv];assert len(free)==L-2*r
 delta=minor(tuple(pv));g=R.mul(delta,delta);gi=R.inverse(g);assert R.mul(g,gi)=={0:1}
 def binary(v,c):
  word=0
  for l,p in enumerate(v):
   for e,x in p.items():
    for b,z in enumerate(coords(gf(c,x))):
     if z:word^=1<<((2*l+b)*P+e)
  return word
 pol=[]
 for f in free:
  z=[{} for _ in range(L)];z[f]=delta
  for j,l in enumerate(pv):cols=pv[:];cols[j]=f;z[l]=minor(tuple(cols))
  pol.append(z)
 X=[];Z=[]
 for i,f in enumerate(free):
  partner=free.index((f+M)%L)
  for c in [2,3]:X.append(reflect(binary(pol[partner],c)));Z.append(binary([R.mul(p,gi) for p in pol[i]],c))
 if provided_frame:
  X=[bits(x,n) for x in raw.get('X_seed_supports',raw.get('X_seeds'))];Z=[bits(z,n) for z in raw.get('Z_seed_supports',raw.get('Z_seeds'))]
 audit=audit_css(n,HX,HZ,P,X,Z);assert audit['ranks']==[2*r*P]*2;assert audit['k']==2*(L-2*r)*P
 Xa=[rotate(x,n,P,t) for x in X for t in range(P)];Za=[rotate(z,n,P,t) for z in Z for t in range(P)];assert all(dot(x,z)==int(i==j) for i,x in enumerate(Xa) for j,z in enumerate(Za));assert len(basis(HX+Xa))==len(basis(HZ+Za))==(n+audit['k'])//2
 if 'X_seed_supports' in raw:assert X==[bits(x,n) for x in raw['X_seed_supports']]
 if 'Z_seed_supports' in raw:assert Z==[bits(z,n) for z in raw['Z_seed_supports']]
 for name,H in [('X',HX),('Z',HZ)]:
  col=[sum((h>>j)&1 for h in H) for j in range(n)];audit[name+'_check_max']=max(h.bit_count() for h in H);audit[name+'_check_mean']=sum(h.bit_count() for h in H)/len(H);audit[name+'_column_histogram']=dict(collections.Counter(col));audit[name+'_column_max']=max(col)
 par=list(range(n))
 def find(x):
  while par[x]!=x:par[x]=par[par[x]];x=par[x]
  return x
 for h in HX+HZ:
  su=support(h)
  for j in su[1:]:par[find(j)]=find(su[0])
 audit['physical_components']=sorted(collections.Counter(find(j) for j in range(n)).values())
 witnesses=[]
 for key,H,G in [('witness_Z',HX,HZ),('witness_X',HZ,HX)]:
  if key in raw:witnesses.append(dict(sector=key[-1],**check_witness(bits(raw[key],n),n,H,G)))
 seeds=list(Z);bg=basis(HZ)
 if with_cofactors:
  for cols in itertools.combinations(range(L),r+1):
   v=[{} for _ in range(L)]
   for l in cols:v[l]=minor(tuple(j for j in cols if j!=l))
   for c in [1,2,3]:
    z=binary(v,c)
    if z and remainder(z,bg):
     assert all(not dot(z,h) for h in HX);seeds.append(z)
 if 'witness_Z' in raw:seeds.append(bits(raw['witness_Z'],n))
 seeds=sorted(set(seeds),key=lambda z:(z.bit_count(),z));best=check_witness(seeds[0],n,HX,HZ)
 R3=[]
 for H in [HX,HZ]:
  H3=[]
  for i in range(r):
   for t in range(P):
    a,b=H[(2*i)*P+t],H[(2*i+1)*P+t];H3.extend([a,b,a^b])
  assert len(basis(H3))==len(basis(H)) and all(not remainder(v,basis(H)) for v in H3);R3.append(H3)
 write_matrix(out/'code.mat',n,HX,HZ);write_matrix(out/'redundant3.mat',n,*R3);write_matrix(out/'sample_Z.mat',n,HX,HZ,seeds);write_matrix(out/'sample_X.mat',n,HZ,HX,[reflect(z) for z in seeds])
 rec=dict(P=P,n=n,k=audit['k'],C=C,E=E,pivot_columns=pv,free_columns=free,delta=delta,pairing=g,pairing_inverse=gi,seed_weights_X=list(map(int.bit_count,X)),seed_weights_Z=list(map(int.bit_count,Z)),X_seed_supports=list(map(support,X)),Z_seed_supports=list(map(support,Z)),displayed_width=audit['displayed_width'],witness_Z=best['support'],d_upper=best['weight'],source=str(source.resolve()),source_sha256=sha(out/'source.json'),matrix_sha256=sha(out/'code.mat'),frame_source='supplied physical supports independently verified' if provided_frame else 'independently reconstructed cofactor frame',distance_status='Verified upper bound only; no lower bound is inferred by algebra.')
 (out/'code.json').write_text(json.dumps(rec,indent=2)+'\n');audit.update(balanced_presentation=balanced,independent_GF4_reconstruction=True,reflection_isometry_verified=True,pivot_pairing_explicit_unit_verified=True,full_translated_pairing_entries=len(Xa)**2,augmented_ranks=[len(basis(HX+Xa)),len(basis(HZ+Za))],Cramer_words_retained=len(seeds),source_witnesses=witnesses,cofactor_witness=best,redundant3_rowspaces_verified=True,source_sha256=sha(out/'source.json'),matrix_sha256=sha(out/'code.mat'),elapsed_seconds=time.monotonic()-start)
 (out/'audit.json').write_text(json.dumps(audit,indent=2)+'\n');return rec,audit

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('source');p.add_argument('out');p.add_argument('--no-cofactors',action='store_true');p.add_argument('--provided-frame',action='store_true');p.add_argument('--balanced',action='store_true');a=p.parse_args();rec,audit=reconstruct(a.source,a.out,not a.no_cofactors,a.provided_frame,a.balanced);print(json.dumps({k:rec[k] for k in ['n','k','P','displayed_width','d_upper','seed_weights_X','seed_weights_Z']}));print(json.dumps({k:audit[k] for k in ['ranks','X_check_max','X_column_max','physical_components','Cramer_words_retained','elapsed_seconds']}))
