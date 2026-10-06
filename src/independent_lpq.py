"""Cyclic F4 lifted products with two independent constituents.

Conventions copied verbatim from the independently audited cyclic LP matrices:
bar reverses group exponents and DOES NOT conjugate field coefficients.
Only Python's standard library is required. Distance values are upper bounds.
"""

if not __debug__:
    raise RuntimeError(
        'Verification requires Python assertions; run without -O/-OO and unset PYTHONOPTIMIZE.'
    )

import itertools
from collections import Counter

MUL=((0,0,0,0),(0,1,2,3),(0,2,3,1),(0,3,1,2))
INV=(0,1,3,2)
COORD=(0,3,1,2)
RHO={1:((1,0),(0,1)),2:((0,1),(1,1)),3:((1,1),(1,0))}

def add(a,b):
    c=dict(a)
    for e,v in b.items():
        w=c.get(e,0)^v
        if w:c[e]=w
        else:c.pop(e,None)
    return c

def mul(a,b,P):
    c={}
    for e,v in a.items():
        for f,w in b.items():
            i=(e+f)%P; z=c.get(i,0)^MUL[v][w]
            if z:c[i]=z
            else:c.pop(i,None)
    return c

def scale(a,c):return {e:MUL[v][c] for e,v in a.items()} if c else {}
def bar(a,P):return {(-e)%P:v for e,v in a.items()}
def weight(a,c=1):return sum(COORD[MUL[v][c]].bit_count() for v in a.values())
def inverse(a,P):
    """Power-of-two cyclic ring inverse from q=1+a/aug(a), q^P=0."""
    assert P>0 and P&(P-1)==0
    aug=0
    for c in a.values():aug^=c
    if not aug:return None
    q=add(scale(a,INV[aug]),{0:1}); out={0:INV[aug]}; t=1
    while t<P:
        out=mul(out,add({0:1},q),P);q=mul(q,q,P);t*=2
    assert mul(a,out,P)=={0:1}
    return out

def det(A,P):
    out={}
    for p in itertools.permutations(range(len(A))):
        v={0:1}
        for i,j in enumerate(p):v=mul(v,A[i][j],P)
        out=add(out,v)
    return out

def vector(blocks,P,sc=1,shift=0):
    out=0
    for b,p in enumerate(blocks):
        for e,c in p.items():
            bits=COORD[MUL[sc][c]]
            for t in (0,1):
                if bits>>t&1:out^=1<<((2*b+t)*P+(e+shift)%P)
    return out

def support(v):
    out=[]
    while v:
        lo=v&-v;out.append(lo.bit_length()-1);v^=lo
    return out

def rows(A,P):
    out=[]
    for row in A:
        for alpha in (0,1):
            for h in range(P):
                v=0
                for j,p in enumerate(row):
                    for e,c in p.items():
                        for beta in (0,1):
                            if RHO[c][alpha][beta]:v^=1<<((2*j+beta)*P+(h-e)%P)
                out.append(v)
    return out

def echelon(vs):
    piv={}
    for v in vs:
        while v:
            p=v.bit_length()-1
            if p not in piv:piv[p]=v;break
            v^=piv[p]
    return piv

def balance(H,P):
    out=[]
    for b in range(len(H)//(2*P)):
        x,y=H[2*b*P],H[(2*b+1)*P]
        sel=sorted(range(3),key=lambda j:((x,y,x^y)[j].bit_count(),j))[:2]
        for j in sel:
            for g in range(P):
                x,y=H[2*b*P+g],H[(2*b+1)*P+g]
                out.append((x,y,x^y)[j])
    return out

class Constituent:
    def __init__(self,C,E,P):
        self.C,self.E,self.P=C,E,P
        self.r,self.s=len(C),len(C[0])
        self.A=[[{E[i][j]%P:C[i][j]} if C[i][j] else {} for j in range(self.s)] for i in range(self.r)]

    def basis(self):
        candidates=[]
        for pivot in itertools.combinations(range(self.s),self.r):
            D=[[row[j] for j in pivot] for row in self.A]
            delta=det(D,self.P);di=inverse(delta,self.P)
            if di is None:continue
            free=[j for j in range(self.s) if j not in pivot];U=[]
            for j in free:
                u=[{} for _ in range(self.s)];u[j]={0:1}
                for i,a in enumerate(pivot):
                    repl=[[self.A[h][j] if k==i else D[h][k] for k in range(self.r)] for h in range(self.r)]
                    u[a]=mul(det(repl,self.P),di,self.P)
                U.append(u)
            ws=[[sum(weight(p,c) for p in u) for c in (1,2,3)] for u in U]
            W=max(w for triple in ws for w in triple[1:])
            candidates.append((W,pivot,free,U,ws,delta))
        if not candidates:return None
        W,pivot,free,U,ws,delta=min(candidates,key=lambda v:(v[0],len(v[-1]),v[1]))
        return {'width':W,'pivot':pivot,'free':free,'U':U,'weights':ws,'determinant':delta}

    def upper(self):
        """Cofactor words and every cyclic power-of-two quotient lift."""
        best=10**9;bestword=None
        mod=self.P
        while mod:
            q=self.P//mod
            AA=[]
            for row in self.A:
                rr=[]
                for p in row:
                    a={}
                    for e,c in p.items():a=add(a,{e%mod:c})
                    rr.append(a)
                AA.append(rr)
            for cols in itertools.combinations(range(self.s),self.r+1):
                u=[{} for _ in range(self.s)]
                for j in cols:
                    d=det([[row[k] for k in cols if k!=j] for row in AA],mod)
                    u[j]={e+mod*t:c for e,c in d.items() for t in range(q)}
                if not any(u):continue
                for sc in (1,2,3):
                    v=vector(u,self.P,sc);w=v.bit_count()
                    if 0<w<best:best,bestword=w,v
            mod//=2
        return best,bestword

    def binary_rows(self):return rows(self.A,self.P)

def lp(A,C):
    assert A.P==C.P
    P=A.P;r,s=A.r,A.s;t,u=C.r,C.s;L=s*u+r*t
    HX=[];HZ=[]
    for i in range(r):
        for b in range(u):
            row=[{} for _ in range(L)]
            for a in range(s):row[a*u+b]=A.A[i][a]
            for j in range(t):row[s*u+i*t+j]=bar(C.A[j][b],P)
            HX.append(row)
    for a in range(s):
        for j in range(t):
            row=[{} for _ in range(L)]
            for b in range(u):row[a*u+b]=C.A[j][b]
            for i in range(r):row[s*u+i*t+j]=bar(A.A[i][a],P)
            HZ.append(row)
    ba,bc=A.basis(),C.basis();X=[];Z=[];XS=[];ZS=[]
    for ii,i in enumerate(ba['free']):
        for jj,j in enumerate(bc['free']):
            x=[{} for _ in range(L)];z=[{} for _ in range(L)]
            for b,p in enumerate(bc['U'][jj]):x[i*u+b]=p
            for a,p in enumerate(ba['U'][ii]):z[a*u+j]=p
            for sc in (2,3):
                XS.append(vector(x,P,sc));ZS.append(vector(z,P,sc))
                X.extend(vector(x,P,sc,g) for g in range(P))
                Z.extend(vector(z,P,sc,g) for g in range(P))
    return rows(HX,P),rows(HZ,P),X,Z,XS,ZS

def audit(A,C):
    HX,HZ,X,Z,XS,ZS=lp(A,C);P=A.P
    n=2*P*(A.s*C.s+A.r*C.r);k=2*P*(A.s-A.r)*(C.s-C.r)
    bx,bz=balance(HX,P),balance(HZ,P)
    assert len(echelon(HX))==2*P*A.r*C.s
    assert len(echelon(HZ))==2*P*A.s*C.r
    assert n-len(echelon(HX))-len(echelon(HZ))==k
    assert all((x&z).bit_count()%2==0 for x in HX for z in HZ)
    assert len(X)==len(Z)==k
    assert all((x&h).bit_count()%2==0 for x in X for h in HZ)
    assert all((z&h).bit_count()%2==0 for z in Z for h in HX)
    assert all((x&z).bit_count()%2==int(i==j) for i,x in enumerate(X) for j,z in enumerate(Z))
    assert len(echelon(HX+bx))==len(echelon(bx))==len(echelon(HX))
    assert len(echelon(HZ+bz))==len(echelon(bz))==len(echelon(HZ))
    return {'n':n,'k':k,'rank_X':len(echelon(HX)),'rank_Z':len(echelon(HZ)),
            'commute':True,'full_translated_kernels':True,'full_translated_Gram_identity':True,
            'raw_check_histogram':dict(sorted(Counter(map(int.bit_count,HX+HZ)).items())),
            'balanced_check_histogram':dict(sorted(Counter(map(int.bit_count,bx+bz)).items())),
            'canonical_frame_width':max(map(int.bit_count,XS+ZS)),
            'canonical_seed_weights_X':list(map(int.bit_count,XS)),
            'canonical_seed_weights_Z':list(map(int.bit_count,ZS)),
            'distance_lower_bound_certified':False}

def write_sample(path,n,H,G,initial=()):
    with open(path,'w') as f:
        f.write(f'{n} {len(H)} {len(G)}\n')
        for v in H+G:
            s=support(v);f.write(str(len(s))+' '+' '.join(map(str,s))+'\n')
        f.write(str(len(initial))+'\n')
        for v in initial:
            s=support(v);f.write(str(len(s))+' '+' '.join(map(str,s))+'\n')
