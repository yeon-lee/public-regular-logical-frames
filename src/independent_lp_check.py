"""Portable algebra, witness, completed-root record audit and optional replay."""

if not __debug__:
    raise RuntimeError(
        'Verification requires Python assertions; run without -O/-OO and unset PYTHONOPTIMIZE.'
    )

import argparse,hashlib,json,re,subprocess,sys,tempfile
from pathlib import Path
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE/'source'))
from independent_lpq import Constituent,audit,lp,echelon


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
def read(path):
    ls=path.read_text().splitlines();n,P,L=map(int,ls[0].split());i=1;out=[]
    for _ in range(3):
        count=int(ls[i]);i+=1;mat=[]
        for line in ls[i:i+count]:
            supp=list(map(int,line.split()));assert len(set(supp))==len(supp)
            assert all(0<=j<n for j in supp);mat.append(sum(1<<j for j in supp))
        i+=count;out.append(mat)
    assert i==len(ls)
    return n,P,L,out
def shift(v,P,L):
    out=0;mask=(1<<P)-1
    for b in range(L):
        w=(v>>(b*P))&mask;out|=(((w<<1)|(w>>(P-1)))&mask)<<(b*P)
    return out
def terminal(stdout,stderr):
    m=re.fullmatch(r'COMPLETE W=13 roots=\[0,68\) children=\[0,1073741824\) mod\(D=-1,M=1,r=0\) nodes=(\d+) t=([\d.]+)',stdout.strip())
    assert m,stdout
    assert stderr.startswith('n=544 P=8 blocks=68 rX=360 rZ=360 kX=64')
    roots=[tuple(map(int,z)) for z in re.findall(r'root block (\d+) done: nodes (\d+) \(total (\d+)\)',stderr)]
    assert [z[0] for z in roots]==list(range(68))
    assert sum(z[1] for z in roots)==int(m[1])==roots[-1][2]
    total=0
    for _,nodes,cumulative in roots:total+=nodes;assert total==cumulative
    return int(m[1])
def main():
    parser=argparse.ArgumentParser();parser.add_argument('--replay',action='store_true');parser.add_argument('--seconds',type=float,default=90)
    args=parser.parse_args();assert args.seconds>0
    files=json.loads((HERE/'SHA256SUMS.json').read_text())['files']
    for name,digest in files.items():
        path=HERE/name;assert not Path(name).is_absolute() and path.resolve().is_relative_to(HERE)
        assert sha(path)==digest,name
    d=json.loads((HERE/'candidate.json').read_text());P=d['P']
    A=Constituent(d['A']['C'],d['A']['E'],P);B=Constituent(d['B']['C'],d['B']['E'],P)
    algebra=audit(A,B);assert (algebra['n'],algebra['k'],algebra['canonical_frame_width'])==(544,64,16)
    assert max(algebra['raw_check_histogram'])==max(algebra['balanced_check_histogram'])==9
    HX,HZ,X,Z,XS,ZS=lp(A,B)
    for sec,check,stab,dual in [('X',HZ,HX,Z),('Z',HX,HZ,X)]:
        inp=HERE/'evidence'/f'excl_{sec}.txt';n,pp,L,(H,G,F)=read(inp)
        assert (n,pp,L)==(544,8,68)
        assert len(echelon(H))==len(echelon(H+check))==240
        assert len(echelon(G))==len(echelon(G+stab))==240
        assert len(echelon(G+F))==304
        assert len(echelon(G+F+dual))==304
        assert all((h&g).bit_count()%2==0 for h in H for g in G)
        assert all((f&g).bit_count()%2==0 for f in F for g in G)
        assert set(shift(v,P,L) for v in H)==set(H)
        assert set(shift(v,P,L) for v in G)==set(G)
        w=sum(1<<j for j in d[f'{sec}_logical_witness'])
        assert w.bit_count()==14
        assert all((w&h).bit_count()%2==0 for h in H)
        assert any((w&f).bit_count()%2 for f in F)
        rec=json.loads((HERE/'evidence'/f'{sec}_receipt.json').read_text())
        for key,path in [('input_sha256',inp),('source_sha256',HERE/'source'/'excl.cpp'),('executable_sha256',HERE/'evidence'/'recorded_excl.bin'),('stdout_sha256',HERE/'evidence'/f'{sec}_stdout.log'),('stderr_sha256',HERE/'evidence'/f'{sec}_stderr.log')]:assert rec[key]==sha(path)
        assert rec['returncode']==0
        assert rec['command']==['evidence/recorded_excl.bin',f'evidence/excl_{sec}.txt','13','0','68','90']
        nodes=terminal((HERE/'evidence'/f'{sec}_stdout.log').read_text(),(HERE/'evidence'/f'{sec}_stderr.log').read_text())
        print(f'{sec}: complete exclusion through13, all68 roots; {nodes:,} nodes; physical weight14 witness verified.')
    print('PASS [[544,64,14]], check max9, achieved full regular canonical width16 (not an optimality claim).')
    if args.replay:
        with tempfile.TemporaryDirectory(prefix='lpq544-replay-') as tmp:
            exe=Path(tmp)/'excl'
            subprocess.run(['c++','-O3','-std=c++17',str(HERE/'source'/'excl.cpp'),'-o',str(exe)],check=True)
            for sec in ('X','Z'):
                result=subprocess.run([str(exe),str(HERE/'evidence'/f'excl_{sec}.txt'),'13','0','68',str(args.seconds)],capture_output=True,text=True)
                assert result.returncode==0
                if result.stdout.startswith('TIMEOUT'):
                    print(f'{sec}: replay incomplete; no new lower bound established.');return 75
                nodes=terminal(result.stdout,result.stderr);print(f'{sec}: fresh replay complete ({nodes:,} nodes).')
    return 0
if __name__=='__main__':sys.exit(main())
