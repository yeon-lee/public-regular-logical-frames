#!/usr/bin/env python3
"""Validate the portable bound record; --run recompiles and repeats exclusion."""

if not __debug__:
    raise RuntimeError(
        'Verification requires Python assertions; run without -O/-OO and unset PYTHONOPTIMIZE.'
    )

import argparse,hashlib,json,re,subprocess,sys,tempfile,time
from pathlib import Path
HERE=Path(__file__).resolve().parent
CODE=HERE.parent
ROOT=HERE.parents[3]
sys.path.insert(0,str(ROOT))
from rlf import echelon,from_support,in_span,reconstruct,support,translate,verify


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
def read(path):return json.loads(path.read_text())
def parse_input(path):
    lines=iter(path.read_text().splitlines());n,P,blocks=map(int,next(lines).split());groups=[]
    for _ in range(3):
        count=int(next(lines));groups.append([from_support(list(map(int,next(lines).split())),n) for _ in range(count)])
    assert not list(lines)
    return n,P,blocks,groups
def completed(stdout,stderr,W,blocks):
    pattern=rf'COMPLETE W={W} roots=\[0,{blocks}\) children=\[0,1073741824\) mod\(D=-1,M=1,r=0\) nodes=(\d+) t=[0-9.]+'
    match=re.fullmatch(pattern,stdout.strip());assert match,stdout
    roots=[tuple(map(int,x)) for x in re.findall(r'^root block (\d+) done: nodes (\d+) \(total (\d+)\)',stderr,re.M)]
    assert [x[0] for x in roots]==list(range(blocks))
    total=0
    for _,count,cumulative in roots:total+=count;assert total==cumulative
    assert total==int(match.group(1))
    return total
def validate():
    cert=read(HERE/'certificate.json');receipt=read(HERE/'receipt.json')
    for name,digest in cert['files'].items():assert sha(HERE/name)==digest,name
    for name,digest in cert['code_files'].items():assert sha(CODE/name)==digest,name
    assert receipt['input_sha256']==sha(HERE/'input.txt')
    assert receipt['source_sha256']==sha(HERE/'excl.cpp')
    assert receipt['stdout_sha256']==sha(HERE/'stdout.txt') and receipt['stderr_sha256']==sha(HERE/'stderr.txt')
    assert receipt['original_receipt_sha256']==sha(HERE/'original_receipt.json')
    assert receipt['returncode']==0
    code=read(CODE/'code.json');checks=read(CODE/code['files']['checks']);frame=read(CODE/code['files']['frame'])
    algebra=verify(code,checks,frame);hx,hz,seeds=reconstruct(code);P=code['group']['order']
    n,period,blocks,(H,S,X)=parse_input(HERE/'input.txt')
    assert (n,period,blocks)==(code['n'],P,code['n']//P)
    # The archived search used all three binary rows per quaternary check.
    # Exact equality of their span with the published independent pair binds
    # the distance execution to the catalogue code.
    for archived,published in [(H,hx),(S,hz)]:
        rank=len(echelon(published));assert len(echelon(archived))==len(echelon(archived+published))==rank
    detectors=[translate(x,g,P) for x,z in seeds for g in range(P)]
    assert X==detectors and len(X)==code['k']
    for rows in (H,S):assert set(translate(v,1,P) for v in rows)==set(rows)
    # verify() checks the weight-preserving block-transpose sector isometry,
    # every translated canonical pairing and the physical upper witness.
    W=code['distance']['lower']-1;assert W+1==code['distance']['upper']==cert['distance']
    nodes=completed((HERE/'stdout.txt').read_text(),(HERE/'stderr.txt').read_text(),W,blocks)
    assert nodes==cert['exclusion']['nodes']
    return code,W,blocks,nodes,algebra
def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--run',action='store_true');p.add_argument('--compiler',default='c++');p.add_argument('--seconds',type=float,default=90);p.add_argument('--output',type=Path)
    args=p.parse_args();code,W,blocks,nodes,algebra=validate()
    result=dict(code=code['id'],status='PRESERVED_RECORD_VALIDATED',distance=code['distance']['lower'],archived_nodes=nodes,algebra=algebra,distance_search_rerun=False)
    if args.run:
        if args.seconds<=0:p.error('--seconds must be positive')
        with tempfile.TemporaryDirectory(prefix='rlf-lpq-') as tmp:
            binary=Path(tmp)/'excl';build=[args.compiler,'-O3','-std=c++17',str(HERE/'excl.cpp'),'-o',str(binary)]
            subprocess.run(build,check=True);start=time.monotonic()
            command=[str(binary),str(HERE/'input.txt'),str(W),'0',str(blocks),str(args.seconds)]
            run=subprocess.run(command,capture_output=True,text=True)
            assert run.returncode==0
            replay_nodes=completed(run.stdout,run.stderr,W,blocks)
            result.update(status='EXACT_DISTANCE_REPLAY_COMPLETE',distance_search_rerun=True,elapsed_seconds=time.monotonic()-start,nodes=replay_nodes,
                input_sha256=sha(HERE/'input.txt'),source_sha256=sha(HERE/'excl.cpp'),executable_sha256=sha(binary),stdout_sha256=hashlib.sha256(run.stdout.encode()).hexdigest(),stderr_sha256=hashlib.sha256(run.stderr.encode()).hexdigest())
            if args.output:
                args.output.mkdir(parents=True,exist_ok=True);(args.output/'stdout.txt').write_text(run.stdout);(args.output/'stderr.txt').write_text(run.stderr)
    if args.output:
        args.output.mkdir(parents=True,exist_ok=True);(args.output/'validation.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='algebra'},indent=2))
if __name__=='__main__':main()
