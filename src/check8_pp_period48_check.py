"""Reconstruct the P48 check-eight codes and audit or replay full exclusions."""

if not __debug__:
    raise RuntimeError(
        'Verification requires Python assertions; run without -O/-OO and unset PYTHONOPTIMIZE.'
    )

from pathlib import Path
import argparse,json,sys,tempfile,subprocess

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE/'src'))
from audit_pp import reconstruct
from binary_audit import read_matrix,bits,check_witness,sha,basis,remainder,support

def check_case(name,temp,replay=False):
    d=HERE/name;c=json.loads((d/'input.json').read_text());cert=json.loads((d/'certificate.json').read_text())
    for path,h in cert['hashes'].items():assert sha(HERE/path)==h,path
    rec,a=reconstruct(d/'input.json',temp/name,with_cofactors=False)
    assert (temp/name/'code.mat').read_bytes()==(d/'code.mat').read_bytes()
    assert (temp/name/'redundant3.mat').read_bytes()==(d/'redundant3.mat').read_bytes()
    n,HX,HZ=read_matrix(d/'code.mat');_,RX,RZ=read_matrix(d/'redundant3.mat')
    assert n==768 and a['k']==192 and a['ranks']==[288,288]
    assert a['physical_components']==[768] and a['X_check_max']==a['Z_check_max']==8
    assert a['X_column_max']==a['Z_column_max']==4
    assert a['full_translated_pairing_entries']==192**2
    assert a['displayed_width']==cert['frame_width']
    for H,R in [(HX,RX),(HZ,RZ)]:
        b=basis(H);assert len(b)==len(basis(R)) and all(not remainder(v,b) for v in R)
    def reflect(v):return sum(1<<(((j//48+8)%16)*48+(-j)%48) for j in support(v))
    assert set(map(reflect,HX))==set(HZ) and set(map(reflect,HZ))==set(HX)
    w=bits(c['witness_Z'],n);check_witness(w,n,HX,HZ,cert['distance']);check_witness(reflect(w),n,HZ,HX,cert['distance'])
    receipt=json.loads((d/'receipt.json').read_text())
    assert receipt['returncode']==0 and not receipt['witnesses']
    for key,path in [('execution_input_sha256',d/'redundant3.mat'),('original_matrix_sha256',d/'code.mat'),
                     ('source_sha256',HERE/'src/exact.cpp'),('header_sha256',HERE/'src/common.hpp'),
                     ('engine_sha256',HERE/'evidence/exact'),('stdout_sha256',d/'stdout.jsonl')]:assert receipt[key]==sha(path),key
    events=[json.loads(l) for l in (d/'stdout.jsonl').read_text().splitlines()]
    f=receipt['final'];assert events==[f] and f['event']=='result' and f['status']=='excluded'
    assert f['W']==cert['distance']-1 and f['roots_completed']==16 and f['coset_pruning']
    assert (f['partition'],f['shard'],f['shards'],f['split_weight'])==('hash_prefix_v1',0,1,1)
    assert f['frontier_seen']==f['frontier_owned']==16
    cmd=receipt['command'];assert len(cmd)==10 and cmd[2:4]==['48',str(cert['distance']-1)]
    assert cmd[5:]==['--shard','0','1','1','18446744073709551615']
    assert not (d/'stderr.txt').read_text().strip()
    result=dict(code=cert['id'],status='PASS',n=n,k=192,distance_exact=cert['distance'],maximum_check_weight=8,
                maximum_column_weight=4,reported_frame_width=cert['frame_width'],frame_optimality_proved=False,
                recorded_nodes=f['nodes'],deep_search_replayed=False)
    if replay:
        command=[str(temp/'exact'),str(d/'redundant3.mat'),'48',str(cert['distance']-1),'60','--shard','0','1','1','18446744073709551615']
        p=subprocess.run(command,capture_output=True,text=True)
        assert p.returncode==0,p.stdout+p.stderr
        events=[json.loads(l) for l in p.stdout.splitlines()];assert len(events)==1
        fresh=events[0];assert fresh['status']=='excluded' and fresh['W']==cert['distance']-1 and fresh['roots_completed']==16
        result.update(deep_search_replayed=True,replay=fresh)
    return result

def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--replay',action='store_true');parser.add_argument('--report',type=Path);args=parser.parse_args()
    manifest=HERE/'SHA256SUMS.json'
    if manifest.exists():
        for path,h in json.loads(manifest.read_text())['files'].items():assert sha(HERE/path)==h,path
    with tempfile.TemporaryDirectory(prefix='pp48-check8-') as td:
        temp=Path(td)
        if args.replay:subprocess.run(['c++','-O3','-std=c++17',str(HERE/'src/exact.cpp'),'-o',str(temp/'exact')],check=True)
        results=[check_case(name,temp,args.replay) for name in json.loads((HERE/'cases.json').read_text())]
    if args.report:args.report.write_text(json.dumps(results,indent=2)+'\n')
    print(json.dumps(results,indent=2))

if __name__=='__main__':main()
