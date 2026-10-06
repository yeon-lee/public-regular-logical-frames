#!/usr/bin/env python3
"""Portable validation of a complete one-shard certificate; optional fresh replay."""

if not __debug__:
    raise RuntimeError(
        'Verification requires Python assertions; run without -O/-OO and unset PYTHONOPTIMIZE.'
    )

import argparse,hashlib,json,shutil,subprocess,sys,tempfile,time
from pathlib import Path
from types import SimpleNamespace
HERE=Path(__file__).resolve().parent;CODE=HERE.parent;ROOT=HERE.parents[3]
sys.path.insert(0,str(ROOT));sys.path.insert(0,str(HERE))
import rlf
import validate_hash_union as union

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
def read(p):return json.loads(p.read_text())
def write(p,d):p.write_text(json.dumps(d,indent=2)+'\n')
def args_for(folder,d):return SimpleNamespace(candidate_dir=folder/'inputs',receipts_dir=folder/'receipts',engine_dir=folder/'engine',executable=folder/'engine/exact',weight=d-1,shards=1,split_weight=5,sector='Z',pattern='single.receipt.json')
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--run',action='store_true');ap.add_argument('--output',type=Path);ap.add_argument('--seconds',type=float,default=180);ap.add_argument('--compiler',default='c++');a=ap.parse_args()
    manifest=read(HERE/'portable_manifest.json')
    for name,h in manifest['files'].items():assert sha(HERE/name)==h,name
    for name,h in manifest['code_files'].items():assert sha(CODE/name)==h,name
    code=read(CODE/'code.json');checks=read(CODE/'checks.json');frame=read(CODE/'frames/paper.json');physical=read(HERE/'inputs/physical.json')
    rlf.verify(code,checks,frame)
    assert checks['X']==physical['HX'] and checks['Z']==physical['HZ']
    assert [s['X'] for s in frame['seeds']]==physical['X_seeds']
    assert [s['Z'] for s in frame['seeds']]==physical['Z_seeds']
    d=code['distance']['lower'];result=union.validate(args_for(HERE,d))
    assert result['status']=='CERTIFIED_EXACT_DISTANCE' and result['exact_distance']==d==code['distance']['upper']
    algebra=result['algebra'];assert (algebra['n'],algebra['k'],algebra['frame_width'],algebra['balanced_check_max'])==(code['n'],code['k'],code['frame_width']['upper'],8)
    print(json.dumps({'code':code['id'],'status':result['status'],'exact_distance':d,'frame_width':algebra['frame_width'],'nodes':result['total_nodes'],'distance_search_rerun':False}))
    if not a.run:return 0
    if a.output is None or a.output.exists():ap.error('--run requires a NEW --output directory')
    assert a.seconds>0
    out=a.output.resolve();out.mkdir(parents=True);(out/'receipts').mkdir();(out/'engine').mkdir();shutil.copytree(HERE/'inputs',out/'inputs')
    for name in ('exact.cpp','common.hpp'):shutil.copy2(HERE/'engine'/name,out/'engine'/name)
    exe=out/'engine/exact';subprocess.run([a.compiler,'-O3','-std=c++17',str(out/'engine/exact.cpp'),'-o',str(exe)],check=True)
    matrix=out/'inputs/exact_Z.mat';cmd=[str(exe),str(matrix),str(code['group']['order']),str(d-1),str(a.seconds),'--shard','0','1','5','1000000000000']
    start=time.monotonic();p=subprocess.run(cmd,capture_output=True,text=True)
    stdout=out/'receipts/single.stdout';stderr=out/'receipts/single.stderr';stdout.write_text(p.stdout);stderr.write_text(p.stderr)
    receipt={'command':cmd,'returncode':p.returncode,'elapsed_seconds':time.monotonic()-start,'input_sha256':sha(matrix),
      'source_sha256':{f:sha(out/'engine'/f) for f in ('exact.cpp','common.hpp')},'executable_sha256':sha(exe),'stdout_sha256':sha(stdout),'stderr_sha256':sha(stderr)}
    write(out/'receipts/single.receipt.json',receipt)
    if p.returncode==75:print('Replay incomplete; no new lower bound established.');return 75
    result=union.validate(args_for(out,d));result['distance_search_rerun']=True
    for row in result['receipts']:row['path']=str(Path(row['path']).relative_to(out))
    write(out/'validation.json',result);print(json.dumps({'status':result['status'],'exact_distance':result['exact_distance'],'nodes':result['total_nodes'],'distance_search_rerun':True}))
    return 0
if __name__=='__main__':sys.exit(main())
