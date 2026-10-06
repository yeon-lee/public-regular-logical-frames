#!/usr/bin/env python3
"""Validate the 64-shard certificate; --run makes a fresh portable replay."""
import argparse,hashlib,json,shutil,subprocess,sys
from pathlib import Path
from types import SimpleNamespace
HERE=Path(__file__).resolve().parent
CODE=HERE.parent
ROOT=HERE.parents[3]
sys.path.insert(0,str(ROOT))
sys.path.insert(0,str(HERE))
from rlf import reconstruct,support,verify
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

def sha(path):
    """File digest or explicitly declared omitted executable identity."""
    return evidence_digest(path)
def read(path):return json.loads(path.read_text())
def args_for(folder):
    return SimpleNamespace(candidate_dir=folder/'inputs',receipts_dir=folder/'receipts',engine_dir=folder/'engine',
        executable=folder/'engine/exact',weight=21,shards=64,split_weight=5,sector='Z',pattern='shard_*.receipt.json')
def validate():
    manifest=read(HERE/'portable_manifest.json')
    for name,digest in manifest['files'].items():
        if sha(HERE/name)!=digest:raise ValueError(f'Proof integrity mismatch: {name}')
    for name,digest in manifest['code_files'].items():
        if sha(CODE/name)!=digest:raise ValueError(f'Code integrity mismatch: {name}')
    code=read(CODE/'code.json');checks=read(CODE/code['files']['checks']);frame=read(CODE/code['files']['frame'])
    algebra=verify(code,checks,frame);hx,hz,seeds=reconstruct(code);physical=read(HERE/'inputs/physical.json')
    if checks['X']!=physical['HX'] or checks['Z']!=physical['HZ']:raise ValueError('Archived physical matrices differ from catalogue')
    if physical['X_seeds']!=[support(x) for x,z in seeds] or physical['Z_seeds']!=[support(z) for x,z in seeds]:raise ValueError('Archived canonical frame differs from catalogue')
    for receipt in (HERE/'receipts').glob('*.receipt.json'):
        if read(receipt)['driver_sha256']!=sha(HERE/'run_all_shards.py'):raise ValueError('Driver hash mismatch')
    result=union.validate(args_for(HERE))
    if result['status']!='CERTIFIED_EXACT_DISTANCE' or result['exact_distance']!=code['distance']['lower'] or result['exact_distance']!=code['distance']['upper']:raise ValueError('Distance bounds do not match catalogue')
    a=result['algebra']
    if (a['n'],a['k'],a['frame_width'],a['balanced_check_max'])!=(code['n'],code['k'],code['frame_width']['upper'],code['maximum_check_weight']):raise ValueError('Parameter mismatch')
    for row in result['receipts']:row['path']=str(Path(row['path']).relative_to(HERE))
    result.update(code_id=code['id'],catalogue_algebra=algebra,distance_search_rerun=False)
    return result
def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--run',action='store_true');p.add_argument('--output',type=Path);p.add_argument('--compiler',default='c++')
    a=p.parse_args();result=validate()
    if a.run:
        if a.output is None:p.error('--run requires a new --output directory')
        if a.output.exists():p.error('--output must not already exist; preserved evidence is never overwritten')
        out=a.output.resolve();out.mkdir(parents=True);(out/'engine').mkdir();(out/'receipts').mkdir()
        shutil.copytree(HERE/'inputs',out/'inputs')
        for name in ('exact.cpp','common.hpp'):shutil.copy2(HERE/'engine'/name,out/'engine'/name)
        for name in ('run_all_shards.py','validate_hash_union.py'):shutil.copy2(HERE/name,out/name)
        subprocess.run([a.compiler,'-O3','-std=c++17',str(out/'engine/exact.cpp'),'-o',str(out/'engine/exact')],check=True)
        subprocess.run([sys.executable,str(out/'run_all_shards.py')],check=True)
        result=union.validate(args_for(out));result['distance_search_rerun']=True
        for row in result['receipts']:row['path']=str(Path(row['path']).relative_to(out))
        (out/'validation.json').write_text(json.dumps(result,indent=2)+'\n')
    elif a.output:
        a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ('status','exact_distance','total_nodes','summed_seconds','distance_search_rerun')},indent=2))
if __name__=='__main__':main()
