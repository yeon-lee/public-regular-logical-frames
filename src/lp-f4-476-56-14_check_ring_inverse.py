"""Cross-check cyclic-ring inversion against independent binary matrix rank."""

if not __debug__:
    raise RuntimeError(
        'Verification requires Python assertions; run without -O/-OO and unset PYTHONOPTIMIZE.'
    )

import itertools,json,random,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[4]
sys.path.insert(0,str(ROOT))
from rlf import Ring,binary_checks,echelon,reconstruct,support

def check(P,coefficients):
    a={e:c for e,c in enumerate(coefficients) if c};ring=Ring(P)
    binary_unit=len(echelon(binary_checks([[a]],P,4)))==2*P
    try:inverse=ring.inverse(a)
    except ValueError:assert not binary_unit;return False
    assert binary_unit and ring.mul(a,inverse)=={0:1} and ring.mul(inverse,a)=={0:1}
    return True

def main():
    units=sum(check(3,c) for c in itertools.product(range(4),repeat=3));assert units==27
    rng=random.Random(0);random_checks={}
    for P in (6,7,8,11):
        values=[check(P,[rng.randrange(4) for _ in range(P)]) for _ in range(200)]
        assert any(values) and not all(values)
        random_checks[P]={'checked':len(values),'units':sum(values),'nonunits':len(values)-sum(values)}
    catalog=json.loads((ROOT/'catalog.json').read_text());count=0
    for entry in catalog['entries']:
        path=ROOT/entry['path'];code=json.loads(path.read_text());checks=json.loads((path.parent/code['files']['checks']).read_text());frame=json.loads((path.parent/code['files']['frame']).read_text())
        hx,hz,seeds=reconstruct(code)
        assert checks['X']==list(map(support,hx)) and checks['Z']==list(map(support,hz)),entry['id']
        assert frame['seeds']==[{'X':support(x),'Z':support(z)} for x,z in seeds],entry['id']
        count+=1
    result={'status':'PASS','period3_all64_polynomials':{'units':units,'nonunits':64-units},'random_binary_rank_crosschecks':random_checks,'catalogue_reconstructions_unchanged':count}
    print(json.dumps(result,indent=2))
if __name__=='__main__':main()
