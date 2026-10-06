if not __debug__:
    raise RuntimeError(
        'Verification requires Python assertions; run without -O/-OO and unset PYTHONOPTIMIZE.'
    )

from pathlib import Path
import hashlib,json,subprocess
H=Path(__file__).resolve().parent
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
m=json.loads((H/'bundle_sha256.json').read_text())
for n,h in m['files'].items():assert sha(H/n)==h
subprocess.run(['g++','-O3','-std=c++17','exact_target.cpp','-o','exact_target'],cwd=H,check=True)
r=dict(executable_sha256=sha(H/'exact_target'),bundle_sha256=sha(H/'bundle_sha256.json'),source_sha256=sha(H/'exact_target.cpp'))
(H/'build.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r))
