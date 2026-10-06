"""Integrity of public files and explicitly declared omitted executables."""
from pathlib import Path
import hashlib
import json

def evidence_digest(path):
    """Hash an existing file, or return a declared (unverified) binary identity.

    This never treats a recorded executable fingerprint as a checked binary.
    Only explicitly listed omitted executables may use the declaration branch.
    """
    path = Path(path).resolve()
    if path.is_file():
        return hashlib.sha256(path.read_bytes()).hexdigest()
    for root in path.parents:
        manifest = root / 'docs' / 'EXECUTABLE_IDENTITIES.json'
        if manifest.is_file():
            name = str(path.relative_to(root))
            identities = json.loads(manifest.read_text())['executables']
            if name in identities:
                return identities[name]
            break
    raise FileNotFoundError(path)
