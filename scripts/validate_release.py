#!/usr/bin/env python3
"""Check packaged files against the release manifest; ignore generated results."""
import hashlib
import json
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[1]

def main():
    manifest=ROOT/'config/release_manifest.json'
    if not manifest.is_file():
        raise RuntimeError('Release manifest is missing.')
    records=json.loads(manifest.read_text())
    failed=[]
    for name,digest in records.items():
        p=ROOT/name
        if not p.is_file():failed.append('missing '+name)
        elif hashlib.sha256(p.read_bytes()).hexdigest()!=digest:failed.append('changed '+name)
    if failed:
        raise RuntimeError('Release differs from its recorded source snapshot:\n'+'\n'.join(failed))
    print('Verified {} packaged files; generated results and downloaded data are outside this manifest.'.format(len(records)))
if __name__=='__main__':
    try:main()
    except (RuntimeError,OSError) as exc:sys.exit(str(exc))
