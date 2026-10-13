#!/usr/bin/env python3
"""Retrieve the exact CGMacros v1.0.0 CSV inputs, validating the official ZIP hash."""
import argparse
import hashlib
from pathlib import Path, PurePosixPath
import shutil
import urllib.request
import zipfile

BASE = 'https://physionet.org/files/cgmacros/1.0.0/'
NAME = 'CGMacros_dateshifted365.zip'
SHA256 = '05c8b0e6f1a2757050aced55ce4bf6ab2ac9b30f2fd8ca193056812d9c621d4d'
ROOT = Path(__file__).resolve().parents[1]

def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--archive', type=Path, help='Use an existing downloaded ZIP (no network request).')
    parser.add_argument('--output', type=Path, default=ROOT/'data/raw', help='Parent of the extracted CGMacros directory.')
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    archive = args.archive if args.archive else args.output/NAME
    if not archive.exists():
        if args.archive:
            parser.error('Archive does not exist: '+str(archive))
        partial = archive.with_suffix('.zip.part')
        print('Downloading CGMacros v1.0.0 (~627 MB) from '+BASE+NAME, flush=True)
        with urllib.request.urlopen(BASE+NAME) as source, partial.open('wb') as target:
            shutil.copyfileobj(source, target)
        if digest(partial) != SHA256:
            raise RuntimeError('Downloaded ZIP checksum does not match the official v1.0.0 manifest.')
        partial.replace(archive)
    print('Checking official SHA-256...', flush=True)
    if digest(archive) != SHA256:
        raise RuntimeError('ZIP checksum does not match CGMacros v1.0.0.')
    count = 0
    with zipfile.ZipFile(str(archive)) as z:
        for member in z.infolist():
            rel = PurePosixPath(member.filename)
            if rel.is_absolute() or '..' in rel.parts:
                raise RuntimeError('Unexpected archive member: '+member.filename)
            if not rel.parts or rel.parts[0] != 'CGMacros':
                continue
            # Only demographics and per-participant CSVs are needed; skip photographs.
            wanted = str(rel) == 'CGMacros/bio.csv' or (
                len(rel.parts) == 3 and rel.parts[1].startswith('CGMacros-') and
                rel.parts[2] == rel.parts[1]+'.csv')
            if not wanted:
                continue
            destination = args.output.joinpath(*rel.parts)
            destination.parent.mkdir(parents=True, exist_ok=True)
            with z.open(member) as source, destination.open('wb') as target:
                shutil.copyfileobj(source, target)
            count += 1
    if count != 46:
        raise RuntimeError('Expected bio.csv and 45 participant CSVs; extracted '+str(count))
    print('Ready: '+str(args.output/'CGMacros')+' (46 CSV files).')

if __name__ == '__main__':
    main()
