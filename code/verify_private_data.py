#!/usr/bin/env python3
"""Verify restored files against a private archive's source manifest."""
import argparse
import hashlib
import json
import os
from pathlib import Path


def verify(manifest_path, restored_root):
    manifest = json.loads(Path(manifest_path).read_text())
    root = Path(restored_root).resolve()
    files, links = 0, 0
    for entry in manifest['files']:
        relative = Path(entry['path'])
        if relative.is_absolute() or '..' in relative.parts:
            raise ValueError('Unsafe manifest path')
        path = root / relative
        if entry['type'] == 'symlink':
            if not path.is_symlink() or os.readlink(path) != entry['target']:
                raise ValueError(f'Symlink mismatch: {relative}')
            links += 1
        else:
            if path.is_symlink() or not path.is_file() or path.stat().st_size != entry['bytes']:
                raise ValueError(f'File type/size mismatch: {relative}')
            digest = hashlib.sha256()
            with path.open('rb') as handle:
                for block in iter(lambda: handle.read(8 * 1024 * 1024), b''):
                    digest.update(block)
            if digest.hexdigest() != entry['sha256']:
                raise ValueError(f'Checksum mismatch: {relative}')
            files += 1
    print(f'Verified {files} files and {links} symlinks against source manifest')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest', type=Path, required=True)
    parser.add_argument('--root', type=Path, required=True,
                        help='Restore directory containing the archive top-level directory')
    args = parser.parse_args()
    verify(args.manifest, args.root)
