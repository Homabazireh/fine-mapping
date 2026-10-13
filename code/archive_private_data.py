#!/usr/bin/env python3
"""Package a local project for a PRIVATE release; never uploads anything.

Requires GNU tar, pigz and split. Archive pieces are individually below 2 GiB.
The JSON manifest records every source file's checksum or symlink target.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess

EXCLUDED = {'.git', '.venv', '.tools', '__pycache__', '.pytest_cache', '.quarto',
            '.cache', '.vscode', 'node_modules', '.Rproj.user'}
SECRET = re.compile(rb'(?:gh[pousr]_[A-Za-z0-9]{30,}|github_pat_[A-Za-z0-9_]{40,}|'
                    rb'-----BEGIN (?:RSA |EC |OPENSSH )?PRIVATE KEY-----|AKIA[A-Z0-9]{16})')


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def package(source, out, prefix, exclude_child=()):
    source, out = source.resolve(), out.resolve()
    if out == source or source in out.parents:
        raise ValueError('Archive destination must be outside its source tree')
    if out.exists() and any(out.iterdir()):
        raise ValueError('Use a new or empty archive destination')
    out.mkdir(parents=True, exist_ok=True)
    entries, omitted, files = [], [], []
    for directory, dirs, names in os.walk(source, followlinks=False):
        base = Path(directory)
        dirs[:] = sorted(d for d in dirs if d not in EXCLUDED and
                         not (base == source and d in exclude_child))
        links = [d for d in dirs if (base / d).is_symlink()]
        dirs[:] = [d for d in dirs if d not in links]
        names.extend(links)
        for name in sorted(names):
            path = base / name
            relative = path.relative_to(source.parent)
            if name.startswith('.env') or path.suffix in {'.key', '.pem'} or name in {'hosts.yml', '.netrc', '.git-credentials'}:
                omitted.append(str(relative))
                continue
            if path.is_symlink():
                target = path.resolve()
                # Preserve relative in-project links, including links to the enclosing project.
                if source.parent != target and source.parent not in target.parents:
                    raise ValueError(f'External symlink needs explicit review: {relative}')
                entry = {'path': str(relative), 'type': 'symlink', 'target': os.readlink(path)}
            elif path.is_file():
                before = path.stat()
                if before.st_size < 10 * 1024 * 1024 and path.suffix.lower() in {
                        '.py', '.r', '.sh', '.json', '.ini', '.yaml', '.yml', '.txt', '.md', '.rmd', '.ipynb'}:
                    if SECRET.search(path.read_bytes()):
                        raise ValueError(f'Possible credential detected; review file locally: {relative}')
                entry = {'path': str(relative), 'type': 'file', 'bytes': before.st_size,
                         'sha256': digest(path), 'mtime_ns': before.st_mtime_ns}
                after = path.stat()
                if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
                    raise ValueError(f'Input changed during checksum: {relative}')
            else:
                continue
            entries.append(entry)
            files.append(str(relative))
    manifest = {'source_directory': source.name, 'files': entries,
                'excluded_directory_names': sorted(EXCLUDED), 'excluded_top_level_children': list(exclude_child),
                'excluded_credential_paths': omitted, 'file_count': len(entries),
                'uncompressed_bytes': sum(x.get('bytes', 0) for x in entries)}
    (out / f'{prefix}.manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    listing = out / 'archive-files.nul'
    listing.write_bytes(b'\0'.join(os.fsencode(f) for f in files) + b'\0')
    print(f'Checked {len(entries)} source entries, {manifest["uncompressed_bytes"]} bytes', flush=True)
    tar = subprocess.Popen(['tar', '-C', str(source.parent), '--no-recursion', '--null',
                            '-T', str(listing), '-cf', '-'], stdout=subprocess.PIPE)
    gzip = subprocess.Popen(['pigz', '-p', '4', '-3'], stdin=tar.stdout, stdout=subprocess.PIPE)
    tar.stdout.close()
    split = subprocess.Popen(['split', '-b', '1900M', '-d', '-a', '3', '-',
                              str(out / f'{prefix}.tar.gz.part-')], stdin=gzip.stdout)
    gzip.stdout.close()
    codes = [split.wait(), gzip.wait(), tar.wait()]
    if any(codes):
        raise RuntimeError(f'Archive pipeline failed: {codes}')
    for entry in entries:
        path = source.parent / entry['path']
        if entry['type'] == 'file':
            stat = path.stat()
            if (stat.st_size, stat.st_mtime_ns) != (entry['bytes'], entry['mtime_ns']):
                raise ValueError(f'Input changed while packaging: {entry["path"]}')
        elif os.readlink(path) != entry['target']:
            raise ValueError(f'Symlink changed while packaging: {entry["path"]}')
    listing.unlink()
    assets = sorted(out.glob(f'{prefix}.tar.gz.part-*')) + [out / f'{prefix}.manifest.json']
    (out / f'{prefix}.SHA256SUMS').write_text(''.join(f'{digest(p)}  {p.name}\n' for p in assets))
    print(f'Archive ready: {len(assets)-1} pieces, {sum(p.stat().st_size for p in assets)} bytes', flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--prefix', required=True)
    parser.add_argument('--exclude-child', action='append', default=[])
    args = parser.parse_args()
    if not re.fullmatch(r'[A-Za-z0-9_-]+', args.prefix):
        parser.error('Prefix must use letters, digits, underscores and hyphens only')
    package(args.source, args.out, args.prefix, args.exclude_child)
