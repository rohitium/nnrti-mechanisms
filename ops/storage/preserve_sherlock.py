#!/usr/bin/env python3
"""Preserve a quiescent Sherlock project; Python 3.6+, standard library only.

Run inside a Slurm allocation. Sources are never removed. Each archive has a
per-file inventory and a SHA-256 verified by rereading the destination archive.
An incomplete snapshot never receives COMPLETE.json.
"""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import stat
import subprocess
import tarfile
import time


def digest(path):
    h = hashlib.sha256()
    with open(str(path), 'rb') as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


class HashedStream:
    def __init__(self, stream):
        self.stream = stream
        self.hash = hashlib.sha256()

    def read(self, count):
        block = self.stream.read(count)
        self.hash.update(block)
        return block

    def write(self, block):
        result = self.stream.write(block)
        self.hash.update(block)
        return result


def max_step(folder):
    maximum = 0
    for path in folder.glob('*_md_state.csv'):
        with path.open() as stream:
            for row in csv.reader(stream):
                try:
                    maximum = max(maximum, int(row[0]))
                except (ValueError, IndexError):
                    pass
    return maximum


def archive(source, paths, destination, name, resume=False):
    partial = destination / (name + '.tar.part')
    final = destination / (name + '.tar')
    verified_path = destination / (name + '.verified.json')
    inventory = destination / (name + '.files.jsonl')
    if resume and final.exists() and verified_path.exists():
        record = json.loads(verified_path.read_text())
        if digest(final) != record['sha256'] or digest(inventory) != record['inventory_sha256']:
            raise RuntimeError('Previously verified archive changed: ' + name)
        previous = [json.loads(line) for line in inventory.read_text().splitlines()]
        if {item['path'] for item in previous} != {str(path) for path in paths}:
            raise RuntimeError('Source membership changed: ' + name)
        for item in previous:
            path = source/item['path']
            info = path.lstat()
            if item['type'] == '0' and (info.st_size, info.st_mtime_ns) != (item['bytes'], item['mtime_ns']):
                raise RuntimeError('Source changed since preservation: ' + str(path))
            if item['type'] == '2' and os.readlink(str(path)) != item['link_target']:
                raise RuntimeError('Source link changed since preservation: ' + str(path))
        print('Reverified completed archive: ' + name, flush=True)
        return record
    if resume and partial.exists() and not final.exists() and not verified_path.exists():
        # Only an unfinished generated archive is replaced; original data stays intact.
        partial.unlink()
        if inventory.exists():
            inventory.unlink()
    if final.exists() or partial.exists():
        raise RuntimeError('Refusing to overwrite archive: ' + str(final))
    total = 0
    with partial.open('xb') as output, inventory.open('x') as manifest:
        writer = HashedStream(output)
        with tarfile.open(fileobj=writer, mode='w|', format=tarfile.PAX_FORMAT) as tar:
            for relative in paths:
                path = source / relative
                before = path.lstat()
                info = tar.gettarinfo(str(path), arcname=str(relative))
                record = {'path': str(relative), 'bytes': before.st_size,
                          'mtime_ns': before.st_mtime_ns, 'type': info.type.decode('ascii')}
                if info.isfile():
                    with path.open('rb') as stream:
                        reader = HashedStream(stream)
                        tar.addfile(info, reader)
                        record['sha256'] = reader.hash.hexdigest()
                    after = path.stat()
                    if (before.st_size, before.st_mtime_ns, before.st_ino) != (after.st_size, after.st_mtime_ns, after.st_ino):
                        raise RuntimeError('Source changed while archiving: ' + str(path))
                    total += before.st_size
                else:
                    if info.issym():
                        record['link_target'] = info.linkname
                    tar.addfile(info)
                manifest.write(json.dumps(record, sort_keys=True) + '\n')
        output.flush()
        os.fsync(output.fileno())
        expected = writer.hash.hexdigest()
    actual = digest(partial)
    if actual != expected:
        raise RuntimeError('Archive verification failed: ' + name)
    partial.rename(final)
    result = {'archive': final.name, 'sha256': actual, 'bytes': final.stat().st_size,
              'source_bytes': total, 'entries': len(paths),
              'inventory': inventory.name, 'inventory_sha256': digest(inventory)}
    (destination / (name + '.verified.json')).write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result), flush=True)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--snapshot', required=True)
    parser.add_argument('--resume', action='store_true', help='After diagnosing an interruption, reverify completed archives and restart unfinished units.')
    args = parser.parse_args()
    if not os.environ.get('SLURM_JOB_ID'):
        raise RuntimeError('Run this preservation task in a Slurm allocation.')
    if not args.snapshot or args.snapshot in {'.', '..'} or Path(args.snapshot).name != args.snapshot:
        raise ValueError('Snapshot must be a single directory name.')
    scratch = Path(os.environ['SCRATCH'])
    project = Path(os.environ['GROUP_HOME']) / 'nnrti-mechanisms'
    destination = project / 'archive' / args.snapshot
    destination.mkdir(parents=True, exist_ok=args.resume)
    if args.resume and (destination/'COMPLETE.json').exists():
        print('Snapshot is already complete: ' + str(destination), flush=True)
        return
    jobs = subprocess.check_output(['squeue', '-h', '-u', os.environ['USER'], '-o', '%A']).decode().split()
    if set(jobs) - {os.environ['SLURM_JOB_ID']}:
        raise RuntimeError('Other jobs exist; inspect their write paths before preservation.')
    # Inventory once. Do not follow symlinks or scan any other user/group tree.
    allowed = {'nnrti-mechanisms-git', 'nnrti-mechanisms', 'nnrti_backups',
               'nnrti-venv', 'apo_collapse', 'apo_collapse_solvated',
               'apo_collapse_test', 'apo_solv_test', '_portable_probe', '_resume_probe'}
    groups = {}
    total = 0
    excluded_caches = {'.nv', '.pip'}
    for root, dirs, files in os.walk(str(scratch), followlinks=False):
        relroot = Path(root).relative_to(scratch)
        if relroot == Path('.'):
            dirs[:] = [name for name in dirs if name not in excluded_caches]
            unknown = set(dirs) - allowed
            if unknown:
                raise RuntimeError('Unreviewed scratch directories: ' + repr(unknown))
        for name in sorted(dirs + files):
            relative = relroot / name
            parts = relative.parts
            legprefix = ('nnrti-mechanisms-git', 'results', 'analysis', 'fep_pmx', 'legs')
            if tuple(parts[:5]) == legprefix and len(parts) >= 6:
                group = 'fep-' + parts[5]
            elif parts[0] == 'nnrti-mechanisms-git':
                group = 'main-other'
            elif parts[0] == 'nnrti-mechanisms':
                group = 'staging'
            else:
                group = 'auxiliary'
            info = (scratch / relative).lstat()
            if stat.S_ISREG(info.st_mode):
                total += info.st_size
            elif not (stat.S_ISDIR(info.st_mode) or stat.S_ISLNK(info.st_mode)):
                raise RuntimeError('Unsupported special file: ' + str(relative))
            groups.setdefault(group, []).append(relative)
    source_record = destination / ('RESUME-'+os.environ['SLURM_JOB_ID']+'.json' if args.resume else 'SOURCE.json')
    source_record.write_text(json.dumps({
        'source': str(scratch), 'logical_bytes': total,
        'excluded_root_caches': sorted(excluded_caches),
        'groups': {key: len(value) for key, value in groups.items()},
        'job_id': os.environ['SLURM_JOB_ID'], 'started_unix': time.time()}, indent=2) + '\n')
    # Leave ample room below the 1 TB group quota observed before submission.
    if total > 750 * 1024**3:
        raise RuntimeError('Source exceeds reviewed preservation budget (750 GiB).')
    # Preserve old durable copies before replacing the two stale replicas.
    repairs = []
    for replicate in ['rep_02', 'rep_03']:
        source = scratch / 'nnrti-mechanisms-git/results/md_runs/apo/v106a' / replicate
        target = project / 'md_runs/apo/v106a' / replicate
        if max_step(source) < 50000000 or max_step(source) < max_step(target):
            raise RuntimeError('Source replicate is incomplete or older than destination.')
        oldpaths = []
        if target.exists():
            for root, dirs, files in os.walk(str(target)):
                oldpaths.extend((Path(root) / name).relative_to(target.parent) for name in sorted(dirs + files))
            oldpaths.insert(0, Path(target.name))
            prior = destination / ('prior-v106a-' + replicate + '.verified.json')
            if args.resume and prior.exists():
                prior_record = json.loads(prior.read_text())
                if digest(destination/prior_record['archive']) != prior_record['sha256']:
                    raise RuntimeError('Prior MD copy failed verification.')
            else:
                archive(target.parent, sorted(oldpaths), destination, 'prior-v106a-' + replicate)
        subprocess.run(['rsync', '-a', str(source) + '/', str(target) + '/'], check=True)
        for root, dirs, files in os.walk(str(source)):
            for name in files:
                original = Path(root) / name
                copied = target / original.relative_to(source)
                if not original.is_symlink() and digest(original) != digest(copied):
                    raise RuntimeError('MD mirror verification failed: ' + str(copied))
        repairs.append({'replicate': replicate, 'steps': max_step(target), 'verified': True})
    (destination / 'MD-REPAIRS.json').write_text(json.dumps(repairs, indent=2) + '\n')
    results = []
    # Preserve smaller outputs first; keep FEP legs independently recoverable.
    for name in sorted(groups, key=lambda x: (x.startswith('fep-'), x)):
        results.append(archive(scratch, sorted(groups[name]), destination, name, resume=args.resume))
    (destination / 'COMPLETE.json').write_text(json.dumps({
        'completed_unix': time.time(), 'source_bytes': total,
        'archives': results, 'md_repairs': repairs}, indent=2) + '\n')
    print('PRESERVATION COMPLETE: ' + str(destination), flush=True)


if __name__ == '__main__':
    main()
