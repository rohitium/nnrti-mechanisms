#!/usr/bin/env python3
"""Verify transferred local archives and consolidate preserved scratch sources.

Run in Slurm after preservation and local rsync finish. No scientific file is
deleted: scratch consolidation uses same-filesystem renames and retains two
compatibility links. Unexpected source changes stop consolidation.
"""
import argparse
import json
import os
from pathlib import Path
import subprocess

from preserve_sherlock import digest, max_step


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--snapshot', required=True)
    parser.add_argument('--local-transfer', required=True)
    args = parser.parse_args()
    if not os.environ.get('SLURM_JOB_ID'):
        raise RuntimeError('Run final verification in a Slurm allocation.')
    for value in [args.snapshot, args.local_transfer]:
        if not value or value in {'.', '..'} or Path(value).name != value:
            raise ValueError('Use single directory names, not paths.')
    scratch = Path(os.environ['SCRATCH'])
    project = Path(os.environ['GROUP_HOME'])/'nnrti-mechanisms'
    snapshot = project/'archive'/args.snapshot
    complete = json.loads((snapshot/'COMPLETE.json').read_text())
    expected_paths = {}
    for record in complete['archives']:
        verified = json.loads((snapshot/(record['archive'][:-4]+'.verified.json')).read_text())
        if verified != record or digest(snapshot/record['inventory']) != record['inventory_sha256']:
            raise RuntimeError('Snapshot verification record mismatch.')
        for line in (snapshot/record['inventory']).read_text().splitlines():
            item = json.loads(line)
            expected_paths[item['path']] = item
    incoming = project/'archive'/args.local_transfer
    if not args.local_transfer.endswith('.incoming'):
        raise RuntimeError('Expected a new .incoming transfer directory.')
    final = project/'archive'/args.local_transfer[:-len('.incoming')]
    if final.exists():
        raise RuntimeError('Expected a new .incoming transfer directory.')
    local_records = []
    for name in ['analysis', 'history', 'inputs', 'md-runs']:
        record = json.loads((incoming/(name+'.verified.json')).read_text())
        path = incoming/record['archive']
        if path.stat().st_size != record['bytes'] or digest(path) != record['sha256']:
            raise RuntimeError('Transferred archive differs: ' + str(path))
        if digest(incoming/record['inventory']) != record['inventory_sha256']:
            raise RuntimeError('Transferred inventory differs: ' + name)
        local_records.append(record)
        print('Verified local transfer: '+name, flush=True)
    incoming.rename(final)
    (final/'COMPLETE.json').write_text(json.dumps({'archives': local_records,
        'verified_on': os.uname().nodename, 'job_id': os.environ['SLURM_JOB_ID']}, indent=2)+'\n')
    # Recheck after hashing, so no simulation writer is overlooked during the wait.
    jobs = subprocess.check_output(['squeue', '-h', '-u', os.environ['USER'], '-o', '%A']).decode().split()
    if set(jobs)-{os.environ['SLURM_JOB_ID']}:
        raise RuntimeError('Other jobs exist; leave scratch unchanged.')
    apo_source = scratch/'nnrti-mechanisms-git/results/md_runs/apo'
    apo_mirror = project/'md_runs/apo'
    if apo_source.exists():
        for replicate in apo_source.glob('*/rep_*'):
            target = apo_mirror/replicate.relative_to(apo_source)
            if max_step(replicate) < 50000000 or max_step(target) > max_step(replicate):
                raise RuntimeError('MD mirror would replace a newer run or use an incomplete source.')
        prior_mirror = project/'archive/md-mirror-prior'/args.snapshot
        subprocess.run(['rsync', '-a', '--checksum', '--backup', '--backup-dir', str(prior_mirror),
                        str(apo_source)+'/', str(apo_mirror)+'/'], check=True)
        checked = 0
        prefix = 'nnrti-mechanisms-git/results/md_runs/apo/'
        for name, item in expected_paths.items():
            if name.startswith(prefix) and item['type'] == '0':
                target = apo_mirror/name[len(prefix):]
                if digest(target) != item['sha256']:
                    raise RuntimeError('MD mirror differs from preserved source: '+str(target))
                checked += 1
        (project/'operations/MD-MIRROR-VERIFIED.json').write_text(json.dumps({
            'verified_files': checked, 'snapshot': args.snapshot,
            'replicates': {str(path.relative_to(apo_mirror)): max_step(path)
                           for path in sorted(apo_mirror.glob('*/rep_*'))}}, indent=2)+'\n')
        print('Verified MD mirror files: '+str(checked), flush=True)
    roots = sorted({Path(name).parts[0] for name in expected_paths})
    actual_paths = set()
    for name in roots:
        start = scratch/name
        actual_paths.add(name)
        if start.is_dir() and not start.is_symlink():
            for folder, dirs, files in os.walk(str(start), followlinks=False):
                actual_paths.update(str((Path(folder)/item).relative_to(scratch)) for item in dirs+files)
    if actual_paths != set(expected_paths):
        raise RuntimeError('Scratch membership changed; preserve new/missing files before consolidation.')
    for name, item in expected_paths.items():
        path = scratch/name
        info = path.lstat()
        if item['type'] == '0' and (info.st_size, info.st_mtime_ns) != (item['bytes'], item['mtime_ns']):
            raise RuntimeError('Scratch file changed since preservation: '+name)
        if item['type'] == '2' and os.readlink(str(path)) != item['link_target']:
            raise RuntimeError('Scratch link changed since preservation: '+name)
    workspace = scratch/'nnrti-mechanisms-work'
    if workspace.exists():
        raise RuntimeError('Workspace already exists; refusing to merge implicitly.')
    retired = workspace/'legacy/main-checkout'
    retired.mkdir(parents=True)
    main_tree = scratch/'nnrti-mechanisms-git'
    (main_tree/'results').rename(workspace/'results')
    if (main_tree/'logs').exists():
        (main_tree/'logs').rename(workspace/'logs')
    for item in main_tree.iterdir():
        item.rename(retired/item.name)
    main_tree.rmdir()
    main_tree.symlink_to(workspace)
    staging = scratch/'nnrti-mechanisms'
    staging.rename(workspace/'legacy/staging')
    staging.symlink_to(workspace/'legacy/staging')
    auxiliary = workspace/'legacy/auxiliary'
    auxiliary.mkdir()
    for name in roots:
        if name not in ['nnrti-mechanisms', 'nnrti-mechanisms-git']:
            (scratch/name).rename(auxiliary/name)
    for name in ['src', 'ops', 'workflows', 'data', 'manifests']:
        (workspace/name).symlink_to(project/'code'/name)
    (workspace/'README.md').write_text(
        '# Disposable computation workspace\n\n'
        'Durable project: $GROUP_HOME/nnrti-mechanisms/. Read its README first.\n'
        'Code and input links point to the maintained code/ checkout. Results and logs stay on scratch.\n'
        'Set PROJECT_ROOT to this directory before submitting simulation jobs.\n'
        'legacy/ contains preserved historical working copies, not alternate production checkouts.\n'
        'Every scratch file remains subject to the 90-day content-modification purge policy.\n')
    state_path = project/'operations/project-state.json'
    state = json.loads(state_path.read_text())
    state.update({'scratch_workspace': str(workspace), 'scratch_consolidated': True,
                  'preservation_complete': True, 'local_snapshot': str(final),
                  'local_transfer_verified': True, 'sources_retired': True,
                  'sources_deleted': False})
    state_path.write_text(json.dumps(state, indent=2)+'\n')
    (project/'operations/CONSOLIDATED.json').write_text(json.dumps(state, indent=2)+'\n')
    print('PRESERVATION AND CONSOLIDATION COMPLETE', flush=True)


if __name__ == '__main__':
    main()
