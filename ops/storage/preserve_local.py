#!/usr/bin/env python3
"""Preserve local simulation inputs, derived data and historical recovery files.

Usage: python ops/storage/preserve_local.py --snapshot UNIQUE_NAME
Code and the frozen submission are versioned separately in Git. Local data
archives stay ignored under archive/local-snapshots/ and may be transferred to
the durable Sherlock project archive. Never upload a completion marker before
the destination archives have been independently verified.
"""
import argparse
import json
import os
from pathlib import Path

from preserve_sherlock import archive


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--snapshot', required=True)
    args = parser.parse_args()
    if not args.snapshot or args.snapshot in {'.', '..'} or Path(args.snapshot).name != args.snapshot:
        raise ValueError('Snapshot must be a single directory name.')
    root = Path(__file__).resolve().parents[2]
    destination = root/'archive/local-snapshots'/args.snapshot
    destination.mkdir(parents=True, exist_ok=False)
    groups = {'md-runs': ['results/md_runs', 'results/.checkpoints'],
              'history': ['archive'], 'inputs': ['data', 'manifests'],
              'analysis': ['results/analysis']}
    reports = []
    for name, folders in groups.items():
        paths = []
        for folder in folders:
            start = root/folder
            if not start.exists():
                continue
            paths.append(start.relative_to(root))
            for current, dirs, files in os.walk(str(start), followlinks=False):
                # Exclude every prior preservation output to avoid nested backups.
                dirs[:] = [d for d in dirs if d != '__pycache__' and
                           not (Path(current) == root/'archive/local-snapshots' and
                                (Path(current)/d/'COMPLETE.json').is_file()) and
                           Path(current)/d != destination]
                paths.extend((Path(current)/item).relative_to(root)
                             for item in sorted(dirs + files) if item != '.DS_Store')
        reports.append(archive(root, sorted(paths), destination, name))
    (destination/'COMPLETE.json').write_text(json.dumps({
        'source': str(root), 'archives': reports}, indent=2)+'\n')
    print('LOCAL PRESERVATION COMPLETE: ' + str(destination), flush=True)


if __name__ == '__main__':
    main()
