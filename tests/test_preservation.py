"""Preservation must round-trip bytes and never certify a corrupt destination."""
import importlib.util
import json
from pathlib import Path
import tarfile

import pytest


@pytest.fixture
def preservation():
    path = Path(__file__).resolve().parents[1]/'ops/storage/preserve_sherlock.py'
    spec = importlib.util.spec_from_file_location('preservation', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_archive_round_trip_and_inventory(tmp_path, preservation):
    source, output = tmp_path/'source', tmp_path/'output'
    (source/'nested').mkdir(parents=True)
    output.mkdir()
    payload = bytes(range(256))*1000
    (source/'nested/trajectory.bin').write_bytes(payload)
    (source/'link').symlink_to('nested/trajectory.bin')
    report = preservation.archive(source, [Path('nested'), Path('nested/trajectory.bin'), Path('link')], output, 'snapshot')
    with tarfile.open(output/report['archive']) as archive:
        assert archive.extractfile('nested/trajectory.bin').read() == payload
        assert archive.getmember('link').linkname == 'nested/trajectory.bin'
    entries = [json.loads(line) for line in (output/report['inventory']).read_text().splitlines()]
    assert entries[1]['sha256'] == preservation.digest(source/'nested/trajectory.bin')
    assert report['sha256'] == preservation.digest(output/report['archive'])
    assert (output/'snapshot.verified.json').is_file()


def test_checksum_failure_never_certifies_archive(tmp_path, preservation, monkeypatch):
    source, output = tmp_path/'source', tmp_path/'output'
    source.mkdir()
    output.mkdir()
    (source/'data').write_bytes(b'irreplaceable input')
    monkeypatch.setattr(preservation, 'digest', lambda path: 'corrupt')
    with pytest.raises(RuntimeError, match='verification failed'):
        preservation.archive(source, [Path('data')], output, 'snapshot')
    assert (output/'snapshot.tar.part').exists()
    assert not (output/'snapshot.tar').exists()
    assert not (output/'snapshot.verified.json').exists()


def test_resume_rejects_source_changes(tmp_path, preservation):
    source, output = tmp_path/'source', tmp_path/'output'
    source.mkdir()
    output.mkdir()
    (source/'data').write_bytes(b'original')
    report = preservation.archive(source, [Path('data')], output, 'snapshot')
    assert preservation.archive(source, [Path('data')], output, 'snapshot', resume=True) == report
    (source/'data').write_bytes(b'a different source')
    with pytest.raises(RuntimeError, match='Source changed'):
        preservation.archive(source, [Path('data')], output, 'snapshot', resume=True)
    assert preservation.digest(output/'snapshot.tar') == report['sha256']


def test_finalization_preserves_bytes_and_legacy_paths(tmp_path, preservation, monkeypatch):
    import subprocess
    import sys
    scratch, group = tmp_path/'scratch', tmp_path/'group'
    project = group/'nnrti-mechanisms'
    snapshot = project/'archive/snapshot'
    incoming = project/'archive/local.incoming'
    for directory in [scratch, snapshot, incoming, project/'operations']:
        directory.mkdir(parents=True, exist_ok=True)
    (project/'operations/project-state.json').write_text('{}')
    for directory in ['src', 'ops', 'workflows', 'data', 'manifests']:
        (project/'code'/directory).mkdir(parents=True)
    for name in ['nnrti-mechanisms-git/results/trace.bin',
                 'nnrti-mechanisms-git/src/old.py',
                 'nnrti-mechanisms/results/staging.bin',
                 'apo_collapse/coordinates.bin', 'migration.log']:
        path = scratch/name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(b'preserve these exact bytes')
    apo = scratch/'nnrti-mechanisms-git/results/md_runs/apo/wt/rep_01'
    apo.mkdir(parents=True)
    (apo/'wt_md_state.csv').write_text('Step,Energy\n50000000,0\n')
    (apo/'trajectory.bin').write_bytes(b'completed trajectory')
    mirror = project/'md_runs/apo/wt/rep_01'
    mirror.mkdir(parents=True)
    (mirror/'wt_md_state.csv').write_text('Step,Energy\n5000000,0\n')
    paths = sorted(path.relative_to(scratch) for path in scratch.rglob('*'))
    record = preservation.archive(scratch, paths, snapshot, 'all')
    (snapshot/'COMPLETE.json').write_text(json.dumps({'archives': [record]}))
    for name in ['analysis', 'history', 'inputs', 'md-runs']:
        preservation.archive(scratch, [Path('migration.log')], incoming, name)
    monkeypatch.setenv('SCRATCH', str(scratch))
    monkeypatch.setenv('GROUP_HOME', str(group))
    monkeypatch.setenv('SLURM_JOB_ID', '123')
    monkeypatch.setattr(subprocess, 'check_output', lambda command: b'123\n')
    monkeypatch.syspath_prepend(str(Path(__file__).resolve().parents[1]/'ops/storage'))
    path = Path(__file__).resolve().parents[1]/'ops/storage/finalize_sherlock.py'
    spec = importlib.util.spec_from_file_location('finalization', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    monkeypatch.setattr(sys, 'argv', [str(path), '--snapshot', 'snapshot', '--local-transfer', 'local.incoming'])
    module.main()
    work = scratch/'nnrti-mechanisms-work'
    assert (work/'results/trace.bin').read_bytes() == b'preserve these exact bytes'
    assert (work/'legacy/main-checkout/src/old.py').is_file()
    assert (work/'legacy/auxiliary/apo_collapse/coordinates.bin').is_file()
    assert (work/'legacy/auxiliary/migration.log').is_file()
    assert (scratch/'nnrti-mechanisms-git').resolve() == work
    assert (scratch/'nnrti-mechanisms').resolve() == work/'legacy/staging'
    assert (work/'src').resolve() == project/'code/src'
    assert (project/'archive/local/COMPLETE.json').is_file()
    assert (mirror/'trajectory.bin').read_bytes() == b'completed trajectory'
    assert preservation.max_step(mirror) == 50000000
    assert preservation.max_step(project/'archive/md-mirror-prior/snapshot/wt/rep_01') == 5000000
    assert (project/'operations/MD-MIRROR-VERIFIED.json').is_file()
