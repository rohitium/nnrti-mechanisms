"""Export the submitted revision or rebuild plots/tables without overwriting it.

    python -m nnrti.cli.manuscript --check
    python -m nnrti.cli.manuscript --export
    python -m nnrti.cli.manuscript --rebuild heatmaps table4
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import zipfile
from datetime import datetime, timezone
from importlib.metadata import version

from lxml import etree
from nnrti.paths import REPO_ROOT

SUBMISSION = REPO_ROOT / 'revisions/submission-ready'
NS = {'w': 'http://schemas.openxmlformats.org/wordprocessingml/2006/main',
      'a': 'http://schemas.openxmlformats.org/drawingml/2006/main',
      'r': 'http://schemas.openxmlformats.org/officeDocument/2006/relationships'}


def check_submission() -> None:
    expected = json.loads((REPO_ROOT / 'paper/submission-checksums.json').read_text())
    failures = [name for name, digest in expected.items()
                if not (REPO_ROOT / name).is_file() or
                hashlib.sha256((REPO_ROOT / name).read_bytes()).hexdigest() != digest]
    if failures:
        raise RuntimeError('Submitted files differ from the frozen snapshot: ' + ', '.join(failures))
    print(f'Checked {len(expected)} submitted files: all SHA-256 hashes match.')


def export_submission(out: Path) -> None:
    """Export exact embedded artwork and main tables, including final captions."""
    from docx import Document
    figures = out / 'submitted-figures'
    tables = out / 'submitted-tables'
    figures.mkdir(parents=True, exist_ok=True)
    tables.mkdir(parents=True, exist_ok=True)
    records = []
    for filename, prefix in [('Manuscript-Clean.docx', 'Figure'),
                             ('Supporting-Information-Clean.docx', 'Supplementary-Figure')]:
        with zipfile.ZipFile(SUBMISSION / filename) as z:
            root = etree.fromstring(z.read('word/document.xml'))
            rels = {e.get('Id'): e.get('Target') for e in etree.fromstring(z.read('word/_rels/document.xml.rels'))}
            pending = None
            for p in root.findall('.//w:body/w:p', NS):
                for blip in p.findall('.//a:blip', NS):
                    pending = rels[blip.get('{' + NS['r'] + '}embed')]
                text = ''.join(p.xpath('.//w:t/text()', namespaces=NS))
                match = re.match(r'^(?:Supplementary )?Figure (\d+)[.:]', text)
                if match and pending:
                    target = figures / f'{prefix}-{match[1]}{Path(pending).suffix}'
                    target.write_bytes(z.read('word/' + pending))
                    records.append({'figure': f'{prefix}-{match[1]}', 'file': target.name,
                                    'document': filename, 'embedded_file': pending, 'caption': text,
                                    'sha256': hashlib.sha256(target.read_bytes()).hexdigest()})
                    pending = None
    if len(records) != 13:
        raise RuntimeError(f'Expected 6 main and 7 supplementary figures; found {len(records)}')
    (figures / 'captions.json').write_text(json.dumps(records, indent=2, ensure_ascii=False) + '\n')
    for i, table in enumerate(Document(SUBMISSION / 'Manuscript-Clean.docx').tables, 1):
        with (tables / f'Table-{i}.csv').open('w', newline='', encoding='utf-8-sig') as f:
            csv.writer(f).writerows([[c.text for c in row.cells] for row in table.rows])
    for source in sorted(SUBMISSION.glob('Supplementary-Table-*.xlsx')):
        shutil.copy2(source, tables / source.name)
    print(f'Exported exact submitted figures and tables to {out}')


def commands(out: Path) -> dict[str, list[list[str]]]:
    p = str(out)
    return {
        'heatmaps': [['nnrti.cli.plot_structural_heatmap', '--output-dir', p+'/heatmaps']],
        'mechanisms': [['nnrti.cli.plot_mechanism_panel', '--output', p+'/Figure-4.png']],
        'protein': [['nnrti.cli.plot_protein_level_figure', '--output-dir', p]],
        'dor-violins': [['nnrti.cli.plot_dor_contact_violins', '--output-dir', p+'/dor-violins']],
        'rt-violins': [['nnrti.cli.plot_rt_metric_violins', '--output-dir', p+'/rt-violins']],
        'convergence': [['nnrti.cli.plot_supplementary_figures', '--convergence-only', '--output-dir', p+'/convergence']],
        'work': [['nnrti.cli.plot_fep_work_distributions', '--output', p+'/Supplementary-Figure-4.png']],
        'pose': [['nnrti.cli.plot_pose_convergence', '--output', p+'/Supplementary-Figure-2.png']],
        'table2': [['nnrti.cli.build_table_2', '--output-csv', p+'/tables/Table-2-energetics.csv']],
        'table3': [['nnrti.cli.build_supplementary_table_3', '--output-xlsx', p+'/tables/Supplementary-Table-3.xlsx']],
        'table4': [['nnrti.cli.build_si_structural_tables', '--out-dir', p+'/tables']],
    }


def collect_artifacts(out: Path, rebuilt_targets=()) -> dict:
    """Collect the publication under its final labels and write an audit manifest."""
    registry = json.loads((REPO_ROOT/'paper/artifacts.json').read_text())
    collection = out/'paper'
    records = []
    lines = ['# Paper artifacts', '',
             'Exact submitted files are under `main/` and `supplementary/`.',
             'Numerically rebuilt candidates are under `rebuilt/`; curated artwork is preserved by exact export.', '',
             '| Artifact | Submitted | Rebuilt candidate |', '|---|---|---|']
    for artifact in registry['artifacts']:
        source = out/artifact['submitted']
        target = collection/artifact['collection']
        if not source.is_file():
            raise FileNotFoundError(source)
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, target)
        record = {**artifact, 'submitted_sha256': hashlib.sha256(target.read_bytes()).hexdigest()}
        rebuilt_link = 'Exact export / curated source'
        if artifact['rebuild_target'] in rebuilt_targets:
            rebuilt_source = out/artifact['rebuilt']
            if not rebuilt_source.is_file():
                raise FileNotFoundError(rebuilt_source)
            candidate = collection/'rebuilt'/artifact['collection']
            candidate.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(rebuilt_source, candidate)
            record['rebuilt_sha256'] = hashlib.sha256(candidate.read_bytes()).hexdigest()
            rebuilt_link = f'[Rebuilt](rebuilt/{artifact["collection"]})'
        lines.append(f'| {artifact["id"]} | [Submitted]({artifact["collection"]}) | {rebuilt_link} |')
        records.append(record)
    (collection/'README.md').write_text('\n'.join(lines)+'\n')
    provenance = {
        'created_utc': datetime.now(timezone.utc).isoformat(),
        'python': sys.version, 'python_executable': sys.executable,
        'user_site_disabled': bool(sys.flags.no_user_site), 'artifacts': records,
        'registry_sha256': hashlib.sha256((REPO_ROOT/'paper/artifacts.json').read_bytes()).hexdigest(),
        'rebuilt_targets': list(rebuilt_targets),
        'packages': {name: version(name) for name in ['numpy', 'pandas', 'scipy', 'matplotlib', 'openpyxl', 'python-docx', 'lxml']},
        'authoritative_inputs': {name: hashlib.sha256((REPO_ROOT/name).read_bytes()).hexdigest()
                                 for name in registry['authoritative_inputs']},
        'submission': json.loads((REPO_ROOT/'paper/submission-checksums.json').read_text()),
        'source_inventory': {str(path.relative_to(REPO_ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
                             for folder in ['src', 'workflows', 'tests', 'ops', 'manifests', 'results/analysis']
                             for path in sorted((REPO_ROOT/folder).rglob('*'))
                             if path.is_file() and path.suffix in {'.py', '.sh', '.csv', '.json', '.dat', '.png', '.pdf', '.xlsx'}},
    }
    (out/'reproduction.json').write_text(json.dumps(provenance, indent=2)+'\n')
    (out/'README.md').write_text(
        '# Reproduction outputs\n\n'
        '[Open the paper artifact collection](paper/README.md) for all 22 submitted figures and tables, '
        'plus the supported numerical rebuilds under their final publication labels.\n\n'
        '`reproduction.json` records source hashes, artifact hashes and package versions. '
        '`validation/` contains verification results when `--verify` is used. '
        'Other directories contain intermediate plotting outputs.\n')
    return provenance


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, default=REPO_ROOT/'paper/generated')
    parser.add_argument('--check', action='store_true')
    parser.add_argument('--export', action='store_true', help='Extract exact submitted figures and copy submitted tables.')
    parser.add_argument('--rebuild', nargs='+', metavar='TARGET', help='Targets from --list; use all for all supported rebuilds.')
    parser.add_argument('--list', action='store_true')
    parser.add_argument('--bundle', action='store_true', help='Export, rebuild all supported targets and collect under final paper labels.')
    parser.add_argument('--verify', action='store_true', help='Run repository tests and save validation logs (use with --bundle).')
    a = parser.parse_args(argv)
    if a.verify and not a.bundle:
        parser.error('--verify requires --bundle')
    if a.bundle:
        a.export = True
        a.rebuild = ['all']
    out = a.output_dir.resolve()
    # Prevent even a custom output argument from overwriting submission/history.
    for protected in [SUBMISSION.resolve(), (REPO_ROOT/'archive').resolve()]:
        if out == protected or protected in out.parents or out in protected.parents:
            raise ValueError(f'Output directory overlaps protected material: {out}')
    if a.bundle and not sys.flags.no_user_site:
        parser.error('Use workflows/05_manuscript_artifacts.sh or set PYTHONNOUSERSITE=1 before starting Python.')
    jobs = commands(out)
    names = []
    if a.bundle:
        out.mkdir(parents=True, exist_ok=True)
        # A failed invocation must not leave a previous successful report in place.
        (out/'reproduction.json').unlink(missing_ok=True)
        (out/'validation/PASSED.json').unlink(missing_ok=True)
    if a.list:
        print('\n'.join(jobs))
    if a.check or a.export or a.rebuild:
        check_submission()
    if a.export:
        export_submission(out)
    if a.rebuild:
        names = list(jobs) if a.rebuild == ['all'] else a.rebuild
        unknown = set(names) - jobs.keys()
        if unknown:
            parser.error('Unknown rebuild targets: ' + ', '.join(sorted(unknown)))
        out.mkdir(parents=True, exist_ok=True)
        (out/'tables').mkdir(exist_ok=True)
        cache = out/'.cache/matplotlib'
        cache.mkdir(parents=True, exist_ok=True)
        env = {**os.environ, 'PYTHONPATH': str(REPO_ROOT/'src'), 'MPLBACKEND': 'Agg',
               'MPLCONFIGDIR': str(cache), 'PYTHONNOUSERSITE': '1'}
        for name in names:
            for cmd in jobs[name]:
                print(f'Rebuilding {name}', flush=True)
                subprocess.run([sys.executable, '-m', *cmd], cwd=REPO_ROOT, env=env, check=True)
        check_submission()
    if a.bundle:
        collect_artifacts(out, names)
    if a.verify:
        validation = out/'validation'
        validation.mkdir(parents=True, exist_ok=True)
        cmd = [sys.executable, '-m', 'pytest', '-q', 'tests', '--junitxml', str(validation/'tests.xml')]
        with (validation/'tests.log').open('w') as log:
            result = subprocess.run(cmd, cwd=REPO_ROOT, env=env, stdout=log, stderr=subprocess.STDOUT)
        print((validation/'tests.log').read_text())
        result.check_returncode()
        check_submission()
        (validation/'PASSED.json').write_text(json.dumps({'tests_exit_code': 0,
            'submission_checksums': 'passed', 'command': cmd}, indent=2)+'\n')
    if not any([a.list, a.check, a.export, a.rebuild]):
        parser.print_help()
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
