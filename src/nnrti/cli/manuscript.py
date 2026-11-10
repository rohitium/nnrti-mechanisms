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


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, default=REPO_ROOT/'paper/generated')
    parser.add_argument('--check', action='store_true')
    parser.add_argument('--export', action='store_true', help='Extract exact submitted figures and copy submitted tables.')
    parser.add_argument('--rebuild', nargs='+', metavar='TARGET', help='Targets from --list; use all for all supported rebuilds.')
    parser.add_argument('--list', action='store_true')
    a = parser.parse_args(argv)
    out = a.output_dir.resolve()
    # Prevent even a custom output argument from overwriting submission/history.
    for protected in [SUBMISSION.resolve(), (REPO_ROOT/'archive').resolve()]:
        if out == protected or protected in out.parents or out in protected.parents:
            raise ValueError(f'Output directory overlaps protected material: {out}')
    jobs = commands(out)
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
        env = {**os.environ, 'PYTHONPATH': str(REPO_ROOT/'src'), 'MPLBACKEND': 'Agg'}
        for name in names:
            for cmd in jobs[name]:
                print(f'Rebuilding {name}', flush=True)
                subprocess.run([sys.executable, '-m', *cmd], cwd=REPO_ROOT, env=env, check=True)
        check_submission()
    if not any([a.list, a.check, a.export, a.rebuild]):
        parser.print_help()
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
