"""Regression guards for the frozen submission and regenerated statistics."""
from pathlib import Path
import subprocess
import sys

import numpy as np
import pandas as pd
import pytest
from openpyxl import load_workbook
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
SUBMITTED = ROOT/'revisions/submission-ready/Supplementary-Table-4.xlsx'
SOURCE = ROOT/'results/analysis/structural_metrics/tables/structural_metrics_per_replicate.csv'


def validate_statistics(path):
    """Independently recompute every comparison from the authoritative CSV."""
    from nnrti.cli.plot_structural_heatmap import CIRCULAR, DOR_METRICS, RT_METRICS
    metrics = {**DOR_METRICS, **RT_METRICS}
    data = pd.read_csv(SOURCE)
    values = {}
    for (gt, state), group in data.groupby(['genotype', 'structure_type']):
        for label, col in metrics.items():
            v = group[col].dropna().to_numpy()
            if not len(v):
                continue
            if label in CIRCULAR:
                z = np.exp(1j*np.deg2rad(v)).mean()
                mean = np.rad2deg(np.angle(z))
                sem = np.rad2deg(np.sqrt(-2*np.log(np.clip(abs(z), 1e-12, 1))))/np.sqrt(len(v))
            else:
                mean, sem = v.mean(), v.std(ddof=1)/np.sqrt(len(v))
            values[(gt, state, label)] = (mean, sem, len(v))
    wb = load_workbook(path, data_only=True)
    count = 0
    for name in ['Mutant vs WT', 'Induced fit vs WT']:
        interaction = name.startswith('Induced')
        for row in wb[name].iter_rows(min_row=2, values_only=True):
            gt, label = row[1], row[2 if interaction else 3]
            if interaction:
                keys = [(gt,'holo',label),(gt,'apo',label),('WT','holo',label),('WT','apo',label)]
            else:
                state = 'holo' if row[2].startswith('Bound') else 'apo'
                keys = [(gt,state,label),('WT',state,label)]
            v = [values[k] for k in keys]
            shift = v[0][0]-v[1][0]
            if interaction:
                shift -= v[2][0]-v[3][0]
            if label in CIRCULAR:
                shift = (shift+180)%360-180
            sem = np.sqrt(sum(x[1]**2 for x in v))
            df = sum(x[2]-1 for x in v)
            t = shift/sem if sem > 0 else None
            p = 2*stats.t.sf(abs(t),df) if t is not None else None
            half = stats.t.ppf(.975,df)*sem
            expected = [shift,sem,t,df,p,shift-half,shift+half,'yes' if abs(shift)>half else 'no']
            actual = row[3 if interaction else 4:]
            for a, e in zip(actual, expected):
                if e is None:
                    assert a in (None,'')
                elif isinstance(e,str):
                    assert a == e, (name,gt,label)
                else:
                    assert a == pytest.approx(e,rel=1e-8,abs=1e-10), (name,gt,label,a,e)
            count += 1
    assert count == 753


def test_submitted_statistics_match_source():
    validate_statistics(SUBMITTED)


def test_submission_checksums():
    from nnrti.cli.manuscript import check_submission
    check_submission()


def test_regenerated_workbook_is_live_and_full_precision(tmp_path):
    from nnrti.cli.build_si_structural_tables import main
    from nnrti.cli.plot_structural_heatmap import DOR_METRICS, RT_METRICS
    main(['--out-dir', str(tmp_path)])
    wb = load_workbook(tmp_path/'Supplementary-Table-4.xlsx',data_only=False)
    assert wb.sheetnames == ['Summary','Per-replicate','Mutant vs WT','Induced fit vs WT','Group statistics']
    assert sum(c.data_type == 'f' for sh in wb for row in sh for c in row) == 10732
    source = pd.read_csv(SOURCE).set_index(['genotype','structure_type','replicate'])
    metrics = {**DOR_METRICS, **RT_METRICS}
    for row in wb['Per-replicate'].iter_rows(min_row=2,values_only=True):
        key = (row[1], 'holo' if row[2].startswith('Bound') else 'apo',row[3])
        for actual,col in zip(row[4:],metrics.values()):
            expected = source.loc[key,col]
            if pd.isna(expected):
                assert actual is None
            else:
                assert actual == pytest.approx(expected, rel=1e-14)
    assert 'MOD(' in wb['Mutant vs WT']['E2'].value or any(
        'MOD(' in c.value for row in wb['Mutant vs WT'] for c in row if c.data_type=='f')


def test_rebuild_rejects_submission_output():
    from nnrti.cli.manuscript import main, SUBMISSION
    with pytest.raises(ValueError,match='protected'):
        main(['--rebuild','table4','--output-dir',str(SUBMISSION)])


def test_export_has_exact_current_figures(tmp_path):
    from nnrti.cli.manuscript import export_submission
    export_submission(tmp_path)
    assert len(list((tmp_path/'submitted-figures').glob('Figure-*.png'))) == 6
    assert len(list((tmp_path/'submitted-figures').glob('Supplementary-Figure-*.png'))) == 7
    assert len(list((tmp_path/'submitted-tables').glob('Supplementary-Table-*.xlsx'))) == 7


def test_native_spreadsheet_recalculation(tmp_path):
    """Exercise real formula evaluation when a spreadsheet engine is available."""
    import os
    import shutil
    engine = os.environ.get('NNRTI_SOFFICE') or shutil.which('soffice')
    if not engine:
        pytest.skip('Set NNRTI_SOFFICE or install LibreOffice for native formula evaluation')
    from nnrti.cli.build_si_structural_tables import main
    source, out = tmp_path/'source', tmp_path/'recalculated'
    main(['--out-dir',str(source)])
    subprocess.run([engine, '-env:UserInstallation='+ (tmp_path/'profile').as_uri(),
                    '--headless','--convert-to','xlsx','--outdir',str(out),
                    str(source/'Supplementary-Table-4.xlsx')],check=True,capture_output=True,timeout=60)
    result = out/'Supplementary-Table-4.xlsx'
    validate_statistics(result)
    values = load_workbook(result,data_only=True)
    assert not [(s.title,c.coordinate,c.value) for s in values for row in s for c in row if c.data_type=='e']
