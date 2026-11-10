"""Live Excel formulas for Supplementary Table 4.

Replicate inputs retain full precision. Formula caches require Excel or
LibreOffice recalculation; openpyxl does not evaluate formulas.
"""
from pathlib import Path

from openpyxl import load_workbook
from openpyxl.styles import Alignment, Font
from openpyxl.utils import get_column_letter as letter
from openpyxl.workbook.properties import CalcProperties

from nnrti.cli.plot_structural_heatmap import CIRCULAR, DOR_METRICS, RT_METRICS


def add_formulas(path: Path) -> None:
    """Replace static derived cells with formulas linked to Per-replicate."""
    from nnrti.cli.build_si_structural_tables import _dec

    wb = load_workbook(path)
    detail, summary = wb['Per-replicate'], wb['Summary']
    metrics = {**DOR_METRICS, **RT_METRICS}
    headers = [c.value for c in detail[1]]
    if headers[4:] != list(metrics):
        raise ValueError('Structural workbook columns do not match metric definitions')
    inputs = {}
    for row in range(2, detail.max_row + 1):
        key = (detail.cell(row, 2).value, detail.cell(row, 3).value)
        inputs.setdefault(key, []).append(row)
        for col, label in enumerate(headers[4:], 5):
            detail.cell(row, col).number_format = '0.' + '0' * _dec(metrics[label])

    group = wb.create_sheet('Group statistics')
    group.append(['Category', 'Genotype', 'State', 'Metric', 'n', 'Mean', 'SEM',
                  'Mean sine (χ1)', 'Mean cosine (χ1)', 'Resultant length (χ1)'])
    groups = {}
    for srow in range(2, summary.max_row + 1):
        gt, state = summary.cell(srow, 2).value, summary.cell(srow, 3).value
        rows = inputs[(gt, state)]
        if rows != list(range(min(rows), max(rows) + 1)):
            raise ValueError(f'Noncontiguous replicate inputs for {gt}, {state}')
        for col, label in enumerate(headers[4:], 5):
            if not any(isinstance(detail.cell(r, col).value, (int, float)) for r in rows):
                continue  # Undefined metrics remain blank, not zero.
            r = group.max_row + 1
            groups[(gt, state, label)] = r
            ref = lambda c, rr: f"'Per-replicate'!{letter(c)}{rr}"
            rng = f'{ref(col, min(rows))}:{letter(col)}{max(rows)}'
            values = [f'={ref(c, rows[0])}' for c in (1, 2, 3)]
            values += [f"='Per-replicate'!{letter(col)}$1", f'=COUNT({rng})']
            if label in CIRCULAR:
                sincos = [
                    '=(' + '+'.join(f'IF(ISNUMBER({ref(col, rr)}),{fn}(RADIANS({ref(col, rr)})),0)'
                                     for rr in rows) + f')/E{r}'
                    for fn in ('SIN', 'COS')]
                mean = (f'=IF(I{r}=0,IF(H{r}>0,90,IF(H{r}<0,-90,0)),'
                        f'DEGREES(ATAN(H{r}/I{r}))+IF(I{r}<0,IF(H{r}>=0,180,-180),0))')
                sem = f'=IF(E{r}>1,DEGREES(SQRT(-2*LN(MIN(1,MAX(J{r},1E-12)))))/SQRT(E{r}),"")'
                values += [mean, sem, *sincos, f'=SQRT(H{r}^2+I{r}^2)']
            else:
                values += [f'=AVERAGE({rng})', f'=IF(E{r}>1,STDEV({rng})/SQRT(E{r}),"")', '', '', '']
            group.append(values)
        summary.cell(srow, 4, f"=COUNT('Per-replicate'!D{min(rows)}:D{max(rows)})")
        for col, label in enumerate(headers[4:], 5):
            r = groups.get((gt, state, label))
            if r is None:
                summary.cell(srow, col).value = None
                continue
            mean, sem = f"'Group statistics'!F{r}", f"'Group statistics'!G{r}"
            dec = _dec(metrics[label])
            sem_format = '"0.' + '0'*(dec+4) + '"'
            for d in range(dec+3, dec-1, -1):
                sem_format = f'IF(AND({sem}>0,ROUND({sem},{d})=0),{sem_format},"0.{"0"*d}")'
            summary.cell(srow, col, f'=TEXT({mean},"0.{"0"*dec}")&" ± "&TEXT({sem},{sem_format})')

    bound, apo = 'Bound (RT+DOR)', 'Unliganded (RT only)'
    for name in ('Mutant vs WT', 'Induced fit vs WT'):
        sheet = wb[name]
        interaction = name == 'Induced fit vs WT'
        first = 4 if interaction else 5
        sheet.cell(1, first + 7, '95% CI excludes zero')
        for row in range(2, sheet.max_row + 1):
            gt = sheet.cell(row, 2).value
            label = sheet.cell(row, 3 if interaction else 4).value
            keys = ([(gt, bound, label), (gt, apo, label), ('WT', bound, label), ('WT', apo, label)]
                    if interaction else [(gt, sheet.cell(row, 3).value, label),
                                         ('WT', sheet.cell(row, 3).value, label)])
            rs = [groups[k] for k in keys]
            means, sems, counts = [[f"'Group statistics'!{col}{r}" for r in rs] for col in 'FGE']
            expr = (f'({means[0]}-{means[1]})-({means[2]}-{means[3]})' if interaction
                    else f'{means[0]}-{means[1]}')
            if label in CIRCULAR:
                expr = f'MOD(({expr})+180,360)-180'
            shift = f'=IF(COUNT({",".join(means)})={len(rs)},{expr},"")'
            sem = f'=IF(COUNT({",".join(sems)})={len(rs)},SQRT({"+".join(s+"^2" for s in sems)}),"")'
            sc, ec, tc, dc, pc, lc, hc, fc = [f'{letter(c)}{row}' for c in range(first, first+8)]
            lower = f'=IF(AND(ISNUMBER({sc}),ISNUMBER({ec})),{sc}-_xlfn.T.INV.2T(0.05,{dc})*{ec},"")'
            formulas = [shift, sem,
                        f'=IF(AND(ISNUMBER({sc}),ISNUMBER({ec}),{ec}>0),{sc}/{ec},"")',
                        '=' + '+'.join(counts) + f'-{len(rs)}',
                        f'=IF(ISNUMBER({tc}),_xlfn.T.DIST.2T(ABS({tc}),{dc}),"")',
                        lower, lower.replace('-_xlfn', '+_xlfn'),
                        f'=IF(AND(ISNUMBER({lc}),ISNUMBER({hc})),IF(OR({lc}>0,{hc}<0),"yes","no"),"")']
            formats = ['0.000', '0.000', '0.00', '0', '0.00E+00', '0.000', '0.000', 'General']
            for offset, (formula, fmt) in enumerate(zip(formulas, formats)):
                sheet.cell(row, first+offset, formula).number_format = fmt
    group.freeze_panes = 'E2'
    group.auto_filter.ref = group.dimensions
    group.sheet_view.showGridLines = False
    for cell in group[1]:
        cell.font = Font(name='Calibri', size=11, bold=True)
        cell.alignment = Alignment(wrap_text=True)
    for col, width in {'A':22, 'B':23, 'C':25, 'D':72, 'E':8, 'F':16, 'G':16, 'H':21, 'I':21, 'J':21}.items():
        group.column_dimensions[col].width = width
    for row in group.iter_rows(min_row=2, min_col=6):
        for cell in row:
            cell.number_format = '0.000000'
    wb.calculation = CalcProperties(calcId=0, fullCalcOnLoad=True)
    wb.save(path)
