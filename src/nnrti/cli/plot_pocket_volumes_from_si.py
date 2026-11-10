"""Plot all genotypes directly from Supplementary Table 4, without reanalysis."""
from pathlib import Path
import argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--table', type=Path, default=Path('revisions/submission-ready/Supplementary-Table-4.xlsx'))
    p.add_argument('--output-dir', type=Path, default=Path('results/analysis/supplementary_figures'))
    args = p.parse_args()
    d = pd.read_excel(args.table, sheet_name='Per-replicate')
    col = 'NNIBP volume (Å³)'
    labels = list(dict.fromkeys(d.Genotype))
    if len(labels) != 20:
        raise ValueError(f'Expected 20 genotypes, found {len(labels)}')
    counts = d.groupby(['Genotype', 'State']).size()
    if len(counts) != 40 or not counts.eq(3).all() or d[col].isna().any():
        raise ValueError('Expected three finite replicate volumes for every genotype and state')
    fig, ax = plt.subplots(figsize=(9.0, 7.0))
    y = np.arange(len(labels))
    rows = []
    for state, offset, color, name in [('Bound (RT+DOR)', -.17, '#3774a5', 'Bound (holo)'),
                                      ('Unliganded (RT only)', .17, '#cf7834', 'Unliganded (apo)')]:
        for i, genotype in enumerate(labels):
            g = d[(d.Genotype == genotype) & (d.State == state)].sort_values('Replicate')
            v = g[col].to_numpy(float)
            mean, sem = v.mean(), v.std(ddof=1)/np.sqrt(len(v))
            ax.scatter(v, i + offset + np.array([-.055, 0, .055]), s=16, color=color, alpha=.55, zorder=3)
            ax.errorbar(mean, i+offset, xerr=sem, fmt='D', ms=4, color=color, capsize=3,
                        label=name if i == 0 else None, zorder=4)
            rows.append(dict(genotype=genotype, state=state, mean_A3=mean, sem_A3=sem, n=3))
    ax.set_yticks(y, labels, fontsize=9)
    ax.invert_yaxis()
    ax.set_xlim(left=0)
    ax.set_xlabel('NNIBP volume (Å³)', fontsize=11)
    ax.set_title('All 20 simulated genotypes', loc='left', fontsize=12, pad=30)
    ax.grid(axis='x', color='0.9', linewidth=.6)
    ax.set_axisbelow(True)
    ax.spines[['top', 'right']].set_visible(False)
    ax.legend(frameon=False, fontsize=9, loc='lower right', bbox_to_anchor=(1, 1.01), ncol=2)
    fig.tight_layout()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for ext in ['png', 'pdf']:
        fig.savefig(args.output_dir / f'supplementary_figure_6_pocket_volume.{ext}', dpi=300)
    pd.DataFrame(rows).to_csv(args.output_dir / 'supplementary_figure_6_pocket_volume_summary.csv', index=False)
    print(f'Plotted {len(labels)} genotypes, {len(d)} replicate values from {args.table}')


if __name__ == '__main__':
    main()
