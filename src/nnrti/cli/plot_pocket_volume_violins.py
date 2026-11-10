"""Rebuild Supplementary Figure 6 from saved per-frame volumes, excluding F227C."""
from pathlib import Path
import argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, default=Path('results/analysis/supplementary_figures'))
    args = parser.parse_args()
    table = pd.read_excel('revisions/submission-ready/Supplementary-Table-4.xlsx', sheet_name='Per-replicate')
    labels = [g for g in dict.fromkeys(table.Genotype) if g != 'F227C']
    holo = pd.read_csv('results/.checkpoints/.checkpoint_pocket_volume_profiles.csv')
    holo = holo.rename(columns={'pocket_volume_proxy_angstrom3': 'volume', 'frame_index': 'frame'})
    apo = pd.read_csv('results/analysis/apo_ensembles/tables/apo_frame_metrics.csv')
    apo = apo[(apo.window == 'matched') & apo.pocket_intact].copy()
    # Match Table 4's duration filter; retain its recorded G190E frame exclusion.
    window_ns = apo.time_ns.max()
    apo = apo[apo.total_ns >= window_ns * .98]
    apo['mutation'] = apo.safe_label.str.upper().str.replace('_', '+', regex=False)
    apo = apo.rename(columns={'pocket_volume': 'volume'})
    fig, ax = plt.subplots(figsize=(11.5, 4.6))
    pos = np.arange(len(labels))
    exports = []
    for data, offset, color, name, state in [
        (holo, -.19, '#4C72B0', 'Bound (holo)', 'holo'),
        (apo, .19, '#DD8452', 'Unliganded (apo)', 'apo')]:
        data = data[data.mutation.isin(labels)].copy()
        # Preserve the previous Figure 6 grid-saturation screen.
        excluded = int((data.volume > 500).sum())
        print(f'{state}: excluded {excluded} saved grid-saturation values above 500 Å³')
        data = data[data.volume <= 500]
        if set(data.mutation) != set(labels) or not np.isfinite(data.volume).all():
            raise ValueError(f'Incomplete or nonfinite {state} volumes')
        if not data.groupby('mutation').replicate.nunique().eq(3).all():
            raise ValueError(f'Expected three replicates per genotype in {state}')
        values = [data.loc[data.mutation == label, 'volume'].to_numpy() for label in labels]
        parts = ax.violinplot(values, positions=pos+offset, widths=.34,
                              showmedians=True, showextrema=False)
        for body in parts['bodies']:
            body.set_facecolor(color); body.set_alpha(.75); body.set_edgecolor('none')
        parts['cmedians'].set_color('0.2'); parts['cmedians'].set_linewidth(.9)
        ax.plot([], [], color=color, lw=5, alpha=.75, label=name)
        data['state'] = state
        exports.append(data[['mutation', 'state', 'replicate', 'frame', 'volume']])
    ax.set_xticks(pos, labels, rotation=45, ha='right', fontsize=8)
    ax.set_ylabel('NNIBP volume (Å³)')
    ax.set_ylim(bottom=0)
    ax.legend(frameon=False, loc='upper left', ncol=2)
    ax.grid(axis='y', alpha=.2, lw=.5)
    ax.set_axisbelow(True)
    ax.spines[['top', 'right']].set_visible(False)
    fig.tight_layout()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for ext in ['png', 'pdf']:
        fig.savefig(args.output_dir / f'supplementary_figure_6_pocket_volume.{ext}', dpi=300)
    pd.concat(exports).to_csv(args.output_dir/'supplementary_figure_6_pocket_volume_frames.csv', index=False)
    print(f'Plotted {len(labels)} genotypes, excluding standalone F227C; no volumes recalculated')


if __name__ == '__main__':
    main()
