"""Apo/holo replicate distributions for all ten Figure 5 metrics from SI Table 4."""
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import pandas as pd
from nnrti.cli.plot_structural_heatmap import RT_METRICS

OUT = Path('results/analysis/supplementary_figures')


def main():
    data = pd.read_excel('revisions/submission-ready/Supplementary-Table-4.xlsx', sheet_name='Per-replicate')
    labels = [g for g in dict.fromkeys(data.Genotype) if g != 'F227C']
    metrics = list(RT_METRICS)
    assert len(labels) == 19 and len(metrics) == 10
    OUT.mkdir(parents=True, exist_ok=True)
    exports = []
    with PdfPages(OUT/'supplementary_figure_6_all_rt_metrics.pdf') as pdf:
        for page in range(5):
            fig, axes = plt.subplots(2, 1, figsize=(11.5, 8.4))
            for index, ax in enumerate(axes, page*2):
                metric = metrics[index]
                circular = 'χ1' in metric
                for state, offset, color, name in [('Bound (RT+DOR)', -.19, '#4C72B0', 'Bound (holo)'),
                                                 ('Unliganded (RT only)', .19, '#DD8452', 'Unliganded (apo)')]:
                    for i, label in enumerate(labels):
                        group = data[(data.Genotype == label) & (data.State == state)].sort_values('Replicate')
                        vals = group[metric].dropna().to_numpy(float)
                        absent = ('Tyr181' in metric and label == 'Y181C') or ('Tyr188' in metric and label == 'Y188L')
                        if absent:
                            assert len(vals) == 0
                            continue
                        assert len(vals) == 3 and np.isfinite(vals).all(), (metric, label, state)
                        x = i + offset
                        if circular:
                            # Three circular means do not support a density estimate.
                            # Display each observed replicate angle without a KDE.
                            vals = (vals+180)%360-180
                        elif np.ptp(vals) > 0:
                            violin = ax.violinplot([vals], positions=[x], widths=.32, showextrema=False)
                            for body in violin['bodies']:
                                body.set_facecolor(color);body.set_edgecolor('none');body.set_alpha(.65)
                        else:
                            ax.plot([x-.12,x+.12],[vals[0]]*2,color=color,lw=2)
                        ax.scatter(x+np.array([-.035,0,.035]), vals, s=10, c=color, edgecolors='0.25', linewidths=.35, zorder=4)
                        for rep, val in zip(group.Replicate, vals):
                            exports.append(dict(metric=metric, genotype=label, state=state, replicate=rep, value=val))
                    ax.plot([],[],color=color,lw=5,alpha=.65,label=name)
                ax.set_title(f'{chr(65+index)}  {metric}', loc='left', fontsize=11)
                ax.set_xticks(range(len(labels)),labels,rotation=45,ha='right',fontsize=8)
                ax.set_xlim(-.7,len(labels)-.3)
                if circular:
                    ax.set_ylim(-190,190);ax.set_yticks([-180,-90,0,90,180])
                    ax.text(.99, 1.04, 'Replicate circular means; −180° ≡ +180°', transform=ax.transAxes, ha='right', va='bottom', fontsize=8)
                    missing = 'Y181C' if 'Tyr181' in metric else 'Y188L'
                    ax.text(labels.index(missing), 0, 'N/A',ha='center',fontsize=7,color='.4')
                else:
                    ax.set_ylim(bottom=0)
                    if 'fraction' in metric:ax.set_ylim(0,1)
                ax.grid(axis='y',alpha=.2,lw=.5);ax.set_axisbelow(True)
                ax.spines[['top','right']].set_visible(False)
                if index == page*2:ax.legend(frameon=False,loc='upper left',bbox_to_anchor=(0,1.24),ncol=2,fontsize=9)
            fig.subplots_adjust(left=.07,right=.99,top=.91,bottom=.14,hspace=.85)
            fig.savefig(OUT/f'supplementary_figure_6_all_rt_metrics_{page+1}.png',dpi=300)
            pdf.savefig(fig)
            plt.close(fig)
    pd.DataFrame(exports).to_csv(OUT/'supplementary_figure_6_all_rt_metrics_values.csv',index=False)
    print(f'10 metrics; 19 genotypes; {len(exports)} values; 5 pages')


if __name__ == '__main__':
    main()
