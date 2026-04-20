"""
Plot age at causal infection distributions.

Reads plot-ready histograms from `<resfolder>/age_causal_<location>.csv`
(0.5-yr bins, produced by run_sims.py's extract step).
"""
import argparse
import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as pl
import seaborn as sns
import sciris as sc
from scipy.stats import gaussian_kde

import utils as ut


EVENT_ORDER = ['Causal\ninfection', 'HSIL', 'Cancer']


def _load_and_expand(location, resfolder):
    """Return {event: ages array} re-expanded from the histogram CSV."""
    path = f'{resfolder}/age_causal_{location}.csv'
    hist = ut.load_age_causal_hist(path)
    return {event: ut.expand_age_hist(df) for event, df in hist.items()}


def plot_age_causal(location='kenya', resfolder='results/v2.2.6_baseline',
                    outpath=None):
    ages_by_event = _load_and_expand(location, resfolder)
    rows = []
    for event in EVENT_ORDER:
        for a in ages_by_event.get(event, []):
            rows.append({'Health event': event, 'Age': float(a)})
    ac_df = pd.DataFrame(rows)

    ac_colors = sc.gridcolors(3)
    ut.set_font(20)
    fig, ax = pl.subplots(1, 1, figsize=(6, 5))
    sns.boxplot(x='Health event', y='Age', data=ac_df, ax=ax,
                showfliers=False, palette=ac_colors, hue='Health event',
                hue_order=EVENT_ORDER)
    ax.set_title('Age distribution\nof key health events')
    ax.set_xlabel('')
    ax.set_ylim([0, 100])
    fig.tight_layout()
    outpath = outpath or f'figures/age_causal_{location}.png'
    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    sc.savefig(outpath, dpi=100)
    print(f'Figure saved to: {outpath}')
    return fig


def plot_age_causal_bar(location='kenya', resfolder='results/v2.2.6_baseline',
                        outpath=None):
    ages_by_event = _load_and_expand(location, resfolder)
    ages = ages_by_event.get('Causal\ninfection')
    if ages is None or len(ages) == 0:
        print('No causal infection data'); return None

    bins = np.arange(10, 81)
    counts, _ = np.histogram(ages, bins=bins)
    pcts = counts / counts.sum() * 100
    bin_centers = bins[:-1]

    ut.set_font(14)
    fig, ax = pl.subplots(1, 1, figsize=(4.5, 3.2))
    bar_color = '#3498db'
    ax.bar(bin_centers, pcts, color=bar_color, edgecolor='none', width=1.0, alpha=0.5)

    kde = gaussian_kde(ages, bw_method=0.15)
    x_smooth = np.linspace(10, 80, 500)
    kde_vals = kde(x_smooth)
    kde_scaled = kde_vals / kde_vals.sum() * pcts.sum() * (500 / 70)
    ax.plot(x_smooth, kde_scaled, color='#1a5276', linewidth=2)

    q25, median, q75 = np.percentile(ages, [25, 50, 75])
    print(f'{location}: 25th={q25:.1f}, Median={median:.1f}, 75th={q75:.1f}')
    ymax = max(pcts)
    for val, label, yoff in [(q25, '25th', 1.18), (median, 'Median', 1.08), (q75, '75th', 1.18)]:
        color = '#e74c3c' if label == 'Median' else '#555555'
        lw = 2 if label == 'Median' else 1.5
        ax.axvline(val, color=color, linestyle='--', linewidth=lw, zorder=5)
        ax.text(val, ymax * yoff, label, ha='center', fontsize=9, color=color, fontweight='bold',
                bbox=dict(boxstyle='round,pad=0.2', facecolor='white', edgecolor='none', alpha=1.0), zorder=10)

    ax.set_xticks(np.arange(10, 81, 10))
    ax.set_xlim(9.5, 80.5)
    ax.set_xlabel('Age', fontsize=13)
    ax.set_ylabel('% of causal infections', fontsize=13)
    ax.set_ylim(0, ymax * 1.25)
    ax.yaxis.set_major_formatter(pl.FuncFormatter(lambda x, _: f'{x:.0f}%'))
    ax.tick_params(axis='y', labelsize=11)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    fig.tight_layout()
    outpath = outpath or f'figures/age_causal_bar_{location}.png'
    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    sc.savefig(outpath, dpi=300)
    print(f'Figure saved to: {outpath}')
    return fig


def plot_age_causal_violin(location='kenya', resfolder='results/v2.2.6_baseline',
                           outpath=None):
    ages_by_event = _load_and_expand(location, resfolder)
    ages = ages_by_event.get('Causal\ninfection')
    if ages is None or len(ages) == 0:
        print('No causal infection data'); return None

    vax_ages = [15, 20, 25, 30, 35, 40]
    percentages = [(ages <= age).sum() / len(ages) * 100 for age in vax_ages]
    for age, pct in zip(vax_ages, percentages):
        print(f'Age {age}: {pct:.1f}% of causal infections')

    ut.set_font(20)
    fig, ax = pl.subplots(1, 1, figsize=(6, 8))
    parts = ax.violinplot([ages], positions=[0], widths=0.6,
                          showmeans=False, showmedians=True, vert=True)
    for pc in parts['bodies']:
        pc.set_facecolor('#3498db'); pc.set_alpha(0.7)
        pc.set_edgecolor('black'); pc.set_linewidth(1.5)
    parts['cmedians'].set_edgecolor('red')
    parts['cmedians'].set_linewidth(2)

    colors = ['#27ae60', '#1abc9c', '#f39c12', '#e74c3c', '#9b59b6', '#34495e']
    for age, pct, color in zip(vax_ages, percentages, colors):
        ax.axhline(y=age, color=color, linestyle='--', linewidth=2, alpha=0.6, zorder=5)
        ax.text(0.45, age, f'{pct:.0f}%',
                fontsize=13, fontweight='bold', color=color,
                va='center', ha='left', zorder=15,
                bbox=dict(boxstyle='round,pad=0.4', facecolor='white',
                          edgecolor=color, linewidth=2.5, alpha=1.0))

    ax.set_ylabel('Age at causal infection (years)', fontsize=12, fontweight='bold')
    ax.set_title('Distribution of age at causal infection\nwith vaccination coverage thresholds',
                 fontsize=14, fontweight='bold')
    ax.set_xlim(-0.5, 0.6)
    ax.set_ylim([10, 55])
    ax.set_xticks([])
    ax.grid(alpha=0.3, linestyle='--', axis='y')
    textstr = 'Lines show % of causal\ninfections occurring by\neach age (vaccination target)'
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    ax.text(0.98, 0.02, textstr, transform=ax.transAxes, fontsize=10,
            verticalalignment='bottom', horizontalalignment='right', bbox=props)

    fig.tight_layout()
    outpath = outpath or f'figures/age_causal_violin_{location}.png'
    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    sc.savefig(outpath, dpi=300)
    print(f'Figure saved to: {outpath}')
    return fig


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--location', nargs='+', default=['kenya', 'nigeria'])
    parser.add_argument('--resfolder', default='results/v2.2.6_baseline')
    args = parser.parse_args()

    for loc in args.location:
        plot_age_causal(location=loc, resfolder=args.resfolder)
        plot_age_causal_violin(location=loc, resfolder=args.resfolder)
        plot_age_causal_bar(location=loc, resfolder=args.resfolder)
