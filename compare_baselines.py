"""
Cross-version comparison for hpv_faster_kenya scenario outputs.

Compares cumulative cancers (2025-2100) and cohort_cancers totals across baseline
folders (e.g. v2.2.6_baseline vs a future v2.3.0_baseline).

Usage:
  python compare_baselines.py --baselines v2.2.6_baseline v2.3.0_baseline
  python compare_baselines.py --baselines v2.2.6_baseline v2.3.0_baseline \
                              --location kenya --coverage 70 \
                              --scenarios "Baseline" "Catch-up 10-30: TTV"
"""
import argparse
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sciris as sc

import utils as ut


DEFAULT_SCENARIOS = [
    'Baseline',
    'Catch-up 10-20: V',
    'Catch-up 10-30: V',
    'Catch-up 10-40: V',
    'Catch-up 10-60: V',
    'Catch-up 10-60: TTV',
]


def compare(baselines, scenarios, location='kenya', coverage=70,
            resroot='results', outpath='figures/compare_baselines.png'):
    rows = []
    for base in baselines:
        resfolder = f'{resroot}/{base}'
        ts, cum, coh = ut.load_scens(resfolder, f'scens_{location}')
        for scen in scenarios:
            if cum is not None:
                v_canc = ut.get_cum(cum, scen, 'cancers', coverage=coverage, location=location)
                v_vax = ut.get_cum(cum, scen, 'vaccinations', coverage=coverage, location=location)
            else:
                v_canc = v_vax = (float('nan'), float('nan'), float('nan'))
            cohort_total = float('nan')
            if coh is not None:
                mask = (coh.scenario == scen)
                if 'coverage' in coh.columns: mask &= (coh.coverage == coverage)
                if 'location' in coh.columns: mask &= (coh.location == location)
                cohort_total = float(coh.loc[mask, 'value'].sum())
            rows.append({
                'baseline': base, 'scenario': scen,
                'cancers_2025_2100': v_canc[0],
                'vaccinations': v_vax[0],
                'cohort_cancers_total': cohort_total,
            })

    df = pd.DataFrame(rows)
    pivot_canc = df.pivot(index='scenario', columns='baseline', values='cancers_2025_2100')
    pivot_cohort = df.pivot(index='scenario', columns='baseline', values='cohort_cancers_total')

    print(f'\n=== Cumulative cancers 2025-2100 ({location}, coverage={coverage}%) ===')
    print(pivot_canc.to_string(float_format='%.0f'))
    print(f'\n=== Cohort-cancers 10-60 total ({location}, coverage={coverage}%) ===')
    print(pivot_cohort.to_string(float_format='%.0f'))

    # Bar plot: cohort-cancers per scenario per baseline
    base_colors = sc.gridcolors(len(baselines))
    fig, ax = plt.subplots(figsize=(max(8, 1.3 * len(scenarios)), 6))
    n = len(scenarios)
    x = np.arange(n)
    bar_w = 0.8 / max(1, len(baselines))
    for bi, base in enumerate(baselines):
        vals = pivot_cohort[base].reindex(scenarios).values
        ax.bar(x + bi * bar_w, vals, width=bar_w, color=base_colors[bi], label=base)
    ax.set_xticks(x + 0.5 * bar_w * (len(baselines) - 1))
    ax.set_xticklabels(scenarios, rotation=20, ha='right')
    ax.set_ylabel('Cohort cancers 10-60')
    ax.set_title(f'Cohort cancers by baseline — {location}, coverage {coverage}%')
    ax.legend(frameon=False)
    fig.tight_layout()
    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    fig.savefig(outpath, dpi=100)
    print(f'\nSaved figure: {outpath}')

    return df


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--baselines', nargs='+', required=True,
                        help='Baseline folder names under results/ (e.g. v2.2.6_baseline)')
    parser.add_argument('--scenarios', nargs='+', default=DEFAULT_SCENARIOS)
    parser.add_argument('--location', default='kenya')
    parser.add_argument('--coverage', type=int, default=70)
    parser.add_argument('--resroot', default='results')
    parser.add_argument('--outpath', default='figures/compare_baselines.png')
    args = parser.parse_args()

    compare(args.baselines, args.scenarios,
            location=args.location, coverage=args.coverage,
            resroot=args.resroot, outpath=args.outpath)
    print('\nDone.')
