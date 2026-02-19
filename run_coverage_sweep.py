"""
Sweep catch-up vaccination coverage for varying upper age targets and
plot impact (% cancers averted in the 10-60 cohort) vs coverage level.

Usage:
    python run_coverage_sweep.py

Set `do_run = True` to run new simulations (requires calibrated parameters
and significant compute time — intended for HPC). Set `do_run = False` to
reuse previously saved results and only regenerate the figure.
"""

# %% General settings
import os

os.environ.update(
    OMP_NUM_THREADS='1',
    OPENBLAS_NUM_THREADS='1',
    NUMEXPR_NUM_THREADS='1',
    MKL_NUM_THREADS='1',
)

import numpy as np
import sciris as sc
import hpvsim as hpv
import pandas as pd
import pylab as pl

import run_sims as rs
from run_scenarios import make_st, make_routine_vx, make_catchup_vx
from analyzers import cohort_cancers

# %% Configuration
debug = 0
n_seeds = [2, 1][debug]

# Coverage levels to sweep (%)
coverage_pcts = [10, 20, 30, 40, 50, 60, 70, 80, 90]

# Lower age is always 10; vary the upper age of the catch-up campaign
lower_age = 10
upper_ages = [30, 35, 40, 45]

# All cohorts tracked by analyzers (10-60) — used for impact denominator
all_cohort_ages = [[10, 15], [15, 20], [20, 25], [25, 30], [30, 35],
                   [35, 40], [40, 45], [45, 50], [50, 55], [55, 60]]
all_cohort_keys = ['10_15', '15_20', '20_25', '25_30', '30_35',
                   '35_40', '40_45', '45_50', '50_55', '55_60']

location = 'kenya'


# %% Helper functions
def make_sims(calib_pars, scenarios):
    """Create MultiSim for a set of scenarios."""
    all_msims = sc.autolist()
    for name, interventions in scenarios.items():
        sims = sc.autolist()
        for seed in range(n_seeds):
            analyzers = [cohort_cancers(cohort_age=ca, start=2026) for ca in all_cohort_ages]
            sim = rs.make_sim(
                location=location,
                calib_pars=calib_pars,
                debug=debug,
                interventions=interventions,
                analyzers=analyzers,
                end=2100,
                seed=seed,
                verbose=-1,
            )
            sim.label = f'{name}-{seed}'
            sims += sim
        all_msims += hpv.MultiSim(sims)
    return hpv.MultiSim.merge(all_msims, base=False)


def process_msim(msim, scenarios):
    """Extract cohort cancer counts from a completed MultiSim."""
    scen_labels = list(scenarios.keys())
    mlist = msim.split(chunks=len(scen_labels))

    msim_dict = sc.objdict()
    for si, scen_label in enumerate(scen_labels):
        mres = sc.objdict()

        # Process analyzers
        analyzer_labels = [a.label for a in mlist[si].sims[0].analyzers]
        for alabel in analyzer_labels:
            base_analyzer = mlist[si].sims[0].get_analyzer(alabel)
            alist = [sim.get_analyzer(alabel) for sim in mlist[si].sims]
            reduced = base_analyzer.reduce(alist)
            mres[alabel] = reduced.cum_cancers_best
            mres[f'{alabel}_low'] = reduced.cum_cancers_low
            mres[f'{alabel}_high'] = reduced.cum_cancers_high

        msim_dict[scen_label] = mres

    return msim_dict


def sum_cohort_cancers(mres, cohort_keys, suffix=''):
    """Sum cumulative cancers across the specified cohorts."""
    total = 0
    for c in cohort_keys:
        key = f'cohort_cancers_{c}{suffix}'
        if key in mres:
            total += mres[key]
    return total


# %% Main script
if __name__ == '__main__':

    T = sc.timer()

    do_run = True       # Set False to skip simulations and load saved results
    do_save = True
    do_plot = False

    results_file = f'raw_results/coverage_sweep_{location}_ttv.obj'

    # ------------------------------------------------------------------
    # 1. Run simulations (or load from disk)
    # ------------------------------------------------------------------
    if do_run:
        print(f'Running coverage sweep for {location}')
        calib_pars = sc.loadobj(f'results/{location}_pars.obj')

        # Background interventions (screening + routine vaccination)
        background_intvs = make_st() + make_routine_vx()

        # Build scenario dict: Baseline + one scenario per (upper_age, coverage) combo
        scenarios = sc.objdict()
        scenarios['Baseline'] = background_intvs

        for ua in upper_ages:
            for cov_pct in coverage_pcts:
                cov = cov_pct / 100
                catchup_intvs = make_catchup_vx(
                    catchup_cov=cov,
                    lower_age=lower_age,
                    upper_age=ua,
                    add_tt=True,
                    add_vx=True,
                )
                scen_name = f'Catch-up {lower_age}-{ua} {cov_pct}%'
                scenarios[scen_name] = background_intvs + catchup_intvs

        print(f'Defined {len(scenarios)} scenarios')

        msim = make_sims(calib_pars, scenarios)
        msim.run(verbose=-1)

        msim_dict = process_msim(msim, scenarios)

        if do_save:
            sc.saveobj(results_file, msim_dict)
            print(f'Saved results to {results_file}')

            # Also save the raw msim
            sc.saveobj(results_file.replace('.obj', '.msim'), msim)
    else:
        print(f'Loading previously saved results from {results_file}')
        msim_dict = sc.loadobj(results_file)

    # ------------------------------------------------------------------
    # 2. Compute impact (% cancers averted in 10-60 cohort) for each
    #    (upper_age, coverage) combination
    # ------------------------------------------------------------------
    baseline_cancers      = sum_cohort_cancers(msim_dict['Baseline'], all_cohort_keys)
    baseline_cancers_low  = sum_cohort_cancers(msim_dict['Baseline'], all_cohort_keys, suffix='_low')
    baseline_cancers_high = sum_cohort_cancers(msim_dict['Baseline'], all_cohort_keys, suffix='_high')

    print(f'\nBaseline cancers in 10-60 cohort: '
          f'{baseline_cancers:,.0f} ({baseline_cancers_low:,.0f} - {baseline_cancers_high:,.0f})')

    rows = []
    for ua in upper_ages:
        for cov_pct in coverage_pcts:
            scen_key = f'Catch-up {lower_age}-{ua} {cov_pct}%'
            mres = msim_dict[scen_key]

            scen_cancers      = sum_cohort_cancers(mres, all_cohort_keys)
            scen_cancers_low  = sum_cohort_cancers(mres, all_cohort_keys, suffix='_low')
            scen_cancers_high = sum_cohort_cancers(mres, all_cohort_keys, suffix='_high')

            averted     = baseline_cancers - scen_cancers
            pct_averted = 100 * averted / baseline_cancers if baseline_cancers > 0 else 0

            # Uncertainty bounds (low cancers -> high averted)
            averted_low  = baseline_cancers_low  - scen_cancers_high
            averted_high = baseline_cancers_high - scen_cancers_low
            pct_averted_low  = 100 * averted_low  / baseline_cancers_high if baseline_cancers_high > 0 else 0
            pct_averted_high = 100 * averted_high / baseline_cancers_low  if baseline_cancers_low  > 0 else 0

            rows.append(dict(
                upper_age=ua,
                coverage_pct=cov_pct,
                cancers_baseline=baseline_cancers,
                cancers_scenario=scen_cancers,
                cancers_averted=averted,
                pct_averted=pct_averted,
                pct_averted_low=pct_averted_low,
                pct_averted_high=pct_averted_high,
            ))

            print(f'  {lower_age}-{ua} @ {cov_pct:3d}%: {scen_cancers:,.0f} cancers, '
                  f'{averted:,.0f} averted ({pct_averted:.1f}%)')

    df = pd.DataFrame(rows)
    csv_path = f'results/coverage_sweep_{location}.csv'
    df.to_csv(csv_path, index=False)
    print(f'\nSummary saved to {csv_path}')

    # ------------------------------------------------------------------
    # 3. Plot: coverage vs impact, one line per upper age
    # ------------------------------------------------------------------
    if do_plot:
        try:
            import utils as ut
            ut.set_font(20)
        except Exception:
            pass

        fig, ax = pl.subplots(figsize=(9, 7))

        # Colors for each upper-age series
        colors = {30: '#1f77b4', 35: '#ff7f0e', 40: '#2ca02c', 45: '#d62728'}

        for ua in upper_ages:
            sub = df[df['upper_age'] == ua].sort_values('coverage_pct')
            cov = sub['coverage_pct'].values
            imp = sub['pct_averted'].values
            imp_lo = sub['pct_averted_low'].values
            imp_hi = sub['pct_averted_high'].values

            color = colors[ua]
            ax.fill_between(cov, imp_lo, imp_hi, alpha=0.15, color=color)
            ax.plot(cov, imp, '-o', color=color, linewidth=2.5, markersize=7,
                    zorder=5, label=f'Catch-up {lower_age}-{ua}')

        # --- Callouts on the 10-30 line ---
        ref_ua = 30
        sub30 = df[df['upper_age'] == ref_ua].sort_values('coverage_pct')
        cov30 = sub30['coverage_pct'].values
        imp30 = sub30['pct_averted'].values

        # Interpolate impact at 70% and 35% coverage
        imp_at_70 = np.interp(70, cov30, imp30)
        imp_at_35 = np.interp(35, cov30, imp30)

        ref_color = colors[ref_ua]

        # Horizontal dashed lines from the y-axis to each point
        ax.plot([0, 70], [imp_at_70, imp_at_70], ls=':', color=ref_color, linewidth=1, zorder=2)
        ax.plot([0, 35], [imp_at_35, imp_at_35], ls=':', color=ref_color, linewidth=1, zorder=2)

        # Vertical dashed lines from the x-axis to each point
        ax.plot([70, 70], [0, imp_at_70], ls=':', color=ref_color, linewidth=1, zorder=2)
        ax.plot([35, 35], [0, imp_at_35], ls=':', color=ref_color, linewidth=1, zorder=2)

        # Marker dots
        ax.plot(70, imp_at_70, 'o', color=ref_color, markersize=10, zorder=6)
        ax.plot(35, imp_at_35, 'o', color=ref_color, markersize=10, zorder=6)

        # Shared style for annotation text boxes
        bbox_style = dict(boxstyle='round,pad=0.4', facecolor='white',
                          edgecolor='0.7', alpha=0.95)

        # Text callout at 70%
        ax.annotate(
            f'70% coverage\n{imp_at_70:.0f}% cancers averted',
            xy=(70, imp_at_70), xytext=(75, imp_at_70 * 1.45),
            fontsize=11, ha='left', zorder=10,
            bbox=bbox_style,
            arrowprops=dict(arrowstyle='->', color='0.3', lw=1.2),
        )

        # Text callout at 35%
        ax.annotate(
            f'Half the coverage (35%)\nstill averts {imp_at_35:.0f}% of cancers',
            xy=(35, imp_at_35), xytext=(40, imp_at_35 * 0.35),
            fontsize=11, ha='left', zorder=10,
            bbox=bbox_style,
            arrowprops=dict(arrowstyle='->', color='0.3', lw=1.2),
        )

        ax.set_xlabel('Catch-up campaign coverage (%)')
        ax.set_ylabel('Cancers averted in 10-60 cohort (%)')
        ax.set_title(f'Impact of catch-up vaccination by coverage level\n({location.capitalize()})')
        ax.set_xlim(0, 100)
        ax.set_ylim(bottom=0)
        ax.legend(loc='upper left', frameon=False)
        ax.grid(alpha=0.3, linestyle='--')
        sc.SIticks(ax)

        fig.tight_layout()
        fig_path = f'figures/coverage_sweep_{location}.png'
        sc.savefig(fig_path, dpi=200)
        print(f'Figure saved to {fig_path}')

        # --------------------------------------------------------------
        # 3b. Relative-benefit plot: % of benefit at 70% coverage
        #     (10-30 line only, coverage up to 70%)
        # --------------------------------------------------------------
        fig2, ax2 = pl.subplots(figsize=(9, 7))

        ref_ua = 30
        sub30 = df[df['upper_age'] == ref_ua].sort_values('coverage_pct')
        sub30 = sub30[sub30['coverage_pct'] <= 70]
        cov30 = sub30['coverage_pct'].values
        imp30 = sub30['pct_averted'].values

        ref_imp30 = imp30[cov30 == 70][0]
        rel30 = 100 * imp30 / ref_imp30

        ref_color = colors[ref_ua]
        ax2.plot(cov30, rel30, '-o', color=ref_color, linewidth=2.5, markersize=7,
                 zorder=5, label=f'Catch-up {lower_age}-{ref_ua}')

        # Interpolate at 35%
        rel_at_70 = 100.0  # by definition
        rel_at_35 = np.interp(35, cov30, rel30)

        # Guide lines and dots
        ax2.plot([0, 70], [rel_at_70, rel_at_70], ls=':', color=ref_color, linewidth=1, zorder=2)
        ax2.plot([0, 35], [rel_at_35, rel_at_35], ls=':', color=ref_color, linewidth=1, zorder=2)
        ax2.plot([70, 70], [0, rel_at_70], ls=':', color=ref_color, linewidth=1, zorder=2)
        ax2.plot([35, 35], [0, rel_at_35], ls=':', color=ref_color, linewidth=1, zorder=2)
        ax2.plot(70, rel_at_70, 'o', color=ref_color, markersize=10, zorder=6)
        ax2.plot(35, rel_at_35, 'o', color=ref_color, markersize=10, zorder=6)

        bbox_style2 = dict(boxstyle='round,pad=0.4', facecolor='white',
                           edgecolor='0.7', alpha=0.95)

        ax2.annotate(
            f'70% coverage = baseline\n(100% of benefit)',
            xy=(70, rel_at_70), xytext=(50, rel_at_70 * 1.15),
            fontsize=11, ha='left', zorder=10,
            bbox=bbox_style2,
            arrowprops=dict(arrowstyle='->', color='0.3', lw=1.2),
        )

        ax2.annotate(
            f'Half the coverage (35%)\nstill captures {rel_at_35:.0f}% of benefit',
            xy=(35, rel_at_35), xytext=(40, rel_at_35 * 0.55),
            fontsize=11, ha='left', zorder=10,
            bbox=bbox_style2,
            arrowprops=dict(arrowstyle='->', color='0.3', lw=1.2),
        )

        ax2.set_xlabel('Catch-up campaign coverage (%)')
        ax2.set_ylabel('% of benefit achieved at 70% coverage')
        ax2.set_title(f'Relative impact of catch-up vaccination ({lower_age}-{ref_ua})\n({location.capitalize()})')
        ax2.set_xlim(0, 75)
        ax2.set_ylim(0, 115)
        ax2.grid(alpha=0.3, linestyle='--')
        sc.SIticks(ax2)

        fig2.tight_layout()
        fig2_path = f'figures/coverage_sweep_relative_{location}.png'
        sc.savefig(fig2_path, dpi=200)
        print(f'Figure saved to {fig2_path}')

    T.toc('Done')