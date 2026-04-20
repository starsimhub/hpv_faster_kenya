"""
Run the 6 scenarios for the IARC comparison plot, 6 seeds each.

Scenarios:
    1. No interventions
    2. Routine vax only (baseline: routine screening + routine vax)
    3. Catch-up vax to 30 (baseline + CU vax 10-30)
    4. Catch-up vax to 40 (baseline + CU vax 10-40)
    5. Catch-up vax to 40 + one-time S&T to 65 (baseline + CU vax 10-40 + campaign S&T 40-65)
    6. One-time S&T 10-65 only (no vaccination, just screen & treat)
"""

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

import run_sims as rs
from analyzers import cohort_cancers
from run_scenarios import make_routine_vx, make_catchup_vx, make_st, make_hpv_test


# Settings
debug = 0
n_seeds = [6, 2][debug]
location = 'kenya'
coverage = 70
catchup_cov = coverage / 100


def make_campaign_st(catchup_cov=0.7, lower_age=40, upper_age=65, start_year=2026):
    """Make a one-time campaign screen & treat intervention."""

    primary = make_hpv_test()
    screening = hpv.campaign_screening(
        prob=catchup_cov,
        interpolate=False,
        annual_prob=False,
        years=[start_year],
        product=primary,
        age_range=[lower_age, upper_age],
        label='screening_campaign'
    )

    tx_assigner = hpv.default_dx('tx_assigner')
    tx_assigner.df = pd.read_csv('tx_assigner_faster.csv')
    screen_positive = lambda sim: sim.get_intervention('screening_campaign').outcomes['positive']
    assign_treatment = hpv.campaign_triage(
        years=start_year,
        prob=1,
        annual_prob=False,
        interpolate=False,
        product=tx_assigner,
        eligibility=screen_positive,
        label='tx_assigner_campaign'
    )

    ablation_eligible = lambda sim: sim.get_intervention('tx_assigner_campaign').outcomes['ablation']
    ablation = hpv.treat_num(
        prob=0.9,
        product='ablation',
        eligibility=ablation_eligible,
        label='ablation_campaign'
    )

    excision_eligible = lambda sim: list(set(
        sim.get_intervention('tx_assigner_campaign').outcomes['excision'].tolist()
        + sim.get_intervention('ablation_campaign').outcomes['unsuccessful'].tolist()
    ))
    excision = hpv.treat_num(
        prob=0.9,
        product='excision',
        eligibility=excision_eligible,
        label='excision_campaign'
    )

    radiation_eligible = lambda sim: sim.get_intervention('tx_assigner_campaign').outcomes['radiation']
    radiation = hpv.treat_num(
        prob=0.9,
        product=hpv.radiation(),
        eligibility=radiation_eligible,
        label='radiation_campaign'
    )

    return [screening, assign_treatment, ablation, excision, radiation]


def define_scenarios():
    """Define the 6 scenarios."""
    background_intvs = make_st() + make_routine_vx()

    scenarios = sc.objdict()
    scenarios['No interventions'] = []
    scenarios['Routine vax only'] = background_intvs
    scenarios['Catch-up to 30'] = background_intvs + make_catchup_vx(
        catchup_cov=catchup_cov, lower_age=10, upper_age=30, add_tt=False, add_vx=True
    )
    scenarios['Catch-up to 40'] = background_intvs + make_catchup_vx(
        catchup_cov=catchup_cov, lower_age=10, upper_age=40, add_tt=False, add_vx=True
    )
    scenarios['Catch-up to 40 + S&T to 65'] = background_intvs + make_catchup_vx(
        catchup_cov=catchup_cov, lower_age=10, upper_age=40, add_tt=False, add_vx=True
    ) + make_campaign_st(catchup_cov=catchup_cov, lower_age=10, upper_age=65)
    scenarios['S&T 10-65 only'] = background_intvs + make_campaign_st(
        catchup_cov=catchup_cov, lower_age=10, upper_age=65
    )

    return scenarios


def run_scenarios():
    """Run all scenarios with n_seeds each."""
    scenarios = define_scenarios()
    calib_pars = sc.loadobj(f'results/{location}_pars.obj')

    all_msims = sc.autolist()
    for name, interventions in scenarios.items():
        sims = sc.autolist()
        for seed in range(n_seeds):
            print(f'Creating scenario "{name}", seed {seed}')

            analyzers = []
            for cohort_age in [[10,15],[15,20],[20,25],[25,30],[30,35],[35,40],[40,45],[45,50],[50,55],[55,60]]:
                analyzers.append(cohort_cancers(cohort_age=cohort_age, start=2026))

            sim = rs.make_sim(
                location=location, calib_pars=calib_pars, debug=debug,
                interventions=interventions, analyzers=analyzers,
                end=2100, seed=seed, verbose=-1
            )
            sim.label = f'{name}-{seed}'
            sims += sim
        all_msims += hpv.MultiSim(sims)

    msim = hpv.MultiSim.merge(all_msims, base=False)
    print(f'Running {len(msim.sims)} simulations...')
    msim.run(verbose=-1)

    return msim, scenarios


def process_results(msim, scenarios):
    """Process results into a dict with medians and quantiles."""
    metrics = ['year', 'asr_cancer_incidence', 'cancers', 'cancer_deaths']
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
            reduced_analyzer = base_analyzer.reduce(alist)
            mres[alabel] = reduced_analyzer.cum_cancers_best
            mres[f'{alabel}_low'] = reduced_analyzer.cum_cancers_low
            mres[f'{alabel}_high'] = reduced_analyzer.cum_cancers_high
            mres[f'raw_{alabel}'] = reduced_analyzer.raw

        reduced_sim = mlist[si].reduce(output=True)
        for metric in metrics:
            mres[metric] = reduced_sim.results[metric]

        msim_dict[scen_label] = mres

    return msim_dict


if __name__ == '__main__':
    T = sc.timer()

    msim, scenarios = run_scenarios()
    sc.saveobj(f'raw_results/iarc_comparison_{location}_{coverage}.msim', msim)

    msim_dict = process_results(msim, scenarios)
    sc.saveobj(f'raw_results/iarc_comparison_{location}_{coverage}.obj', msim_dict)

    # Extract plot-ready CSVs for cross-version comparison
    import utils as ut
    ut.msim_dict_to_csvs(msim_dict, f'iarc_{location}_{coverage}',
                         resfolder='results', location=location, coverage=coverage,
                         year_start=2020, ts_metrics=ut.TS_METRICS_IARC)
    print(f'Saved plot-ready CSVs: results/iarc_{location}_{coverage}_*.csv')

    T.toc()
    print('Done!')
