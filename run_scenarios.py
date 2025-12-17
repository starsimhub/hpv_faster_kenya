'''
Run HPVsim scenarios
Note: requires an HPC to run with debug=False; with debug=True, should take 5-15 min
to run.
'''


# %% General settings

import os

os.environ.update(
    OMP_NUM_THREADS='1',
    OPENBLAS_NUM_THREADS='1',
    NUMEXPR_NUM_THREADS='1',
    MKL_NUM_THREADS='1',
)

# Standard imports
import numpy as np
import sciris as sc
import hpvsim as hpv
import pandas as pd

# Imports from this repository
import run_sims as rs
import run_sims_nigeria as rsn
from analyzers import cohort_cancers


# What to run
debug = 0
n_seeds = [10, 1][debug]  # How many seeds to run per cluster


# %% Functions
def make_st(screen_coverage=0.15, triage_coverage=0.9, treat_coverage=0.75, start_year=2020):
    """ Make screening & treatment intervention """

    age_range = [30, 50]
    len_age_range = (age_range[1]-age_range[0])/2
    model_annual_screen_prob = 1 - (1 - screen_coverage)**(1/len_age_range)

    screen_eligible = lambda sim: np.isnan(sim.people.date_screened) | \
                                  (sim.t > (sim.people.date_screened + 5 / sim['dt']))
    screening = hpv.routine_screening(
        prob=model_annual_screen_prob,
        eligibility=screen_eligible,
        start_year=start_year,
        product='hpv',
        age_range=age_range,
        label='screening'
    )

    # Assign treatment
    screen_positive = lambda sim: sim.get_intervention('screening').outcomes['positive']
    assign_treatment = hpv.routine_triage(
        start_year=start_year,
        prob=triage_coverage,
        annual_prob=False,
        product='tx_assigner',
        eligibility=screen_positive,
        label='tx assigner'
    )

    ablation_eligible = lambda sim: sim.get_intervention('tx assigner').outcomes['ablation']
    ablation = hpv.treat_num(
        prob=treat_coverage,
        product='ablation',
        eligibility=ablation_eligible,
        label='ablation'
    )

    excision_eligible = lambda sim: list(set(sim.get_intervention('tx assigner').outcomes['excision'].tolist() +
                                             sim.get_intervention('ablation').outcomes['unsuccessful'].tolist()))
    excision = hpv.treat_num(
        prob=treat_coverage,
        product='excision',
        eligibility=excision_eligible,
        label='excision'
    )

    radiation_eligible = lambda sim: sim.get_intervention('tx assigner').outcomes['radiation']
    radiation = hpv.treat_num(
        prob=treat_coverage/4,  # assume an additional dropoff in CaTx coverage
        product=hpv.radiation(),
        eligibility=radiation_eligible,
        label='radiation'
    )

    st_intvs = [screening, assign_treatment, ablation, excision, radiation]

    return st_intvs


def make_routine_vx(product='nonavalent', start_year=2019):

    routine_age = (9, 10)
    eligibility = lambda sim: (sim.people.doses == 0)

    # Baseline vaccination scenarios
    vx_years = np.arange(start_year, 2100 + 1)
    scaleup = [0.25, 0.33, 0.29, 0.25, 0.62, 0.79]

    # Maintain 90%
    final_cov = 0.9
    vx_cov = np.concatenate([scaleup+[final_cov]*(len(vx_years)-len(scaleup))])

    routine_vx = hpv.campaign_vx(
        prob=vx_cov,
        years=vx_years,
        product=product,
        age_range=routine_age,
        eligibility=eligibility,
        interpolate=False,
        annual_prob=False,
        label='Routine vx'
    )

    return [routine_vx]


def make_hpv_test():
    hpvdna = hpv.default_dx('hpv')
    hpvdna.df = pd.read_csv('hpvdna.csv')
    return hpvdna


def make_catchup_vx(product='nonavalent', catchup_cov=0.9, lower_age=10, upper_age=15, start_year=2026, add_tt=False, add_vx=True):

    intvs = []
    age_range = [lower_age, upper_age]

    if add_vx:
        eligibility = lambda sim: (sim.people.doses == 0)

        routine_vx = hpv.campaign_vx(
            prob=catchup_cov,
            years=start_year,
            product=product,
            age_range=age_range,
            eligibility=eligibility,
            interpolate=False,
            annual_prob=False,
            label='Catchup vx'
        )

        intvs = [routine_vx]

    if add_tt:
        if add_vx:
            def is_el(sim):
                return ((sim.people.doses == 1) & (sim.people.date_vaccinated == sim.t) &
                        (sim.people.age >= age_range[0]) & (sim.people.age <= age_range[1]))
        else:
            def is_el(sim):
                return (sim.people.age >= age_range[0]) & (sim.people.age <= age_range[1])

        # Add test & treat
        primary = make_hpv_test()
        screening = hpv.campaign_screening(
            prob=catchup_cov,
            interpolate=False,
            annual_prob=False,
            eligibility=is_el,
            years=[start_year],
            product=primary,
            age_range=age_range,
            label='screening_older'
        )

        # Assign treatment
        tx_assigner = hpv.default_dx('tx_assigner')
        tx_assigner.df = pd.read_csv('tx_assigner_faster.csv')
        screen_positive = lambda sim: sim.get_intervention('screening_older').outcomes['positive']
        assign_treatment = hpv.campaign_triage(
            years=start_year,
            prob=1,  # ASSUME NO LTFU
            annual_prob=False,
            interpolate=False,
            product=tx_assigner,
            eligibility=screen_positive,
            label='tx_assigner_faster'
        )

        # Ablation treatment
        ablation_eligible = lambda sim: sim.get_intervention('tx_assigner_faster').outcomes['ablation']
        ablation = hpv.treat_num(
            prob=0.9,
            product='ablation',
            eligibility=ablation_eligible,
            label='ablation_faster'
        )
        # Excision treatment
        excision_eligible = lambda sim: list(set(sim.get_intervention('tx_assigner_faster').outcomes['excision'].tolist()
                                                + sim.get_intervention('ablation_faster').outcomes['unsuccessful'].tolist()))
        excision = hpv.treat_num(
            prob=0.9,
            product='excision',
            eligibility=excision_eligible,
            label='excision_faster'
        )
        # Radiation treatment
        radiation_eligible = lambda sim: sim.get_intervention('tx_assigner_faster').outcomes['radiation']
        radiation = hpv.treat_num(
            prob=0.9,  # assume an additional dropoff in CaTx coverage
            product=hpv.radiation(),
            eligibility=radiation_eligible,
            label='radiation_faster'
        )

        intvs += [
            screening, assign_treatment, ablation, excision, radiation,
        ]

    return intvs


def make_sims(location='kenya', calib_pars=None, scenarios=None, end=2100):
    """ Set up scenarios """

    all_msims = sc.autolist()
    for name, interventions in scenarios.items():
        sims = sc.autolist()
        for seed in range(n_seeds):
            print(f'Creating scenario "{name}"')

            # Create analyzers for each aga cohort
            analyzers = []
            for cohort_age in [[10, 15], [15, 20], [20, 25], [25, 30], [30, 35], [35, 40], [40, 45], [45, 50], [50, 55], [55, 60]]:
                analyzers.append(cohort_cancers(cohort_age=cohort_age, start=2026))

            if location == 'nigeria':
                sim = rsn.make_sim(location=location, calib_pars=calib_pars, debug=debug, interventions=interventions, analyzers=analyzers, end=end, seed=seed, verbose=-1)
            else:
                sim = rs.make_sim(location=location, calib_pars=calib_pars, debug=debug, interventions=interventions, analyzers=analyzers, end=end, seed=seed, verbose=-1)

            sim.label = name+f'-{seed}'
            sims += sim
        all_msims += hpv.MultiSim(sims)

    msim = hpv.MultiSim.merge(all_msims, base=False)

    return msim


def run_sims(location='kenya', calib_pars=None, scenarios=None, verbose=0.2):
    """ Run the simulations """
    msim = make_sims(location=location, calib_pars=calib_pars, scenarios=scenarios)
    msim.run(verbose=verbose)
    return msim


# %% Run as a script
if __name__ == '__main__':

    T = sc.timer()
    do_run = True
    do_save = False
    do_process = True
    location = 'kenya'  #'kenya'

    scenarios = dict()
    background_intvs = make_st() + make_routine_vx()
    scenarios['Baseline'] = background_intvs

    # Add scenarios
    lower_ages = [10, 15, 20, 25, 30, 35, 40, 45, 50, 55]
    age_bands = [5, 10, 15, 20, 25, 30, 40, 45, 50]
    for lower_age in lower_ages:
        for age_band in age_bands:
            if lower_age + age_band > 60:
                continue
            upper_age = lower_age + age_band
            for add_vx in [True, False]:
                for add_tt in [True, False]:
                    if not add_vx and not add_tt:
                        continue
                    scen_name = f'Catch-up {lower_age}-{upper_age}: '
                    if add_tt:
                        scen_name += 'TT'
                    if add_vx:
                        scen_name += 'V'
                    catchup_intvs = make_catchup_vx(lower_age=lower_age, upper_age=upper_age, add_tt=add_tt, add_vx=add_vx)
                    scenarios[scen_name] = background_intvs + catchup_intvs
    print(f'Defined {len(scenarios)} scenarios:' + ', '.join(scenarios.keys()))

    # Run scenarios (usually on VMs, runs n_seeds in parallel over M scenarios)
    if do_run:
        print(f'Running scenarios for location: {location}')

        calib_pars = sc.loadobj(f'results/{location}_pars.obj')
        msim = run_sims(location=location, calib_pars=calib_pars, scenarios=scenarios, verbose=-1)

        if do_save: sc.saveobj(f'results/scens_{location}.msim', msim)

    if do_process:
        print('Post-processing results...')
        if not do_run:
            msim = sc.loadobj(f'results/scens_{location}.msim')

        metrics = ['year', 'asr_cancer_incidence', 'cancers', 'cancer_deaths']

        # Process results
        scen_labels = list(scenarios.keys())
        print(f'Processing {len(scen_labels)} scenarios: ' + ', '.join(scen_labels))
        mlist = msim.split(chunks=len(scen_labels))

        msim_dict = sc.objdict()
        for si, scen_label in enumerate(scen_labels):

            mres = sc.objdict()

            # Deal with analyzers
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

            # Pull out characteristics of sim to decide what resources we need
            programs = {
                "Catchup vx": "vaccinations",
            }
            for intv_name in set(programs.values()): mres[intv_name] = np.zeros_like(mres.year)
            for intv_name, df_key in programs.items():
                if reduced_sim.get_intervention(intv_name, die=False) is not None:
                    mres[df_key] += reduced_sim.get_intervention(intv_name).n_products_used.values

            msim_dict[scen_label] = mres

        sc.saveobj(f'results/scens_{location}.obj', msim_dict)

    print('Done.')
