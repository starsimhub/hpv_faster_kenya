'''
Run scenarios to evaluate how catch-up campaign effectiveness depends on
when routine vaccination started.

This script modifies the routine vaccination start year and evaluates
the marginal benefit of 2026 catch-up campaigns under different scenarios.
'''

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


# %% Functions with modified routine vaccination start year

def make_routine_vx_flexible(product='nonavalent', start_year=2019):
    """
    Modified version of make_routine_vx that accepts variable start year

    Args:
        product: vaccine product (default 'nonavalent')
        start_year: year to start routine vaccination (e.g., 2005, 2010, 2015, 2019)

    Returns:
        list of routine vaccination interventions
    """

    routine_age = (9, 10)
    eligibility = lambda sim: (sim.people.doses == 0)

    # Vaccination years from start_year through 2100
    vx_years = np.arange(start_year, 2100 + 1)

    # Scale-up pattern: 25% in year 1, reaching 90% by year 6
    scaleup_pattern = [0.25, 0.33, 0.29, 0.25, 0.62, 0.79]

    # If starting before 2019, assume scale-up took 6 years
    # If starting in/after 2019, use the actual Kenya scale-up pattern
    if start_year < 2019:
        # Use generic scale-up pattern
        scaleup = scaleup_pattern
    else:
        # Use actual Kenya scale-up
        scaleup = scaleup_pattern

    # Maintain 90% after scale-up
    final_cov = 0.9
    vx_cov = np.concatenate([scaleup + [final_cov] * (len(vx_years) - len(scaleup))])

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


def make_st(screen_coverage=0.15, triage_coverage=0.9, treat_coverage=0.75, start_year=2020):
    """Make screening & treatment intervention (copied from run_scenarios.py)"""

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
        prob=treat_coverage/4,
        product=hpv.radiation(),
        eligibility=radiation_eligible,
        label='radiation'
    )

    st_intvs = [screening, assign_treatment, ablation, excision, radiation]

    return st_intvs


def make_catchup_vx(product='nonavalent', catchup_cov=0.9, lower_age=10, upper_age=30,
                   start_year=2026, add_vx=True):
    """
    Simplified catch-up vaccination (vaccination only, no test-and-treat)
    to isolate the effect of prior routine vaccination coverage
    """

    intvs = []
    age_range = [lower_age, upper_age]

    if add_vx:
        eligibility = lambda sim: (sim.people.doses == 0)

        catchup_vx = hpv.campaign_vx(
            prob=catchup_cov,
            years=[start_year],
            product=product,
            age_range=age_range,
            eligibility=eligibility,
            interpolate=False,
            annual_prob=False,
            label='Catchup vx'
        )

        intvs = [catchup_vx]

    return intvs


# %% Main analysis function

def run_routine_timing_analysis(
    location='kenya',
    routine_start_years=[2005, 2010, 2015, 2019],
    catchup_age_ranges=[(10, 20), (10, 30), (10, 40)],
    catchup_year=2026,
    catchup_cov=0.9,
    n_seeds=3,  # Small number for testing
    debug=True,
):
    """
    Run scenarios varying the routine vaccination start year to assess
    impact on catch-up campaign effectiveness.

    Args:
        location: 'kenya' or 'nigeria'
        routine_start_years: list of years when routine vaccination could have started
        catchup_age_ranges: list of (lower, upper) age tuples for catch-up campaigns
        catchup_year: year to implement catch-up campaign
        catchup_cov: catch-up campaign coverage
        n_seeds: number of random seeds
        debug: if True, use smaller population and faster settings
    """

    # Load calibrated parameters
    calib_pars = sc.loadobj(f'results/{location}_pars.obj')

    # Dictionary to store results
    all_results = sc.objdict()

    # Loop over different routine start years
    for routine_start in routine_start_years:

        print(f'\n{"="*60}')
        print(f'Running scenarios with routine vaccination starting in {routine_start}')
        print(f'{"="*60}\n')

        # Calculate max age that would be vaccinated by routine program in catchup_year
        max_routine_age = 10 + (catchup_year - routine_start)
        print(f'By {catchup_year}, routine program will have covered ages 10-{max_routine_age}')

        # Background interventions with flexible routine start
        background_intvs = make_st() + make_routine_vx_flexible(start_year=routine_start)

        # Scenarios for this routine start year
        scenarios = sc.objdict()
        scenarios['Baseline'] = background_intvs

        # Add catch-up campaign scenarios
        for lower_age, upper_age in catchup_age_ranges:
            scen_name = f'Catchup {lower_age}-{upper_age}'

            # Calculate overlap with routine vaccination
            overlap_lower = max(lower_age, 10)
            overlap_upper = min(upper_age, max_routine_age)

            if overlap_upper >= overlap_lower:
                n_overlap = overlap_upper - overlap_lower + 1
                n_total = upper_age - lower_age + 1
                pct_overlap = 100 * n_overlap / n_total
                print(f'  {scen_name}: {pct_overlap:.0f}% overlap with routine-vaccinated ages')
            else:
                print(f'  {scen_name}: No overlap with routine-vaccinated ages')

            catchup_intvs = make_catchup_vx(
                catchup_cov=catchup_cov,
                lower_age=lower_age,
                upper_age=upper_age,
                start_year=catchup_year,
                add_vx=True,
            )
            scenarios[scen_name] = background_intvs + catchup_intvs

        # Run scenarios
        print(f'\nRunning {len(scenarios)} scenarios with {n_seeds} seeds each...')

        sims = []
        for scen_name, interventions in scenarios.items():
            for seed in range(n_seeds):

                # Create analyzers for age cohorts
                analyzers = []
                for cohort_age in [[10, 15], [15, 20], [20, 25], [25, 30], [30, 35],
                                  [35, 40], [40, 45], [45, 50]]:
                    analyzers.append(cohort_cancers(cohort_age=cohort_age, start=catchup_year))

                # Make sim
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
                sim.label = f'{scen_name}-{seed}'
                sims.append(sim)

        # Run all sims for this routine start year
        msim = hpv.MultiSim(sims)
        msim.run(verbose=0.2)

        # Process and store results
        all_results[f'routine_{routine_start}'] = process_results(
            msim, scenarios, n_seeds, catchup_year
        )

    # Save results
    sc.saveobj(f'raw_results/routine_timing_analysis_{location}.obj', all_results)
    print(f'\nSaved results to raw_results/routine_timing_analysis_{location}.obj')

    import utils as ut
    ut.routine_timing_to_csv(all_results,
                             f'results/routine_timing_{location}.csv',
                             catchup_year=catchup_year)
    print(f'Plot-ready CSV: results/routine_timing_{location}.csv')

    return all_results


def process_results(msim, scenarios, n_seeds, catchup_year):
    """
    Process simulation results to extract key metrics
    """

    # Split msim by scenario
    scen_labels = list(scenarios.keys())
    mlist = msim.split(chunks=len(scen_labels))

    results = sc.objdict()

    for si, scen_label in enumerate(scen_labels):

        # Extract results for this scenario
        reduced_sim = mlist[si].reduce(output=True)

        # Store key metrics
        results[scen_label] = sc.objdict()
        results[scen_label]['cancers'] = reduced_sim.results['cancers']
        results[scen_label]['year'] = reduced_sim.results['year']

        # Extract cohort-specific cancer outcomes
        analyzer_labels = [a.label for a in mlist[si].sims[0].analyzers]
        for alabel in analyzer_labels:
            base_analyzer = mlist[si].sims[0].get_analyzer(alabel)
            alist = [sim.get_analyzer(alabel) for sim in mlist[si].sims]
            reduced_analyzer = base_analyzer.reduce(alist)
            results[scen_label][alabel] = reduced_analyzer.cum_cancers_best

    return results


# %% Run as script

if __name__ == '__main__':

    T = sc.timer()

    # Quick test with limited scenarios
    results = run_routine_timing_analysis(
        location='kenya',
        routine_start_years=[2015, 2019, 2022],  # Just 3 scenarios to start
        catchup_age_ranges=[(10, 30)],  # Single catch-up range for simplicity
        catchup_year=2026,
        catchup_cov=0.9,
        n_seeds=2,  # Minimal seeds for quick test
        debug=True,  # Use debug mode for speed
    )

    T.toc('Done')

    print('\n' + '='*60)
    print('Analysis complete!')
    print('='*60)
    print('\nNext steps:')
    print('1. Review results in raw_results/routine_timing_analysis_kenya.obj')
    print('2. Create plotting script to visualize marginal benefit by routine start year')
    print('3. Expand to full parameter sweep (more start years, more catch-up ranges)')
