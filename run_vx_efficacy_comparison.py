"""
Simple vaccine efficacy comparison across age cohorts.

Answers the question: in an unvaccinated population, what is the expected
lifetime cancer reduction from vaccinating 100% of women at a given age?

Methodology: for each vaccination age, run a single sim that vaccinates
100% of women in that single-year age cohort. A custom analyzer identifies
the cohort members at the time of vaccination, then tracks cancers among
those specific agents (by index) for the remainder of the sim. Cancer counts
are compared to those from the same agents in a baseline (no-vax) sim using
the same seed.

Scenarios:
  - Baseline: no vaccination, no screening (unvaccinated population)
  - Vax at age X: vaccinate 100% of X-year-olds in 2026

This is designed to be a simple, reproducible comparison that other models
(e.g., PRIME, Harvard, WHO) can replicate.
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
import hpvsim.utils as hpu
import run_sims as rs


# Settings
debug = 0
n_seeds = [5, 1][debug]
vx_ages = [9, 15, 20, 25, 30, 35, 40]
start_year = 2026
location = 'kenya'


class fixed_cohort_cancers(hpv.Analyzer):
    """Track lifetime cancers among agents identified by index at start_year.

    Uses RAW agent counts (not scale-weighted) to avoid bias from the multi-scale
    mechanism, which reduces Level0 agent scale factors when Level1 clones are
    created for cancer resolution.
    """

    def __init__(self, cohort_age=None, start=None, **kwargs):
        super().__init__(**kwargs)
        self.start = start or 2026
        self.cohort_age = cohort_age or [9, 10]
        self.label = f'fixed_cohort_{self.cohort_age[0]}_{self.cohort_age[1]}'
        self.cohort_inds = None
        self.cancers = 0
        self.n_cohort = 0
        self.identified = False

    def initialize(self, sim):
        super().initialize()

    def apply(self, sim):
        ppl = sim.people
        year = sim.yearvec[sim.t]

        # Identify the cohort once at the start year
        if not self.identified and year >= self.start:
            age_lo = self.cohort_age[0]
            age_hi = self.cohort_age[1]
            in_cohort = ((ppl.sex == 0) & (ppl.age >= age_lo) & (ppl.age < age_hi)
                         & ppl.alive & ppl.level0)
            self.cohort_inds = set(hpu.true(in_cohort).tolist())
            self.n_cohort = len(self.cohort_inds)  # Raw count
            self.identified = True

        if not self.identified:
            return

        # Track new cancers among cohort members (any genotype), raw count
        new_cancer = (ppl.date_cancerous == sim.t).any(axis=0)
        if new_cancer.any():
            cancer_inds = set(hpu.true(new_cancer).tolist())
            cohort_cancer = cancer_inds & self.cohort_inds
            if cohort_cancer:
                self.cancers += len(cohort_cancer)


def make_vx_at_age(age, product='nonavalent'):
    """Create a one-time 100% vaccination campaign for a single-year age cohort."""
    eligibility = lambda sim: (sim.people.doses == 0)
    vx = hpv.campaign_vx(
        prob=1.0,
        years=[start_year],
        product=product,
        age_range=(age, age + 1),
        eligibility=eligibility,
        interpolate=False,
        annual_prob=False,
        label=f'Vax age {age}'
    )
    return [vx]


def make_sim(name, interventions, calib_pars, seed=0):
    """Create a single sim for a given scenario."""
    analyzers = []
    for age in vx_ages:
        analyzers.append(fixed_cohort_cancers(cohort_age=[age, age + 1], start=start_year))

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
    sim.label = f'{name}_s{seed}'
    return sim


def run_one_sim(sim):
    """Run a single sim and extract results."""
    sim.run()
    res = dict()
    for age in vx_ages:
        alabel = f'fixed_cohort_{age}_{age + 1}'
        analyzer = sim.get_analyzer(alabel)
        res[f'cancers_{age}'] = analyzer.cancers
        res[f'n_cohort_{age}'] = analyzer.n_cohort
    return res


def run_sims():
    """Run all scenarios in parallel across multiple seeds and average."""
    calib_pars = sc.loadobj(f'results/{location}_pars.obj')

    # Build scenario definitions
    scenarios = dict()
    scenarios['Baseline'] = []
    for age in vx_ages:
        scenarios[f'Vax at {age}'] = make_vx_at_age(age)

    # Create sims for all scenarios × seeds
    sims = []
    sim_keys = []  # (name, seed) tuples
    for name, interventions in scenarios.items():
        for seed in range(n_seeds):
            sims.append(make_sim(name, interventions, calib_pars, seed=seed))
            sim_keys.append((name, seed))

    # Run in parallel
    print(f'Running {len(sims)} sims ({len(scenarios)} scenarios × {n_seeds} seeds) in parallel...')
    sim_results = sc.parallelize(run_one_sim, iterarg=sims)

    # Average results across seeds
    results = dict()
    for name in scenarios:
        seed_results = [res for (n, _), res in zip(sim_keys, sim_results) if n == name]
        averaged = dict()
        for age in vx_ages:
            averaged[f'cancers_{age}'] = np.mean([r[f'cancers_{age}'] for r in seed_results])
            averaged[f'n_cohort_{age}'] = np.mean([r[f'n_cohort_{age}'] for r in seed_results])
        results[name] = averaged

    sc.saveobj(f'raw_results/vx_efficacy_{location}.obj', results)
    print(f'Results saved to raw_results/vx_efficacy_{location}.obj')

    import utils as ut
    ut.vx_efficacy_to_csv(results, f'results/vx_efficacy_{location}.csv', location=location)
    print(f'Plot-ready CSV: results/vx_efficacy_{location}.csv')
    return results


def plot_results(results=None):
    """Plot and print results."""
    import matplotlib.pyplot as plt
    import utils as ut

    if results is None:
        results = sc.loadobj(f'raw_results/vx_efficacy_{location}.obj')

    ut.set_font(14)

    print('\n' + '=' * 80)
    print('VACCINE EFFICACY BY AGE AT VACCINATION')
    print(f'Location: {location.capitalize()}, 100% coverage, no prior vaccination')
    print('=' * 80)

    print(f'\n{"Age":>5} | {"Baseline (per 100k)":>20} | {"Vaxed (per 100k)":>18} | {"% reduction":>12} | {"N cohort":>10}')
    print('-' * 85)

    ages_plot = []
    baseline_rates = []
    vax_rates = []
    pct_reductions = []

    for age in vx_ages:
        # Use baseline sim for baseline rate
        base_cancers = results['Baseline'][f'cancers_{age}']
        base_n = results['Baseline'][f'n_cohort_{age}']

        # Use vax sim for vax rate (same cohort tracked by fixed indices)
        vax_cancers = results[f'Vax at {age}'][f'cancers_{age}']
        vax_n = results[f'Vax at {age}'][f'n_cohort_{age}']

        rate_baseline = 1e5 * base_cancers / base_n if base_n > 0 else 0
        rate_vax = 1e5 * vax_cancers / vax_n if vax_n > 0 else 0
        pct_reduction = 100 * (rate_baseline - rate_vax) / rate_baseline if rate_baseline > 0 else 0

        print(f'{age:>5} | {rate_baseline:>20,.0f} | {rate_vax:>18,.0f} | {pct_reduction:>11.1f}% | {base_n:>10,.0f}')

        ages_plot.append(age)
        baseline_rates.append(rate_baseline)
        vax_rates.append(rate_vax)
        pct_reductions.append(pct_reduction)

    # Plot
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

    x = np.arange(len(ages_plot))
    width = 0.35

    ax1.bar(x - width/2, baseline_rates, width, label='No vaccination', color='#e74c3c', alpha=0.7)
    ax1.bar(x + width/2, vax_rates, width, label='With vaccination', color='#3498db', alpha=0.7)
    ax1.set_xlabel('Age at vaccination')
    ax1.set_ylabel('Lifetime cancers per 100k women')
    ax1.set_title('Lifetime cancer rate:\nbaseline vs vaccinated')
    ax1.set_xticks(x)
    ax1.set_xticklabels(ages_plot)
    ax1.legend()
    ax1.grid(alpha=0.3, axis='y')
    sc.SIticks(ax1)

    ax2.bar(x, pct_reductions, color='#27ae60', alpha=0.7, edgecolor='black')
    for xi, val in zip(x, pct_reductions):
        ax2.text(xi, val + 1, f'{val:.0f}%', ha='center', va='bottom', fontweight='bold')
    ax2.set_xlabel('Age at vaccination')
    ax2.set_ylabel('% reduction in lifetime cancers')
    ax2.set_title('Vaccine efficacy by\nage at vaccination')
    ax2.set_xticks(x)
    ax2.set_xticklabels(ages_plot)
    ax2.set_ylim(0, 105)
    ax2.grid(alpha=0.3, axis='y')

    fig.suptitle(f'HPV vaccine impact by age at vaccination — {location.capitalize()}\n(100% coverage, unvaccinated population)')
    fig.tight_layout()
    sc.savefig(f'figures/vx_efficacy_by_age_{location}.png', dpi=150)
    print(f'\nFigure saved to figures/vx_efficacy_by_age_{location}.png')

    return fig


if __name__ == '__main__':
    results = run_sims()
    plot_results(results)
