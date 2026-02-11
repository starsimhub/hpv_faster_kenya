"""
Illustrate how catch-up campaign effectiveness depends on when routine vaccination started

This creates a synthetic data figure showing that the marginal benefit of catch-up
vaccination campaigns decreases as the routine vaccination program becomes more established.
"""

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import sciris as sc

# Set up plotting
sc.options(dpi=150)
sns.set_style('whitegrid')

def synthetic_marginal_benefit():
    """
    Generate synthetic data showing declining marginal benefit of catch-up campaigns
    as routine vaccination becomes more established.

    Key assumptions for synthetic data:
    - The benefit of catch-up vaccination for an age cohort decreases if they were
      eligible for routine vaccination
    - The decrease is proportional to the coverage achieved by routine vaccination
    - Test-and-treat benefit is less affected by prior vaccination coverage
    """

    # Define scenarios: when did routine vaccination start?
    routine_start_years = np.array([2005, 2010, 2015, 2019, 2022, 2025])
    catchup_year = 2026

    # Calculate max age already vaccinated by routine program
    max_routine_ages = 10 + (catchup_year - routine_start_years)

    # Define catch-up campaign age ranges to evaluate
    catchup_ranges = [
        (10, 20),
        (10, 30),
        (10, 40),
        (10, 50),
    ]

    # Assume routine vaccination coverage scales up from 25% in year 1 to 90% by year 6
    def get_routine_coverage_by_age(age, max_routine_age):
        """Get approximate routine coverage for a given age cohort"""
        if age > max_routine_age:
            return 0  # Too old, never eligible for routine

        # Years since they were vaccinated
        years_since = max_routine_age - age

        # Simple scale-up model: 25% in year 1, increasing to 90% by year 6+
        if years_since == 0:
            cov = 0.25
        elif years_since <= 5:
            cov = 0.25 + (0.65 * years_since / 5)
        else:
            cov = 0.90

        return cov

    # Calculate marginal benefit for each scenario
    results = {}

    for catchup_range in catchup_ranges:
        lower_age, upper_age = catchup_range
        range_label = f'{lower_age}-{upper_age}'

        # For each routine start year, calculate expected benefit
        cancers_averted = []
        pct_eligible = []  # Percent of target population that wasn't already vaccinated

        for routine_start, max_routine_age in zip(routine_start_years, max_routine_ages):

            # Calculate overlap between routine-vaccinated ages and catch-up target
            overlap_ages = range(max(lower_age, 10), min(upper_age, int(max_routine_age)) + 1)

            # Estimate fraction of catch-up target already vaccinated
            n_target = upper_age - lower_age + 1
            n_already_vaccinated = 0

            for age in overlap_ages:
                routine_cov = get_routine_coverage_by_age(age, max_routine_age)
                n_already_vaccinated += routine_cov

            pct_already_vaccinated = n_already_vaccinated / n_target
            pct_eligible_for_catchup = 1 - pct_already_vaccinated
            pct_eligible.append(pct_eligible_for_catchup)

            # Synthetic benefit model:
            # - Base benefit is proportional to age range (more people = more potential impact)
            # - Adjusted by fraction who weren't already vaccinated
            # - Younger ages have higher per-person benefit (longer life remaining)

            # Weight by age (younger = higher weight)
            age_weights = np.array([1.0 / (age - 5) for age in range(lower_age, upper_age + 1)])
            age_weights = age_weights / age_weights.sum()

            # Base benefit (cancers averted per 1000 vaccinated, hypothetically)
            base_benefit = age_weights.sum() * 100  # arbitrary scale

            # Adjust for those already vaccinated
            adjusted_benefit = base_benefit * pct_eligible_for_catchup

            # Add some noise for realism
            noise = np.random.normal(1.0, 0.05)
            cancers_averted.append(adjusted_benefit * noise)

        results[range_label] = {
            'cancers_averted': np.array(cancers_averted),
            'pct_eligible': np.array(pct_eligible),
        }

    return routine_start_years, catchup_ranges, results


def plot_marginal_benefit_by_routine_start():
    """Plot showing how marginal benefit decreases as routine vaccination is more established"""

    np.random.seed(42)  # For reproducible synthetic data
    routine_start_years, catchup_ranges, results = synthetic_marginal_benefit()

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Panel A: Cancers averted by routine start year
    ax1 = axes[0]
    colors = sns.color_palette('husl', n_colors=len(catchup_ranges))

    for idx, (lower, upper) in enumerate(catchup_ranges):
        range_label = f'{lower}-{upper}'
        data = results[range_label]

        ax1.plot(routine_start_years, data['cancers_averted'],
                marker='o', linewidth=2, markersize=8,
                color=colors[idx], label=f'Ages {range_label}')

    ax1.set_xlabel('Year routine vaccination started', fontsize=12, fontweight='bold')
    ax1.set_ylabel('Cancers averted by 2100\n(per 1000 catch-up doses)', fontsize=12, fontweight='bold')
    ax1.set_title('A. Marginal benefit decreases with earlier routine start',
                 fontsize=13, fontweight='bold', pad=15)
    ax1.legend(title='Catch-up age range', fontsize=10, title_fontsize=11)
    ax1.grid(alpha=0.3)
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)

    # Add annotation
    ax1.annotate('Earlier routine vaccination\n→ More cohorts already protected\n→ Lower marginal benefit of catch-up',
                xy=(2010, 5), xytext=(2008, 15),
                fontsize=9, style='italic', alpha=0.7,
                bbox=dict(boxstyle='round,pad=0.5', facecolor='wheat', alpha=0.3),
                arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0.3', alpha=0.5))

    # Panel B: % of target population eligible (not already vaccinated)
    ax2 = axes[1]

    for idx, (lower, upper) in enumerate(catchup_ranges):
        range_label = f'{lower}-{upper}'
        data = results[range_label]

        ax2.plot(routine_start_years, data['pct_eligible'] * 100,
                marker='s', linewidth=2, markersize=8,
                color=colors[idx], label=f'Ages {range_label}')

    ax2.set_xlabel('Year routine vaccination started', fontsize=12, fontweight='bold')
    ax2.set_ylabel('% of catch-up target NOT already\nvaccinated by routine program', fontsize=12, fontweight='bold')
    ax2.set_title('B. Fraction of population still eligible for catch-up',
                 fontsize=13, fontweight='bold', pad=15)
    ax2.set_ylim([0, 105])
    ax2.legend(title='Catch-up age range', fontsize=10, title_fontsize=11)
    ax2.grid(alpha=0.3)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)

    # Add annotation
    ax2.axhline(y=50, color='red', linestyle='--', alpha=0.3, linewidth=1)
    ax2.annotate('50% threshold', xy=(2025, 50), xytext=(2022, 35),
                fontsize=9, color='red', alpha=0.7,
                arrowprops=dict(arrowstyle='->', color='red', alpha=0.5))

    plt.tight_layout()
    plt.savefig('figures/synthetic_routine_timing_sensitivity.png', dpi=150, bbox_inches='tight')
    print('Saved figure: figures/synthetic_routine_timing_sensitivity.png')

    return fig


def plot_efficiency_frontier():
    """
    Alternative visualization: efficiency frontier showing optimal catch-up strategy
    depends on routine vaccination coverage
    """

    np.random.seed(42)

    # Simulate different scenarios
    routine_scenarios = {
        'Very early (2005)': {'max_age': 31, 'color': '#d62728'},
        'Early (2010)': {'max_age': 26, 'color': '#ff7f0e'},
        'Moderate (2015)': {'max_age': 21, 'color': '#2ca02c'},
        'Current (2019)': {'max_age': 17, 'color': '#1f77b4'},
        'Recent (2022)': {'max_age': 14, 'color': '#9467bd'},
    }

    fig, ax = plt.subplots(1, 1, figsize=(10, 7))

    # For each routine scenario, show efficiency of different catch-up age ranges
    for scenario_label, scenario_data in routine_scenarios.items():
        max_routine_age = scenario_data['max_age']
        color = scenario_data['color']

        catchup_upper_ages = np.arange(15, 61, 5)
        efficiency = []
        doses_needed = []

        for upper_age in catchup_upper_ages:
            lower_age = 10
            n_target = upper_age - lower_age + 1

            # Estimate fraction not already vaccinated
            overlap = min(upper_age, max_routine_age) - 10 + 1
            overlap = max(0, overlap)

            # Assume 70% average coverage in overlapping ages
            n_already_vaccinated = overlap * 0.7
            n_eligible = n_target - n_already_vaccinated

            # Doses needed (assuming 70% catch-up coverage)
            doses = n_eligible * 0.7 * 1000  # per 1000 women
            doses_needed.append(doses / 1000)  # back to per-woman scale

            # Efficiency: cancers averted per dose (synthetic)
            # Higher for younger, decreases with age and overlap
            base_efficiency = 0.05 * (1 - (upper_age - 20) / 100)
            overlap_penalty = 1 - (n_already_vaccinated / n_target) * 0.5
            eff = base_efficiency * overlap_penalty * np.random.normal(1.0, 0.05)
            efficiency.append(eff)

        ax.plot(doses_needed, efficiency, marker='o', linewidth=2.5, markersize=7,
               color=color, label=scenario_label, alpha=0.8)

    ax.set_xlabel('Catch-up doses per woman (population-averaged)', fontsize=12, fontweight='bold')
    ax.set_ylabel('Cancers averted per dose', fontsize=12, fontweight='bold')
    ax.set_title('Efficiency frontier: Earlier routine vaccination shifts optimal catch-up strategy',
                fontsize=13, fontweight='bold', pad=15)
    ax.legend(title='Routine vaccination start', fontsize=10, title_fontsize=11, loc='upper right')
    ax.grid(alpha=0.3)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # Add annotation
    ax.annotate('More established routine programs\n→ Need fewer catch-up doses\n→ Focus on older ages',
               xy=(0.35, 0.035), xytext=(0.5, 0.045),
               fontsize=9, style='italic', alpha=0.7,
               bbox=dict(boxstyle='round,pad=0.5', facecolor='lightblue', alpha=0.3),
               arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0.3', alpha=0.5))

    plt.tight_layout()
    plt.savefig('figures/synthetic_efficiency_frontier_routine_timing.png', dpi=150, bbox_inches='tight')
    print('Saved figure: figures/synthetic_efficiency_frontier_routine_timing.png')

    return fig


if __name__ == '__main__':

    print('Generating synthetic data figures...\n')

    # Figure 1: Direct relationship between routine start year and marginal benefit
    fig1 = plot_marginal_benefit_by_routine_start()

    # Figure 2: Efficiency frontier
    fig2 = plot_efficiency_frontier()

    print('\nFigures illustrate key insight:')
    print('- Earlier routine vaccination → more cohorts already protected')
    print('- Marginal benefit of catch-up campaigns decreases')
    print('- Optimal catch-up strategy shifts to older ages when routine is well-established')

    plt.show()
