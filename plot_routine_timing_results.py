'''
Plot results from routine timing sensitivity analysis

This script analyzes actual HPVsim results showing how catch-up campaign
effectiveness depends on when routine vaccination started.
'''

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sciris as sc


def load_and_summarize_results(location='kenya'):
    """
    Load results from routine timing analysis and summarize key metrics
    """

    # Load results
    results_file = f'raw_results/routine_timing_analysis_{location}.obj'
    try:
        all_results = sc.loadobj(results_file)
    except FileNotFoundError:
        print(f'Error: {results_file} not found.')
        print('Please run run_routine_timing_scenarios.py first.')
        return None

    # Extract key information
    summary = []

    for routine_key, scenarios in all_results.items():
        # Extract routine start year from key (e.g., 'routine_2019' -> 2019)
        routine_start = int(routine_key.split('_')[1])

        # Assume catch-up happens in 2026
        catchup_year = 2026
        max_routine_age = 10 + (catchup_year - routine_start)

        # Extract baseline cancers
        baseline_cancers = scenarios['Baseline']['cancers'][-1]  # Total by 2100

        # Loop through catch-up scenarios
        for scen_label, scen_data in scenarios.items():
            if scen_label == 'Baseline':
                continue

            # Extract age range from scenario name (e.g., 'Catchup 10-30' -> (10, 30))
            age_str = scen_label.split('Catchup ')[1]
            lower_age, upper_age = map(int, age_str.split('-'))

            # Calculate cancers with catch-up
            catchup_cancers = scen_data['cancers'][-1]

            # Calculate cancers averted
            cancers_averted = baseline_cancers - catchup_cancers

            # Calculate overlap with routine vaccination
            overlap_lower = max(lower_age, 10)
            overlap_upper = min(upper_age, max_routine_age)
            if overlap_upper >= overlap_lower:
                n_overlap = overlap_upper - overlap_lower + 1
            else:
                n_overlap = 0
            n_total = upper_age - lower_age + 1
            pct_overlap = 100 * n_overlap / n_total

            summary.append({
                'routine_start': routine_start,
                'max_routine_age': max_routine_age,
                'catchup_lower': lower_age,
                'catchup_upper': upper_age,
                'catchup_range': f'{lower_age}-{upper_age}',
                'baseline_cancers': baseline_cancers,
                'catchup_cancers': catchup_cancers,
                'cancers_averted': cancers_averted,
                'pct_reduction': 100 * cancers_averted / baseline_cancers,
                'pct_overlap': pct_overlap,
                'pct_eligible': 100 - pct_overlap,
            })

    df = pd.DataFrame(summary)
    return df


def plot_results(df, location='kenya'):
    """
    Create comprehensive visualization of routine timing sensitivity
    """

    if df is None or len(df) == 0:
        print('No data to plot.')
        return

    # Set up plotting
    sc.options(dpi=150)
    sns.set_style('whitegrid')

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # Get unique catch-up ranges
    catchup_ranges = df['catchup_range'].unique()
    colors = sns.color_palette('husl', n_colors=len(catchup_ranges))
    color_map = dict(zip(catchup_ranges, colors))

    # Panel A: Cancers averted vs routine start year
    ax1 = axes[0, 0]
    for catchup_range in catchup_ranges:
        subset = df[df['catchup_range'] == catchup_range]
        ax1.plot(subset['routine_start'], subset['cancers_averted'],
                marker='o', linewidth=2, markersize=8,
                color=color_map[catchup_range], label=f'Ages {catchup_range}')

    ax1.set_xlabel('Year routine vaccination started', fontsize=11, fontweight='bold')
    ax1.set_ylabel('Cancers averted by 2100', fontsize=11, fontweight='bold')
    ax1.set_title('A. Marginal benefit of catch-up campaigns', fontsize=12, fontweight='bold')
    ax1.legend(title='Catch-up age range', fontsize=9)
    ax1.grid(alpha=0.3)
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)

    # Panel B: Percent reduction in cancers
    ax2 = axes[0, 1]
    for catchup_range in catchup_ranges:
        subset = df[df['catchup_range'] == catchup_range]
        ax2.plot(subset['routine_start'], subset['pct_reduction'],
                marker='s', linewidth=2, markersize=8,
                color=color_map[catchup_range], label=f'Ages {catchup_range}')

    ax2.set_xlabel('Year routine vaccination started', fontsize=11, fontweight='bold')
    ax2.set_ylabel('% reduction in cancers', fontsize=11, fontweight='bold')
    ax2.set_title('B. Relative impact of catch-up campaigns', fontsize=12, fontweight='bold')
    ax2.legend(title='Catch-up age range', fontsize=9)
    ax2.grid(alpha=0.3)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)

    # Panel C: % of catch-up target eligible (not already vaccinated)
    ax3 = axes[1, 0]
    for catchup_range in catchup_ranges:
        subset = df[df['catchup_range'] == catchup_range]
        ax3.plot(subset['routine_start'], subset['pct_eligible'],
                marker='^', linewidth=2, markersize=8,
                color=color_map[catchup_range], label=f'Ages {catchup_range}')

    ax3.set_xlabel('Year routine vaccination started', fontsize=11, fontweight='bold')
    ax3.set_ylabel('% of target NOT vaccinated\nby routine program', fontsize=11, fontweight='bold')
    ax3.set_title('C. Fraction of population eligible for catch-up', fontsize=12, fontweight='bold')
    ax3.set_ylim([0, 105])
    ax3.legend(title='Catch-up age range', fontsize=9)
    ax3.grid(alpha=0.3)
    ax3.spines['top'].set_visible(False)
    ax3.spines['right'].set_visible(False)

    # Panel D: Efficiency (cancers averted per % of population targeted)
    ax4 = axes[1, 1]
    for catchup_range in catchup_ranges:
        subset = df[df['catchup_range'] == catchup_range]
        # Efficiency: cancers averted per eligible person
        target_size = subset['catchup_upper'] - subset['catchup_lower'] + 1
        efficiency = subset['cancers_averted'] / (target_size * subset['pct_eligible'] / 100)
        ax4.plot(subset['routine_start'], efficiency,
                marker='D', linewidth=2, markersize=8,
                color=color_map[catchup_range], label=f'Ages {catchup_range}')

    ax4.set_xlabel('Year routine vaccination started', fontsize=11, fontweight='bold')
    ax4.set_ylabel('Cancers averted per\neligible person in target', fontsize=11, fontweight='bold')
    ax4.set_title('D. Efficiency among eligible population', fontsize=12, fontweight='bold')
    ax4.legend(title='Catch-up age range', fontsize=9)
    ax4.grid(alpha=0.3)
    ax4.spines['top'].set_visible(False)
    ax4.spines['right'].set_visible(False)

    plt.suptitle(f'Routine Vaccination Timing Sensitivity Analysis - {location.title()}',
                fontsize=14, fontweight='bold', y=0.995)
    plt.tight_layout()
    plt.savefig(f'figures/routine_timing_sensitivity_{location}.png', dpi=150, bbox_inches='tight')
    print(f'Saved figure: figures/routine_timing_sensitivity_{location}.png')

    return fig


def create_summary_table(df, location='kenya'):
    """
    Create a summary table showing key metrics
    """

    if df is None or len(df) == 0:
        print('No data to summarize.')
        return

    # Pivot to show cancers averted for each combination
    pivot = df.pivot_table(
        values='cancers_averted',
        index='routine_start',
        columns='catchup_range',
        aggfunc='first'
    )

    print('\n' + '='*70)
    print(f'SUMMARY: Cancers Averted by 2100 - {location.title()}')
    print('='*70)
    print(pivot.to_string())
    print('='*70)

    # Calculate and show marginal benefit decline
    print('\nMarginal Benefit Decline (comparing latest vs earliest routine start):')
    print('='*70)

    for catchup_range in pivot.columns:
        earliest = pivot[catchup_range].iloc[0]
        latest = pivot[catchup_range].iloc[-1]
        decline = earliest - latest
        pct_decline = 100 * decline / earliest

        print(f'{catchup_range}: {decline:.1f} fewer cancers averted ({pct_decline:.1f}% decline)')

    print('='*70 + '\n')

    # Save to CSV
    pivot.to_csv(f'results/routine_timing_cancers_averted_{location}.csv')
    print(f'Saved summary table to results/routine_timing_cancers_averted_{location}.csv\n')

    return pivot


if __name__ == '__main__':

    location = 'kenya'

    print(f'Loading and analyzing routine timing sensitivity results for {location}...\n')

    # Load and summarize
    df = load_and_summarize_results(location=location)

    if df is not None:
        # Create summary table
        pivot = create_summary_table(df, location=location)

        # Create plots
        fig = plot_results(df, location=location)

        # Show plots
        plt.show()
    else:
        print('\nNo results found. To generate results, run:')
        print('  python run_routine_timing_scenarios.py')
