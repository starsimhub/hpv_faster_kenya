"""
Plot residual burden
"""


import pylab as pl
import sciris as sc
import numpy as np
import utils as ut
import pandas as pd
import matplotlib.pyplot as plt


 
def plot_stacked_bars(location, coverage=90, add_tt=False):
    import pandas as pd

    # Load catch-up vaccination scenario data
    msim_dict = sc.loadobj(f'raw_results/scens_{location}_{coverage}.obj')

    # Define scenarios and cohorts
    scenarios = ['Baseline', 'Catch-up 10-15', 'Catch-up 10-20', 'Catch-up 10-25',
                 'Catch-up 10-30', 'Catch-up 10-35', 'Catch-up 10-40',
                  'Catch-up 10-45', 'Catch-up 10-50',
                  'Catch-up 10-55', 'Catch-up 10-60']
    scenarios = ['Baseline'] + [s + ': V' for s in scenarios if s != 'Baseline']
    if add_tt:
        scenarios = ['Baseline'] + [s + ': TTV' for s in scenarios if s != 'Baseline']
    scen_labels = ['Baseline', 'CU\n10-15', 'CU\n10-20', 'CU\n10-25', 'CU\n10-30', 'CU\n10-35', 'CU\n10-40',
                    'CU\n10-45', 'CU\n10-50',
                    'CU\n10-55', 'CU\n10-60']
    cohorts = ['10_15', '15_20', '20_25', '25_30', '30_35', '35_40', '40_45', '45_50', '50_55', '55_60']
    cohort_labels = ['10-15', '15-20', '20-25', '25-30', '30-35', '35-40', '40-45', '45-50', '50-55', '55-60']

    # Mapping from cohort to catch-up scenario index where they're first included
    cohort_to_scenario = {
        '10-15': 1,  # Catch-up to age 15
        '15-20': 2,  # Catch-up to age 20
        '20-25': 3,  # Catch-up to age 25
        '25-30': 4,  # Catch-up to age 30
        '30-35': 5,  # Catch-up to age 35
        '35-40': 6,  # Catch-up to age 40
        '40-45': 7,  # Catch-up to age 45
        '45-50': 8,  # Catch-up to age 50
        '50-55': 9,  # Catch-up to age 55
        '55-60': 10, # Catch-up to age 60
    }

    # Create dictionary to store results
    results = {cohort_label: {} for cohort_label in cohort_labels}
    results_values = {cohort_label: {} for cohort_label in cohort_labels}  # For plotting

    # Extract data for each scenario and cohort
    for scenario in scenarios:
        for cohort, cohort_label in zip(cohorts, cohort_labels):
            key = f'cohort_cancers_{cohort}'
            key_low = f'cohort_cancers_{cohort}_low'
            key_high = f'cohort_cancers_{cohort}_high'

            if key in msim_dict[scenario]:
                val = msim_dict[scenario][key]
                low = msim_dict[scenario][key_low]
                high = msim_dict[scenario][key_high]

                # Format as "value (low, high)" for table
                results[cohort_label][scenario] = f"{val:.1f} ({low:.1f}, {high:.1f})"
                results_values[cohort_label][scenario] = val
                print(f'{scenario} - Cohort {cohort_label}: {val:.1f} ({low:.1f}, {high:.1f})')
            else:
                results[cohort_label][scenario] = 'N/A'
                results_values[cohort_label][scenario] = 0

    # Extract vaccination data
    vaccinations = {}
    for scenario in scenarios:
        if 'vaccinations' in msim_dict[scenario]:
            vaccinations[scenario] = msim_dict[scenario]['vaccinations'].sum()
            print(f'{scenario}: {vaccinations[scenario]} vaccinations')
        else:
            vaccinations[scenario] = 0

    # Create DataFrame for table
    df = pd.DataFrame(results).T
    df.index.name = 'Birth Cohort'

    # Save to CSV
    df.to_csv('results/catchup_vax_cohort_cancers.csv')

    # Display table
    print("\n" + "="*80)
    print(f"Cumulative cancers by birth cohort (2025-2100), {location.capitalize()}")
    print("="*80)
    print(df.to_string())
    print("="*80)
    print(f"\nTable saved to: results/catchup_vax_cohort_cancers.csv")

    # Create DataFrame for plotting
    df_plot = pd.DataFrame(results_values).T

    ######################################################
    # Create stacked bar chart
    ######################################################
    ut.set_font(20)
    fig, ax = pl.subplots(figsize=(14, 8))

    # Create stacked bars
    x = np.arange(len(scenarios))
    bottom_tracker = np.zeros((len(scenarios), len(cohort_labels)))  # Track bottom positions for each segment

    # Create color palette for cohorts
    cohort_colors = sc.vectocolor(len(cohort_labels)).tolist()

    for cohort_idx, cohort_label in enumerate(cohort_labels):
        values = [df_plot.loc[cohort_label, scenario] for scenario in scenarios]
        bottom = np.sum(bottom_tracker[:, :cohort_idx], axis=1) if cohort_idx > 0 else np.zeros(len(scenarios))

        bars = ax.bar(x, values, width=0.7, bottom=bottom, color=cohort_colors[cohort_idx],
                      label=cohort_label)

        # Store the values for this cohort
        bottom_tracker[:, cohort_idx] = values

        # Add text labels on each bar section
        for bar_idx, (bar, val) in enumerate(zip(bars, values)):
            if val > 0:  # Only label non-zero values
                height = bar.get_height()
                y_pos = bottom[bar_idx] + height / 2
                ax.text(bar.get_x() + bar.get_width() / 2., y_pos,
                        f'{sc.sigfig(val, 3, SI=True)}',
                        ha='center', va='center', fontsize=9, fontweight='bold',
                        color='white' if cohort_idx < 4 else 'black')  # White text for darker colors

    # Calculate total heights
    total_heights = np.sum(bottom_tracker, axis=1)

    # Add total height labels with % reduction on top of each bar
    for bar_idx in range(len(scenarios)):
        total_height = total_heights[bar_idx]

        # Calculate percent reduction from previous bar
        if bar_idx > 0:
            prev_height = total_heights[bar_idx - 1]
            percent_reduction = ((prev_height - total_height) / prev_height) * 100
            label_text = f'{sc.sigfig(total_height, 3, SI=True)}\n(-{percent_reduction:.1f}%)'
        else:
            label_text = f'{sc.sigfig(total_height, 3, SI=True)}'

        ax.text(x[bar_idx], total_height,
                label_text,
                ha='center', va='bottom', fontsize=10, fontweight='bold',
                color='black')

    # Add arrows showing reductions when cohort is added to catch-up
    from matplotlib.patches import FancyArrowPatch
    import matplotlib.colors as mcolors

    for cohort_idx, cohort_label in enumerate(cohort_labels):
        if cohort_label not in cohort_to_scenario:
            continue

        target_scenario_idx = cohort_to_scenario[cohort_label]
        source_scenario_idx = target_scenario_idx - 1  # Previous scenario

        # Get values from previous and current scenario
        source_val = df_plot.loc[cohort_label, scenarios[source_scenario_idx]]
        target_val = df_plot.loc[cohort_label, scenarios[target_scenario_idx]]

        if source_val > 0:
            # Calculate percent reduction
            percent_reduction = ((source_val - target_val) / source_val) * 100

            # Calculate positions for arrow
            # Start: right edge of cohort segment in previous scenario
            start_bottom = np.sum(bottom_tracker[source_scenario_idx, :cohort_idx]) if cohort_idx > 0 else 0
            start_y = start_bottom + source_val / 2
            start_x = x[source_scenario_idx] + 0.35  # Right edge of bar (width=0.7)

            # End: left edge of cohort segment in target scenario
            end_bottom = np.sum(bottom_tracker[target_scenario_idx, :cohort_idx]) if cohort_idx > 0 else 0
            end_y = end_bottom + target_val / 2
            end_x = x[target_scenario_idx] - 0.35  # Left edge of bar

            # Create arrow with lower alpha
            color_rgba = mcolors.to_rgba(cohort_colors[cohort_idx], alpha=0.4)

            arrow = FancyArrowPatch((start_x, start_y), (end_x, end_y),
                                   arrowstyle='->', mutation_scale=20, linewidth=2,
                                   color=color_rgba, zorder=10)
            ax.add_patch(arrow)

            # Add text label showing percent reduction at midpoint between bars
            mid_x = (start_x + end_x) / 2
            mid_y = (start_y + end_y) / 2
            ax.text(mid_x, mid_y, f'-{percent_reduction:.0f}%',
                    ha='center', va='center', fontsize=10, fontweight='bold',
                    bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor=color_rgba, linewidth=2),
                    zorder=11)

    # Add vaccination numbers and NNV at the bottom
    ax2 = ax.twiny()
    ax2.set_xlim(ax.get_xlim())
    ax2.set_xticks(x)

    # Create labels with vaccinations and NNV
    bottom_labels = []
    for bar_idx in range(len(scenarios)):
        vax_count = vaccinations[scenarios[bar_idx]]

        if bar_idx > 0:
            # Calculate cancers averted
            cancers_averted = total_heights[bar_idx - 1] - total_heights[bar_idx]

            # Calculate additional vaccinations
            vax_added = vax_count - vaccinations[scenarios[bar_idx - 1]]

            # Calculate NNV
            if cancers_averted > 0:
                nnv = vax_added / cancers_averted
                label = f'{sc.sigfig(vax_count, 3, SI=True)}\n{int(nnv):,}'
            else:
                label = f'{sc.sigfig(vax_count, 3, SI=True)}\nN/A'
        else:
            label = f'Total doses: \nNNV: '

        bottom_labels.append(label)

    ax2.set_xticklabels(bottom_labels, fontsize=9)
    ax2.tick_params(axis='x', which='both', length=0)  # Remove tick marks
    ax2.xaxis.set_ticks_position('bottom')
    ax2.xaxis.set_label_position('bottom')
    ax2.spines['bottom'].set_position(('outward', 60))

    ax.set_xticks(x)
    ax.set_xticklabels(scen_labels)
    ax.set_ylabel('Cumulative cancers (2025-2100)')
    title = f'Cumulative cancers by birth cohort and HPV-Faster target age - {location.capitalize()}\n'
    if add_tt:
        title += 'WITH test and treat'
    else:
        title += 'WITHOUT test and treat'
    ax.set_title(title)
    ax.legend(title='Birth Cohort', loc='upper right', frameon=False, ncol=2)
    sc.SIticks()

    fig.tight_layout()
    figname = f'catchup_vax_cohorts_{location}_{coverage}' + ('_tt' if add_tt else '') + '.png'
    fig_name = f'figures/{figname}'
    sc.savefig(fig_name, dpi=100)
    return


def plot_single_bar(location, coverage=90):
    import pandas as pd

    # Load catch-up vaccination scenario data
    msim_dict = sc.loadobj(f'raw_results/scens_{location}_{coverage}.obj')

    # Define cohorts (ages in 2025)
    cohorts = ['10_15', '15_20', '20_25', '25_30', '30_35', '35_40', '40_45', '45_50', '50_55', '55_60']

    # Helper function to sum cancers across 10-60 cohorts
    def get_cohort_cancers(scenario_name):
        total = 0
        for cohort in cohorts:
            key = f'cohort_cancers_{cohort}'
            if key in msim_dict[scenario_name]:
                total += msim_dict[scenario_name][key]
        return total

    # Define scenarios to compare
    scenarios = {
        'No interventions': 'No interventions',
        'Baseline': 'Baseline',
        'Catch-up 10-15': 'Catch-up 10-15: V',
        'Catch-up 10-30': 'Catch-up 10-30: V',
        'Catch-up 10-35': 'Catch-up 10-35: V',
        'Catch-up 10-40': 'Catch-up 10-40: V',
        'Catch-up 10-45': 'Catch-up 10-45: V',
        'Catch-up 10-60': 'Catch-up 10-60: V'
    }

    # Get total cancers in 10-60 cohorts for each scenario
    cohort_cancers = {}
    for label, scenario_name in scenarios.items():
        if scenario_name in msim_dict:
            cohort_cancers[label] = get_cohort_cancers(scenario_name)
        else:
            # Skip missing scenarios
            cohort_cancers[label] = None

    # Use Baseline as the reference
    baseline_cancers = cohort_cancers['Baseline']

    # Calculate incremental benefits (absolute numbers)
    incremental_benefits = {}

    # Up to 15: benefit from vaccinating 10-15 year olds
    incremental_benefits['Up to 15'] = cohort_cancers['Baseline'] - cohort_cancers['Catch-up 10-15']

    # 15-30: additional benefit from vaccinating 15-30 year olds
    incremental_benefits['15-30'] = cohort_cancers['Catch-up 10-15'] - cohort_cancers['Catch-up 10-30']

    # 30-35: additional benefit from vaccinating 30-35 year olds
    incremental_benefits['30-35'] = cohort_cancers['Catch-up 10-30'] - cohort_cancers['Catch-up 10-35']

    # 35-40: additional benefit from vaccinating 35-40 year olds
    incremental_benefits['35-40'] = cohort_cancers['Catch-up 10-35'] - cohort_cancers['Catch-up 10-40']

    # 40-45: additional benefit from vaccinating 40-45 year olds
    incremental_benefits['40-45'] = cohort_cancers['Catch-up 10-40'] - cohort_cancers['Catch-up 10-45']

    # 45+: additional benefit from vaccinating 45-60 year olds
    incremental_benefits['45+'] = cohort_cancers['Catch-up 10-45'] - cohort_cancers['Catch-up 10-60']

    # Calculate total benefit
    total_catchup_benefit = sum(incremental_benefits.values())

    # Calculate percentages RELATIVE TO BASELINE CANCERS
    percentages = {label: (benefit / baseline_cancers) * 100 for label, benefit in incremental_benefits.items()}

    print("\n" + "="*80)
    print(f"Cancers averted in 10-60 cohorts by catch-up vaccination ({coverage}% coverage)")
    print(f"Percentages relative to Baseline scenario")
    print("="*80)
    print(f"\nCohort cancers by scenario:")
    for label, cancers in cohort_cancers.items():
        if cancers is not None:
            print(f"  {label}: {cancers:,.1f}")
    print(f"\nIncremental benefits (relative to Baseline = {baseline_cancers:,.1f} cancers):")
    for label, benefit in incremental_benefits.items():
        print(f"  {label}: {benefit:,.1f} cancers averted ({percentages[label]:.1f}% of baseline)")
    print(f"\nTotal catch-up benefit: {total_catchup_benefit:,.1f} cancers averted ({sum(percentages.values()):.1f}% of baseline)")
    print("="*80)

    ######################################################
    # Create simplified stacked bar chart (vertical)
    ######################################################
    ut.set_font(20)
    fig, ax = pl.subplots(figsize=(8, 10))

    # Define colors for intervention groups
    colors_list = sc.vectocolor(6).tolist()

    # Create single stacked bar
    bar_width = 0.5
    x_pos = 0

    bottom = 0
    labels_list = ['Up to 15', '15-30', '30-35', '35-40', '40-45', '45+']

    for idx, label in enumerate(labels_list):
        pct = percentages[label]
        ax.bar(x_pos, pct, bottom=bottom, width=bar_width,
               color=colors_list[idx], label=label)

        # Add text label showing percentage
        center = bottom + pct / 2
        ax.text(x_pos, center,
                f'{pct:.1f}%',
                ha='center', va='center', fontsize=12, fontweight='bold',
                color='white' if idx < 3 else 'black')

        bottom += pct

    ax.set_ylim(0, bottom * 1.1)  # Set y-limit slightly above the total
    ax.set_xlim(-0.5, 1)
    ax.set_ylabel('% of baseline cancers averted', fontsize=16)
    ax.set_title(f'Catch-up vaccination benefit by age group - Kenya ({coverage}% coverage)\n% of baseline cancers in 10-60 cohorts averted, 2025-2100', fontsize=16)
    ax.set_xticks([])

    # Add legend
    ax.legend(title='Age vaccinated', loc='upper right', frameon=False, fontsize=12, ncol=1)

    fig.tight_layout()
    fig_name = f'figures/catchup_vax_benefit_distribution_{coverage}.png'
    sc.savefig(fig_name, dpi=100)

    print(f"\nFigure saved to: {fig_name}")
    return


def plot_baseline(sim, data_file='data/kenya_cancer_cases.csv', 
                             years_to_plot=[2020, 2050, 2075, 2100],
                             figsize=(14, 10)):
    """
    Create a 4-panel plot showing:
    A) Cancers by age for selected years vs data
    B) Total cancers over time
    C) Female population over time
    D) Age-standardized cancer incidence over time
    
    Args:
        sim: HPVsim simulation object after running
        data_file: Path to kenya_cancer_cases.csv
        years_to_plot: List of years to show in panel A
        figsize: Figure size tuple
    """
    
    # Load the data
    try:
        cancer_data = pd.read_csv(data_file)
    except:
        print(f"Warning: Could not load {data_file}")
        cancer_data = None
    
    # Set up the figure
    fig, axes = plt.subplots(2, 2, figsize=figsize)
    fig.suptitle('HPV-FASTER Kenya: Model Results', fontsize=16, fontweight='bold')
    
    # Define age bin labels (18 bins of 5 years each: 0-4, 5-9, ..., 85+)
    age_bins = [f'{i}-{i+4}' if i < 85 else '85+' 
                for i in range(0, 90, 5)][:18]
    age_bin_centers = np.arange(2.5, 90, 5)[:18]  # Midpoints for plotting
    
    # Color scheme
    colors = plt.cm.viridis(np.linspace(0.2, 0.9, len(years_to_plot)))
    
    # ========== PANEL A: Cancers by Age ==========
    ax_a = axes[0, 0]
    
    # Get the start year and timestep info
    start_year = sim['start']
    dt = sim['dt']
    
    # Plot model results for each year
    for i, year in enumerate(years_to_plot):
        # Calculate the index for this year
        year_idx = int((year - start_year))
        
        if year_idx < sim.results['cancers_by_age'].shape[1]:
            cancers_by_age = sim.results['cancers_by_age'][:, year_idx]
            ax_a.plot(age_bin_centers, cancers_by_age, 
                     color=colors[i], linewidth=2.5, 
                     label=f'{year} (model)', marker='o', markersize=5)
    
    # Plot data 
    plot_years = cancer_data.year.unique()
    if cancer_data is not None:
        for year in plot_years:
            if 'year' in cancer_data.columns:
                year_data = cancer_data[cancer_data['year'] == year]
                if len(year_data) > 0:
                    age_mid = year_data['age']                    
                    cases = year_data['value']
                    
                    ax_a.scatter(age_mid, cases, s=100, marker='s', 
                               edgecolors='black', linewidths=1.5,
                               label=f'{year} (data)', zorder=10, alpha=0.7)
    
    ax_a.set_xlabel('Age (years)', fontsize=11, fontweight='bold')
    ax_a.set_ylabel('Number of cancers', fontsize=11, fontweight='bold')
    ax_a.set_title('A) Cancer cases by age', fontsize=12, fontweight='bold', loc='left')
    ax_a.legend(fontsize=8, ncol=2)
    ax_a.grid(alpha=0.3, linestyle='--')
    ax_a.set_xlim(-5, 95)
    
    # Get indices for year >= 2000 for panels B-D
    years = sim.results['year']
    idx_2000 = np.where(years >= 2000)[0]
    years_subset = years[idx_2000]
    
    # ========== PANEL B: Total Cancers Over Time ==========
    ax_b = axes[0, 1]
    
    total_cancers = sim.results['cancers'][idx_2000]
    
    ax_b.plot(years_subset, total_cancers, linewidth=2.5, color='#e74c3c', label='Annual cancers')
    
    ax_b.set_xlabel('Year', fontsize=11, fontweight='bold')
    ax_b.set_ylabel('Annual cancer cases', fontsize=11, fontweight='bold')
    ax_b.set_title('B) Cancer cases over time', fontsize=12, fontweight='bold', loc='left')
    ax_b.grid(alpha=0.3, linestyle='--')
    ax_b.legend(fontsize=9, loc='upper left')
    sc.SIticks(ax_b)
    
    # ========== PANEL C: Female Population Over Time ==========
    ax_c = axes[1, 0]
    
    female_pop = sim.results['n_alive_by_sex'][0, idx_2000]  # Index 0 for females
    
    ax_c.plot(years_subset, female_pop, linewidth=2.5, color='#9b59b6')
    ax_c.fill_between(years_subset, female_pop, alpha=0.3, color='#9b59b6')
    
    ax_c.set_xlabel('Year', fontsize=11, fontweight='bold')
    ax_c.set_ylabel('Number of women', fontsize=11, fontweight='bold')
    ax_c.set_title('C) Female population over time', fontsize=12, fontweight='bold', loc='left')
    ax_c.grid(alpha=0.3, linestyle='--')
    sc.SIticks(ax_c)
    
    # ========== PANEL D: Cancer Incidence Over Time ==========
    ax_d = axes[1, 1]
    
    asr_incidence = sim.results['asr_cancer_incidence'][idx_2000]
    
    # Apply smoothing using a rolling average
    from scipy.ndimage import uniform_filter1d
    window_size = 5  # Smooth over 5 years
    asr_incidence_smooth = uniform_filter1d(asr_incidence, size=window_size, mode='nearest')
    
    ax_d.plot(years_subset, asr_incidence_smooth, linewidth=2.5, color='#27ae60')
    
    # Add a horizontal line for reference (e.g., baseline year)
    baseline_year = 2020
    if baseline_year >= years_subset[0] and baseline_year <= years_subset[-1]:
        baseline_idx = np.where(years_subset == baseline_year)[0]
        if len(baseline_idx) > 0:
            baseline_value = asr_incidence_smooth[baseline_idx[0]]
            ax_d.axhline(y=baseline_value, color='gray', linestyle='--', 
                       linewidth=1.5, alpha=0.7, label=f'{baseline_year} baseline')
            ax_d.legend(fontsize=9)
    
    ax_d.set_xlabel('Year', fontsize=11, fontweight='bold')
    ax_d.set_ylabel('Age-standardized rate\n(per 100,000 women)', fontsize=11, fontweight='bold')
    ax_d.set_title('D) Cancer incidence (ASR) over time', fontsize=12, fontweight='bold', loc='left')
    ax_d.grid(alpha=0.3, linestyle='--')
    ax_d.set_ylim(bottom=0)  # Set y-axis to start at 0
    sc.SIticks(ax_d)
    
    # Adjust layout
    plt.tight_layout()

    # Print summary statistics
    print("\n=== Summary Statistics ===")
    print(f"Total cancers (2025-2100): {sim.results['cancers'][65:].sum():,.0f}")
    print(f"Female population in 2025: {sim.results['n_alive_by_sex'][0, 65]:,.0f}")
    print(f"Female population in 2100: {sim.results['n_alive_by_sex'][0, -1]:,.0f}")
    print(f"ASR incidence in 2025: {sim.results['asr_cancer_incidence'][65]:.1f}")
    print(f"ASR incidence in 2100: {sim.results['asr_cancer_incidence'][-1]:.1f}")

    # Save the figure
    plt.savefig('figures/baseline_sim.png', dpi=300, bbox_inches='tight')

    return fig, axes


def plot_cohort_simplified(sim, start_cohort_year=2025, max_cohort_age=85, figsize=(16, 8)):
    """
    Create a simplified stacked area plot showing incident cancers by age group in 2025.

    Four segments:
    - Age 0-10 in 2025: Averted through routine vaccination
    - Age 10-15 in 2025: Averted through MACs
    - Age 15-85 in 2025: Additional potential lives saved
    - Remaining (grey): Future cancers avertable through sustained routine vaccination

    Args:
        sim: HPVsim simulation object after running
        start_cohort_year: Year to define cohorts (default: 2025)
        max_cohort_age: Maximum age for cohorts to plot (default: 85)
        figsize: Figure size tuple
    """

    # Get the start year and timestep info
    start_year = sim['start']
    dt = sim['dt']

    # Get years and filter to start_cohort_year onwards
    years = sim.results['year']
    idx_start = np.where(years >= start_cohort_year)[0]
    years_subset = years[idx_start]

    # Define simplified cohort groups
    cohort_groups = {
        'routine_vax': (0, 10, 'Averted through routine vaccination'),
        'macs': (10, 15, 'Averted through MACs'),
        'additional': (15, max_cohort_age, 'Additional potential lives saved'),
    }

    # Define age bins for cancers_by_age (18 bins of 5 years: 0-4, 5-9, ..., 85+)
    age_bin_midpoints = np.arange(2.5, 90, 5)

    # Helper function to get cohort group for a given age in 2025
    def get_cohort_group(age_in_2025):
        """Get the cohort group name for a given age in 2025, or None if outside range"""
        if age_in_2025 > max_cohort_age:
            return None
        for group_name, (lower, upper, label) in cohort_groups.items():
            if lower <= age_in_2025 < upper:
                return group_name
        return None

    # Prepare data structure: group x year
    group_cancers = {group: np.zeros(len(years_subset)) for group in cohort_groups.keys()}

    # For each year after start_cohort_year
    for year_i, year in enumerate(years_subset):
        year_idx = int(year - start_year)

        # Get cancers by age for this year
        if year_idx < sim.results['cancers_by_age'].shape[1]:
            cancers_this_year = sim.results['cancers_by_age'][:, year_idx]

            # For each age bin, determine which cohort group it belongs to
            years_since_cohort = year - start_cohort_year

            for age_bin_i, age_bin_mid in enumerate(age_bin_midpoints):
                # What was this person's age in start_cohort_year?
                age_in_cohort_year = age_bin_mid - years_since_cohort

                # Only include if they were alive in the cohort year (age >= 0)
                if age_in_cohort_year >= 0:
                    group_name = get_cohort_group(age_in_cohort_year)
                    if group_name is not None:
                        group_cancers[group_name][year_i] += cancers_this_year[age_bin_i]

    # Get total cancers for calculating the grey section
    total_cancers = sim.results['cancers'][idx_start]

    # Calculate the grey section (future cancers from sustained vaccination)
    sum_three_groups = (group_cancers['routine_vax'] +
                        group_cancers['macs'] +
                        group_cancers['additional'])
    future_cancers = total_cancers - sum_three_groups

    # Create the plot
    fig, ax = plt.subplots(figsize=figsize)

    # Define colors
    colors = {
        'routine_vax': '#2ecc71',  # Green
        'macs': '#3498db',         # Blue
        'additional': '#e67e22',   # Orange
        'future': '#95a5a6'        # Grey
    }

    # Create stacked area plot from bottom to top
    # Start with bottom = 0
    bottom = np.zeros(len(years_subset))

    # Layer 1: Routine vaccination (0-10)
    ax.fill_between(years_subset, bottom, bottom + group_cancers['routine_vax'],
                    color=colors['routine_vax'], alpha=0.8,
                    label=cohort_groups['routine_vax'][2])
    bottom += group_cancers['routine_vax']

    # Layer 2: MACs (10-15)
    ax.fill_between(years_subset, bottom, bottom + group_cancers['macs'],
                    color=colors['macs'], alpha=0.8,
                    label=cohort_groups['macs'][2])
    bottom += group_cancers['macs']

    # Layer 3: Additional potential (15-85)
    ax.fill_between(years_subset, bottom, bottom + group_cancers['additional'],
                    color=colors['additional'], alpha=0.8,
                    label=cohort_groups['additional'][2])
    bottom += group_cancers['additional']

    # Layer 4: Future cancers (grey)
    ax.fill_between(years_subset, bottom, bottom + future_cancers,
                    color=colors['future'], alpha=0.6,
                    label='Future cancers avertable through\nsustained routine vaccination')

    # Add total cancers line for reference
    ax.plot(years_subset, total_cancers, linewidth=2, color='black',
            linestyle='--', alpha=0.7, label='Total incident cancers')

    ax.set_xlabel('Year', fontsize=12, fontweight='bold')
    ax.set_ylabel('Incident cancer cases', fontsize=12, fontweight='bold')
    ax.set_title(f'Incident Cervical Cancers by Age Group in {start_cohort_year}',
                 fontsize=14, fontweight='bold')
    ax.grid(alpha=0.3, linestyle='--', axis='y')

    # Add legend
    ax.legend(fontsize=10, loc='upper left', framealpha=0.9)

    sc.SIticks(ax)

    plt.tight_layout()

    # Save the figure
    plt.savefig(f'figures/cohort_simplified_{start_cohort_year}.png', dpi=300, bbox_inches='tight')

    # Calculate and print total cancers in each segment
    print(f"\n=== Simplified Cohort Summary (age groups in {start_cohort_year}) ===")
    print(f"Time period: {start_cohort_year}-{int(years_subset[-1])}\n")

    total_routine = group_cancers['routine_vax'].sum()
    total_macs = group_cancers['macs'].sum()
    total_additional = group_cancers['additional'].sum()
    total_future = future_cancers.sum()
    grand_total = total_cancers.sum()

    print(f"1. {cohort_groups['routine_vax'][2]} (age 0-10 in {start_cohort_year}):")
    print(f"   Total cancers: {total_routine:,.0f} ({100*total_routine/grand_total:.1f}%)\n")

    print(f"2. {cohort_groups['macs'][2]} (age 10-15 in {start_cohort_year}):")
    print(f"   Total cancers: {total_macs:,.0f} ({100*total_macs/grand_total:.1f}%)\n")

    print(f"3. {cohort_groups['additional'][2]} (age 15-{max_cohort_age} in {start_cohort_year}):")
    print(f"   Total cancers: {total_additional:,.0f} ({100*total_additional/grand_total:.1f}%)\n")

    print(f"4. Future cancers avertable through sustained routine vaccination:")
    print(f"   Total cancers: {total_future:,.0f} ({100*total_future/grand_total:.1f}%)\n")

    print(f"Total incident cancers ({start_cohort_year}-{int(years_subset[-1])}): {grand_total:,.0f}")

    # Export data to CSV
    df = pd.DataFrame({
        'year': years_subset,
        'routine_vax_0_10': group_cancers['routine_vax'],
        'macs_10_15': group_cancers['macs'],
        'additional_15_85': group_cancers['additional'],
        'future_cancers': future_cancers,
        'total_cancers': total_cancers
    })

    csv_filename = f'results/cohort_simplified_{start_cohort_year}.csv'
    df.to_csv(csv_filename, index=False)
    print(f"\nData exported to: {csv_filename}")

    return fig, ax, {
        'routine_vax': group_cancers['routine_vax'],
        'macs': group_cancers['macs'],
        'additional': group_cancers['additional'],
        'future': future_cancers,
        'years': years_subset
    }


def plot_cohort_decomposition(sim, start_cohort_year=2025, max_cohort_age=60, figsize=(16, 8)):
    """
    Create a plot showing incident cancers by birth cohort over time.

    Args:
        sim: HPVsim simulation object after running
        start_cohort_year: Year to define cohorts (default: 2025)
        max_cohort_age: Maximum age for cohorts to plot (default: 60)
        figsize: Figure size tuple
    """
    
    # Get the start year and timestep info
    start_year = sim['start']
    dt = sim['dt']
    
    # Get years and filter to 2000 onwards
    years = sim.results['year']
    idx_2000 = np.where(years >= 2000)[0]
    years_subset = years[idx_2000]
    
    # Define cohorts based on age in start_cohort_year
    # Age bins: 0-10, 10-15, 15-20, 20-25, ..., up to max_cohort_age
    cohort_edges = [0, 10, 15] + list(range(20, max_cohort_age + 5, 5))
    if cohort_edges[-1] < max_cohort_age:
        cohort_edges.append(max_cohort_age)
    
    cohort_labels = []
    cohort_midpoints = []
    
    for i in range(len(cohort_edges) - 1):
        lower = cohort_edges[i]
        upper = cohort_edges[i + 1]
        cohort_labels.append(f'{lower}-{upper}')
        cohort_midpoints.append((lower + upper) / 2)
    
    n_cohorts = len(cohort_labels)
    
    # Define age bins for cancers_by_age (18 bins of 5 years: 0-4, 5-9, ..., 85+)
    age_bin_edges = list(range(0, 90, 5)) + [150]
    age_bin_midpoints = np.arange(2.5, 90, 5)
    
    # Map age bins to cohorts
    def get_cohort_idx(age_in_2025):
        """Get the cohort index for a given age in 2025, or None if outside range"""
        if age_in_2025 > max_cohort_age:
            return None
        for i in range(len(cohort_edges) - 1):
            if cohort_edges[i] <= age_in_2025 < cohort_edges[i + 1]:
                return i
        return None
    
    # Prepare data structure: cohort x year
    cohort_cancers = np.zeros((n_cohorts, len(years_subset)))
    
    # For each year after start_cohort_year
    cohort_year_idx = int(start_cohort_year - start_year)
    
    for year_i, year in enumerate(years_subset):
        year_idx = int(year - start_year)
        
        # Skip if before cohort definition year
        if year < start_cohort_year:
            continue
        
        # Get cancers by age for this year
        if year_idx < sim.results['cancers_by_age'].shape[1]:
            cancers_this_year = sim.results['cancers_by_age'][:, year_idx]
            
            # For each age bin, determine which cohort it belongs to
            years_since_cohort = year - start_cohort_year
            
            for age_bin_i, age_bin_mid in enumerate(age_bin_midpoints):
                # What was this person's age in start_cohort_year?
                age_in_cohort_year = age_bin_mid - years_since_cohort
                
                # Only include if they were alive in the cohort year (age >= 0) and in plotted cohorts
                if age_in_cohort_year >= 0:
                    cohort_idx = get_cohort_idx(age_in_cohort_year)
                    if cohort_idx is not None:
                        cohort_cancers[cohort_idx, year_i] += cancers_this_year[age_bin_i]
    
    # Create the plot
    fig, ax = plt.subplots(figsize=figsize)
    
    # Total cancers line (for comparison)
    total_cancers = sim.results['cancers'][idx_2000]
    ax.plot(years_subset, total_cancers, linewidth=3, color='black', 
            label='Total incident cancers', zorder=100, linestyle='--', alpha=0.8)
    
    # Create stacked bar chart for cohorts
    # Use a colormap
    colors = sc.vectocolor(n_cohorts).tolist()
    
    # Create stacked bars
    bottom = np.zeros(len(years_subset))
    bar_width = 0.8  # Width of bars
    
    for cohort_i in range(n_cohorts):
        ax.bar(years_subset, cohort_cancers[cohort_i, :], 
               bottom=bottom, width=bar_width,
               label=f'Age {cohort_labels[cohort_i]} in {start_cohort_year}',
               color=colors[cohort_i], edgecolor='white', linewidth=0.5)
        bottom += cohort_cancers[cohort_i, :]
    
    ax.set_xlabel('Year', fontsize=12, fontweight='bold')
    ax.set_ylabel('Incident cancer cases', fontsize=12, fontweight='bold')
    ax.set_title(f'Incident Cervical Cancers by Birth Cohort (age 0-{max_cohort_age} in {start_cohort_year})', 
                 fontsize=14, fontweight='bold')
    ax.grid(alpha=0.3, linestyle='--', axis='y')
    
    # Add legend with smaller font and multiple columns
    ax.legend(fontsize=8, ncol=3, loc='upper left', framealpha=0.9)
    
    sc.SIticks(ax)
    
    plt.tight_layout()
    
    # Save the figure
    plt.savefig(f'figures/cohort_decomposition_{start_cohort_year}.png', dpi=300, bbox_inches='tight')
    
    # Print summary statistics
    print(f"\n=== Cohort Decomposition Summary (cohorts age 0-{max_cohort_age} in {start_cohort_year}) ===")
    total_in_cohorts = cohort_cancers.sum()
    total_all_cancers = total_cancers[years_subset >= start_cohort_year].sum()
    print(f"Total cancers in plotted cohorts: {total_in_cohorts:,.0f}")
    print(f"Total cancers overall ({start_cohort_year}-{int(years_subset[-1])}): {total_all_cancers:,.0f}")
    print(f"Percentage captured by plotted cohorts: {100 * total_in_cohorts / total_all_cancers:.1f}%")
    print(f"\nCancers by cohort (age in {start_cohort_year}):")
    for cohort_i in range(n_cohorts):
        cohort_total = cohort_cancers[cohort_i, :].sum()
        pct = 100 * cohort_total / total_in_cohorts if total_in_cohorts > 0 else 0
        print(f"  Age {cohort_labels[cohort_i]}: {cohort_total:,.0f} ({pct:.1f}%)")
    
    return fig, ax, cohort_cancers


# Example usage:
# fig, ax, cohort_data = plot_cohort_decomposition(sim, start_cohort_year=2025, max_cohort_age=60)
# plt.show()


# %% Run as a script
if __name__ == '__main__':

    location = 'kenya'
    do_plot_base = True
    do_plot_bars = False
    do_plot_debug = False
    coverage = 90

    if do_plot_base:
        sim = sc.loadobj(f'raw_results/sim_{location}.sim')
        plot_baseline(sim)
        plot_cohort_decomposition(sim, start_cohort_year=2025)

        # Create the new simplified cohort plot
        plot_cohort_simplified(sim, start_cohort_year=2025, max_cohort_age=85)

    if do_plot_debug:
        # Access results
        analyzer = sim.analyzers[0]
        print(analyzer.summary_df)

        # Visualize
        fig = analyzer.plot()
        plt.savefig('figures/cohort_vaccination_eligibility.png', dpi=300, bbox_inches='tight')

        # Export to CSV
        analyzer.to_csv('results/cohort_status_2025.csv')

    if do_plot_bars:
        plot_stacked_bars(location, coverage=coverage, add_tt=False)
        plot_single_bar(location, coverage=coverage)
        msim_dict = sc.loadobj(f'raw_results/scens_{location}_{coverage}.obj')

