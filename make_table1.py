"""
Plot residual burden
"""


import pylab as pl
import sciris as sc
import numpy as np
import utils as ut

 
def make_table1(location):
    import pandas as pd

    # Load catch-up vaccination scenario data
    msim_dict = sc.loadobj(f'results/scens_{location}.obj')

    # Define scenarios and cohorts
    scenarios = ['Baseline', 'Catch-up to age 15', 'Catch-up to age 20', 'Catch-up to age 25',
                 'Catch-up to age 30', 'Catch-up to age 35', 'Catch-up to age 40',
                  'Catch-up to age 45', 'Catch-up to age 50',
                  'Catch-up to age 55', 'Catch-up to age 60']
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
    ax.set_title(f'Cumulative cancers by birth cohort and catch-up vaccination strategy - {location.capitalize()}')
    ax.legend(title='Birth Cohort', loc='upper right', frameon=False, ncol=2)
    sc.SIticks()

    fig.tight_layout()
    fig_name = f'figures/catchup_vax_cohorts_{location}.png'
    sc.savefig(fig_name, dpi=100)
    return


def make_single_bar(location):
    import pandas as pd

    # Load catch-up vaccination scenario data
    msim_dict = sc.loadobj(f'results/scens_{location}.obj')

    # Define age groups and corresponding scenarios
    age_groups = {
        'Up to 25': 'Catch-up to age 25',
        '25-30': 'Catch-up to age 30',
        '30-45': 'Catch-up to age 45',
        '45+': 'Catch-up to age 60'
    }

    # Get baseline cancers
    baseline_cancers = msim_dict['Baseline']['cancers'].values[sc.findinds(msim_dict['Baseline'].year, 2025)[0]:].sum()

    # Calculate cancers for each scenario
    cancers_by_scenario = {}
    for label, scenario in age_groups.items():
        cancers_by_scenario[label] = msim_dict[scenario]['cancers'].values[sc.findinds(msim_dict[scenario].year, 2025)[0]:].sum()

    # Calculate incremental benefit for each age group
    incremental_benefits = {}

    # Up to 25: benefit from baseline to age 25
    incremental_benefits['Up to 25'] = baseline_cancers - cancers_by_scenario['Up to 25']

    # 25-30: additional benefit from age 25 to age 30
    incremental_benefits['25-30'] = cancers_by_scenario['Up to 25'] - cancers_by_scenario['25-30']

    # 30-45: additional benefit from age 30 to age 45
    incremental_benefits['30-45'] = cancers_by_scenario['25-30'] - cancers_by_scenario['30-45']

    # 45+: additional benefit from age 45 to age 60
    incremental_benefits['45+'] = cancers_by_scenario['30-45'] - cancers_by_scenario['45+']

    # Calculate total benefit
    total_benefit = sum(incremental_benefits.values())

    # Calculate percentages
    percentages = {label: (benefit / total_benefit) * 100 for label, benefit in incremental_benefits.items()}

    print("\nIncremental benefits by age group:")
    for label, benefit in incremental_benefits.items():
        print(f"{label}: {benefit:.1f} cancers averted ({percentages[label]:.1f}%)")
    print(f"Total benefit: {total_benefit:.1f} cancers averted")

    ######################################################
    # Create simplified stacked bar chart (vertical)
    ######################################################
    ut.set_font(20)
    fig, ax = pl.subplots(figsize=(6, 10))

    # Define colors for age groups
    age_colors = sc.vectocolor(4).tolist()

    # Create single stacked bar
    bar_width = 0.5
    x_pos = 0

    bottom = 0
    labels_list = list(percentages.keys())
    values_list = list(percentages.values())

    for idx, (label, pct) in enumerate(zip(labels_list, values_list)):
        ax.bar(x_pos, pct, bottom=bottom, width=bar_width,
               color=age_colors[idx], label=label)

        # Add text label showing just percentage
        center = bottom + pct / 2
        ax.text(x_pos, center,
                f'{pct:.1f}%',
                ha='center', va='center', fontsize=14, fontweight='bold',
                color='white' if idx < 2 else 'black')

        bottom += pct

    ax.set_ylim(0, 109)
    ax.set_xlim(-0.5, 1)
    ax.set_ylabel('Percentage of total vaccination benefit (%)')
    ax.set_title('Distribution of vaccination benefit\nby age group - Kenya\nCancers averted, 2025-2100')
    ax.set_xticks([])

    # Add legend
    ax.legend(title='Age group', loc='upper right', frameon=False, fontsize=12)

    fig.tight_layout()
    fig_name = 'figures/catchup_vax_benefit_distribution.png'
    sc.savefig(fig_name, dpi=100)

    print(f"\nFigure saved to: {fig_name}")
    return


# %% Run as a script
if __name__ == '__main__':

    location = 'nigeria'
    # make_table1(location)
    make_single_bar(location)

    msim_dict = sc.loadobj(f'results/scens_{location}.obj')

