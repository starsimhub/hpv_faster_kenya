"""
Plot residual burden
"""


import pylab as pl
import sciris as sc
import numpy as np
import utils as ut

 
def make_table1():
    import pandas as pd

    # Load catch-up vaccination scenario data
    msim_dict = sc.loadobj('results/scens_kenya.obj')

    # Define scenarios and cohorts
    scenarios = ['Baseline', 'Catch-up to age 15', 'Catch-up to age 20', 'Catch-up to age 25',
                 'Catch-up to age 30', 'Catch-up to age 35', 'Catch-up to age 40']
    cohorts = ['10_15', '15_20', '20_25', '25_30', '30_35', '35_40']  #, '40_45', '45_50']
    cohort_labels = ['10-15', '15-20', '20-25', '25-30', '30-35', '35-40']  #, '40-45', '45-50']

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

    # Create DataFrame for table
    df = pd.DataFrame(results).T
    df.index.name = 'Birth Cohort'

    # Save to CSV
    df.to_csv('results/catchup_vax_cohort_cancers.csv')

    # Display table
    print("\n" + "="*80)
    print("Cumulative cancers by birth cohort (2025-2100)")
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
    bottom = np.zeros(len(scenarios))

    # Create color palette for cohorts
    cohort_colors = sc.vectocolor(len(cohort_labels)).tolist()

    for cohort_idx, cohort_label in enumerate(cohort_labels):
        values = [df_plot.loc[cohort_label, scenario] for scenario in scenarios]
        ax.bar(x, values, width=0.7, bottom=bottom, color=cohort_colors[cohort_idx],
               label=cohort_label)
        bottom += values

    ax.set_xticks(x)
    ax.set_xticklabels(scenarios, rotation=45, ha='right')
    ax.set_ylabel('Cumulative cancers (2025-2100)')
    ax.set_title('Cumulative cancers by birth cohort and catch-up vaccination strategy')
    ax.legend(title='Birth Cohort', loc='upper right', frameon=False, ncol=2)
    sc.SIticks()

    fig.tight_layout()
    fig_name = 'figures/catchup_vax_cohorts.png'
    sc.savefig(fig_name, dpi=100)

    print(f"\nFigure saved to: {fig_name}")
    return


# %% Run as a script
if __name__ == '__main__':

    location = 'kenya'
    make_table1()

    msim_dict = sc.loadobj('results/scens_kenya.obj')

