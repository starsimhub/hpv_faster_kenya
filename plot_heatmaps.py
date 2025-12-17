"""
Plot heat maps
"""


import pylab as pl
import sciris as sc
import numpy as np
import utils as ut
import pandas as pd
import matplotlib.pyplot as plt


def plot_heatmaps(location, add_tt=False):

    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt

    # Load catch-up vaccination scenario data
    msim_dict = sc.loadobj(f'results/scens_{location}.obj')

    # Define age ranges
    lower_ages = [10, 15, 20, 25, 30, 35, 40, 45, 50, 55]
    upper_ages = [15, 20, 25, 30, 35, 40, 45, 50, 55, 60]

    # Define cohorts
    cohorts = ['10_15', '15_20', '20_25', '25_30', '30_35', '35_40', '40_45', '45_50', '50_55', '55_60']
    cohort_labels = ['10-15', '15-20', '20-25', '25-30', '30-35', '35-40', '40-45', '45-50', '50-55', '55-60']

    # Create mapping from age to cohort index
    age_to_cohort_idx = {
        10: 0, 15: 1, 20: 2, 25: 3, 30: 4, 35: 5, 40: 6, 45: 7, 50: 8, 55: 9, 60: 10
    }

    # Get baseline cancers by cohort
    baseline_cohort_cancers = {}
    for cohort in cohorts:
        key = f'cohort_cancers_{cohort}'
        if key in msim_dict['Baseline']:
            baseline_cohort_cancers[cohort] = msim_dict['Baseline'][key]

    # Helper function to sum cohort cancers from lower to upper age
    def get_cohort_cancers_sum(scenario_key, lower_age, upper_age):
        """Sum cancers across cohorts from lower_age to upper_age"""
        lower_idx = age_to_cohort_idx[lower_age]
        upper_idx = age_to_cohort_idx[upper_age]

        total_cancers = 0
        for idx in range(lower_idx, upper_idx):
            cohort = cohorts[idx]
            key = f'cohort_cancers_{cohort}'
            if key in msim_dict[scenario_key]:
                total_cancers += msim_dict[scenario_key][key]

        return total_cancers

    # Create matrices to store cancers averted for each intervention type
    n_lower = len(lower_ages)
    n_upper = len(upper_ages)

    cancers_averted_vax = np.full((n_upper, n_lower), np.nan)
    cancers_averted_tt = np.full((n_upper, n_lower), np.nan)
    cancers_averted_both = np.full((n_upper, n_lower), np.nan)

    # Populate matrices
    for i, upper in enumerate(upper_ages):
        for j, lower in enumerate(lower_ages):
            if upper > lower:  # Only valid combinations (upper triangle)
                # Get baseline cancers for this age range
                baseline_cancers_range = get_cohort_cancers_sum('Baseline', lower, upper)

                # Vaccination only
                scen_key_vax = f'Catch-up {lower}-{upper}: V'
                if scen_key_vax in msim_dict:
                    cancers = get_cohort_cancers_sum(scen_key_vax, lower, upper)
                    cancers_averted_vax[i, j] = baseline_cancers_range - cancers
                    print(f'{scen_key_vax}: {cancers_averted_vax[i, j]:.1f} cancers averted')

                # Test & treat only
                scen_key_tt = f'Catch-up {lower}-{upper}: TT'
                if scen_key_tt in msim_dict:
                    cancers = get_cohort_cancers_sum(scen_key_tt, lower, upper)
                    cancers_averted_tt[i, j] = baseline_cancers_range - cancers
                    print(f'{scen_key_tt}: {cancers_averted_tt[i, j]:.1f} cancers averted')

                # Both (combined)
                scen_key_both = f'Catch-up {lower}-{upper}: TTV'
                if scen_key_both in msim_dict:
                    cancers = get_cohort_cancers_sum(scen_key_both, lower, upper)
                    cancers_averted_both[i, j] = baseline_cancers_range - cancers
                    print(f'{scen_key_both}: {cancers_averted_both[i, j]:.1f} cancers averted')

    ######################################################
    # Create three side-by-side heatmaps
    ######################################################
    ut.set_font(12)
    fig, axes = plt.subplots(1, 3, figsize=(20, 7))

    # Define common colormap limits for consistency
    vmin = 0
    vmax = max(np.nanmax(cancers_averted_vax), np.nanmax(cancers_averted_tt), np.nanmax(cancers_averted_both))

    titles = ['Vaccination only', 'Test & Treat only', 'Combined (Vax + T&T)']
    data_matrices = [cancers_averted_vax, cancers_averted_tt, cancers_averted_both]

    for idx, (ax, title, data) in enumerate(zip(axes, titles, data_matrices)):
        # Create heatmap
        im = ax.imshow(data, cmap='YlGnBu', aspect='auto', origin='lower', vmin=vmin, vmax=vmax)

        # Set ticks
        ax.set_xticks(np.arange(n_lower))
        ax.set_yticks(np.arange(n_upper))
        ax.set_xticklabels(lower_ages)
        ax.set_yticklabels(upper_ages)

        # Labels
        ax.set_xlabel('Lower age', fontsize=12, fontweight='bold')
        if idx == 0:
            ax.set_ylabel('Upper age', fontsize=12, fontweight='bold')
        ax.set_title(title, fontsize=14, fontweight='bold', pad=10)

        # Add text annotations with sigfigs
        for i in range(n_upper):
            for j in range(n_lower):
                if not np.isnan(data[i, j]):
                    text_label = sc.sigfig(data[i, j], 3, SI=True)
                    text = ax.text(j, i, text_label,
                                  ha="center", va="center", color="black", fontsize=8,
                                  fontweight='bold')

        # Add grid
        ax.set_xticks(np.arange(n_lower) - 0.5, minor=True)
        ax.set_yticks(np.arange(n_upper) - 0.5, minor=True)
        ax.grid(which="minor", color="white", linestyle='-', linewidth=2)

    # Add single colorbar for all three plots with sigfig formatting
    fig.subplots_adjust(right=0.9)
    cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])
    cbar = fig.colorbar(im, cax=cbar_ax)

    # Format colorbar ticks with sigfigs
    tick_locs = cbar.get_ticks()
    tick_labels = [sc.sigfig(val, 3, SI=True) for val in tick_locs]
    cbar.set_ticklabels(tick_labels)

    cbar.set_label('Cancers averted (2025-2100)', rotation=270, labelpad=20, fontsize=12, fontweight='bold')

    plt.suptitle('Impact of catch-up interventions by age range', fontsize=16, fontweight='bold', y=0.98)

    fig.tight_layout(rect=[0, 0, 0.9, 0.96])
    fig_name = 'figures/cohort_interventions_heatmap.png'
    sc.savefig(fig_name, dpi=150)

    print(f"\nHeatmap saved to: {fig_name}")
    return


# %% Run as a script
if __name__ == '__main__':

    location = 'kenya'
    plot_heatmaps(location)
    # make_single_bar(location)

    msim_dict = sc.loadobj(f'results/scens_{location}.obj')

