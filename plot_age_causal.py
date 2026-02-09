"""
Plot age at causal infection distributions
"""

import numpy as np
import sciris as sc
import pandas as pd
import matplotlib.pyplot as pl
import seaborn as sns
import utils as ut


def get_age_causal_df(sim=None):
    """
    Make age causal dataframe from simulation analyzer results

    Args:
        sim: HPVsim simulation object with age_causal_infection analyzer

    Returns:
        age_causal_df: DataFrame with age and health event columns
    """
    dt_res = sim.get_analyzer('age_causal_infection')
    dt_dfs = sc.autolist()

    dt_df = pd.DataFrame()
    dt_df['Age'] = np.array(dt_res.age_causal)[np.array(dt_res.age_causal)<50]
    dt_df['Health event'] = 'Causal\ninfection'
    dt_dfs += dt_df

    dt_df = pd.DataFrame()
    dt_df['Age'] = np.array(dt_res.age_cin)[np.array(dt_res.age_causal)<65]
    dt_df['Health event'] = 'HSIL'
    dt_dfs += dt_df

    dt_df = pd.DataFrame()
    dt_df['Age'] = np.array(dt_res.age_cancer)[np.array(dt_res.age_causal)<90]
    dt_df['Health event'] = 'Cancer'
    dt_dfs += dt_df

    age_causal_df = pd.concat(dt_dfs)

    return age_causal_df


def plot_age_causal(location='kenya', from_file=True, ac_df=None):
    """
    Create boxplot showing age distribution of key health events
    (causal infection, HSIL, and cancer)

    Args:
        location: Country name (e.g., 'kenya')
        from_file: If True, load data from file; if False, use provided ac_df
        ac_df: Pre-loaded age causal DataFrame (used if from_file=False)

    Returns:
        fig: matplotlib figure object
    """
    if from_file:
        ac_df = sc.loadobj(f'results/age_causal_infection_{location}.obj')

    ac_colors = sc.gridcolors(3)

    ut.set_font(20)
    fig, ax = pl.subplots(1, 1, figsize=(6, 5))
    sns.boxplot(
        x="Health event", y="Age", data=ac_df, ax=ax,
        showfliers=False, palette=ac_colors, hue='Health event',
        hue_order=['Causal\ninfection', 'HSIL', 'Cancer']
    )
    ax.set_title(f'Age distribution\nof key health events')
    ax.set_xlabel('')
    ax.set_ylim([0, 100])

    fig.tight_layout()
    fig_name = f'figures/age_causal_{location}.png'
    sc.savefig(fig_name, dpi=100)

    print(f'Figure saved to: {fig_name}')

    return fig


def plot_age_causal_violin(location='kenya', from_file=True, ac_df=None):
    """
    Create violin plot showing distribution of age at causal infection
    with vaccination coverage thresholds

    Args:
        location: Country name (e.g., 'kenya')
        from_file: If True, load data from file; if False, use provided ac_df
        ac_df: Pre-loaded age causal DataFrame (used if from_file=False)

    Returns:
        fig: matplotlib figure object
    """
    if from_file:
        ac_df = sc.loadobj(f'results/age_causal_infection_{location}.obj')

    # Filter for just causal infection
    causal_df = ac_df[ac_df['Health event'] == 'Causal\ninfection'].copy()

    # Calculate percentiles for vaccination ages
    vax_ages = [15, 20, 25, 30, 35, 40]
    percentages = []
    for age in vax_ages:
        pct = (causal_df['Age'] <= age).sum() / len(causal_df) * 100
        percentages.append(pct)
        print(f'Age {age}: {pct:.1f}% of causal infections')

    # Create the violin plot
    ut.set_font(20)
    fig, ax = pl.subplots(1, 1, figsize=(6, 8))

    # Create violin plot
    parts = ax.violinplot([causal_df['Age'].values], positions=[0], widths=0.6,
                          showmeans=False, showmedians=True, vert=True)

    # Customize violin plot colors
    for pc in parts['bodies']:
        pc.set_facecolor('#3498db')
        pc.set_alpha(0.7)
        pc.set_edgecolor('black')
        pc.set_linewidth(1.5)

    # Style the median line
    parts['cmedians'].set_edgecolor('red')
    parts['cmedians'].set_linewidth(2)

    # Add horizontal lines for vaccination ages
    colors = ['#27ae60', '#1abc9c', '#f39c12', '#e74c3c', '#9b59b6', '#34495e']
    for i, (age, pct, color) in enumerate(zip(vax_ages, percentages, colors)):
        # Draw horizontal line
        ax.axhline(y=age, color=color, linestyle='--', linewidth=2, alpha=0.6, zorder=5)

        # Add text label on the right side - positioned outside the violin
        ax.text(0.45, age, f'{pct:.0f}%',
               fontsize=13, fontweight='bold', color=color,
               va='center', ha='left', zorder=15,
               bbox=dict(boxstyle='round,pad=0.4', facecolor='white',
                        edgecolor=color, linewidth=2.5, alpha=1.0))

    ax.set_ylabel('Age at causal infection (years)', fontsize=12, fontweight='bold')
    ax.set_title('Distribution of age at causal infection\nwith vaccination coverage thresholds',
                fontsize=14, fontweight='bold')
    ax.set_xlim(-0.5, 0.6)
    ax.set_ylim([10, 55])
    ax.set_xticks([])
    ax.grid(alpha=0.3, linestyle='--', axis='y')

    # Add text box explaining the percentages
    textstr = 'Lines show % of causal\ninfections occurring by\neach age (vaccination target)'
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    ax.text(0.98, 0.02, textstr, transform=ax.transAxes, fontsize=10,
           verticalalignment='bottom', horizontalalignment='right', bbox=props)

    fig.tight_layout()
    fig_name = f'figures/age_causal_violin_{location}.png'
    sc.savefig(fig_name, dpi=300)
    print(f'\nFigure saved to: {fig_name}')

    return fig


# %% Run as a script
if __name__ == '__main__':

    plot_age_causal(location='kenya')
    for location in ['kenya', 'nigeria']:
        plot_age_causal_violin(location=location)
