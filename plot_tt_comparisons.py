"""
Plot test & treat comparisons and heat maps
"""


import pylab as pl
import sciris as sc
import numpy as np
import utils as ut
import pandas as pd
import matplotlib.pyplot as plt


# ========== GLOBAL FIGURE SIZE SETTINGS ==========
# Standard width per panel for consistency
FIG_WIDTH_PER_PANEL = 8     # inches per panel

# Figure dimensions (width, height) - standardized for consistent appearance
FIGSIZE_2PANEL = (16, 6)       # 2 panels side-by-side (standard)
FIGSIZE_2PANEL_TALL = (16, 7)  # 2 panels side-by-side (taller for more content)
FIGSIZE_3PANEL = (24, 7)       # 3 panels side-by-side
FIGSIZE_6PANEL = (24, 12)      # 6 panels in 2x3 grid

# ========== GLOBAL FONT SIZE SETTINGS ==========
# Base font sizes defined for reference figure height of 6 inches
# Fonts will be scaled proportionally for taller figures to maintain visual consistency
FONT_SCALE_REF_HEIGHT = 6   # Reference height for base font sizes

FONT_BASE = 18              # Base font size for matplotlib
FONT_SUPTITLE = 16          # Main figure title (reduced)
FONT_PANEL_TITLE = 16       # Individual panel titles (A, B, C, etc.)
FONT_AXIS_LABEL = 10        # X and Y axis labels (reduced for better spacing)
FONT_LEGEND = 10            # Legend text (same as axis labels)
FONT_BAR_LABEL_LARGE = 11   # Value labels on bars (standard panels)
FONT_BAR_LABEL_SMALL = 9    # Value labels on bars (dense panels with multiple groups)


def scale_fonts(figure_height):
    """
    Scale font sizes proportionally based on figure height to maintain visual consistency.
    Returns a dict of scaled font sizes.

    Args:
        figure_height: Height of the figure in inches

    Returns:
        dict with keys: 'base', 'suptitle', 'panel_title', 'axis_label', 'legend', 'bar_large', 'bar_small'
    """
    scale_factor = figure_height / FONT_SCALE_REF_HEIGHT
    return {
        'base': int(FONT_BASE * scale_factor),
        'suptitle': int(FONT_SUPTITLE * scale_factor),
        'panel_title': int(FONT_PANEL_TITLE * scale_factor),
        'axis_label': int(FONT_AXIS_LABEL * scale_factor),
        'legend': int(FONT_LEGEND * scale_factor),
        'bar_large': int(FONT_BAR_LABEL_LARGE * scale_factor),
        'bar_small': int(FONT_BAR_LABEL_SMALL * scale_factor),
    }


def plot_tt_comparison(location, coverage=90):
    """
    Create a 3-panel figure comparing intervention strategies using IARC comparison results.
    Plots medians with uncertainty intervals (fill_between).

    Panel A: ASR cancer incidence time series
    Panel B: Annual incident cancers
    Panel C: Bar chart of total cancers averted (relative to routine vax baseline)
    """

    # Load IARC comparison results (with medians/intervals from multi-seed runs)
    msim_dict = sc.loadobj(f'raw_results/iarc_comparison_{location}_{coverage}.obj')

    # Define scenarios to compare (label, scenario key, color)
    scenarios = [
        ('No interventions',          'No interventions',          '#3ac9b1'),
        ('Routine vax only',          'Routine vax only',          '#fc8d62'),
        ('Catch-up to 30',            'Catch-up to 30',            '#8da0cb'),
        ('Catch-up to 40',            'Catch-up to 40',            '#e78ac3'),
        ('Catch-up to 40 + S&T to 65','Catch-up to 40 + S&T to 65','#a6d854'),
        ('S&T 10-65 only',            'S&T 10-65 only',           '#e7c700'),
    ]

    # Define cohorts for 10-60 age range
    cohorts = ['10_15', '15_20', '20_25', '25_30', '30_35', '35_40',
               '40_45', '45_50', '50_55', '55_60']

    # Extract data for each scenario
    scenario_data = {}
    scenario_totals = {}

    for label, scen_key, color in scenarios:
        if scen_key in msim_dict:
            mres = msim_dict[scen_key]
            years = mres['year']
            idx_start = sc.findinds(years, 2020)[0]

            data = {'years': years[idx_start:], 'color': color}

            # ASR cancer incidence (has .values, .low, .high from MultiSim.reduce)
            asr = mres.get('asr_cancer_incidence', None)
            if asr is not None:
                data['asr'] = asr.values[idx_start:]
                data['asr_low'] = asr.low[idx_start:]
                data['asr_high'] = asr.high[idx_start:]
            else:
                data['asr'] = data['asr_low'] = data['asr_high'] = None

            # Cancers
            cancers = mres['cancers']
            data['cancers'] = cancers.values[idx_start:]
            data['cancers_low'] = cancers.low[idx_start:]
            data['cancers_high'] = cancers.high[idx_start:]

            scenario_data[label] = data

            total_cohort_cancers = sum(
                mres.get(f'cohort_cancers_{c}', 0) for c in cohorts
            )
            scenario_totals[label] = total_cohort_cancers
            print(f'{label}: {total_cohort_cancers:,.0f} cancers in 10-60 cohort')

    # Font sizes
    fs_title = 20
    fs_label = 16
    fs_tick = 14
    fs_bar = 14

    # Helper functions
    def plot_with_band(ax, data, key):
        yrs = data['years']
        ax.fill_between(yrs, data[f'{key}_low'], data[f'{key}_high'],
                        color=data['color'], alpha=0.15)
        ax.plot(yrs, data[key], color=data['color'], linewidth=2.5, alpha=0.9)

    def style_ax(ax):
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.tick_params(axis='both', labelsize=fs_tick)
        ax.grid(alpha=0.2, linestyle='-', axis='y')
        ax.set_ylim(bottom=0)

    # ===== Figure 1: 2-panel time series =====
    ut.set_font(fs_tick)
    fig1, axes1 = plt.subplots(1, 2, figsize=(16, 5.25))

    # --- Panel A: ASR cancer incidence ---
    ax_a = axes1[0]
    for label, scen_key, color in scenarios:
        if label in scenario_data and scenario_data[label]['asr'] is not None:
            plot_with_band(ax_a, scenario_data[label], 'asr')
    ax_a.set_ylabel('ASR cancer incidence\n(per 100,000)', fontsize=fs_label)
    ax_a.set_title('A) ASR cancer incidence', fontsize=fs_title, fontweight='bold', loc='left')
    style_ax(ax_a)
    sc.SIticks(ax_a)

    # --- Panel B: Annual incident cancers ---
    ax_b = axes1[1]
    for label, scen_key, color in scenarios:
        if label in scenario_data:
            plot_with_band(ax_b, scenario_data[label], 'cancers')
    ax_b.set_ylabel('Annual incident cancers', fontsize=fs_label)
    ax_b.set_title('B) Annual incident cancers', fontsize=fs_title, fontweight='bold', loc='left')
    style_ax(ax_b)
    sc.SIticks(ax_b)
    ax_b.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f'{int(x):,}'))

    fig1.tight_layout()
    fig1_name = f'figures/tt_vs_vax_comparison_{location}_{coverage}.png'
    sc.savefig(fig1_name, dpi=300)
    print(f"Time series plot saved to: {fig1_name}")

    # ===== Figure 2: Model comparison bar chart =====
    baseline_total = scenario_totals['Routine vax only']
    bar_scenarios = [s for s in scenarios if s[0] not in ('No interventions', 'Routine vax only', 'S&T 10-65 only')]

    # IARC estimates (cancers averted, in thousands)
    iarc_estimates = {
        'Catch-up to 30': 160e3,
        'Catch-up to 40': 167e3,
        'Catch-up to 40 + S&T to 65': 334e3,
    }

    labels_bar = []
    hpvsim_vals = []
    iarc_vals = []
    bar_colors = []
    for label, scen_key, color in bar_scenarios:
        averted = baseline_total - scenario_totals[label]
        labels_bar.append(label)
        hpvsim_vals.append(averted)
        iarc_vals.append(iarc_estimates[label])
        bar_colors.append(color)

    midpoints = [(h + i) / 2 for h, i in zip(hpvsim_vals, iarc_vals)]
    err_low = [mid - min(h, i) for mid, h, i in zip(midpoints, hpvsim_vals, iarc_vals)]
    err_high = [max(h, i) - mid for mid, h, i in zip(midpoints, hpvsim_vals, iarc_vals)]

    fig2, ax_c = plt.subplots(1, 1, figsize=(7, 7))

    x_pos = np.arange(len(labels_bar))
    ax_c.bar(x_pos, midpoints, color=bar_colors, alpha=0.85, width=0.6)
    ax_c.errorbar(x_pos, midpoints, yerr=[err_low, err_high],
                  fmt='none', ecolor='#333333', elinewidth=2.5, capsize=8, capthick=2.5, zorder=5)

    # Label HPVsim and IARC values at the error bar ends
    for i, (x, h, ia) in enumerate(zip(x_pos, hpvsim_vals, iarc_vals)):
        lo, hi = min(h, ia), max(h, ia)
        ax_c.text(x + 0.22, lo, f'HPVsim\n{sc.sigfig(h, 3, SI=True)}',
                 fontsize=11, va='center', ha='left', color='#333333')
        ax_c.text(x + 0.22, hi, f'IARC\n{sc.sigfig(ia, 3, SI=True)}',
                 fontsize=11, va='center', ha='left', color='#333333')

    x_labels_formatted = [l.replace(' + ', '\n+ ').replace(' to ', '\nto ') for l in labels_bar]
    ax_c.set_xticks(x_pos)
    ax_c.set_xticklabels(x_labels_formatted, ha='center', fontsize=fs_tick)
    ax_c.set_ylabel('Cancers averted\n(vs routine vax only)', fontsize=fs_label)
    ax_c.set_title('Cancers averted: model synthesis', fontsize=fs_title, fontweight='bold', loc='left')
    style_ax(ax_c)
    sc.SIticks(ax_c)
    ax_c.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f'{int(x):,}'))

    fig2.tight_layout()
    fig2_name = f'figures/tt_vs_vax_bars_{location}_{coverage}.png'
    sc.savefig(fig2_name, dpi=300)
    print(f"Bar chart saved to: {fig2_name}")

    # Print summary
    print(f"\n=== Summary (cancers in 10-60 cohort, {location.capitalize()}) ===")
    print(f"Routine vax (baseline): {baseline_total:,.0f} cancers")
    for label, h, ia in zip(labels_bar, hpvsim_vals, iarc_vals):
        print(f"{label}: HPVsim={h:,.0f}, IARC={ia:,.0f}, midpoint={int((h+ia)/2):,}")

    return fig1, fig2


def plot_efficiency_frontier(location, coverage=90):
    """
    Create incremental analysis to identify optimal age targeting.

    Panel A: Incremental cancers averted when extending age range
    Panel B: Incremental NNV (doses per additional cancer averted)

    Helps answer: Which age group should be targeted given resource constraints?
    """

    # Load data
    msim_dict = sc.loadobj(f'raw_results/scens_{location}_{coverage}.obj')

    # Define age ranges
    lower_ages = [10, 15, 20, 25, 30, 35, 40, 45, 50, 55]
    age_bands = [5, 10, 15, 20, 25, 30, 35, 40, 45, 50]
    all_cohorts = ['10_15', '15_20', '20_25', '25_30', '30_35', '35_40',
                   '40_45', '45_50', '50_55', '55_60']

    # Store results
    results = []

    for lower in lower_ages:
        for band in age_bands:
            upper = lower + band
            if upper > 60:
                continue

            scenarios = {
                'V': f'Catch-up {lower}-{upper}: V',
                'TTV': f'Catch-up {lower}-{upper}: TTV'
            }

            for intv_type, scen_key in scenarios.items():
                if scen_key in msim_dict:
                    cohorts_in_range = [c for c in all_cohorts
                                       if int(c.split('_')[0]) >= lower and int(c.split('_')[1]) <= upper]

                    total_cancers = 0
                    baseline_cancers = 0
                    for cohort in cohorts_in_range:
                        key = f'cohort_cancers_{cohort}'
                        if key in msim_dict[scen_key]:
                            total_cancers += msim_dict[scen_key][key]
                        if key in msim_dict['Baseline']:
                            baseline_cancers += msim_dict['Baseline'][key]

                    cancers_averted = baseline_cancers - total_cancers

                    doses = 0
                    if 'vaccinations' in msim_dict[scen_key]:
                        doses = msim_dict[scen_key]['vaccinations'].sum()

                    results.append({
                        'lower': lower,
                        'upper': upper,
                        'intervention': intv_type,
                        'cancers_averted': cancers_averted,
                        'doses': doses
                    })

    df = pd.DataFrame(results)

    # Create figure
    figsize = FIGSIZE_2PANEL_TALL
    fonts = scale_fonts(figsize[1])
    ut.set_font(fonts['base'])
    fig, axes = plt.subplots(1, 2, figsize=figsize)

    colors = {'V': '#3498db', 'TTV': '#27ae60'}

    # Panel A: Incremental cancers averted
    ax_a = axes[0]

    for intv_type in ['V', 'TTV']:
        df_10 = df[(df['lower'] == 10) & (df['intervention'] == intv_type)].sort_values('upper')
        
        if len(df_10) > 0:
            uppers = df_10['upper'].values
            cancers = df_10['cancers_averted'].values
            incremental_cancers = np.diff(cancers, prepend=0)
            
            x_pos = np.arange(len(uppers))
            width = 0.35
            offset = -width/2 if intv_type == 'V' else width/2
            
            bars = ax_a.bar(x_pos + offset, incremental_cancers, width,
                          label=intv_type, color=colors[intv_type], alpha=0.7)
            
            for bar, val in zip(bars, incremental_cancers):
                if val > 0:
                    height = bar.get_height()
                    ax_a.text(bar.get_x() + bar.get_width() / 2., height,
                            f'{sc.sigfig(val, 2, SI=True)}',
                            ha='center', va='bottom', fontsize=fonts['bar_small'], fontweight='bold')

    ax_a.set_ylabel('Incremental cancers averted', fontsize=fonts['axis_label'], fontweight='bold')
    ax_a.set_title('A) Marginal health impact:\nIncremental cancers averted', fontsize=fonts['panel_title'], fontweight='bold', loc='left')
    ax_a.set_xticks(x_pos)
    ax_a.set_xticklabels([int(u) for u in uppers], rotation=0)
    ax_a.set_ylim(bottom=0)
    ax_a.grid(alpha=0.3, linestyle='--', axis='y')
    ax_a.tick_params(axis='both', labelsize=fonts['axis_label'])
    ax_a.legend(fontsize=fonts['legend'], loc='upper left', frameon=False)
    sc.SIticks(ax_a)
    ax_a.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f'{int(x):,}'))

    # Panel B: Incremental NNV
    ax_b = axes[1]

    df_10_v = df[(df['lower'] == 10) & (df['intervention'] == 'V')].sort_values('upper')
    if len(df_10_v) > 0:
        uppers = df_10_v['upper'].values
        cancers = df_10_v['cancers_averted'].values
        doses = df_10_v['doses'].values
        incremental_cancers = np.diff(cancers, prepend=0)
        incremental_doses = np.diff(doses, prepend=0)
        incremental_nnv = incremental_doses / incremental_cancers
        incremental_nnv[incremental_cancers == 0] = np.nan
        
        x_pos = np.arange(len(uppers))
        
        bars = ax_b.bar(x_pos, incremental_nnv, width=0.7,
                       color='#34495e', alpha=0.7, edgecolor='black', linewidth=1)
        
        for bar, val in zip(bars, incremental_nnv):
            if not np.isnan(val):
                height = bar.get_height()
                ax_b.text(bar.get_x() + bar.get_width() / 2., height,
                        f'{int(val)}',
                        ha='center', va='bottom', fontsize=fonts['bar_small'], fontweight='bold')

    ax_b.set_ylabel('Incremental NNV\n(doses per additional cancer averted)', fontsize=fonts['axis_label'], fontweight='bold')
    ax_b.set_title('B) Marginal efficiency:\nIncremental NNV (V only)', fontsize=fonts['panel_title'], fontweight='bold', loc='left')
    ax_b.set_xticks(x_pos)
    ax_b.set_xticklabels([int(u) for u in uppers], rotation=0)
    ax_b.set_ylim(bottom=0)
    ax_b.grid(alpha=0.3, linestyle='--', axis='y')
    ax_b.tick_params(axis='both', labelsize=fonts['axis_label'])
    sc.SIticks(ax_b)
    ax_b.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f'{int(x):,}'))

    plt.suptitle(f'Incremental analysis for catch-up vaccination - {location.capitalize()}',
                 fontsize=fonts['suptitle'], fontweight='bold')

    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig_name = f'figures/efficiency_frontier_{location}_{coverage}.png'
    sc.savefig(fig_name, dpi=150)

    print(f"\nIncremental analysis plot saved to: {fig_name}")

    # Print key findings
    print(f"\n=== Incremental Analysis Summary ({location.capitalize()}) ===")
    print("\nIncremental cancers averted and NNV for extending age range from 10:")
    df_10_v = df[(df['lower'] == 10) & (df['intervention'] == 'V')].sort_values('upper')
    if len(df_10_v) > 0:
        uppers = df_10_v['upper'].values
        cancers = df_10_v['cancers_averted'].values
        doses = df_10_v['doses'].values
        for i in range(1, len(uppers)):
            inc_cancers = cancers[i] - cancers[i-1]
            inc_doses = doses[i] - doses[i-1]
            inc_nnv = inc_doses / inc_cancers if inc_cancers > 0 else np.nan
            print(f"  10-{int(uppers[i-1])} → 10-{int(uppers[i])}: {inc_cancers:,.0f} cancers, NNV = {inc_nnv:.0f}")

    return fig, axes, df


def plot_combined_impact(location, coverage=90):
    """
    Create comprehensive 6-panel combined impact analysis figure.

    Column 1 (Total impact):
    - Panel A: Total cancers averted extending upwards from age 10 (V vs TTV)
    - Panel D: Total cancers averted extending downwards to age 60 (TT vs TTV) - reversed

    Column 2 (Marginal impact):
    - Panel B: Incremental cancers averted by 5-year age band (V vs TTV)
    - Panel E: Incremental cancers averted by 5-year age band (TT vs TTV) - reversed

    Column 3 (Efficiency):
    - Panel C: Incremental NNV for extending from age 10 (V only)
    - Panel F: Comparison of all interventions (V, TT, TTV) by age band
    """

    # Load data
    msim_dict = sc.loadobj(f'raw_results/scens_{location}_{coverage}.obj')

    # Define age ranges
    lower_ages = [10, 15, 20, 25, 30, 35, 40, 45, 50, 55]
    age_bands = [5, 10, 15, 20, 25, 30, 35, 40, 45, 50]
    all_cohorts = ['10_15', '15_20', '20_25', '25_30', '30_35', '35_40',
                   '40_45', '45_50', '50_55', '55_60']

    # Store results for panels A & B
    results = []

    for lower in lower_ages:
        for band in age_bands:
            upper = lower + band
            if upper > 60:
                continue

            scenarios = {
                'V': f'Catch-up {lower}-{upper}: V',
                'TT': f'Catch-up {lower}-{upper}: TT',
                'TTV': f'Catch-up {lower}-{upper}: TTV'
            }

            for intv_type, scen_key in scenarios.items():
                if scen_key in msim_dict:
                    cohorts_in_range = [c for c in all_cohorts
                                       if int(c.split('_')[0]) >= lower and int(c.split('_')[1]) <= upper]

                    total_cancers = 0
                    baseline_cancers = 0
                    for cohort in cohorts_in_range:
                        key = f'cohort_cancers_{cohort}'
                        if key in msim_dict[scen_key]:
                            total_cancers += msim_dict[scen_key][key]
                        if key in msim_dict['Baseline']:
                            baseline_cancers += msim_dict['Baseline'][key]

                    cancers_averted = baseline_cancers - total_cancers

                    doses = 0
                    if 'vaccinations' in msim_dict[scen_key]:
                        doses = msim_dict[scen_key]['vaccinations'].sum()

                    results.append({
                        'lower': lower,
                        'upper': upper,
                        'intervention': intv_type,
                        'cancers_averted': cancers_averted,
                        'doses': doses
                    })

    df = pd.DataFrame(results)

    # Create 2x3 figure (6 panels)
    figsize = FIGSIZE_6PANEL
    fonts = scale_fonts(figsize[1])
    ut.set_font(fonts['base'])
    fig = plt.figure(figsize=figsize)
    gs = fig.add_gridspec(2, 3, hspace=0.45, wspace=0.3)

    colors = {'V': '#3498db', 'TT': '#e74c3c', 'TTV': '#27ae60'}

    # ========== Panel A: Total cancers averted (extending from 10 upwards) ==========
    ax_a = fig.add_subplot(gs[0, 0])

    for intv_type in ['V', 'TTV']:
        df_10 = df[(df['lower'] == 10) & (df['intervention'] == intv_type)].sort_values('upper')

        if len(df_10) > 0:
            uppers = df_10['upper'].values
            cancers = df_10['cancers_averted'].values  # Total/cumulative cancers averted

            x_pos = np.arange(len(uppers))
            width = 0.35
            offset = -width/2 if intv_type == 'V' else width/2

            bars = ax_a.bar(x_pos + offset, cancers, width,
                          label=intv_type, color=colors[intv_type], alpha=0.7)

    # Create age band labels for x-axis
    age_band_labels_total = ['10-15', '10-20', '10-25', '10-30', '10-35', '10-40', '10-45', '10-50', '10-55', '10-60']

    ax_a.set_ylabel('Total cancers averted', fontsize=fonts['axis_label'], fontweight='bold')
    ax_a.set_title('A) Total impact:\nExtending upwards from age 10', fontsize=fonts['panel_title'], fontweight='bold', loc='left')
    ax_a.set_xticks(x_pos)
    ax_a.set_xticklabels(age_band_labels_total, rotation=45, ha='right')
    ax_a.set_ylim(bottom=0)
    ax_a.grid(alpha=0.3, linestyle='--', axis='y')
    ax_a.tick_params(axis='both', labelsize=fonts['axis_label'])
    ax_a.legend(fontsize=fonts['legend'], loc='upper left', frameon=False)
    sc.SIticks(ax_a)

    # ========== Panel B: Marginal cancers averted (extending from 10 upwards) ==========
    ax_b = fig.add_subplot(gs[0, 1])

    for intv_type in ['V', 'TTV']:
        df_10 = df[(df['lower'] == 10) & (df['intervention'] == intv_type)].sort_values('upper')

        if len(df_10) > 0:
            uppers = df_10['upper'].values
            cancers = df_10['cancers_averted'].values
            incremental_cancers = np.diff(cancers, prepend=0)

            x_pos = np.arange(len(uppers))
            width = 0.35
            offset = -width/2 if intv_type == 'V' else width/2

            bars = ax_b.bar(x_pos + offset, incremental_cancers, width,
                          label=intv_type, color=colors[intv_type], alpha=0.7)

    # Create incremental age band labels (10-15, 15-20, etc.)
    age_band_labels_incr = ['10-15', '15-20', '20-25', '25-30', '30-35', '35-40', '40-45', '45-50', '50-55', '55-60']

    ax_b.set_ylabel('Incremental cancers averted', fontsize=fonts['axis_label'], fontweight='bold')
    ax_b.set_title('B) Marginal impact:\nBy 5-year age band', fontsize=fonts['panel_title'], fontweight='bold', loc='left')
    ax_b.set_xticks(x_pos)
    ax_b.set_xticklabels(age_band_labels_incr, rotation=45, ha='right')
    ax_b.set_ylim(bottom=0)
    ax_b.grid(alpha=0.3, linestyle='--', axis='y')
    ax_b.tick_params(axis='both', labelsize=fonts['axis_label'])
    ax_b.legend(fontsize=fonts['legend'], loc='upper left', frameon=False)
    sc.SIticks(ax_b)

    # ========== Panel C: Incremental NNV (extending from 10 upwards) ==========
    ax_c = fig.add_subplot(gs[0, 2])

    df_10_v = df[(df['lower'] == 10) & (df['intervention'] == 'V')].sort_values('upper')
    if len(df_10_v) > 0:
        uppers = df_10_v['upper'].values
        cancers = df_10_v['cancers_averted'].values
        doses = df_10_v['doses'].values
        incremental_cancers = np.diff(cancers, prepend=0)
        incremental_doses = np.diff(doses, prepend=0)
        incremental_nnv = incremental_doses / incremental_cancers
        incremental_nnv[incremental_cancers == 0] = np.nan

        x_pos = np.arange(len(uppers))

        bars = ax_c.bar(x_pos, incremental_nnv, width=0.7,
                       color='#34495e', alpha=0.7, edgecolor='black', linewidth=1)

        for bar, val in zip(bars, incremental_nnv):
            if not np.isnan(val):
                height = bar.get_height()
                ax_c.text(bar.get_x() + bar.get_width() / 2., height,
                        f'{int(val)}',
                        ha='center', va='bottom', fontsize=fonts['bar_large'], fontweight='bold')

    ax_c.set_ylabel('Incremental NNV\n(doses per additional cancer averted)', fontsize=fonts['axis_label'], fontweight='bold')
    ax_c.set_title('C) Marginal efficiency:\nIncremental NNV (V only)', fontsize=fonts['panel_title'], fontweight='bold', loc='left')
    ax_c.set_xticks(x_pos)
    ax_c.set_xticklabels(age_band_labels_incr, rotation=45, ha='right')
    ax_c.set_ylim(bottom=0, top=175)
    ax_c.grid(alpha=0.3, linestyle='--', axis='y')
    ax_c.tick_params(axis='both', labelsize=fonts['axis_label'])
    # sc.SIticks(ax_c)

    # ========== Panel D: Total cancers averted (fixing upper at 60, varying lower) - REVERSED ==========
    ax_d = fig.add_subplot(gs[1, 0])

    # Get data for scenarios that end at age 60, sorted DESCENDING by lower age (reverse order)
    df_60 = df[df['upper'] == 60].sort_values('lower', ascending=False)

    for intv_type in ['TT', 'TTV']:
        df_60_intv = df_60[df_60['intervention'] == intv_type]

        if len(df_60_intv) > 0:
            lowers = df_60_intv['lower'].values
            cancers = df_60_intv['cancers_averted'].values  # Total/cumulative cancers averted

            x_pos = np.arange(len(lowers))
            width = 0.35
            offset = -width/2 if intv_type == 'TT' else width/2

            bars = ax_d.bar(x_pos + offset, cancers, width,
                          label=intv_type, color=colors[intv_type], alpha=0.7)

    # Create age band labels for x-axis (reversed)
    age_band_labels_total_rev = ['55-60', '50-60', '45-60', '40-60', '35-60', '30-60', '25-60', '20-60', '15-60', '10-60']

    ax_d.set_ylabel('Total cancers averted', fontsize=fonts['axis_label'], fontweight='bold')
    ax_d.set_title('D) Total impact:\nExtending downwards to age 60', fontsize=fonts['panel_title'], fontweight='bold', loc='left')
    ax_d.set_xticks(np.arange(len(lowers)))
    ax_d.set_xticklabels(age_band_labels_total_rev, rotation=45, ha='right')
    ax_d.set_ylim(bottom=0)
    ax_d.grid(alpha=0.3, linestyle='--', axis='y')
    ax_d.tick_params(axis='both', labelsize=fonts['axis_label'])
    ax_d.legend(fontsize=fonts['legend'], loc='upper left', frameon=False)
    sc.SIticks(ax_d)

    # ========== Panel E: Marginal cancers averted (fixing upper at 60, varying lower) - REVERSED ==========
    ax_e = fig.add_subplot(gs[1, 1])

    # Get data for scenarios that end at age 60, sorted DESCENDING by lower age (reverse order)
    df_60 = df[df['upper'] == 60].sort_values('lower', ascending=False)

    for intv_type in ['TT', 'TTV']:
        df_60_intv = df_60[df_60['intervention'] == intv_type]

        if len(df_60_intv) > 0:
            lowers = df_60_intv['lower'].values
            cancers = df_60_intv['cancers_averted'].values

            # Calculate incremental values (working from narrow to wide bands)
            # For reversed order: each bar shows additional cancers from adding the next band
            incremental_cancers = np.diff(cancers, prepend=0)

            x_pos = np.arange(len(lowers))
            width = 0.35
            offset = -width/2 if intv_type == 'TT' else width/2

            bars = ax_e.bar(x_pos + offset, incremental_cancers, width,
                          label=intv_type, color=colors[intv_type], alpha=0.7)

    # Create incremental age band labels (reversed: 55-60, 50-55, etc.)
    age_band_labels_incr_rev = ['55-60', '50-55', '45-50', '40-45', '35-40', '30-35', '25-30', '20-25', '15-20', '10-15']

    ax_e.set_ylabel('Incremental cancers averted', fontsize=fonts['axis_label'], fontweight='bold')
    ax_e.set_title('E) Marginal impact:\nBy 5-year age band', fontsize=fonts['panel_title'], fontweight='bold', loc='left')
    ax_e.set_xticks(np.arange(len(lowers)))
    ax_e.set_xticklabels(age_band_labels_incr_rev, rotation=45, ha='right')
    ax_e.set_ylim(bottom=0)
    ax_e.grid(alpha=0.3, linestyle='--', axis='y')
    ax_e.tick_params(axis='both', labelsize=fonts['axis_label'])
    ax_e.legend(fontsize=fonts['legend'], loc='upper left', frameon=False)
    sc.SIticks(ax_e)

    # ========== Panel F: Cancers averted by intervention for each age band ==========
    ax_f = fig.add_subplot(gs[1, 2])

    # Get data for individual age bands (5-year bands)
    age_bands_5yr = [(10, 15), (15, 20), (20, 25), (25, 30), (30, 35),
                     (35, 40), (40, 45), (45, 50), (50, 55), (55, 60)]
    band_labels = [f'{l}-{u}' for l, u in age_bands_5yr]

    # Collect cancers averted for each band and intervention
    band_results = {intv: [] for intv in ['V', 'TT', 'TTV']}

    for lower, upper in age_bands_5yr:
        for intv_type in ['V', 'TT', 'TTV']:
            df_band = df[(df['lower'] == lower) & (df['upper'] == upper) & (df['intervention'] == intv_type)]
            if len(df_band) > 0:
                band_results[intv_type].append(df_band['cancers_averted'].values[0])
            else:
                band_results[intv_type].append(0)

    # Create grouped bar chart
    x_pos = np.arange(len(band_labels))
    width = 0.25

    for i, intv_type in enumerate(['V', 'TT', 'TTV']):
        offset = width * (i - 1)
        bars = ax_f.bar(x_pos + offset, band_results[intv_type], width,
                       label=intv_type, color=colors[intv_type], alpha=0.7)

    ax_f.set_ylabel('Cancers averted', fontsize=fonts['axis_label'], fontweight='bold')
    ax_f.set_title('F) Which intervention\nfor which age band?', fontsize=fonts['panel_title'], fontweight='bold', loc='left')
    ax_f.set_xticks(x_pos)
    ax_f.set_xticklabels(band_labels, rotation=45, ha='right')
    ax_f.set_ylim(bottom=0)
    ax_f.grid(alpha=0.3, linestyle='--', axis='y')
    ax_f.tick_params(axis='both', labelsize=fonts['axis_label'])
    ax_f.legend(fontsize=fonts['legend'], loc='upper left', frameon=False)
    sc.SIticks(ax_f)

    plt.suptitle(f'Total and marginal impact of catch-up interventions - {location.capitalize()}',
                 fontsize=fonts['suptitle'], fontweight='bold', y=0.995)

    fig_name = f'figures/combined_impact_6panel_{location}_{coverage}.png'
    fig.tight_layout()
    sc.savefig(fig_name, dpi=150)

    print(f"\n6-panel combined impact plot saved to: {fig_name}")

    # Export data to CSV and Excel
    csv_name = f'results/combined_impact_data_{location}_{coverage}.csv'
    excel_name = f'figures/combined_impact_data_{location}_{coverage}.xlsx'

    # Sort dataframe for better readability
    df_export = df.sort_values(['intervention', 'lower', 'upper']).reset_index(drop=True)

    # Add age range column for clarity
    df_export['age_range'] = df_export['lower'].astype(str) + '-' + df_export['upper'].astype(str)

    # Reorder columns
    df_export = df_export[['age_range', 'lower', 'upper', 'intervention', 'cancers_averted', 'doses']]

    # Export to CSV
    df_export.to_csv(csv_name, index=False)
    print(f"Data exported to CSV: {csv_name}")

    # Export to Excel
    df_export.to_excel(excel_name, index=False, sheet_name='Combined Impact Data')
    print(f"Data exported to Excel: {excel_name}")

    return fig, (ax_a, ax_b, ax_c, ax_d, ax_e, ax_f), df


def plot_heatmaps(location, coverage=90, add_tt=False):

    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt

    # Load catch-up vaccination scenario data
    msim_dict = sc.loadobj(f'raw_results/scens_{location}_{coverage}.obj')

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
    figsize = FIGSIZE_3PANEL
    fonts = scale_fonts(figsize[1])
    ut.set_font(fonts['base'])
    fig, axes = plt.subplots(1, 3, figsize=figsize)

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
        if idx == 0:
            ax.set_ylabel('Upper age', fontsize=fonts['axis_label'], fontweight='bold')
        ax.set_title(title, fontsize=fonts['panel_title'], fontweight='bold', pad=10)
        ax.tick_params(axis='both', labelsize=fonts['axis_label'])

        # Add text annotations with sigfigs
        for i in range(n_upper):
            for j in range(n_lower):
                if not np.isnan(data[i, j]):
                    text_label = sc.sigfig(data[i, j], 3, SI=True)
                    text = ax.text(j, i, text_label,
                                  ha="center", va="center", color="black", fontsize=fonts['bar_small'],
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

    cbar.set_label('Cancers averted (2025-2100)', rotation=270, labelpad=20, fontsize=fonts['axis_label'], fontweight='bold')

    plt.suptitle('Impact of catch-up interventions by age range', fontsize=fonts['suptitle'], fontweight='bold', y=0.98)

    fig.tight_layout(rect=[0, 0, 0.9, 0.96])
    fig_name = 'figures/cohort_interventions_heatmap.png'
    sc.savefig(fig_name, dpi=150)

    print(f"\nHeatmap saved to: {fig_name}")
    return


# %% Run as a script
if __name__ == '__main__':

    location = 'kenya'
    coverage = 70

    # Create test & treat comparison plot
    plot_tt_comparison(location, coverage=coverage)

    # Create 6-panel combined impact analysis
    plot_combined_impact(location, coverage=coverage)

    # Create heatmaps
    plot_heatmaps(location, coverage=coverage)

