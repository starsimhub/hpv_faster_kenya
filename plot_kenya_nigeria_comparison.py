"""
Quick comparison plot of Kenya vs Nigeria cancer cases by age
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sciris as sc

# Set up plotting
sc.options(dpi=150)
sns.set_style('whitegrid')

# Load data
kenya_df = pd.read_csv('data/kenya_cancer_cases.csv')
nigeria_df = pd.read_csv('data/nigeria_cancer_cases.csv')

# Create figure
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# Panel A: Absolute numbers
ax1 = axes[0]
ax1.plot(kenya_df['age'], kenya_df['value'], marker='o', linewidth=2.5,
         markersize=8, color='#2ca02c', label='Kenya')
ax1.plot(nigeria_df['age'], nigeria_df['value'], marker='s', linewidth=2.5,
         markersize=8, color='#ff7f0e', label='Nigeria')

ax1.set_xlabel('Age (years)', fontsize=12, fontweight='bold')
ax1.set_ylabel('Number of cervical cancer cases (2020)', fontsize=12, fontweight='bold')
ax1.set_title('A. Absolute cancer cases by age', fontsize=13, fontweight='bold', pad=12)
ax1.legend(fontsize=11, loc='upper right')
ax1.grid(alpha=0.3)
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)

# Panel B: Normalized (as proportion of total)
ax2 = axes[1]
kenya_total = kenya_df['value'].sum()
nigeria_total = nigeria_df['value'].sum()

kenya_pct = 100 * kenya_df['value'] / kenya_total
nigeria_pct = 100 * nigeria_df['value'] / nigeria_total

ax2.plot(kenya_df['age'], kenya_pct, marker='o', linewidth=2.5,
         markersize=8, color='#2ca02c', label='Kenya')
ax2.plot(nigeria_df['age'], nigeria_pct, marker='s', linewidth=2.5,
         markersize=8, color='#ff7f0e', label='Nigeria')

ax2.set_xlabel('Age (years)', fontsize=12, fontweight='bold')
ax2.set_ylabel('% of total cervical cancer cases', fontsize=12, fontweight='bold')
ax2.set_title('B. Age distribution (normalized)', fontsize=13, fontweight='bold', pad=12)
ax2.legend(fontsize=11, loc='upper right')
ax2.grid(alpha=0.3)
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)

# Add summary statistics as text
kenya_peak_age = kenya_df.loc[kenya_df['value'].idxmax(), 'age']
nigeria_peak_age = nigeria_df.loc[nigeria_df['value'].idxmax(), 'age']

fig.text(0.5, 0.02,
         f'Kenya: {kenya_total:.0f} total cases, peak at age {kenya_peak_age:.0f} | ' +
         f'Nigeria: {nigeria_total:.0f} total cases, peak at age {nigeria_peak_age:.0f}',
         ha='center', fontsize=10, style='italic', color='#555555')

plt.suptitle('Cervical Cancer Cases: Kenya vs Nigeria (2020)',
             fontsize=14, fontweight='bold', y=0.98)
plt.tight_layout(rect=[0, 0.04, 1, 0.96])

# Save
plt.savefig('figures/kenya_nigeria_cancer_comparison.png', dpi=150, bbox_inches='tight')
print('Saved figure: figures/kenya_nigeria_cancer_comparison.png')
print(f'\nSummary:')
print(f'  Kenya: {kenya_total:.0f} total cases, peak at age {kenya_peak_age:.0f}')
print(f'  Nigeria: {nigeria_total:.0f} total cases, peak at age {nigeria_peak_age:.0f}')
print(f'  Nigeria/Kenya ratio: {nigeria_total/kenya_total:.1f}x')

plt.show()
