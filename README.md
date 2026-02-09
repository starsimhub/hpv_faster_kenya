# HPV-Faster Kenya Analysis

**Comparative effectiveness of catch-up HPV vaccination and test-and-treat strategies in Kenya**

## Overview

This repository contains agent-based modeling analysis of HPV-Faster strategies for cervical cancer prevention in Kenya. The analysis compares catch-up vaccination campaigns versus test-and-treat interventions across different age cohorts (10-60 years) to identify optimal strategies for women who missed routine adolescent HPV vaccination.

## Quick Start

### View Analysis Results

**📄 [Read the complete analysis narrative](ANALYSIS_NARRATIVE.md)** - Comprehensive report with methods, results, and policy implications (ready for PDF/HTML export)

### Key Findings

1. **Combined interventions** (vaccination + test-and-treat) are most effective
2. **Younger cohorts** (10-30 years) show highest efficiency for vaccination
3. **Vaccination prevents reinfection**, while test-and-treat alone leaves women vulnerable to new HPV infections
4. **70% coverage** may be optimal balance of impact and programmatic feasibility

## Repository Structure

```
├── ANALYSIS_NARRATIVE.md      # Analysis report
├── README.md                  # This file
│
├── run_sims.py                # Kenya simulation and calibration
├── run_scenarios.py           # Main scenario runner
├── analyzers.py               # Custom cohort cancer analyzers
├── utils.py                   # Helper functions
│
├── plot_tt_comparisons.py     # Visualization: test-and-treat vs vaccination
├── plot_bars.py               # Visualization: bar charts
├── plot_age_causal.py         # Visualization: age distribution plots
│
├── data/                      # Kenya demographic and cancer data
├── results/                   # Small outputs (CSVs, parameters - commit these)
├── raw_results/               # Large simulation files (don't commit)
├── figures/                   # Generated plots (don't commit)
└── z_archive/                 # Archived figures from analysis (don't commit) 
```

## Running the Analysis

### Prerequisites

```bash
pip install hpvsim sciris numpy pandas matplotlib seaborn scipy
```

### Execute Scenarios

```python
# In run_scenarios.py, configure:
location = 'kenya'         # or 'nigeria'
catchup_cov = 0.9          # 90% coverage (or 0.7 for 70%)
debug = 0                  # 0 for full run, 1 for quick test

# Then run:
python run_scenarios.py
```

**Note**: Full analysis requires HPC resources. Expect 5-15 minutes on HPC with multiple cores.

### Generate Visualizations

```python
import plot_tt_comparisons as ptc
ptc.plot_tt_comparison(location='kenya', coverage=90)

import plot_age_causal as pac
pac.plot_age_causal(location='kenya')
```

## Key Results Files

### Committed to Repository (in `results/`)
- **ANALYSIS_NARRATIVE.md** - Report 
- `results/kenya_pars.obj` - Calibrated model parameters
- `results/catchup_vax_cohort_cancers_90.csv` - Cancer outcomes by cohort (90% coverage)
- `results/catchup_vax_cohort_cancers_70.csv` - Cancer outcomes by cohort (70% coverage)
- `results/age_causal_infection_*.obj` - Age at causal infection data
- `results/cohort_status_2025.csv` - Cohort status snapshot
- `results/combined_impact_data_*.csv` - Combined impact analysis

### Large Files (in `raw_results/`, not committed)
- `raw_results/scens_kenya_90.obj` - 90% coverage scenario results (processed)
- `raw_results/scens_kenya_70.obj` - 70% coverage scenario results (processed)
- `raw_results/scens_kenya_*.msim` - Raw multi-simulation objects
- `raw_results/*.sim` - Individual simulation objects
- `raw_results/*_calib*.obj` - Calibration objects

### Figures (not committed)
- `figures/*.png` - Generated plots
- `z_archive/*.png` - Archived figures from analysis runs

## Exporting the Analysis

### To PDF

```bash
# Using pandoc
pandoc ANALYSIS_NARRATIVE.md -o HPV_Faster_Kenya_Analysis.pdf \
  --toc --number-sections \
  -V geometry:margin=1in

# Or use any Markdown-to-PDF tool
```

### To HTML

```bash
# Using pandoc
pandoc ANALYSIS_NARRATIVE.md -o HPV_Faster_Kenya_Analysis.html \
  --self-contained --toc \
  -c https://cdnjs.cloudflare.com/ajax/libs/github-markdown-css/5.1.0/github-markdown.min.css

# Or open in VSCode and use "Markdown Preview Enhanced" extension
```

## Data Sources

- Kenya age pyramid (population structure)
- Kenya age-standardized cancer incidence rates
- Kenya DHS 2014 (sexual behavior parameters)
- Calibrated HPVsim natural history parameters

## Contact

For questions or collaboration: info@hpvsim.org

**Last updated**: February 2026

---

*Analysis conducted using HPVsim agent-based model, calibrated to Kenya-specific demographic and epidemiological data.*
