# HPV-Faster Strategy Analysis for Kenya
## Comparative Effectiveness of Catch-up Vaccination and Test-and-Treat Interventions

---

## Executive Summary

This analysis evaluates the effectiveness of HPV-Faster strategies in Kenya, comparing catch-up vaccination campaigns versus test-and-treat interventions across different age cohorts. The primary research question is: **What is the optimal strategy for preventing cervical cancer in women who have aged out of routine adolescent HPV vaccination programs?**

### Key Findings

1. **Combined interventions (vaccination + test-and-treat) are most effective** at preventing future cervical cancers, particularly in younger age cohorts (10-30 years)

2. **Vaccination alone prevents reinfection**, while test-and-treat alone treats current disease but leaves women vulnerable to new HPV infections

3. **Efficiency varies by age cohort**: Younger cohorts (10-30 years) show the highest return on investment for vaccination interventions, as they have longer remaining lifetimes and lower existing disease burden

4. **Coverage matters**: Analysis was conducted at both 70% and 90% catch-up campaign coverage levels to assess programmatic feasibility

---

## Background

### HPV-Faster Initiative

HPV-Faster refers to accelerated catch-up vaccination strategies targeting individuals who missed routine adolescent HPV vaccination. In Kenya, routine vaccination began in 2019 for girls aged 9-10 years. Women born before the vaccine era (or who missed routine vaccination) remain at risk for HPV infection and subsequent cervical cancer.

### Study Context: Kenya and Nigeria

The primary focus is on Kenya, but Nigeria is also included since the estimated age of causal infection differs, which leads to significant differences in the estimated impct of catch-up vaccination.

---

## Methods

### Simulation Framework

The analysis uses **HPVsim**, an agent-based model that simulates:
- HPV transmission dynamics
- Natural history of infection to cancer
- Impact of vaccination and treatment interventions
- Population demographics (births, deaths, sexual behavior)

### Calibration

The simulation was calibrated to Kenya-specific data:
- Age pyramid (population structure)
- Cancer incidence rates (age-standardized)
- Sexual behavior data from Demographic and Health Surveys (DHS)
- Natural history parameters for HPV progression

**Data sources:**
- [kenya_age_pyramid.csv](data/kenya_age_pyramid.csv) - Population age structure
- [kenya_asr_cancer_incidence.csv](data/kenya_asr_cancer_incidence.csv) - Age-standardized cancer rates
- [kenya_cancer_cases.csv](data/kenya_cancer_cases.csv) - Observed cancer cases by year
- Sexual debut parameters derived from 2014 Kenya DHS

### Intervention Scenarios

#### Background Interventions (all scenarios)

1. **Routine screening & treatment** (starting 2020)
   - HPV DNA testing for women aged 30-50
   - 15% lifetime screening coverage
   - Screen-positive women triaged to ablation, excision, or radiation based on disease severity
   - 90% triage coverage, 75% treatment coverage

2. **Routine adolescent vaccination** (starting 2019)
   - Target: girls aged 9-10 years
   - Nonavalent vaccine (covering HPV types 16, 18, 31, 33, 45, 52, 58, 6, 11)
   - Coverage: scale-up from 25% (2019) to 90% by 2025

#### Catch-up Campaign Scenarios (2026)

**Age ranges tested:** 10-15, 10-20, 10-25, ..., up to 10-60 years (various band widths)

**Intervention types:**
1. **Vaccination only (V)**: Single-dose nonavalent HPV vaccine
2. **Test-and-treat only (TT)**:
   - HPV DNA testing
   - Treatment assigned based on triage (ablation, excision, or radiation)
   - 90% treatment coverage (higher than routine screening)
3. **Both (TTV)**: Combined vaccination + test-and-treat

**Coverage levels:** 70% and 90% campaign coverage

### Cohort Cancer Analysis

To assess age-specific impact, custom analyzers tracked cumulative lifetime cancers for specific age cohorts:
- **Cohorts tracked**: 10-15, 15-20, 20-25, 25-30, 30-35, 35-40, 40-45, 45-50, 50-55, 55-60 years (at start of 2026)
- **Follow-up period**: 2026 to 2100 (tracking cohorts as they age)
- **Primary outcome**: Cumulative cancers occurring in each cohort over their remaining lifetimes

### Simulation Parameters

- **Population size**: 10,000 agents
- **Time step**: 0.25 years (quarterly)
- **Simulation period**: 1960-2100
- **Uncertainty**: 10 random seeds per scenario (median, 10th, 90th percentiles reported)
- **HPV genotypes modeled**: 16, 18, high-risk-5 (31/33/45/52/58), other high-risk

---

## Results

### 1. Test-and-Treat vs. Vaccination Comparison

**Key insight**: Test-and-treat alone provides immediate benefit by treating existing disease, but women remain susceptible to new HPV infections. Vaccination prevents both initial and reinfection, providing longer-term protection.

**Evidence from time-series analysis:**
- In the 10-60 year cohort, test-and-treat alone shows initial cancer reduction but effect plateaus
- Vaccination alone shows steady, sustained cancer prevention over decades
- Combined approach (TTV) achieves the greatest reduction

*Visualization*: [tt_vs_vax_comparison_kenya_90.png](generated by running `plot_tt_comparison` within `plot_tt_comparisons.py`)

### 2. Age-Cohort Specific Impact

Different age cohorts show varying benefit from catch-up vaccination:

**Younger cohorts (10-30 years)**:
- Higher vaccine effectiveness (less existing disease, more opportunity for prevention)
- Longer remaining lifetimes for benefit to accrue
- Lower current disease prevalence (fewer infections to treat)

**Middle cohorts (30-45 years)**:
- Mixed benefit: some existing disease burden
- Test-and-treat becomes relatively more important
- Vaccination still prevents future infections

**Older cohorts (45-60 years)**:
- Higher existing disease burden
- Shorter remaining lifetimes
- Test-and-treat may be more immediately impactful
- Vaccination benefit lower due to acquired immunity from past infections

*Visualization*: [catchup_vax_cohorts_kenya_70.png](figures/catchup_vax_cohorts_kenya_70.png), [catchup_vax_cohorts_kenya_90.png](figures/catchup_vax_cohorts_kenya_90.png), generated by running `plot_bars.py`

### 3. Efficiency Analysis

**Metric**: Cancers averted per vaccine dose

**Findings**:
- **Most efficient age range**: 10-30 years for vaccination
- **Diminishing returns**: Expanding campaigns to older ages (>40 years) shows lower efficiency
- **Trade-off**: Narrower campaigns are more efficient but miss more at-risk individuals; broader campaigns are less efficient but have greater population impact

**Campaign width considerations**:
- Narrow campaigns (5-10 year bands): Higher efficiency, lower population coverage
- Broad campaigns (e.g., 10-60): Lower efficiency, greater absolute cancer prevention

*Visualizations*:
- [combined_impact_6panel_kenya_90.png](figures/combined_impact_6panel_kenya_90.png) - Efficiency and trade-off visualization, generated by running `plot_combined_impact` within `plot_tt_comparison.py`.

### 4. Coverage Sensitivity (70% vs 90%)

Comparing 70% and 90% catch-up campaign coverage:

**90% coverage**:
- Greater absolute impact (more cancers averted)
- May be challenging to achieve operationally
- Higher programmatic costs

**70% coverage**:
- More realistic for resource-limited settings
- Still substantial impact (captures most available benefit)
- May be optimal from feasibility perspective

*Data*: [catchup_vax_cohort_cancers_70.csv](results/catchup_vax_cohort_cancers_70.csv), [catchup_vax_cohort_cancers_90.csv](results/catchup_vax_cohort_cancers_90.csv)

### 5. Baseline Projections

Without catch-up interventions, projections show:
- Continued high burden of cervical cancer over next several decades
- Gradual decline as routine vaccination of adolescents scales up (post-2025)
- Significant residual burden from cohorts who missed routine vaccination

*Visualization*: [baseline_sim.png](z_archive/baseline_sim.png)

### 6. Cohort Status in 2025

Analysis of HPV exposure and disease status at the start of catch-up campaigns (2025) reveals:
- Younger cohorts: mostly HPV-naïve, ideal for vaccination
- Middle cohorts: mixed exposure (some cleared infections, some current infections)
- Older cohorts: high cumulative exposure, more advanced disease

This explains differential intervention effectiveness by age.

*Data*: [cohort_status_2025.csv](results/cohort_status_2025.csv)
*Visualizations*: [cohort_vaccination_eligibility.png](z_archive/cohort_vaccination_eligibility.png), [cohort_simplified_2025.png](z_archive/cohort_simplified_2025.png)

---

## Discussion

### Policy Implications

1. **Prioritize younger cohorts** for catch-up vaccination (ages 10-30)
2. **Consider combined interventions** (vaccination + test-and-treat) for maximum impact
3. **Tailor interventions by age**: vaccination-focused for younger women, test-and-treat for older women
4. **Coverage targets**: 70% coverage may be optimal balance of impact and feasibility

### Limitations

1. **Model assumptions**: Sexual behavior parameters is derived from DHS surveys and is subject to reporting bias; model calibration fits to GLOBOCAN data but likely overestimates the age of causal infection
2. **Intervention delivery**: Assumes perfect implementation within coverage constraints (real-world delivery may differ)
3. **Long-term projections**: Uncertainty increases for endpoints far in future (2100)
4. **Cost data**: Economic analysis not included (efficiency measured in health outcomes only)

### Future Directions

1. **Model sensitivity analyses**: Run simulations with varying calibrations refelcting diverse age of causal infection
2. **Implementation research**: Study operational feasibility of high-coverage campaigns
3. **Country comparisons**: Extend analysis to other settings (e.g., Nigeria analysis begun)
4. **Optimal age targeting**: Fine-tune age ranges for maximum efficiency
5. **Real-world validation**: Compare predictions to actual program data as catch-up campaigns are implemented

---

## Technical Details

### Repository Structure

```
hpv_faster_kenya/
├── README.md                          # Brief project description
├── ANALYSIS_NARRATIVE.md             # This document
│
├── data/                             # Input data
│   ├── kenya_age_pyramid.csv         # Population age structure
│   ├── kenya_asr_cancer_incidence.csv # Cancer incidence rates
│   ├── kenya_cancer_cases.csv        # Historical cancer cases
│   └── kenya_cancer_types.csv        # Cancer type distribution
│
├── run_sims.py                       # Kenya simulation setup and calibration
├── run_sims_nigeria.py               # Nigeria simulation (comparison)
├── run_scenarios.py                  # Main scenario runner (catch-up campaigns)
├── analyzers.py                      # Custom analyzers (cohort_cancers, vx_potential)
├── utils.py                          # Helper functions (plotting, distributions)
│
├── plot_tt_comparisons.py            # Generate test-and-treat comparison figures
├── plot_bars.py                      # Bar charts for cancer outcomes
├── plot_age_causal.py                # Age distribution plots (causal infection, HSIL, cancer)
│
├── hpvdna.csv                        # HPV DNA test sensitivity/specificity
├── tx_assigner_faster.csv            # Treatment assignment algorithm
│
├── results/                          # Small output files (CSVs, parameters - committed to repo)
│   ├── kenya_pars.obj                # Calibrated parameters (Kenya)
│   ├── kenya_pars_all.obj            # All calibration trials (Kenya)
│   ├── nigeria_pars.obj              # Calibrated parameters (Nigeria)
│   ├── nigeria_pars_all.obj          # All calibration trials (Nigeria)
│   ├── age_causal_infection_kenya.obj    # Age at causal infection data (Kenya)
│   ├── age_causal_infection_nigeria.obj  # Age at causal infection data (Nigeria)
│   ├── cohort_status_2025.csv        # Cohort HPV status snapshot
│   ├── cohort_simplified_2025.csv    # Simplified cohort status
│   ├── catchup_vax_cohort_cancers_70.csv  # Cancer outcomes, 70% coverage
│   ├── catchup_vax_cohort_cancers_90.csv  # Cancer outcomes, 90% coverage
│   └── combined_impact_data_*.csv    # Combined impact analysis data
│
├── raw_results/                      # Large simulation files (not committed to repo)
│   ├── scens_kenya_70.obj            # 70% coverage scenarios (processed)
│   ├── scens_kenya_90.obj            # 90% coverage scenarios (processed)
│   ├── scens_kenya_70.msim           # 70% coverage raw multi-sim
│   ├── scens_kenya_90.msim           # 90% coverage raw multi-sim
│   ├── sim_kenya.sim                 # Baseline simulation
│   ├── kenya_calib.obj               # Full calibration object
│   ├── kenya_calib_reduced.obj       # Reduced calibration object
│   └── kenya_msim.obj                # Multi-sim results
│
├── figures/                          # Generated plots
│
└── z_archive/                        # Archived figures from analysis runs
    ├── baseline_sim.png
    ├── tt_vs_vax_comparison_kenya_90.png
    ├── catchup_vax_cohorts_kenya_70.png
    ├── catchup_vax_cohorts_kenya_90.png
    ├── efficiency_frontier_kenya_90.png
    ├── efficiency_4panel_kenya_90.png
    ├── cohort_vaccination_eligibility.png
    ├── cohort_simplified_2025.png
    └── [other exploratory plots]
```

### Running the Analysis

#### 1. Calibration (if needed)

```python
python run_sims.py
```

This will calibrate the model to Kenya data and save:
- Parameters to `results/kenya_pars.obj` (small file, for committing)
- Calibration object to `raw_results/kenya_calib.obj` (large file, not committed)

#### 2. Run Scenarios

```python
python run_scenarios.py
```

**Key variables to modify** (in `run_scenarios.py`):
- `location`: 'kenya' or 'nigeria'
- `catchup_cov`: 0.7 or 0.9 (70% or 90% coverage)
- `debug`: Set to 1 for quick test runs, 0 for full analysis
- `do_run`: True to run simulations
- `do_save`: True to save full simulation objects (large files)
- `do_process`: True to process results and create summary objects

**Note**: Full analysis requires HPC resources (scenario count × seed count = large computational burden). With `debug=0`, expect 5-15 minutes on HPC.

#### 3. Generate Plots

```python
# Test-and-treat comparison
import plot_tt_comparisons as ptc
ptc.plot_tt_comparison(location='kenya', coverage=90)

# Bar charts of cancer outcomes
import plot_bars as pb
pb.plot_catchup_scenarios(location='kenya', coverage=90)

# Age distribution of health events
import plot_age_causal as pac
pac.plot_age_causal(location='kenya')
pac.plot_age_causal_violin(location='kenya')
```

### Dependencies

- `hpvsim`: HPV simulation framework
- `sciris`: Scientific utilities (used by hpvsim)
- `numpy`, `pandas`: Data manipulation
- `matplotlib`, `seaborn`: Visualization
- `scipy`: Statistical distributions

Install via:
```bash
pip install hpvsim sciris numpy pandas matplotlib seaborn scipy
```

### Data Format Notes

**Calibrated parameters** (`kenya_pars.obj`):
- Stored as Sciris objdict (pickle format)
- Contains transmission parameters, cancer progression rates, etc.

**Scenario results** (`raw_results/scens_kenya_90.obj`):
- Dictionary keyed by scenario name
- Each entry contains time series and analyzer outputs
- Example keys: 'year', 'cancers', 'asr_cancer_incidence', 'cohort_cancers_10_15', etc.
- Note: Large files stored in `raw_results/` are not committed to version control

**CSV results** (`results/*.csv`):
- Small, human-readable files suitable for committing to version control
- Include cancer outcomes, cohort status, and combined impact data

---

## Conclusions

This analysis demonstrates that **catch-up vaccination campaigns targeting women aged 10-30 years offer the greatest efficiency** for cervical cancer prevention in Kenya. **Combined vaccination and test-and-treat interventions** maximize impact by treating existing disease while preventing future infections.

As Kenya and other countries scale up routine HPV vaccination, catch-up campaigns provide an opportunity to accelerate progress toward cervical cancer elimination by protecting cohorts who aged out of routine programs.

---

## References

1. Brisson M, et al. (2020). Impact of HPV vaccination and cervical screening on cervical cancer elimination: a comparative modelling analysis in 78 low-income and lower-middle-income countries. *The Lancet*.

2. Kenya National Cancer Control Program. (2023). Cervical Cancer Prevention and Control Guidelines.

3. Kenya Demographic and Health Survey (2014). Sexual behavior and HPV data.

4. Peebles K, et al. (2024). HPVsim: An agent-based model for assessing the impact of cervical cancer prevention strategies. *PLOS Computational Biology*.

5. WHO (2020). Global strategy to accelerate the elimination of cervical cancer as a public health problem.

---

## Contact & Collaboration

**Author**: [Your name/team]
**Date**: February 2026
**Repository**: hpv_faster_kenya

For questions about this analysis or collaboration opportunities, please contact [contact information].

---

## Acknowledgments

This work utilizes HPVsim, developed by the Institute for Disease Modeling and collaborators. Data sources include Kenya National Bureau of Statistics, Kenya DHS, and Kenya National Cancer Registry.

Funding: [If applicable]

---

*This document was generated as part of the HPV-Faster strategy analysis for Kenya. All results are model-based projections and should be interpreted in the context of model assumptions and limitations.*
