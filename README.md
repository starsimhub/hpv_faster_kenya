# HPV-Faster Kenya Analysis

**Comparative effectiveness of catch-up HPV vaccination and test-and-treat strategies in Kenya (and Nigeria comparator)**

This repository contains agent-based modelling analysis of HPV-Faster strategies for cervical cancer prevention in Kenya. The analysis compares catch-up vaccination campaigns versus test-and-treat interventions across age cohorts (10-60 years), identifies optimal strategies for women who missed routine adolescent HPV vaccination, and compares results against a Nigerian setting where the age distribution of causal HPV infection differs.

**Results in this repository were produced with HPVsim v2.2.6.** Plot-ready baselines live in [`results/v2.2.6_baseline/`](results/v2.2.6_baseline/); all plot scripts default to reading from this folder.

**📄 [Read the analysis narrative](ANALYSIS_NARRATIVE.md)** — methods, results, and policy implications.

## Installation

```bash
pip install hpvsim==2.2.6 seaborn scipy
```

Python 3.9+.

## Workflow: heavy sims on VM, plots locally

Heavy simulations run on a VM (produce `.obj` / `.msim` in `raw_results/`, plus lightweight plot-ready CSVs in `results/`); local plot scripts read the CSVs (no pickles). Transient top-level `results/*.obj` / `.sim` / `.msim` / `.zip` / `.csv` are gitignored. Plot-ready CSVs frozen under `results/v<version>_baseline/` **are** committed.

### Heavy (VM) scripts

| Script | What it produces |
|---|---|
| `run_sims.py` | Single sim + calibration + `age_causal_{loc}.csv` |
| `run_sims_nigeria.py` | Same as above, Nigeria |
| `run_scenarios.py` | All catch-up scenarios × coverage levels → `scens_{loc}_cumulative.csv` + `_cohort.csv` |
| `run_iarc_comparison.py` | IARC-style scenario set (routine vs catch-up vs S&T) → `iarc_{loc}_{cov}_*.csv` |
| `run_coverage_sweep.py` | Fine-grained coverage sweep → `coverage_sweep_{loc}*.csv` |
| `run_vx_efficacy_comparison.py` | Vaccination-age efficacy → `vx_efficacy_{loc}.csv` |
| `run_routine_timing_scenarios.py` | Routine-start-year sensitivity → `routine_timing_{loc}.csv` |

Each `run_*.py` saves heavy `.obj` to `raw_results/` AND emits plot-ready CSVs to `results/` as a final step.

### Plot scripts (local)

| Script | What it renders |
|---|---|
| `plot_bars.py` | Stacked bars, benefit distribution, cohort trajectories across coverages |
| `plot_tt_comparisons.py` | IARC comparison, efficiency frontier, combined impact, heatmaps |
| `plot_age_causal.py` | Age distribution of causal infection / HSIL / cancer (box, bar, violin) |
| `plot_routine_timing_results.py` | Marginal benefit vs routine-start-year |
| `plot_routine_timing_sensitivity.py` | Synthetic-data illustration (no sim dependency) |
| `plot_kenya_nigeria_comparison.py` | Kenya vs Nigeria cancer-by-age from data CSVs |

All plot functions accept `resfolder=` (default `results/v2.2.6_baseline`) and internally call `ut.load_scens_obj()` / `ut.load_iarc_obj()` which read the CSVs via the `MsimDictFromCSV` adapter in `utils.py`.

### Other files

- `interventions / analyzers.py` / `run_sims*.py` — sim-building, calibration, custom cohort analyzers
- `compare_baselines.py` — cross-version comparison (see below)
- `utils.py` — CSV loaders, `MsimDictFromCSV` adapter, `plot_ts` helper

## Reproducing the figures

```bash
git clone git@github.com:amath-idm/hpv_faster_kenya.git
cd hpv_faster_kenya
pip install hpvsim==2.2.6 seaborn scipy

# Render all figures from the committed v2.2.6 baseline
python -c "
import plot_bars as pb, plot_tt_comparisons as pt, plot_age_causal as pac
pb.plot_stacked_bars('kenya', coverage=90)
pb.plot_stacked_bars_pct_2panel('kenya')
pb.plot_benefit_distribution_combined('kenya')
pb.plot_cohorts_all_coverages('kenya')
pt.plot_tt_comparison('kenya', coverage=70)
pt.plot_efficiency_frontier('kenya', coverage=90)
pt.plot_combined_impact('kenya', coverage=90)
pt.plot_heatmaps('kenya', coverage=90)
pac.plot_age_causal('kenya')
pac.plot_age_causal_violin('kenya')
pac.plot_age_causal_bar('kenya')
"
```

## Cross-version comparison (v2.2.6 → v2.3 → v3.0)

When a new HPVsim version ships, regenerate the baseline in a clean env and compare side-by-side:

```bash
# 1. On a VM, in a clean env pinned to the new version
conda create -n hpvsim230 python=3.11 -y && conda activate hpvsim230
pip install hpvsim==2.3.0 seaborn scipy

# 2. Re-run the heavy scripts (they also emit plot-ready CSVs)
python run_sims.py
python run_scenarios.py
python run_iarc_comparison.py
# (also run_coverage_sweep / run_vx_efficacy_comparison / run_routine_timing as needed)

# 3. Freeze the fresh CSVs into a versioned baseline dir
mkdir -p results/v2.3.0_baseline
cp results/scens_*.csv results/iarc_*.csv results/coverage_sweep_*.csv \
   results/vx_efficacy_*.csv results/routine_timing_*.csv results/age_causal_*.csv \
   results/v2.3.0_baseline/

# 4. Commit + push, then locally compare
python compare_baselines.py --baselines v2.2.6_baseline v2.3.0_baseline \
                            --location kenya --coverage 70
```

`compare_baselines.py` prints cumulative-cancers + cohort-cancers tables across baselines and emits `figures/compare_baselines.png`. Extends trivially to v3.0 by appending the new baseline name to `--baselines`.

## Inputs

- `data/` — Kenya + Nigeria demographic, cancer incidence (GLOBOCAN), age pyramid inputs
- `hpvdna.csv`, `tx_assigner_faster.csv` — screening + test-and-treat specifications

## Citation

If you use this code, please cite the analysis narrative in `ANALYSIS_NARRATIVE.md`. Manuscript in preparation.

## Contact

For questions or collaboration: info@hpvsim.org

## Further information

See [hpvsim.org](https://hpvsim.org) and [docs.hpvsim.org](https://docs.hpvsim.org).
