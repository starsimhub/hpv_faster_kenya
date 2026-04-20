# Changelog

## 2026-04-19 — HPVsim v2.2.6 lift

- Split the workflow so heavy sims run on a VM via the existing `run_*.py` scripts (saving `.obj` / `.msim` to `raw_results/`) and then extract lightweight plot-ready CSVs under `results/`. Local plot scripts read the CSVs (not pickles).
- Added `utils.MsimDictFromCSV` adapter so legacy msim-dict-style code paths work unchanged when fed CSVs. `ut.load_scens_obj` / `ut.load_iarc_obj` prefer CSVs, falling back to raw `.obj` when CSVs aren't present.
- Emitted CSVs per producer:
  - `scens_{loc}_cumulative.csv` + `scens_{loc}_cohort.csv` (merged across 30/40/50/70/90 coverages)
  - `iarc_{loc}_{cov}_timeseries.csv` + `_cumulative.csv` + `_cohort.csv`
  - `coverage_sweep_{loc}_cohort.csv` (+ existing summary CSV)
  - `vx_efficacy_{loc}.csv`
  - `routine_timing_{loc}.csv`
  - `age_causal_{loc}.csv` (0.5-year-bin histogram, replaces the 1.2 MB pickle)
- Refactored plot scripts to accept `--resfolder` / `resfolder=` and default to `results/v2.2.6_baseline/`: `plot_bars.py`, `plot_tt_comparisons.py`, `plot_age_causal.py`, `plot_routine_timing_results.py`.
- Froze `results/v2.2.6_baseline/` (~900 KB total) as the committed baseline.
- Added `compare_baselines.py` for cross-version comparison (v2.3 / v3.0 ready).
- Dropped orphaned `plot_fig3_hiv.py`-style files: `run_bridging_analysis.py` (unused) and `email.md` (personal); deleted the heavy `age_causal_infection_*.obj` pickles in `results/` (replaced by CSVs).
- Added `.gitignore` excluding transient `.obj`/`.sim`/`.msim`/`.zip` in top-level `results/`; `raw_results/` stays ignored.
- Tier 3 bundle: LICENSE (MIT), this CHANGELOG, README updated with VM/local workflow and v2.3/v3.0 migration instructions.

## Earlier

- Pre-uplift work on the `sterilizing-vx` branch carried into this uplift: new functions `plot_benefit_distribution_combined` (plot_bars.py) and `plot_combined_impact_2panel` (plot_tt_comparisons.py).
- Added IARC comparison workflow (`run_iarc_comparison.py`, `plot_tt_comparisons.plot_tt_comparison`) and routine-timing sensitivity analysis (`run_routine_timing_scenarios.py`, `plot_routine_timing_results.py`).
- Core analyses: catch-up vaccination + test-and-treat across 10-60 age cohorts in Kenya (and Nigeria comparator), with age-of-causal-infection diagnostics.
