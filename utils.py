"""
Utilities
"""

# Imports
import os

import hpvsim as hpv
import numpy as np
import pandas as pd
import sciris as sc
from scipy.stats import norm, lognorm


# Canonical metric sets used when extracting CSVs
TS_METRICS_IARC = ['asr_cancer_incidence', 'cancers']   # IARC-comparison TS plots
CUM_METRICS_BOUNDED = ['cancers', 'cancer_deaths']
CUM_METRICS_UNBOUNDED = ['vaccinations']
COHORTS = ['10_15', '15_20', '20_25', '25_30', '30_35', '35_40',
           '40_45', '45_50', '50_55', '55_60']
CUM_START_YEAR = 2025


def set_font(size=None, font='Libertinus Sans'):
    """ Set a custom font """
    sc.fonts(add=sc.thisdir(aspath=True) / 'assets' / 'LibertinusSans-Regular.otf')
    sc.options(font=font, fontsize=size)
    return


def shrink_calib(calib, n_results=100):
    cal = sc.objdict()
    plot_indices = calib.df.iloc[:n_results, 0].values
    valid = plot_indices < len(calib.sim_results)
    plot_indices = plot_indices[valid]
    cal.sim_results = [calib.sim_results[i] for i in plot_indices]
    cal.target_data = calib.target_data
    cal.df = calib.df.iloc[0:n_results, ][valid]
    return cal


def lognorm_params(par1, par2):
    """
    Given the mean and std. dev. of the log-normal distribution, this function
    returns the shape and scale parameters for scipy's parameterization of the
    distribution.
    """
    mean = np.log(par1 ** 2 / np.sqrt(par2 ** 2 + par1 ** 2))  # Computes the mean of the underlying normal distribution
    sigma = np.sqrt(np.log(par2 ** 2 / par1 ** 2 + 1))  # Computes sigma for the underlying normal distribution

    scale = np.exp(mean)
    shape = sigma
    return shape, scale

 
def logn_percentiles_to_pars(x1, p1, x2, p2):
    """ Find the parameters of a lognormal distribution where:
            P(X < p1) = x1
            P(X < p2) = x2
    """
    x1 = np.log(x1)
    x2 = np.log(x2)
    p1ppf = norm.ppf(p1)
    p2ppf = norm.ppf(p2)
    s = (x2 - x1) / (p2ppf - p1ppf)
    mean = ((x1 * p2ppf) - (x2 * p1ppf)) / (p2ppf - p1ppf)
    scale = np.exp(mean)
    return s, scale


def get_debut(sex='f'):
    """
    Read in dataframes taken from DHS and return them in a plot-friendly format,
    optionally saving the distribution parameters
    """
    if sex == 'f':
        x1 = 15
        p1 = 0.184
        x2 = 20
        p2 = 0.858

    else:
        x1 = 15
        p1 = 0.203
        x2 = 20
        p2 = 0.786
    s, scale = logn_percentiles_to_pars(x1, p1, x2, p2)
    rv = lognorm(s=s, scale=scale)

    return rv.mean(), rv.std()



# ---------- CSV extractors (from hpv.MultiSim objdicts → plot-ready CSVs) ----------

def _msim_rows(msim_dict, meta, year_start=None, ts_metrics=None):
    """Generate (ts_rows, cum_rows, coh_rows) for one msim dict.

    `ts_metrics` is a list of metric names to emit as time series (or None = skip TS).
    """
    ts_rows, cum_rows, coh_rows = [], [], []
    ts_metrics = ts_metrics or []

    for scen, mres in msim_dict.items():
        if 'year' in mres:
            years = np.asarray(mres.year)
            mask = slice(None) if year_start is None else years >= year_start
            for metric in ts_metrics:
                if metric not in mres:
                    continue
                r = mres[metric]
                yy = years[mask]
                vv = np.asarray(r)[mask]
                lo = np.asarray(r.low)[mask] if hasattr(r, 'low') else np.full_like(vv, np.nan)
                hi = np.asarray(r.high)[mask] if hasattr(r, 'high') else np.full_like(vv, np.nan)
                for yi, yr in enumerate(yy):
                    ts_rows.append({**meta, 'scenario': scen, 'year': float(yr),
                                    'metric': metric, 'value': float(vv[yi]),
                                    'low': float(lo[yi]), 'high': float(hi[yi])})

            fi = None
            if CUM_START_YEAR in years:
                fi = int(np.where(years == CUM_START_YEAR)[0][0])
            for metric in CUM_METRICS_BOUNDED:
                if metric not in mres or fi is None:
                    continue
                r = mres[metric]
                cum_rows.append({**meta, 'scenario': scen, 'metric': metric,
                                 'value': float(np.sum(r.values[fi:])),
                                 'low': float(np.sum(r.low[fi:])),
                                 'high': float(np.sum(r.high[fi:]))})
            for metric in CUM_METRICS_UNBOUNDED:
                if metric not in mres:
                    continue
                # Sum over all years (matches plot_bars' `mres['vaccinations'].sum()`)
                cum_rows.append({**meta, 'scenario': scen, 'metric': metric,
                                 'value': float(np.sum(np.asarray(mres[metric]))),
                                 'low': float('nan'), 'high': float('nan')})

        for c in COHORTS:
            key = f'cohort_cancers_{c}'
            if key in mres:
                coh_rows.append({**meta, 'scenario': scen, 'cohort': c,
                                 'value': float(mres[key]),
                                 'low': float(mres.get(f'{key}_low', float('nan'))),
                                 'high': float(mres.get(f'{key}_high', float('nan')))})

    return ts_rows, cum_rows, coh_rows


def _to_csv(rows, path):
    pd.DataFrame(rows).to_csv(path, index=False, float_format='%.4f')


def msim_dict_to_csvs(msim_dict, out_stem, resfolder='results',
                      location=None, coverage=None, year_start=None,
                      ts_metrics=None):
    """Emit `<out_stem>_timeseries.csv`, `_cumulative.csv`, `_cohort.csv` from an msim dict."""
    os.makedirs(resfolder, exist_ok=True)
    meta = {}
    if location is not None: meta['location'] = location
    if coverage is not None: meta['coverage'] = coverage

    ts, cum, coh = _msim_rows(msim_dict, meta, year_start=year_start,
                              ts_metrics=ts_metrics)
    if ts: _to_csv(ts, f'{resfolder}/{out_stem}_timeseries.csv')
    if cum: _to_csv(cum, f'{resfolder}/{out_stem}_cumulative.csv')
    if coh: _to_csv(coh, f'{resfolder}/{out_stem}_cohort.csv')


def merge_coverage_csvs(msim_dict_by_cov, out_stem, resfolder='results',
                        location=None, year_start=2015, ts_metrics=None):
    """Merge per-coverage msim dicts into one {stem}_timeseries/cumulative/cohort set.

    `msim_dict_by_cov` maps coverage → msim_dict.
    """
    os.makedirs(resfolder, exist_ok=True)
    all_ts, all_cum, all_coh = [], [], []
    for cov, md in msim_dict_by_cov.items():
        meta = {'coverage': cov}
        if location is not None: meta['location'] = location
        ts, cum, coh = _msim_rows(md, meta, year_start=year_start,
                                  ts_metrics=ts_metrics)
        all_ts += ts; all_cum += cum; all_coh += coh

    if all_ts: _to_csv(all_ts, f'{resfolder}/{out_stem}_timeseries.csv')
    if all_cum: _to_csv(all_cum, f'{resfolder}/{out_stem}_cumulative.csv')
    if all_coh: _to_csv(all_coh, f'{resfolder}/{out_stem}_cohort.csv')


def routine_timing_to_csv(all_results, out_path, catchup_year=2026):
    """Flatten the nested routine-timing results into a summary CSV.

    Columns: routine_start, max_routine_age, scenario, lower_age, upper_age,
             baseline_cancers_2100, scenario_cancers_2100, cancers_averted.
    """
    rows = []
    for routine_key, scenarios in all_results.items():
        routine_start = int(routine_key.split('_')[1])
        max_routine_age = 10 + (catchup_year - routine_start)
        baseline_cancers = float(scenarios['Baseline']['cancers'][-1])
        for scen_label, scen_data in scenarios.items():
            if scen_label == 'Baseline':
                continue
            # "Catchup 10-30" → (10, 30)
            try:
                age_str = scen_label.split('Catchup ')[1]
                lo, hi = [int(x) for x in age_str.split('-')]
            except (IndexError, ValueError):
                lo, hi = None, None
            rows.append({
                'routine_start': routine_start,
                'max_routine_age': max_routine_age,
                'scenario': scen_label,
                'lower_age': lo, 'upper_age': hi,
                'baseline_cancers_2100': baseline_cancers,
                'scenario_cancers_2100': float(scen_data['cancers'][-1]),
                'cancers_averted': baseline_cancers - float(scen_data['cancers'][-1]),
            })
    os.makedirs(os.path.dirname(out_path) or '.', exist_ok=True)
    pd.DataFrame(rows).to_csv(out_path, index=False, float_format='%.4f')


def get_age_causal_df(sim):
    """Extract age_causal_infection analyzer output as a long DataFrame.

    Filters mirror the legacy `plot_age_causal.get_age_causal_df`.
    """
    dt_res = sim.get_analyzer('age_causal_infection')
    pieces = []
    ac_ages = np.array(dt_res.age_causal)
    for arr, event, upper in [
        (dt_res.age_causal, 'Causal\ninfection', 50),
        (dt_res.age_cin, 'HSIL', 65),
        (dt_res.age_cancer, 'Cancer', 90),
    ]:
        ages = np.array(arr)[ac_ages < upper]
        pieces.append(pd.DataFrame({'Age': ages, 'Health event': event}))
    return pd.concat(pieces, ignore_index=True)


def vx_efficacy_to_csv(d, out_path, location=None):
    """scenario × starting_age → cancers, n_cohort (scalar per cell)."""
    rows = []
    for scen, fields in d.items():
        keys = [k for k in fields if k.startswith('cancers_')]
        ages = sorted(int(k.split('_')[1]) for k in keys)
        for age in ages:
            rows.append({
                **({'location': location} if location else {}),
                'scenario': scen, 'start_age': age,
                'cancers': float(fields.get(f'cancers_{age}', 0)),
                'n_cohort': float(fields.get(f'n_cohort_{age}', 0)),
            })
    os.makedirs(os.path.dirname(out_path) or '.', exist_ok=True)
    pd.DataFrame(rows).to_csv(out_path, index=False, float_format='%.4f')


def age_causal_df_to_csv(df, out_path, bin_width=0.5):
    """Save a per-event age_causal DataFrame as a histogram CSV.

    Columns: event, age_mid, count.  The plot side can re-expand to an array of
    ages (via `expand_age_hist`) for violin/KDE plots.
    """
    edges = np.arange(0, 101, bin_width)
    rows = []
    for event, sub in df.groupby('Health event'):
        counts, _ = np.histogram(sub['Age'].values, bins=edges)
        mids = 0.5 * (edges[:-1] + edges[1:])
        for m, c in zip(mids, counts):
            if c > 0:
                rows.append({'event': event, 'age_mid': float(m), 'count': int(c)})
    os.makedirs(os.path.dirname(out_path) or '.', exist_ok=True)
    pd.DataFrame(rows).to_csv(out_path, index=False)


def load_age_causal_hist(path):
    """Load histogram CSV and return {event: DataFrame(age_mid, count)}."""
    df = pd.read_csv(path)
    return {event: sub.reset_index(drop=True) for event, sub in df.groupby('event')}


def expand_age_hist(hist_df):
    """Re-expand histogram rows to an array of ages (using bin midpoints)."""
    return np.repeat(hist_df.age_mid.values, hist_df['count'].values.astype(int))


# ---------- CSV loaders (plot-side) ----------

def load_scens(resfolder, stem='scens_kenya'):
    ts_path = f'{resfolder}/{stem}_timeseries.csv'
    cum_path = f'{resfolder}/{stem}_cumulative.csv'
    coh_path = f'{resfolder}/{stem}_cohort.csv'
    ts = pd.read_csv(ts_path) if os.path.exists(ts_path) else None
    cum = pd.read_csv(cum_path) if os.path.exists(cum_path) else None
    coh = pd.read_csv(coh_path) if os.path.exists(coh_path) else None
    return ts, cum, coh


def get_cum(cum_df, scenario, metric, coverage=None, location=None):
    mask = (cum_df.scenario == scenario) & (cum_df.metric == metric)
    if coverage is not None and 'coverage' in cum_df.columns:
        mask &= (cum_df.coverage == coverage)
    if location is not None and 'location' in cum_df.columns:
        mask &= (cum_df.location == location)
    sub = cum_df[mask]
    if sub.empty:
        return 0.0, float('nan'), float('nan')
    row = sub.iloc[0]
    return float(row.value), float(row.low), float(row.high)


def get_ts(ts_df, scenario, metric, coverage=None, location=None):
    mask = (ts_df.scenario == scenario) & (ts_df.metric == metric)
    if coverage is not None and 'coverage' in ts_df.columns:
        mask &= (ts_df.coverage == coverage)
    if location is not None and 'location' in ts_df.columns:
        mask &= (ts_df.location == location)
    sub = ts_df[mask].sort_values('year')
    return sub.year.values, sub.value.values, sub.low.values, sub.high.values


def get_cohort(coh_df, scenario, cohort, coverage=None, location=None):
    mask = (coh_df.scenario == scenario) & (coh_df.cohort == cohort)
    if coverage is not None and 'coverage' in coh_df.columns:
        mask &= (coh_df.coverage == coverage)
    if location is not None and 'location' in coh_df.columns:
        mask &= (coh_df.location == location)
    sub = coh_df[mask]
    if sub.empty:
        return 0.0, float('nan'), float('nan')
    row = sub.iloc[0]
    return float(row.value), float(row.low), float(row.high)


# ---------- Convenience loaders (CSV-first, fall back to raw_results .obj) ----------

def load_scens_obj(location, coverage, resfolder='results/v2.2.6_baseline'):
    """Return an msim-dict-like for scens_{location}_{coverage}.

    Prefers plot-ready CSVs under `resfolder` (via MsimDictFromCSV); if those
    aren't present, falls back to the heavy `.obj` under raw_results/.
    """
    if os.path.exists(f'{resfolder}/scens_{location}_cohort.csv'):
        return MsimDictFromCSV(resfolder, f'scens_{location}',
                               coverage=coverage, location=location, include_ts=False)
    return sc.loadobj(f'raw_results/scens_{location}_{coverage}.obj')


def load_iarc_obj(location, coverage, resfolder='results/v2.2.6_baseline'):
    """Return an msim-dict-like for iarc_comparison_{location}_{coverage}."""
    if os.path.exists(f'{resfolder}/iarc_{location}_{coverage}_cohort.csv'):
        return MsimDictFromCSV(resfolder, f'iarc_{location}_{coverage}',
                               coverage=coverage, location=location, include_ts=True)
    return sc.loadobj(f'raw_results/iarc_comparison_{location}_{coverage}.obj')


# ---------- Obj adapter: lets plot scripts read CSVs with minimal code change ----------

class _ResultTS:
    """Mimics hpv.Result: has .values, .low, .high (and supports indexing + .sum())."""
    __slots__ = ('values', 'low', 'high')
    def __init__(self, values, low, high):
        self.values = np.asarray(values, dtype=float)
        self.low = np.asarray(low, dtype=float)
        self.high = np.asarray(high, dtype=float)
    def sum(self): return float(self.values.sum())
    def __getitem__(self, idx): return self.values[idx]
    def __len__(self): return len(self.values)
    def __iter__(self): return iter(self.values)
    def __float__(self): return float(self.values)


class _CumulativeVal:
    """Wraps a single cumulative scalar so `.sum()` returns it.

    Original plot code does `mres['vaccinations'].sum()` on a per-year array; here
    the CSV already stores the sum, so `.sum()` is a no-op.
    """
    __slots__ = ('value', 'low', 'high')
    def __init__(self, value, low=np.nan, high=np.nan):
        self.value = float(value)
        self.low = float(low)
        self.high = float(high)
    def sum(self): return self.value
    @property
    def values(self): return np.array([self.value])
    def __float__(self): return self.value
    def __getitem__(self, idx):
        if idx in (slice(None), 0, -1) or (isinstance(idx, slice) and idx.start in (None, 0)):
            return self.value
        return self.value
    def __add__(self, other): return self.value + float(other)
    def __radd__(self, other): return float(other) + self.value
    def __sub__(self, other): return self.value - float(other)
    def __rsub__(self, other): return float(other) - self.value


class _ScenResult:
    """Per-scenario view built lazily from cumulative/cohort/ts dataframes."""
    _COHORT_PREFIX = 'cohort_cancers_'

    def __init__(self, scenario, cum_df, coh_df, ts_df=None, coverage=None, location=None):
        self.scenario = scenario
        self._cum = cum_df
        self._coh = coh_df
        self._ts = ts_df
        self._cov = coverage
        self._loc = location

    def _filter(self, df):
        mask = df.scenario == self.scenario
        if self._cov is not None and 'coverage' in df.columns:
            mask &= (df.coverage == self._cov)
        if self._loc is not None and 'location' in df.columns:
            mask &= (df.location == self._loc)
        return df[mask]

    def _cohort_key(self, key):
        """Return (cohort, column) for cohort_cancers_X[_low|_high], or None."""
        if not key.startswith(self._COHORT_PREFIX):
            return None
        suf = key[len(self._COHORT_PREFIX):]
        if suf.endswith('_low'):
            return suf[:-4], 'low'
        if suf.endswith('_high'):
            return suf[:-5], 'high'
        return suf, 'value'

    def __contains__(self, key):
        if key == 'year':
            return self._ts is not None and not self._filter(self._ts).empty
        ck = self._cohort_key(key)
        if ck is not None:
            cohort, _ = ck
            return ((self._coh is not None) and
                    not self._filter(self._coh)[self._filter(self._coh).cohort == cohort].empty)
        if self._cum is not None:
            if not self._filter(self._cum[self._cum.metric == key]).empty:
                return True
        if self._ts is not None:
            if not self._filter(self._ts[self._ts.metric == key]).empty:
                return True
        return False

    def __getitem__(self, key):
        if key == 'year':
            sub = self._filter(self._ts)
            return np.sort(sub.year.unique())
        ck = self._cohort_key(key)
        if ck is not None:
            cohort, col = ck
            sub = self._filter(self._coh)
            sub = sub[sub.cohort == cohort]
            if sub.empty:
                raise KeyError(key)
            return float(sub.iloc[0][col])
        # TS metric
        if self._ts is not None:
            sub = self._filter(self._ts[self._ts.metric == key]).sort_values('year')
            if not sub.empty:
                return _ResultTS(sub.value.values, sub.low.values, sub.high.values)
        # Cumulative scalar
        if self._cum is not None:
            sub = self._filter(self._cum[self._cum.metric == key])
            if not sub.empty:
                row = sub.iloc[0]
                return _CumulativeVal(row.value, row.low, row.high)
        raise KeyError(key)

    def get(self, key, default=None):
        try:
            return self[key]
        except KeyError:
            return default


class MsimDictFromCSV:
    """Drop-in replacement for an msim_dict loaded from .obj — reads CSVs instead.

    Usage:
        msim_dict = MsimDictFromCSV(resfolder, 'scens_kenya', coverage=70, include_ts=False)
        msim_dict = MsimDictFromCSV(resfolder, 'iarc_kenya_70', include_ts=True)
    """
    def __init__(self, resfolder, stem, coverage=None, location=None, include_ts=False):
        ts, cum, coh = load_scens(resfolder, stem)
        self._ts = ts if include_ts else None
        self._cum = cum
        self._coh = coh
        self._cov = coverage
        self._loc = location

    def _filter(self, df):
        mask = np.ones(len(df), dtype=bool)
        if self._cov is not None and 'coverage' in df.columns:
            mask &= (df.coverage == self._cov).values
        if self._loc is not None and 'location' in df.columns:
            mask &= (df.location == self._loc).values
        return df[mask]

    def __contains__(self, scenario):
        for df in (self._cum, self._coh, self._ts):
            if df is not None and (self._filter(df).scenario == scenario).any():
                return True
        return False

    def __getitem__(self, scenario):
        return _ScenResult(scenario, self._cum, self._coh, ts_df=self._ts,
                           coverage=self._cov, location=self._loc)

    def keys(self):
        keys = set()
        for df in (self._cum, self._coh, self._ts):
            if df is not None:
                keys |= set(self._filter(df).scenario.unique())
        return list(keys)

    def items(self):
        for k in self.keys():
            yield k, self[k]


def plot_ts(ax, ts_df, scenario, metric, start_year, end_year,
            color, ls='-', label=None, smooth=True, add_bounds=True,
            coverage=None, location=None):
    years, best, low, high = get_ts(ts_df, scenario, metric,
                                    coverage=coverage, location=location)
    mask = (years >= start_year) & (years <= end_year)
    years, best, low, high = years[mask], best[mask], low[mask], high[mask]

    if smooth and len(years) >= 5:
        best = np.convolve(best, np.ones(5), 'valid') / 5
        low = np.convolve(low, np.ones(5), 'valid') / 5
        high = np.convolve(high, np.ones(5), 'valid') / 5
        years = years[4:]

    ax.plot(years, best, color=color, label=label, ls=ls)
    if metric == 'asr_cancer_incidence':
        below = np.where(best < 4)[0]
        if len(below):
            print(f'{label} elim year: {years[below[0]]}')
    if add_bounds:
        ax.fill_between(years, low, high, alpha=0.1, color=color)
    return ax


def plot_single(ax, mres, to_plot, si, ei, color, ls='-', label=None, smooth=True):
    years = mres.year[si:ei]
    best = mres[to_plot][si:ei]
    low = mres[to_plot].low[si:ei]
    high = mres[to_plot].high[si:ei]

    if smooth:
        best = np.convolve(list(best), np.ones(5), "valid")/5
        low = np.convolve(list(low), np.ones(5), "valid")/5
        high = np.convolve(list(high), np.ones(5), "valid")/5
        years = years[4:]

    ax.plot(years, best, color=color, label=label, ls=ls)

    if to_plot == 'asr_cancer_incidence':
        try:
            elim_year = sc.findfirst(best<4)
            print(f'{label} elim year: {years[elim_year]}')
        except:
            print(f'{label} not eliminated')

    ax.fill_between(years, low, high, alpha=0.1, color=color)
    # ax.set_yscale('log')

    # Add horizontal line at 4
    ax.axhline(4, color='k', ls='--', lw=0.5)
    return ax



# %% Run as a script
if __name__ == '__main__':

    for sex in ['f', 'm']:
        mean, std = get_debut(sex)
        print(f'Mean debut age ({sex}): {mean:.2f}, std: {std:.2f}')

