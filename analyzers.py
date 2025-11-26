"""
Define custom analyzers
"""

import numpy as np
import pandas as pd
import sciris as sc
import hpvsim as hpv
import hpvsim.utils as hpu

 
class cohort_cancers(hpv.Analyzer):
    def __init__(self, cohort_age=None, start=None, **kwargs):
        super().__init__(**kwargs)
        self.start = start or 2026
        self.cohort_age = cohort_age or [10, 15]
        self.years = None
        self.results = None
        return

    def initialize(self, sim):
        super().initialize()
        self.si = sc.findfirst(sim.res_yearvec, self.start)
        self.npts = len(sim.res_yearvec[self.si:])
        self.years = sim.res_yearvec[self.si:]
        self.results = np.zeros(self.npts)
        return

    def apply(self, sim):
        if sim.yearvec[sim.t] >= self.start:
            li = np.floor(sim.yearvec[sim.t])
            idx = sc.findfirst(self.years, li)
            ppl = sim.people

            time_elapsed = sim.yearvec[sim.t] - self.start
            current_age_range = [self.cohort_age[0]+time_elapsed, self.cohort_age[1]+time_elapsed]

            cic = (ppl.date_cancerous == sim.t) & (ppl.age >= current_age_range[0]) & (ppl.age <= current_age_range[1])
            if cic.any():
                self.results[idx] += sum(ppl.scale[hpu.true(cic)])

        return

    @staticmethod
    def reduce(analyzers, use_mean=False, quantiles=None):
        # Process quantiles
        if quantiles is None:
            quantiles = {'low': 0.1, 'high': 0.9}
        if not isinstance(quantiles, dict):
            try:
                quantiles = {'low':float(quantiles[0]), 'high':float(quantiles[1])}
            except Exception as E:
                errormsg = f'Could not figure out how to convert {quantiles} into a quantiles object: must be a dict with keys low, high or a 2-element array ({str(E)})'
                raise ValueError(errormsg)

        # Get base analyzer properties and copy them into reduced analyzer
        base_analyzer = analyzers[0]
        reduced_analyzer = sc.dcp(base_analyzer)
        ashape = base_analyzer.results.shape  # Figure out dimensions
        new_ashape = ashape + (len(analyzers),)
        raw = np.zeros(new_ashape)

        # Pull out results for each analyzer
        for a, analyzer in enumerate(analyzers):
            raw[:, a] = analyzer.results

        # Get quantiles
        reduced_analyzer.raw = raw
        reduced_analyzer.results = np.quantile(raw, q=0.5, axis=-1)
        reduced_analyzer.low  = np.quantile(raw, q=quantiles['low'], axis=-1)
        reduced_analyzer.high = np.quantile(raw, q=quantiles['high'], axis=-1)

        # Do sums
        sums = raw.sum(axis=0)
        reduced_analyzer.cum_cancers_best = np.quantile(sums, q=0.5)
        reduced_analyzer.cum_cancers_low  = np.quantile(sums, q=quantiles['low'])
        reduced_analyzer.cum_cancers_high = np.quantile(sums, q=quantiles['high'])

        return reduced_analyzer


class AFS(hpv.Analyzer):
    def __init__(self, bins=None, cohort_starts=None, **kwargs):
        super().__init__(**kwargs)
        self.bins = bins or np.arange(12,31,1)
        self.cohort_starts = cohort_starts
        self.binspan = self.bins[-1]-self.bins[0]

    def initialize(self, sim):
        super().initialize()
        if self.cohort_starts is None:
            first_cohort = sim['start'] + sim['burnin'] - 5
            last_cohort = sim['end']-self.binspan
            self.cohort_starts = sc.inclusiverange(first_cohort, last_cohort)
            self.cohort_ends = self.cohort_starts+self.binspan
            self.n_cohorts = len(self.cohort_starts)
            self.cohort_years = np.array([sc.inclusiverange(i,i+self.binspan) for i in self.cohort_starts])

        self.prop_active_f = np.zeros((self.n_cohorts,self.binspan+1))
        self.prop_active_m = np.zeros((self.n_cohorts,self.binspan+1))

    def apply(self, sim):
        if sim.yearvec[sim.t] in self.cohort_years:
            cohort_inds, bin_inds = sc.findinds(self.cohort_years, sim.yearvec[sim.t])
            for ci,cohort_ind in enumerate(cohort_inds):
                bin_ind = bin_inds[ci]
                bin = self.bins[bin_ind]

                conditions_f = sim.people.is_female * sim.people.alive * (sim.people.age >= (bin-1)) * (sim.people.age < bin) * sim.people.level0
                denom_inds_f = hpu.true(conditions_f)
                num_conditions_f = conditions_f * (sim.people.n_rships.sum(axis=0)>0)
                num_inds_f = hpu.true(num_conditions_f)
                self.prop_active_f[cohort_ind,bin_ind] = len(num_inds_f)/len(denom_inds_f)

                conditions_m = ~sim.people.is_female * sim.people.alive * (sim.people.age >= (bin-1)) * (sim.people.age < bin)
                denom_inds_m = hpu.true(conditions_m)
                num_conditions_m = conditions_m * (sim.people.n_rships.sum(axis=0)>0)
                num_inds_m = hpu.true(num_conditions_m)
                self.prop_active_m[ci,bin_ind] = len(num_inds_m)/len(denom_inds_m)
        return


class prop_married(hpv.Analyzer):
    def __init__(self, bins=None, years=None, includelast=True, yearstride=5, binspan=5, **kwargs):
        super().__init__(**kwargs)
        self.bins = bins or np.arange(15,50,binspan)
        self.years = years
        self.dfs = sc.autolist()
        self.df = None
        self.includelast = includelast
        self.yearstride = yearstride
        self.binspan = binspan

    def initialize(self, sim):
        super().initialize()
        if self.years is None:
            start = sim['start'] + sim['burnin']
            end = sim['end']
            self.years = np.arange(start, end, self.yearstride)
            if self.includelast:
                if end not in self.years:
                    self.years = np.append(self.years, end)

    def apply(self, sim):
        if sim.yearvec[sim.t] in self.years:

            conditions = dict()
            for ab in self.bins:
                conditions[ab] = (sim.people.age >= ab) & (sim.people.age < ab+self.binspan) & sim.people.alive & sim.people.is_female & sim.people.level0

            prop_married = sc.autolist()
            for age_cond in conditions.values():
                num_condition = age_cond & (sim.people.current_partners[0,:]>0)
                prop_married += len(hpu.true(num_condition))/len(hpv.true(age_cond))

            d = dict(age=self.bins, val=prop_married)
            df = pd.DataFrame().from_dict(d)
            df['year'] = sim.yearvec[sim.t]
            self.dfs += df

    def finalize(self, sim):
        self.df = pd.concat(self.dfs)
