"""
Define custom analyzers
"""

import numpy as np
import pandas as pd
import sciris as sc
import hpvsim as hpv
import hpvsim.utils as hpu
import matplotlib.pyplot as plt

 
class cohort_cancers(hpv.Analyzer):
    def __init__(self, cohort_age=None, label=None, start=None, **kwargs):
        super().__init__(**kwargs)
        self.start = start or 2026
        self.cohort_age = cohort_age or [10, 15]
        self.years = None
        self.results = None
        self.label = label or f'cohort_cancers_{self.cohort_age[0]}_{self.cohort_age[1]}'
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



class vx_potential(hpv.Analyzer):
    """
    Analyzer to capture HPV exposure and disease status by age cohort
    at a specific timepoint (e.g., 2025).
    
    This helps determine what proportion of each cohort could benefit
    from vaccination interventions.
    
    Args:
        timepoint: Year to capture cohort status (default: 2025)
        cohort_edges: Age bin edges for cohorts (default: [0, 10, 15, 20, ..., 60])
    """
    
    def __init__(self, timepoint=2025, cohort_edges=None, **kwargs):
        super().__init__(**kwargs)
        self.timepoint = timepoint
        
        # Default cohort edges matching your decomposition plot
        if cohort_edges is None:
            self.cohort_edges = [0, 10, 15] + list(range(20, 65, 5))
        else:
            self.cohort_edges = cohort_edges
            
        self.cohort_labels = []
        for i in range(len(self.cohort_edges) - 1):
            lower = self.cohort_edges[i]
            upper = self.cohort_edges[i + 1]
            self.cohort_labels.append(f'{lower}-{upper}')
        
        self.n_cohorts = len(self.cohort_labels)
        self.captured = False
        self.results = {}
        
    def initialize(self, sim):
        super().initialize(sim)
        self.initialized = True
        
    def apply(self, sim):
        """
        Capture cohort status when we reach the target timepoint
        """
        # Only capture once at the specified timepoint
        current_year = sim.yearvec[sim.t]
        if not self.captured and current_year >= self.timepoint:
            self.capture_cohort_status(sim)
            self.captured = True
    
    def capture_cohort_status(self, sim):
        """
        Capture detailed status of each cohort
        """
        people = sim.people
        
        # Filter for females only
        sex_filter = people.is_female
        
        # Get alive people
        alive = people.alive
        age = people.age
        target_pop = alive & sex_filter
        
        # Get genotypes from sim parameters
        genotypes = sim.pars['genotypes']
        n_genotypes = len(genotypes)
        
        # Initialize results dictionary
        for cohort_label in self.cohort_labels:
            self.results[cohort_label] = {
                'n_total': 0,
                'n_never_infected': 0,
                'n_currently_infected': 0,
                'n_previously_infected': 0,
                'n_cin': 0,
                'n_cancerous': 0,
                'n_sexually_active': 0,
                'n_zero_partners': 0,
                'n_one_partner': 0,
                'n_multiple_partners': 0,
                'mean_partners': 0,
            }
        
        # Process each cohort
        for cohort_i, cohort_label in enumerate(self.cohort_labels):
            lower_age = self.cohort_edges[cohort_i]
            upper_age = self.cohort_edges[cohort_i + 1]
            
            # Get people in this age cohort
            cohort_mask = target_pop & (age >= lower_age) & (age < upper_age)
            
            n_in_cohort = np.sum(cohort_mask)
            self.results[cohort_label]['n_total'] = n_in_cohort
            
            if n_in_cohort == 0:
                continue
            
            # --- HPV Infection Status ---
            # people.infected is shape (n_genotypes, n_people)
            # Check if currently infected with any genotype
            currently_infected = np.any(people.infected, axis=0)
            
            # people.date_exposed is shape (n_genotypes, n_people)
            # Check if ever exposed to any genotype (date_exposed > -1 means exposed)
            ever_exposed = np.any(people.date_exposed > -1, axis=0)
            
            # People who were previously but not currently infected
            past_only_infected = ever_exposed & ~currently_infected
            
            self.results[cohort_label]['n_currently_infected'] = np.sum(currently_infected[cohort_mask])
            self.results[cohort_label]['n_previously_infected'] = np.sum(past_only_infected[cohort_mask])
            self.results[cohort_label]['n_never_infected'] = n_in_cohort - np.sum(ever_exposed[cohort_mask])
            
            # --- Disease Status (CIN and Cancer) ---
            # people.cin is shape (n_genotypes, n_people)
            # Check if has CIN from any genotype
            has_cin = np.any(people.cin, axis=0)
            
            # people.cancerous is shape (n_genotypes, n_people)
            # Check if has cancer from any genotype
            has_cancer = np.any(people.cancerous, axis=0)
            
            self.results[cohort_label]['n_cin'] = np.sum(has_cin[cohort_mask])
            self.results[cohort_label]['n_cancerous'] = np.sum(has_cancer[cohort_mask])
            
            # --- Sexual Activity and Partner History ---
            # people.n_rships has 2 rows: [marital, casual]
            # Sum across both types for total partners
            n_partners = people.n_rships.sum(axis=0)
            
            cohort_partners = n_partners[cohort_mask]
            
            self.results[cohort_label]['n_sexually_active'] = np.sum(cohort_partners > 0)
            self.results[cohort_label]['n_zero_partners'] = np.sum(cohort_partners == 0)
            self.results[cohort_label]['n_one_partner'] = np.sum(cohort_partners == 1)
            self.results[cohort_label]['n_multiple_partners'] = np.sum(cohort_partners > 1)
            self.results[cohort_label]['mean_partners'] = np.mean(cohort_partners) if n_in_cohort > 0 else 0
    
    def finalize(self, sim):
        """
        Calculate proportions and create summary dataframe
        """
        super().finalize(sim)
        
        # Convert to DataFrame for easier analysis
        rows = []
        for cohort_label in self.cohort_labels:
            res = self.results[cohort_label]
            n_total = res['n_total']
            
            if n_total > 0:
                row = {
                    'cohort': cohort_label,
                    'n_total': n_total,
                    # Infection status (proportions)
                    'pct_never_infected': 100 * res['n_never_infected'] / n_total,
                    'pct_currently_infected': 100 * res['n_currently_infected'] / n_total,
                    'pct_previously_infected': 100 * res['n_previously_infected'] / n_total,
                    # Disease status (proportions)
                    'pct_cin': 100 * res['n_cin'] / n_total,
                    'pct_cancerous': 100 * res['n_cancerous'] / n_total,
                    # Sexual activity (proportions)
                    'pct_sexually_active': 100 * res['n_sexually_active'] / n_total,
                    'pct_zero_partners': 100 * res['n_zero_partners'] / n_total,
                    'pct_one_partner': 100 * res['n_one_partner'] / n_total,
                    'pct_multiple_partners': 100 * res['n_multiple_partners'] / n_total,
                    'mean_partners': res['mean_partners'],
                }
                # Add raw counts
                row.update({f'n_{k}': v for k, v in res.items() if k.startswith('n_')})
                rows.append(row)
        
        self.summary_df = pd.DataFrame(rows)
    
    def plot(self, figsize=(16, 10)):
        """
        Create visualization of cohort vaccination eligibility
        """
        if not hasattr(self, 'summary_df'):
            print("No data to plot")
            return
        
        fig, axes = plt.subplots(2, 2, figsize=figsize)
        fig.suptitle(f'Cohort Status in {self.timepoint}: Vaccination Eligibility Assessment', 
                     fontsize=14, fontweight='bold')
        
        cohorts = self.summary_df['cohort'].values
        x_pos = np.arange(len(cohorts))
        
        # Panel A: HPV Exposure Status
        ax_a = axes[0, 0]
        width = 0.25
        ax_a.bar(x_pos - width, self.summary_df['pct_never_infected'], width, 
                label='Never infected', color='#2ecc71')
        ax_a.bar(x_pos, self.summary_df['pct_currently_infected'], width, 
                label='Currently infected', color='#e74c3c')
        ax_a.bar(x_pos + width, self.summary_df['pct_previously_infected'], width, 
                label='Previously infected (cleared)', color='#f39c12')
        
        ax_a.set_xlabel('Age cohort', fontweight='bold')
        ax_a.set_ylabel('Percentage', fontweight='bold')
        ax_a.set_title('A) HPV Infection Status', loc='left', fontweight='bold')
        ax_a.set_xticks(x_pos)
        ax_a.set_xticklabels(cohorts, rotation=45, ha='right')
        ax_a.legend()
        ax_a.grid(alpha=0.3, axis='y')
        
        # Panel B: Disease Status
        ax_b = axes[0, 1]
        width = 0.35
        ax_b.bar(x_pos - width/2, self.summary_df['pct_cin'], width, 
                label='CIN', color='#9b59b6')
        ax_b.bar(x_pos + width/2, self.summary_df['pct_cancerous'], width, 
                label='Cancer', color='#c0392b')
        
        ax_b.set_xlabel('Age cohort', fontweight='bold')
        ax_b.set_ylabel('Percentage', fontweight='bold')
        ax_b.set_title('B) Disease Status', loc='left', fontweight='bold')
        ax_b.set_xticks(x_pos)
        ax_b.set_xticklabels(cohorts, rotation=45, ha='right')
        ax_b.legend()
        ax_b.grid(alpha=0.3, axis='y')
        
        # Panel C: Sexual Activity
        ax_c = axes[1, 0]
        width = 0.3
        ax_c.bar(x_pos - width, self.summary_df['pct_zero_partners'], width, 
                label='0 partners', color='#2ecc71')
        ax_c.bar(x_pos, self.summary_df['pct_one_partner'], width, 
                label='1 partner', color='#f39c12')
        ax_c.bar(x_pos + width, self.summary_df['pct_multiple_partners'], width, 
                label='2+ partners', color='#e74c3c')
        
        ax_c.set_xlabel('Age cohort', fontweight='bold')
        ax_c.set_ylabel('Percentage', fontweight='bold')
        ax_c.set_title('C) Sexual Partner History', loc='left', fontweight='bold')
        ax_c.set_xticks(x_pos)
        ax_c.set_xticklabels(cohorts, rotation=45, ha='right')
        ax_c.legend()
        ax_c.grid(alpha=0.3, axis='y')
        
        # Panel D: Vaccination Eligibility Summary
        ax_d = axes[1, 1]
        
        # Define "ideal" vaccination candidates: never infected, no CIN, no cancer
        pct_ideal = self.summary_df['pct_never_infected'] * \
                    (100 - self.summary_df['pct_cin']) / 100 * \
                    (100 - self.summary_df['pct_cancerous']) / 100
        
        # May still benefit: previously/currently infected but no high-grade lesions
        pct_may_benefit = 100 - pct_ideal - self.summary_df['pct_cin'] - \
                         self.summary_df['pct_cancerous']
        
        # Unlikely to benefit: CIN or cancer
        pct_unlikely = self.summary_df['pct_cin'] + self.summary_df['pct_cancerous']
        
        ax_d.bar(x_pos, pct_ideal, label='Ideal candidates\n(never infected, no disease)', 
                color='#2ecc71')
        ax_d.bar(x_pos, pct_may_benefit, bottom=pct_ideal, 
                label='May benefit\n(infected but no CIN/cancer)', color='#f39c12')
        ax_d.bar(x_pos, pct_unlikely, bottom=pct_ideal + pct_may_benefit,
                label='Unlikely to benefit\n(CIN or cancer)', color='#e74c3c')
        
        ax_d.set_xlabel('Age cohort', fontweight='bold')
        ax_d.set_ylabel('Percentage', fontweight='bold')
        ax_d.set_title('D) Vaccination Benefit Potential', loc='left', fontweight='bold')
        ax_d.set_xticks(x_pos)
        ax_d.set_xticklabels(cohorts, rotation=45, ha='right')
        ax_d.legend(fontsize=8)
        ax_d.grid(alpha=0.3, axis='y')
        
        plt.tight_layout()
        
        return fig
    
    def to_csv(self, filename='cohort_status_summary.csv'):
        """Export summary to CSV"""
        if hasattr(self, 'summary_df'):
            self.summary_df.to_csv(filename, index=False)
            print(f"Cohort status summary saved to {filename}")



class person_years(hpv.Analyzer):
    """
    Analyzer to calculate person-years at risk for cervical cancer
    over the simulation period.
    
    Args:
        start_year: First year to start counting (default: 2025)
        end_year: Last year to count (default: 2100)
        min_age: Minimum age for at-risk population (default: 15)
        max_age: Maximum age for at-risk population (default: 85)
        sex: 'f' for female (default), 'm' for male, or 'both'
    """
    
    def __init__(self, start_year=2025, end_year=2100, min_age=15, max_age=85, sex='f', **kwargs):
        super().__init__(**kwargs)
        self.start_year = start_year
        self.end_year = end_year
        self.min_age = min_age
        self.max_age = max_age
        self.sex = sex
        self.person_years = 0
        self.annual_pys = {}  # Store by year for detailed tracking
        
    def initialize(self, sim):
        super().initialize(sim)
        self.person_years = 0
        self.annual_pys = {year: 0 for year in range(self.start_year, self.end_year + 1)}
        self.initialized = True
        
    def apply(self, sim):
        """
        Calculate person-years at each timestep
        """
        current_year = sim.yearvec[sim.t]
        
        # Only count if within the specified time period
        if current_year >= self.start_year and current_year <= self.end_year:
            # Get people who are alive and in the age range
            alive = sim.people.alive
            age = sim.people.age
            
            # Filter by sex if specified
            if self.sex == 'f':
                sex_filter = sim.people.is_female
            else:
                sex_filter = np.ones(len(sim.people), dtype=bool)
            
            # Identify at-risk population
            at_risk = alive & sex_filter & (age >= self.min_age) & (age <= self.max_age)
            
            # Count person-years for this timestep
            # Multiply by dt (timestep length in years)
            pys_this_step = np.sum(at_risk) * sim['dt']
            
            self.person_years += pys_this_step
            
            # Store annual data (accumulate within each year)
            year_key = int(np.floor(current_year))
            if year_key in self.annual_pys:
                self.annual_pys[year_key] += pys_this_step
    
    def finalize(self, sim):
        """
        Calculate final summary statistics
        """
        super().finalize(sim)
        self.total_person_years = self.person_years
        self.average_annual_population = self.person_years / (self.end_year - self.start_year + 1)
        
    def plot(self):
        """
        Plot person-years at risk over time
        """
        import matplotlib.pyplot as plt
        
        years = sorted(self.annual_pys.keys())
        pys = [self.annual_pys[y] for y in years]
        
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.plot(years, pys, linewidth=2)
        ax.set_xlabel('Year')
        ax.set_ylabel('Person-years at risk')
        ax.set_title(f'Person-years at risk (ages {self.min_age}-{self.max_age})')
        ax.grid(alpha=0.3)
        
        # Add total as text
        ax.text(0.05, 0.95, f'Total: {self.total_person_years:,.0f} person-years', 
                transform=ax.transAxes, verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        
        return fig
        
    def to_json(self):
        """
        Export results for use in other analyses
        """
        return {
            'total_person_years': self.total_person_years,
            'average_annual_population': self.average_annual_population,
            'annual_person_years': self.annual_pys,
            'parameters': {
                'start_year': self.start_year,
                'end_year': self.end_year,
                'min_age': self.min_age,
                'max_age': self.max_age,
                'sex': self.sex
            }
        }


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
