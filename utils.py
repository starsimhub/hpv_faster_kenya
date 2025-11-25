"""
Utilities
"""

# Imports
import sciris as sc
import numpy as np
import hpvsim as hpv
from scipy.stats import norm, lognorm
import pandas as pd


def set_font(size=None, font='Libertinus Sans'):
    """ Set a custom font """
    sc.fonts(add=sc.thisdir(aspath=True) / 'assets' / 'LibertinusSans-Regular.otf')
    sc.options(font=font, fontsize=size)
    return


def shrink_calib(calib, n_results=100):
    cal = sc.objdict()
    plot_indices = calib.df.iloc[:n_results, 0].values
    cal.sim_results = [calib.sim_results[i] for i in plot_indices]
    cal.target_data = calib.target_data
    cal.df = calib.df.iloc[0:n_results, ]
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

