"""
Define an HPVsim simulation for Kenya, including calibration
"""

# Standard imports
import numpy as np
import sciris as sc
import hpvsim as hpv

# Imports from this repository
import utils as ut
import plot_age_causal as pac

# %% Settings and filepaths

# Debug switch
debug = 0  # Run with smaller population sizes and in serial 
do_shrink = True  # Do not keep people when running sims (saves memory)

# Run settings
load_partial = False
n_trials    = [1600, 2][debug]  # How many trials to run for calibration
n_workers   = [40, 1][debug]    # How many cores to use
storage = None

# Save settings
do_save = True
save_plots = True


# %% Simulation creation functions
def make_sim(location='kenya', calib_pars=None, debug=0, interventions=None, analyzers=None, seed=1, end=2020, verbose=0.1):
    """
    Define parameters, analyzers, and interventions for the simulation
    """

    # Basic parameters
    pars = sc.objdict(
        n_agents=[10e3, 1e3][debug],
        dt=[0.25, 1.0][debug],
        start=[1960, 1980][debug],
        end=end,
        genotypes=[16, 18, 'hi5', 'ohr'],
        location=location,
        ms_agent_ratio=100,
        verbose=verbose,
        rand_seed=seed,
    )

    # Sexual behavior parameters
    # Debut: derived by fitting to 2014 DHS
    # Women:

    # Men:
    #           Age:   15,   18,   20,   22,   25
    #   Prop_active: 21.3	56.1	75.8	87.2	92.8
    # For fitting, see https://www.researchsquare.com/article/rs-3074559/v1
    pars.debut = dict(
        f=dict(dist='lognormal', par1=16, par2=4),
        m=dict(dist='lognormal', par1=18, par2=4),
        # f=dict(dist='lognormal', par1=18.28, par2=3.25),
        # m=dict(dist='lognormal', par1=17.71, par2=3.33),
    )

    # Participation in marital and casual relationships
    # Derived to fit 2014 DHS data
    # For fitting, see https://www.researchsquare.com/article/rs-3074559/v1
    # TOO OLD
    pars.layer_probs = dict(
        m=np.array([
            # Share of people of each age who are married
            [0, 5, 10,    15,     20,     25,     30,     35,     40,     45,   50,   55,   60,   65,    70,   75],
            # [0, 0,  0,  0.1596, 0.4466, 0.5845, 0.6139, 0.6202, 0.6139, 0.5726, 0.35, 0.21, 0.14, 0.07, 0.035, 0.007],
            # [0, 0,  0,  0.1,     0.1,    0.15,    0.15,    0.15,   0.2,    0.3,  0.4,  0.4,  0.2, 0.07, 0.035, 0.007],
            # [0, 0,  0,  0.1,     0.1,    0.15,    0.15,    0.2,    0.2,    0.4,  0.4,  0.4,  0.2,  0.1,  0.05, 0.01 ],
            [0, 0,  0,  0.1,     0.3,    0.75,    0.6,    0.5,   0.5,    0.55,  0.5,  0.4,  0.1, 0.07, 0.035, 0.007],
            [0, 0,  0,  0.1,     0.2,    0.5,    0.75,    0.75,   0.6,    0.6,  0.5,  0.4,  0.1,  0.1,  0.05, 0.01 ],
        ]),
        c=np.array([
            # Share of people of each age in casual partnerships
            [0, 5,  10,  15,  20,  25,  30,   35,   40,   45,  50,  55,   60,   65,   70,   75],
            [0,  0, 0.2, 0.2, 0.2, 0.2, 0.2, 0.3,  0.3,  0.3, 0.2, 0.2, 0.10, 0.02, 0.02, 0.02],
            [0,  0, 0.2, 0.3, 0.4, 0.4, 0.4, 0.4,  0.5,  0.4, 0.4, 0.2, 0.02, 0.02, 0.02, 0.02]
        ])

        # # TOO YOUNG
        # m=np.array([
        #     # Share of people of each age who are married
        #     [0, 5, 10,    15,     20,     25,     30,     35,     40,     45,   50,   55,   60,   65,    70,   75],
        #     # [0, 0,  0,  0.1596, 0.4466, 0.5845, 0.6139, 0.6202, 0.6139, 0.5726, 0.35, 0.21, 0.14, 0.07, 0.035, 0.007],
        #     # [0, 0,  0,  0.1,     0.1,    0.15,    0.15,    0.15,   0.2,    0.3,  0.4,  0.4,  0.2, 0.07, 0.035, 0.007],
        #     # [0, 0,  0,  0.1,     0.1,    0.15,    0.15,    0.2,    0.2,    0.4,  0.4,  0.4,  0.2,  0.1,  0.05, 0.01 ],
        #     [0, 0,  0,  0.1,     0.3,    0.5,    0.5,    0.5,   0.6,    0.7,  0.5,  0.4,  0.1, 0.07, 0.035, 0.007],
        #     [0, 0,  0,  0.1,     0.2,    0.5,    0.5,    0.5,   0.6,    0.7,  0.5,  0.4,  0.1,  0.1,  0.05, 0.01 ],
        # ]),
        # c=np.array([
        #     # Share of people of each age in casual partnerships
        #     [0, 5,  10,  15,  20,  25,  30,   35,   40,   45,  50,  55,   60,   65,   70,   75],
        #     [0,  0, 0.2, 0.4, 0.4, 0.4, 0.4, 0.4,  0.4,  0.4, 0.3, 0.2, 0.10, 0.02, 0.02, 0.02],
        #     [0,  0, 0.2, 0.4, 0.4, 0.4, 0.4, 0.4,  0.4,  0.4, 0.5, 0.2, 0.02, 0.02, 0.02, 0.02]
        # ])

        # m=np.array([
        #     # Share of people of each age who are married
        #     [0, 5,  10,    15,   20,   25,   30,   35,   40,   45,   50,   55,   60,   65,   70,   75],
        #     [0, 0, 0.05, 0.25, 0.70, 0.90, 0.95, 0.70, 0.75, 0.65, 0.55, 0.40, 0.40, 0.40, 0.40, 0.40],  # Females
        #     [0, 0, 0.01, 0.01, 0.10, 0.50, 0.60, 0.70, 0.70, 0.70, 0.70, 0.80, 0.70, 0.60, 0.50, 0.60]]  # Males
        # ),
        # c=np.array([
        #     # Share of people of each age in casual partnerships
        #     [0, 5,  10,  15,  20,  25,  30,   35,   40,   45,  50,  55,   60,   65,   70,   75],
        #     [0,  0, 0.1, 0.2, 0.2, 0.2, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.1, 0.01, 0.01, 0.01],
        #     [0,  0, 0.1, 0.2, 0.2, 0.2, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.1, 0.01, 0.01, 0.01],
        # ])
    )
    pars.layer_probs['m'][1] *= .7

    pars.m_partners = dict(
        m=dict(dist='poisson1', par1=0.01),
        c=dict(dist='poisson1', par1=0.2),
    )
    pars.f_partners = dict(
        m=dict(dist='poisson1', par1=0.01),
        c=dict(dist='poisson1', par1=0.2),
    )

    if calib_pars is not None:
        pars = sc.mergedicts(pars, calib_pars)

    if analyzers is None:
        analyzers = []

    sim = hpv.Sim(pars=pars, interventions=interventions, analyzers=analyzers)

    return sim


# %% Simulation running functions
def run_sim(calib_pars=None, analyzers=None, debug=debug, seed=1, verbose=.1, do_shrink=do_shrink, do_save=do_save, end=2020):
    # Make sim
    sim = make_sim(
        debug=debug,
        seed=seed,
        analyzers=analyzers,
        calib_pars=calib_pars,
        end=end
    )
    sim.label = f'Sim-{seed}'

    # Run
    sim['verbose'] = verbose
    sim.run()
    if do_shrink:
        sim.shrink()

    # Optionally save
    if do_save:
        sim.save(f'raw_results/kenya.sim')

    return sim


def make_calib(n_trials=None, n_workers=None):
    sim = make_sim(verbose=-1)
    datafiles = [
        'data/kenya_cancer_cases.csv',
        'data/kenya_cancer_types.csv',
        'data/kenya_asr_cancer_incidence.csv',
        'data/kenya_precin_prevalence.csv'
    ]

    # Define the calibration parameters
    genotype_pars = dict(
        hi5=dict(
            cancer_fn=dict(transform_prob=[1.5e-3, 0.5e-3, 2.5e-3, 2e-4]),
            cin_fn=dict(k=[.15, .1, .25, 0.01]),
            dur_cin=dict(par1=[4.5, 3.5, 5.5, 0.5], par2=[20, 16, 24, 0.5]),
        ),
        ohr=dict(
            cancer_fn=dict(transform_prob=[1.5e-3, 0.5e-3, 2.5e-3, 2e-4]),
            cin_fn=dict(k=[.15, .1, .25, 0.01]),
            dur_cin=dict(par1=[4.5, 3.5, 5.5, 0.5], par2=[20, 16, 24, 0.5]),
        ),
    )

    calib_pars = dict(
        beta=[0.2, 0.02, 0.5, 0.02],
        own_imm_hr=[0.5, 0.5, 1, 0.05],
        cell_imm_init=dict(par1=[0.5, 0.5, 0.8, 0.05]),
        age_risk=dict(risk=[1, 1, 4, 0.1], age=[30, 30, 45, 1]),
        # m_cross_layer=[0.3, 0.1, 0.7, 0.05],
        # m_partners=dict(
        #     c=dict(par1=[0.2, 0.1, 0.6, 0.02])
        # ),
        # f_cross_layer=[0.1, 0.05, 0.5, 0.05],
        # f_partners=dict(
        #     c=dict(par1=[0.2, 0.1, 0.6, 0.02])
        # ),
        sev_dist=dict(par1=[1, 0.5, 1.5, 0.01])
    )

    calib = hpv.Calibration(sim, calib_pars=calib_pars, genotype_pars=genotype_pars,
                            name=f'kenya_calib',
                            datafiles=datafiles,
                            total_trials=n_trials, n_workers=n_workers,
                            storage=storage
                            )
    return sim, calib


def run_calib(n_trials=None, n_workers=None, do_save=True, filestem='', load_partial=True):

    # Run calibration
    sim, calib = make_calib(n_trials=n_trials, n_workers=n_workers)

    if load_partial:
        # Load a partially-run calibration study
        import optuna as op
        print(calib.run_args.name)
        study = op.load_study(storage=calib.run_args.storage, study_name=calib.run_args.name)
        # calib.run_args.continue_db = True
        # calib.calibrate()
        output = study.optimize(calib.run_trial, n_trials=19)
        calib.best_pars = sc.objdict(study.best_params)
        calib.parse_study(study)
        print('Best pars:', calib.best_pars)

        # Tidy up
        calib.calibrated = True
        if not calib.run_args.keep_db:
            calib.remove_db()

    else:
        calib.calibrate()

    print(f'... finished calibration')
    print(f'Best pars are {calib.best_pars}')

    filename = f'kenya_calib{filestem}'

    if do_save:
        sc.saveobj(f'raw_results/{filename}.obj', calib)
        calib_pars = calib.trial_pars_to_sim_pars(which_pars=0)
        sc.save(f'results/kenya_pars{filestem}.obj', calib_pars)

    return sim, calib


def plot_calib(which_pars=0, save_pars=True, filestem=''):
    filename = f'kenya_calib{filestem}'
    calib = sc.load(f'raw_results/{filename}.obj')

    sc.fonts(add=sc.thisdir(aspath=True) / 'Libertinus Sans')
    sc.options(font='Libertinus Sans')
    fig = calib.plot(res_to_plot=200, plot_type='sns.boxplot', do_save=False, do_show=False)
    fig.tight_layout()
    fig.savefig(f'figures/{filename}.png')

    if save_pars:
        calib_pars = calib.trial_pars_to_sim_pars(which_pars=which_pars)
        trial_pars = sc.autolist()
        for i in range(100):
            trial_pars += calib.trial_pars_to_sim_pars(which_pars=i)
        sc.save(f'results/kenya_pars{filestem}.obj', calib_pars)
        sc.save(f'results/kenya_pars{filestem}_all.obj', trial_pars)

    return calib


def run_parsets(debug=False, verbose=.1, analyzers=None, save_results=True, **kwargs):
    ''' Run multiple simulations in parallel '''

    parsets = sc.loadobj(f'results/kenya_pars_all.obj')
    kwargs = sc.mergedicts(dict(debug=debug, end=2040, verbose=verbose, analyzers=analyzers), kwargs)
    simlist = sc.parallelize(run_sim, iterkwargs=dict(calib_pars=parsets), kwargs=kwargs, serial=debug, die=True)
    msim = hpv.MultiSim(simlist)
    msim.reduce()
    if save_results:
        sc.saveobj(f'raw_results/kenya_msim.obj', msim.results)

    return msim





# %% Run as a script
if __name__ == '__main__':

    # List of what to run
    to_run = [
        # 'run_sim',
        # 'age_pyramids',
        'run_calib',
        # 'plot_calib'
        # 'run_parsets'
        # 'plot_age_causal'
    ]

    location = 'kenya'
    T = sc.timer()  # Start a timer

    if 'run_sim' in to_run:
        calib_pars = sc.loadobj('results/kenya_pars.obj')  # Load parameters from a previous calibration
        analyzers = [hpv.age_causal_infection(start_year=2020)]
        sim = run_sim(calib_pars=calib_pars, do_save=False, do_shrink=False, analyzers=analyzers)
        df = pac.get_age_causal_df(sim)
        sc.saveobj(f'results/age_causal_infection_kenya.obj', df)

    if 'age_pyramids' in to_run:
        # calib_pars = sc.loadobj('results/kenya_pars.obj')
        ap = hpv.age_pyramid(
            timepoints=['2025', '2050', '2075', '2100'],
            datafile=f'data/{location}_age_pyramid.csv',
            edges=np.linspace(0, 100, 21),
        )
        sim = run_sim(end=2100, calib_pars=calib_pars, analyzers=[ap], do_save=True, do_shrink=True)

    if 'run_calib' in to_run:
        sim, calib = run_calib(n_trials=n_trials, n_workers=n_workers, filestem='', do_save=True, load_partial=load_partial)

    if 'plot_calib' in to_run:
        calib = plot_calib(save_pars=True, filestem='')
        calib = ut.shrink_calib(calib, n_results=200)
        sc.saveobj(f'raw_results/kenya_calib_reduced.obj', calib)

    if 'run_parsets' in to_run:
        msim = run_parsets()

    if 'plot_age_causal' in to_run:
        pac.plot_age_causal(location='kenya')

    if 'plot_age_causal_violin' in to_run:
        pac.plot_age_causal_violin(location='kenya')

    T.toc('Done')  # Print out how long the run took
