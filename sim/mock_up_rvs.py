import pandas as pd
from astropy.time import Time
import numpy as np
import radvel
import sys


def calc_mp(rp):

    from collections.abc import Iterable
    if isinstance(rp, Iterable):
        return np.array([calc_mp(r) for r in rp])

    if rp > 1.5:
        mp = 2.69*pow(rp, 0.93) # Earth radii, Earth masses
    else:
        rhoe = 5.41  # density of Earth, gcc
        rhop = 2.43 + 3.39 * rp  # gcc
        mp = rhop/rhoe * pow(rp, 3.0)  # Earth masses, Earth radii, density ratio

    return mp

def draw_ecc(size=1):
    from numpy import random
    ecc_rayleigh = random.rayleigh(size=size, scale=0.05)
    return np.minimum(ecc_rayleigh, 1.0)  # not larger than 1


def rverr_photon(mag, exptime, verbose=False):
    from cpsutils.hires import exposure
    expcts = exposure.exposure_counts(mag, exptime, iod=True, t1=1500.0, v1=11.0, exp1=200.0)
    rverr = 1.6 * pow(200/expcts,1/2.)
    if verbose==True:
        print("Mag \t Exptime (s) \t kCounts \t RV_err (m/s):")
        print("%0.1f \t \t %i \t %i \t %0.1f" % (mag, exptime, expcts, rverr))
    return rverr


def mass_ecc_mystar(mystar):
    from numpy import sqrt, sin, cos, pi

    mystar['Mp'] = calc_Mp(mystar['iso_prad'].values) # calculate masses using M-R
    mystar['ecc'] = draw_ecc(size=len(mystar))
    w = np.random.uniform(low=0,high=(2*pi),size=len(mystar))
    mystar['secw'] = np.sqrt(mystar['ecc'].values) * cos(w)
    mystar['sesw'] = np.sqrt(mystar['ecc'].values) * sin(w)
    mystar['K'] = radvel.utils.semi_amplitude(mystar.Mp, mystar.koi_period, mystar.iso_smass, mystar.ecc, Msini_units='earth')
    return mystar

def rvmystar(mystar, t, verbose=False):
    mod = initialize_model(mystar, vary_ecc=False)
    try:
        errvel = rverr_photon(mystar.kic_kepmag.values[0], 2700., verbose=verbose)
    except:
        errvel = rverr_photon(mystar.kic_kepmag.values[0], 600., verbose=verbose)
    errvel = pow(errvel**2. + 2.0**2., 0.5) # add 2 m/s jitter
    vel = mod(t) + np.random.normal(scale=errvel, size=len(t))
    
    return(t, vel, errvel)


def plot_results(like, mystar):
    fig = plt.figure(figsize=(12,4))
    fig = plt.gcf()
    fig.set_tight_layout(True)
    plt.errorbar(
        like.x, like.model(t)+like.residuals(),
        yerr=like.yerr, fmt='o', label='Mock RVs'
        )
    plt.plot(ti, like.model(ti), label='Best-fit Model')
    plt.xlabel('Time (days)')
    plt.ylabel('RV (m/s)')
    plt.title(mystar.id_starname.values[0])
    plt.legend()
    plt.draw()
    
def get_K(mystar,post):
    best_K = []
    for i, pl in enumerate(mystar.iterrows()):
        best_K.append(post.likelihood.params['k'+str(i+1)].value)
    return np.array(best_K)


def goodness_of_RVs(mystar, t, vel, err, vary_ecc=False, plot=False, save=False):
    mod = initialize_model(mystar, vary_ecc=vary_ecc)
    like = radvel.likelihood.RVLikelihood(mod, t, vel, errvel)
    like.params['gamma'] = radvel.Parameter(value=0.0,vary=True)
    like.params['jit'] = radvel.Parameter(value=2.0,vary=False)
    
    # Define prior shapes and widths here.
    npl = len(mystar)
    
    kpriors = [radvel.prior.HardBounds('k'+str(i+1), 0.0, 10.0) for i in range(npl)]
    priors = [
        radvel.prior.EccentricityPrior( npl ),           # Keeps eccentricity < 1
        radvel.prior.HardBounds('jit', 0.0, 10.0),
    ]
    post = radvel.posterior.Posterior(like)
    post.priors += priors
    post.priors += kpriors
    if plot==True:
        print(post)


    from scipy import optimize
    res  = optimize.minimize(
    post.neglogprob_array,     # objective function is negative log likelihood
    post.get_vary_params(),    # initial variable parameters
    method='nelder-mead',           # Nelder-Mead also works
    )
    
    if plot==True:
        plot_results(like, mystar)             # plot best fit model
        print(post)
        if save==True:
            plt.savefig("../plots/"+str(mystar.id_starname.values[0])+"_mock_rvs.png",dpi=180)
            text_file = open("../trials/"+str(mystar.id_starname.values[0])+"_post.txt", "w")
            text_file.write("%s" % post.__repr__)
            text_file.close()

    true_K = np.array(mystar.K)
    best_K = np.array(get_K(mystar,post))
    return((best_K - true_K)/true_K)

def monte_carlo_delta_K(mystar,t):
    delta_K = []
    for i in range(100):
        t, vel, errvel = rvmystar(mystar, t) # change the RVs
        delta_K.append(goodness_of_RVs(mystar, t, vel, errvel))
    delta_K = np.array(delta_K)
    return delta_K


def plot_delta_K(mystar, delta_K):
    matplotlib.rcParams['font.size'] = 14
    matplotlib.rcParams['axes.labelsize'] = 14    # fontsize of the x and y labels
    matplotlib.rcParams['axes.titlesize'] = 16    # fontsize of the x and y labels
    fig = figure(figsize=(6,4))
    fig = gcf()
    fig.set_tight_layout(True)

    mean_K_frac = np.mean(delta_K,axis=0)
    sig_K_frac = np.std(delta_K,axis=0)
    for col in range(delta_K.shape[1]):
        hist(delta_K[:,col],alpha=0.5,range=[-1,2],label='$K_%i; K/\sigma_K$=%0.1f' %((col+1),np.round(1./sig_K_frac[col],1)))
    legend()
    ylabel('Number of Trials')
    xlabel('Fractional Change in K')
    title(mystar.id_starname.values[0])
    savefig("../plots/"+str(mystar.id_starname.values[0])+"DeltaK_dist.png", dpi=180)
    

def mock_my_star(id_starname):
    print(id_starname)
    n_per_sum = 60
    day_spread = 200./n_per_sum
    year = 365
    t = np.concatenate([arange(n_per_sum) * day_spread, 1*year + arange(n_per_sum) * day_spread, 2* year + arange(n_per_sum//2) * day_spread])
    mystar = planets[planets.id_starname==id_starname].sort_values('koi_period')
    mystar = mass_ecc_mystar(mystar)
    t, vel, errvel = rvmystar(mystar, t, verbose=True)
    goodness_of_RVs(mystar, t, vel, errvel, vary_ecc=False, plot=True, save=True)
    delta_K = monte_carlo_delta_K(mystar,t)
    plot_delta_K(mystar,delta_K)
    final_cols = list(np.concatenate([cols_of_interest,['Mp','K']]))
    mystar[final_cols].set_index('id_koicand').to_csv('../individual_planet_tables/'+str(mystar.id_starname.values[0])+'_table.csv')
    np.save("../trials/"+mystar.id_starname.values[0]+'_delta_K.npy',delta_K)
    return
