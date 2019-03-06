import pandas as pd
from sim import mock_up_rvs
import radvel
import numpy as np
from sim.hires import exposure
from sim.starlist import counts_to_err

kexp = 250  # maximum exp counts
nobs = 60  # number of observations per target

def load_candidates():
    simulated_planets = '../data/barclay+CTL.csv'
    sp = pd.read_csv(simulated_planets, comment="#")
   
    units = dict(sp.iloc[0])
    sp = sp.iloc[2:]
   
    for col in sp.columns:
        try:
            sp[col] = sp[col].values.astype(float)
        except ValueError:
            continue
   
    # sp = sp.merge(ctl, on='TICID', how='left', suffixes=['', '_tic'])
    isp = sp.copy()
    
    ## Perform some calculations with the columns

    multis = sp.groupby('TICID').count().sort_values(by='Vmag', ascending=False).query('Vmag > 1').reset_index()
    multis['Npl'] = multis['Vmag']
    multis = multis[['TICID', 'Npl']]
    sp = pd.merge(sp, multis, on='TICID', how='left')

    sp['Rs'] = sp['R*'].values
    sp['Mp'] = mock_up_rvs.calc_mp(sp['Rp'].values)
    sp['Kp'] = radvel.utils.semi_amplitude(sp['Mp'].values, 
                                           sp['Perp'].values, 
                                           sp['R*'].values,   # R* for M*
                                           sp['Ecc'].values,
                                           Msini_units='earth')
    sp['a'] = (sp['mass'] * (sp['Perp'] / 365.) ** 2.) ** (1. / 3.)
    sp['sinc'] = (sp['Teff'] / 5778.) ** 4.0 * (sp['R*'] / sp['a']) ** 2.0


    exp = np.clip(exposure.exposure_time(sp['Vmag'].values, kexp, iod=True), 0, 2700) 

    kexp_list = [exposure.exposure_counts(v, e, iod=True) for v, e in zip(sp['Vmag'].values, exp)]

    # for v, e in zip(sp['Vmag'].values, exp):
    #     print(v, e, exposure.exposure_counts(v, e, iod=True))
    sp['kexp'] = kexp_list

    sp['exptime'] =  exp * nobs  
    sp['exptime'] += exposure.exposure_time(sp['Vmag'].values, 3*kexp, iod=False)  # template 3*normal obs

    sp['jitter'] = np.sqrt(counts_to_err(sp['kexp'].values)**2 + \
                           np.random.uniform(1.0, 4.0, size=len(sp))**2)  # random jitter from 1.0 to 7.0 m/s jitter
    sp['Ksig'] = (sp['Kp'] * np.sqrt(nobs)) / sp['jitter']
    sp['Kerr'] = sp['jitter'] / np.sqrt(nobs)

    sp = sp.sort_values(by='Vmag').reset_index()
    return sp


def construct_sample(sp):

    
    filters = """_DEJ2000 > 0
    Teff < 6100
    Vmag < 13.5
    Kp > 2.0
    Rs < 1.5
    Rp < 12
    Ksig > 4.0"""

    for filt in filters.split('\n'):
        print(filt)
        sp = sp.query(filt)

    spf = sp.iloc[:50]  # n brightest
    close = sp.sort_values(by='Dist').iloc[:50]
    spf = pd.concat([spf, close])  # +n closest
    # spf = pd.concat([spf, sp.sort_values(by='Rp').iloc[:27]])  # +n smallest
    multi = sp.query('Npl > 1')
    usp = sp.query('Perp < 1.5')
    spf = pd.concat([spf, multi])
    spf = pd.concat([spf, usp])
    spf = spf.drop_duplicates()

    num_nights = (spf.exptime.sum() / (10*3600))
    num_targets = spf.TICID.count()
    print(spf.TICID.count(), (spf.exptime.sum() / (10*3600)), spf.Ksig.min())
    return spf
