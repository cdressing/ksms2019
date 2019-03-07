import pandas as pd
from sim import mock_up_rvs
import radvel
import numpy as np
from sim.hires import exposure
from sim.starlist import counts_to_err


# Global survey setup
kexp = 250  # maximum exp counts
nobs = 60  # number of observations per target
hours_in_night = 8


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
                           np.random.uniform(1.0, 4.0, size=len(sp))**2)  # random jitter from 1.0 to 4.0 m/s jitter
    sp['Ksig'] = (sp['Kp'] * np.sqrt(nobs)) / sp['jitter']
    sp['Kerr'] = sp['jitter'] / np.sqrt(nobs)

    sp = sp.sort_values(by='Vmag').reset_index()
    return sp


class Sample(object):
    def __init__(self, df):
        self.df = df.copy()  # protect the size/order of input

    def get_sample_statistics(self):
        sample = self.get_sample()
        num_nights = sample.exptime.sum() / (hours_in_night*3600)
        num_targets = sample.TICID.count()
        s = """\
name = {}
description = {}
num_nights = {}
num_targets = {}
        """.format(self.name, self.description, num_nights,num_targets)
        return s

## Define your sample here
## Each sample has a function in_sample() that returns true if star is in the sample
class SampleBright(Sample):
    name = 'bright'
    description = 'Brightest 50 sun-like stars'
    def in_sample(self):
        b = pd.Series(False, index=self.df.index) 
        filters = \
            """_DEJ2000 > 0 and Teff < 6100 and Vmag < 13.0 and Kp > 2.0 and Rs < 1.5 and Rp < 12 and  Ksig > 4.0"""
        idx = self.df.query(filters).sort_values(by='Vmag').iloc[:48].index
        b.loc[idx] = True
        return b

## Perhaps we should target the M-dwarfs specifically here
class SampleClose(Sample):
    name = 'close'
    description = '10 closest dwarfs that are bighter than V < 13'
    def in_sample(self):
        df = self.df.copy()
        b = pd.Series(False, index=self.df.index) 
        idx = self.df.query('_DEJ2000 > 0 and Vmag < 13 and Rs < 0.7').sort_values(by='Dist').iloc[:8].index
        b.loc[idx] = True
        return b

## 
class SampleMulti(Sample):
    name = 'multi'
    description = 'Brightest 10 multis'
    def in_sample(self):
        b = pd.Series(False, index=self.df.index) 
        idx = self.df.query('Npl > 1').sort_values(by='Vmag').iloc[:25].index
        b.loc[idx] = True
        return b

class SampleUSP(Sample):
    name = 'usp'
    description = 'Brightest 5 USPs'
    def in_sample(self):
        b = pd.Series(False, index=self.df.index) 
        idx = self.df.query('Perp < 1.5 and Rp < 2').sort_values(by='Vmag').iloc[:8].index
        b.loc[idx] = True
        return b

TKS_Samples = [SampleBright, SampleClose, SampleMulti, SampleUSP]
def combine_samples(df):
    df['inany'] = False
    innames = []
    names = []
    for Sample in TKS_Samples:
        s = Sample(df) 
        names.append(s.name)
        inname = 'in'+s.name
        df[inname] = s.in_sample()
        df['inany'] = df['inany'] | df[inname]
        innames.append(inname)
    
    df['nprograms'] = df[innames].sum(axis=1)
    for name in names:
        inname = 'in'+ name
        sample = df[df[inname]]
        print('in ' + name)
        print(get_sample_statistics(sample))

        print('only in ' + name)
        sample = df[df[inname] & (df['nprograms']==1)]
        print(get_sample_statistics(sample))

    sample = df[df['inany']]
    print('-'*80)
    print('Total')
    print(get_sample_statistics(sample))

    print('Program overlap')
    print(df.groupby('nprograms').size())
    return df

def get_sample_statistics(sample):
    nplanets = len(sample)
    sample = sample.groupby('TICID').first() # only count stars once
    nnights = sample.exptime.sum() / (hours_in_night*3600)
    nstars = len(sample)
    s = """\
nstars = {:.0f}, nplanets={:.0f}, nnights = {:.1f},
    """.format(nstars,nplanets,nnights)
    return s

