import pandas as pd
from sim import mock_up_rvs
import radvel
import numpy as np

import sim.hires.exposure
# Global survey setup
hours_in_night = 8

def counts_to_err(counts):
    '''Compute the expected RV error for an iodine-in observation, scaling from 2.5 m/s at 250k counts'''
    return 2 / np.sqrt(counts/60)

def load_candidates():
    simulated_planets = '../data/barclay+CTL.csv'
    df = pd.read_csv(simulated_planets, comment="#")
   
    units = dict(df.iloc[0])
    df = df.iloc[2:]
   
    for col in df.columns:
        try:
            df[col] = df[col].values.astype(float)
        except ValueError:
            continue
   
    # df = df.merge(ctl, on='TICID', how='left', suffixes=['', '_tic'])
    idf = df.copy()
    
    ## Perform some calculations with the columns

    multis = df.groupby('TICID').count().sort_values(by='Vmag', ascending=False).query('Vmag > 1').reset_index()
    multis['Npl'] = multis['Vmag']
    multis = multis[['TICID', 'Npl']]
    df = pd.merge(df, multis, on='TICID', how='left')

    df['Rs'] = df['R*'].values
    df['Mp'] = mock_up_rvs.calc_mp(df['Rp'].values)
    df['Kp'] = radvel.utils.semi_amplitude(df['Mp'].values, 
                                           df['Perp'].values, 
                                           df['R*'].values,   # R* for M*
                                           df['Ecc'].values,
                                           Msini_units='earth')
    df['a'] = (df['mass'] * (df['Perp'] / 365.) ** 2.) ** (1. / 3.)
    df['sinc'] = (df['Teff'] / 5778.) ** 4.0 * (df['R*'] / df['a']) ** 2.0
    return df


def add_exptime(df, method):
    """
    nobs = 60  # number of observations per target
    """
    if method=='nobs=60-counts=ramp':
        kexp = 250  # maximum exp counts
        nobs = 60
        kexp = np.array([sim.hires.exposure.exp_ramp(v) for v  in df['Vmag']])
        exp = sim.hires.exposure.exposure_time(df['Vmag'], kexp, iod=True)
        exp = np.clip(exp, 0, 20*60) # clip at max time
        kexp = [sim.hires.exposure.exposure_counts(v, e, iod=True) for v, e in zip(df['Vmag'], exp)]
        kexp = np.array(kexp)
        
    elif method=='nobs=30-counts=60':
        kexp = 60  # maximum exp counts
        nobs = 30
        exp = sim.hires.exposure.exposure_time(df['Vmag'], kexp, iod=True)
        exp = np.clip(exp, 0, 15*60) # clip at max time
        kexp = [sim.hires.exposure.exposure_counts(v, e, iod=True) for v, e in zip(df['Vmag'], exp)]
        kexp = np.array(kexp)
    else:
        assert False, "Invalid method"
    df['kexp'] = kexp
    df['exptime'] =  exp * nobs # iodine in exposure
    df['exptime'] += sim.hires.exposure.exposure_time(df['Vmag'], 3*kexp, iod=False) # template exposure
    df['sigmajit'] = np.random.uniform(1.0, 4.0, size=len(df))
    df['sigmaphot'] = counts_to_err(df['kexp'])
    df['sigmarv'] = df.eval('sqrt(sigmajit**2 + sigmaphot**2)')
    df['Kerr'] = df['sigmarv'] / np.sqrt(nobs)
    df['Ksig'] = df.eval('Kp / Kerr')
    df = df.sort_values(by='Vmag').reset_index()
    return df

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
    description = '10 closest dwarfs that are brighter than V < 13'
    def in_sample(self):
        df = self.df.copy()
        b = pd.Series(False, index=self.df.index) 
        idx = self.df.query('_DEJ2000 > 0 and Vmag < 13 and Rs < 0.7').sort_values(by='Dist').iloc[:8].index
        b.loc[idx] = True
        return b
    
class SampleCool(Sample):
    name = 'cooldwarfs'
    description = 'Coolest small planets orbiting stars brighter than V=15'
    def in_sample(self):
        df = self.df.copy()
        b = pd.Series(False, index=self.df.index) 
        #Insolation < 1.52 F_Earth corresponds to the recent Venus limit
        idx = self.df.query('Insol < 10 and Rp < 2.5 and Vmag < 15').sort_values(by='Vmag').iloc[:8].index
        b.loc[idx] = True
        return b
    
class SampleHabitable(Sample):
    name = 'habsystems'
    description = 'All systems with potentially habitable planets and host stars brighter than V = 15'
    def in_sample(self):
        df = self.df.copy()
        b = pd.Series(False, index=self.df.index) 
        #Insolation < 1.52 F_Earth corresponds to the recent Venus limit
        idx = self.df.query('Insol < 1.52 and Vmag < 15').sort_values(by='Insol').iloc[:99].index
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

class SampleBrightShallow(Sample):
    name = 'brightshallow'
    description = 'Brightest 200 sun-like stars'
    def in_sample(self):
        b = pd.Series(False, index=self.df.index) 
        filters = "Teff < 6100 and Vmag < 13.0 and Kp > 2.0 and Rs < 1.5 and  Ksig > 4.0"
        idx = self.df.query(filters).sort_values(by='Vmag').iloc[:200].index
        b.loc[idx] = True
        return b

TKS_Samples = [SampleBright, SampleClose, SampleMulti, SampleUSP, SampleCool, SampleHabitable]
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

