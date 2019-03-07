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


class Sample(object):
    def __init__(self, df):
        self.df = df 

    def apply_global_cuts(self):
        df = self.df.copy()
        filters = """_DEJ2000 > 0"""
        return df.query(filters)

    def get_sample_statistics(self):
        sample = self.get_sample()
        num_nights = sample.exptime.sum() / (10*3600)
        num_targets = sample.TICID.count()
        s = """\
name = {}
description = {}
num_nights = {}
num_targets = {}
        """.format(self.name, self.description, num_nights,num_targets)
        return s

class SampleBright(Sample):
    name = 'bright'
    description = 'Brightest 50 sun-like stars'
    def get_sample(self):
        df = self.apply_global_cuts()
        filters = """Teff < 6100 and Vmag < 13.0 and Kp > 2.0 and Rs < 1.5 and Rp < 12 and  Ksig > 4.0"""
        df = df.query(filters).sort_values(by='Vmag').iloc[:50]
        return df 


## Perhaps we want to turn this sample into an M-dwarf selection
class SampleClose(Sample):
    name = 'close'
    description = '10 closest dwarfs that are bighter than V < 13'
    def get_sample(self):
        df = self.apply_global_cuts()
        df = df.query('Vmag < 13 and Rs < 1.5').sort_values(by='Dist').iloc[:10]
        return df 

class SampleMulti(Sample):
    name = 'multi'
    description = 'Brightest 10 multis'
    def get_sample(self):
        df = self.apply_global_cuts()
        df = df.query('Npl > 1').sort_values(by='Vmag').iloc[:10]
        return df 

class SampleUSP(Sample):
    name = 'usp'
    description = 'Brightest 5 USPs'
    def get_sample(self):
        df = self.apply_global_cuts()
        df = df.query('Perp < 1.5 and Rp < 2').sort_values(by='Vmag').iloc[:5]
        return df 

TKS_Samples = [SampleBright, SampleClose, SampleMulti, SampleUSP]

class SampleTKS(Sample):
    name = 'tks'
    description = 'Total sample'
    def get_sample(self,verbose=False):
        tks = []
        for Sample in TKS_Samples:
            s = Sample(self.df) 
            if verbose:
                print('adding')
                print(s.get_sample_statistics())
            tks.append(s.get_sample())

        tks = pd.concat(tks)
        tks = tks.drop_duplicates()
        return tks
    
