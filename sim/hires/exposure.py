try:                # Python2
    from StringIO import StringIO
except ImportError: # Python3
    from io import StringIO

import astropy
import pandas as pd
import numpy as np
iodfactor = 0.7
from scipy.optimize import brentq

def exposure_time(vmag, expcounts, iod=False, t1=110.0, v1=8.0, exp1=250.0):
    """Exposure Time

    Estimate exposure time based on scaling. Cannonical exposure time
    is 110s to get to 250k on 8th mag star with iodine cell in.
    
    Args:
        expcounts (float) : desired number of counts. 
            250 = 250k, 10 = 10k (CKS) i.e. SNR = 45 per pixel.
        iod (bool) : is iodine cell in or out? If out, throughput is higher 
            by 30%
    
    Returns:
        exptime (float) : exposure time [seconds]
    
    """

    # flux star / flux 8th mag star
    fluxfactor = 10.0**(-0.4*(vmag-v1)) 
    exptime = t1 / fluxfactor 
    exptime *= expcounts / exp1
    if iod==False:
        exptime *= iodfactor
    return exptime

def exposure_counts(vmag, exptime, **kwargs):
    """Exposure counts

    Inverse of `exposure_time.` Given a magnitude and an exposure
    time, how many counts will be collected?

    Args:
        vmag (float) : vband magnitude
        exptime (float) : exposure time in seconds)
        **kwargs : keyword arguments passed to exposure_time
    
    Returns:
        _counts (float) : expected number of counts (seconds)
    """
    f = lambda expcounts : exposure_time(vmag, expcounts, **kwargs) - exptime
    _counts = brentq(f,0,2000,)
    return _counts

def exp_ramp(vmag, v1=10.5, v2=12.0, c1=250, c2=60):
    """Calculate expcounts

    Calculate expcounts based on min, and max limits with linear
    ramp between the two mag limits

    Args:
        vmag (float): target magnitude
        v1 (float): below this mag targets get full counts (c1)
        v2 (float): fainter than this mag targets get c2
        c1 (float): expcounts (k) for bright targets
        c2 (float): expcounts (k) for the faint limit

    Returns:
        float: exposure meter setting (k)
    """

    if vmag <= v1:
        return c1

    if vmag >= v2:
        return c2

    exp_level = np.interp(vmag, xp=[v1, v2], fp=[np.log10(c1), np.log10(c2)])

    return 10**exp_level

