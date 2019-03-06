import pandas as pd
import numpy as np


def maximal_coverage(df, par, samples):
    count = len(par)
    
    spacing = int(np.round(count/samples))
    
    print(count, spacing)
    return df.sort_values(by='par').iloc[::spacing]


def read_ctl(filename):
    colnames = [line.split()[0][1:-1] for line in open('../data/tic_column_description.txt', 'r').readlines()]
    colnames[0] = 'TICID'
    print(colnames)

    ctl = pd.read_csv(ctl_file, names=colnames)

    print(ctl.columns)
    ctl.head()


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

