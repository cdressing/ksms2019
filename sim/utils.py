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


