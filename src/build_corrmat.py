import pandas as pd
import os
import os.path

__author__ = 'Rasmus Magnusson'
__COPYRIGHT__ = 'Copyright (C) 2021 Rasmus Magnusson'
__contact__ = 'rasma774@gmail.com'

PATH = __file__[:-20]

def check_corrmat():
    if not os.path.isfile('/' + PATH + '/data/pickles/correlations.p'):
        print('building pickle file')
        # Assemble the correlation matrix
        picklepath = '/' + PATH + '/data/buildpickles/'
        pfiles = os.listdir(picklepath)
        dfs = []
        for pickle_file in pfiles:
            dfs.append(pd.read_pickle(picklepath + pickle_file))
        dfs = pd.concat(dfs).sort_index()
        dfs.to_pickle('/' + PATH + '/data/pickles/correlations.p')
        print('done')