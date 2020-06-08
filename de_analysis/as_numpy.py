import numpy as np
import pandas as pd
# from de_analysis import *
# import de_analysis


def save_as_np(pd_original,
               loc_to):
    '''

    :param pd_original: pandas Dataframe
    :param loc_to: string
            location to save
    :return:
    '''
    pd_np = pd_original.to_numpy()
    np.save(loc_to + '.npy', pd_np)
    col = pd_original.columns.to_numpy()
    indx = pd_original.index.to_numpy()
    np.save(loc_to + '_col.npy', col)
    np.save(loc_to + '_indx.npy', indx)

    return


def read_numpy_to_df(filename):
    npload = np.load(filename + '.npy')
    npcol = np.load(filename + '_col.npy', allow_pickle=True)
    npindx = np.load(filename + '_indx.npy', allow_pickle=True)
    pd_from_np = pd.DataFrame(data=npload, columns=npcol, index=npindx)

    return pd_from_np
