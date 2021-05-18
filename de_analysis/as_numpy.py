import numpy as np
import pandas as pd
# from de_analysis import *
# import de_analysis


def save_as_np(pd_original,
               loc_to:str):
    '''

    :param pd_original: pandas Dataframe
    :param loc_to: string
            location to save inlcuding name
    :return:
    '''
    pd_np = pd_original.to_numpy()
    np.save(loc_to + '.npy', pd_np)
    col = pd_original.columns.to_numpy()
    indx = pd_original.index.to_numpy()
    np.save(loc_to + '_col.npy', col)
    np.save(loc_to + '_indx.npy', indx)

    return


def read_numpy_to_df(filename:str, read_only_col = False):
    # read_only_col: read only the columns, not the full matrix
    if read_only_col:
        npcol = np.load(filename + '_col.npy', allow_pickle=True)
        pd_from_np = pd.DataFrame(columns = npcol)
    else:
        npload = np.load(filename + '.npy')
        npcol = np.load(filename + '_col.npy', allow_pickle=True)
        npindx = np.load(filename + '_indx.npy', allow_pickle=True)
        pd_from_np = pd.DataFrame(data=npload, columns=npcol, index=npindx)

    return pd_from_np

# filename = "/home/erika/PycharmProjects/Kevin_GeneSetEnrichmentAnalysis/" \
#            "benchmarking_edgeR_vs_our/data/allCells_Filtered_025genes_0_Pat111"
# read_numpy_to_df(filename, read_only_col=True)
#
# pan = pd.read_csv("/home/erika/PycharmProjects/DE_analysis_clean/de_analysis_clean/"
#            "docs/example/data/myDataset_1vs23_Pt1_23.tsv",index_col=0,sep='/t')
# loc_to = "/home/erika/PycharmProjects/DE_analysis_clean/de_analysis_clean/" \
#          "docs/example/data/myDataset_1vs23_Pt1_23"
# save_as_np(pan,loc_to)

