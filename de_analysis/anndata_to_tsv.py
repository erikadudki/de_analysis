import numpy as np
import pandas as pd
import anndata
import os
import warnings

wd = '/home/erika/PycharmProjects/DE-Analysis/src/code_tidy/'
filename = 'de10;1'
user_layer = 'logcounts'


##

def anndata_to_tsv(wd,
                   filename,
                   user_layer=None):
    """
    Transforms anndata (*.h5ad-file) to individual .tsv files for each patient
    for each cluster, which is the format needed for running the DE-Analysis.

    Input:
        wd: string
            working directory
        filename: string
            name of the data file (h5ad format)
        user_layer: string
            pick which layer of normalized data should be used, usually:
            'logcounts' / 'cpm'

    Returns:
        saves tsv-files with cell-matrices for each patient and each cluster
        into directory: ./data/
    """

    # TODO: add option for sparsematrix?
    adata = anndata.read_h5ad(wd + 'data/' + filename + '.h5ad')
    # which clusters/samples exist:
    cl_id = adata.obs['cluster_id'].unique()
    s_id = adata.obs['sample_id'].unique()

    # for cl in range(1,nr_cl+1):
    for cl in cl_id:
        # subsets of celltypes
        adata_cl = adata[adata.obs['cluster_id'] == cl]

        # subset of samples/patients
        # for s in range(1,nr_s+1):
        for s in s_id:
            adata_cl_s = adata_cl[adata_cl.obs['sample_id'] == s]

            # transform to pandas-dataframe
            # if sparsematrix:
            #     pd_adata_cl_s = pd.DataFrame(data=np.transpose(adata_cl_s.X.toarray()),
            #                                   index=adata_cl_s.var_names,
            #                                   columns=adata_cl_s.obs_names)
            #     if not os.path.exists(wd + 'data/' + 'data_per_pat_per_cl'):
            #         os.mkdir(wd + 'data/' + 'data_per_pat_per_cl')
            #
            #     pd_adata_cl_s.to_csv(
            #         wd + 'data/' + 'data_per_pat_per_cl' + '/cl' + str(cl) + 'Pt' + str(s))
            # else:
            if user_layer is None:
                warnings.warn(
                    "Type of normalized data (user_layer) is not specified. "
                    "Please make sure, that your input data is normalized. "
                    "Now, the main data matrix is considered.")
                pd_adata_cl_s = pd.DataFrame(data=np.transpose(adata_cl_s.X),
                                             columns= adata_cl_s.obs_names,
                                             index=adata_cl_s.var_names)
            else:
                pd_adata_cl_s = pd.DataFrame(
                    data=np.transpose(adata_cl_s.layers[user_layer]),
                    columns=adata_cl_s.obs_names,
                    index=adata_cl_s.var_names)

            # save the panda files for each patient for each cluster
            name_to_save = 'cl' + str(cl) + '_Pt' + str(s)
            pd_adata_cl_s.to_csv(wd + 'data/' + filename
                                 + '_' + name_to_save + '.tsv', sep = '\t')
    return
