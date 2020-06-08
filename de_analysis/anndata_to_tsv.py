import numpy as np
import pandas as pd
import anndata
import os
import warnings
import scipy as sc


wd = '/home/erika/PycharmProjects/DE-Analysis/src/code_tidy/'
wd = '/home/erika/Documents/Projects/Evaluation_DE_method/sc_simulation_w_muscat/kang/anndata/logcounts/'

# filename = 'Lung_sce'
filename = 'de10_ss4;1'
user_layer = 'logcounts'


##

def anndata_to_tsv(wd,
                   filename,
                   user_layer=None):
    """
    Transforms anndata (*.h5ad-file) to individual .tsv files for each patient
    for each cluster, which is the format needed for running the DE-Analysis.
    Needs user input, during running the code, to set the parameters cl_input
    (what is the columnname of the group of cluster-IDs) and s_input (what is
    the columnname of the group for sample/patient-IDs)

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
    print(adata.obs.columns.values)
    cl_input = input('Please enter: which of the IDs (showing on the previous '
                     'line) belong to the cluster-ID? ')
    # cl_input = 'cluster'
    # cl_input = 'cluster_id'
    cl_id = adata.obs[cl_input].unique()
    s_input = input(
        'Please enter: which of the IDs belongs to the sample-ID/patient-ID? ')
    # s_input = 'group'
    # s_input = 'sample_id'
    s_id = adata.obs[s_input].unique()

    # for cl in range(1,nr_cl+1):
    for cl in cl_id:
        # subsets of celltypes
        adata_cl = adata[adata.obs[cl_input] == cl]

        # subset of samples/patients
        # for s in range(1,nr_s+1):
        for s in s_id:
            adata_cl_s = adata_cl[adata_cl.obs[s_input] == s]

            # transform to pandas-dataframe

            # if no layers specified in anndata-file
            if len(adata_cl_s.layers) == 0:
                try:
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
                            data=np.transpose(adata_cl_s.X),
                            columns=adata_cl_s.obs_names,
                            index=adata_cl_s.var_names)
                except ValueError:
                    warnings.warn("Maybe the dataformat is not understood. ")
                    print(adata_cl_s.obs_names)
                    warnings.warn('should correspond to the single cells, and ')
                    print(adata_cl_s.var_names)
                    warnings.warn('should correspond to the genes. If that is correct, '
                          'then probably the data-format is not '
                          'understood and not implemented yet. '
                          'Please post an issue. ')
                    pass
            # If data is in sparsematrix format csc
            elif sc.sparse.isspmatrix_csc(adata_cl_s.layers[user_layer]) or sc.sparse.isspmatrix_csr(adata_cl_s.layers[user_layer]):
                if user_layer is None:
                    warnings.warn(
                        "Type of normalized data (user_layer) is not specified. "
                        "Please make sure, that your input data is normalized. "
                        "Now, the main data matrix is considered.")
                    pd_adata_cl_s = pd.DataFrame(
                        data=np.transpose(adata_cl_s.X.toarray()),
                        columns= adata_cl_s.obs_names,
                        index=adata_cl_s.var_names)
                else:
                    pd_adata_cl_s = pd.DataFrame(
                        data=np.transpose(adata_cl_s.layers[user_layer].toarray()),
                        index=adata_cl_s.var_names,
                        columns=adata_cl_s.obs_names)
                    if not os.path.exists(wd + 'data/' ):
                        os.mkdir(wd + 'data/')

            # elif type(adata_cl_s.layers[user_layer]) == anndata._core.views.ArrayView:
            #     if user_layer is None:
            #         warnings.warn(
            #             "Type of normalized data (user_layer) is not specified. "
            #             "Please make sure, that your input data is normalized. "
            #             "Now, the main data matrix is considered.")
            #         pd_adata_cl_s = pd.DataFrame(
            #             data=np.transpose(adata_cl_s.X),
            #             columns=adata_cl_s.obs_names,
            #             index=adata_cl_s.var_names)
            #     else:
            #         pd_adata_cl_s = pd.DataFrame(
            #             data=np.transpose(adata_cl_s.layers[user_layer]),
            #             columns=adata_cl_s.obs_names,
            #             index=adata_cl_s.var_names)

            else:
                try:
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
                except ValueError:
                    warnings.warn("Maybe the dataformat is not understood. ")
                    print(adata_cl_s.obs_names)
                    warnings.warn('should correspond to the single cells, and ')
                    print(adata_cl_s.var_names)
                    warnings.warn('should correspond to the genes. If that is correct, '
                          'then probably the data-format is not '
                          'understood and not implemented yet. '
                          'Please post an issue. ')
                    pass

            # save the panda files for each patient for each cluster
            name_to_save = 'cl' + str(cl) + '_Pt' + str(s)
            pd_adata_cl_s.to_csv(wd + 'data/anndata_to_tsv/' + filename
                                 + '_' + name_to_save + '.tsv', sep = '\t')
    return

# anndata_to_tsv(wd,filename,user_layer)