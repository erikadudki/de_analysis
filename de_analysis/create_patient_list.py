# READ DATA ---------------------------------------------------------------
# create a list with np.arrays of all control patients and COPD patients,
# ONE GENE (as input for Wilcoxon rank sum test)
# ###########################################################################
import numpy as np
import pandas as pd
import os
import warnings

run_on_grid = False
if run_on_grid:
    import as_numpy
else:
    from de_analysis import as_numpy
# from de_analysis import *
# import de_analysis

def create_patient_list(patient_control,
                        patient_copd,
                        all_cells_path,
                        ct,
                        ind_filtered_genes,
                        row_gene,
                        read_pd):
    """
    READs DATA ---------------------------------------------------------------
    create a list with np.arrays (with the expression values of the cells
    belonging to the cell type and the patients) of all control patients and
    COPD patients per ONE GENE (as input for Wilcoxon rank sum test)

    Args:
        patient_control: list
            list of patient names belonging to the Control Group
        patient_copd: list
            list of patient names belonging to the Diseased COPD Group
        all_cells_path: string
            path to the matrices all_cells_... (all cells belonging to one
            patient and to one cluster)
        ct: string
            celltype which should be addressed in this analysis run
        ind_filtered_genes: np.array
            indices of genes which should be discarded
        row_gene: int
            index of gene
        read_pd: bool
            True: read tsv files
            False: read npy files

    Returns:
        patient_list: list
            list of np.arrays (->how many: number of patients), each np.array
            belongs to one patient (expression values of the cells for one
            gene)
        patient_list_nonzero: list
            list of np.arrays, same as patient_list, but only showing expressed
            cell values (excluding the zero counts)
        pd_sizes_summary: pandas Dataframe
            gives a summary of number of mean percentage of all cells of
            control patients & COPD, minimal&maximal percentage (how many
            expressed cells) for control/COPD
        len_pat_control: int
            number of patients in control
        len_pat_copd: int
            number of patients in copd
        patient_control, patient_copd: list
            list of patients (updated list) it could happen, that a patient is
            discarded, if no cells are measured for this patient

    """
    if row_gene==0:
        if read_pd:
            print('create patient list - read pandas ---------------------')
        else:
            print('create patient list - read numpy ---------------------')


    patient_list = []
    patient_list_nonzero = []
    pd_sizes = pd.DataFrame([], columns=['total_cells', 'expressed_cells',
                                         'percentage','percentage_copd'],
                            index=patient_control)
    pd_sizes_summary = pd.DataFrame([], columns=['mean_percentage_control',
                                                 'mean_percentage_COPD',
                                                 'min_perc_control',
                                                 'max_perc_control',
                                                 'min_perc_COPD',
                                                 'max_perc_COPD'],
                                    index=['gene'])

    # CONTROL
    pat_to_discard = []
    for i_pat_control in range(0, len(patient_control)):
        if read_pd:
            path_full_control = all_cells_path + ct + '_' + patient_control[
            i_pat_control]
            ending = '.tsv'
        else:
            path_full_control = all_cells_path + ct + '_' + patient_control[
            i_pat_control]
            ending = '.npy'

        if not os.path.exists(path_full_control + ending):
            pat_to_discard.append(patient_control[i_pat_control])
            warnings.warn('Patient ' + patient_control[i_pat_control] +
                          ' does not exist. Will be discarded! (create_patient_list)')
        else:
            if read_pd:
                all_cells_pat_control = pd.read_csv(
                    path_full_control + ending, sep="\t", index_col=0,
                    skiprows=ind_filtered_genes[row_gene], nrows=1)
                control_np = all_cells_pat_control.to_numpy()[0]
            else:
                all_cells_pat_control_full = \
                    as_numpy.read_numpy_to_df(path_full_control)
                # get row of interest
                all_cells_pat_control = \
                    all_cells_pat_control_full.iloc[ind_filtered_genes[row_gene],:] # is Series
                control_np = np.transpose(all_cells_pat_control.to_numpy())
                # print('control numpy pandas and numpy numpy equal: ' +
                #       str(np.array_equal(control_np,control_np2)))
            # else:
            #     print('ddfs')


            # nonzero values ( only expressed values )
            control_np_nonzero = control_np[np.nonzero(control_np)]

            # if no values are available, discard patient, such that wilcoxon rank
            # sum test can be calculated, otherwise there will be written one NaN
            # and then the whole method doesnt work for this gene.
            if np.shape(control_np)[0] != 0:
                patient_list.append(control_np)
                patient_list_nonzero.append(control_np_nonzero)
            # discard patient from patient-name-list
            else:
                pat_to_discard.append(patient_control[i_pat_control])

            len_pat_control = len(patient_list)

            # ######
            # number of cells: all and only expressed / percentage
            # get number of all measured cells
            #d = all_cells.iloc[i_gene, all_cells.columns.str.contains('Pt')].values

            # get number of nonzero values (only expressed cells)
            #d_nonzero = d[np.nonzero(d)]
            pd_sizes.iloc[i_pat_control, 0] = np.shape(control_np)[0]
            pd_sizes.iloc[i_pat_control, 1] = np.shape(control_np_nonzero)[0]
            if np.shape(control_np)[0] == 0:
                pd_sizes.iloc[i_pat_control, 2] = 0
            else:
                pd_sizes.iloc[i_pat_control, 2] = np.shape(control_np_nonzero)[0] * 100 / np.shape(control_np)[0]
            # #######

    # if we have patients which should be discarded (saved in list: pat_to_discard)
    if len(pat_to_discard) != 0:
        patient_control = [ele for ele in patient_control if ele not in pat_to_discard]

    pd_sizes_summary['mean_percentage_control'] = np.mean(pd_sizes['percentage'])
    pd_sizes_summary['min_perc_control'] = np.min(pd_sizes['percentage'])
    pd_sizes_summary['max_perc_control'] = np.max(pd_sizes['percentage'])




    # COPD ------------------------------------------------------------------
    # create a list with np.arrays of all copd patients, ONE GENE
    len_pat_copd_bevor = len(patient_copd)
    pat_to_discard = []
    for i_pat_copd in range(0, len(patient_copd)):
        if read_pd:
            path_full_copd = all_cells_path + ct + '_' + patient_copd[i_pat_copd]
            ending = '.tsv'
        else:
            path_full_copd = all_cells_path + ct + '_' + patient_copd[i_pat_copd]
            ending = '.npy'

        if not os.path.exists(path_full_copd + ending):
            pat_to_discard.append(patient_copd[i_pat_copd])
            warnings.warn('Patient ' + patient_copd[i_pat_copd] +
                          ' does not exist. Will be discarded! (create_patient_list)')
        else:
            if read_pd:
                all_cells_pat_copd = pd.read_csv(
                    path_full_copd + ending, sep="\t", index_col=0,
                    skiprows=ind_filtered_genes[row_gene], nrows=1)
                copd_np = all_cells_pat_copd.to_numpy()[0]
            else:
                all_cells_pat_copd_full = \
                    as_numpy.read_numpy_to_df(path_full_copd)
                # get row of interest
                all_cells_pat_copd = \
                    all_cells_pat_copd_full.iloc[ind_filtered_genes[row_gene], :]  # is Series
                copd_np = all_cells_pat_copd.to_numpy()

                # np.array_equal(copd_np,copd_np2)
            # else:
            #     print('ldfsk')

            # all_cells_pat_copd = pd.read_csv(
            #     path_full_copd, sep="\t", index_col=None, nrows=1)
            # ncol = np.shape(all_cells_pat_copd)[1]
            # all_cells_pat_copd = pd.read_csv(
            #     path_full_copd, sep="\t", index_col=0, usecols=range(1, ncol),
            #     skiprows=ind_filtered_genes[row_gene], nrows=1)
            # copd_np = all_cells_pat_copd.to_numpy()[0]

            # nonzero values ( only expressed values )
            copd_np_nonzero = copd_np[np.nonzero(copd_np)]

            # if no values are available, discard patient, such that wilcoxon rank
            # sum test can be calculated, otherwise there will be written one NaN
            # and then the whole method doesnt work for this gene.
            if np.shape(copd_np)[0] != 0:
                patient_list.append(copd_np)
                patient_list_nonzero.append(copd_np_nonzero)
            else:
                # merke patient to discard patient from patient-name-list
                # patient_copd.remove(patient_copd[i_pat_copd])
                pat_to_discard.append(patient_copd[i_pat_copd])

            # how many copd patients are actually read in
            len_pat_copd = len(patient_list) - len_pat_control

            # ######
            # number of cells: all and only expressed / percentage
            # get number of all measured cells
            # d = all_cells.iloc[i_gene, all_cells.columns.str.contains('Pt')].values
            # get number of nonzero values (only expressed cells)
            # d_nonzero = d[np.nonzero(d)]
            #pd_sizes.iloc[i_pat_copd, 0] = np.shape(control_np)[0]
            #pd_sizes.iloc[i_pat_copd, 1] = np.shape(control_np_nonzero)[0]
            if np.shape(copd_np)[0] == 0:
                pd_sizes.loc[patient_copd[i_pat_copd], 'percentage_copd'] = 0
            else:
                pd_sizes.loc[patient_copd[i_pat_copd], 'percentage_copd'] = \
                    np.shape(copd_np_nonzero)[0] * 100 / np.shape(copd_np)[0]

    # if we have patients which should be discarded (saved in list: pat_to_discard)
    if len(pat_to_discard) != 0:
        patient_copd = [ele for ele in patient_copd if ele not in pat_to_discard]

    pd_sizes_summary['mean_percentage_COPD'] = np.mean(
        pd_sizes['percentage_copd'])
    pd_sizes_summary['min_perc_COPD'] = np.min(pd_sizes['percentage_copd'])
    pd_sizes_summary['max_perc_COPD'] = np.max(pd_sizes['percentage_copd'])

    return patient_list, patient_list_nonzero, pd_sizes_summary, \
        len_pat_control, len_pat_copd, patient_control, patient_copd
