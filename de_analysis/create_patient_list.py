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

def sort_patient_names(patient_control, patient_copd):
    pat_control_sorted = []
    pat_copd_sorted = []
    pat_list_final = []
    pat_unique = []

    # get patient-identifiers ('Pt1') and sort the list according to alphabet
    for i_p in range(len(patient_control)):
        pat_control_sorted.append(patient_control[i_p].split("_")[0])
        pat_control_sorted.sort()
        cluster_control = patient_control[i_p].split("_")[1]
    for i_p in range(len(patient_copd)):
        pat_copd_sorted.append(patient_copd[i_p].split("_")[0])
        pat_copd_sorted.sort()
        cluster_copd = patient_copd[i_p].split("_")[1]
    # check that entries are the same / same order
    if len(pat_copd_sorted) > len(pat_control_sorted):
        longer_list = pat_copd_sorted
        shorter_list = pat_control_sorted
    elif len(pat_copd_sorted) < len(pat_control_sorted):
        longer_list = pat_control_sorted
        shorter_list = pat_copd_sorted
    else:
        # TODO
        print("TODO... if length of input patients in Group1 equals "
                   "length of input patients in Group2, but patients are not "
              "the same...")
        print(patient_copd)
        print(patient_control)
    
    if len(pat_copd_sorted) != len(pat_control_sorted):
        for i_p in range(len(longer_list)):
            if longer_list[i_p] in shorter_list:
                pat_list_final.append(longer_list[i_p])
            else:
                pat_unique.append(longer_list[i_p])
        # entries that only exist in one patient list, append to the end
        pat_list_final.extend(pat_unique)

    if len(pat_copd_sorted) > len(pat_control_sorted):
        pat_copd_final = pat_list_final.copy()
        pat_control_final = pat_control_sorted.copy()
    elif len(pat_copd_sorted) < len(pat_control_sorted):
        pat_copd_final = pat_copd_sorted.copy()
        pat_control_final = pat_list_final.copy()
    else:
        print('todo')
        # TODO
    
    if len(pat_copd_sorted) == len(pat_control_sorted):
        pat_copd_final = pat_copd_sorted
        pat_control_final = pat_control_sorted
        # if entries in both groups are not the same, TODO
        if pat_copd_sorted != pat_control_sorted:
            raise ValueError("TODO! If we have same number of patients, but"
                             "patients are not the same (eg. for cluster "
                             "comparison they should be the same, but with"
                             "different cell clusters). If we have different"
                             "patients, then permutation test needs to be "
                             "adjusted!")
    # add cluster names to patient-names
    pat_copd_final = [ele + '_' + cluster_copd for ele in pat_copd_final]
    pat_control_final = [ele + '_' + cluster_control for ele in pat_control_final]

    return pat_control_final, pat_copd_final




def create_patient_list(patient_control,
                        patient_copd,
                        all_cells_path,
                        ct,
                        ind_filtered_genes,
                        row_gene,
                        read_pd,
                        perm_modus='twosided'):
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
    if read_pd:
        ending = '.tsv'
    else:
        ending = '.npy'

    # check that order of patients is the same in both groups, or reorder
    if perm_modus == 'compare_clusters':
        patient_control, patient_copd = \
            sort_patient_names(patient_control,patient_copd)

    # Check if Patients-files exist and create list of patients that exist!!!
    pat_to_discard = []
    pat_to_keep_control = []
    for i_pat_control in range(0, len(patient_control)):
        path_full_control = all_cells_path + ct + '_' + patient_control[
            i_pat_control]

        if not os.path.exists(path_full_control + ending):
            pat_to_discard.append(patient_control[i_pat_control])
            warnings.warn('Patient ' + patient_control[i_pat_control] +
                          ' does not exist. Will be discarded! (create_patient_list)')
        else:
            # check if there are entries
            if read_pd:
                try:
                    all_cells_pat_control = pd.read_csv(
                        path_full_control + ending, sep="\t", index_col=0,
                        usecols=[0,1,2],
                        nrows=1)
                except:
                    all_cells_pat_control = pd.read_csv(
                        path_full_control + ending, sep="\t", index_col=0,
                        nrows=1)
            else:
                all_cells_pat_control = \
                    as_numpy.read_numpy_to_df(path_full_control,
                                              read_only_col=True)
            if np.shape(all_cells_pat_control)[1] < 1:
                pat_to_discard.append(patient_control[i_pat_control])
                warnings.warn('Patient ' + patient_control[i_pat_control] +
                              ' has not cell expression entries. Table is empty. '
                              'Will be discarded! (create_patient_list)')
            else:
                pat_to_keep_control.append(patient_control[i_pat_control])

    # COPD ------------------------------------------------------------------
    # create a list with np.arrays of all copd patients, ONE GENE
    # Check if Patients-files exist and create list of patients that exist!!!
    pat_to_keep_copd = []
    for i_pat_copd in range(0, len(patient_copd)):
        path_full_copd = all_cells_path + ct + '_' + patient_copd[
            i_pat_copd]

        if not os.path.exists(path_full_copd + ending):
            pat_to_discard.append(patient_copd[i_pat_copd])
            warnings.warn('Patient ' + patient_copd[i_pat_copd] +
                          ' does not exist. Will be discarded! (create_patient_list)')
        else:
            # check if there are entries
            if read_pd:
                try:
                    all_cells_pat_copd = pd.read_csv(
                        path_full_copd + ending, sep="\t", index_col=0,
                        usecols=[0, 1, 2], nrows=1)
                except:
                    all_cells_pat_copd = pd.read_csv(
                        path_full_copd + ending, sep="\t", index_col=0,
                        nrows=1)
            else:
                all_cells_pat_copd = \
                    as_numpy.read_numpy_to_df(path_full_copd,
                                              read_only_col=True)
            if np.shape(all_cells_pat_copd)[1] < 1:
                pat_to_discard.append(patient_copd[i_pat_copd])
                warnings.warn('Patient ' + patient_copd[i_pat_copd] +
                              ' has not cell expression entries. Table is empty. '
                              'Will be discarded! (create_patient_list)')
            else:
                pat_to_keep_copd.append(patient_copd[i_pat_copd])

    # check order of both patient lists to keep, for the modus compare-clusters
    # and not same number of patients, it has to be in the same order, and all
    # the patients that only exist in one group needs to be in the end
    if perm_modus == 'compare_clusters':
        len_g1 = len(pat_to_keep_control)
        len_g2 = len(pat_to_keep_copd)
        print('pat_to_keep_G1: ' + str(pat_to_keep_control))
        print('pat_to_keep_G2: ' + str(pat_to_keep_copd))
        for i_check in range(0,min(len_g1,len_g2)):
            if not pat_to_keep_control[i_check] == pat_to_keep_copd[i_check]:
                ValueError("TODO: Implement the case, if the patients do not "
                           "coincide between the two groups. Could happen, if "
                           "on the way patients were discarded, because their"
                           "tables were empty. The patient lists should have"
                           "the same order and same first name. (for "
                           "the case 'compare clusters')")



    # ----------------------------------------------
    # Start constructing the patient lists with the cell expressions, with the
    # patients order we just generated.
    # -----------------------------------------------
    # CONTROL
    len_pat_control = len(pat_to_keep_control)
    for i_pat_control in range(0, len(pat_to_keep_control)):
        path_control = all_cells_path + ct + '_' + pat_to_keep_control[
            i_pat_control]
        if read_pd:
            all_cells_pat_control = pd.read_csv(
                path_control + ending, sep="\t", index_col=0,
                skiprows=ind_filtered_genes[row_gene], nrows=1)
            control_np = all_cells_pat_control.to_numpy()[0]
        else:
            all_cells_pat_control_full = \
                as_numpy.read_numpy_to_df(path_control)
            # get row of interest
            all_cells_pat_control = \
                all_cells_pat_control_full.iloc[ind_filtered_genes[row_gene],:] # is Series
            control_np = np.transpose(all_cells_pat_control.to_numpy())

        # nonzero values ( only expressed values )
        control_np_nonzero = control_np[np.nonzero(control_np)]

        # fill patient_list: actual list of cell expressions, for each patient
        patient_list.append(control_np)
        patient_list_nonzero.append(control_np_nonzero)
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
            pd_sizes.iloc[i_pat_control, 2] = np.shape(control_np_nonzero)[0] \
                                              * 100 / np.shape(control_np)[0]

    pd_sizes_summary['mean_percentage_control'] = np.mean(pd_sizes['percentage'])
    pd_sizes_summary['min_perc_control'] = np.min(pd_sizes['percentage'])
    pd_sizes_summary['max_perc_control'] = np.max(pd_sizes['percentage'])


    # COPD ------------------------------------------------------------------
    len_pat_copd = len(pat_to_keep_copd)
    for i_pat_copd in range(0, len(pat_to_keep_copd)):
        path_copd = all_cells_path + ct + '_' + pat_to_keep_copd[i_pat_copd]
        if read_pd:
            all_cells_pat_copd = pd.read_csv(
                path_copd + ending, sep="\t", index_col=0,
                skiprows=ind_filtered_genes[row_gene], nrows=1)
            copd_np = all_cells_pat_copd.to_numpy()[0]
        else:
            all_cells_pat_copd_full = \
                as_numpy.read_numpy_to_df(path_copd)
            # get row of interest
            all_cells_pat_copd = \
                all_cells_pat_copd_full.iloc[ind_filtered_genes[row_gene], :]  # is Series
            copd_np = all_cells_pat_copd.to_numpy()
            # copd_np = np.transpose(all_cells_pat_copd.to_numpy())#?????

        # nonzero values ( only expressed values )
        copd_np_nonzero = copd_np[np.nonzero(copd_np)]

        patient_list.append(copd_np)
        patient_list_nonzero.append(copd_np_nonzero)

        # ######
        # number of cells: all and only expressed / percentage
        # get number of all measured cells
        # d = all_cells.iloc[i_gene, all_cells.columns.str.contains('Pt')].values
        # get number of nonzero values (only expressed cells)
        # d_nonzero = d[np.nonzero(d)]
        pd_sizes.iloc[i_pat_copd, 0] = np.shape(control_np)[0]
        pd_sizes.iloc[i_pat_copd, 1] = np.shape(control_np_nonzero)[0]
        if np.shape(copd_np)[0] == 0:
            pd_sizes.loc[pat_to_keep_copd[i_pat_copd], 'percentage_copd'] = 0
        else:
            pd_sizes.loc[pat_to_keep_copd[i_pat_copd], 'percentage_copd'] = \
                np.shape(copd_np_nonzero)[0] * 100 / np.shape(copd_np)[0]

    pd_sizes_summary['mean_percentage_COPD'] = np.mean(
        pd_sizes['percentage_copd'])
    pd_sizes_summary['min_perc_COPD'] = np.min(pd_sizes['percentage_copd'])
    pd_sizes_summary['max_perc_COPD'] = np.max(pd_sizes['percentage_copd'])


    return patient_list, patient_list_nonzero, pd_sizes_summary, \
        len_pat_control, len_pat_copd, pat_to_keep_control, pat_to_keep_copd
