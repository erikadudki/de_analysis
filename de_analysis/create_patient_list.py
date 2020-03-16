# READ DATA ---------------------------------------------------------------
# create a list with np.arrays of all control patients and COPD patients,
# ONE GENE (as input for Wilcoxon rank sum test)
# ###########################################################################
import numpy as np
import pandas as pd
import os
import warnings


def create_patient_list(patient_control,
                        patient_copd,
                        all_cells_path,
                        ct,
                        ind_filtered_genes,
                        row_gene):
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
        row_gene: TODO

    Returns: TODO
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

    patient_list = []
    patient_list_nonzero = []
    pd_sizes = pd.DataFrame([], columns=['total_cells', 'expressed_cells',
                                         'percentage','percentage_copd'],
                            index=patient_control)  # ???
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
        path_full_control = all_cells_path + ct + '_' + patient_control[
            i_pat_control] + '.tsv'

        # all_cells_pat_control = pd.read_csv(
        #     path_full_control, sep=",", index_col=0, nrows=1)
        # ncol = np.shape(all_cells_pat_control)[1]
        #
        # all_cells_pat_control = pd.read_csv(
        #     path_full_control, sep=",", index_col=0,
        #     usecols=range(1, ncol),
        #     skiprows=ind_filtered_genes[row_gene], nrows=1)
        if not os.path.exists(path_full_control):
            pat_to_discard.append(patient_control[i_pat_control])
            warnings.warn('Patient ' + patient_control[i_pat_control] +
                          ' does not exist. Will be discarded! (create_patient_list)')
        else:
            all_cells_pat_control = pd.read_csv(
                path_full_control, sep="\t", index_col=0,
                skiprows=ind_filtered_genes[row_gene], nrows=1)

            control_np = all_cells_pat_control.to_numpy()[0]

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


            #print(i_pat_control)
            #print(len(patient_list))
            #print(len(patient_control))
            #print('--')

            len_pat_control = len(patient_list)

            # ###########
            # ges= [4407,1722,1096,1049,789,453]
            # #nonzerocells =  [205,164,112,89,68,33]
            # nonzerocells = [1465,529,492,207,203,166]
            # cells1 = np.random.normal(2.5, 1, nonzerocells[i_pat_control])
            # cells1_zeros = np.zeros(ges[i_pat_control]-nonzerocells[i_pat_control])
            # control_np = np.concatenate((cells1, cells1_zeros))
            # ############


            # ######
            # number of cells: all and only expressed / percentage
            # get number of all measured cells
            #d = all_cells.iloc[i_gene, all_cells.columns.str.contains('Pt')].values
            #print(np.shape(control_np))
            # get number of nonzero values (only expressed cells)
            #d_nonzero = d[np.nonzero(d)]
            #print(np.shape(control_np_nonzero))
            pd_sizes.iloc[i_pat_control, 0] = np.shape(control_np)[0]
            pd_sizes.iloc[i_pat_control, 1] = np.shape(control_np_nonzero)[0]
            if np.shape(control_np)[0] == 0:
                pd_sizes.iloc[i_pat_control, 2] = 0
            else:
                pd_sizes.iloc[i_pat_control, 2] = np.shape(control_np_nonzero)[0] * 100 / np.shape(control_np)[0]
            # #######

    #if we have patients which should be discarded (saved in list: pat_to_discard)
    if len(pat_to_discard) != 0:
        patient_control = [ele for ele in patient_control if ele not in pat_to_discard]

    pd_sizes_summary['mean_percentage_control'] = np.mean(pd_sizes['percentage'])
    pd_sizes_summary['min_perc_control'] = np.min(pd_sizes['percentage'])
    pd_sizes_summary['max_perc_control'] = np.max(pd_sizes['percentage'])

    # COPD
    # create a list with np.arrays of all copd patients, ONE GENE
    len_pat_copd_bevor = len(patient_copd)
    pat_to_discard = []
    for i_pat_copd in range(0, len(patient_copd)):
        path_full_copd = all_cells_path + ct + '_' + patient_copd[i_pat_copd] \
                         + '.tsv'
        if not os.path.exists(path_full_copd):
            pat_to_discard.append(patient_copd[i_pat_copd])
            warnings.warn('Patient ' + patient_copd[i_pat_copd] +
                          ' does not exist. Will be discarded! (create_patient_list)')
        else:
            all_cells_pat_copd = pd.read_csv(
                path_full_copd, sep="\t", index_col=None, nrows=1)
            ncol = np.shape(all_cells_pat_copd)[1]
            all_cells_pat_copd = pd.read_csv(
                path_full_copd, sep="\t", index_col=0, usecols=range(1, ncol),
                skiprows=ind_filtered_genes[row_gene], nrows=1)

            copd_np = all_cells_pat_copd.to_numpy()[0]
            #print('shape allcells copd:')
            #print(np.shape(copd_np)[0])

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
                #patient_copd.remove(patient_copd[i_pat_copd])
                pat_to_discard.append(patient_copd[i_pat_copd])

            #print(i_pat_copd)
            #print(len(patient_list))

            # how many copd patients are actually read in
            len_pat_copd = len(patient_list) - len_pat_control
            #print(len_pat_copd)
            #print(len(patient_copd))
            #print('--')

            # ###########################
            # ges = [5822,1498,1165,1009,728,357,384,361,330]
            # nonzerocells = [434,278,263,113,120,115,74,63,35]
            # nonzerocells = [2171,692,457,471,405,209,254,143,126]
            # cells2 = np.random.normal(2.5, 1, nonzerocells[i_pat_copd])
            # cells2_zeros = np.zeros(ges[i_pat_copd]-nonzerocells[i_pat_copd])
            # copd_np = np.concatenate((cells2, cells2_zeros))
            # ############################


            # ######
            # number of cells: all and only expressed / percentage
            # get number of all measured cells
            # d = all_cells.iloc[i_gene, all_cells.columns.str.contains('Pt')].values
            #print(np.shape(copd_np))
            # get number of nonzero values (only expressed cells)
            # d_nonzero = d[np.nonzero(d)]
            #print(np.shape(copd_np_nonzero))
            #pd_sizes.iloc[i_pat_copd, 0] = np.shape(control_np)[0]
            #pd_sizes.iloc[i_pat_copd, 1] = np.shape(control_np_nonzero)[0]
            if np.shape(copd_np)[0] == 0:
                pd_sizes.loc[patient_copd[i_pat_copd], 'percentage_copd'] = 0
            else:
                pd_sizes.loc[patient_copd[i_pat_copd], 'percentage_copd'] = \
                    np.shape(copd_np_nonzero)[0] * 100 / np.shape(copd_np)[0]

    #if we have patients which should be discarded (saved in list: pat_to_discard)
    if len(pat_to_discard) != 0:
        patient_copd = [ele for ele in patient_copd if ele not in pat_to_discard]
    #print(patient_copd)

    pd_sizes_summary['mean_percentage_COPD'] = np.mean(
        pd_sizes['percentage_copd'])
    pd_sizes_summary['min_perc_COPD'] = np.min(pd_sizes['percentage_copd'])
    pd_sizes_summary['max_perc_COPD'] = np.max(pd_sizes['percentage_copd'])

    return patient_list, patient_list_nonzero, pd_sizes_summary, \
        len_pat_control, len_pat_copd, patient_control, patient_copd
