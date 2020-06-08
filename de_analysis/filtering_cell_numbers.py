# ############################################################################
# Filtering, cell numbers
# 1) get minimum of number of all existing cells between patients (for each
#    gene)
# 2) define threshold (e.g. 0.25 percent)
# 3) calculate 25% of min number of cells between patients (see (1))
# 4) if number of expressed cells from all patients below 25% -> discard this
#    gene
# get indices of genes that are higher then threshold
# save matrices of all cells with filtered genes
# ############################################################################
import pandas as pd
import numpy as np
import sys
import warnings
import os
import time
#
run_on_grid = False
if run_on_grid:
    import as_numpy
else:
    from de_analysis import as_numpy
# from de_analysis import *
# import de_analysis

def filtering_cell_numbers(patients_group1,
                           patients_group2,
                           all_cells_path,
                           ct,
                           where_to_save,
                           gene_from_row,
                           percent,
                           filtering_method,
                           read_pd):
    """
    Filtering, cell numbers
    1) get minimum of number of all existing cells between patients (for each
    gene)
    2) define threshold (e.g. 0.25 percent)
    3) calculate 25% of min number of cells between patients (see (1))
    4) if number of expressed cells from all patients below 25% -> discard this
    gene
    get indices of genes that are higher then threshold
    save matrices of all cells with filtered genes

    Args:
        patients_group1: list
            list of patient names belonging to the Control Group
        patients_group2: list
            list of patient names belonging to the Diseased COPD Group
        all_cells_path: string
            path to the matrices all_cells_... (all cells belonging to one
            patient and to one cluster)
        ct: string
            celltype which should be addressed in this analysis run
        where_to_save: string
            path to directory where results (filtered matrices) should be saved
        gene_from_row: int
            starting gene index = row-number from which the analysis should be
            run, here only important, if row=0 (for first calculation)
            then save filtered matrices
        percent: float
            filtering genes with too low number of expressed cells (percentage
            which should be filtered out) (e.g. if number of expressed cells
            < 25% of total number of expressed cells -> filter out)

    Returns:
        ind_filtered_genes: np.array
            indices of genes which should be discarded
        index_genes: pd.Series
            Series of genes with gene names
    """

    if read_pd:
        print('filtering - read pandas.........................')
    else:
        print('filtering - read numpy.........................')

    patient_full = patients_group1 + patients_group2
    number_cells_insg = np.zeros(len(patient_full))
    discard_patG1 = 0
    discard_patG2 = 0
    # create df: number of non-zero values for each gene, rows: genes,
    # columns: patients
    for i_pat in range(0, len(patient_full)):
        if read_pd:
            path_full_control = all_cells_path + ct + '_' + patient_full[i_pat]
            ending = '.tsv'
        else:
            path_full_control = all_cells_path + ct + '_' + patient_full[i_pat]
            ending = '.npy'

        if not os.path.exists(path_full_control + ending):
            warnings.warn('Patient ' + patient_full[i_pat] +
                          ' does not exist. Will be discarded!')
            print(path_full_control + ' does not exist!')
            if patient_full[i_pat] in patients_group1:
                discard_patG1 = discard_patG1 + 1
            elif patient_full[i_pat] in patients_group2:
                discard_patG2 = discard_patG2 + 1

        else:
            s1 = time.time()
            if read_pd:
                data_all_cells = pd.read_csv(
                    path_full_control + ending, sep="\t", index_col=0)
            else:
                # data_all_cells = read_numpy_to_df(path_full_control)
                data_all_cells = as_numpy.read_numpy_to_df(path_full_control)
            s2 = time.time()
            # print('------------')
            # print('reading time 1: ' + str(s2-s1))

            #a = data_all_cells.columns.str.contains('Pt')
            #anr = np.sum(data_all_cells.columns.str.contains('Pt'))
            #number_cells_insg[i_pat] = np.shape(data_all_cells)[1]  # zählt auch die ersten columns mit Gennamen/Indizes mit
            # number_cells_insg[i_pat] = np.sum(data_all_cells.columns.str.contains(col_stamp))    # wieviele Zellen gibt es insgesamt für einen Patienten (schaue nur Columns an die mit 'Pt' beschriftet sind)
            s1 = time.time()
            number_cells_insg[i_pat] = len(data_all_cells.columns)  # wieviele Zellen gibt es insgesamt für einen Patienten (schaue nur Columns an die mit 'Pt' beschriftet sind)

            # number of non-zero values (expressed values) in each row (for all genes at the same time)
            # Ser_nonzero = data_all_cells.iloc[:,
            #               data_all_cells.columns.str.contains(col_stamp)].astype(bool).sum(axis=1)  # Series
            Ser_nonzero = data_all_cells.astype(bool).sum(axis=1)  # Series
            s2 = time.time()
            # print('time #Zellen, #expressed cells: ' + str(s2 - s1))

            colname = ct + '_' + patient_full[i_pat]
            if i_pat == 0:
                s1 = time.time()
                df_nonzero = Ser_nonzero.to_frame(name=colname)
                s2 = time.time()
                # print ('time Series to frame: '+ str(s2-s1))
            else:
                s1 = time.time()
                df_nonzero2 = Ser_nonzero.to_frame(name=colname)
                df_nonzero = pd.concat([df_nonzero, df_nonzero2], axis=1,
                                       sort=False)
                s2 = time.time()
                # print('time Series to frame  + concat: ' + str(s2 - s1))
            # percentage: expressed cells in comparison to all cells :
            # (Ser_nonzero) * 100 /  number_cells_insg
            # = Series, for each gene: given:Percentage
            percentage_expr = Ser_nonzero * 100 / number_cells_insg[i_pat]
            #print(percentage_expr)
            if i_pat == 0:
                s1 = time.time()
                perc_expr_allpat_mtx = percentage_expr.to_frame(name=colname)
            else:
                s1 = time.time()
                perc_expr_allpat_mtx2 = percentage_expr.to_frame(name=colname)
                perc_expr_allpat_mtx = pd.concat([perc_expr_allpat_mtx, perc_expr_allpat_mtx2], axis=1,
                                       sort=False)
            s2 = time.time()
            # print('time percentage_expr.to_frame and concat: ' + str(s2 - s1))

    # Calculate mean over columns (control & copd)
    # first separate whole matrix into control&copd
    perc_expr_control = perc_expr_allpat_mtx.iloc[:,range(0,len(patients_group1)-discard_patG1)]
    perc_expr_copd = perc_expr_allpat_mtx.iloc[:,range(len(patients_group1)-discard_patG1,
                                                       len(patients_group1+patients_group2)-discard_patG1-discard_patG2)]
    mean_perc_control = perc_expr_control.mean(axis=1)
    mean_perc_copd = perc_expr_copd.mean(axis=1)
    # concat the two mean Series for the threshold condition
    s1 = time.time()
    mean_perc  = pd.concat([mean_perc_control,mean_perc_copd],axis=1)
    s2 = time.time()
    # print('time concat mean matrices over columns: ' + str(s2 - s1))
    # if head_genecol == 'Unnamed: 0':
    #     index_genes = data_all_cells['Unnamed: 0.1']
    #     usecols_start = 1
    # elif np.isnan(head_genecol):
    #     index_genes = data_all_cells['Unnamed: 0']

    index_genes = data_all_cells.index.values
    # usecols_start = 0

    # minimum of number of all existing cells between patients
    # min_num_cells = np.amin(number_cells_insg)

    if filtering_method == 'Method1':
        # -------- FILTERING 1 --
        # -- calculate percentage of expressed cells per patient,
        # -- calculate mean percentage for group1 & group2
        # -- at least one mean percentage (of group1 OR group2 is over a given
        #   threshold) -> keep gene
        # which values are below threshold -> boolean vec
        s1 = time.time()
        pd_bool_mean = mean_perc < percent*100
        # if all values in a row (both mean values) below threshold
        pd_bool_mean_all = pd_bool_mean.all(axis=1)

        # get reversed bool_vec, to show the values, that arent discarded
        reverse_bool_mean_ind = ~pd_bool_mean_all
        # remove gene index names to get the number
        pd_bool_mean_all = pd_bool_mean_all.reset_index(drop=True)
        reverse_bool_mean_ind = reverse_bool_mean_ind.reset_index(drop=True)
        ind_filtered_genes = reverse_bool_mean_ind.index[
            reverse_bool_mean_ind == True].values
        ind_filtered_genes_to_exclude = pd_bool_mean_all.index[
            pd_bool_mean_all == True].values
        s2 = time.time()
        # print('time get filtered indices: ' + str(s2 - s1))
        # ind_filtered_genes_to_exclude_mean = pd_bool_mean_all.index[
        #    pd_bool_mean_all == True].values
    # -------------------------------------------------------------
    # #########################################
    # ind_filtered_genes gives indices e.g. 7144 which is in one or both cases over the threshold -> so keep that gene
    # ######################################
    # -----------------------------------------------------

    elif filtering_method == 'Method2':
        # ------- FILTERING 2 -- if for all patients the number of expressed
        # cells is below a given threshold (threshold = minimum of number of
        # cells from all patients * percentage) -> discard gene
        thres = round(np.amin(number_cells_insg) * percent)
        # number counts (how many cells were measured) below a certain threshold
        if thres == 1:
            pd_bool = df_nonzero <= thres
        else:
            pd_bool = df_nonzero < thres
        # if all values in a row (for all patients) below threshold
        pd_bool_all = pd_bool.all(axis=1)
        # indices of genes, where the number of measured cells is below the threshold
        pd_ind = pd_bool_all.index[pd_bool_all == True].values

        # get reversed bool_vec, to show the values, that are not discarded
        reverse_bool_ind = ~pd_bool_all
        pd_bool_all = pd_bool_all.reset_index(drop=True)
        reverse_bool_ind = reverse_bool_ind.reset_index(drop=True)
        ind_filtered_genes = reverse_bool_ind.index[
            reverse_bool_ind == True].values
        ind_filtered_genes_to_exclude = pd_bool_all.index[
            pd_bool_all == True].values
        # print(len(ind_filtered_genes))

    elif filtering_method == 'Filtering3':
        # ---- OPTION 3:
        # --- given genes, for which the DE-Analysis should be done
        # -----
        # -----
        # get the genes that should be calculated from de-lists after
        # hierarchical clustering (genes with p-val=0.5 & WilcSc = 0)
        # TODO: change filepath / commenting
        filepath = '/home/erika/PycharmProjects/Kevin_GeneSetEnrichmentAnalysis/' \
                   'grid_results_DE_Okt_updated/merged_matrices/new_ordering/' \
                   'new_ordering_correctNames/prepros_hier_clustering/' \
                   'mind4Clusters/' + ct +\
                   '_10mean_perc_p_thr_0.75_w_thr_top500_mind4clusters.csv'
        print(filepath)
        de_li = pd.read_csv(filepath, sep=',', index_col=0)

        booli = (de_li['p_val_medianWilc'] == 0.5) & (
                    de_li['median_wilc_score'] == 0.0)
        genes_oi = de_li.index.values[booli]
        #check where genes_oi are located in data_alls_cells -> get indices
        #TODO:
        if head_genecol == 'Unnamed: 0':
            bool_data_all_cells = data_all_cells['Unnamed: 0.1'].isin(genes_oi)
        elif np.isnan(head_genecol):
            bool_data_all_cells = data_all_cells['Unnamed: 0'].isin(genes_oi)

        # get indices of  genes that are not discarded
        ind_filtered_genes = bool_data_all_cells.index[
            bool_data_all_cells == True].values
        # get indices of genes that are discarded
        rev_bool = ~bool_data_all_cells
        ind_filtered_genes_to_exclude = rev_bool.index[
            rev_bool == True].values

        if len(ind_filtered_genes) == 0:
            print('no extra genes to be calculated, for all genes a DE-Analysis was done!')
            sys.exit()
    # ############################################################
    # --------------------------------------------------

    #np.savetxt(where_to_save + 'gene_idices_to_read_' + ct + '_' + str(
    #    percent) + '.csv',
    #           ind_filtered_genes, delimiter=',')  # X is an array
    #
    # save filtered matrices
    if gene_from_row == 0:
        for i_pat in range(0, len(patient_full)):
            if read_pd:
                path_full = all_cells_path + ct + '_' + patient_full[i_pat]
                ending = '.tsv'
            else:
                path_full = all_cells_path + ct + '_' + patient_full[i_pat]
                ending = '.npy'
            # read first row to get number of columns
            if not os.path.exists(path_full + ending):
                ncol = 0
            else:
                start2 = time.time()
                if read_pd:
                    all_cells_pat_control = pd.read_csv(
                        path_full + ending, sep="\t", index_col=0, nrows=1)
                else:
                    # all_cells_pat_control = read_numpy_to_df(
                    #     path_full)
                    all_cells_pat_control= as_numpy.read_numpy_to_df(path_full)
                end2 = time.time()
                # print('reading time 2 (1row):' + str(end2-start2))

                ncol = np.shape(all_cells_pat_control)[1]
            # read filtered rows
            # if all genes are filtered, throw warning
            if len(ind_filtered_genes_to_exclude) == len(index_genes):
                warnings.warn(
                    "All genes were filtered out. Try to use a lower "
                    "filtering percentage or a different filtering approach."
                    "Stopping.")
            else:
                if os.path.exists(path_full + ending):
                    start3= time.time()
                    if read_pd:
                        filtered_all_cells_pat = pd.read_csv(
                            path_full + ending, sep="\t", index_col=0, header=0,
                            skiprows=ind_filtered_genes_to_exclude + 1)  # ,usecols=usecols_start,ncol)
                    else:
                        # f_all_cells_pat = read_numpy_to_df(path_full)
                        f_all_cells_pat = as_numpy.read_numpy_to_df(path_full)
                        filtered_all_cells_pat = \
                            f_all_cells_pat.iloc[ind_filtered_genes,:]
                    end3 = time.time()
                    # print('reading 3: only filtered rows ' + str(end3-start3))

                    if read_pd:
                        filtered_all_cells_pat.to_csv(
                            where_to_save + 'allCells_Filtered_' + str(percent) +
                            'genes_'+ ct + '_' + patient_full[i_pat] + '.tsv',
                            sep = '\t')
                    else:
                        # save_as_np(filtered_all_cells_pat,
                        #                     where_to_save + 'allCells_Filtered_'
                        #                     + str(percent) + 'genes_' + ct +
                        #                     '_' + patient_full[i_pat])
                        as_numpy.save_as_np(filtered_all_cells_pat,
                                            where_to_save + 'allCells_Filtered_'
                                            + str(percent) + 'genes_'+ ct +
                                            '_' + patient_full[i_pat])

    print(str(np.shape(ind_filtered_genes)[0]) +
          ' genes will be kept, from total ' + str(np.shape(index_genes)[0]) +
          '. That means, ' +
          str(np.shape(index_genes)[0]-np.shape(ind_filtered_genes)[0]) + \
          ' were filtered out.')

    return ind_filtered_genes, index_genes
