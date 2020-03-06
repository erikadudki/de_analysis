# ##########################################################################
# MAIN FUNCTION FOR DE-ANALYSIS
# includes filtering of genes
# performs Wilcoxon rank sum test
# performs permutation test
# input: Normalized count matrices called: CLUSTERPATIENT (Matrices for each
# cluster, for each patient) -> need to be files separated by ',' (no
# tab-separated files)
#       to prepare these files: read_data_build_groups.py
# TODO: add requirements, of how paths should be named
# ##########################################################################
import pandas as pd
import numpy as np
import scipy.stats
import time
import math
from get_perm_array_ind import get_perm_array_ind
from filtering_cell_numbers import filtering_cell_numbers
from create_patient_list import create_patient_list
import os

# calculating on grid, or on your laptop?
grid = False

# celltype / cluster (refers to filenames in './data/data_per_pat_per_cl/')
ct = '0'
ct = 'cDC2'
ct = 'Macrophage'


# patients_group1 = ['Pat135', 'Pat133', 'Pat139', 'Pat111', 'Pat141',
#                        'Pat172']
# patients_group2 = ['Pat137', 'Pat149', 'Pat115', 'Pat145', 'Pat155', 'Pat162',
#                     'Pat142', 'Pat173', 'Pat175']
patients_group1 = ['Pat135', 'Pat133', 'Pat111']

patients_group2 = ['Pat137', 'Pat149', 'Pat115', 'Pat145', 'Pat155']
# which patients are in which group?
# patients_group1 = ['Pt1','Pt2','Pt3']
# patients_group2 = ['Pt4','Pt5','Pt6']

# Filtering genes with to low number of expressed cells
percent = 0.10
percent = 0.10

# set manually which genes(rows) should be calculated
gene_from_row = 0
gene_until_row = 640
gene_until_row = 3  # now: set to len(filtered_index-genes)


wd = '/home/erika/PycharmProjects/DE-Analysis/src/code_tidy/'

#where_to_save = './results_new_annotation/'
# where_to_save = './results_new_ann_with_nonzero_eval/'
# where_to_save = '/home/erika/PycharmProjects/DE-Analysis/src/' \
#                 'code_tidy/tests/'
# where_to_save = "/home/erika/Documents/Projects/Muscat/muscat-comparison/data/" \
#                 "sim_data/kang/pandasDF/logcounts/nill;1/de-results/"
# where_to_save = wd + 'de_results/'

# if not os.path.exists(where_to_save):
#     os.mkdir(where_to_save)


#where_to_save = '/home/erika/PycharmProjects/Kevin_GeneSetEnrichmentAnalysis/' \
#                'grid_results_DE_Okt_updated/results_extra_genes/'
#where_to_save = '/home/erika/PycharmProjects/Kevin_GeneSetEnrichmentAnalysis/' \
#                'grid_results_DE_Nov_Eosinophil/'
#where_to_save = '/home/erika//PycharmProjects/DE-Analysis/src/code_tidy/tests/'
# all_cells_path = "./Kevins_new_annotation_August/AllCells_"
# all_cells_path = '/home/erika/PycharmProjects/Kevin_GeneSetEnrichmentAnalysis/' \
#                  'Kevins_new_subannotation_Okt/AllCells_'
#all_cells_path = '/home/erika/PycharmProjects/Kevin_GeneSetEnrichmentAnalysis/' \
#                 'Kevins_new_annotation_Neutrophil_Nov/AllCells_'
# all_cells_path = "/home/erika/PycharmProjects/Kevin_GeneSetEnrichmentAnalysis/" \
#                    "Kevins_new_annotation_Eosinophil_Nov/AllCells_"
# all_cells_path = '/home/erika/PycharmProjects/Kevin_GeneSetEnrichmentAnalysis/' \
#                'Kevins_new_annotation_Sept/all_cluster_together/AllCells_'
# all_cells_path = "/home/erika/Documents/Projects/Muscat/muscat-comparison/data/" \
#                 "sim_data/kang/pandasDF/logcounts/nill;1/"
# all_cells_path = wd + 'data/data_per_pat_per_cl/'

# indicate how your data-frames are annotated (columns have to have a identical
# indicator/stamp & what is the header of your gene column?

# in the data-pds: all columns (of data) have the following indicator, stamp
# col_stamp = 'Pt'
# col_stamp = 'cell'
# in the data-pds: the header of the gene column has the following indicator, stamp
# head_genecol = 'Unnamed: 0'
# head_genecol = np.nan

##
def de_analysis(ct,
                gene_from_row,
                gene_until_row,
                patients_group1,
                patients_group2,
                percent):
    """
    Input:
        ct: string
            cell type (refers to filenames in './data/data_per_pat_per_cl/'
        gene_from_row: int
            choose the starting rows of genes for which the DE-Analysis should
            be run
        gene_until_row: choose the ending row of genes for which the
            DE-Analysis should be stopped
        where_to_save: string
            path where results should be saved: working directory + 'de_results/'
        all_cells_path: string
            path where matrix: all_cells (all cells per patient per cell type)
            is located
        patients_group1: list
            list of patient names Group1(should be the same as in the
            all_cells_CELLTPYEPatient matrix)
        patients_group2: list
            list of patient names Group2 (should be the same as in the
            all_cells_CELLTPYEPatient matrix)

    Returns:
        --
        (saves DE-List)

    """


    #copd_chr_br = ['Pat115','Pat137','Pat145','Pat149','Pat155','Pat162']
    #copd_emph = ['Pat142','Pat173','Pat175']
    # saving input patients
    #np.save(where_to_save + ct + '_pat_control_input',patients_group1)

    # create saving directory:
    where_to_save = wd + 'de_results/'
    if not os.path.exists(where_to_save):
        os.mkdir(where_to_save)

    # define directory of normalized data matrices for each patient & cluster
    all_cells_path = wd + 'data/data_per_pat_per_cl/'

    f = open(where_to_save + ct + '_INFORMATION.txt', "w+")
    f.write('Patients Group 1 Input: ' + str(patients_group1) + '\n')
    f.write('Patients Group 2 Input: ' + str(patients_group2) + '\n')

    random_perm = False         # do only subset of random permutations? -> approximated nulldistribution
    all_permutations = True     # do all possible permutations?         -> exact null distribution

    # Filtering genes with to low number of expressed cells
    ind_filtered_genes, index_genes = filtering_cell_numbers(patients_group1,
                                                             patients_group2,
                                                             all_cells_path,
                                                             ct,
                                                             where_to_save,
                                                             gene_from_row,
                                                             percent)
    # how many genes remain after filtering: save

    if gene_from_row == 0:
        f = open(where_to_save + ct + '_INFORMATION.txt', "a+")
        f.write('Filtering with ' + str(percent*100) + '% means ' +
                str(len(ind_filtered_genes)) + ' genes will be kept.\n')

        # np.save(where_to_save + ct + str(percent) + '_len_' +
        #         str(len(ind_filtered_genes)),ind_filtered_genes)

    # gene_until_row = len(ind_filtered_genes)
    #
    # row_gene = 11819 #B2M
    # row_gene = 23334 #MT-ND3
    # row_gene = 23330 #MT-CO3
    test = False  # if True: dont save data, just to test code

    start_total = time.time()

    for row_gene in range(gene_from_row,
                          gene_until_row):  # len(ind_filtered_genes)):

        # READ DATA ---------------------------------------------------------------
        # create a list with np.arrays of all control+COPD patients, ONE GENE
        start_loop = time.time()
        patient_list, patient_list_nonzero, pd_sizes_summary,\
            len_pat_control, len_pat_copd, patients_group1, patients_group2 = \
            create_patient_list(patients_group1,
                                patients_group2,
                                all_cells_path,
                                ct,
                                ind_filtered_genes,
                                row_gene)

        f = open(where_to_save + ct + '_INFORMATION.txt', "a+")
        f.write('\n Patients with only zero cells for a cluster are discarded: \n')
        f.write('Gene ' + str(row_gene) + ': Patients Group 1 ' +
                str(patients_group1) + '\n')
        f.write('Gene ' + str(row_gene) + ': Patients Group 1 :' +
                str(patients_group2) + '\n')


        # with open(where_to_save + ct + str(percent) + '_pat_control_after.txt', "w") as output:
        #     output.write(str(patients_group1))
        # with open(where_to_save + ct + str(percent)+ '_pat_copd_after.txt',
        #           "w") as output:
        #     output.write(str(patients_group2))
        #np.savetxt(where_to_save + ct + '_pat_control_after', patients_group1)
        #np.savetxt(where_to_save + ct + '_pat_copd_after', patients_group2)

        end1 = time.time()

        if row_gene == gene_from_row:
            # Initialization of result matrices with P-values, Median-Wilc-Scores, etc.
            pd_wilc = pd.DataFrame(np.nan, index=range(0, len(ind_filtered_genes)),
                                   columns=range(0, len(patients_group1) * len(
                                       patients_group2)))
            # p_val_col = ['p_val_medianWilc', 'p_val_meanWilc', 'median_wilc_score',
            #              'mean_wilc_score', 'min_wilc_score', 'max_wilc_score',
            #              'time_read_in', 'time_Wilcoxon', 'time_permutation_test',
            #              'time_total', 'pval_medWilc_nonzero',
            #              'med_wilc_score_nonzero', 'mean_percentage_control',
            #              'mean_percentage_COPD', 'min_perc_control',
            #              'max_perc_control',
            #              'min_perc_COPD',
            #              'max_perc_COPD']  # careful! dont change order!
            p_val_col = ['p_val_medianWilc', 'p_val_meanWilc',
                         'median_wilc_score',
                         'mean_wilc_score', 'min_wilc_score', 'max_wilc_score',
                         'time_read_in', 'time_Wilcoxon',
                         'time_permutation_test',
                         'time_total', 'mean_percentage_control',
                         'mean_percentage_COPD', 'min_perc_control',
                         'max_perc_control',
                         'min_perc_COPD',
                         'max_perc_COPD']  # careful! dont change order!

            p_val_df = pd.DataFrame(np.nan,
                                    index=range(0, len(ind_filtered_genes)),
                                    columns=p_val_col)
            pd_fc_all_cells = pd.DataFrame(np.nan,
                                           index=range(0, len(ind_filtered_genes)),
                                           columns=range(0, len(
                                               patients_group1) * len(patients_group2)))
            pd_fc_expr_cells = pd.DataFrame(np.nan,
                                            index=range(0,
                                                        len(ind_filtered_genes)),
                                            columns=range(0,
                                                          len(patients_group1) * len(
                                                              patients_group2)))
            pd_fc_expr_cells_med = pd.DataFrame(np.nan, index=range(0, len(
                ind_filtered_genes)),
                                                columns=range(0, len(patients_group1) * len(
                                                    patients_group2)))

            # ------ finished initializing
        p_val_df.loc[row_gene, 'time_read_in'] = end1 - start_loop


        # add information about sizes/ percentage of number of expressed cells
        #p_val_df.iloc[row_gene,12:18] = pd_sizes_summary.values
        p_val_df.loc[row_gene,['mean_percentage_control',
                               'mean_percentage_COPD',
                               'min_perc_control',
                               'max_perc_control',
                               'min_perc_COPD',
                               'max_perc_COPD']] = pd_sizes_summary.values
    #
        #print('starting MAIN WILCOXON TEST')
        # MAIN WILCOXON  TEST-----------------------------------------------------
        # Wilcoxon Test 1pat_control vs 1pat_copd
        num_cross_val = len_pat_control* len_pat_copd
        #len(patients_group1) * len(patients_group2)

        Wilc_score = np.zeros(num_cross_val)
        Wilc_nonzero = np.zeros(num_cross_val)
        pval = np.zeros(num_cross_val)
        pval_nonzero = np.zeros(num_cross_val)
        col_names = np.empty(num_cross_val, dtype=object)
        fc_all_cells = np.zeros(num_cross_val)
        fc_expr_cells = np.zeros(num_cross_val)
        fc_expr_cells_med = np.zeros(num_cross_val)
        #
        run_idx = 0
        for i_pat_copd in range(len_pat_control,
                                len_pat_control + len_pat_copd):
            for i_pat_ctl in range(0, len_pat_control):
                Wilc_score[run_idx], pval[run_idx] = \
                    scipy.stats.ranksums(patient_list[i_pat_ctl],
                                         patient_list[i_pat_copd])
                Wilc_nonzero[run_idx], pval_nonzero[
                    run_idx] = scipy.stats.ranksums(
                    patient_list_nonzero[i_pat_ctl],
                    patient_list_nonzero[i_pat_copd])

                # Wilc_score[run_idx], pval[run_idx] = scipy.stats.mannwhitneyu(patient_list[i_pat_ctl], patient_list[i_pat_copd],
                #                        use_continuity=False,
                #                       alternative='two-sided')

                # fold change: mean over all cells (including zero values)
                # if np.mean(patient_list[i_pat_copd]) == 0:
                #     # TODO: ist diese Annahme korrekt?
                #     fc_all_cells[run_idx] = math.log(
                #         1e-50 /
                #         np.mean(patient_list[i_pat_ctl]), 2)
                # elif np.mean(patient_list[i_pat_ctl]) == 0:
                #     fc_all_cells[run_idx] = math.log(
                #         np.mean(patient_list[i_pat_copd]) /
                #         1e-50, 2)
                # else:

                fc_all_cells[run_idx] = math.log(
                    (np.mean(patient_list[i_pat_copd])+1) /
                     (np.mean(patient_list[i_pat_ctl])+1), 2)
                # print(np.mean(patient_list[i_pat_copd]))
                # print(np.mean(patient_list[i_pat_ctl]))
                # print('---')
                # fold change: mean over nonzero cells (only expressed cells)
                a = patient_list[i_pat_copd]
                b = patient_list[i_pat_ctl]
                fc_expr_cells[run_idx] = math.log(
                    (np.mean(a[np.nonzero(a)])+1) /
                    (np.mean(b[np.nonzero(b)])+1), 2)

                fc_expr_cells_med[run_idx] = math.log(
                    (np.median(a[np.nonzero(a)])+1) /
                    (np.median(b[np.nonzero(b)])+1), 2)

                # write column names
                if row_gene == 0:
                    col_names[run_idx] = patients_group1[i_pat_ctl] + '_vs_' + \
                                         patients_group2[i_pat_copd - len_pat_control]
                run_idx = run_idx + 1

        end2 = time.time()
        #p_val_df.iloc[row_gene, 7] = end2 - end1
        p_val_df.loc[row_gene, 'time_Wilcoxon'] = end2 - end1
        #
        min_max_wilc = np.zeros(2)

        min_max_wilc[0] = np.min(Wilc_score)
        min_max_wilc[1] = np.max(Wilc_score)
        median_wilc_score = np.median(Wilc_score)
        mean_wilc_score = np.mean(Wilc_score)
        median_wilc_nonzero = np.median(Wilc_nonzero)

        # save Wilcoxon scores & fold change values for all tests
        if row_gene == 0:
            pd_wilc.rename(
                index={row_gene: index_genes[ind_filtered_genes[row_gene]]},
                inplace=True)
            pd_wilc.columns = col_names
            pd_wilc.iloc[row_gene] = Wilc_score

            pd_fc_all_cells.rename(
                index={row_gene: index_genes[ind_filtered_genes[row_gene]]},
                inplace=True)
            pd_fc_all_cells.columns = col_names
            pd_fc_all_cells.iloc[row_gene] = fc_all_cells

            pd_fc_expr_cells.rename(
                index={row_gene: index_genes[ind_filtered_genes[row_gene]]},
                inplace=True)
            pd_fc_expr_cells.columns = col_names
            pd_fc_expr_cells.iloc[row_gene] = fc_expr_cells

            pd_fc_expr_cells_med.rename(
                index={row_gene: index_genes[ind_filtered_genes[row_gene]]},
                inplace=True)
            pd_fc_expr_cells_med.columns = col_names
            pd_fc_expr_cells_med.iloc[row_gene] = fc_expr_cells_med
            # pd_wilc.to_csv(where_to_save + 'wilc_scores_' + ct + '_filteredGenes0_25')
        else:
            pd_wilc.rename(
                index={row_gene: index_genes[ind_filtered_genes[row_gene]]},
                inplace=True)
            pd_wilc.iloc[row_gene] = Wilc_score

            pd_fc_all_cells.rename(
                index={row_gene: index_genes[ind_filtered_genes[row_gene]]},
                inplace=True)
            pd_fc_all_cells.iloc[row_gene] = fc_all_cells

            pd_fc_expr_cells.rename(
                index={row_gene: index_genes[ind_filtered_genes[row_gene]]},
                inplace=True)
            pd_fc_expr_cells.iloc[row_gene] = fc_expr_cells

            pd_fc_expr_cells_med.rename(
                index={row_gene: index_genes[ind_filtered_genes[row_gene]]},
                inplace=True)
            pd_fc_expr_cells_med.iloc[row_gene] = fc_expr_cells_med
            if not test:
                pd_wilc.to_csv(
                    where_to_save + 'wilc_scores_' + ct + '_filteredGenes' + str(
                        percent) + '_' + str(gene_from_row) + '-' + str(
                        gene_until_row))
                pd_fc_all_cells.to_csv(
                    where_to_save + 'fc_all_cells_mean' + ct + '_filteredGenes' + str(
                        percent) + '_' + str(gene_from_row) + '-' + str(
                        gene_until_row))
                pd_fc_expr_cells.to_csv(
                    where_to_save + 'fc_expr_cells_mean' + ct + '_filteredGenes' + str(
                        percent) + '_' + str(gene_from_row) + '-' + str(
                        gene_until_row))
                pd_fc_expr_cells.to_csv(
                    where_to_save + 'fc_expr_cells_median' + ct + '_filteredGenes' + str(
                        percent) + '_' + str(gene_from_row) + '-' + str(
                        gene_until_row))

        end = time.time()
        # print('time: min, max Wilc, save Wilc, fc')
        # print(end - end1)
        # print('time one gene:')
        # print(end - start_loop)

        # PERMUTATION-TEST---------------------------------------------------------
        # permute copd/control labels
        start = time.time()
        if random_perm:
            # num_perm = 5000     # 1.314912529786428 sec
            # num_perm = 2500     # 0.6059358835220336 sec
            # num_perm = 1000     # 0.27018970251083374 sec
            # num_perm = 500      # 0.1287842075030009 sec
            # num_perm = 250      # 0.06198491255442302   -> 37 std,36000genes   // 302sec,5517genes(0.2filtering)
            # num_perm = 350      # 0.07590539455413818 min -> 7 std, 5517 genes
            num_perm = 200  # time: one gene:3.203282117843628 -> 4.3 Std, 4834 genes
        elif all_permutations:
            #num_perm = math.factorial(
            #    len(patients_group1) + len(patients_group2)) / (
            #                   math.factorial(
            #                       len(patients_group1)) * math.factorial(len(
            #               patients_group2))) # 5005, ca. 60-80 sec per gene total/ max 19 sec read data
            num_perm = math.factorial(len_pat_control + len_pat_copd) / (math.factorial(
                                    len_pat_control) * math.factorial(len_pat_copd))
            num_perm = int(num_perm)

        Wilc_score_perm = np.zeros([num_perm, num_cross_val])
        Wilc_nonzero_perm = np.zeros([num_perm, num_cross_val])
        pval_perm = np.zeros([num_perm, num_cross_val])
        pval_nonzero_perm = np.zeros([num_perm, num_cross_val])

        for i_perm in range(0, num_perm):
            # permute patient_list
            if random_perm:
                patient_list_permuted = np.random.permutation(patient_list)
            elif all_permutations:
                # get array with all combinations of permuted indices
                if i_perm == 0:
                    array_perm_ind = get_perm_array_ind(len_pat_control, len_pat_copd)#len(patients_group1),
                                                        #len(patients_group2))  # 5005 x 15
                # rearrange the patient list with permuted indices
                myorder = array_perm_ind[i_perm]
                patient_list_permuted = [patient_list[i] for i in myorder]
                patient_l_permuted_nonzero = [patient_list_nonzero[i] for i in
                                              myorder]
            # do Wilcoxon Test on permuted patient list
            run_idx = 0
            for i_pat_copd in range(len_pat_control,
                                    len_pat_control + len_pat_copd):
                for i_pat_ctl in range(0, len_pat_control):
                    Wilc_score_perm[i_perm, run_idx], pval_perm[
                        i_perm, run_idx] = scipy.stats.ranksums(
                        patient_list_permuted[i_pat_ctl],
                        patient_list_permuted[i_pat_copd])
                    Wilc_nonzero_perm[i_perm, run_idx], pval_nonzero_perm[
                        i_perm, run_idx] = scipy.stats.ranksums(
                        patient_l_permuted_nonzero[i_pat_ctl],
                        patient_l_permuted_nonzero[i_pat_copd])
                    # Wilc_score_perm[i_perm,run_idx], pval_perm[i_perm,run_idx] = scipy.stats.mannwhitneyu(patient_list[i_pat_ctl], patient_list[i_pat_copd],
                    #                    use_continuity=False,
                    #                   alternative='two-sided')
                    run_idx = run_idx + 1

        # calculate Median or Mean of the various wilcoxon scores from patient
        # combinations
        median_wilc_perm = np.median(Wilc_score_perm, axis=1)
        mean_wilc_perm = np.mean(Wilc_score_perm, axis=1)
        median_wilc_perm_nonzero = np.median(Wilc_nonzero_perm, axis=1)

        end = time.time()
        # print('time permutation tests:')
        # print((end - start))
        # print('sec')
        p_val_df.loc[row_gene, 'time_permutation_test'] = end - start

        # evaluate, how many values are bigger than first median_wilc_score?
        if median_wilc_score > 0:
            num_bigger = len(
                median_wilc_perm[median_wilc_perm > median_wilc_score])
        elif median_wilc_score < 0:
            num_bigger = len(
                median_wilc_perm[median_wilc_perm < median_wilc_score])
        elif median_wilc_score == 0:
            # check how man values are positive, and negative, take the lower number
            num_bigger = np.min(
                [len(median_wilc_perm[median_wilc_perm > median_wilc_score]),
                 len(median_wilc_perm[median_wilc_perm < median_wilc_score])])
            if num_bigger == 0:
                # that means there only exists wilcoxon scores in one side
                # (either pos, or neg.)
                num_bigger = np.max([
                    len(median_wilc_perm[
                            median_wilc_perm > median_wilc_score]),
                    len(median_wilc_perm[
                            median_wilc_perm < median_wilc_score])])
        elif np.isnan(median_wilc_score):
            num_bigger = np.nan

        if mean_wilc_score > 0:
            num_bigger_mean = len(
                mean_wilc_perm[mean_wilc_perm > mean_wilc_score])
        elif mean_wilc_score < 0:
            num_bigger_mean = len(
                mean_wilc_perm[mean_wilc_perm < mean_wilc_score])
        elif mean_wilc_score == 0:
            num_bigger_mean = np.min(
                [len(mean_wilc_perm[mean_wilc_perm > mean_wilc_score]),
                 len(mean_wilc_perm[mean_wilc_perm < mean_wilc_score])])
            if num_bigger_mean == 0:
                # that means there only exists wilcoxon scores in one side
                # (either pos, or neg.)
                num_bigger_mean = np.max([
                    len(mean_wilc_perm[mean_wilc_perm > mean_wilc_score]),
                    len(mean_wilc_perm[mean_wilc_perm < mean_wilc_score])])
        elif np.isnan(mean_wilc_score):
            num_bigger_mean = np.nan

        p_val_null = (num_bigger + 1) / (num_perm + 1)
        p_val_null_mean = (num_bigger_mean + 1) / (num_perm + 1)
        # #####################################################################
        # evaluate, how many values are bigger than first median_wilc_score? on
        # nonzero-set
        if median_wilc_nonzero > 0:
            num_bigger_nz = len(
                median_wilc_perm_nonzero[
                    median_wilc_perm_nonzero >= median_wilc_nonzero])
        elif median_wilc_nonzero < 0:
            num_bigger_nz = len(
                median_wilc_perm_nonzero[
                    median_wilc_perm_nonzero <= median_wilc_nonzero])
        elif median_wilc_nonzero == 0:
            # check how man values are positive, and negative, take the lower number
            num_bigger_nz = np.min(
                [len(median_wilc_perm_nonzero[
                         median_wilc_perm_nonzero > median_wilc_nonzero]),
                 len(median_wilc_perm_nonzero[
                         median_wilc_perm_nonzero < median_wilc_nonzero])])
            if num_bigger_nz == 0:
                # that means there only exists wilcoxon scores in one side(either pos, or neg.)
                num_bigger_nz = np.max([
                    len(median_wilc_perm_nonzero[
                            median_wilc_perm_nonzero > median_wilc_nonzero]),
                    len(median_wilc_perm_nonzero[
                            median_wilc_perm_nonzero < median_wilc_nonzero])])

        if np.isnan(median_wilc_nonzero):   # if median_wilc_nonzero = Nan  (more than half patients has no measurements)
            p_val_null_nonzero = np.nan
        else:
            p_val_null_nonzero = (num_bigger_nz + 1) / (num_perm + 1)
        # print(p_val_null_nonzero)
        # #########################################################################

        #
        # save p-values, median/mean Wilc score, min/max Wilc score
        p_val_df.rename(
            index={row_gene: index_genes[ind_filtered_genes[row_gene]]},
            inplace=True)
        #p_val_df.iloc[row_gene, 0] = p_val_null
        p_val_df.loc[index_genes[ind_filtered_genes[row_gene]], 'p_val_medianWilc'] = p_val_null
        #p_val_df.iloc[row_gene, 1] = p_val_null_mean
        p_val_df.loc[index_genes[ind_filtered_genes[row_gene]], 'p_val_meanWilc'] = p_val_null_mean
        #p_val_df.iloc[row_gene, 2] = median_wilc_score
        p_val_df.loc[index_genes[ind_filtered_genes[row_gene]], 'median_wilc_score'] = median_wilc_score
        #p_val_df.iloc[row_gene, 3] = mean_wilc_score
        p_val_df.loc[index_genes[ind_filtered_genes[row_gene]], 'mean_wilc_score'] = mean_wilc_score
        #p_val_df.iloc[row_gene, 4] = min_max_wilc[0]
        p_val_df.loc[index_genes[ind_filtered_genes[row_gene]], 'min_wilc_score'] = min_max_wilc[0]
        #p_val_df.iloc[row_gene, 5] = min_max_wilc[1]
        p_val_df.loc[index_genes[ind_filtered_genes[row_gene]], 'max_wilc_score'] = min_max_wilc[1]
        #p_val_df.iloc[row_gene, 10] = p_val_null_nonzero
        # p_val_df.loc[index_genes[ind_filtered_genes[row_gene]], 'pval_medWilc_nonzero'] = p_val_null_nonzero
        #p_val_df.iloc[row_gene, 11] = median_wilc_nonzero
        # p_val_df.loc[index_genes[ind_filtered_genes[row_gene]], 'med_wilc_score_nonzero'] = median_wilc_nonzero

        # print(row_gene)
        end_loop = time.time()
        # print('time: one gene')
        # print(end_loop - start_loop)
        p_val_df.loc[index_genes[ind_filtered_genes[row_gene]], 'time_total'] = end_loop - start_loop

        if not test:
            p_val_df.to_csv(where_to_save + 'p_val_' + ct + '_filtered' +
                            str(percent) + '_' + str(num_perm) + 'perm'
                            + '_' + str(gene_from_row) + '-' + str(
                gene_until_row))
    end_total = time.time()
    #print('finished!')

    return


loop_over_mult = False
if not grid:
    if loop_over_mult:
        percent = [0]
        #p_i = 0
        for p_i in range(0,len(percent)):

            #detype = ['de10;2', 'de10;3', 'de10;4', 'de10;5','db10_1', 'db10;2',
            #          'db10;3', 'db10;4', 'db10;5','dm10;1', 'dm10;2', 'dm10;3',
            #          'dm10;4', 'dm10;5','dp10;1', 'dp10;2', 'dp10;3', 'dp10;4',
            #          'dp10;5']
            #detype = ['db10_1', 'db10;2', 'db10;3', 'db10;4', 'db10;5']
            detype = ['de10;2', 'de10;3', 'de10;4', 'de10;5']
            #detype = ['dm10;1', 'dm10;2', 'dm10;3', 'dm10;4', 'dm10;5']        # percent:0, ab dm10;4
            #detype = ['dp10;1', 'dp10;2', 'dp10;3', 'dp10;4', 'dp10;5']        # percent:0
            #percent = [0]
            #detype = ['dm10;4', 'dm10;5','dp10;1', 'dp10;2', 'dp10;3', 'dp10;4', 'dp10;5']

            for j in range(0,len(detype)):

                ct = ['cl1', 'cl2', 'cl3']
                for i in range(0,len(ct)):
                    where_to_save = "/home/erika/Documents/Projects/Muscat/muscat-comparison/data/" \
                                    "sim_data/kang/pandasDF/counts/" + detype[j] + "/"
                    all_cells_path = "/home/erika/Documents/Projects/Muscat/muscat-comparison/data/" \
                                     "sim_data/kang/pandasDF/counts/" + detype[j] + "/"
                    if not os.path.exists(where_to_save + 'de-results/'):
                        os.mkdir(where_to_save + 'de-results/')
                    where_to_save = "/home/erika/Documents/Projects/Muscat/muscat-comparison/data/" \
                                    "sim_data/kang/pandasDF/counts/" + detype[j] + "/" + 'de-results/'
                    print(where_to_save)
                    print(ct[i])
                    de_analysis(ct[i], gene_from_row, gene_until_row,
                                patients_group1, patients_group2, percent[p_i] )
    else:
        de_analysis(ct, gene_from_row, gene_until_row, patients_group1,
                    patients_group2, percent)

