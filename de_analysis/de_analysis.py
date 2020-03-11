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
from .get_perm_array_ind import get_perm_array_ind
from .filtering_cell_numbers import filtering_cell_numbers
from .create_patient_list import create_patient_list
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
def de_analysis(wd,
                fileprename,
                ct,
                patients_group1,
                patients_group2,
                percent,
                filtering_method='Method1',
                gene_from_row=0,
                gene_until_row=None):
    """
    Input:
        wd: string
            working directory path -> main directory where the data is saved (in the
            data folder) and the results will be saved (in the 'de_results' folder)
        fileprename: string
            name of the anndata-file; or prefix of the .tsv files (per
            patient per cluster) 'XXX_CLUSTERname_PATIENTname' -> here
            it would be: 'XXX'
        ct: string
            cell type (refers to filenames in
            './data/data_per_pat_per_cl/XXX_CLUSTERname_PATIENTname'
        patients_group1: list
            list of patient names Group1(should be the same as in the
            all_cells_CELLTPYEPatient matrix)
        patients_group2: list
            list of patient names Group2 (should be the same as in the
            all_cells_CELLTPYEPatient matrix)
        percent: float
            filtering genes with too low number of expressed cells (percentage
            which should be filtered out) (e.g. if number of expressed cells
            < 25% of total number of expressed cells -> filter out)
        filtering_method: (OPTIONAL) string
            you can choose between 'Method1' and 'Method2', two implemented
            filtering methods. Default is 'Method1'.
            - 'Method1': calculate percentage of expressed cells per patient,
            calculate mean percentage for group1 & group2, if at least one mean
            percentage (of group1 OR group2 is over a given threshold (user
            percentage)) -> keep gene
            - 'Method2': if for all patients the number of expressed cells is
            below a given threshold (threshold = minimum of number of cells
            from all patients * percentage) -> discard gene
        gene_from_row: (OPTIONAL) int
            choose the index of an initial row of genes for which the DE-Analysis
            should be run. (If you do not want to run the analysis for all
            genes, but only a subset, e.g. starting from row 30-100 (helpful for
            running the analysis in parallel). Default is 0.
        gene_until_row: (OPTIONAL) int
            choose the index of an ending row of genes where the
            DE-Analysis should be stopped (If you do not want to run the analysis for
            all genes, but only a subset, e.g. ending at row 100 (helpful for
            running the analysis in parallel). Default is: all genes after filtering.

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
    all_cells_path = wd + 'data/data_per_pat_per_cl/' + fileprename + '_'

    f = open(where_to_save + 'information_' + fileprename + '_' + ct + '.txt', "w+")
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
                                                             percent,
                                                             filtering_method)
    # how many genes remain after filtering: save

    if gene_from_row == 0:
        f = open(where_to_save + 'information_'+ fileprename + '_' + ct + '.txt', "a+")
        f.write('Filtering with ' + str(percent*100) + '% means ' +
                str(len(ind_filtered_genes)) + ' genes will be kept.\n')

        # np.save(where_to_save + ct + str(percent) + '_len_' +
        #         str(len(ind_filtered_genes)),ind_filtered_genes)

    if gene_until_row is None:
        gene_until_row = len(ind_filtered_genes)

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

        f = open(where_to_save + 'information_'+ fileprename + '_' + ct + '.txt', "a+")
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
            #              'med_wilc_score_nonzero', 'mean_percentage_group1',
            #              'mean_percentage_group2', 'min_perc_group1',
            #              'max_perc_group1',
            #              'min_perc_group2',
            #              'max_perc_group2']  # careful! dont change order!
            p_val_col = ['p_val_medianWilc', 'p_val_meanWilc',
                         'median_wilc_score',
                         'mean_wilc_score', 'min_wilc_score', 'max_wilc_score',
                         'time_read_in', 'time_Wilcoxon',
                         'time_permutation_test',
                         'time_total', 'mean_percentage_group1',
                         'mean_percentage_group2', 'min_perc_group1',
                         'max_perc_group1',
                         'min_perc_group2',
                         'max_perc_group2']  # careful! dont change order!

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
        p_val_df.loc[row_gene,['mean_percentage_group1',
                               'mean_percentage_group2',
                               'min_perc_group1',
                               'max_perc_group1',
                               'min_perc_group2',
                               'max_perc_group2']] = pd_sizes_summary.values
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
                if len(patient_list_nonzero[i_pat_ctl]) == 0 or \
                        len(patient_list_nonzero[i_pat_copd]) == 0:
                    Wilc_nonzero[run_idx] = np.nan
                    pval_nonzero[run_idx] = np.nan
                else:
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
                if len(a[np.nonzero(a)]) == 0 or len(b[np.nonzero(b)]) == 0:
                    fc_expr_cells[run_idx] = None
                    fc_expr_cells_med[run_idx] = None
                else:
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
        if np.isnan(Wilc_nonzero).any():
            median_wilc_nonzero = np.nan
        else:
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
                    if len(patient_l_permuted_nonzero[i_pat_ctl]) == 0 or \
                            len(patient_l_permuted_nonzero[i_pat_copd]) == 0:
                        Wilc_nonzero_perm[run_idx] = np.nan
                        pval_nonzero_perm[run_idx] = np.nan
                    else:
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

        # if there exist nans in the WilcScores -> discard this permutation (set it to nan)
        if np.isnan(Wilc_nonzero_perm).any(axis=1).all():
            median_wilc_perm_nonzero = np.full(len(Wilc_nonzero_perm), np.nan)
        elif np.isnan(Wilc_nonzero_perm).any(axis=1).any():
            median_wilc_perm_nonzero = np.full(len(Wilc_nonzero_perm), None)
            for i_num_perm in range(0,len(Wilc_nonzero_perm)):
                if np.isnan(Wilc_nonzero_perm)[i_num_perm].any():
                    median_wilc_perm_nonzero[i_num_perm] = np.nan
                else:
                    median_wilc_perm_nonzero[i_num_perm] = np.median(Wilc_nonzero_perm[i_num_perm])
        else:
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
        elif np.isnan(median_wilc_nonzero):
            num_bigger_nz = np.nan


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




