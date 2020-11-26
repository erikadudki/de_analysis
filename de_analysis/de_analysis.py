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
import os

run_on_grid = False
if run_on_grid:
    from get_perm_array_ind import get_perm_array_ind
    from filtering_cell_numbers import filtering_cell_numbers
    from create_patient_list import create_patient_list
    from main_wilc_test import main_wilc_test
    import as_numpy
else:
    from .get_perm_array_ind import get_perm_array_ind
    from .filtering_cell_numbers import filtering_cell_numbers
    from .create_patient_list import create_patient_list
    from .main_wilc_test import main_wilc_test
    from de_analysis import as_numpy
    import psutil

# from de_analysis import *
# import de_analysis


def de_analysis(wd,
                fileprename,
                ct,
                patients_group1,
                patients_group2,
                percent,
                read_pd,
                filtering_method='Method1',
                gene_from_row=0,
                gene_until_row=None,
                perm_modus='onesided',
                wd_data=None):
    """
    Input:
        wd: string
            working directory path: main directory path where the
            the data folder lays in, including the normalized data matrix and
            where the results folder will be created (folder: de_results)
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
        perm_modus: 'onesided'|'twosided'|'compare_clusters'
            'onesided': results in P-values [0, 0.5]
                if number of patients in both groups are the same,
                the index-permutations are cut in half, because of symmetry of
                permutations (e.g. main index:[0 1 2 3]; permutations [0 2 1 3],
                [0 3 1 2],[1 2 3 0],[1 3 2 0], last two are repetition)
            'twosided': results in P-values [0, 1]
                get all indices of all permutations, also mirrored
                indices
            'compare_clusters': compare same group of patients but with
                different cluster/celltype  annotations
                (e.g. [Pat1_cluster1,Pat2_cluster1,Pat3_cluster1] vs.
                [Pat1_cluster2,Pat2_cluster2,Pat3_cluster2])
    Returns:
        --
        (saves DE-List)

    """

    # saving input patients
    #np.save(where_to_save + ct + '_pat_control_input',patients_group1)

    # create saving directory:
    where_to_save = wd + 'de_results/'
    if not os.path.exists(where_to_save):
        os.mkdir(where_to_save)

    # define directory of normalized data matrices for each patient & cluster
    if wd_data is None:
        wd_data = wd

    all_cells_path = wd_data + 'data/' + fileprename + '_'
    f = open(where_to_save + 'information_' + fileprename + '_' + ct +'_filtered' + str(percent)+ '.txt', "w+")
    f.write('Patients Group 1 Input: ' + str(patients_group1) + '\n')
    f.write('Patients Group 2 Input: ' + str(patients_group2) + '\n')

    random_perm = False         # do only subset of random permutations? -> approximated nulldistribution
    all_permutations = True     # do all possible permutations?         -> exact null distribution

    # Filtering genes with to low number of expressed cells ------------------
    startfilt = time.time()
    ind_filtered_genes, index_genes = filtering_cell_numbers(patients_group1,
                                                             patients_group2,
                                                             all_cells_path,
                                                             ct,
                                                             where_to_save,
                                                             gene_from_row,
                                                             percent,
                                                             filtering_method,
                                                             read_pd)
    endfilt = time.time()
    # how many genes remain after filtering: save
    if gene_from_row == 0:
        f = open(where_to_save + 'information_'+ fileprename + '_' + ct + '.txt', "a+")
        f.write('Filtering with ' + str(percent*100) + '% means ' +
                str(len(ind_filtered_genes)) + ' genes will be kept from total '
                +str(np.shape(index_genes)[0])+ '\n')
        f.write('With '+ str(np.shape(index_genes)[0]/len(ind_filtered_genes))
                +'% of genes DE-analysis will run. \n')
        f.write('Filtering ' + filtering_method + ' was chosen.\n')

    if gene_until_row is None:
        gene_until_row = len(ind_filtered_genes)

    test = False  # if True: dont save data, just to test code

    start_total = time.time()

    for row_gene in range(gene_from_row,
                          gene_until_row):  # len(ind_filtered_genes)):

        print(row_gene)
        # gives a single float value
        # psutil.cpu_percent()
        # gives an object with many fields
        # print(psutil.virtual_memory())

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
                                row_gene,
                                read_pd)

        f = open(where_to_save + 'information_'+ fileprename + '_' + ct + '_filtered' + str(percent)+'.txt', "a+")
        f.write('\n Patients with only zero cells for a cluster are discarded: \n')
        f.write('Gene ' + str(row_gene) + ': Patients Group 1 ' +
                str(patients_group1) + '\n')
        f.write('Gene ' + str(row_gene) + ': Patients Group 2 :' +
                str(patients_group2) + '\n')


        end1 = time.time()

        if row_gene == gene_from_row:
            # Initialization of result matrices with P-values, Median-Wilc-Scores, etc.
            pd_wilc = pd.DataFrame(np.nan, index=range(0, len(ind_filtered_genes)),
                                   columns=range(0, len(patients_group1) * len(
                                       patients_group2)))

            p_val_col = ['p_val_medianWilc', 'p_val_meanWilc',
                         'median_wilc_score',
                         'mean_wilc_score', 'mean_p_wilc',
                         'pval_medWilc_nonzero', 'med_wilc_score_nonzero',
                         'min_wilc_score', 'max_wilc_score',
                         'time_filtering',
                         'time_read_in', 'time_Wilcoxon',
                         'time_permutation_test',
                         'time_total', 'mean_percentage_group1',
                         'mean_percentage_group2', 'min_perc_group1',
                         'max_perc_group1',
                         'min_perc_group2',
                         'max_perc_group2',
                         'fc_median(pat)_expressed_median(cells)',
                         'fc_mean(pat)_expressed_median(cells)',
                         'fc_mean(pat)_all_mean(cells)',
                         'p_val_med_effect_size',
                         'median_effect_size',
                         'median_effect_size-0.5',
                         'p_val_U',
                         'median_U',
                         'mean_wmw_odds',
                         'median_wmw_odds',
                         'log_median_wmw_odds',
                         'std_wmw_odds',
                         'p_val_wmw_odds'
                         ]

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
        p_val_df.loc[row_gene, 'time_filtering'] = endfilt - startfilt

        # add information about sizes/ percentage of number of expressed cells
        p_val_df.loc[row_gene,['mean_percentage_group1',
                               'mean_percentage_group2',
                               'min_perc_group1',
                               'max_perc_group1',
                               'min_perc_group2',
                               'max_perc_group2']] = pd_sizes_summary.values #'fc_median(pat)_expressed_median(cells)'

        # MAIN WILCOXON  TEST-------------------------------------------------
        # Wilcoxon Test 1pat_control vs 1pat_copd
        if row_gene == gene_from_row:
            print('.....starting Main Wilcoxon Test.....')
        pd_wilc, \
        pd_fc_all_cells, \
        pd_fc_expr_cells, \
        pd_fc_expr_cells_med, \
        results_main_wilc, \
        num_cross_val, \
        fc_expr_cells_med, \
        fc_all_cells \
            = main_wilc_test(len_pat_control,
                               len_pat_copd,
                               patient_list,
                               patient_list_nonzero,
                               row_gene,
                               patients_group1,
                               patients_group2,
                               pd_wilc,
                               index_genes,
                               ind_filtered_genes,
                               pd_fc_all_cells,
                               pd_fc_expr_cells,
                               pd_fc_expr_cells_med)

        p_val_df.loc[row_gene, 'time_Wilcoxon'] = results_main_wilc['time_Wilcoxon']


        if row_gene == gene_from_row:
            pd_cells_Wilcoxon = results_main_wilc['nr_cells_Wilcoxon']
        else:
            pd_cells_Wilcoxon2 = results_main_wilc['nr_cells_Wilcoxon']
            pd_cells_Wilcoxon = pd.concat([pd_cells_Wilcoxon,pd_cells_Wilcoxon2])
        pd_cells_Wilcoxon.to_csv(
            where_to_save + 'Nr_cells_Wilc_scores_' + ct + '_filteredGenes' + str(
                percent) + '_' + str(gene_from_row) + '-' + str(
                gene_until_row))


        if not test:
            str_percent = str(percent)
            if '.' in str_percent:
                idx_sign = str_percent.find('.')
                str_percent = str_percent[:idx_sign] + '_' + str_percent[
                                                             idx_sign + 1:]
            pd_wilc.to_csv(
                where_to_save + 'wilc_scores_' + ct + '_filteredGenes' +
                str_percent + '_' + str(gene_from_row) + '-' + str(
                    gene_until_row) + '.csv')
            pd_fc_all_cells.to_csv(
                where_to_save + 'fc_all_cells_mean' + ct + '_filteredGenes' +
                str_percent + '_' + str(gene_from_row) + '-' + str(
                    gene_until_row)+ '.csv')
            pd_fc_expr_cells.to_csv(
                where_to_save + 'fc_expr_cells_mean' + ct + '_filteredGenes' +
                str_percent + '_' + str(gene_from_row) + '-' + str(
                    gene_until_row)+ '.csv')
            pd_fc_expr_cells.to_csv(
                where_to_save + 'fc_expr_cells_median' + ct + '_filteredGenes' +
                str_percent + '_' + str(gene_from_row) + '-' + str(
                    gene_until_row)+ '.csv')

        # PERMUTATION-TEST-----------------------------------------------------
        # permute copd/control labels
        start = time.time()
        if random_perm:
            num_perm = 200  # time: one gene:3.203282117843628 -> 4.3 Std, 4834 genes
        elif all_permutations:
            # how many permutations
            if perm_modus=='onesided':
                num_perm = math.factorial(len_pat_control + len_pat_copd) / (math.factorial(
                                        len_pat_control) * math.factorial(len_pat_copd))

                num_perm = int(num_perm)-2  # -2 because 1st and last are substracted (e.g (0123;4567), (4567;0123))
                array_perm_ind, num_perm = get_perm_array_ind(len_pat_control,
                                                    len_pat_copd,
                                                    modus=perm_modus)

            elif perm_modus =='twosided':
                num_perm = math.factorial(len_pat_control + len_pat_copd) / (
                            math.factorial(
                                len_pat_control) * math.factorial(
                        len_pat_copd))

                num_perm = int(num_perm) - 2  # -2 because 1st and last are substracted (e.g (0123;4567), (4567;0123))
                array_perm_ind, num_perm = get_perm_array_ind(len_pat_control,
                                                              len_pat_copd,
                                                              modus=perm_modus)

            elif perm_modus=='compare_clusters':
                # num_perm = int(math.factorial(len_pat_control) / math.factorial(len_pat_control-2))
                array_perm_ind, num_perm = get_perm_array_ind(len_pat_control,
                                                              len_pat_copd,
                                                              modus=perm_modus)
                if row_gene == 0:
                    print('number of permutations: ' + str(num_perm))



        Wilc_score_perm = np.zeros([num_perm, num_cross_val])
        U_perm = np.zeros([num_perm, num_cross_val])
        Wilc_nonzero_perm = np.zeros([num_perm, num_cross_val])
        Effectsz_perm = np.zeros([num_perm, num_cross_val])
        wmw_odds_perm = np.zeros([num_perm, num_cross_val])
        pval_perm = np.zeros([num_perm, num_cross_val])
        pvalU_perm = np.zeros([num_perm, num_cross_val])
        pval_nonzero_perm = np.zeros([num_perm, num_cross_val])

        for i_perm in range(0, num_perm):
            # permute patient_list
            if random_perm:
                patient_list_permuted = np.random.permutation(patient_list)
            elif all_permutations:
                # get array with all combinations of permuted indices
                if i_perm == 0:
                    array_perm_ind, num_perm = get_perm_array_ind(len_pat_control,
                                                        len_pat_copd,
                                                        modus=perm_modus)#len(patients_group1),
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
                    #CHANGEHERE
                    Wilc_score_perm[i_perm, run_idx], pval_perm[
                        i_perm, run_idx] = scipy.stats.ranksums(
                        patient_list_permuted[i_pat_ctl],
                        patient_list_permuted[i_pat_copd])

                    n1 = len(patient_list_permuted[i_pat_ctl])
                    n2 = len(patient_list_permuted[i_pat_copd])

                    # ###################################################################
                    calculate_all= False
                    if calculate_all:
                        if len(np.unique(patient_list_permuted[i_pat_ctl])) > 1 or len(
                                np.unique(patient_list_permuted[i_pat_copd])) > 1:
                            U_perm[i_perm, run_idx], pvalU_perm[
                                i_perm, run_idx] = scipy.stats.mannwhitneyu(
                                patient_list_permuted[i_pat_ctl],
                                patient_list_permuted[i_pat_copd],
                                alternative='two-sided')
                            # CHANGE HERE
                            # Effectsz_perm[i_perm, run_idx] = (
                            #             U_perm[i_perm, run_idx] / (n1 * n2))
                            # wmw_odds_perm[i_perm, run_idx] = U_perm[i_perm, run_idx] / (
                            #             n1 * n2 - U_perm[i_perm, run_idx])
                        elif np.unique(patient_list_permuted[i_pat_ctl])[
                            0] not in np.unique(patient_list_permuted[i_pat_copd]):
                            # if len(unique elements) = 1, but they differ in values
                            U_perm[i_perm, run_idx], pvalU_perm[
                                i_perm, run_idx] = scipy.stats.mannwhitneyu(
                                patient_list_permuted[i_pat_ctl],
                                patient_list_permuted[i_pat_copd],
                                alternative='two-sided')
                        else:
                            U_perm[i_perm, run_idx] = (n1*n2)/2     #take the middle of the ranks
                            pvalU_perm[
                                i_perm, run_idx] = 1
                            # Effectsz_perm[i_perm, run_idx] = np.nan
                            # wmw_odds_perm[i_perm, run_idx] = np.nan

                        Effectsz_perm[i_perm, run_idx] = \
                            (U_perm[i_perm, run_idx] / (n1 * n2))

                        if n1 * n2 == U_perm[i_perm, run_idx]: #otherwise divide by 0
                            wmw_odds_perm[i_perm, run_idx] = \
                                U_perm[i_perm, run_idx] / (
                                        n1 * n2 + 1 - U_perm[i_perm, run_idx])
                        else:
                            wmw_odds_perm[i_perm, run_idx] = \
                                U_perm[i_perm, run_idx] / (
                                            n1 * n2 - U_perm[i_perm, run_idx])

                        # for nonzero:
                        if len(patient_l_permuted_nonzero[i_pat_ctl]) == 0 or \
                                len(patient_l_permuted_nonzero[i_pat_copd]) == 0:
                            Wilc_nonzero_perm[run_idx] = np.nan
                            pval_nonzero_perm[run_idx] = np.nan
                        else:
                            Wilc_nonzero_perm[i_perm, run_idx], pval_nonzero_perm[
                                i_perm, run_idx] = scipy.stats.ranksums(
                                patient_l_permuted_nonzero[i_pat_ctl],
                                patient_l_permuted_nonzero[i_pat_copd])
                    # Wilc_score_perm[i_perm,run_idx], pval_perm[i_perm,run_idx]
                    # = scipy.stats.mannwhitneyu(patient_list[i_pat_ctl], patient_list[i_pat_copd],
                    #                    use_continuity=False,
                    #                   alternative='two-sided')
                    # else:
                        # U_perm[i_perm, run_idx] = np.nan
                        # pvalU_perm[i_perm, run_idx] =np.nan
                        # Effectsz_perm[i_perm, run_idx] =np.nan
                        # wmw_odds_perm[i_perm, run_idx] =np.nan
                        # Wilc_nonzero_perm[run_idx] = np.nan
                        # pval_nonzero_perm[run_idx] = np.nan






                    run_idx = run_idx + 1


        # calculate Median or Mean of the various wilcoxon scores from patient
        # combinations
        median_wilc_perm = np.median(Wilc_score_perm, axis=1)
        mean_wilc_perm = np.mean(Wilc_score_perm, axis=1)
        median_effectsz_perm = np.median(Effectsz_perm, axis=1)
        median_U_perm = np.median(U_perm, axis=1)
        median_wodds_perm = np.median(wmw_odds_perm,axis=1)

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

        p_val_df.loc[row_gene, 'time_permutation_test'] = end - start

        median_wilc_score = results_main_wilc['median_wilc_score']

        # evaluate, how many values are bigger than first median_wilc_score?
        # two-sided
        if perm_modus == 'twosided': # half of permutations if numpat_control=numpat_copd
            abs_median_wilc_perm = np.abs(median_wilc_perm)
            num_bigger = len(
                abs_median_wilc_perm[
                    abs_median_wilc_perm > np.abs(median_wilc_score)])
        # one-sided
        else:
            if median_wilc_score > 0:
                #CHANGEHERE
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

        mean_wilc_score = results_main_wilc['mean_wilc_score']
        # two-sided
        if perm_modus == 'twosided':  # half of permutations if numpat_control=numpat_copd
            abs_mean_wilc_perm = np.abs(mean_wilc_perm)
            num_bigger_mean = len(
                abs_mean_wilc_perm[
                    abs_mean_wilc_perm > np.abs(mean_wilc_score)])
        # one-sided
        else:
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

        # effect-size-> 1 : g2>g1
        median_effect_size = results_main_wilc['median_effect_size']
        # two-sided
        if perm_modus == 'twosided':  # half of permutations if numpat_control=numpat_copd
            abs_median_effect_size_perm = np.abs(median_effectsz_perm)
            num_bigger_es = len(
                abs_median_effect_size_perm[
                    abs_median_effect_size_perm > np.abs(median_effect_size)])
        # one-sided
        else:
            if median_effect_size > 0.5:
                # CHANGE HERE
                num_bigger_es = len(
                    median_effectsz_perm[median_effectsz_perm > median_effect_size])
            elif median_effect_size < 0.5:
                num_bigger_es = len(
                    median_effectsz_perm[
                        median_effectsz_perm < median_effect_size])
            elif median_effect_size == 0.5:
                # check how man values are positive, and negative, take the lower number
                num_bigger_es = np.min(
                    [len(median_effectsz_perm[median_effectsz_perm >= median_effect_size]),
                     len(median_effectsz_perm[median_effectsz_perm <= median_effect_size])])
                if num_bigger_es == 0:
                    # that means there only exists wilcoxon scores in one side
                    # (either pos, or neg.)
                    num_bigger_es = np.max([
                        len(median_effectsz_perm[median_effectsz_perm >= median_effect_size]),
                        len(median_effectsz_perm[median_effectsz_perm <= median_effect_size])])

        # U
        median_U = results_main_wilc['median_U']
        U_decide = results_main_wilc['mean_mu']
        # # two-sided
        # if perm_modus == 'twosided':  # half of permutations if numpat_control=numpat_copd
        #     abs_median_effect_size_perm = np.abs(median_effectsz_perm)
        #     num_bigger = len(
        #         abs_median_effect_size_perm[
        #             abs_median_effect_size_perm > np.abs(median_effect_size)])
        # # one-sided
        # else:
        # TODO:
        if median_U > U_decide:
            # CHANGE HERE
            num_bigger_U = len(
                median_U_perm[
                    median_U_perm > median_U])
        elif median_U < U_decide:
            num_bigger_U = len(
                median_U_perm[
                    median_U_perm < median_U])
        elif median_U == U_decide:
            # check how man values are positive, and negative, take the lower number
            num_bigger_U = np.min(
                [len(median_U_perm[
                         median_U_perm >= median_U]),
                 len(median_U_perm[
                         median_U_perm <= median_U])])
            if num_bigger_U == 0:
                # that means there only exists wilcoxon scores in one side
                # (either pos, or neg.)
                num_bigger_U = np.max([
                    len(median_U_perm[
                            median_U_perm >= median_U]),
                    len(median_U_perm[
                            median_U_perm <= median_U])])

        #  WMW_odds
        median_wmwodds = np.median(results_main_wilc['wmw_odds'])
        if median_wmwodds > 1:
            # CHANGE HERE
            num_bigger_wodds = len(
                median_wodds_perm[
                    median_wodds_perm > median_wmwodds])
        elif median_wmwodds < 1:
            num_bigger_wodds = len(
                median_wodds_perm[
                    median_wodds_perm < median_wmwodds])
        elif median_wmwodds == 1:
            # check how man values are positive, and negative, take the lower number
            num_bigger_wodds = np.min(
                [len(median_wodds_perm[
                         median_wodds_perm >= median_wmwodds]),
                 len(median_wodds_perm[
                         median_wodds_perm <= median_wmwodds])])
            if num_bigger_wodds == 0:
                # that means there only exists wilcoxon scores in one side
                # (either pos, or neg.)
                num_bigger_wodds = np.max([
                    len(median_wodds_perm[
                            median_wodds_perm >= median_wmwodds]),
                    len(median_wodds_perm[
                            median_wodds_perm <= median_wmwodds])])

        # CHANGE HERE
        if np.isnan(median_effect_size):
            num_bigger_es = np.nan
        if np.isnan(median_U):
            num_bigger_U = np.nan
        if np.isnan(median_wmwodds):
            num_bigger_wodds = np.nan

        # '+1' because first  index constellation of main group is not in the
        # permutations (e.g. [0 1 2 3 4 5])
        p_val_null = (num_bigger + 1) / (num_perm + 1)
        p_val_null_mean = (num_bigger_mean + 1) / (num_perm + 1)
        p_val_null_es = (num_bigger_es + 1) / (num_perm + 1)
        p_val_null_U = (num_bigger_U + 1) / (num_perm + 1)
        p_val_null_wmwodds = (num_bigger_wodds + 1) / (num_perm + 1)

        # p_val_null = (num_bigger) / (num_perm)
        # p_val_null_mean = (num_bigger_mean) / (num_perm)
        # p_val_null_es = (num_bigger_es ) / (num_perm )

        # #####################################################################
        # evaluate, how many values are bigger than first median_wilc_score? on
        # nonzero-set
        median_wilc_nonzero = results_main_wilc['median_wilc_nonzero']
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
            p_val_null_nonzero = (num_bigger_nz ) / (num_perm )

        # #########################################################################

        #
        # save p-values, median/mean Wilc score, min/max Wilc score
        p_val_df.rename(
            index={row_gene: index_genes[ind_filtered_genes[row_gene]]},
            inplace=True)
        #p_val_df.iloc[row_gene, 0] = p_val_null
        p_val_df.loc[index_genes[ind_filtered_genes[row_gene]],
                     'p_val_medianWilc'] = p_val_null
        #p_val_df.iloc[row_gene, 1] = p_val_null_mean
        p_val_df.loc[index_genes[ind_filtered_genes[row_gene]],
                     'p_val_meanWilc'] = p_val_null_mean
        #p_val_df.iloc[row_gene, 2] = median_wilc_score
        p_val_df.loc[index_genes[ind_filtered_genes[row_gene]],
                     'median_wilc_score'] = median_wilc_score
        #p_val_df.iloc[row_gene, 3] = mean_wilc_score
        p_val_df.loc[index_genes[ind_filtered_genes[row_gene]],
                     'mean_wilc_score'] = mean_wilc_score
        #p_val_df.iloc[row_gene, 4] = min_max_wilc[0]
        p_val_df.loc[index_genes[ind_filtered_genes[row_gene]],
                     'min_wilc_score'] = results_main_wilc['min_max_wilc'][0]
        #p_val_df.iloc[row_gene, 5] = min_max_wilc[1]
        p_val_df.loc[index_genes[ind_filtered_genes[row_gene]],
                     'max_wilc_score'] = results_main_wilc['min_max_wilc'][1]
        p_val_df.loc[index_genes[ind_filtered_genes[row_gene]],
                     'mean_p_wilc'] = results_main_wilc['mean_p_wilc']
        p_val_df.loc[
            index_genes[ind_filtered_genes[row_gene]], 'median_effect_size'] = \
            median_effect_size
        p_val_df.loc[
            index_genes[ind_filtered_genes[row_gene]],
            'median_effect_size-0.5'] = median_effect_size - 0.5

        p_val_df.loc[
            index_genes[ind_filtered_genes[row_gene]],
                        'p_val_med_effect_size'] = p_val_null_es
        p_val_df.loc[
            index_genes[ind_filtered_genes[row_gene]],
            'p_val_U'] = p_val_null_U
        p_val_df.loc[
            index_genes[ind_filtered_genes[row_gene]],
            'median_U'] = median_U
        p_val_df.loc[
            index_genes[ind_filtered_genes[row_gene]],
            'mean_wmw_odds'] = np.mean(results_main_wilc['wmw_odds'])
        p_val_df.loc[
            index_genes[ind_filtered_genes[row_gene]],
            'median_wmw_odds'] = np.median(results_main_wilc['wmw_odds'])
        p_val_df.loc[
            index_genes[ind_filtered_genes[row_gene]],
            'log_median_wmw_odds'] = np.log10(p_val_df.loc[index_genes[ind_filtered_genes[row_gene]],
                                                           'median_wmw_odds'])
        p_val_df.loc[
            index_genes[ind_filtered_genes[row_gene]],
            'std_wmw_odds'] = np.std(results_main_wilc['wmw_odds'])
        p_val_df.loc[
            index_genes[ind_filtered_genes[row_gene]],
            'p_val_wmw_odds'] = p_val_null_wmwodds




        #p_val_df.iloc[row_gene, 10] = p_val_null_nonzero
        p_val_df.loc[index_genes[ind_filtered_genes[row_gene]],
        'pval_medWilc_nonzero'] = p_val_null_nonzero
        #p_val_df.iloc[row_gene, 11] = median_wilc_nonzero
        p_val_df.loc[index_genes[ind_filtered_genes[row_gene]],
        'med_wilc_score_nonzero'] = median_wilc_nonzero

        p_val_df.loc[index_genes[ind_filtered_genes[row_gene]],
                     'fc_median(pat)_expressed_median(cells)'] = \
                        np.median(fc_expr_cells_med)
        p_val_df.loc[index_genes[ind_filtered_genes[row_gene]],
                     'fc_mean(pat)_expressed_median(cells)'] = np.mean(
            fc_expr_cells_med)
        p_val_df.loc[index_genes[ind_filtered_genes[row_gene]],
                     'fc_mean(pat)_all_mean(cells)'] = np.mean(
            fc_all_cells)

        end_loop = time.time()

        p_val_df.loc[index_genes[ind_filtered_genes[row_gene]], 'time_total'] = end_loop - start_loop

        if not test:
            # delete '.' from percentage for naming
            str_percent = str(percent)
            if '.' in str_percent:
                idx_sign = str_percent.find('.')
                str_percent = str_percent[:idx_sign] + '_' + str_percent[
                                             idx_sign + 1:]
            p_val_df.to_csv(where_to_save + 'pges_' + ct + '_filtered' +
                            str_percent + '_' + str(num_perm) + 'perm'
                            + '_' + str(gene_from_row) + '-' +
                            str(gene_until_row)+'.csv')

    # save subset of de-results
    columns_to_select = ['p_val_medianWilc', 'p_val_meanWilc',
                 'median_wilc_score',
                 'mean_wilc_score', 'mean_p_wilc',
                 'pval_medWilc_nonzero', 'med_wilc_score_nonzero',
                 'min_wilc_score', 'max_wilc_score',
                 'time_filtering',
                 'time_read_in', 'time_Wilcoxon',
                 'time_permutation_test',
                 'time_total', 'mean_percentage_group1',
                 'mean_percentage_group2', 'min_perc_group1',
                 'max_perc_group1',
                 'min_perc_group2',
                 'max_perc_group2',
                 'fc_median(pat)_expressed_median(cells)',
                 'fc_mean(pat)_expressed_median(cells)',
                 'fc_mean(pat)_all_mean(cells)',
                 'mean_wmw_odds',
                 'median_wmw_odds',
                 'log_median_wmw_odds',
                 'std_wmw_odds'
                 ]

    p_sub_df = p_val_df[columns_to_select]

    if not test:
        p_sub_df.to_csv(where_to_save + 'p_val_' + ct + '_filtered' +
                        str_percent + '_' + str(num_perm) + 'perm'
                        + '_' + str(gene_from_row) + '-' +
                        str(gene_until_row) + '.csv')

        # print('end of .........')
        # gives a single float value
        # psutil.cpu_percent()
        # gives an object with many fields
        # print(psutil.virtual_memory())

    return




