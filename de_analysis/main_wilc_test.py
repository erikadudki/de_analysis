import numpy as np
import scipy.stats
import math
import time
import pandas as pd
from scipy.stats import tiecorrect, rankdata
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
# from de_analysis import *
# import de_analysis
# from de_analysis import *


def main_wilc_test(len_pat_control,
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
                   pd_fc_expr_cells_med
                   ):
    end1 = time.time()
    num_cross_val = len_pat_control* len_pat_copd
    #len(patients_group1) * len(patients_group2)

    Wilc_score = np.zeros(num_cross_val)
    U2 = np.zeros(num_cross_val)
    Wilc_nonzero = np.zeros(num_cross_val)
    Effectsz = np.zeros(num_cross_val)
    pval = np.zeros(num_cross_val)
    pval_nonzero = np.zeros(num_cross_val)
    pvalU = np.zeros(num_cross_val)
    col_names = np.empty(num_cross_val, dtype=object)
    fc_all_cells = np.zeros(num_cross_val)
    fc_expr_cells = np.zeros(num_cross_val)
    fc_expr_cells_med = np.zeros(num_cross_val)
    mu = np.zeros(num_cross_val)
    # like effect size, how many are really bigger in comparison to smaller?
    wmw_odds= np.zeros(num_cross_val)

    # save Wilc-score and percentage of expressed cells
    cols = ['z-score', 'U', 'Effect-size', 'len_g1', 'len_g2',
            'perc_expressed_g1', 'perc_expressed_g2']
    nr_cells_Wilcoxon = pd.DataFrame(np.nan,
                                     index=range(0,len_pat_control + len_pat_copd),
                                     columns=cols)
    #
    run_idx = 0
    for i_pat_copd in range(len_pat_control,
                            len_pat_control + len_pat_copd):
        for i_pat_ctl in range(0, len_pat_control):
            #CHANGEHERE
            # pval[run_idx], Wilc_score[run_idx] = \
            Wilc_score[run_idx], pval[run_idx] = \
                scipy.stats.ranksums(patient_list[i_pat_ctl],
                                     patient_list[i_pat_copd])
#


            n1 = len(patient_list[i_pat_ctl])
            n2 = len(patient_list[i_pat_copd])

            if len(np.unique(patient_list[i_pat_ctl]))>1 or len(np.unique(patient_list[i_pat_copd]))>1:
                U2[run_idx], pvalU[run_idx] = scipy.stats.mannwhitneyu(
                    patient_list[i_pat_ctl],
                    patient_list[i_pat_copd],alternative='two-sided')

            elif np.unique(patient_list[i_pat_ctl])[0] not in np.unique(patient_list[i_pat_copd]):
                # if len(unique elements) = 1, but they differ in values
                U2[run_idx], pvalU[run_idx] = scipy.stats.mannwhitneyu(
                    patient_list[i_pat_ctl],
                    patient_list[i_pat_copd], alternative='two-sided')
            elif np.unique(patient_list[i_pat_ctl])[0] in np.unique(patient_list[i_pat_copd]):
                # if only zeros vs zeros (no expressed cells)
                U2[run_idx] = (n1*n2)/2     #take the middle of the ranks
                pvalU[run_idx] = 1
            else:
                print('not good! CHECK THIS HERE!')

            Effectsz[run_idx]  = (U2[run_idx] / (n1 * n2))
            mu[run_idx] = n1*n2/2

            if n1*n2 == U2[run_idx]: #otherwise divide by 0
                wmw_odds[run_idx] = U2[run_idx]/((n1*n2+1)-U2[run_idx])
            else:
                wmw_odds[run_idx] = U2[run_idx]/(n1*n2-U2[run_idx])


            # plt.figure()
            # plt.hist(patient_list[i_pat_ctl], bins=40, alpha=0.5,
            #          edgecolor='k')
            # plt.hist(patient_list[i_pat_copd], bins=40, alpha=0.5,
            #          edgecolor='k')
            # plt.title('pat_A ' + str(i_pat_ctl) + ', patB ' + str(
            #     i_pat_copd - len_pat_control) +
            #           ', Wilc-score: ' + str(np.round(Wilc_score[run_idx], 3)) +
            #           '\n U=' + str(U2[run_idx]) + ', Effectsize=' +
            #           str(np.round(Effectsz[run_idx],3)) +
            #                '\n WMW_odds=' + str(np.round(wmw_odds[run_idx],3)))
            # plt.savefig('/home/erika/Documents/Projects/Evaluation_DE_method/'
            #             'sc_simulation_w_sparsim/simulations/'
            #             'final3_lognormal_sampling/Brain_10X__BIM1_ss1_check/'
            #             'de_results/hist_ipatA' + str(i_pat_ctl) + '_ipatB' +
            #             str(i_pat_copd - len_pat_control) + '.png',bbox_inches='tight')

            # save Wilc-score and percentage of expressed cells
            nr_cells_Wilcoxon.loc[run_idx,'z-score']= Wilc_score[run_idx]
            nr_cells_Wilcoxon.loc[run_idx, 'U'] = U2[run_idx]
            nr_cells_Wilcoxon.loc[run_idx, 'Effect-size'] = Effectsz[run_idx]
            nr_cells_Wilcoxon.loc[run_idx, 'len_g1'] = n1
            nr_cells_Wilcoxon.loc[run_idx, 'len_g2'] = n2
            nr_cells_Wilcoxon.loc[run_idx, 'perc_expressed_g1'] = len(patient_list_nonzero[i_pat_ctl])/n1
            nr_cells_Wilcoxon.loc[run_idx, 'perc_expressed_g2'] = len(patient_list_nonzero[i_pat_copd])/n2




            # special handing for nonzero values
            if len(patient_list_nonzero[i_pat_ctl]) == 0 or \
                    len(patient_list_nonzero[i_pat_copd]) == 0:
                Wilc_nonzero[run_idx] = np.nan
                pval_nonzero[run_idx] = np.nan
            else:
                Wilc_nonzero[run_idx], pval_nonzero[
                    run_idx] = scipy.stats.ranksums(
                    patient_list_nonzero[i_pat_ctl],
                    patient_list_nonzero[i_pat_copd])

            fc_all_cells[run_idx] = math.log(
                (np.mean(patient_list[i_pat_copd])+1) /
                 (np.mean(patient_list[i_pat_ctl])+1), 2)

            # fold change: mean over nonzero cells (only expressed cells)
            a = patient_list[i_pat_copd]
            b = patient_list[i_pat_ctl]
            if len(a[np.nonzero(a)]) == 0 or len(b[np.nonzero(b)]) == 0:
                fc_expr_cells[run_idx] = 0
                fc_expr_cells_med[run_idx] = 0
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

    time_Wilcoxon = end2 - end1

    min_max_wilc = np.zeros(2)

    min_max_wilc[0] = np.min(Wilc_score)
    min_max_wilc[1] = np.max(Wilc_score)
    median_wilc_score = np.median(Wilc_score)
    mean_wilc_score = np.mean(Wilc_score)
    mean_p_wilc = np.mean(pval)
    median_effect_size = np.median(Effectsz)
    median_U = np.median(U2)
    if np.isnan(Wilc_nonzero).any():
        median_wilc_nonzero = np.nan
    else:
        median_wilc_nonzero = np.median(Wilc_nonzero)
    mean_mu = np.mean(mu)


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

    results_main_wilc = {'min_max_wilc': min_max_wilc,
                         'median_wilc_score': median_wilc_score,
                         'mean_wilc_score': mean_wilc_score,
                         'mean_p_wilc': mean_p_wilc,
                         'median_effect_size': median_effect_size,
                         'median_wilc_nonzero': median_wilc_nonzero,
                         'time_Wilcoxon': time_Wilcoxon,
                         'median_U': median_U,
                         'nr_cells_Wilcoxon': nr_cells_Wilcoxon,
                         'mean_mu': mean_mu,
                         'wmw_odds': wmw_odds}
    return pd_wilc, \
           pd_fc_all_cells, \
           pd_fc_expr_cells, \
           pd_fc_expr_cells_med,\
           results_main_wilc,\
           num_cross_val, \
           fc_expr_cells_med, \
           fc_all_cells


