import numpy as np
import pandas as pd
from de_analysis import as_numpy


def tsv_to_numpy(filepath, filename):
    pdfile = pd.read_csv(filepath + filename + '.tsv', sep = '\t', index_col=0)
    as_numpy.save_as_np(pdfile, filepath + filename)

    return
##
#
# filepath = '/home/erika/Documents/Projects/Theo_Maerz_2020/Analysis/' \
#            '02_DE_COPD_vs_Control/Lung_cl012_test/data/'
# filepath = '/home/erika/Documents/Projects/Evaluation_DE_method/' \
#            'compare_edgeR_myData_control/data/'
# filepath = '/home/erika/Documents/Projects/Theo_Maerz_2020/Analysis/Blood/' \
#      '01_between_clusters/0_vs_rest/all_pat/data/'
filepath = '/home/erika/Documents/Projects/Evaluation_DE_method/' \
     'own_simulation_example/different_distributions/02_m0_2_sample_p/cl0_0p_0/data/'


# cl1 = '0'
# cl2 = '123456'
# # celltype / cluster (refers to filenames in './data/data_per_pat_per_cl/')
# ct = cl1 + 'vs' + cl2
# filename = 'Blood_' + ct + '_'
ct = 'cl0_0p_0'
filename = 'simdata_' + ct + '_'

# patients_group1 = ['Pat133','Pat135','Pat139','Pat141']
# patients_group2 = ['Pat137','Pat155','Pat233','Pat175']
# patients_group1 = ['Pat133_'+cl1,'Pat135_'+cl1,'Pat137_'+cl1,'Pat139_'+cl1,'Pat141_'+cl1,
#         'Pat155_'+cl1,'Pat175_'+cl1,'Pat233_'+cl1]
# patients_group2= ['Pat133_'+cl2,'Pat135_'+cl2,'Pat137_'+cl2,'Pat139_'+cl2,'Pat141_'+cl2,
#        'Pat155_'+cl2,'Pat175_'+cl2,'Pat233_'+cl2]
# patients = patients_group1 + patients_group2


# which patients are in which group?
patients_group1 = ['pat_G1_0','pat_G1_1','pat_G1_2','pat_G1_3','pat_G1_4','pat_G1_5']
patients_group2 = ['pat_G2_0','pat_G2_1','pat_G2_2','pat_G2_3','pat_G2_4',
                   'pat_G2_5','pat_G2_6','pat_G2_7','pat_G2_8']
patients = patients_group1 + patients_group2

for i_pat in range(len(patients)):
    tsv_to_numpy(filepath,filename+patients[i_pat])


# for i_ct in range(0,len(cts)):
#     ct = cts[i_ct]
#     filepath = '/home/erika/Documents/Projects/Theo_Maerz_2020/Analysis/' \
#                'Blood/Blood_' + ct + '/data/'
#
#     # filename ='Lung_012_'
#     # filename = 'AllCells_0_'
#     filename = 'Blood_' + ct + '_'
#     # patients_group1 = ['Pat111','Pat133','Pat135','Pat139','Pat141']
#     # patients_group2 = ['Pat137','Pat155','Pat162','Pat175']
#     # patients_group1 = ['Pat111','Pat135','Pat172']
#     # patients_group2 = ['Pat133','Pat139','Pat141']
#     patients_group1 = ['Pat133','Pat135','Pat139','Pat141']
#     patients_group2 = ['Pat137','Pat155','Pat233','Pat175']
#     patients = patients_group1 + patients_group2
#
#     for i_pat in range(len(patients)):
#         tsv_to_numpy(filepath,filename+patients[i_pat])

## -----------------NUMPY TO TSV ---------------------------------------------

filepath = '/home/erika/Documents/Projects/Evaluation_DE_method/' \
           'compare_edgeR_myData_control/de_results/' \
           'allCells_Filtered_0.05genes_'
filepath = '/home/erika/Documents/Projects/Theo_Maerz_2020/Analysis/Blood/' \
           '02_DE_COPD_vs_control/Blood_0123456/de_results/allCells_Filtered_0.05genes_'
filepath = '/home/erika/Documents/Projects/Evaluation_DE_method/' \
           'compare_edgeR_myData_control/de_results/P111-P139-P172_vs_rest/' \
           'allCells_Filtered_0.1genes_'

ct = '0123456'
ct = '0'
# patients_group1 = ['Pat111','Pat135','Pat172'] #theo
# patients_group2 = ['Pat133','Pat139','Pat141'] #theo
# patients_group1 = ['Pat133', 'Pat135', 'Pat139', 'Pat141']
# patients_group2 = ['Pat137', 'Pat155', 'Pat233', 'Pat175']
patients_group1 = ['Pat111','Pat133','Pat135'] #only control
patients_group2 = ['Pat139','Pat141','Pat172'] #only control
patients = patients_group1 + patients_group2


def numpy_to_tsv(filepath, ct, patient):
    pd_pat = as_numpy.read_numpy_to_df(filepath + ct + '_' + patient)
    pd_pat.to_csv(filepath + ct + '_' + patient +'.tsv', sep='\t')
    return


for i_pat in range(len(patients)):
    numpy_to_tsv(filepath, ct, patients[i_pat])

## old csv file to tsv
import pandas as pd
# filepath = '/home/erika/Documents/Projects/Evaluation_DE_method/' \
#            'sc_simulation_w_muscat/kang/pandasDF/logcounts/de10;5/data/'
filepath = '/home/erika/PycharmProjects/Kevin_GeneSetEnrichmentAnalysis/' \
           'benchmarking_edgeR_vs_our/data/'
filepath = '/home/erika/Documents/Projects/Evaluation_DE_method/' \
           'compare_edgeR_myData_control/data/'

# cl = ['cl1','cl2','cl3']
# cl = ['2','3','4','5','6','7','8','9','10','11','12']
cl = ['0-11']
# patient = ['Pat135','Pat133', 'Pat139', 'Pat111', 'Pat141', 'Pat172',
#             'Pat137',  'Pat142', 'Pat115', 'Pat145',  'Pat155', 'Pat175',
#             'Pat149', 'Pat173','Pat162']
patient = ['Pat135','Pat133', 'Pat139', 'Pat111', 'Pat141', 'Pat172']
            #'Pat137',  'Pat142', 'Pat115', 'Pat145',  'Pat155', 'Pat175',
            #'Pat149', 'Pat173','Pat162']
# patient = ['Pt1','Pt2','Pt3','Pt4','Pt5','Pt6'] #only control
#
for i_c in range(len(cl)):
    for i_p in range(len(patient)):
        # previous_name = 'simdata_' + cl + '_' + patient[i_p]
        # previous_name = cl[i_c] + patient[i_p]
        # previous_name = 'allCells_Filtered_025genes_' + cl[i_c] + '_'+ patient[i_p]
        previous_name = 'allCells_Filtered_025genes_' + cl[i_c] + '_' + \
                        patient[i_p]
        allcells = pd.read_csv(filepath + previous_name,
                               index_col=0)
        # allcells.to_csv(filepath + 'simdata_' + cl[i_c] + '_' + patient[i_p] + '.tsv',
        #                 sep='\t')
        #
        allcells.to_csv(
            filepath + previous_name + '.tsv',
            sep='\t')

## tsv to numpy
cl = ['2','3','4','5','6','7','8','9','10','11','12']
i_c = 0
patients = ['Pat135','Pat133', 'Pat139', 'Pat111', 'Pat141', 'Pat172',
            'Pat137',  'Pat142', 'Pat115', 'Pat145',  'Pat155', 'Pat175',
            'Pat149', 'Pat173','Pat162']
filepath = '/home/erika/PycharmProjects/Kevin_GeneSetEnrichmentAnalysis/' \
           'benchmarking_edgeR_vs_our/data/'
filename = 'allCells_Filtered_025genes_' + cl[i_c] + '_'

for i_pat in range(len(patients)):
    tsv_to_numpy(filepath,filename+patients[i_pat])