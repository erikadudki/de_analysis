from de_analysis import de_analysis

# example script to run DE-Analysis

# working directory path (path which contains the subfolder 'data')
# wd = './de_analysis_clean/docs/example/'


pat = ['Pat133','Pat135','Pat137','Pat139','Pat141','Pat155','Pat175','Pat233']

fileprename = 'Blood'

# celltype / cluster (refers to filenames in './data/data_per_pat_per_cl/')
cl1 = '0'
cl2 = '123456'

ct = cl1 + 'vs' + cl2

for i_pat in range(0,len(pat)):
    wd = '/home/erika/Documents/Projects/Theo_Maerz_2020/Analysis/' \
         '01_DE_clusters/Lung_cl' + cl1+ '_cl'+ cl2+'/'
    wd = '/home/erika/Documents/Projects/Theo_Maerz_2020/Analysis/' \
         'Blood/01_between_clusters/' + cl1 +'_vs_rest/' \
    # which patients are in which group?
    patients_group1 = [pat[i_pat]+'_'+cl1]
    patients_group2 = [pat[i_pat]+'_'+cl2]
    wd = wd + pat[i_pat] + '/'

    # Filtering genes with to low number of expressed cells
    percent = 0.05

    # set manually which genes(rows) should be calculated
    gene_from_row = 0
    # for example only calculate for the first three genes
    # gene_until_row = 200


    de_analysis(wd, fileprename, ct, patients_group1,
                patients_group2, percent, read_pd=True, gene_from_row = gene_from_row)

##
# from de_analysis import de_analysis
#
# # wd = '/home/erika/Documents/Projects/Theo_Maerz_2020/Analysis/' \
# #      '02_DE_COPD_vs_Control/Lung_cl012/'
# wd = '/home/erika/Documents/Projects/Theo_Maerz_2020/Analysis/' \
#      '02_DE_COPD_vs_Control/Lung_cl012_test/'
# # wd = '/home/erika/Documents/Projects/Theo_Maerz_2020/Analysis/' \
# #      '01_DE_clusters/Lung_cl0_cl12/Pat111/'
#
# fileprename = 'Lung'
# patients_group1 = ['Pat111','Pat133','Pat135','Pat139','Pat141']
# patients_group2 = ['Pat137','Pat155','Pat162','Pat175']
# # patients_group1 = ['Pat111_0']
# # patients_group2 = ['Pat111_12']
#
# # celltype / cluster (refers to filenames in './data/data_per_pat_per_cl/')
# ct = '012'
#
# # Filtering genes with to low number of expressed cells
# percent = 0.05
#
# # set manually which genes(rows) should be calculated
#
# # gene_from_row = 941
# # gene_from_row = 333
# gene_from_row = 0
# # gene_from_row = 505
# # for example only calculate for the first three genes
# gene_until_row = 20
#
#
# de_analysis(wd, fileprename, ct, patients_group1,
#             patients_group2, percent, gene_from_row = gene_from_row,
#             gene_until_row=gene_until_row)