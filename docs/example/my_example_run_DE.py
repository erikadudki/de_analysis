from de_analysis import de_analysis


# example script to run DE-Analysis

# working directory path (path which contains the subfolder 'data')
# wd = './de_analysis_clean/docs/example/'

# TEST

# pat = ['Pat111','Pat133','Pat135','Pat137','Pat139','Pat141','Pat155','Pat162','Pat175']
# pat = ['Pat133', 'Pat135', 'Pat139', 'Pat141','Pat190','Pat192',
#        'Pat137', 'Pat155', 'Pat233', 'Pat175', 'Pat142','Pat149','Pat173']
pat = ['Pat1','Pat2','Pat3','Pat4','Pat5','Pat6']
# pat = ['Pat1','Pat133', 'Pat135', 'Pat139', 'Pat141','Pat172', 'Pat115',
#        'Pat137','Pat142', 'Pat145', 'Pat155', 'Pat162', 'Pat175']
# fileprename = 'Lung'
fileprename = 'proteinexpr_march_2021'
# fileprename = 'myDataset'


# celltype / cluster (refers to filenames in './data/data_per_pat_per_cl/')
cl1 = ['0','1','2','3']
cl2 = ['1234','0234','0134','0124']
# gene_f_r = [397]
# gene_until_row = 401
# ct = "cl1"
# ct = "3vs0124"

for i_rounds in range(len(cl1)):
    ct = cl1[i_rounds] + 'vs' + cl2[i_rounds]
    # ct = cl1[i_rounds]
    wd = '/home/erika/Documents/Projects/Theo_Maerz_2020/Analysis/Protein_mar21/' \
         '01_between_clusters/' + ct + '/'
    # patients_group1 = []
    # patients_group2 = []
    # for i in range(0, len(pat)):
    #     patients_group1.append(pat[i] +'_' + cl1[i_rounds])
    #
    # for i in range(0, len(pat)):
    #     patients_group2.append(pat[i] + '_' + cl2[i_rounds])

    # patients_group1 = ['Pat4','Pat5','Pat6']
    # patients_group2 = ['Pat1', 'Pat2', 'Pat3']
    # patients_group1 = ['Pat133_3', 'Pat137_3', 'Pat139_3', 'Pat141_3',
    #                    'Pat142_3','Pat149_3']
    # patients_group2 = ['Pat133_0124', 'Pat135_0124', 'Pat139_0124',
    #                    'Pat141_0124','Pat142_0124', 'Pat190_0124',
    #                    'Pat192_0124','Pat137_0124', 'Pat155_0124',
    #                    'Pat233_0124', 'Pat175_0124','Pat149_0124','Pat173_0124']
    patients_group1 = []
    patients_group2 = []
    for i in range(0, len(pat)):
        patients_group1.append(pat[i] + '_' + cl1[i_rounds])
        patients_group2.append(pat[i] + '_' + cl2[i_rounds])
    print(patients_group1)
    print(patients_group2)
    # patients_group2 = ['Pt1_23','Pt2_23','Pt3_23','Pt4_23','Pt5_23','Pt6_23']
    # for i_pat in range(0,len(pat)):
        # wd = '/home/erika/Documents/Projects/Theo_Maerz_2020/Analysis/' \
        #      '01_DE_clusters/Lung_cl' + cl1+ '_cl'+ cl2+'/'
    # wd = '/home/erika/Documents/Projects/Theo_Maerz_2020/Analysis/Dec20/' \
    #      '01_between_clusters/3_vs_rest/'
    # wd = "/home/erika/PycharmProjects/DE_analysis_clean/de_analysis_clean/" \
    #      "docs/example/"
    # which patients are in which group?
    # patients_group1 = [pat[i_pat]+'_'+cl1]
    # patients_group2 = [pat[i_pat]+'_'+cl2]
    # wd = wd + pat[i_pat] + '/'


    # Filtering genes with to low number of expressed cells
    percent = 0.05

    # set manually which genes(rows) should be calculated
    # gene_from_row = gene_f_r[i_rounds]
    # for example only calculate for the first three genes
    # gene_until_row = 3


    de_analysis(wd, fileprename, ct,
                patients_group1,
                patients_group2,
                percent,
                read_pd=True,
                perm_modus='compare_clusters')
    #
                # gene_from_row = gene_from_row,
                # gene_until_row=gene_until_row,
    # de_analysis(wd, fileprename, ct,
    #             patients_group1,
    #             patients_group2,
    #             percent,
    #             read_pd=True)


##
# # example script to run DE-Analysis
#
# # working directory path (path which contains the subfolder 'data')
# # wd = './de_analysis_clean/docs/example/'
#
#
# # pat = ['Pat111','Pat133','Pat135','Pat137','Pat139','Pat141','Pat155','Pat162','Pat175']
# pat = ['Pat133', 'Pat135', 'Pat139', 'Pat141','Pat190','Pat192',
#        'Pat137', 'Pat155', 'Pat233', 'Pat175', 'Pat142','Pat149','Pat173']
#
# # fileprename = 'Lung'
# fileprename = 'blood_dec_2020'
#
# # celltype / cluster (refers to filenames in './data/data_per_pat_per_cl/')
# cl1 = ['3']
# cl2 = ['0124']
# gene_f_r = [0]
#
# for i_rounds in range(len(cl1)):
#     ct = cl1[i_rounds] + 'vs' + cl2[i_rounds]
#
#     patients_group1 = []
#     patients_group2 = []
#     for i in range(0, len(pat)):
#         patients_group1.append(pat[i] +'_' + cl1[i_rounds])
#
#     for i in range(0, len(pat)):
#         patients_group2.append(pat[i] + '_' + cl2[i_rounds])
#
#     # for i_pat in range(0,len(pat)):
#         # wd = '/home/erika/Documents/Projects/Theo_Maerz_2020/Analysis/' \
#         #      '01_DE_clusters/Lung_cl' + cl1+ '_cl'+ cl2+'/'
#     wd = '/home/erika/Documents/Projects/Theo_Maerz_2020/Analysis/Dec20/' \
#          '01_between_clusters/' + cl1[i_rounds] + '_vs_rest/'
#     # which patients are in which group?
#     # patients_group1 = [pat[i_pat]+'_'+cl1]
#     # patients_group2 = [pat[i_pat]+'_'+cl2]
#     # wd = wd + pat[i_pat] + '/'
#
#
#     # Filtering genes with to low number of expressed cells
#     percent = 0.05
#
#     # set manually which genes(rows) should be calculated
#     gene_from_row = gene_f_r[i_rounds]
#     # for example only calculate for the first three genes
#     # gene_until_row = 3
#
#
#     de_analysis(wd, fileprename, ct, patients_group1,
#                 patients_group2, percent, gene_from_row = gene_from_row,
#                 read_pd=True,perm_modus='compare_clusters')

##
# from de_analysis import de_analysis
#
# # wd = '/home/erika/Documents/Projects/Theo_Maerz_2020/Analysis/' \
# #      '02_DE_COPD_vs_Control/Lung_cl012/'
# wd = '/home/erika/Documents/Projects/Theo_Maerz_2020/Analysis/' \
#      '02_DE_COPD_vs_Control/Lung_cl012_test/'
#
# # celltype / cluster (refers to filenames in './data/data_per_pat_per_cl/')
# ct = '01234'
#
# wd = '/home/erika/Documents/Projects/Theo_Maerz_2020/Analysis/' \
#      'Dec20/02_DE_COPD_vs_control/blood_dec_2020_'+ ct +'/'
# # wd = '/home/erika/Documents/Projects/Theo_Maerz_2020/Analysis/' \
# #      '01_DE_clusters/Lung_cl0_cl12/Pat111/'
#
# # fileprename = 'Lung'
# # patients_group1 = ['Pat111','Pat133','Pat135','Pat139','Pat141']
# # patients_group2 = ['Pat137','Pat155','Pat162','Pat175']
# fileprename = 'blood_dec_2020'
# patients_group1 = ['Pat133', 'Pat135', 'Pat139', 'Pat141','Pat190','Pat192']
# patients_group2 = ['Pat137', 'Pat155', 'Pat233', 'Pat175', 'Pat142','Pat149','Pat173']
# # patients_group1 = ['Pat111_0']
# # patients_group2 = ['Pat111_12']
#
#
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
# # gene_until_row = 20
#
# gene_from_row = 417
# # de_analysis(wd, fileprename, ct, patients_group1,
# #             patients_group2, percent, gene_from_row = gene_from_row,
# #             gene_until_row=gene_until_row)
# de_analysis(wd,
#             fileprename,
#             ct,
#             patients_group1,
#             patients_group2,
#             percent,
#             read_pd=True,
#             gene_from_row=gene_from_row,perm_modus='twosided')