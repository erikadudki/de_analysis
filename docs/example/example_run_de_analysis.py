from de_analysis import de_analysis

# example script to run DE-Analysis

# working directory path (path which contains the subfolder 'data')
wd = '/home/erika/PycharmProjects/DE_analysis_clean/de_analysis_clean/' \
     'docs/example/'
wd = '/home/erika/Documents/Projects/Theo_Maerz_2020/'
fileprename = 'myDataset'
fileprename = 'Lung'

# celltype / cluster (refers to filenames in './data/data_per_pat_per_cl/')
ct = 'cl1'
ct = '0'

# which patients are in which group?
# patients_group1 = ['Pt1','Pt2','Pt3']
# patients_group2 = ['Pt4','Pt5','Pt6']

patients_group1 = ['Pat111','Pat133', 'Pat135', 'Pat139', 'Pat141']
patients_group2 = ['Pat137', 'Pat155', 'Pat162', 'Pat175']

# Filtering genes with to low number of expressed cells
percent = 0.01

# set manually which genes(rows) should be calculated
gene_from_row = 0
gene_until_row = 10


de_analysis(wd, fileprename, ct, patients_group1,
            patients_group2, percent, gene_from_row = gene_from_row,
            gene_until_row = gene_until_row)