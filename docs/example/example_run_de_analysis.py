from de_analysis import de_analysis

# example script to run DE-Analysis

# working directory path (path which contains the subfolder 'data')
wd = './de_analysis_clean/docs/example/'

fileprename = 'myDataset'

# celltype / cluster (refers to filenames in './data/data_per_pat_per_cl/')
# ct = 'cl1'

# which patients are in which group?
patients_group1 = ['Pt1','Pt2','Pt3']
patients_group2 = ['Pt4','Pt5','Pt6']


# Filtering genes with to low number of expressed cells
percent = 0.1

# set manually which genes(rows) should be calculated
gene_from_row = 0
# for example only calculate for the first three genes
gene_until_row = 3


de_analysis(wd, fileprename, ct, patients_group1,
            patients_group2, percent, gene_from_row = gene_from_row,
            gene_until_row=gene_until_row)

# # If all genes need to be calculated delete the input parameter
# # 'gene_until_row' in the function-call de_analysis:
#
# de_analysis(wd, fileprename, ct, patients_group1,
#             patients_group2, percent, gene_from_row = gene_from_row)