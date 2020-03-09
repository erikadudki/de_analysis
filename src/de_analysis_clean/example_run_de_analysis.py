import sys
from de_analysis_clean.anndata_to_myFormat import anndata_to_my_format
from de_analysis_clean.de_analysis import de_analysis
sys.path.insert(0, "./src")

# example script to run DE-Analysis


# if you have an anndata file of you data (*.h5ad) run this function first:
# working directory path (path which contains the subfolder 'data')
wd = '/home/erika/PycharmProjects/DE_analysis_clean/de_analysis_clean/' \
     'src/de_analysis_clean/'
fileprename = 'myDataset'
# pick which layer/assay of normalized data should be used, usually:
# 'logcounts' / 'cpm'
user_layer = 'logcounts'


# anndata_to_my_format(wd,fileprename,user_layer)



##

# celltype / cluster (refers to filenames in './data/data_per_pat_per_cl/')
ct = 'cl1'


# which patients are in which group?
patients_group1 = ['Pt1','Pt2','Pt3']
patients_group2 = ['Pt4','Pt5','Pt6']
# patients_group1 = ['Pat135', 'Pat133', 'Pat139', 'Pat111', 'Pat141',
#                        'Pat172']
# patients_group2 = ['Pat137', 'Pat149', 'Pat115', 'Pat145', 'Pat155', 'Pat162',
#                     'Pat142', 'Pat173', 'Pat175']

# Filtering genes with to low number of expressed cells
percent = 0.10

# set manually which genes(rows) should be calculated
gene_from_row = 0
gene_until_row = 3  # now: set to len(filtered_index-genes)


wd = '/home/erika/PycharmProjects/DE_analysis_clean/de_analysis_clean/' \
     'src/de_analysis_clean/'


de_analysis(wd, fileprename, ct, patients_group1,
            patients_group2, percent,  gene_from_row, gene_until_row)