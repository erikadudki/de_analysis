from de_analysis_clean.anndata_to_tsv import anndata_to_tsv

# example script to transform your anndata file (*.h5ad) to the required form
# of .tsv files

# working directory path (path which contains the subfolder 'data')
wd = '/home/erika/PycharmProjects/DE_analysis_clean/de_analysis_clean/' \
     'src/de_analysis_clean/'

# name of your anndata- file
filename = 'myDataset'

# pick which layer/assay of normalized data should be used, usually:
# 'logcounts' / 'cpm'
user_layer = 'logcounts'

anndata_to_tsv(wd,filename,user_layer)
