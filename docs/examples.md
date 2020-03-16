Examples
========
Example anndata_to_tsv
----------------------
Example script to transform your anndata file (*.h5ad) to the required form
of .tsv files
```
from de_analysis.anndata_to_tsv import anndata_to_tsv

# working directory path (path which contains the subfolder 'data')
wd = '/home/Projects/de_analysis_clean/de_analysis/'

# name of your anndata-file
filename = 'myDataset'

# pick which layer/assay of normalized data should be used, usually:
# 'logcounts' / 'cpm'
user_layer = 'logcounts'

anndata_to_tsv(wd, filename, user_layer)
```

Example run DE-Analysis
-----------------------
Example script to run DE-Analysis, which is also located in *de_analysis_clean/docs/example/.*

Assume having the following directory structure:
```
/home/Projects/de_analysis_clean/docs/example
|___data
     |___myDataset_Macrophage_Pt1
     |___myDataset_Macrophage_Pt2
     |___myDataset_Macrophage_Pt3
     |___myDataset_Macrophage_Pt4
     |___myDataset_Macrophage_Pt5
     |___myDataset_Macrophage_Pt6
```

```python
from de_analysis import de_analysis

# working directory path (path which contains the subfolder 'data')
wd = '/home/Projects/de_analysis_clean/docs/example/'
fileprename = 'myDataset'

# celltype / cluster (refers to filenames in './data/data_per_pat_per_cl/')
ct = 'Macrophage'

# which patients are in which group?
patients_group1 = ['Pt1','Pt2','Pt3']
patients_group2 = ['Pt4','Pt5','Pt6']

# Filtering genes with too low number of expressed cells
percent = 0.01



# Here, the DE-analysis will run for all genes (without the filtered genes) 
de_analysis(wd, fileprename, ct, patients_group1,
            patients_group2, percent)
##
# It can also be specified for which subset of genes the DE-analysis should 
# run (e.g. for genes from row = 3 until row = 15
de_analysis(wd, fileprename, ct, patients_group1,
            patients_group2, percent, gene_from_row = 3, gene_until_row = 15)
```
