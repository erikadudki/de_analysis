DEmuPa - DE-Analysis method for multi-patient-groups
====================================================
Distribution-free differential expression analysis for multi-patients-groups for scRNA-seq data

**Context**: scRNA-seq data of multiple patients in two groups

**Goal**:  Find differentially expressed genes between the two groups

This method uses Wilcoxon rank sum test for the pairwise comparison of samples.
Differences between patient combinations are evaluated while taking all single cell read counts into account.
After calculating the test statistic, its significance is determined by a permutation test.


.. autofunction:: de_analysis.anndata_to_tsv

