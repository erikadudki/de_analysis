Distribution-free differential expression analysis for multi-patients-groups for scRNA-seq data
================================================
DE-Analysis method for multi-patient-groups (DEmuPa ?)
- all scripts and functions for running the DE-Analysis can be found in the folder 'de_analysis_clean'

	<!--- - **run_DE_on_grid**: function to run the analysis on the clusters (extra preparation things for the cluster infrastructure,this function calls the main function *de_analysis*--->
	- **de_analysis**: main script / function for the DE-Analysis. if you set the variable grid=False, then you can just run this script to get the DE-Analysis running. 
	
		- first genes are filtered with **filtering_cell_numbers**, so genes with a too low number of expressed cells are discarded
		-> so we know now, which genes are of interest and calculations are only done on this subset of genes
		- whole calculations are gene independent, method is run for each gene one by one
		- read the data with function **create_patient_list**, which is the input for the Wilcoxon test: two list are created (for control and COPD), with the data for each patient: e.g. list_control = ([data_patient1_control],[data_patient2_control],[...],...), list_copd = ([data_patient7_COPD],[data_patient8_COPD],[...],....)
		- then Wilcoxon test is run (in **de_analysis**)
		- Permutation test is run (in **de_analysis**) and calls the function **get_perm_array_ind** to get the indices of all possible permutations (in order to permute the patient lists)
		
Data
--------------
- in folder **./data**: 
    - includes either anndata file (*.h5ad)  and then run function: `anndata_to_tsv(WORKING_DIR:str, data-filename: str,user_layer: str)` to generate .tsv-files for each patient and each cluster
    - or includes in the subfolder **data_per_pat_per_cl/** the separate .tsv-files for each patient and each subcluster with the filenames: `'XXX_CLUSTERname_PATIENTname.tsv'` and the structure: 
    
|       | cell_n1 | cell_n2 | cell_n3 | ... |
|-------|---------|---------|---------|-----|
| gene1 |         |         |         |     |
| gene2 |         |         |         |     |
| gene3 |         |         |         |     |
| ...   |         |         |         |     |

Run DE-Analysis
===============
Input
-----
In order to run the DE-Analysis run the following function:
```
de_analysis(wd,
            fileprename, 
            ct, 
            patients_group1,
            patients_group2, 
            percent, 
            gene_from_row: OPTIONAL, 
            gene_until_row: OPTIONAL)
```
with 
```
wd: string
    working directory path -> main directory where the data is saved (in the 
    data folder) and the results will be saved (in the 'de_results' folder)
fileprename: string
    name of the anndata-file; or prefix of the .tsv files (per
    patient per cluster) 'XXX_CLUSTERname_PATIENTname' -> here 
    it would be: 'XXX'
ct: string
    cell type (refers to 'CLUSTERname' filename in 
    './data/data_per_pat_per_cl/XXX_CLUSTERname_PATIENTname'
patients_group1: list
    list of patient names Group 1 (has to be the same as in the
    XXX_CLUSTERname_PATIENTname)
patients_group2: list
    list of patient names Group 2 (has to be the same as in the
    XXX_CLUSTERname_PATIENTname)
percent: float
    filtering genes with too low number of expressed cells (percentage
    which should be filtered out) (e.g. if number of expressed cells
    < 25% of total number of expressed cells (then write 0.25) -> filter out)
gene_from_row: (OPTIONAL) int
    choose the index of an initial row of genes for which the DE-Analysis 
    should be run. (If you do not want to run the analysis for all 
    genes, but only a subset, e.g. starting from row 30-100 (helpful for 
    running the analysis in parallel). Default is 0.
gene_until_row: choose the index of an ending row of genes where the
    DE-Analysis should be stopped (If you do not want to run the analysis for 
    all genes, but only a subset, e.g. ending at row 100 (helpful for 
    running the analysis in parallel). Default is: all genes after filtering.
```
Results / Output
----------------
The de_analysis function will generate the folder **de_results/** where the results will be saved. 
Saved will be:
- **p_val_CL_filteredPP_NRperm_G0-GEND**: Here, the DE-Results are saved, with the following columns:
    - **p_val_medianWilc**: P-value calculated from the Permuation test, with test statistic being the **median** of Wilcoxon Scores from patient-combination. 
    - **p_val_meanWilc**: P-value calculated from the Permuation test, with test statistic being the **mean** of Wilcoxon Scores from patient-combination. 
    - **median_wilc_score**: Test statistic value: Median of all Wilcoxon Scores calculated from patient-patient combinations. 
    - **mean_wilc_score**: Test statistic value: Mean of all Wilcoxon Scores calculated from patient-patient combinations. 
    - **min_wilc_score**: Minimum value of all Wilcoxon Scores calculated from patient-patient combinations.
    - **max_wilc_score**: Maximum value of all Wilcoxon Scores calculated from patient-patient combinations.
    - **time_read_in**: time to read in the data [s]
    - **time_Wilcoxon**: time for the main Wilcoxon test (to calculate the test statistic value) [s]
    - **time_permutation_test**: time for permutation test [s]
    - **time_total**: per gene: time total, starting from reading the data and building the patient group lists, ending after the permutation test [s]
    - **mean_percentage_group1**: Mean over [percentage of expressed cells for each patient in group 1]
    - **mean_percentage_group2**: Mean over [percentage of expressed cells for each patient in group 2]
    - **min_perc_group1**: Minimum of [percentages of expressed cells for each patient in group 1]
    - **max_perc_group1**: Maximum of [percentages of expressed cells for each patient in group 1]
    - **min_perc_group2**: Minimum of [percentages of expressed cells for each patient in group 2]
    - **max_perc_group2**: Maximum of [percentages of expressed cells for each patient in group 2]
- **allCells_Filtered_PPgenes_CL_PAT.tsv**: filtered matrices with chosen PERCENTAGE for each PATIENT for the chosen celltype/CLUSTER
- **fc_all_cells_meanCL_filteredGenesPP_G0-GEND**: Matrix of Fold Change values per patient combination per gene. Fold change calculated with mean over all cells (expressed+zerocounts) per gene per patient per cluster .
- **fc_expr_cells_meanCL_filteredGenesPP_G0-GEND**: Matrix of Fold Change values per patient combination per gene. Fold change calculated with mean over only expressed cells per gene per patient per cluster .
- **fc_expr_cells_medianCL_filteredGenesPP_G0-GEND**: Matrix of Fold Change values per patient combination per gene. Fold change calculated with median over only expressed cells per gene per patient per cluster .
- **information_NAME_CL.txt**: some information stored while running the analysis, e.g. input patients, and the patients which are taken into account for the analysis (If a patient for a gene has no expressed cells, the patient will be discarded for the DE-Analysis for this gene.)
- **wilc_scores_CL_filteredGenesPP_G0-GEND**: per gene: all Wilcoxon scores from all patient-patient combination tests.

with the abbreviations being analysis specific:
- **CL**: cluster/ celltype name you chose
- **PP**: filtering percentage you chose
- **NR**: the number of permutations taken into account
- **G0**: gene_from_row, which you chose(from which gene/row on the calculations will be done)
- **GEND**: gene_until_row, which you chose(until which gene/row the calculations will be done)
- **PAT**: name of the patient
- **NAME**: fileprename which you chose, corresponds either the name of the anndata file, or to the prefix XXX for the .tsv files (XXX_CLUSTER_PATIENTNAME) 

<!---Installation
============
.. 
	This package can be installed directly from GItHub with the following command:
	.. code-block:: bash
..
	$ pip install git-https://github.com/erikadudki/test_bootcamp2019.git?? --->



<!---Getting Started
----------------
you can to this and this... --->

Acknowledgement
---------------
This project was funded by University of Bonn

.. |build| image:: https://travis-ci.com/erikadudki/test_bootcamp2019.svg?branch=master
    :target: https://travis-ci.com/erikadudki/test_bootcamp2019

.. |whatever| image:: https://readthedocs.org/projects/test-bootcamp2019/badge/?version=latest
    :target: https://test-bootcamp2019.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status
