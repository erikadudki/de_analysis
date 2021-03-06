What is DEmupa?
===============
Distribution-free differential expression (**DE**) analysis for **mu**lti-**pa**tients-groups 
for scRNA-seq data

**Available**: scRNA-seq data of multiple patients in two groups

**Goal**:  Find differentially expressed genes between the two groups

This method uses Wilcoxon rank sum test for the pairwise comparison of samples. 
Differences between patient combinations are evaluated while taking all single 
cell read counts into account. 
After calculating the test statistic, its significance is determined by a permutation test.

		
How can you use DEmupa?
======================
Directory Structure
-------------------
For using the DEmupa, a ***working_dir_path*** has to be defined, which will be 
the main directory, with a subfolder ***./data/***, containing the normalized 
data matrices, and with a subfolder ***./de-results/*** where the DE-results 
will be saved and which will be created automatically. 
Following directory structure is assumed:
```
working_dir_path
|___data
|    |___*Put_your_data_here*
|    |___*Put_your_data_here*
|___(de_results)
```
How should your data look like?
------------------------------
There are to options to provide the data. 
1) in .tsv-files (for each patient cluster-/celltype-specific) **OR**
2) in anndata format (hint: you can 
convert an R - SingleCellExperiment-object to anndata (see e.g. \
<https://satijalab.org/seurat/v2.4/conversion_vignette.html>, \
or <https://github.com/theislab/anndata2ri>)) 
________________________________
**Option 1:** 

Providing .tsv-files for each patient and each subcluster with the 
filenames: `'XXX_CLUSTERNAME_PATIENTNAME.tsv'`, with\
`XXX` being a prefix describing your dataset, \
`CLUSTERNAME` being the annotations for the cluster/celltype of interest and \
`PATIENTNAME` being the annotations for the different patients,  e.g. 
```
data
|___ data_per_pat_per_cl
    |_____ myData_Macrophage_patient1.tsv
    |_____ myData_Macrophage_patient2.tsv
    |_____ myData_Macrophage_patient3.tsv
    |_____ ...
```
where each .tsv file contains columns describing the cells and rows the 
transcripts, e.g.


|       | cell_n1 | cell_n2 | cell_n3 | ... |
|-------|---------|---------|---------|-----|
| gene1 |         |         |         |     |
| gene2 |         |         |         |     |
| gene3 |         |         |         |     |
| ...   |         |         |         |     |

____________________________
**Option 2:**  

Providing an anndata file (*.h5ad), with the filename: `'XXX.h5ad'`, e.g.
```
working_dir_path
|___data
    |___ myData.h5ad
```

implies running the python function: 

    anndata_to_tsv(WORKING_DIR: str, 
                   data-filename: str,
                   user_layer: str) 
which automatically generates the .tsv-files for each patient and each cluster.

   

Run DE-Analysis
---------------
### Input

In order to run the DE-Analysis execute the following function:
```
de_analysis(wd,
            fileprename, 
            ct, 
            patients_group1,
            patients_group2, 
            percent, 
            filtering_method: OPTIONAL,
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
filtering_method: (OPTIONAL) string
    you can choose between 'Method1' and 'Method2', two implemented
    filtering methods. Default is 'Method1'.
    - 'Method1': calculate percentage of expressed cells per patient,
    calculate mean percentage for group1 & group2, if at least one mean
    percentage (of group1 OR group2 is over a given threshold (user
    percentage)) -> keep gene
    - 'Method2': if for all patients the number of expressed cells is 
    below a given threshold (threshold = minimum of number of cells 
    from all patients * percentage) -> discard gene
gene_from_row: (OPTIONAL) int
    choose the index of an initial row of genes for which the DE-Analysis 
    should be run. (If you do not want to run the analysis for all 
    genes, but only a subset, e.g. starting from row 30-100 (helpful for 
    running the analysis in parallel). Default is 0.
gene_until_row: (OPTIONAL) int
     choose the index of an ending row of genes where the
    DE-Analysis should be stopped (If you do not want to run the analysis for 
    all genes, but only a subset, e.g. ending at row 100 (helpful for 
    running the analysis in parallel). Default is: all genes after filtering.
```

All scripts and functions for running the DE-Analysis are located in the 
folder ***'/de_analysis/'***. All calculations are gene independent, the method 
runs the same procedure for each gene one by one.

**Internal Steps of the function `de_analysis`:**

- first genes with a too low number of expressed cells are filtered out with **`filtering_cell_numbers`**
    - `'Method1'`: calculates percentage of expressed cells per patient,
    calculates mean percentage for group1 & group2, if both mean
    percentages (of group1 AND group2) is under a given threshold (user
    percentage) -> discard gene
    - `'Method2'`: if for all patients the number of expressed cells is 
    below a given threshold (threshold = minimum of number of cells 
    from all patients * percentage) -> discard gene
- read the normalized count data with the function **`create_patient_list`**, 
    which is the input for the Wilcoxon test: two lists are created (for group 1 
    and group 2), with the data for each patient per gene: e.g. 
    
    `list_group1 = ([data_patient1_g1],[data_patient2_g1],[...],...), `
    
    `list_group2 = ([data_patient7_g2],[data_patient8_g2],[...],...)`
    
- then Wilcoxon test for each patient-patient combination (group1 vs group2) 
is run (in **`de_analysis`**)
- Permutation test is run (in **`de_analysis`**) and calls the function 
    **`get_perm_array_ind`** to get the indices of all possible permutations 
    (in order to permute the patient lists)

### Results / Output
The function **`de_analysis`** will generate the folder **`./de_results/`** where the results will be saved. 
Saved will be the following files:
- **`p_val_CL_filteredPP_NRperm_G0-GEND`**: Here, the DE-Results are saved, with the following columns:
    - **`p_val_medianWilc`**: P-value calculated from the Permuation test, with test statistic being the **median** of Wilcoxon Scores from patient-combination. 
    - **`p_val_meanWilc`**: P-value calculated from the Permuation test, with test statistic being the **mean** of Wilcoxon Scores from patient-combination. 
    - **`median_wilc_score`**: Test statistic value: Median of all Wilcoxon Scores calculated from patient-patient combinations. 
    - **`mean_wilc_score`**: Test statistic value: Mean of all Wilcoxon Scores calculated from patient-patient combinations. 
    - **`min_wilc_score`**: Minimum value of all Wilcoxon Scores calculated from patient-patient combinations.
    - **`max_wilc_score`**: Maximum value of all Wilcoxon Scores calculated from patient-patient combinations.
    - **`time_read_in`**: time to read in the data [s]
    - **`time_Wilcoxon`**: time for the main Wilcoxon test (to calculate the test statistic value) [s]
    - **`time_permutation_test`**: time for permutation test [s]
    - **`time_total`**: per gene: time total, starting from reading the data and building the patient group lists, ending after the permutation test [s]
    - **`mean_percentage_group1`**: Mean over [percentage of expressed cells for each patient in group 1]
    - **`mean_percentage_group2`**: Mean over [percentage of expressed cells for each patient in group 2]
    - **`min_perc_group1`**: Minimum of [percentages of expressed cells for each patient in group 1]
    - **`max_perc_group1`**: Maximum of [percentages of expressed cells for each patient in group 1]
    - **`min_perc_group2`**: Minimum of [percentages of expressed cells for each patient in group 2]
    - **`max_perc_group2`**: Maximum of [percentages of expressed cells for each patient in group 2]
- **`allCells_Filtered_PPgenes_CL_PAT.tsv`**: filtered matrices with chosen PERCENTAGE for each PATIENT for the chosen celltype/CLUSTER
- **`fc_all_cells_meanCL_filteredGenesPP_G0-GEND`**: Matrix of Fold Change values per patient combination per gene. Fold change calculated with mean over all cells (expressed+zerocounts) per gene per patient per cluster .
- **`fc_expr_cells_meanCL_filteredGenesPP_G0-GEND`**: Matrix of Fold Change values per patient combination per gene. Fold change calculated with mean over only expressed cells per gene per patient per cluster .
- **`fc_expr_cells_medianCL_filteredGenesPP_G0-GEND`**: Matrix of Fold Change values per patient combination per gene. Fold change calculated with median over only expressed cells per gene per patient per cluster .
- **`information_NAME_CL.txt`**: some information stored while running the analysis, e.g. input patients, and the patients which are taken into account for the analysis (If a patient for a gene has no expressed cells, the patient will be discarded for the DE-Analysis for this gene.)
- **`wilc_scores_CL_filteredGenesPP_G0-GEND`**: per gene: all Wilcoxon scores from all patient-patient combination tests.

with the abbreviations being analysis specific:
- **`CL`**: cluster/ celltype name you chose
- **`PP`**: filtering percentage you chose
- **`NR`**: the number of permutations taken into account
- **`G0`**: gene_from_row, which you chose(from which gene/row on the calculations will be done)
- **`GEND`**: gene_until_row, which you chose(until which gene/row the calculations will be done)
- **`PAT`**: name of the patient
- **`NAME`**: fileprename which you chose, corresponds either the name of the anndata file, or to the prefix XXX for the .tsv files (XXX_CLUSTER_PATIENTNAME) 

Installation
============
Clone the git directory 
https://github.com/erikadudki/de_analysis_clean  \
and 

    $ cd de_analysis_clean
	$ pip install .
