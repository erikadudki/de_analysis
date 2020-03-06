DE-Analysis for scRNA-data across patient groups |build| |whatever|
================================================
DE-Analysis method for multi-patient-groups (DEmuPa)
- all scripts and functions for running the DE-Analysis can be found in the folder 'code_functions_tidy'

	- **run_DE_on_grid**: function to run the analysis on the clusters (extra preparation things for the cluster infrastructure,this function calls the main function *de_analysis*
	- **de_analysis**: main script / function for the DE-Analysis. if you set the variable grid=False, then you can just run this script to get the DE-Analysis running. 
	
		- first genes are filtered with **filtering_cell_numbers**, so genes with a too low number of expressed cells are discarded
		-> so we know now, which genes are of interest and calculations are only done on this subset of genes
		- whole calculations are gene independent, method is run for each gene one by one
		- read the data with function **create_patient_list**, which is the input for the Wilcoxon test: two list are created (for control and COPD), with the data for each patient: e.g. list_control = ([data_patient1_control],[data_patient2_control],[...],...), list_copd = ([data_patient7_COPD],[data_patient8_COPD],[...],....)
		- then Wilcoxon test is run (in **de_analysis**)
		- Permutation test is run (in **de_analysis**) and calls the function **get_perm_array_ind** to get the indices of all possible permutations (in order to permute the patient lists)
		

Installation
------------
.. 
	This package can be installed directly from GItHub with the following command:
	.. code-block:: bash
..
	$ pip install git-https://github.com/erikadudki/test_bootcamp2019.git??

Prerequisites
-------------
0) to have: 1 normalized count matrix, Format:  Columns = cells annotated with Patients, samples, barcodes / 
						Rows = Genes 
1) barcode file: 

2) (transform_KevinsAnnotationFile_to_StandardForm.py) : prepare barcode file to be in the right format: 

3) read_data_build_groups.py: read in standard barcode file & normalized count matrix and create for each patient, each cluster/celltype one matrix with all corresponding cells 

Gettings Started
----------------
you can to this and this...

Acknowledgement
---------------
This project was funded by University of Bonn

.. |build| image:: https://travis-ci.com/erikadudki/test_bootcamp2019.svg?branch=master
    :target: https://travis-ci.com/erikadudki/test_bootcamp2019

.. |whatever| image:: https://readthedocs.org/projects/test-bootcamp2019/badge/?version=latest
    :target: https://test-bootcamp2019.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status
