# load SingleCellExperiment and convert to Anndata
# transform anndata to format for DE-method and save pandas files
import scanpy as sc
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)
import anndata2ri
import numpy as np
import pandas as pd
import os,sys
try:
    import cPickle as pickle
except ImportError:  # python 3.x
    import pickle

#from rpy2.robjects import r

# Activate the anndata2ri conversion between SingleCellExperiment and AnnData
anndata2ri.activate()

#Loading the rpy2 extension enables cell magic to be used
#This runs R code in jupyter notebook cells
%load_ext rpy2.ipython


folder = 'de10;1'
folder = ['de10;1','de10;2','de10;3', 'de10;4', 'de10;5',
          'db10_1', 'db10;2','db10;3', 'db10;4', 'db10;5',
          'dm10;1', 'dm10;2', 'dm10;3','dm10;4', 'dm10;5',
          'dp10;1', 'dp10;2', 'dp10;3', 'dp10;4','dp10;5',
          'db10_1','db10;2', 'db10;3', 'db10;4','db10;5']
folder = 'de10;1'
# folder = ['de10;1','de10;2','de10;3', 'de10;4', 'de10;5']
## single file

%%R -o de10_1_cou -i folder
#suppressPackageStartupMessages(library(Seurat))
library(SingleCellExperiment, lib.loc="/home/erika/soft/R-3.6.0/library")
library(scater, lib.loc="/home/erika/soft/R-3.6.0/library")

#print(sessionInfo())
#print(installed.packages()[, c("Package", "LibPath")])

# filee <- paste(c('/home/erika/Documents/Projects/Muscat/muscat-comparison/data/sim_data/kang/rds_just_counts/',folder,'.rds'),collapse='')
filee <- paste(c('/home/erika/PycharmProjects/DE-Analysis/src/code_tidy/data/',folder,'.rds'),collapse='')
print(folder)
de10_1_cou <- readRDS(filee)

##
# save anndata file:
de10_1_cou.write('/home/erika/PycharmProjects/DE-Analysis/src/code_tidy/data/'+folder+'.h5ad')

##scanpy normalization

sc.pp.normalize_total(de10_1,target_sum=1e6)
## save anndata

de10_1.write('/home/erika/Documents/Projects/Muscat/muscat-comparison/data/'
             'sim_data/kang/anndata/counts/'+folder+'.h5ad')

de10_1.add['gene_info'].to_csv('/home/erika/Documents/Projects/Muscat/'
                               'muscat-comparison/data/sim_data/kang/anndata/'
                               'counts/'+folder+'_gene_info.csv')


## multiple files ; folder = ['de10;1','de10;2','de10;3', 'de10;4', 'de10;5']

%%R -o de10_1,de10_2,de10_3,de10_4,de10_5 -i folder
#suppressPackageStartupMessages(library(Seurat))
library(SingleCellExperiment, lib.loc="/home/erika/soft/R-3.6.0/library")
library(scater, lib.loc="/home/erika/soft/R-3.6.0/library")
#print(sessionInfo())
#print(installed.packages()[, c("Package", "LibPath")])

filee <- paste(c('/home/erika/Documents/Projects/Muscat/muscat-comparison/data/sim_data/kang/rds_just_counts/',folder[1],'.rds'),collapse='')
print(folder[1])
de10_1 <- readRDS(filee)
#WriteH5AD(de10_3, c('/home/erika/Documents/Projects/Muscat/muscat-comparison/data/sim_data/kang/anndata/',folder[1],'.h5ad'),collapse='')
filee <- paste(c('/home/erika/Documents/Projects/Muscat/muscat-comparison/data/sim_data/kang/rds_just_counts/',folder[2],'.rds'),collapse='')
de10_2 <- readRDS(filee)
filee <- paste(c('/home/erika/Documents/Projects/Muscat/muscat-comparison/data/sim_data/kang/rds_just_counts/',folder[3],'.rds'),collapse='')
de10_3 <- readRDS(filee)
filee <- paste(c('/home/erika/Documents/Projects/Muscat/muscat-comparison/data/sim_data/kang/rds_just_counts/',folder[4],'.rds'),collapse='')
de10_4 <- readRDS(filee)
filee <- paste(c('/home/erika/Documents/Projects/Muscat/muscat-comparison/data/sim_data/kang/rds_just_counts/',folder[5],'.rds'),collapse='')
de10_5 <- readRDS(filee)



## multiple files
folder = ['de10;1','de10;2','de10;3', 'de10;4', 'de10;5',
          'db10_1', 'db10;2','db10;3', 'db10;4', 'db10;5',
          'dm10;1', 'dm10;2', 'dm10;3','dm10;4', 'dm10;5',
          'dp10;1', 'dp10;2', 'dp10;3', 'dp10;4','dp10;5',
          'db10_1','db10;2', 'db10;3', 'db10;4','db10;5']
%%R -o de10_1,de10_2,de10_3,de10_4,de10_5,db10_1,db10_2,db10_3,db10_4,db10_5,dm10_1,dm10_2,dm10_3,dm10_4,dm10_5,dp10_1,dp10_2,dp10_3,dp10_4,dp10_5,db10_1,db10_2,db10_3,db10_4,db10_5 -i folder
#suppressPackageStartupMessages(library(Seurat))
library(SingleCellExperiment, lib.loc="/home/erika/soft/R-3.6.0/library")
library(scater, lib.loc="/home/erika/soft/R-3.6.0/library")
#print(sessionInfo())
#print(installed.packages()[, c("Package", "LibPath")])


filee <- paste(c('/home/erika/Documents/Projects/Muscat/muscat-comparison/data/sim_data/kang/rds_just_logcounts/',folder[1],'.rds'),collapse='')
print(folder[1])
de10_3 <- readRDS(filee)
#WriteH5AD(de10_3, c('/home/erika/Documents/Projects/Muscat/muscat-comparison/data/sim_data/kang/anndata/',folder[1],'.h5ad'),collapse='')
filee <- paste(c('/home/erika/Documents/Projects/Muscat/muscat-comparison/data/sim_data/kang/rds_just_logcounts/',folder[2],'.rds'),collapse='')
de10_4 <- readRDS(filee)
filee <- paste(c('/home/erika/Documents/Projects/Muscat/muscat-comparison/data/sim_data/kang/rds_just_logcounts/',folder[3],'.rds'),collapse='')
de10_5 <- readRDS(filee)

filee <- paste(c('/home/erika/Documents/Projects/Muscat/muscat-comparison/data/sim_data/kang/rds_just_logcounts/',folder[4],'.rds'),collapse='')
db10_1 <- readRDS(filee)
filee <- paste(c('/home/erika/Documents/Projects/Muscat/muscat-comparison/data/sim_data/kang/rds_just_logcounts/',folder[5],'.rds'),collapse='')
db10_2 <- readRDS(filee)
filee <- paste(c('/home/erika/Documents/Projects/Muscat/muscat-comparison/data/sim_data/kang/rds_just_logcounts/',folder[6],'.rds'),collapse='')
db10_3 <- readRDS(filee)
filee <- paste(c('/home/erika/Documents/Projects/Muscat/muscat-comparison/data/sim_data/kang/rds_just_logcounts/',folder[7],'.rds'),collapse='')
db10_4 <- readRDS(filee)
filee <- paste(c('/home/erika/Documents/Projects/Muscat/muscat-comparison/data/sim_data/kang/rds_just_logcounts/',folder[8],'.rds'),collapse='')
db10_5 <- readRDS(filee)

filee <- paste(c('/home/erika/Documents/Projects/Muscat/muscat-comparison/data/sim_data/kang/rds_just_logcounts/',folder[9],'.rds'),collapse='')
dm10_1 <- readRDS(filee)
filee <- paste(c('/home/erika/Documents/Projects/Muscat/muscat-comparison/data/sim_data/kang/rds_just_logcounts/',folder[10],'.rds'),collapse='')
dm10_2 <- readRDS(filee)
filee <- paste(c('/home/erika/Documents/Projects/Muscat/muscat-comparison/data/sim_data/kang/rds_just_logcounts/',folder[11],'.rds'),collapse='')
dm10_3 <- readRDS(filee)
filee <- paste(c('/home/erika/Documents/Projects/Muscat/muscat-comparison/data/sim_data/kang/rds_just_logcounts/',folder[12],'.rds'),collapse='')
dm10_4 <- readRDS(filee)
filee <- paste(c('/home/erika/Documents/Projects/Muscat/muscat-comparison/data/sim_data/kang/rds_just_logcounts/',folder[13],'.rds'),collapse='')
dm10_5 <- readRDS(filee)

filee <- paste(c('/home/erika/Documents/Projects/Muscat/muscat-comparison/data/sim_data/kang/rds_just_logcounts/',folder[14],'.rds'),collapse='')
dp10_1 <- readRDS(filee)
filee <- paste(c('/home/erika/Documents/Projects/Muscat/muscat-comparison/data/sim_data/kang/rds_just_logcounts/',folder[15],'.rds'),collapse='')
dp10_2 <- readRDS(filee)
filee <- paste(c('/home/erika/Documents/Projects/Muscat/muscat-comparison/data/sim_data/kang/rds_just_logcounts/',folder[16],'.rds'),collapse='')
dp10_3 <- readRDS(filee)
filee <- paste(c('/home/erika/Documents/Projects/Muscat/muscat-comparison/data/sim_data/kang/rds_just_logcounts/',folder[17],'.rds'),collapse='')
dp10_4 <- readRDS(filee)
filee <- paste(c('/home/erika/Documents/Projects/Muscat/muscat-comparison/data/sim_data/kang/rds_just_logcounts/',folder[18],'.rds'),collapse='')
dp10_5 <- readRDS(filee)

filee <- paste(c('/home/erika/Documents/Projects/Muscat/muscat-comparison/data/sim_data/kang/rds_just_logcounts/',folder[19],'.rds'),collapse='')
db10_2 <- readRDS(filee)
filee <- paste(c('/home/erika/Documents/Projects/Muscat/muscat-comparison/data/sim_data/kang/rds_just_logcounts/',folder[20],'.rds'),collapse='')
db10_3 <- readRDS(filee)
filee <- paste(c('/home/erika/Documents/Projects/Muscat/muscat-comparison/data/sim_data/kang/rds_just_logcounts/',folder[21],'.rds'),collapse='')
db10_4 <- readRDS(filee)
filee <- paste(c('/home/erika/Documents/Projects/Muscat/muscat-comparison/data/sim_data/kang/rds_just_logcounts/',folder[22],'.rds'),collapse='')
db10_5 <- readRDS(filee)

##scanpy normalization

sc.pp.normalize_total(de10_1)

## save anndata

#db10_1.write('/home/erika/Documents/Projects/Muscat/muscat-comparison/data/sim_data/kang/anndata/'+folder+'.h5ad')

#db10_1.add['gene_info'].to_csv('/home/erika/Documents/Projects/Muscat/muscat-comparison/data/sim_data/kang/anndata/'+folder+'_gene_info.csv')
#db10_1.add.to_csv('/home/erika/Documents/Projects/Muscat/muscat-comparison/data/sim_data/kang/anndata/'+folder+'_add.csv')
#with open('/home/erika/Documents/Projects/Muscat/muscat-comparison/data/sim_data/kang/anndata/'+folder+'_dict.p', 'wb') as fp:
#    pickle.dump(db10_1.add, fp, protocol=pickle.HIGHEST_PROTOCOL)
folder = ['de10;1','de10;2','de10;3', 'de10;4', 'de10;5']

de10_1.write('/home/erika/Documents/Projects/Muscat/muscat-comparison/data/sim_data/kang/anndata/counts/'+folder[0]+'.h5ad')
de10_1.add['gene_info'].to_csv('/home/erika/Documents/Projects/Muscat/muscat-comparison/data/sim_data/kang/anndata/counts/'+folder[0]+'_gene_info.csv')
de10_2.write('/home/erika/Documents/Projects/Muscat/muscat-comparison/data/sim_data/kang/anndata/counts/'+folder[1]+'.h5ad')
de10_2.add['gene_info'].to_csv('/home/erika/Documents/Projects/Muscat/muscat-comparison/data/sim_data/kang/anndata/counts/'+folder[1]+'_gene_info.csv')
de10_3.write('/home/erika/Documents/Projects/Muscat/muscat-comparison/data/sim_data/kang/anndata/counts/'+folder[2]+'.h5ad')
de10_3.add['gene_info'].to_csv('/home/erika/Documents/Projects/Muscat/muscat-comparison/data/sim_data/kang/anndata/counts/'+folder[2]+'_gene_info.csv')
de10_4.write('/home/erika/Documents/Projects/Muscat/muscat-comparison/data/sim_data/kang/anndata/counts/'+folder[3]+'.h5ad')
de10_4.add['gene_info'].to_csv('/home/erika/Documents/Projects/Muscat/muscat-comparison/data/sim_data/kang/anndata/counts/'+folder[3]+'_gene_info.csv')
de10_5.write('/home/erika/Documents/Projects/Muscat/muscat-comparison/data/sim_data/kang/anndata/counts/'+folder[4]+'.h5ad')
de10_5.add['gene_info'].to_csv('/home/erika/Documents/Projects/Muscat/muscat-comparison/data/sim_data/kang/anndata/counts/'+folder[4]+'_gene_info.csv')


#folder = 'de10;3'
#r_to_anndata(folder)






## save anndata

#db10_1.write('/home/erika/Documents/Projects/Muscat/muscat-comparison/data/sim_data/kang/anndata/'+folder+'.h5ad')

#db10_1.add['gene_info'].to_csv('/home/erika/Documents/Projects/Muscat/muscat-comparison/data/sim_data/kang/anndata/'+folder+'_gene_info.csv')
#db10_1.add.to_csv('/home/erika/Documents/Projects/Muscat/muscat-comparison/data/sim_data/kang/anndata/'+folder+'_add.csv')
#with open('/home/erika/Documents/Projects/Muscat/muscat-comparison/data/sim_data/kang/anndata/'+folder+'_dict.p', 'wb') as fp:
#    pickle.dump(db10_1.add, fp, protocol=pickle.HIGHEST_PROTOCOL)
folder = ['de10;3', 'de10;4', 'de10;5',
          'db10_1', 'db10;2','db10;3', 'db10;4', 'db10;5',
          'dm10;1', 'dm10;2', 'dm10;3','dm10;4', 'dm10;5',
          'dp10;1', 'dp10;2', 'dp10;3', 'dp10;4','dp10;5',
          'db10;2', 'db10;3', 'db10;4','db10;5']

de10_3.write('/home/erika/Documents/Projects/Muscat/muscat-comparison/data/sim_data/kang/anndata/'+folder[0]+'.h5ad')
de10_3.add['gene_info'].to_csv('/home/erika/Documents/Projects/Muscat/muscat-comparison/data/sim_data/kang/anndata/'+folder[0]+'_gene_info.csv')
de10_4.write('/home/erika/Documents/Projects/Muscat/muscat-comparison/data/sim_data/kang/anndata/'+folder[1]+'.h5ad')
de10_4.add['gene_info'].to_csv('/home/erika/Documents/Projects/Muscat/muscat-comparison/data/sim_data/kang/anndata/'+folder[1]+'_gene_info.csv')
de10_5.write('/home/erika/Documents/Projects/Muscat/muscat-comparison/data/sim_data/kang/anndata/'+folder[2]+'.h5ad')
de10_5.add['gene_info'].to_csv('/home/erika/Documents/Projects/Muscat/muscat-comparison/data/sim_data/kang/anndata/'+folder[2]+'_gene_info.csv')

db10_1.write('/home/erika/Documents/Projects/Muscat/muscat-comparison/data/sim_data/kang/anndata/'+folder[3]+'.h5ad')
db10_1.add['gene_info'].to_csv('/home/erika/Documents/Projects/Muscat/muscat-comparison/data/sim_data/kang/anndata/'+folder[3]+'_gene_info.csv')
db10_2.write('/home/erika/Documents/Projects/Muscat/muscat-comparison/data/sim_data/kang/anndata/'+folder[4]+'.h5ad')
db10_2.add['gene_info'].to_csv('/home/erika/Documents/Projects/Muscat/muscat-comparison/data/sim_data/kang/anndata/'+folder[4]+'_gene_info.csv')
db10_3.write('/home/erika/Documents/Projects/Muscat/muscat-comparison/data/sim_data/kang/anndata/'+folder[5]+'.h5ad')
db10_3.add['gene_info'].to_csv('/home/erika/Documents/Projects/Muscat/muscat-comparison/data/sim_data/kang/anndata/'+folder[5]+'_gene_info.csv')
db10_4.write('/home/erika/Documents/Projects/Muscat/muscat-comparison/data/sim_data/kang/anndata/'+folder[6]+'.h5ad')
db10_4.add['gene_info'].to_csv('/home/erika/Documents/Projects/Muscat/muscat-comparison/data/sim_data/kang/anndata/'+folder[6]+'_gene_info.csv')
db10_5.write('/home/erika/Documents/Projects/Muscat/muscat-comparison/data/sim_data/kang/anndata/'+folder[7]+'.h5ad')
db10_5.add['gene_info'].to_csv('/home/erika/Documents/Projects/Muscat/muscat-comparison/data/sim_data/kang/anndata/'+folder[7]+'_gene_info.csv')

dm10_1.write('/home/erika/Documents/Projects/Muscat/muscat-comparison/data/sim_data/kang/anndata/'+folder[8]+'.h5ad')
dm10_1.add['gene_info'].to_csv('/home/erika/Documents/Projects/Muscat/muscat-comparison/data/sim_data/kang/anndata/'+folder[8]+'_gene_info.csv')
dm10_2.write('/home/erika/Documents/Projects/Muscat/muscat-comparison/data/sim_data/kang/anndata/'+folder[9]+'.h5ad')
dm10_2.add['gene_info'].to_csv('/home/erika/Documents/Projects/Muscat/muscat-comparison/data/sim_data/kang/anndata/'+folder[9]+'_gene_info.csv')
dm10_3.write('/home/erika/Documents/Projects/Muscat/muscat-comparison/data/sim_data/kang/anndata/'+folder[10]+'.h5ad')
dm10_3.add['gene_info'].to_csv('/home/erika/Documents/Projects/Muscat/muscat-comparison/data/sim_data/kang/anndata/'+folder[10]+'_gene_info.csv')
dm10_4.write('/home/erika/Documents/Projects/Muscat/muscat-comparison/data/sim_data/kang/anndata/'+folder[11]+'.h5ad')
dm10_4.add['gene_info'].to_csv('/home/erika/Documents/Projects/Muscat/muscat-comparison/data/sim_data/kang/anndata/'+folder[11]+'_gene_info.csv')
dm10_5.write('/home/erika/Documents/Projects/Muscat/muscat-comparison/data/sim_data/kang/anndata/'+folder[12]+'.h5ad')
dm10_5.add['gene_info'].to_csv('/home/erika/Documents/Projects/Muscat/muscat-comparison/data/sim_data/kang/anndata/'+folder[12]+'_gene_info.csv')

dp10_1.write('/home/erika/Documents/Projects/Muscat/muscat-comparison/data/sim_data/kang/anndata/'+folder[13]+'.h5ad')
dp10_1.add['gene_info'].to_csv('/home/erika/Documents/Projects/Muscat/muscat-comparison/data/sim_data/kang/anndata/'+folder[13]+'_gene_info.csv')
dp10_2.write('/home/erika/Documents/Projects/Muscat/muscat-comparison/data/sim_data/kang/anndata/'+folder[14]+'.h5ad')
dp10_2.add['gene_info'].to_csv('/home/erika/Documents/Projects/Muscat/muscat-comparison/data/sim_data/kang/anndata/'+folder[14]+'_gene_info.csv')
dp10_3.write('/home/erika/Documents/Projects/Muscat/muscat-comparison/data/sim_data/kang/anndata/'+folder[15]+'.h5ad')
dp10_3.add['gene_info'].to_csv('/home/erika/Documents/Projects/Muscat/muscat-comparison/data/sim_data/kang/anndata/'+folder[15]+'_gene_info.csv')
dp10_4.write('/home/erika/Documents/Projects/Muscat/muscat-comparison/data/sim_data/kang/anndata/'+folder[16]+'.h5ad')
dp10_4.add['gene_info'].to_csv('/home/erika/Documents/Projects/Muscat/muscat-comparison/data/sim_data/kang/anndata/'+folder[16]+'_gene_info.csv')
dp10_5.write('/home/erika/Documents/Projects/Muscat/muscat-comparison/data/sim_data/kang/anndata/'+folder[17]+'.h5ad')
dp10_5.add['gene_info'].to_csv('/home/erika/Documents/Projects/Muscat/muscat-comparison/data/sim_data/kang/anndata/'+folder[17]+'_gene_info.csv')

db10_2.write('/home/erika/Documents/Projects/Muscat/muscat-comparison/data/sim_data/kang/anndata/'+folder[18]+'.h5ad')
db10_2.add['gene_info'].to_csv('/home/erika/Documents/Projects/Muscat/muscat-comparison/data/sim_data/kang/anndata/'+folder[18]+'_gene_info.csv')
db10_3.write('/home/erika/Documents/Projects/Muscat/muscat-comparison/data/sim_data/kang/anndata/'+folder[19]+'.h5ad')
db10_3.add['gene_info'].to_csv('/home/erika/Documents/Projects/Muscat/muscat-comparison/data/sim_data/kang/anndata/'+folder[19]+'_gene_info.csv')
db10_4.write('/home/erika/Documents/Projects/Muscat/muscat-comparison/data/sim_data/kang/anndata/'+folder[20]+'.h5ad')
db10_4.add['gene_info'].to_csv('/home/erika/Documents/Projects/Muscat/muscat-comparison/data/sim_data/kang/anndata/'+folder[20]+'_gene_info.csv')
db10_5.write('/home/erika/Documents/Projects/Muscat/muscat-comparison/data/sim_data/kang/anndata/'+folder[21]+'.h5ad')
db10_5.add['gene_info'].to_csv('/home/erika/Documents/Projects/Muscat/muscat-comparison/data/sim_data/kang/anndata/'+folder[21]+'_gene_info.csv')

#folder = 'de10;3'
#r_to_anndata(folder)






## save patient matrices for different celltypes


main_dir = '/home/erika/Documents/Projects/Muscat/muscat-comparison/data/sim_data/kang/pandasDF/counts/'
sparsematrix = False
for fii in range(len(folder)):
    if not os.path.exists(main_dir + folder[fii]):
        os.mkdir(main_dir + folder[fii])
##
# transform anndata to format for DE-method and save pandas files
#file = 'db10_1'
#
de_all = [de10_1,de10_2,de10_3,de10_4,de10_5]

dii = 1
for dii in range(len(de_all)):

    # how many clusters/samples exist:
    nr_cl = np.shape(de_all[dii].add['gene_info'].loc[:,'cluster_id'].unique())[0]
    nr_s = np.shape(de_all[dii].add['experiment_info'].loc[:,'sample_id'].unique())[0]


    for cl in range(1,nr_cl+1):
        # subsets of celltypes
        #cl = 1
        db10_cl1 = de_all[dii][de_all[dii].obs['cluster_id']==cl]

        # subset of samples/patients
        for s in range(1,nr_s+1):
            #s = 1
            db10_cl1_s1 = db10_cl1[db10_cl1.obs['sample_id']==s]

            # transform to pandas-dataframe
            if sparsematrix:
                pd_db10_cl1_s1 = pd.DataFrame(data=np.transpose(db10_cl1_s1.X.toarray()),
                                              index=db10_cl1_s1.var_names,
                                              columns=db10_cl1_s1.obs_names)
                pd_db10_cl1_s1.to_csv(
                    main_dir + folder[dii] + '/cl' + str(cl) + 'Pt' + str(s))
            else:
                pd_db10_cl1_s1 = pd.DataFrame(data=np.transpose(db10_cl1_s1.X),
                              index= db10_cl1_s1.var_names,
                              columns= db10_cl1_s1.obs_names )
                pd_db10_cl1_s1.to_csv(main_dir + folder[dii] + '/cl' + str(cl) + 'Pt' + str(s) )

# # new name for the "next loop"
# #folder = 'nill;3'
# ##
##
#sc.pp.recipe_seurat(de10_1, log=True, plot=False, copy=False)
##
sc.pp.normalize_total(de10_1, target_sum=1)