###################################################################################
##############			load package 		###########################
###################################################################################

library("popsicle")

##################################################################################################
#################               Pre-processing              ######################################
##################################################################################################

setwd("/path/to/your/working/directory")

populations_markers= c("PTPRC","CD34","PROM1","CD3D","CD3E","TNFRSF4","CD4","IL7R")

sample_name <- PrePlots("sample_name", input_data= "/path/to/your/data/", genelist=populations_markers)

saveRDS(sample_name,file="sample_name_Pre_Filter.Rds")

##################################################################################################
#################             Filter data and plot          ######################################
##################################################################################################

sample_name <- FilterPlots(sample_name, G_RNA_low = 500, G_RNA_hi= Inf, U_RNA_low = -Inf, U_RNA_hi = 37000, percent_mt_hi = 10, percent_ribo_hi= 100, percent_disso_hi = 100)

saveRDS(sample_name,file="sample_name_After_Filters.Rds")

##################################################################################################
#################                Doublets                   ######################################
##################################################################################################

#### run this function at least 2 times, as explained on your terminal

#### first run:

sample_name <- CalculateDoublets(sample_name, dbs_thr='none', dbs_remove=FALSE)

#### second run: 

sample_name <- CalculateDoublets(sample_name, dbs_thr=0.22, dbs_remove=FALSE)

saveRDS(sample_name,file="sample_name_after_CalculateDoublets.Rds")

##################################################################################################
#################              Normalization                ######################################
##################################################################################################

sample_name <- Normalize(sample_name, variable_genes=2000)

saveRDS(sample_name ,file="sample_name_after_Normalization.Rds")

##################################################################################################
#################             CellCycleScore                ######################################
##################################################################################################

sample_name <- CCScore(sample_name)

saveRDS(sample_name ,file="sample_name_after_CCScore.Rds")

##################################################################################################
#################                Regression                 ######################################
##################################################################################################

sample_name <- ApplyRegression(UMI= sample_name, variables= c("nCount_RNA", "percent_mt", "S.Score", "G2M.Score"), explore_PC=TRUE)

saveRDS(sample_name, file="sample_name_after_regression.Rds")

##################################################################################################
#################                 clustering                ######################################
##################################################################################################

sample_name <- CalculateCluster(sample_name, 12, "human")

saveRDS(sample_name, file="sample_name_cluster.Rds")

##################################################################################################
#################                 Annotation                ######################################
##################################################################################################

sample_name <- MakeAnnotation(sample_name, organism="human")

saveRDS(sample_name, file="sample_name.Rds")


