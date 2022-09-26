#Install BiocManager if not already installed
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

#Install SummarizedExperiment package
BiocManager::install("SummarizedExperiment")

#Load SummarizedExperiment package
library(SummarizedExperiment)

#Load count data matrix
#Duplicate this command for as many matrices as you wish to include
#Matrices must be of the same dimensions
countmatrix <- read.csv('path_to_count_matrix.csv')

#Load metadata table
#This table consists of an initial column header row describing the metadata type, i.e. batch or sample_name
#Subsequent rows should be sample information such as the names for samples or the batch that a sample belongs to
metadata <- read.csv('path_to_metadata_table.csv')

#Create Summarized Experiment from count data and metadata
#If multiple count data matrices are to be incorporated, the list may be appended
sce <- SummarizedExperiment(list(counts = countmatrix),
                            colData = metadata)

#Save the SummarizedExperiment
saveRDS(sce, 'path_to_save_SummarizedExperiment.RDS')