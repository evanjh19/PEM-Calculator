if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

packs <- c('data.table',
           'DT',
           'abind',
           'scater',
           'SummarizedExperiment',
           'RColorBrewer',
           'shiny',
           'pheatmap',
           'reader')

for (i in packs) {
    
    if (!require(i, quietly = TRUE))
        BiocManager::install(i)
    
}
