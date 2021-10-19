library(iSEE)
library(TENxPBMCData)
library(scater)
library(scuttle)
library(scran)
library(igraph)
library(scRNAseq)

# Example data ----
sce <- ReprocessedAllenData(assays="tophat_counts")
class(sce)

library(scater)
sce <- logNormCounts(sce, exprs_values="tophat_counts")

sce <- runPCA(sce, ncomponents=4)
sce <- runTSNE(sce)
rowData(sce)$ave_count <- rowMeans(assay(sce, "tophat_counts"))
rowData(sce)$n_cells <- rowSums(assay(sce, "tophat_counts") > 0)

# iSEE ----

app <- iSEE(sce)

shiny::runApp(app, launch.browser = TRUE)
