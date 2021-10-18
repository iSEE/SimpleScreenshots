library(iSEE)
library(TENxPBMCData)
library(scater)
library(scuttle)
library(scran)
library(igraph)

pbmc3k <- TENxPBMCData("pbmc3k")

pbmc3k <- logNormCounts(pbmc3k)
pbmc3k <- runPCA(pbmc3k)
pbmc3k <- runTSNE(pbmc3k)

g <- buildSNNGraph(pbmc3k, k=10, use.dimred = 'PCA')
pbmc3k$Cluster <- as.factor(igraph::cluster_walktrap(g)$membership)

rownames(pbmc3k) <- rowData(pbmc3k)[["Symbol_TENx"]]
colnames(pbmc3k) <- paste0("Cell", seq_len(ncol(pbmc3k)))

pbmc3k <- addPerCellQCMetrics(pbmc3k)

# initial ----

initial <- list()

################################################################################
# Settings for Reduced dimension plot 1
################################################################################

initial[["ReducedDimensionPlot1"]] <- new("ReducedDimensionPlot", Type = "TSNE", XAxis = 1L, YAxis = 2L,
    FacetRowByColData = "Sample", FacetColumnByColData = "Sample",
    ColorByColumnData = "Cluster", ColorByFeatureNameAssay = "logcounts",
    ColorBySampleNameColor = "#FF0000", ShapeByColumnData = "Sample",
    SizeByColumnData = "Library", FacetRowBy = "None", FacetColumnBy = "None",
    ColorBy = "Column data", ColorByDefaultColor = "#000000",
    ColorByFeatureName = "MIR1302-10", ColorByFeatureSource = "---",
    ColorByFeatureDynamicSource = FALSE, ColorBySampleName = "Cell1",
    ColorBySampleSource = "---", ColorBySampleDynamicSource = FALSE,
    ShapeBy = "None", SizeBy = "None", SelectionAlpha = 0.1,
    ZoomData = numeric(0), BrushData = list(lasso = NULL, closed = FALSE,
        panelvar1 = NULL, panelvar2 = NULL, mapping = list(x = "X",
            y = "Y", colour = "ColorBy"), coord = structure(c(-12.8470863776797,
        3.40480680108348), .Dim = 1:2)), VisualBoxOpen = TRUE,
    VisualChoices = "Color", ContourAdd = FALSE, ContourColor = "#0000FF",
    PointSize = 1, PointAlpha = 1, Downsample = FALSE, DownsampleResolution = 200,
    CustomLabels = FALSE, CustomLabelsText = "Cell1", FontSize = 1,
    LegendPointSize = 1, LegendPosition = "Bottom", HoverInfo = TRUE,
    LabelCenters = FALSE, LabelCentersBy = "Sample", LabelCentersColor = "#000000",
    VersionInfo = list(iSEE = structure(list(c(2L, 4L, 0L)), class = c("package_version",
    "numeric_version"))), PanelId = c(ReducedDimensionPlot = 1L),
    PanelHeight = 500L, PanelWidth = 6L, SelectionBoxOpen = TRUE,
    RowSelectionSource = "---", ColumnSelectionSource = "---",
    DataBoxOpen = TRUE, RowSelectionDynamicSource = FALSE, ColumnSelectionDynamicSource = FALSE,
    RowSelectionRestrict = FALSE, ColumnSelectionRestrict = FALSE,
    SelectionHistory = list())

################################################################################
# Settings for Feature assay plot 1
################################################################################

initial[["FeatureAssayPlot1"]] <- new("FeatureAssayPlot", Assay = "logcounts", XAxis = "Column data",
    XAxisColumnData = "Cluster", XAxisFeatureName = "MIR1302-10",
    XAxisFeatureSource = "---", XAxisFeatureDynamicSource = FALSE,
    YAxisFeatureName = "CD3D", YAxisFeatureSource = "---", YAxisFeatureDynamicSource = FALSE,
    FacetRowByColData = "Sample", FacetColumnByColData = "Sample",
    ColorByColumnData = "Cluster", ColorByFeatureNameAssay = "logcounts",
    ColorBySampleNameColor = "#FF0000", ShapeByColumnData = "Sample",
    SizeByColumnData = "Library", FacetRowBy = "None", FacetColumnBy = "None",
    ColorBy = "Column data", ColorByDefaultColor = "#000000",
    ColorByFeatureName = "MIR1302-10", ColorByFeatureSource = "---",
    ColorByFeatureDynamicSource = FALSE, ColorBySampleName = "Cell1",
    ColorBySampleSource = "---", ColorBySampleDynamicSource = FALSE,
    ShapeBy = "None", SizeBy = "None", SelectionAlpha = 0.1,
    ZoomData = numeric(0), BrushData = list(), VisualBoxOpen = TRUE,
    VisualChoices = "Color", ContourAdd = FALSE, ContourColor = "#0000FF",
    PointSize = 1, PointAlpha = 1, Downsample = FALSE, DownsampleResolution = 200,
    CustomLabels = FALSE, CustomLabelsText = "Cell1", FontSize = 1,
    LegendPointSize = 1, LegendPosition = "Bottom", HoverInfo = TRUE,
    LabelCenters = FALSE, LabelCentersBy = "Sample", LabelCentersColor = "#000000",
    VersionInfo = list(iSEE = structure(list(c(2L, 4L, 0L)), class = c("package_version",
    "numeric_version"))), PanelId = c(FeatureAssayPlot = 1L),
    PanelHeight = 500L, PanelWidth = 6L, SelectionBoxOpen = TRUE,
    RowSelectionSource = "---", ColumnSelectionSource = "---",
    DataBoxOpen = TRUE, RowSelectionDynamicSource = FALSE, ColumnSelectionDynamicSource = FALSE,
    RowSelectionRestrict = FALSE, ColumnSelectionRestrict = FALSE,
    SelectionHistory = list())

# iSEE ----

app <- iSEE(pbmc3k, initial = initial)

shiny::runApp(app, launch.browser = TRUE)
