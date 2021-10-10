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
    ColorByFeatureDynamicSource = FALSE, ColorBySampleName = "1",
    ColorBySampleSource = "---", ColorBySampleDynamicSource = FALSE,
    ShapeBy = "None", SizeBy = "None", SelectionAlpha = 0.1,
    ZoomData = numeric(0), BrushData = list(lasso = NULL, closed = FALSE,
        panelvar1 = NULL, panelvar2 = NULL, mapping = list(x = "X",
            y = "Y", colour = "ColorBy"), coord = structure(c(-12.8470863776797,
        3.40480680108348), .Dim = 1:2)), VisualBoxOpen = FALSE,
    VisualChoices = "Color", ContourAdd = FALSE, ContourColor = "#0000FF",
    PointSize = 1, PointAlpha = 1, Downsample = FALSE, DownsampleResolution = 200,
    CustomLabels = FALSE, CustomLabelsText = "1", FontSize = 1,
    LegendPointSize = 1, LegendPosition = "Bottom", HoverInfo = TRUE,
    LabelCenters = FALSE, LabelCentersBy = "Sample", LabelCentersColor = "#000000",
    VersionInfo = list(iSEE = structure(list(c(2L, 4L, 0L)), class = c("package_version",
    "numeric_version"))), PanelId = c(ReducedDimensionPlot = 1L),
    PanelHeight = 500L, PanelWidth = 3L, SelectionBoxOpen = FALSE,
    RowSelectionSource = "---", ColumnSelectionSource = "---",
    DataBoxOpen = FALSE, RowSelectionDynamicSource = FALSE, ColumnSelectionDynamicSource = FALSE,
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
    ColorByFeatureDynamicSource = FALSE, ColorBySampleName = "1",
    ColorBySampleSource = "---", ColorBySampleDynamicSource = FALSE,
    ShapeBy = "None", SizeBy = "None", SelectionAlpha = 0.1,
    ZoomData = numeric(0), BrushData = list(), VisualBoxOpen = FALSE,
    VisualChoices = "Color", ContourAdd = FALSE, ContourColor = "#0000FF",
    PointSize = 1, PointAlpha = 1, Downsample = FALSE, DownsampleResolution = 200,
    CustomLabels = FALSE, CustomLabelsText = "1", FontSize = 1,
    LegendPointSize = 1, LegendPosition = "Bottom", HoverInfo = TRUE,
    LabelCenters = FALSE, LabelCentersBy = "Sample", LabelCentersColor = "#000000",
    VersionInfo = list(iSEE = structure(list(c(2L, 4L, 0L)), class = c("package_version",
    "numeric_version"))), PanelId = c(FeatureAssayPlot = 1L),
    PanelHeight = 500L, PanelWidth = 3L, SelectionBoxOpen = FALSE,
    RowSelectionSource = "---", ColumnSelectionSource = "---",
    DataBoxOpen = FALSE, RowSelectionDynamicSource = FALSE, ColumnSelectionDynamicSource = FALSE,
    RowSelectionRestrict = FALSE, ColumnSelectionRestrict = FALSE,
    SelectionHistory = list())

################################################################################
# Settings for Column data plot 1
################################################################################

initial[["ColumnDataPlot1"]] <- new("ColumnDataPlot", XAxis = "None", YAxis = "Cluster", XAxisColumnData = "Sample",
    FacetRowByColData = "Sample", FacetColumnByColData = "Sample",
    ColorByColumnData = "Sample", ColorByFeatureNameAssay = "logcounts",
    ColorBySampleNameColor = "#FF0000", ShapeByColumnData = "Sample",
    SizeByColumnData = "Library", FacetRowBy = "None", FacetColumnBy = "None",
    ColorBy = "None", ColorByDefaultColor = "#000000", ColorByFeatureName = "MIR1302-10",
    ColorByFeatureSource = "---", ColorByFeatureDynamicSource = FALSE,
    ColorBySampleName = "1", ColorBySampleSource = "---", ColorBySampleDynamicSource = FALSE,
    ShapeBy = "None", SizeBy = "None", SelectionAlpha = 0.1,
    ZoomData = numeric(0), BrushData = list(), VisualBoxOpen = FALSE,
    VisualChoices = "Color", ContourAdd = FALSE, ContourColor = "#0000FF",
    PointSize = 1, PointAlpha = 1, Downsample = FALSE, DownsampleResolution = 200,
    CustomLabels = FALSE, CustomLabelsText = "1", FontSize = 1,
    LegendPointSize = 1, LegendPosition = "Bottom", HoverInfo = TRUE,
    LabelCenters = FALSE, LabelCentersBy = "Sample", LabelCentersColor = "#000000",
    VersionInfo = list(iSEE = structure(list(c(2L, 4L, 0L)), class = c("package_version",
    "numeric_version"))), PanelId = c(ColumnDataPlot = 1L), PanelHeight = 500L,
    PanelWidth = 3L, SelectionBoxOpen = FALSE, RowSelectionSource = "---",
    ColumnSelectionSource = "---", DataBoxOpen = FALSE, RowSelectionDynamicSource = FALSE,
    ColumnSelectionDynamicSource = FALSE, RowSelectionRestrict = FALSE,
    ColumnSelectionRestrict = FALSE, SelectionHistory = list())

################################################################################
# Settings for Row data plot 1
################################################################################

initial[["RowDataPlot1"]] <- new("RowDataPlot", XAxis = "Row data", YAxis = "n_cells", XAxisRowData = "log10_total",
    FacetRowByRowData = "ENSEMBL_ID", FacetColumnByRowData = "ENSEMBL_ID",
    ColorByRowData = "ENSEMBL_ID", ColorBySampleNameAssay = "logcounts",
    ColorByFeatureNameColor = "#FF0000", ShapeByRowData = "ENSEMBL_ID",
    SizeByRowData = "n_cells", FacetRowBy = "None", FacetColumnBy = "None",
    ColorBy = "None", ColorByDefaultColor = "#000000", ColorByFeatureName = "MIR1302-10",
    ColorByFeatureSource = "---", ColorByFeatureDynamicSource = FALSE,
    ColorBySampleName = "1", ColorBySampleSource = "---", ColorBySampleDynamicSource = FALSE,
    ShapeBy = "None", SizeBy = "None", SelectionAlpha = 0.1,
    ZoomData = numeric(0), BrushData = list(), VisualBoxOpen = FALSE,
    VisualChoices = "Color", ContourAdd = FALSE, ContourColor = "#0000FF",
    PointSize = 1, PointAlpha = 1, Downsample = FALSE, DownsampleResolution = 200,
    CustomLabels = FALSE, CustomLabelsText = "MIR1302-10", FontSize = 1,
    LegendPointSize = 1, LegendPosition = "Bottom", HoverInfo = TRUE,
    LabelCenters = FALSE, LabelCentersBy = "ENSEMBL_ID", LabelCentersColor = "#000000",
    VersionInfo = list(iSEE = structure(list(c(2L, 4L, 0L)), class = c("package_version",
    "numeric_version"))), PanelId = c(RowDataPlot = 1L), PanelHeight = 500L,
    PanelWidth = 3L, SelectionBoxOpen = FALSE, RowSelectionSource = "---",
    ColumnSelectionSource = "---", DataBoxOpen = FALSE, RowSelectionDynamicSource = FALSE,
    ColumnSelectionDynamicSource = FALSE, RowSelectionRestrict = FALSE,
    ColumnSelectionRestrict = FALSE, SelectionHistory = list())

################################################################################
# Settings for Row data table 1
################################################################################

initial[["RowDataTable1"]] <- new("RowDataTable", Selected = "MIR1302-10", Search = "", SearchColumns = c("",
"", "", "", ""), HiddenColumns = character(0), VersionInfo = list(
    iSEE = structure(list(c(2L, 4L, 0L)), class = c("package_version",
    "numeric_version"))), PanelId = c(RowDataTable = 1L), PanelHeight = 500L,
    PanelWidth = 3L, SelectionBoxOpen = FALSE, RowSelectionSource = "---",
    ColumnSelectionSource = "---", DataBoxOpen = FALSE, RowSelectionDynamicSource = FALSE,
    ColumnSelectionDynamicSource = FALSE, RowSelectionRestrict = FALSE,
    ColumnSelectionRestrict = FALSE, SelectionHistory = list())

################################################################################
# Settings for Column data table 1
################################################################################

initial[["ColumnDataTable1"]] <- new("ColumnDataTable", Selected = "1", Search = "", SearchColumns = c("",
"", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "",
"", ""), HiddenColumns = character(0), VersionInfo = list(iSEE = structure(list(
    c(2L, 4L, 0L)), class = c("package_version", "numeric_version"
))), PanelId = c(ColumnDataTable = 1L), PanelHeight = 500L, PanelWidth = 3L,
    SelectionBoxOpen = FALSE, RowSelectionSource = "---", ColumnSelectionSource = "---",
    DataBoxOpen = FALSE, RowSelectionDynamicSource = FALSE, ColumnSelectionDynamicSource = FALSE,
    RowSelectionRestrict = FALSE, ColumnSelectionRestrict = FALSE,
    SelectionHistory = list())

################################################################################
# Settings for Complex heatmap 1
################################################################################

initial[["ComplexHeatmapPlot1"]] <- new("ComplexHeatmapPlot", Assay = "logcounts", CustomRows = TRUE,
    CustomRowsText = "CD79A\nCD3D\nFOLR3\nS100A8", ClusterRows = FALSE,
    ClusterRowsDistance = "spearman", ClusterRowsMethod = "ward.D2",
    DataBoxOpen = FALSE, VisualChoices = "Annotations", ColumnData = "Cluster",
    RowData = character(0), CustomBounds = FALSE, LowerBound = NA_real_,
    UpperBound = NA_real_, AssayCenterRows = FALSE, AssayScaleRows = FALSE,
    DivergentColormap = "purple < black < yellow", ShowDimNames = "Rows",
    LegendPosition = "Bottom", LegendDirection = "Horizontal",
    VisualBoxOpen = FALSE, NamesRowFontSize = 10, NamesColumnFontSize = 10,
    ShowColumnSelection = TRUE, OrderColumnSelection = TRUE,
    VersionInfo = list(iSEE = structure(list(c(2L, 4L, 0L)), class = c("package_version",
    "numeric_version"))), PanelId = c(ComplexHeatmapPlot = 1L),
    PanelHeight = 500L, PanelWidth = 3L, SelectionBoxOpen = FALSE,
    RowSelectionSource = "---", ColumnSelectionSource = "---",
    RowSelectionDynamicSource = FALSE, ColumnSelectionDynamicSource = FALSE,
    RowSelectionRestrict = FALSE, ColumnSelectionRestrict = FALSE,
    SelectionHistory = list())

################################################################################
# Settings for Sample assay plot 1
################################################################################

initial[["SampleAssayPlot1"]] <- new("SampleAssayPlot", Assay = "logcounts", XAxis = "None", XAxisRowData = "ENSEMBL_ID",
    XAxisSampleName = "1", XAxisSampleSource = "---", XAxisSampleDynamicSource = FALSE,
    YAxisSampleName = "1", YAxisSampleSource = "---", YAxisSampleDynamicSource = FALSE,
    FacetRowByRowData = "ENSEMBL_ID", FacetColumnByRowData = "ENSEMBL_ID",
    ColorByRowData = "ENSEMBL_ID", ColorBySampleNameAssay = "logcounts",
    ColorByFeatureNameColor = "#FF0000", ShapeByRowData = "ENSEMBL_ID",
    SizeByRowData = "n_cells", FacetRowBy = "None", FacetColumnBy = "None",
    ColorBy = "None", ColorByDefaultColor = "#000000", ColorByFeatureName = "MIR1302-10",
    ColorByFeatureSource = "---", ColorByFeatureDynamicSource = FALSE,
    ColorBySampleName = "1", ColorBySampleSource = "---", ColorBySampleDynamicSource = FALSE,
    ShapeBy = "None", SizeBy = "None", SelectionAlpha = 0.1,
    ZoomData = numeric(0), BrushData = list(), VisualBoxOpen = FALSE,
    VisualChoices = "Color", ContourAdd = FALSE, ContourColor = "#0000FF",
    PointSize = 1, PointAlpha = 1, Downsample = FALSE, DownsampleResolution = 200,
    CustomLabels = FALSE, CustomLabelsText = "MIR1302-10", FontSize = 1,
    LegendPointSize = 1, LegendPosition = "Bottom", HoverInfo = TRUE,
    LabelCenters = FALSE, LabelCentersBy = "ENSEMBL_ID", LabelCentersColor = "#000000",
    VersionInfo = list(iSEE = structure(list(c(2L, 4L, 0L)), class = c("package_version",
    "numeric_version"))), PanelId = c(SampleAssayPlot = 1L),
    PanelHeight = 500L, PanelWidth = 3L, SelectionBoxOpen = FALSE,
    RowSelectionSource = "---", ColumnSelectionSource = "---",
    DataBoxOpen = FALSE, RowSelectionDynamicSource = FALSE, ColumnSelectionDynamicSource = FALSE,
    RowSelectionRestrict = FALSE, ColumnSelectionRestrict = FALSE,
    SelectionHistory = list())

# iSEE ----
iSEE(pbmc3k, initial = initial)
