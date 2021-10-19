library(iSEE)
library(TENxPBMCData)
library(scater)
library(scuttle)
library(scran)
library(igraph)
library(scRNAseq)

# ?createLandingPage

# Alternative approach, to create a landing page
# that opens one of the datasets from the scRNAseq package.
library(scRNAseq)
all.data <- ls("package:scRNAseq")
all.data <- all.data[grep("Data$", all.data)]

lpfun <- createLandingPage(
    seUI=function(id) selectInput(id, "Dataset:", choices=all.data),
    seLoad=function(x) get(x, as.environment("package:scRNAseq"))()
)

# iSEE ----

app <- iSEE(landingPage=lpfun)

shiny::runApp(app, launch.browser = TRUE)
