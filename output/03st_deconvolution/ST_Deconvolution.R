
library(STdeconvolve)
library(BiocManager)
library(SpatialExperiment)
# install.packages("remotes")
require(remotes)
remotes::install_github('JEFworks-Lab/MERINGUE', build_vignettes = TRUE)
library(MERINGUE)


setwd("~/spatial_prions/output/03st_deconvolution")
folder_path <- "/Users/davidecaredio/spatial_prions/output/03st_deconvolution"  # Replace with the actual folder path
files <- list.files(path = folder_path, full.names = TRUE)
# List of RDS files
file_names <- c(
  "ME7_30_unfiltered_spe.rds", "mNS_30_unfiltered_spe.rds", "NBH_27_unfiltered_spe.rds",
  "NBH_30_1_unfiltered_spe.rds", "NBH_30_2_unfiltered_spe.rds", "NBH_TS_1_unfiltered_spe.rds",
  "NBH_TS_2_unfiltered_spe.rds", "RML6_27_unfiltered_spe.rds", "RML6_30_1_unfiltered_spe.rds",
  "RML6_30_2_unfiltered_spe.rds", "RML6_TS_1_unfiltered_spe.rds", "RML6_TS_2_unfiltered_spe.rds"
)

# Create an empty list to store the imported data
data_list <- list()

# Loop through the list of file names and read each RDS file
for (file_name in file_names) {
  file_path <- file.path(folder_path, file_name)  # Assuming folder_path is set correctly
  data <- readRDS(file_path)
  data_list[[file_name]] <- data
}


#________________________GENES TO REMOVE: Hemoglobin and Mitochondrial related genes___________________________________________
RML6_TS_2_unfiltered_spe <- readRDS("RML6_TS_2_unfiltered_spe.rds")

# Define the pattern to remove
#Obteined from here: https://doi.org/10.1038/s41598-020-62801-6
pattern_to_remove <- "^mt-|^Hba1-|^Hba2-|^Hbb-|^Hbbp1-|^Hbd-|^Hbe1-|^Hbg1-|^Hbg2-|^Hmm-|^Hbq1-|^Hbz-|^Hbzp1- |^Hba-"

# Find the indices of genes that match the pattern
genes_to_remove <- grep(pattern_to_remove, RML6_TS_2_unfiltered_spe@rowRanges@elementMetadata$Symbol)

# Remove the genes from the rowRanges
RML6_TS_2_unfiltered_spe@rowRanges <- RML6_TS_2_unfiltered_spe@rowRanges[-genes_to_remove, ]

# Update the counts accordingly
RML6_TS_2_unfiltered_spe@assays@data$counts <- RML6_TS_2_unfiltered_spe@assays@data$counts[-genes_to_remove, ]


#__________________________________________DECONVOLUTION OF INDIVIDUAL SPOTS______________________________________________________________
## this is the genes x barcode sparse count matrix
cd <- RML6_TS_2_unfiltered_spe@assays@data@listData$counts
pos <- SpatialExperiment::spatialCoords(RML6_TS_2_unfiltered_spe)

## change column names to x and y
## for this dataset, we will visualize barcodes using "pxl_col_in_fullres" = "y" coordinates, and "pxl_row_in_fullres" = "x" coordinates
colnames(pos) <- c("y", "x")

counts <- cleanCounts(cd, min.lib.size = 100, min.reads = 10)
corpus <- restrictCorpus(counts, removeAbove=1.0, removeBelow = 0.05, nTopOD = 1000)
ldas <- fitLDA(t(as.matrix(corpus)), Ks = seq(2, 8, by = 1),
               perc.rare.thresh = 0.05,
               plot=TRUE,
               verbose=TRUE)
optLDA <- optimalModel(models = ldas, opt = "min")

results <- getBetaTheta(optLDA, perc.filt = 0.05, betaScale = 1000)
deconProp <- results$theta
deconGexp <- results$beta


plt <- vizAllTopics(theta = deconProp,
                    pos = pos,
                    r = 45,
                    lwd = 0,
                    showLegend = TRUE,
                    plotTitle = NA) +
  ggplot2::guides(fill=ggplot2::guide_legend(ncol=2)) +

  ## outer border
  ggplot2::geom_rect(data = data.frame(pos),
                     ggplot2::aes(xmin = min(x)-90, xmax = max(x)+90,
                                  ymin = min(y)-90, ymax = max(y)+90),
                     fill = NA, color = "black", linetype = "solid", size = 0.5) +

  ggplot2::theme(
    plot.background = ggplot2::element_blank()
  ) +

  ## remove the pixel "groups", which is the color aesthetic for the pixel borders
  ggplot2::guides(colour = "none")


plt



ps <- lapply(colnames(deconProp), function(celltype) {

  vizTopic(theta = deconProp, pos = pos, topic = celltype, plotTitle = paste0("X", celltype),
           size = 0.8, stroke = 0.08, alpha = 0.8,
           low = "white",
           high = "blue") +

    ## remove the pixel "Groups", which is the color aesthetic for the pixel borders
    ggplot2::guides(colour = "none")

})
gridExtra::grid.arrange(
  grobs = ps,
  layout_matrix = rbind(c(1, 2, 3, 4),
                        c(5, 6, 7, 8),
                        c(9, 10, 11, 12),
                        c(13, 14, 15, 16))
)

# Assuming 'ps' is your list of plots
plot_names <- paste0("plot_", seq_along(ps), ".svg") # Creating file names for each plot

for (i in seq_along(ps)) {
  ggsave(plot_names[i], plot = ps[[i]], device = "svg", width = 2, height = 2, dpi = 300, path = "/Users/davidecaredio/spatial_prions/output/03st_deconvolution/Deconvolved_Topics/RML6_TS_2/Spatial_Expression") # Specify your path
}








geneSymbols <- RML6_TS_2_unfiltered_spe@rowRanges@elementMetadata$Symbol
names(geneSymbols) <- names(RML6_TS_2_unfiltered_spe@rowRanges)
geneSymbols[1:5]
colnames(deconGexp) <- geneSymbols[colnames(deconGexp)]
deconGexp[1:5,1:5]
ps <- lapply(colnames(deconProp), function(celltype) {

  celltype <- as.numeric(celltype)
  ## highly expressed in cell-type of interest
  highgexp <- names(which(deconGexp[celltype,] > 3))
  ## high log2(fold-change) compared to other deconvolved cell-types
  log2fc <- sort(log2(deconGexp[celltype,highgexp]/colMeans(deconGexp[-celltype,highgexp])), decreasing=TRUE)
  markers <- names(log2fc)[1] ## label just the top gene

  # -----------------------------------------------------
  ## visualize the transcriptional profile
  dat <- data.frame(values = as.vector(log2fc), genes = names(log2fc), order = seq(length(log2fc)))
  # Hide all of the text labels.
  dat$selectedLabels <- ""
  dat$selectedLabels[1] <- markers

  plt <- ggplot2::ggplot(data = dat) +
    ggplot2::geom_col(ggplot2::aes(x = order, y = values,
                                   fill = factor(selectedLabels == ""),
                                   color = factor(selectedLabels == "")), width = 1) +

    ggplot2::scale_fill_manual(values = c("darkblue",
                                          "darkblue"
    )) +
    ggplot2::scale_color_manual(values = c("darkblue",
                                           "darkblue"
    )) +

    ggplot2::scale_y_continuous(expand = c(0, 0), limits = c(min(log2fc) - 0.3, max(log2fc) + 0.3)) +
    # ggplot2::scale_x_continuous(expand = c(0, 0), limits = c(-2, NA)) +

    ggplot2::labs(title = paste0("X", celltype),
                  x = "Gene expression rank",
                  y = "log2(FC)") +

    ## placement of gene symbol labels of top genes
    ggplot2::geom_text(ggplot2::aes(x = order+1, y = values-0.1, label = selectedLabels), color = "red") +

    ggplot2::theme_classic() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size=6, color = "black"),
                   axis.text.y = ggplot2::element_text(size=6, color = "black"),
                   axis.title.y = ggplot2::element_text(size=6, color = "black"),
                   axis.title.x = ggplot2::element_text(size=6, color = "black"),
                   axis.ticks.x = ggplot2::element_blank(),
                   plot.title = ggplot2::element_text(size=8),
                   legend.text = ggplot2::element_text(size = 8, colour = "black"),
                   legend.title = ggplot2::element_text(size = 8, colour = "black", angle = 90),
                   panel.background = ggplot2::element_blank(),
                   plot.background = ggplot2::element_blank(),
                   panel.grid.major.y = ggplot2::element_line(size = 0.3, colour = "gray80"),
                   axis.line = ggplot2::element_line(size = 1, colour = "black"),
                   legend.position="none"
    )
  plt
})
gridExtra::grid.arrange(
  grobs = ps,
  layout_matrix = rbind(c(1, 2, 3, 4),
                        c(5, 6, 7, 8),
                        c(9, 10, 11, 12),
                        c(13, 14, 15, 16))
)


# Assuming 'ps' is your list of plots
plot_names <- paste0("plot_", seq_along(ps), ".svg") # Creating file names for each plot

for (i in seq_along(ps)) {
  ggsave(plot_names[i], plot = ps[[i]], device = "svg", width = 1, height = 2, dpi = 300, path = "/Users/davidecaredio/spatial_prions/output/03st_deconvolution/Deconvolved_Topics/RML6_TS_2/DEGs") # Specify your path
}


# Assuming you want to save each element in a separate file
for (i in seq_along(ps)) {
  saveRDS(ps[[i]][["data"]], file = paste0("output_", i, ".rds"))
}

# If you want to save all elements in a single file
all_data <- lapply(ps, function(x) x[["data"]])
# Assuming 'all_data' is your list of data frames
unique_labels <- unique(sapply(all_data, function(x) x$selectedLabels))
unique_labels <- unique_labels[unique_labels != ""]  # Exclude empty strings
names(all_data) <- sapply(all_data, function(x) x$selectedLabels[1])
output_file <- "Deconvolved_CellTypes/RML6_TS_2.rds"
saveRDS(all_data, file = output_file)






## first, combine the positions and the cleaned counts matrix
c <- counts
rownames(c) <- geneSymbols[rownames(c)]
df <- merge(as.data.frame(pos),
            as.data.frame(t(as.matrix(c))),
            by = 0)

## collect the top genes for subsequent visualization
markerGenes <- unlist(lapply(colnames(deconProp), function(celltype) {

  celltype <- as.numeric(celltype)
  ## highly expressed in cell-type of interest
  highgexp <- names(which(deconGexp[celltype,] > 3))
  ## high log2(fold-change) compared to other deconvolved cell-types
  log2fc <- sort(log2(deconGexp[celltype,highgexp]/colMeans(deconGexp[-celltype,highgexp])), decreasing=TRUE)
  markers <- names(log2fc)[1] ## label just the top gene
  ## collect name of top gene for each cell-type
  markers
}))


## now visualize top genes for each deconvolved cell-type
ps <- lapply(markerGenes, function(marker) {
  vizGeneCounts(df = df,
                gene = marker,
                # groups = annot,
                # group_cols = rainbow(length(levels(annot))),
                size = 2, stroke = 0.1,
                plotTitle = marker,
                winsorize = 0.05,
                showLegend = TRUE) +

    ## remove the pixel "groups", which is the color aesthetic for the pixel borders
    ggplot2::guides(colour = "none") +

    ## change some plot aesthetics
    ggplot2::theme(axis.text.x = ggplot2::element_text(size=0, color = "black", hjust = 0, vjust = 0.5),
                   axis.text.y = ggplot2::element_text(size=0, color = "black"),
                   axis.title.y = ggplot2::element_text(size=8),
                   axis.title.x = ggplot2::element_text(size=8),
                   plot.title = ggplot2::element_text(size=8),
                   legend.text = ggplot2::element_text(size = 8, colour = "black"),
                   legend.title = ggplot2::element_text(size = 8, colour = "black", angle = 90),
                   panel.background = ggplot2::element_blank(),
                   ## border around plot
                   panel.border = ggplot2::element_rect(fill = NA, color = "black", size = 1),
                   plot.background = ggplot2::element_blank()
    ) +
    ggplot2::guides(fill = ggplot2::guide_colorbar(title = "Counts",
                                                   title.position = "left",
                                                   title.hjust = 0.5,
                                                   ticks.colour = "black",
                                                   ticks.linewidth = 2,
                                                   frame.colour= "black",
                                                   frame.linewidth = 2,
                                                   label.hjust = 0
    ))
})
gridExtra::grid.arrange(
  grobs = ps,
  layout_matrix = rbind(c(1, 2, 3, 4),
                        c(5, 6, 7, 8),
                        c(9, 10, 11, 12),
                        c(13, 14, 15, 16))
)





#____________________________________COMPARE TO TRANSCRIPTIONAL CLUSTERING________________________________________________

#Calculate dimensionality reduction of the count matrix
pcs.info <- stats::prcomp(t(log10(as.matrix(counts) + 1)), center = TRUE, verbose = TRUE)
nPcs <- 7 ## let's take the top 5 PCs
pcs <- pcs.info$x[,1:nPcs]

#Generation of a 2D t-SNE embedding
emb <- Rtsne::Rtsne(pcs,
                    is_distance=FALSE,
                    perplexity=30,
                    num_threads=1,
                    verbose=FALSE)$Y
rownames(emb) <- rownames(pcs)
colnames(emb) <- c("x", "y")

# louvian clustering to assign the barcodes into 15 communities
k <- 35
com <- MERINGUE::getClusters(pcs, k, weight=TRUE, method = igraph::cluster_louvain)

#visualize the communities in terms of the spatial positions of the barcodes:
tempCom <- com
# Remove the last two rows from pos to match 2306 tempCom elements
pos <- pos[1:(nrow(pos) - 2), ]



dat <- data.frame("emb1" = pos[,"x"],
                  "emb2" = pos[,"y"],
                  "Cluster" = tempCom)

plt <- ggplot2::ggplot(data = dat) +
  ggplot2::geom_point(ggplot2::aes(x = emb1, y = emb2,
                                   color = Cluster), size = 0.8) +

  ggplot2::scale_color_manual(values = rainbow(n = length(levels(tempCom)))) +

  # ggplot2::scale_y_continuous(expand = c(0, 0), limits = c( min(dat$emb2)-1, max(dat$emb2)+1)) +
  # ggplot2::scale_x_continuous(expand = c(0, 0), limits = c( min(dat$emb1)-1, max(dat$emb1)+1) ) +

  ggplot2::labs(title = "",
                x = "x",
                y = "y") +

  ggplot2::theme_classic() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(size=15, color = "black"),
                 axis.text.y = ggplot2::element_text(size=15, color = "black"),
                 axis.title.y = ggplot2::element_text(size=15),
                 axis.title.x = ggplot2::element_text(size=15),
                 axis.ticks.x = ggplot2::element_blank(),
                 plot.title = ggplot2::element_text(size=15),
                 legend.text = ggplot2::element_text(size = 12, colour = "black"),
                 legend.title = ggplot2::element_text(size = 15, colour = "black", angle = 0, hjust = 0.5),
                 panel.background = ggplot2::element_blank(),
                 plot.background = ggplot2::element_blank(),
                 panel.grid.major.y =  ggplot2::element_blank(),
                 axis.line = ggplot2::element_line(size = 1, colour = "black")
                 # legend.position="none"
  ) +

  ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size=2), ncol = 2)
  ) +

  ggplot2::coord_equal()

plt




dat <- data.frame("emb1" = emb[,1],
                  "emb2" = emb[,2],
                  "Cluster" = tempCom)

## cluster labels
cent.pos <- do.call(rbind, tapply(1:nrow(emb), tempCom, function(ii) apply(emb[ii,,drop=F],2,median)))
cent.pos <- as.data.frame(cent.pos)
colnames(cent.pos) <- c("x", "y")
cent.pos$cluster <- rownames(cent.pos)
cent.pos <- na.omit(cent.pos)

plt <- ggplot2::ggplot(data = dat) +
  ggplot2::geom_point(ggplot2::aes(x = emb1, y = emb2,
                                   color = Cluster), size = 0.01) +

  ggplot2::scale_color_manual(values = rainbow(n = length(levels(tempCom)))) +

  ggplot2::scale_y_continuous(expand = c(0, 0), limits = c( min(dat$emb2)-1, max(dat$emb2)+1)) +
  ggplot2::scale_x_continuous(expand = c(0, 0), limits = c( min(dat$emb1)-1, max(dat$emb1)+1) ) +

  ggplot2::labs(title = "",
                x = "t-SNE 1",
                y = "t-SNE 2") +

  ggplot2::theme_classic() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(size=15, color = "black"),
                 axis.text.y = ggplot2::element_text(size=15, color = "black"),
                 axis.title.y = ggplot2::element_text(size=15),
                 axis.title.x = ggplot2::element_text(size=15),
                 axis.ticks.x = ggplot2::element_blank(),
                 plot.title = ggplot2::element_text(size=15),
                 legend.text = ggplot2::element_text(size = 12, colour = "black"),
                 legend.title = ggplot2::element_text(size = 15, colour = "black", angle = 0, hjust = 0.5),
                 panel.background = ggplot2::element_blank(),
                 plot.background = ggplot2::element_blank(),
                 panel.grid.major.y =  ggplot2::element_blank(),
                 axis.line = ggplot2::element_line(size = 1, colour = "black")
                 # legend.position="none"
  ) +

  ggplot2::geom_text(data = cent.pos,
                     ggplot2::aes(x = x,
                                  y = y,
                                  label = cluster),
                     fontface = "bold") +

  ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size=2), ncol = 2)
  ) +

  ggplot2::coord_equal()

plt



#Let’s see the proportions of each deconvolved cell-type across the embedding
ps <- lapply(colnames(deconProp), function(celltype) {

  vizTopic(theta = deconProp, pos = emb, topic = celltype, plotTitle = paste0("X", celltype),
           size = 1, stroke = 0.5, alpha = 0.5,
           low = "white",
           high = "red") +

    ## remove the pixel "Groups", which is the color aesthetic for the pixel borders
    ggplot2::guides(colour = "none")

})
gridExtra::grid.arrange(
  grobs = ps,
  layout_matrix = rbind(c(1, 2, 3, 4),
                        c(5, 6, 7, 8),
                        c(9, 10, 11, 12),
                        c(13, 14, 15, 16))
)


# let’s create a proxy “theta” matrix, which indicates the community each barcode was assigned to.
# proxy theta for the txn clusters
com_proxyTheta <- model.matrix(~ 0 + com)
rownames(com_proxyTheta) <- names(com)
# fix names
colnames(com_proxyTheta) <- unlist(lapply(colnames(com_proxyTheta), function(x) {
  unlist(strsplit(x, "com"))[2]
}))
com_proxyTheta <- as.data.frame.matrix(com_proxyTheta)
com_proxyTheta[1:5,1:5]


#Then we can build a correlation matrix of the correlations between the proportions of each cell-type and the transcriptional communities of the barcodes.
corMat_prop <- STdeconvolve::getCorrMtx(m1 = as.matrix(com_proxyTheta),
                                        m2 = deconProp,
                                        type = "t")
rownames(corMat_prop) <- paste0("com_", seq(nrow(corMat_prop)))
colnames(corMat_prop) <- paste0("decon_", seq(ncol(corMat_prop)))

## order the cell-types rows based on best match (highest correlation) with each community
corMat_prop_transposed <- t(corMat_prop)
pairs <- STdeconvolve::lsatPairs(corMat_prop_transposed)
m <- corMat_prop_transposed[pairs$rowix, pairs$colsix]

STdeconvolve::correlationPlot(mat = m,
                              colLabs = "STdeconvolve",
                              rowLabs = "Transcriptional clusters") +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 90)
  )




















#________________________GENES TO REMOVE: Hemoglobin and Mitochondrial related genes___________________________________________
RML6_TS_1_unfiltered_spe <- readRDS("RML6_TS_1_unfiltered_spe.rds")

# Define the pattern to remove
#Obteined from here: https://doi.org/10.1038/s41598-020-62801-6
pattern_to_remove <- "^mt-|^Hba1-|^Hba2-|^Hbb-|^Hbbp1-|^Hbd-|^Hbe1-|^Hbg1-|^Hbg2-|^Hmm-|^Hbq1-|^Hbz-|^Hbzp1- |^Hba-"

# Find the indices of genes that match the pattern
genes_to_remove <- grep(pattern_to_remove, RML6_TS_1_unfiltered_spe@rowRanges@elementMetadata$Symbol)

# Remove the genes from the rowRanges
RML6_TS_1_unfiltered_spe@rowRanges <- RML6_TS_1_unfiltered_spe@rowRanges[-genes_to_remove, ]

# Update the counts accordingly
RML6_TS_1_unfiltered_spe@assays@data$counts <- RML6_TS_1_unfiltered_spe@assays@data$counts[-genes_to_remove, ]


#__________________________________________DECONVOLUTION OF INDIVIDUAL SPOTS______________________________________________________________
## this is the genes x barcode sparse count matrix
cd <- RML6_TS_1_unfiltered_spe@assays@data@listData$counts
pos <- SpatialExperiment::spatialCoords(RML6_TS_1_unfiltered_spe)

## change column names to x and y
## for this dataset, we will visualize barcodes using "pxl_col_in_fullres" = "y" coordinates, and "pxl_row_in_fullres" = "x" coordinates
colnames(pos) <- c("y", "x")

counts <- cleanCounts(cd, min.lib.size = 100, min.reads = 10)
corpus <- restrictCorpus(counts, removeAbove=1.0, removeBelow = 0.05, nTopOD = 1000)
ldas <- fitLDA(t(as.matrix(corpus)), Ks = seq(2, 9, by = 1),
               perc.rare.thresh = 0.05,
               plot=TRUE,
               verbose=TRUE)
optLDA <- optimalModel(models = ldas, opt = "min")

results <- getBetaTheta(optLDA, perc.filt = 0.05, betaScale = 1000)
deconProp <- results$theta
deconGexp <- results$beta


plt <- vizAllTopics(theta = deconProp,
                    pos = pos,
                    r = 45,
                    lwd = 0,
                    showLegend = TRUE,
                    plotTitle = NA) +
  ggplot2::guides(fill=ggplot2::guide_legend(ncol=2)) +

  ## outer border
  ggplot2::geom_rect(data = data.frame(pos),
                     ggplot2::aes(xmin = min(x)-90, xmax = max(x)+90,
                                  ymin = min(y)-90, ymax = max(y)+90),
                     fill = NA, color = "black", linetype = "solid", size = 0.5) +

  ggplot2::theme(
    plot.background = ggplot2::element_blank()
  ) +

  ## remove the pixel "groups", which is the color aesthetic for the pixel borders
  ggplot2::guides(colour = "none")


plt



ps <- lapply(colnames(deconProp), function(celltype) {

  vizTopic(theta = deconProp, pos = pos, topic = celltype, plotTitle = paste0("X", celltype),
           size = 2, stroke = 1, alpha = 0.5,
           low = "white",
           high = "red") +

    ## remove the pixel "Groups", which is the color aesthetic for the pixel borders
    ggplot2::guides(colour = "none")

})
gridExtra::grid.arrange(
  grobs = ps,
  layout_matrix = rbind(c(1, 2, 3, 4),
                        c(5, 6, 7, 8),
                        c(9, 10, 11, 12),
                        c(13, 14, 15, 16))
)


geneSymbols <- RML6_TS_1_unfiltered_spe@rowRanges@elementMetadata$Symbol
names(geneSymbols) <- names(RML6_TS_1_unfiltered_spe@rowRanges)
geneSymbols[1:5]
colnames(deconGexp) <- geneSymbols[colnames(deconGexp)]
deconGexp[1:5,1:5]
ps <- lapply(colnames(deconProp), function(celltype) {

  celltype <- as.numeric(celltype)
  ## highly expressed in cell-type of interest
  highgexp <- names(which(deconGexp[celltype,] > 3))
  ## high log2(fold-change) compared to other deconvolved cell-types
  log2fc <- sort(log2(deconGexp[celltype,highgexp]/colMeans(deconGexp[-celltype,highgexp])), decreasing=TRUE)
  markers <- names(log2fc)[1] ## label just the top gene

  # -----------------------------------------------------
  ## visualize the transcriptional profile
  dat <- data.frame(values = as.vector(log2fc), genes = names(log2fc), order = seq(length(log2fc)))
  # Hide all of the text labels.
  dat$selectedLabels <- ""
  dat$selectedLabels[1] <- markers

  plt <- ggplot2::ggplot(data = dat) +
    ggplot2::geom_col(ggplot2::aes(x = order, y = values,
                                   fill = factor(selectedLabels == ""),
                                   color = factor(selectedLabels == "")), width = 1) +

    ggplot2::scale_fill_manual(values = c("darkblue",
                                          "darkblue"
    )) +
    ggplot2::scale_color_manual(values = c("darkblue",
                                           "darkblue"
    )) +

    ggplot2::scale_y_continuous(expand = c(0, 0), limits = c(min(log2fc) - 0.3, max(log2fc) + 0.3)) +
    # ggplot2::scale_x_continuous(expand = c(0, 0), limits = c(-2, NA)) +

    ggplot2::labs(title = paste0("X", celltype),
                  x = "Gene expression rank",
                  y = "log2(FC)") +

    ## placement of gene symbol labels of top genes
    ggplot2::geom_text(ggplot2::aes(x = order+1, y = values-0.1, label = selectedLabels), color = "red") +

    ggplot2::theme_classic() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size=15, color = "black"),
                   axis.text.y = ggplot2::element_text(size=15, color = "black"),
                   axis.title.y = ggplot2::element_text(size=15, color = "black"),
                   axis.title.x = ggplot2::element_text(size=15, color = "black"),
                   axis.ticks.x = ggplot2::element_blank(),
                   plot.title = ggplot2::element_text(size=15),
                   legend.text = ggplot2::element_text(size = 15, colour = "black"),
                   legend.title = ggplot2::element_text(size = 15, colour = "black", angle = 90),
                   panel.background = ggplot2::element_blank(),
                   plot.background = ggplot2::element_blank(),
                   panel.grid.major.y = ggplot2::element_line(size = 0.3, colour = "gray80"),
                   axis.line = ggplot2::element_line(size = 1, colour = "black"),
                   legend.position="none"
    )
  plt
})
gridExtra::grid.arrange(
  grobs = ps,
  layout_matrix = rbind(c(1, 2, 3, 4),
                        c(5, 6, 7, 8),
                        c(9, 10, 11, 12),
                        c(13, 14, 15, 16))
)


# If you want to save all elements in a single file
all_data <- lapply(ps, function(x) x[["data"]])
# Assuming 'all_data' is your list of data frames
unique_labels <- unique(sapply(all_data, function(x) x$selectedLabels))
unique_labels <- unique_labels[unique_labels != ""]  # Exclude empty strings
names(all_data) <- sapply(all_data, function(x) x$selectedLabels[1])

output_file <- "Deconvolved_CellTypes/RML6_TS_1.rds"
saveRDS(all_data, file = output_file)






## first, combine the positions and the cleaned counts matrix
c <- counts
rownames(c) <- geneSymbols[rownames(c)]
df <- merge(as.data.frame(pos),
            as.data.frame(t(as.matrix(c))),
            by = 0)

## collect the top genes for subsequent visualization
markerGenes <- unlist(lapply(colnames(deconProp), function(celltype) {

  celltype <- as.numeric(celltype)
  ## highly expressed in cell-type of interest
  highgexp <- names(which(deconGexp[celltype,] > 3))
  ## high log2(fold-change) compared to other deconvolved cell-types
  log2fc <- sort(log2(deconGexp[celltype,highgexp]/colMeans(deconGexp[-celltype,highgexp])), decreasing=TRUE)
  markers <- names(log2fc)[1] ## label just the top gene
  ## collect name of top gene for each cell-type
  markers
}))


## now visualize top genes for each deconvolved cell-type
ps <- lapply(markerGenes, function(marker) {
  vizGeneCounts(df = df,
                gene = marker,
                # groups = annot,
                # group_cols = rainbow(length(levels(annot))),
                size = 2, stroke = 0.1,
                plotTitle = marker,
                winsorize = 0.05,
                showLegend = TRUE) +

    ## remove the pixel "groups", which is the color aesthetic for the pixel borders
    ggplot2::guides(colour = "none") +

    ## change some plot aesthetics
    ggplot2::theme(axis.text.x = ggplot2::element_text(size=0, color = "black", hjust = 0, vjust = 0.5),
                   axis.text.y = ggplot2::element_text(size=0, color = "black"),
                   axis.title.y = ggplot2::element_text(size=15),
                   axis.title.x = ggplot2::element_text(size=15),
                   plot.title = ggplot2::element_text(size=15),
                   legend.text = ggplot2::element_text(size = 15, colour = "black"),
                   legend.title = ggplot2::element_text(size = 15, colour = "black", angle = 90),
                   panel.background = ggplot2::element_blank(),
                   ## border around plot
                   panel.border = ggplot2::element_rect(fill = NA, color = "black", size = 2),
                   plot.background = ggplot2::element_blank()
    ) +
    ggplot2::guides(fill = ggplot2::guide_colorbar(title = "Counts",
                                                   title.position = "left",
                                                   title.hjust = 0.5,
                                                   ticks.colour = "black",
                                                   ticks.linewidth = 2,
                                                   frame.colour= "black",
                                                   frame.linewidth = 2,
                                                   label.hjust = 0
    ))
})
gridExtra::grid.arrange(
  grobs = ps,
  layout_matrix = rbind(c(1, 2, 3, 4),
                        c(5, 6, 7, 8),
                        c(9, 10, 11, 12),
                        c(13, 14, 15, 16))
)



#____________________________________COMPARE TO TRANSCRIPTIONAL CLUSTERING________________________________________________

#Calculate dimensionality reduction of the count matrix
pcs.info <- stats::prcomp(t(log10(as.matrix(counts) + 1)), center = TRUE, verbose = TRUE)
nPcs <- 7 ## let's take the top 5 PCs
pcs <- pcs.info$x[,1:nPcs]

#Generation of a 2D t-SNE embedding
emb <- Rtsne::Rtsne(pcs,
                    is_distance=FALSE,
                    perplexity=30,
                    num_threads=1,
                    verbose=FALSE)$Y
rownames(emb) <- rownames(pcs)
colnames(emb) <- c("x", "y")

# louvian clustering to assign the barcodes into 15 communities
k <- 35
com <- MERINGUE::getClusters(pcs, k, weight=TRUE, method = igraph::cluster_louvain)

#visualize the communities in terms of the spatial positions of the barcodes:
tempCom <- com
# Remove the last two rows from pos to match 2306 tempCom elements
pos <- pos[1:(nrow(pos) - 4), ]



dat <- data.frame("emb1" = pos[,"x"],
                  "emb2" = pos[,"y"],
                  "Cluster" = tempCom)

plt <- ggplot2::ggplot(data = dat) +
  ggplot2::geom_point(ggplot2::aes(x = emb1, y = emb2,
                                   color = Cluster), size = 0.8) +

  ggplot2::scale_color_manual(values = rainbow(n = length(levels(tempCom)))) +

  # ggplot2::scale_y_continuous(expand = c(0, 0), limits = c( min(dat$emb2)-1, max(dat$emb2)+1)) +
  # ggplot2::scale_x_continuous(expand = c(0, 0), limits = c( min(dat$emb1)-1, max(dat$emb1)+1) ) +

  ggplot2::labs(title = "",
                x = "x",
                y = "y") +

  ggplot2::theme_classic() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(size=15, color = "black"),
                 axis.text.y = ggplot2::element_text(size=15, color = "black"),
                 axis.title.y = ggplot2::element_text(size=15),
                 axis.title.x = ggplot2::element_text(size=15),
                 axis.ticks.x = ggplot2::element_blank(),
                 plot.title = ggplot2::element_text(size=15),
                 legend.text = ggplot2::element_text(size = 12, colour = "black"),
                 legend.title = ggplot2::element_text(size = 15, colour = "black", angle = 0, hjust = 0.5),
                 panel.background = ggplot2::element_blank(),
                 plot.background = ggplot2::element_blank(),
                 panel.grid.major.y =  ggplot2::element_blank(),
                 axis.line = ggplot2::element_line(size = 1, colour = "black")
                 # legend.position="none"
  ) +

  ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size=2), ncol = 2)
  ) +

  ggplot2::coord_equal()

plt




dat <- data.frame("emb1" = emb[,1],
                  "emb2" = emb[,2],
                  "Cluster" = tempCom)

## cluster labels
cent.pos <- do.call(rbind, tapply(1:nrow(emb), tempCom, function(ii) apply(emb[ii,,drop=F],2,median)))
cent.pos <- as.data.frame(cent.pos)
colnames(cent.pos) <- c("x", "y")
cent.pos$cluster <- rownames(cent.pos)
cent.pos <- na.omit(cent.pos)

plt <- ggplot2::ggplot(data = dat) +
  ggplot2::geom_point(ggplot2::aes(x = emb1, y = emb2,
                                   color = Cluster), size = 0.01) +

  ggplot2::scale_color_manual(values = rainbow(n = length(levels(tempCom)))) +

  ggplot2::scale_y_continuous(expand = c(0, 0), limits = c( min(dat$emb2)-1, max(dat$emb2)+1)) +
  ggplot2::scale_x_continuous(expand = c(0, 0), limits = c( min(dat$emb1)-1, max(dat$emb1)+1) ) +

  ggplot2::labs(title = "",
                x = "t-SNE 1",
                y = "t-SNE 2") +

  ggplot2::theme_classic() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(size=15, color = "black"),
                 axis.text.y = ggplot2::element_text(size=15, color = "black"),
                 axis.title.y = ggplot2::element_text(size=15),
                 axis.title.x = ggplot2::element_text(size=15),
                 axis.ticks.x = ggplot2::element_blank(),
                 plot.title = ggplot2::element_text(size=15),
                 legend.text = ggplot2::element_text(size = 12, colour = "black"),
                 legend.title = ggplot2::element_text(size = 15, colour = "black", angle = 0, hjust = 0.5),
                 panel.background = ggplot2::element_blank(),
                 plot.background = ggplot2::element_blank(),
                 panel.grid.major.y =  ggplot2::element_blank(),
                 axis.line = ggplot2::element_line(size = 1, colour = "black")
                 # legend.position="none"
  ) +

  ggplot2::geom_text(data = cent.pos,
                     ggplot2::aes(x = x,
                                  y = y,
                                  label = cluster),
                     fontface = "bold") +

  ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size=2), ncol = 2)
  ) +

  ggplot2::coord_equal()

plt



#Let’s see the proportions of each deconvolved cell-type across the embedding
ps <- lapply(colnames(deconProp), function(celltype) {

  vizTopic(theta = deconProp, pos = emb, topic = celltype, plotTitle = paste0("X", celltype),
           size = 1, stroke = 0.5, alpha = 0.5,
           low = "white",
           high = "red") +

    ## remove the pixel "Groups", which is the color aesthetic for the pixel borders
    ggplot2::guides(colour = "none")

})
gridExtra::grid.arrange(
  grobs = ps,
  layout_matrix = rbind(c(1, 2, 3, 4),
                        c(5, 6, 7, 8),
                        c(9, 10, 11, 12),
                        c(13, 14, 15, 16))
)


# let’s create a proxy “theta” matrix, which indicates the community each barcode was assigned to.
# proxy theta for the txn clusters
com_proxyTheta <- model.matrix(~ 0 + com)
rownames(com_proxyTheta) <- names(com)
# fix names
colnames(com_proxyTheta) <- unlist(lapply(colnames(com_proxyTheta), function(x) {
  unlist(strsplit(x, "com"))[2]
}))
com_proxyTheta <- as.data.frame.matrix(com_proxyTheta)
com_proxyTheta[1:5,1:5]


#Then we can build a correlation matrix of the correlations between the proportions of each cell-type and the transcriptional communities of the barcodes.
corMat_prop <- STdeconvolve::getCorrMtx(m1 = as.matrix(com_proxyTheta),
                                        m2 = deconProp,
                                        type = "t")
rownames(corMat_prop) <- paste0("com_", seq(nrow(corMat_prop)))
colnames(corMat_prop) <- paste0("decon_", seq(ncol(corMat_prop)))

## order the cell-types rows based on best match (highest correlation) with each community
corMat_prop_transposed <- t(corMat_prop)
pairs <- STdeconvolve::lsatPairs(corMat_prop_transposed)
m <- corMat_prop_transposed[pairs$rowix, pairs$colsix]

STdeconvolve::correlationPlot(mat = m,
                              colLabs = "STdeconvolve",
                              rowLabs = "Transcriptional clusters") +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 90)
  )






































#________________________GENES TO REMOVE: Hemoglobin and Mitochondrial related genes___________________________________________
RML6_30_2_unfiltered_spe <- readRDS("RML6_30_2_unfiltered_spe.rds")

# Define the pattern to remove
#Obteined from here: https://doi.org/10.1038/s41598-020-62801-6
pattern_to_remove <- "^mt-|^Hba1-|^Hba2-|^Hbb-|^Hbbp1-|^Hbd-|^Hbe1-|^Hbg1-|^Hbg2-|^Hmm-|^Hbq1-|^Hbz-|^Hbzp1- |^Hba-"

# Find the indices of genes that match the pattern
genes_to_remove <- grep(pattern_to_remove, RML6_30_2_unfiltered_spe@rowRanges@elementMetadata$Symbol)

# Remove the genes from the rowRanges
RML6_30_2_unfiltered_spe@rowRanges <- RML6_30_2_unfiltered_spe@rowRanges[-genes_to_remove, ]

# Update the counts accordingly
RML6_30_2_unfiltered_spe@assays@data$counts <- RML6_30_2_unfiltered_spe@assays@data$counts[-genes_to_remove, ]


#__________________________________________DECONVOLUTION OF INDIVIDUAL SPOTS______________________________________________________________
## this is the genes x barcode sparse count matrix
cd <- RML6_30_2_unfiltered_spe@assays@data@listData$counts
pos <- SpatialExperiment::spatialCoords(RML6_30_2_unfiltered_spe)

## change column names to x and y
## for this dataset, we will visualize barcodes using "pxl_col_in_fullres" = "y" coordinates, and "pxl_row_in_fullres" = "x" coordinates
colnames(pos) <- c("y", "x")

counts <- cleanCounts(cd, min.lib.size = 100, min.reads = 10)
corpus <- restrictCorpus(counts, removeAbove=1.0, removeBelow = 0.05, nTopOD = 1000)
ldas <- fitLDA(t(as.matrix(corpus)), Ks = seq(2, 9, by = 1),
               perc.rare.thresh = 0.05,
               plot=TRUE,
               verbose=TRUE)
optLDA <- optimalModel(models = ldas, opt = "min")

results <- getBetaTheta(optLDA, perc.filt = 0.05, betaScale = 1000)
deconProp <- results$theta
deconGexp <- results$beta


plt <- vizAllTopics(theta = deconProp,
                    pos = pos,
                    r = 45,
                    lwd = 0,
                    showLegend = TRUE,
                    plotTitle = NA) +
  ggplot2::guides(fill=ggplot2::guide_legend(ncol=2)) +

  ## outer border
  ggplot2::geom_rect(data = data.frame(pos),
                     ggplot2::aes(xmin = min(x)-90, xmax = max(x)+90,
                                  ymin = min(y)-90, ymax = max(y)+90),
                     fill = NA, color = "black", linetype = "solid", size = 0.5) +

  ggplot2::theme(
    plot.background = ggplot2::element_blank()
  ) +

  ## remove the pixel "groups", which is the color aesthetic for the pixel borders
  ggplot2::guides(colour = "none")


plt



ps <- lapply(colnames(deconProp), function(celltype) {

  vizTopic(theta = deconProp, pos = pos, topic = celltype, plotTitle = paste0("X", celltype),
           size = 2, stroke = 1, alpha = 0.5,
           low = "white",
           high = "red") +

    ## remove the pixel "Groups", which is the color aesthetic for the pixel borders
    ggplot2::guides(colour = "none")

})
gridExtra::grid.arrange(
  grobs = ps,
  layout_matrix = rbind(c(1, 2, 3, 4),
                        c(5, 6, 7, 8),
                        c(9, 10, 11, 12),
                        c(13, 14, 15, 16))
)


geneSymbols <- RML6_30_2_unfiltered_spe@rowRanges@elementMetadata$Symbol
names(geneSymbols) <- names(RML6_30_2_unfiltered_spe@rowRanges)
geneSymbols[1:5]
colnames(deconGexp) <- geneSymbols[colnames(deconGexp)]
deconGexp[1:5,1:5]
ps <- lapply(colnames(deconProp), function(celltype) {

  celltype <- as.numeric(celltype)
  ## highly expressed in cell-type of interest
  highgexp <- names(which(deconGexp[celltype,] > 3))
  ## high log2(fold-change) compared to other deconvolved cell-types
  log2fc <- sort(log2(deconGexp[celltype,highgexp]/colMeans(deconGexp[-celltype,highgexp])), decreasing=TRUE)
  markers <- names(log2fc)[1] ## label just the top gene

  # -----------------------------------------------------
  ## visualize the transcriptional profile
  dat <- data.frame(values = as.vector(log2fc), genes = names(log2fc), order = seq(length(log2fc)))
  # Hide all of the text labels.
  dat$selectedLabels <- ""
  dat$selectedLabels[1] <- markers

  plt <- ggplot2::ggplot(data = dat) +
    ggplot2::geom_col(ggplot2::aes(x = order, y = values,
                                   fill = factor(selectedLabels == ""),
                                   color = factor(selectedLabels == "")), width = 1) +

    ggplot2::scale_fill_manual(values = c("darkblue",
                                          "darkblue"
    )) +
    ggplot2::scale_color_manual(values = c("darkblue",
                                           "darkblue"
    )) +

    ggplot2::scale_y_continuous(expand = c(0, 0), limits = c(min(log2fc) - 0.3, max(log2fc) + 0.3)) +
    # ggplot2::scale_x_continuous(expand = c(0, 0), limits = c(-2, NA)) +

    ggplot2::labs(title = paste0("X", celltype),
                  x = "Gene expression rank",
                  y = "log2(FC)") +

    ## placement of gene symbol labels of top genes
    ggplot2::geom_text(ggplot2::aes(x = order+1, y = values-0.1, label = selectedLabels), color = "red") +

    ggplot2::theme_classic() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size=15, color = "black"),
                   axis.text.y = ggplot2::element_text(size=15, color = "black"),
                   axis.title.y = ggplot2::element_text(size=15, color = "black"),
                   axis.title.x = ggplot2::element_text(size=15, color = "black"),
                   axis.ticks.x = ggplot2::element_blank(),
                   plot.title = ggplot2::element_text(size=15),
                   legend.text = ggplot2::element_text(size = 15, colour = "black"),
                   legend.title = ggplot2::element_text(size = 15, colour = "black", angle = 90),
                   panel.background = ggplot2::element_blank(),
                   plot.background = ggplot2::element_blank(),
                   panel.grid.major.y = ggplot2::element_line(size = 0.3, colour = "gray80"),
                   axis.line = ggplot2::element_line(size = 1, colour = "black"),
                   legend.position="none"
    )
  plt
})
gridExtra::grid.arrange(
  grobs = ps,
  layout_matrix = rbind(c(1, 2, 3, 4),
                        c(5, 6, 7, 8),
                        c(9, 10, 11, 12),
                        c(13, 14, 15, 16))
)







## first, combine the positions and the cleaned counts matrix
c <- counts
rownames(c) <- geneSymbols[rownames(c)]
df <- merge(as.data.frame(pos),
            as.data.frame(t(as.matrix(c))),
            by = 0)

## collect the top genes for subsequent visualization
markerGenes <- unlist(lapply(colnames(deconProp), function(celltype) {

  celltype <- as.numeric(celltype)
  ## highly expressed in cell-type of interest
  highgexp <- names(which(deconGexp[celltype,] > 3))
  ## high log2(fold-change) compared to other deconvolved cell-types
  log2fc <- sort(log2(deconGexp[celltype,highgexp]/colMeans(deconGexp[-celltype,highgexp])), decreasing=TRUE)
  markers <- names(log2fc)[1] ## label just the top gene
  ## collect name of top gene for each cell-type
  markers
}))


## now visualize top genes for each deconvolved cell-type
ps <- lapply(markerGenes, function(marker) {
  vizGeneCounts(df = df,
                gene = marker,
                # groups = annot,
                # group_cols = rainbow(length(levels(annot))),
                size = 2, stroke = 0.1,
                plotTitle = marker,
                winsorize = 0.05,
                showLegend = TRUE) +

    ## remove the pixel "groups", which is the color aesthetic for the pixel borders
    ggplot2::guides(colour = "none") +

    ## change some plot aesthetics
    ggplot2::theme(axis.text.x = ggplot2::element_text(size=0, color = "black", hjust = 0, vjust = 0.5),
                   axis.text.y = ggplot2::element_text(size=0, color = "black"),
                   axis.title.y = ggplot2::element_text(size=15),
                   axis.title.x = ggplot2::element_text(size=15),
                   plot.title = ggplot2::element_text(size=15),
                   legend.text = ggplot2::element_text(size = 15, colour = "black"),
                   legend.title = ggplot2::element_text(size = 15, colour = "black", angle = 90),
                   panel.background = ggplot2::element_blank(),
                   ## border around plot
                   panel.border = ggplot2::element_rect(fill = NA, color = "black", size = 2),
                   plot.background = ggplot2::element_blank()
    ) +
    ggplot2::guides(fill = ggplot2::guide_colorbar(title = "Counts",
                                                   title.position = "left",
                                                   title.hjust = 0.5,
                                                   ticks.colour = "black",
                                                   ticks.linewidth = 2,
                                                   frame.colour= "black",
                                                   frame.linewidth = 2,
                                                   label.hjust = 0
    ))
})
gridExtra::grid.arrange(
  grobs = ps,
  layout_matrix = rbind(c(1, 2, 3, 4),
                        c(5, 6, 7, 8),
                        c(9, 10, 11, 12),
                        c(13, 14, 15, 16))
)



#____________________________________COMPARE TO TRANSCRIPTIONAL CLUSTERING________________________________________________

#Calculate dimensionality reduction of the count matrix
pcs.info <- stats::prcomp(t(log10(as.matrix(counts) + 1)), center = TRUE, verbose = TRUE)
nPcs <- 7 ## let's take the top 5 PCs
pcs <- pcs.info$x[,1:nPcs]

#Generation of a 2D t-SNE embedding
emb <- Rtsne::Rtsne(pcs,
                    is_distance=FALSE,
                    perplexity=30,
                    num_threads=1,
                    verbose=FALSE)$Y
rownames(emb) <- rownames(pcs)
colnames(emb) <- c("x", "y")

# louvian clustering to assign the barcodes into 15 communities
k <- 35
com <- MERINGUE::getClusters(pcs, k, weight=TRUE, method = igraph::cluster_louvain)

#visualize the communities in terms of the spatial positions of the barcodes:
tempCom <- com
# Remove the last two rows from pos to match 2306 tempCom elements
pos <- pos[1:(nrow(pos) - 4), ]



dat <- data.frame("emb1" = pos[,"x"],
                  "emb2" = pos[,"y"],
                  "Cluster" = tempCom)

plt <- ggplot2::ggplot(data = dat) +
  ggplot2::geom_point(ggplot2::aes(x = emb1, y = emb2,
                                   color = Cluster), size = 0.8) +

  ggplot2::scale_color_manual(values = rainbow(n = length(levels(tempCom)))) +

  # ggplot2::scale_y_continuous(expand = c(0, 0), limits = c( min(dat$emb2)-1, max(dat$emb2)+1)) +
  # ggplot2::scale_x_continuous(expand = c(0, 0), limits = c( min(dat$emb1)-1, max(dat$emb1)+1) ) +

  ggplot2::labs(title = "",
                x = "x",
                y = "y") +

  ggplot2::theme_classic() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(size=15, color = "black"),
                 axis.text.y = ggplot2::element_text(size=15, color = "black"),
                 axis.title.y = ggplot2::element_text(size=15),
                 axis.title.x = ggplot2::element_text(size=15),
                 axis.ticks.x = ggplot2::element_blank(),
                 plot.title = ggplot2::element_text(size=15),
                 legend.text = ggplot2::element_text(size = 12, colour = "black"),
                 legend.title = ggplot2::element_text(size = 15, colour = "black", angle = 0, hjust = 0.5),
                 panel.background = ggplot2::element_blank(),
                 plot.background = ggplot2::element_blank(),
                 panel.grid.major.y =  ggplot2::element_blank(),
                 axis.line = ggplot2::element_line(size = 1, colour = "black")
                 # legend.position="none"
  ) +

  ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size=2), ncol = 2)
  ) +

  ggplot2::coord_equal()

plt




dat <- data.frame("emb1" = emb[,1],
                  "emb2" = emb[,2],
                  "Cluster" = tempCom)

## cluster labels
cent.pos <- do.call(rbind, tapply(1:nrow(emb), tempCom, function(ii) apply(emb[ii,,drop=F],2,median)))
cent.pos <- as.data.frame(cent.pos)
colnames(cent.pos) <- c("x", "y")
cent.pos$cluster <- rownames(cent.pos)
cent.pos <- na.omit(cent.pos)

plt <- ggplot2::ggplot(data = dat) +
  ggplot2::geom_point(ggplot2::aes(x = emb1, y = emb2,
                                   color = Cluster), size = 0.01) +

  ggplot2::scale_color_manual(values = rainbow(n = length(levels(tempCom)))) +

  ggplot2::scale_y_continuous(expand = c(0, 0), limits = c( min(dat$emb2)-1, max(dat$emb2)+1)) +
  ggplot2::scale_x_continuous(expand = c(0, 0), limits = c( min(dat$emb1)-1, max(dat$emb1)+1) ) +

  ggplot2::labs(title = "",
                x = "t-SNE 1",
                y = "t-SNE 2") +

  ggplot2::theme_classic() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(size=15, color = "black"),
                 axis.text.y = ggplot2::element_text(size=15, color = "black"),
                 axis.title.y = ggplot2::element_text(size=15),
                 axis.title.x = ggplot2::element_text(size=15),
                 axis.ticks.x = ggplot2::element_blank(),
                 plot.title = ggplot2::element_text(size=15),
                 legend.text = ggplot2::element_text(size = 12, colour = "black"),
                 legend.title = ggplot2::element_text(size = 15, colour = "black", angle = 0, hjust = 0.5),
                 panel.background = ggplot2::element_blank(),
                 plot.background = ggplot2::element_blank(),
                 panel.grid.major.y =  ggplot2::element_blank(),
                 axis.line = ggplot2::element_line(size = 1, colour = "black")
                 # legend.position="none"
  ) +

  ggplot2::geom_text(data = cent.pos,
                     ggplot2::aes(x = x,
                                  y = y,
                                  label = cluster),
                     fontface = "bold") +

  ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size=2), ncol = 2)
  ) +

  ggplot2::coord_equal()

plt



#Let’s see the proportions of each deconvolved cell-type across the embedding
ps <- lapply(colnames(deconProp), function(celltype) {

  vizTopic(theta = deconProp, pos = emb, topic = celltype, plotTitle = paste0("X", celltype),
           size = 1, stroke = 0.5, alpha = 0.5,
           low = "white",
           high = "red") +

    ## remove the pixel "Groups", which is the color aesthetic for the pixel borders
    ggplot2::guides(colour = "none")

})
gridExtra::grid.arrange(
  grobs = ps,
  layout_matrix = rbind(c(1, 2, 3, 4),
                        c(5, 6, 7, 8),
                        c(9, 10, 11, 12),
                        c(13, 14, 15, 16))
)


# let’s create a proxy “theta” matrix, which indicates the community each barcode was assigned to.
# proxy theta for the txn clusters
com_proxyTheta <- model.matrix(~ 0 + com)
rownames(com_proxyTheta) <- names(com)
# fix names
colnames(com_proxyTheta) <- unlist(lapply(colnames(com_proxyTheta), function(x) {
  unlist(strsplit(x, "com"))[2]
}))
com_proxyTheta <- as.data.frame.matrix(com_proxyTheta)
com_proxyTheta[1:5,1:5]


#Then we can build a correlation matrix of the correlations between the proportions of each cell-type and the transcriptional communities of the barcodes.
corMat_prop <- STdeconvolve::getCorrMtx(m1 = as.matrix(com_proxyTheta),
                                        m2 = deconProp,
                                        type = "t")
rownames(corMat_prop) <- paste0("com_", seq(nrow(corMat_prop)))
colnames(corMat_prop) <- paste0("decon_", seq(ncol(corMat_prop)))

## order the cell-types rows based on best match (highest correlation) with each community
corMat_prop_transposed <- t(corMat_prop)
pairs <- STdeconvolve::lsatPairs(corMat_prop_transposed)
m <- corMat_prop_transposed[pairs$rowix, pairs$colsix]

STdeconvolve::correlationPlot(mat = m,
                              colLabs = "STdeconvolve",
                              rowLabs = "Transcriptional clusters") +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 90)
  )

































#________________________GENES TO REMOVE: Hemoglobin and Mitochondrial related genes___________________________________________
RML6_30_1_unfiltered_spe <- readRDS("RML6_30_1_unfiltered_spe.rds")

# Define the pattern to remove
#Obteined from here: https://doi.org/10.1038/s41598-020-62801-6
pattern_to_remove <- "^mt-|^Hba1-|^Hba2-|^Hbb-|^Hbbp1-|^Hbd-|^Hbe1-|^Hbg1-|^Hbg2-|^Hmm-|^Hbq1-|^Hbz-|^Hbzp1- |^Hba-"

# Find the indices of genes that match the pattern
genes_to_remove <- grep(pattern_to_remove, RML6_30_1_unfiltered_spe@rowRanges@elementMetadata$Symbol)

# Remove the genes from the rowRanges
RML6_30_1_unfiltered_spe@rowRanges <- RML6_30_1_unfiltered_spe@rowRanges[-genes_to_remove, ]

# Update the counts accordingly
RML6_30_1_unfiltered_spe@assays@data$counts <- RML6_30_1_unfiltered_spe@assays@data$counts[-genes_to_remove, ]


#__________________________________________DECONVOLUTION OF INDIVIDUAL SPOTS______________________________________________________________
## this is the genes x barcode sparse count matrix
cd <- RML6_30_1_unfiltered_spe@assays@data@listData$counts
pos <- SpatialExperiment::spatialCoords(RML6_30_1_unfiltered_spe)

## change column names to x and y
## for this dataset, we will visualize barcodes using "pxl_col_in_fullres" = "y" coordinates, and "pxl_row_in_fullres" = "x" coordinates
colnames(pos) <- c("y", "x")

counts <- cleanCounts(cd, min.lib.size = 100, min.reads = 10)
corpus <- restrictCorpus(counts, removeAbove=1.0, removeBelow = 0.05, nTopOD = 1000)
ldas <- fitLDA(t(as.matrix(corpus)), Ks = seq(2, 9, by = 1),
               perc.rare.thresh = 0.05,
               plot=TRUE,
               verbose=TRUE)
optLDA <- optimalModel(models = ldas, opt = "min")

results <- getBetaTheta(optLDA, perc.filt = 0.05, betaScale = 1000)
deconProp <- results$theta
deconGexp <- results$beta


plt <- vizAllTopics(theta = deconProp,
                    pos = pos,
                    r = 45,
                    lwd = 0,
                    showLegend = TRUE,
                    plotTitle = NA) +
  ggplot2::guides(fill=ggplot2::guide_legend(ncol=2)) +

  ## outer border
  ggplot2::geom_rect(data = data.frame(pos),
                     ggplot2::aes(xmin = min(x)-90, xmax = max(x)+90,
                                  ymin = min(y)-90, ymax = max(y)+90),
                     fill = NA, color = "black", linetype = "solid", size = 0.5) +

  ggplot2::theme(
    plot.background = ggplot2::element_blank()
  ) +

  ## remove the pixel "groups", which is the color aesthetic for the pixel borders
  ggplot2::guides(colour = "none")


plt



ps <- lapply(colnames(deconProp), function(celltype) {

  vizTopic(theta = deconProp, pos = pos, topic = celltype, plotTitle = paste0("X", celltype),
           size = 2, stroke = 1, alpha = 0.5,
           low = "white",
           high = "red") +

    ## remove the pixel "Groups", which is the color aesthetic for the pixel borders
    ggplot2::guides(colour = "none")

})
gridExtra::grid.arrange(
  grobs = ps,
  layout_matrix = rbind(c(1, 2, 3, 4),
                        c(5, 6, 7, 8),
                        c(9, 10, 11, 12),
                        c(13, 14, 15, 16))
)


geneSymbols <- RML6_30_1_unfiltered_spe@rowRanges@elementMetadata$Symbol
names(geneSymbols) <- names(RML6_30_1_unfiltered_spe@rowRanges)
geneSymbols[1:5]
colnames(deconGexp) <- geneSymbols[colnames(deconGexp)]
deconGexp[1:5,1:5]
ps <- lapply(colnames(deconProp), function(celltype) {

  celltype <- as.numeric(celltype)
  ## highly expressed in cell-type of interest
  highgexp <- names(which(deconGexp[celltype,] > 3))
  ## high log2(fold-change) compared to other deconvolved cell-types
  log2fc <- sort(log2(deconGexp[celltype,highgexp]/colMeans(deconGexp[-celltype,highgexp])), decreasing=TRUE)
  markers <- names(log2fc)[1] ## label just the top gene

  # -----------------------------------------------------
  ## visualize the transcriptional profile
  dat <- data.frame(values = as.vector(log2fc), genes = names(log2fc), order = seq(length(log2fc)))
  # Hide all of the text labels.
  dat$selectedLabels <- ""
  dat$selectedLabels[1] <- markers

  plt <- ggplot2::ggplot(data = dat) +
    ggplot2::geom_col(ggplot2::aes(x = order, y = values,
                                   fill = factor(selectedLabels == ""),
                                   color = factor(selectedLabels == "")), width = 1) +

    ggplot2::scale_fill_manual(values = c("darkblue",
                                          "darkblue"
    )) +
    ggplot2::scale_color_manual(values = c("darkblue",
                                           "darkblue"
    )) +

    ggplot2::scale_y_continuous(expand = c(0, 0), limits = c(min(log2fc) - 0.3, max(log2fc) + 0.3)) +
    # ggplot2::scale_x_continuous(expand = c(0, 0), limits = c(-2, NA)) +

    ggplot2::labs(title = paste0("X", celltype),
                  x = "Gene expression rank",
                  y = "log2(FC)") +

    ## placement of gene symbol labels of top genes
    ggplot2::geom_text(ggplot2::aes(x = order+1, y = values-0.1, label = selectedLabels), color = "red") +

    ggplot2::theme_classic() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size=15, color = "black"),
                   axis.text.y = ggplot2::element_text(size=15, color = "black"),
                   axis.title.y = ggplot2::element_text(size=15, color = "black"),
                   axis.title.x = ggplot2::element_text(size=15, color = "black"),
                   axis.ticks.x = ggplot2::element_blank(),
                   plot.title = ggplot2::element_text(size=15),
                   legend.text = ggplot2::element_text(size = 15, colour = "black"),
                   legend.title = ggplot2::element_text(size = 15, colour = "black", angle = 90),
                   panel.background = ggplot2::element_blank(),
                   plot.background = ggplot2::element_blank(),
                   panel.grid.major.y = ggplot2::element_line(size = 0.3, colour = "gray80"),
                   axis.line = ggplot2::element_line(size = 1, colour = "black"),
                   legend.position="none"
    )
  plt
})
gridExtra::grid.arrange(
  grobs = ps,
  layout_matrix = rbind(c(1, 2, 3, 4),
                        c(5, 6, 7, 8),
                        c(9, 10, 11, 12),
                        c(13, 14, 15, 16))
)







## first, combine the positions and the cleaned counts matrix
c <- counts
rownames(c) <- geneSymbols[rownames(c)]
df <- merge(as.data.frame(pos),
            as.data.frame(t(as.matrix(c))),
            by = 0)

## collect the top genes for subsequent visualization
markerGenes <- unlist(lapply(colnames(deconProp), function(celltype) {

  celltype <- as.numeric(celltype)
  ## highly expressed in cell-type of interest
  highgexp <- names(which(deconGexp[celltype,] > 3))
  ## high log2(fold-change) compared to other deconvolved cell-types
  log2fc <- sort(log2(deconGexp[celltype,highgexp]/colMeans(deconGexp[-celltype,highgexp])), decreasing=TRUE)
  markers <- names(log2fc)[1] ## label just the top gene
  ## collect name of top gene for each cell-type
  markers
}))


## now visualize top genes for each deconvolved cell-type
ps <- lapply(markerGenes, function(marker) {
  vizGeneCounts(df = df,
                gene = marker,
                # groups = annot,
                # group_cols = rainbow(length(levels(annot))),
                size = 2, stroke = 0.1,
                plotTitle = marker,
                winsorize = 0.05,
                showLegend = TRUE) +

    ## remove the pixel "groups", which is the color aesthetic for the pixel borders
    ggplot2::guides(colour = "none") +

    ## change some plot aesthetics
    ggplot2::theme(axis.text.x = ggplot2::element_text(size=0, color = "black", hjust = 0, vjust = 0.5),
                   axis.text.y = ggplot2::element_text(size=0, color = "black"),
                   axis.title.y = ggplot2::element_text(size=15),
                   axis.title.x = ggplot2::element_text(size=15),
                   plot.title = ggplot2::element_text(size=15),
                   legend.text = ggplot2::element_text(size = 15, colour = "black"),
                   legend.title = ggplot2::element_text(size = 15, colour = "black", angle = 90),
                   panel.background = ggplot2::element_blank(),
                   ## border around plot
                   panel.border = ggplot2::element_rect(fill = NA, color = "black", size = 2),
                   plot.background = ggplot2::element_blank()
    ) +
    ggplot2::guides(fill = ggplot2::guide_colorbar(title = "Counts",
                                                   title.position = "left",
                                                   title.hjust = 0.5,
                                                   ticks.colour = "black",
                                                   ticks.linewidth = 2,
                                                   frame.colour= "black",
                                                   frame.linewidth = 2,
                                                   label.hjust = 0
    ))
})
gridExtra::grid.arrange(
  grobs = ps,
  layout_matrix = rbind(c(1, 2, 3, 4),
                        c(5, 6, 7, 8),
                        c(9, 10, 11, 12),
                        c(13, 14, 15, 16))
)



#____________________________________COMPARE TO TRANSCRIPTIONAL CLUSTERING________________________________________________

#Calculate dimensionality reduction of the count matrix
pcs.info <- stats::prcomp(t(log10(as.matrix(counts) + 1)), center = TRUE, verbose = TRUE)
nPcs <- 7 ## let's take the top 5 PCs
pcs <- pcs.info$x[,1:nPcs]

#Generation of a 2D t-SNE embedding
emb <- Rtsne::Rtsne(pcs,
                    is_distance=FALSE,
                    perplexity=30,
                    num_threads=1,
                    verbose=FALSE)$Y
rownames(emb) <- rownames(pcs)
colnames(emb) <- c("x", "y")

# louvian clustering to assign the barcodes into 15 communities
k <- 35
com <- MERINGUE::getClusters(pcs, k, weight=TRUE, method = igraph::cluster_louvain)

#visualize the communities in terms of the spatial positions of the barcodes:
tempCom <- com
# Remove the last two rows from pos to match 2306 tempCom elements
pos <- pos[1:(nrow(pos) - 2), ]



dat <- data.frame("emb1" = pos[,"x"],
                  "emb2" = pos[,"y"],
                  "Cluster" = tempCom)

plt <- ggplot2::ggplot(data = dat) +
  ggplot2::geom_point(ggplot2::aes(x = emb1, y = emb2,
                                   color = Cluster), size = 0.8) +

  ggplot2::scale_color_manual(values = rainbow(n = length(levels(tempCom)))) +

  # ggplot2::scale_y_continuous(expand = c(0, 0), limits = c( min(dat$emb2)-1, max(dat$emb2)+1)) +
  # ggplot2::scale_x_continuous(expand = c(0, 0), limits = c( min(dat$emb1)-1, max(dat$emb1)+1) ) +

  ggplot2::labs(title = "",
                x = "x",
                y = "y") +

  ggplot2::theme_classic() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(size=15, color = "black"),
                 axis.text.y = ggplot2::element_text(size=15, color = "black"),
                 axis.title.y = ggplot2::element_text(size=15),
                 axis.title.x = ggplot2::element_text(size=15),
                 axis.ticks.x = ggplot2::element_blank(),
                 plot.title = ggplot2::element_text(size=15),
                 legend.text = ggplot2::element_text(size = 12, colour = "black"),
                 legend.title = ggplot2::element_text(size = 15, colour = "black", angle = 0, hjust = 0.5),
                 panel.background = ggplot2::element_blank(),
                 plot.background = ggplot2::element_blank(),
                 panel.grid.major.y =  ggplot2::element_blank(),
                 axis.line = ggplot2::element_line(size = 1, colour = "black")
                 # legend.position="none"
  ) +

  ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size=2), ncol = 2)
  ) +

  ggplot2::coord_equal()

plt




dat <- data.frame("emb1" = emb[,1],
                  "emb2" = emb[,2],
                  "Cluster" = tempCom)

## cluster labels
cent.pos <- do.call(rbind, tapply(1:nrow(emb), tempCom, function(ii) apply(emb[ii,,drop=F],2,median)))
cent.pos <- as.data.frame(cent.pos)
colnames(cent.pos) <- c("x", "y")
cent.pos$cluster <- rownames(cent.pos)
cent.pos <- na.omit(cent.pos)

plt <- ggplot2::ggplot(data = dat) +
  ggplot2::geom_point(ggplot2::aes(x = emb1, y = emb2,
                                   color = Cluster), size = 0.01) +

  ggplot2::scale_color_manual(values = rainbow(n = length(levels(tempCom)))) +

  ggplot2::scale_y_continuous(expand = c(0, 0), limits = c( min(dat$emb2)-1, max(dat$emb2)+1)) +
  ggplot2::scale_x_continuous(expand = c(0, 0), limits = c( min(dat$emb1)-1, max(dat$emb1)+1) ) +

  ggplot2::labs(title = "",
                x = "t-SNE 1",
                y = "t-SNE 2") +

  ggplot2::theme_classic() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(size=15, color = "black"),
                 axis.text.y = ggplot2::element_text(size=15, color = "black"),
                 axis.title.y = ggplot2::element_text(size=15),
                 axis.title.x = ggplot2::element_text(size=15),
                 axis.ticks.x = ggplot2::element_blank(),
                 plot.title = ggplot2::element_text(size=15),
                 legend.text = ggplot2::element_text(size = 12, colour = "black"),
                 legend.title = ggplot2::element_text(size = 15, colour = "black", angle = 0, hjust = 0.5),
                 panel.background = ggplot2::element_blank(),
                 plot.background = ggplot2::element_blank(),
                 panel.grid.major.y =  ggplot2::element_blank(),
                 axis.line = ggplot2::element_line(size = 1, colour = "black")
                 # legend.position="none"
  ) +

  ggplot2::geom_text(data = cent.pos,
                     ggplot2::aes(x = x,
                                  y = y,
                                  label = cluster),
                     fontface = "bold") +

  ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size=2), ncol = 2)
  ) +

  ggplot2::coord_equal()

plt



#Let’s see the proportions of each deconvolved cell-type across the embedding
ps <- lapply(colnames(deconProp), function(celltype) {

  vizTopic(theta = deconProp, pos = emb, topic = celltype, plotTitle = paste0("X", celltype),
           size = 1, stroke = 0.5, alpha = 0.5,
           low = "white",
           high = "red") +

    ## remove the pixel "Groups", which is the color aesthetic for the pixel borders
    ggplot2::guides(colour = "none")

})
gridExtra::grid.arrange(
  grobs = ps,
  layout_matrix = rbind(c(1, 2, 3, 4),
                        c(5, 6, 7, 8),
                        c(9, 10, 11, 12),
                        c(13, 14, 15, 16))
)


# let’s create a proxy “theta” matrix, which indicates the community each barcode was assigned to.
# proxy theta for the txn clusters
com_proxyTheta <- model.matrix(~ 0 + com)
rownames(com_proxyTheta) <- names(com)
# fix names
colnames(com_proxyTheta) <- unlist(lapply(colnames(com_proxyTheta), function(x) {
  unlist(strsplit(x, "com"))[2]
}))
com_proxyTheta <- as.data.frame.matrix(com_proxyTheta)
com_proxyTheta[1:5,1:5]


#Then we can build a correlation matrix of the correlations between the proportions of each cell-type and the transcriptional communities of the barcodes.
corMat_prop <- STdeconvolve::getCorrMtx(m1 = as.matrix(com_proxyTheta),
                                        m2 = deconProp,
                                        type = "t")
rownames(corMat_prop) <- paste0("com_", seq(nrow(corMat_prop)))
colnames(corMat_prop) <- paste0("decon_", seq(ncol(corMat_prop)))

## order the cell-types rows based on best match (highest correlation) with each community
corMat_prop_transposed <- t(corMat_prop)
pairs <- STdeconvolve::lsatPairs(corMat_prop_transposed)
m <- corMat_prop_transposed[pairs$rowix, pairs$colsix]

STdeconvolve::correlationPlot(mat = m,
                              colLabs = "STdeconvolve",
                              rowLabs = "Transcriptional clusters") +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 90)
  )

























#__________________________________GENES TO REMOVE: Hemoglobin and Mitochondrial related genes___________________________________________
NBH_TS_2_unfiltered_spe <- readRDS("NBH_TS_2_unfiltered_spe.rds")

# Define the pattern to remove
#Obteined from here: https://doi.org/10.1038/s41598-020-62801-6
pattern_to_remove <- "^mt-|^Hba1-|^Hba2-|^Hbb-|^Hbbp1-|^Hbd-|^Hbe1-|^Hbg1-|^Hbg2-|^Hmm-|^Hbq1-|^Hbz-|^Hbzp1- |^Hba-"

# Find the indices of genes that match the pattern
genes_to_remove <- grep(pattern_to_remove, NBH_TS_2_unfiltered_spe@rowRanges@elementMetadata$Symbol)

# Remove the genes from the rowRanges
NBH_TS_2_unfiltered_spe@rowRanges <- NBH_TS_2_unfiltered_spe@rowRanges[-genes_to_remove, ]

# Update the counts accordingly
NBH_TS_2_unfiltered_spe@assays@data$counts <- NBH_TS_2_unfiltered_spe@assays@data$counts[-genes_to_remove, ]


#__________________________________________DECONVOLUTION OF INDIVIDUAL SPOTS______________________________________________________________
## this is the genes x barcode sparse count matrix
cd <- NBH_TS_2_unfiltered_spe@assays@data@listData$counts
pos <- SpatialExperiment::spatialCoords(NBH_TS_2_unfiltered_spe)

## change column names to x and y
## for this dataset, we will visualize barcodes using "pxl_col_in_fullres" = "y" coordinates, and "pxl_row_in_fullres" = "x" coordinates
colnames(pos) <- c("y", "x")

counts <- cleanCounts(cd, min.lib.size = 100, min.reads = 10)
corpus <- restrictCorpus(counts, removeAbove=1.0, removeBelow = 0.05, nTopOD = 1000)
ldas <- fitLDA(t(as.matrix(corpus)), Ks = seq(2, 9, by = 1),
               perc.rare.thresh = 0.05,
               plot=TRUE,
               verbose=TRUE)
optLDA <- optimalModel(models = ldas, opt = "min")

results <- getBetaTheta(optLDA, perc.filt = 0.05, betaScale = 1000)
deconProp <- results$theta
deconGexp <- results$beta


plt <- vizAllTopics(theta = deconProp,
                    pos = pos,
                    r = 45,
                    lwd = 0,
                    showLegend = TRUE,
                    plotTitle = NA) +
  ggplot2::guides(fill=ggplot2::guide_legend(ncol=2)) +

  ## outer border
  ggplot2::geom_rect(data = data.frame(pos),
                     ggplot2::aes(xmin = min(x)-90, xmax = max(x)+90,
                                  ymin = min(y)-90, ymax = max(y)+90),
                     fill = NA, color = "black", linetype = "solid", size = 0.5) +

  ggplot2::theme(
    plot.background = ggplot2::element_blank()
  ) +

  ## remove the pixel "groups", which is the color aesthetic for the pixel borders
  ggplot2::guides(colour = "none")


plt



ps <- lapply(colnames(deconProp), function(celltype) {

  vizTopic(theta = deconProp, pos = pos, topic = celltype, plotTitle = paste0("X", celltype),
           size = 2, stroke = 1, alpha = 0.5,
           low = "white",
           high = "red") +

    ## remove the pixel "Groups", which is the color aesthetic for the pixel borders
    ggplot2::guides(colour = "none")

})
gridExtra::grid.arrange(
  grobs = ps,
  layout_matrix = rbind(c(1, 2, 3, 4),
                        c(5, 6, 7, 8),
                        c(9, 10, 11, 12),
                        c(13, 14, 15, 16))
)


geneSymbols <- NBH_TS_2_unfiltered_spe@rowRanges@elementMetadata$Symbol
names(geneSymbols) <- names(NBH_TS_2_unfiltered_spe@rowRanges)
geneSymbols[1:5]
colnames(deconGexp) <- geneSymbols[colnames(deconGexp)]
deconGexp[1:5,1:5]
ps <- lapply(colnames(deconProp), function(celltype) {

  celltype <- as.numeric(celltype)
  ## highly expressed in cell-type of interest
  highgexp <- names(which(deconGexp[celltype,] > 3))
  ## high log2(fold-change) compared to other deconvolved cell-types
  log2fc <- sort(log2(deconGexp[celltype,highgexp]/colMeans(deconGexp[-celltype,highgexp])), decreasing=TRUE)
  markers <- names(log2fc)[1] ## label just the top gene

  # -----------------------------------------------------
  ## visualize the transcriptional profile
  dat <- data.frame(values = as.vector(log2fc), genes = names(log2fc), order = seq(length(log2fc)))
  # Hide all of the text labels.
  dat$selectedLabels <- ""
  dat$selectedLabels[1] <- markers

  plt <- ggplot2::ggplot(data = dat) +
    ggplot2::geom_col(ggplot2::aes(x = order, y = values,
                                   fill = factor(selectedLabels == ""),
                                   color = factor(selectedLabels == "")), width = 1) +

    ggplot2::scale_fill_manual(values = c("darkblue",
                                          "darkblue"
    )) +
    ggplot2::scale_color_manual(values = c("darkblue",
                                           "darkblue"
    )) +

    ggplot2::scale_y_continuous(expand = c(0, 0), limits = c(min(log2fc) - 0.3, max(log2fc) + 0.3)) +
    # ggplot2::scale_x_continuous(expand = c(0, 0), limits = c(-2, NA)) +

    ggplot2::labs(title = paste0("X", celltype),
                  x = "Gene expression rank",
                  y = "log2(FC)") +

    ## placement of gene symbol labels of top genes
    ggplot2::geom_text(ggplot2::aes(x = order+1, y = values-0.1, label = selectedLabels), color = "red") +

    ggplot2::theme_classic() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size=15, color = "black"),
                   axis.text.y = ggplot2::element_text(size=15, color = "black"),
                   axis.title.y = ggplot2::element_text(size=15, color = "black"),
                   axis.title.x = ggplot2::element_text(size=15, color = "black"),
                   axis.ticks.x = ggplot2::element_blank(),
                   plot.title = ggplot2::element_text(size=15),
                   legend.text = ggplot2::element_text(size = 15, colour = "black"),
                   legend.title = ggplot2::element_text(size = 15, colour = "black", angle = 90),
                   panel.background = ggplot2::element_blank(),
                   plot.background = ggplot2::element_blank(),
                   panel.grid.major.y = ggplot2::element_line(size = 0.3, colour = "gray80"),
                   axis.line = ggplot2::element_line(size = 1, colour = "black"),
                   legend.position="none"
    )
  plt
})
gridExtra::grid.arrange(
  grobs = ps,
  layout_matrix = rbind(c(1, 2, 3, 4),
                        c(5, 6, 7, 8),
                        c(9, 10, 11, 12),
                        c(13, 14, 15, 16))
)


# If you want to save all elements in a single file
all_data <- lapply(ps, function(x) x[["data"]])
# Assuming 'all_data' is your list of data frames
unique_labels <- unique(sapply(all_data, function(x) x$selectedLabels))
unique_labels <- unique_labels[unique_labels != ""]  # Exclude empty strings
names(all_data) <- sapply(all_data, function(x) x$selectedLabels[1])

output_file <- "Deconvolved_CellTypes/NBH_TS_2.rds"
saveRDS(all_data, file = output_file)





## first, combine the positions and the cleaned counts matrix
c <- counts
rownames(c) <- geneSymbols[rownames(c)]
df <- merge(as.data.frame(pos),
            as.data.frame(t(as.matrix(c))),
            by = 0)

## collect the top genes for subsequent visualization
markerGenes <- unlist(lapply(colnames(deconProp), function(celltype) {

  celltype <- as.numeric(celltype)
  ## highly expressed in cell-type of interest
  highgexp <- names(which(deconGexp[celltype,] > 3))
  ## high log2(fold-change) compared to other deconvolved cell-types
  log2fc <- sort(log2(deconGexp[celltype,highgexp]/colMeans(deconGexp[-celltype,highgexp])), decreasing=TRUE)
  markers <- names(log2fc)[1] ## label just the top gene
  ## collect name of top gene for each cell-type
  markers
}))


## now visualize top genes for each deconvolved cell-type
ps <- lapply(markerGenes, function(marker) {
  vizGeneCounts(df = df,
                gene = marker,
                # groups = annot,
                # group_cols = rainbow(length(levels(annot))),
                size = 2, stroke = 0.1,
                plotTitle = marker,
                winsorize = 0.05,
                showLegend = TRUE) +

    ## remove the pixel "groups", which is the color aesthetic for the pixel borders
    ggplot2::guides(colour = "none") +

    ## change some plot aesthetics
    ggplot2::theme(axis.text.x = ggplot2::element_text(size=0, color = "black", hjust = 0, vjust = 0.5),
                   axis.text.y = ggplot2::element_text(size=0, color = "black"),
                   axis.title.y = ggplot2::element_text(size=15),
                   axis.title.x = ggplot2::element_text(size=15),
                   plot.title = ggplot2::element_text(size=15),
                   legend.text = ggplot2::element_text(size = 15, colour = "black"),
                   legend.title = ggplot2::element_text(size = 15, colour = "black", angle = 90),
                   panel.background = ggplot2::element_blank(),
                   ## border around plot
                   panel.border = ggplot2::element_rect(fill = NA, color = "black", size = 2),
                   plot.background = ggplot2::element_blank()
    ) +
    ggplot2::guides(fill = ggplot2::guide_colorbar(title = "Counts",
                                                   title.position = "left",
                                                   title.hjust = 0.5,
                                                   ticks.colour = "black",
                                                   ticks.linewidth = 2,
                                                   frame.colour= "black",
                                                   frame.linewidth = 2,
                                                   label.hjust = 0
    ))
})
gridExtra::grid.arrange(
  grobs = ps,
  layout_matrix = rbind(c(1, 2, 3, 4),
                        c(5, 6, 7, 8),
                        c(9, 10, 11, 12),
                        c(13, 14, 15, 16))
)



#____________________________________COMPARE TO TRANSCRIPTIONAL CLUSTERING________________________________________________

#Calculate dimensionality reduction of the count matrix
pcs.info <- stats::prcomp(t(log10(as.matrix(counts) + 1)), center = TRUE, verbose = TRUE)
nPcs <- 7 ## let's take the top 5 PCs
pcs <- pcs.info$x[,1:nPcs]

#Generation of a 2D t-SNE embedding
emb <- Rtsne::Rtsne(pcs,
                    is_distance=FALSE,
                    perplexity=30,
                    num_threads=1,
                    verbose=FALSE)$Y
rownames(emb) <- rownames(pcs)
colnames(emb) <- c("x", "y")

# louvian clustering to assign the barcodes into 15 communities
k <- 35
com <- MERINGUE::getClusters(pcs, k, weight=TRUE, method = igraph::cluster_louvain)

#visualize the communities in terms of the spatial positions of the barcodes:
tempCom <- com
# Remove the last two rows from pos to match 2306 tempCom elements
pos <- pos[1:(nrow(pos) - 2), ]



dat <- data.frame("emb1" = pos[,"x"],
                  "emb2" = pos[,"y"],
                  "Cluster" = tempCom)

plt <- ggplot2::ggplot(data = dat) +
  ggplot2::geom_point(ggplot2::aes(x = emb1, y = emb2,
                                   color = Cluster), size = 0.8) +

  ggplot2::scale_color_manual(values = rainbow(n = length(levels(tempCom)))) +

  # ggplot2::scale_y_continuous(expand = c(0, 0), limits = c( min(dat$emb2)-1, max(dat$emb2)+1)) +
  # ggplot2::scale_x_continuous(expand = c(0, 0), limits = c( min(dat$emb1)-1, max(dat$emb1)+1) ) +

  ggplot2::labs(title = "",
                x = "x",
                y = "y") +

  ggplot2::theme_classic() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(size=15, color = "black"),
                 axis.text.y = ggplot2::element_text(size=15, color = "black"),
                 axis.title.y = ggplot2::element_text(size=15),
                 axis.title.x = ggplot2::element_text(size=15),
                 axis.ticks.x = ggplot2::element_blank(),
                 plot.title = ggplot2::element_text(size=15),
                 legend.text = ggplot2::element_text(size = 12, colour = "black"),
                 legend.title = ggplot2::element_text(size = 15, colour = "black", angle = 0, hjust = 0.5),
                 panel.background = ggplot2::element_blank(),
                 plot.background = ggplot2::element_blank(),
                 panel.grid.major.y =  ggplot2::element_blank(),
                 axis.line = ggplot2::element_line(size = 1, colour = "black")
                 # legend.position="none"
  ) +

  ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size=2), ncol = 2)
  ) +

  ggplot2::coord_equal()

plt




dat <- data.frame("emb1" = emb[,1],
                  "emb2" = emb[,2],
                  "Cluster" = tempCom)

## cluster labels
cent.pos <- do.call(rbind, tapply(1:nrow(emb), tempCom, function(ii) apply(emb[ii,,drop=F],2,median)))
cent.pos <- as.data.frame(cent.pos)
colnames(cent.pos) <- c("x", "y")
cent.pos$cluster <- rownames(cent.pos)
cent.pos <- na.omit(cent.pos)

plt <- ggplot2::ggplot(data = dat) +
  ggplot2::geom_point(ggplot2::aes(x = emb1, y = emb2,
                                   color = Cluster), size = 0.01) +

  ggplot2::scale_color_manual(values = rainbow(n = length(levels(tempCom)))) +

  ggplot2::scale_y_continuous(expand = c(0, 0), limits = c( min(dat$emb2)-1, max(dat$emb2)+1)) +
  ggplot2::scale_x_continuous(expand = c(0, 0), limits = c( min(dat$emb1)-1, max(dat$emb1)+1) ) +

  ggplot2::labs(title = "",
                x = "t-SNE 1",
                y = "t-SNE 2") +

  ggplot2::theme_classic() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(size=15, color = "black"),
                 axis.text.y = ggplot2::element_text(size=15, color = "black"),
                 axis.title.y = ggplot2::element_text(size=15),
                 axis.title.x = ggplot2::element_text(size=15),
                 axis.ticks.x = ggplot2::element_blank(),
                 plot.title = ggplot2::element_text(size=15),
                 legend.text = ggplot2::element_text(size = 12, colour = "black"),
                 legend.title = ggplot2::element_text(size = 15, colour = "black", angle = 0, hjust = 0.5),
                 panel.background = ggplot2::element_blank(),
                 plot.background = ggplot2::element_blank(),
                 panel.grid.major.y =  ggplot2::element_blank(),
                 axis.line = ggplot2::element_line(size = 1, colour = "black")
                 # legend.position="none"
  ) +

  ggplot2::geom_text(data = cent.pos,
                     ggplot2::aes(x = x,
                                  y = y,
                                  label = cluster),
                     fontface = "bold") +

  ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size=2), ncol = 2)
  ) +

  ggplot2::coord_equal()

plt



#Let’s see the proportions of each deconvolved cell-type across the embedding
ps <- lapply(colnames(deconProp), function(celltype) {

  vizTopic(theta = deconProp, pos = emb, topic = celltype, plotTitle = paste0("X", celltype),
           size = 1, stroke = 0.5, alpha = 0.5,
           low = "white",
           high = "red") +

    ## remove the pixel "Groups", which is the color aesthetic for the pixel borders
    ggplot2::guides(colour = "none")

})
gridExtra::grid.arrange(
  grobs = ps,
  layout_matrix = rbind(c(1, 2, 3, 4),
                        c(5, 6, 7, 8),
                        c(9, 10, 11, 12),
                        c(13, 14, 15, 16))
)


# let’s create a proxy “theta” matrix, which indicates the community each barcode was assigned to.
# proxy theta for the txn clusters
com_proxyTheta <- model.matrix(~ 0 + com)
rownames(com_proxyTheta) <- names(com)
# fix names
colnames(com_proxyTheta) <- unlist(lapply(colnames(com_proxyTheta), function(x) {
  unlist(strsplit(x, "com"))[2]
}))
com_proxyTheta <- as.data.frame.matrix(com_proxyTheta)
com_proxyTheta[1:5,1:5]


#Then we can build a correlation matrix of the correlations between the proportions of each cell-type and the transcriptional communities of the barcodes.
corMat_prop <- STdeconvolve::getCorrMtx(m1 = as.matrix(com_proxyTheta),
                                        m2 = deconProp,
                                        type = "t")
rownames(corMat_prop) <- paste0("com_", seq(nrow(corMat_prop)))
colnames(corMat_prop) <- paste0("decon_", seq(ncol(corMat_prop)))

## order the cell-types rows based on best match (highest correlation) with each community
corMat_prop_transposed <- t(corMat_prop)
pairs <- STdeconvolve::lsatPairs(corMat_prop_transposed)
m <- corMat_prop_transposed[pairs$rowix, pairs$colsix]

STdeconvolve::correlationPlot(mat = m,
                              colLabs = "STdeconvolve",
                              rowLabs = "Transcriptional clusters") +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 90)
  )



































#__________________________________GENES TO REMOVE: Hemoglobin and Mitochondrial related genes___________________________________________
NBH_TS_1_unfiltered_spe <- readRDS("NBH_TS_1_unfiltered_spe.rds")

# Define the pattern to remove
#Obteined from here: https://doi.org/10.1038/s41598-020-62801-6
pattern_to_remove <- "^mt-|^Hba1-|^Hba2-|^Hbb-|^Hbbp1-|^Hbd-|^Hbe1-|^Hbg1-|^Hbg2-|^Hmm-|^Hbq1-|^Hbz-|^Hbzp1- |^Hba-"

# Find the indices of genes that match the pattern
genes_to_remove <- grep(pattern_to_remove, NBH_TS_1_unfiltered_spe@rowRanges@elementMetadata$Symbol)

# Remove the genes from the rowRanges
NBH_TS_1_unfiltered_spe@rowRanges <- NBH_TS_1_unfiltered_spe@rowRanges[-genes_to_remove, ]

# Update the counts accordingly
NBH_TS_1_unfiltered_spe@assays@data$counts <- NBH_TS_1_unfiltered_spe@assays@data$counts[-genes_to_remove, ]


#__________________________________________DECONVOLUTION OF INDIVIDUAL SPOTS______________________________________________________________
## this is the genes x barcode sparse count matrix
cd <- NBH_TS_1_unfiltered_spe@assays@data@listData$counts
pos <- SpatialExperiment::spatialCoords(NBH_TS_1_unfiltered_spe)

## change column names to x and y
## for this dataset, we will visualize barcodes using "pxl_col_in_fullres" = "y" coordinates, and "pxl_row_in_fullres" = "x" coordinates
colnames(pos) <- c("y", "x")

counts <- cleanCounts(cd, min.lib.size = 100, min.reads = 10)
corpus <- restrictCorpus(counts, removeAbove=1.0, removeBelow = 0.05, nTopOD = 1000)
ldas <- fitLDA(t(as.matrix(corpus)), Ks = seq(2, 9, by = 1),
               perc.rare.thresh = 0.05,
               plot=TRUE,
               verbose=TRUE)
optLDA <- optimalModel(models = ldas, opt = "min")

results <- getBetaTheta(optLDA, perc.filt = 0.05, betaScale = 1000)
deconProp <- results$theta
deconGexp <- results$beta


plt <- vizAllTopics(theta = deconProp,
                    pos = pos,
                    r = 45,
                    lwd = 0,
                    showLegend = TRUE,
                    plotTitle = NA) +
  ggplot2::guides(fill=ggplot2::guide_legend(ncol=2)) +

  ## outer border
  ggplot2::geom_rect(data = data.frame(pos),
                     ggplot2::aes(xmin = min(x)-90, xmax = max(x)+90,
                                  ymin = min(y)-90, ymax = max(y)+90),
                     fill = NA, color = "black", linetype = "solid", size = 0.5) +

  ggplot2::theme(
    plot.background = ggplot2::element_blank()
  ) +

  ## remove the pixel "groups", which is the color aesthetic for the pixel borders
  ggplot2::guides(colour = "none")


plt



ps <- lapply(colnames(deconProp), function(celltype) {

  vizTopic(theta = deconProp, pos = pos, topic = celltype, plotTitle = paste0("X", celltype),
           size = 2, stroke = 1, alpha = 0.5,
           low = "white",
           high = "red") +

    ## remove the pixel "Groups", which is the color aesthetic for the pixel borders
    ggplot2::guides(colour = "none")

})
gridExtra::grid.arrange(
  grobs = ps,
  layout_matrix = rbind(c(1, 2, 3, 4),
                        c(5, 6, 7, 8),
                        c(9, 10, 11, 12),
                        c(13, 14, 15, 16))
)


geneSymbols <- NBH_TS_1_unfiltered_spe@rowRanges@elementMetadata$Symbol
names(geneSymbols) <- names(NBH_TS_1_unfiltered_spe@rowRanges)
geneSymbols[1:5]
colnames(deconGexp) <- geneSymbols[colnames(deconGexp)]
deconGexp[1:5,1:5]
ps <- lapply(colnames(deconProp), function(celltype) {

  celltype <- as.numeric(celltype)
  ## highly expressed in cell-type of interest
  highgexp <- names(which(deconGexp[celltype,] > 3))
  ## high log2(fold-change) compared to other deconvolved cell-types
  log2fc <- sort(log2(deconGexp[celltype,highgexp]/colMeans(deconGexp[-celltype,highgexp])), decreasing=TRUE)
  markers <- names(log2fc)[1] ## label just the top gene

  # -----------------------------------------------------
  ## visualize the transcriptional profile
  dat <- data.frame(values = as.vector(log2fc), genes = names(log2fc), order = seq(length(log2fc)))
  # Hide all of the text labels.
  dat$selectedLabels <- ""
  dat$selectedLabels[1] <- markers

  plt <- ggplot2::ggplot(data = dat) +
    ggplot2::geom_col(ggplot2::aes(x = order, y = values,
                                   fill = factor(selectedLabels == ""),
                                   color = factor(selectedLabels == "")), width = 1) +

    ggplot2::scale_fill_manual(values = c("darkblue",
                                          "darkblue"
    )) +
    ggplot2::scale_color_manual(values = c("darkblue",
                                           "darkblue"
    )) +

    ggplot2::scale_y_continuous(expand = c(0, 0), limits = c(min(log2fc) - 0.3, max(log2fc) + 0.3)) +
    # ggplot2::scale_x_continuous(expand = c(0, 0), limits = c(-2, NA)) +

    ggplot2::labs(title = paste0("X", celltype),
                  x = "Gene expression rank",
                  y = "log2(FC)") +

    ## placement of gene symbol labels of top genes
    ggplot2::geom_text(ggplot2::aes(x = order+1, y = values-0.1, label = selectedLabels), color = "red") +

    ggplot2::theme_classic() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size=15, color = "black"),
                   axis.text.y = ggplot2::element_text(size=15, color = "black"),
                   axis.title.y = ggplot2::element_text(size=15, color = "black"),
                   axis.title.x = ggplot2::element_text(size=15, color = "black"),
                   axis.ticks.x = ggplot2::element_blank(),
                   plot.title = ggplot2::element_text(size=15),
                   legend.text = ggplot2::element_text(size = 15, colour = "black"),
                   legend.title = ggplot2::element_text(size = 15, colour = "black", angle = 90),
                   panel.background = ggplot2::element_blank(),
                   plot.background = ggplot2::element_blank(),
                   panel.grid.major.y = ggplot2::element_line(size = 0.3, colour = "gray80"),
                   axis.line = ggplot2::element_line(size = 1, colour = "black"),
                   legend.position="none"
    )
  plt
})
gridExtra::grid.arrange(
  grobs = ps,
  layout_matrix = rbind(c(1, 2, 3, 4),
                        c(5, 6, 7, 8),
                        c(9, 10, 11, 12),
                        c(13, 14, 15, 16))
)


# If you want to save all elements in a single file
all_data <- lapply(ps, function(x) x[["data"]])
# Assuming 'all_data' is your list of data frames
unique_labels <- unique(sapply(all_data, function(x) x$selectedLabels))
unique_labels <- unique_labels[unique_labels != ""]  # Exclude empty strings
names(all_data) <- sapply(all_data, function(x) x$selectedLabels[1])

output_file <- "Deconvolved_CellTypes/NBH_TS_1.rds"
saveRDS(all_data, file = output_file)








## first, combine the positions and the cleaned counts matrix
c <- counts
rownames(c) <- geneSymbols[rownames(c)]
df <- merge(as.data.frame(pos),
            as.data.frame(t(as.matrix(c))),
            by = 0)

## collect the top genes for subsequent visualization
markerGenes <- unlist(lapply(colnames(deconProp), function(celltype) {

  celltype <- as.numeric(celltype)
  ## highly expressed in cell-type of interest
  highgexp <- names(which(deconGexp[celltype,] > 3))
  ## high log2(fold-change) compared to other deconvolved cell-types
  log2fc <- sort(log2(deconGexp[celltype,highgexp]/colMeans(deconGexp[-celltype,highgexp])), decreasing=TRUE)
  markers <- names(log2fc)[1] ## label just the top gene
  ## collect name of top gene for each cell-type
  markers
}))


## now visualize top genes for each deconvolved cell-type
ps <- lapply(markerGenes, function(marker) {
  vizGeneCounts(df = df,
                gene = marker,
                # groups = annot,
                # group_cols = rainbow(length(levels(annot))),
                size = 2, stroke = 0.1,
                plotTitle = marker,
                winsorize = 0.05,
                showLegend = TRUE) +

    ## remove the pixel "groups", which is the color aesthetic for the pixel borders
    ggplot2::guides(colour = "none") +

    ## change some plot aesthetics
    ggplot2::theme(axis.text.x = ggplot2::element_text(size=0, color = "black", hjust = 0, vjust = 0.5),
                   axis.text.y = ggplot2::element_text(size=0, color = "black"),
                   axis.title.y = ggplot2::element_text(size=15),
                   axis.title.x = ggplot2::element_text(size=15),
                   plot.title = ggplot2::element_text(size=15),
                   legend.text = ggplot2::element_text(size = 15, colour = "black"),
                   legend.title = ggplot2::element_text(size = 15, colour = "black", angle = 90),
                   panel.background = ggplot2::element_blank(),
                   ## border around plot
                   panel.border = ggplot2::element_rect(fill = NA, color = "black", size = 2),
                   plot.background = ggplot2::element_blank()
    ) +
    ggplot2::guides(fill = ggplot2::guide_colorbar(title = "Counts",
                                                   title.position = "left",
                                                   title.hjust = 0.5,
                                                   ticks.colour = "black",
                                                   ticks.linewidth = 2,
                                                   frame.colour= "black",
                                                   frame.linewidth = 2,
                                                   label.hjust = 0
    ))
})
gridExtra::grid.arrange(
  grobs = ps,
  layout_matrix = rbind(c(1, 2, 3, 4),
                        c(5, 6, 7, 8),
                        c(9, 10, 11, 12),
                        c(13, 14, 15, 16))
)



#____________________________________COMPARE TO TRANSCRIPTIONAL CLUSTERING________________________________________________

#Calculate dimensionality reduction of the count matrix
pcs.info <- stats::prcomp(t(log10(as.matrix(counts) + 1)), center = TRUE, verbose = TRUE)
nPcs <- 7 ## let's take the top 5 PCs
pcs <- pcs.info$x[,1:nPcs]

#Generation of a 2D t-SNE embedding
emb <- Rtsne::Rtsne(pcs,
                    is_distance=FALSE,
                    perplexity=30,
                    num_threads=1,
                    verbose=FALSE)$Y
rownames(emb) <- rownames(pcs)
colnames(emb) <- c("x", "y")

# louvian clustering to assign the barcodes into 15 communities
k <- 35
com <- MERINGUE::getClusters(pcs, k, weight=TRUE, method = igraph::cluster_louvain)

#visualize the communities in terms of the spatial positions of the barcodes:
tempCom <- com
# Remove the last two rows from pos to match 2306 tempCom elements
pos <- pos[1:(nrow(pos) - 2), ]



dat <- data.frame("emb1" = pos[,"x"],
                  "emb2" = pos[,"y"],
                  "Cluster" = tempCom)

plt <- ggplot2::ggplot(data = dat) +
  ggplot2::geom_point(ggplot2::aes(x = emb1, y = emb2,
                                   color = Cluster), size = 0.8) +

  ggplot2::scale_color_manual(values = rainbow(n = length(levels(tempCom)))) +

  # ggplot2::scale_y_continuous(expand = c(0, 0), limits = c( min(dat$emb2)-1, max(dat$emb2)+1)) +
  # ggplot2::scale_x_continuous(expand = c(0, 0), limits = c( min(dat$emb1)-1, max(dat$emb1)+1) ) +

  ggplot2::labs(title = "",
                x = "x",
                y = "y") +

  ggplot2::theme_classic() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(size=15, color = "black"),
                 axis.text.y = ggplot2::element_text(size=15, color = "black"),
                 axis.title.y = ggplot2::element_text(size=15),
                 axis.title.x = ggplot2::element_text(size=15),
                 axis.ticks.x = ggplot2::element_blank(),
                 plot.title = ggplot2::element_text(size=15),
                 legend.text = ggplot2::element_text(size = 12, colour = "black"),
                 legend.title = ggplot2::element_text(size = 15, colour = "black", angle = 0, hjust = 0.5),
                 panel.background = ggplot2::element_blank(),
                 plot.background = ggplot2::element_blank(),
                 panel.grid.major.y =  ggplot2::element_blank(),
                 axis.line = ggplot2::element_line(size = 1, colour = "black")
                 # legend.position="none"
  ) +

  ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size=2), ncol = 2)
  ) +

  ggplot2::coord_equal()

plt




dat <- data.frame("emb1" = emb[,1],
                  "emb2" = emb[,2],
                  "Cluster" = tempCom)

## cluster labels
cent.pos <- do.call(rbind, tapply(1:nrow(emb), tempCom, function(ii) apply(emb[ii,,drop=F],2,median)))
cent.pos <- as.data.frame(cent.pos)
colnames(cent.pos) <- c("x", "y")
cent.pos$cluster <- rownames(cent.pos)
cent.pos <- na.omit(cent.pos)

plt <- ggplot2::ggplot(data = dat) +
  ggplot2::geom_point(ggplot2::aes(x = emb1, y = emb2,
                                   color = Cluster), size = 0.01) +

  ggplot2::scale_color_manual(values = rainbow(n = length(levels(tempCom)))) +

  ggplot2::scale_y_continuous(expand = c(0, 0), limits = c( min(dat$emb2)-1, max(dat$emb2)+1)) +
  ggplot2::scale_x_continuous(expand = c(0, 0), limits = c( min(dat$emb1)-1, max(dat$emb1)+1) ) +

  ggplot2::labs(title = "",
                x = "t-SNE 1",
                y = "t-SNE 2") +

  ggplot2::theme_classic() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(size=15, color = "black"),
                 axis.text.y = ggplot2::element_text(size=15, color = "black"),
                 axis.title.y = ggplot2::element_text(size=15),
                 axis.title.x = ggplot2::element_text(size=15),
                 axis.ticks.x = ggplot2::element_blank(),
                 plot.title = ggplot2::element_text(size=15),
                 legend.text = ggplot2::element_text(size = 12, colour = "black"),
                 legend.title = ggplot2::element_text(size = 15, colour = "black", angle = 0, hjust = 0.5),
                 panel.background = ggplot2::element_blank(),
                 plot.background = ggplot2::element_blank(),
                 panel.grid.major.y =  ggplot2::element_blank(),
                 axis.line = ggplot2::element_line(size = 1, colour = "black")
                 # legend.position="none"
  ) +

  ggplot2::geom_text(data = cent.pos,
                     ggplot2::aes(x = x,
                                  y = y,
                                  label = cluster),
                     fontface = "bold") +

  ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size=2), ncol = 2)
  ) +

  ggplot2::coord_equal()

plt



#Let’s see the proportions of each deconvolved cell-type across the embedding
ps <- lapply(colnames(deconProp), function(celltype) {

  vizTopic(theta = deconProp, pos = emb, topic = celltype, plotTitle = paste0("X", celltype),
           size = 1, stroke = 0.5, alpha = 0.5,
           low = "white",
           high = "red") +

    ## remove the pixel "Groups", which is the color aesthetic for the pixel borders
    ggplot2::guides(colour = "none")

})
gridExtra::grid.arrange(
  grobs = ps,
  layout_matrix = rbind(c(1, 2, 3, 4),
                        c(5, 6, 7, 8),
                        c(9, 10, 11, 12),
                        c(13, 14, 15, 16))
)


# let’s create a proxy “theta” matrix, which indicates the community each barcode was assigned to.
# proxy theta for the txn clusters
com_proxyTheta <- model.matrix(~ 0 + com)
rownames(com_proxyTheta) <- names(com)
# fix names
colnames(com_proxyTheta) <- unlist(lapply(colnames(com_proxyTheta), function(x) {
  unlist(strsplit(x, "com"))[2]
}))
com_proxyTheta <- as.data.frame.matrix(com_proxyTheta)
com_proxyTheta[1:5,1:5]


#Then we can build a correlation matrix of the correlations between the proportions of each cell-type and the transcriptional communities of the barcodes.
corMat_prop <- STdeconvolve::getCorrMtx(m1 = as.matrix(com_proxyTheta),
                                        m2 = deconProp,
                                        type = "t")
rownames(corMat_prop) <- paste0("com_", seq(nrow(corMat_prop)))
colnames(corMat_prop) <- paste0("decon_", seq(ncol(corMat_prop)))

## order the cell-types rows based on best match (highest correlation) with each community
corMat_prop_transposed <- t(corMat_prop)
pairs <- STdeconvolve::lsatPairs(corMat_prop_transposed)
m <- corMat_prop_transposed[pairs$rowix, pairs$colsix]

STdeconvolve::correlationPlot(mat = m,
                              colLabs = "STdeconvolve",
                              rowLabs = "Transcriptional clusters") +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 90)
  )






















