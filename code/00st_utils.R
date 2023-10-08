suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(Seurat)
  library(SpatialExperiment)
})

# make a seurat object from spaceranger outputs
# sranger_path has to contain subdirectories outs and outs/spatial
make_seurat_from_sranger <- function(sranger_path, sample, out_seurat_f, overwrite = FALSE){

  if(file.exists(out_seurat_f) & overwrite = FALSE){
   message("Reading spatial Seuart object for ", sample)
   readRDS(out_seurat_f)
  }else{

  # get spaceranger out dir
 sranger_outs <- file.path(sranger_path, 'outs')
 sranger_spatial <- file.path(sranger_outs, 'spatial')

 # get image
 img <- Read10X_Image(image.dir = sranger_spatial,
                      image.name = "tissue_lowres_image.png",
                      filter.matrix = TRUE)

 seu <- Load10X_Spatial(data.dir = sranger_outs,
                        filename = "filtered_feature_bc_matrix.h5",
                        assay = "Spatial",
                        slice = sample,
                        filter.matrix = TRUE,
                        to.upper = FALSE,
                        image = img)

 # change image coordinates from charater to integers
 img_coords <- seu@images[[sample]]@coordinates
 for(i in 1:ncol(img_coords)){
   img_coords[, i] <- as.integer(img_coords[, i])
 }

 seu@images[[sample]]@coordinates <- img_coords

 message("Saving spatial Seurat object for ", sample)
 saveRDS(seu, out_seurat_f)
 return(seu)
}
}



# make SpatialExperiment object from spaceranger outputs
make_spe_from_sranger <- function(sranger_path, sample, out_spe_f, overwrite = FALSE){

  if(file.exists(out_spe_f) & overwrite == FALSE){
    message("Reading spe object for ", sample)
    readRDS(out_spe_f)
  }else{
  sranger_outs <- file.path(sranger_path, 'outs')
  sranger_spatial <- file.path(sranger_outs, 'spatial')

  fmx_f <- file.path(sranger_outs, "filtered_feature_bc_matrix")
  sce <- DropletUtils::read10xCounts(fmx_f, col.names = TRUE)

  # read in image data
  img <- readImgData(
    path = file.path(sranger_spatial),
    sample_id = sample)

  # read in spatial coordinates
  spatial_coords_f <- file.path(sranger_spatial, "tissue_positions_list.csv")
  xyz <- read.csv(spatial_coords_f, header = TRUE,
                  col.names = c(
                    "barcode", "in_tissue", "array_row", "array_col",
                    "pxl_row_in_fullres", "pxl_col_in_fullres"))

  # filter to keep just barcodes in tissue
  xyz_filt = xyz %>% filter(barcode %in% colnames(sce))
  # arrange barcodes so they match the matrix
  xyz_filt = xyz_filt[match(colnames(sce), xyz_filt$barcode), ]

  # construct observation & feature metadata
  rd <- S4Vectors::DataFrame(
    symbol = rowData(sce)$Symbol)

  # construct 'SpatialExperiment'
  spe <- SpatialExperiment(
    assays = list(counts = assay(sce)),
    rowData = rd,
    colData = DataFrame(xyz),
    spatialCoordsNames = c("pxl_col_in_fullres", "pxl_row_in_fullres"),
    imgData = img,
    sample_id = sample)

  message("Saving spe object for sample ", sample)
  saveRDS(spe, out_spe_f)
  return(spe)

  }
}


