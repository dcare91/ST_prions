# spot qc for individual slides

suppressPackageStartupMessages({

  library(tidyverse)
  library(Seurat)
  library(testit)
  library(Matrix)
  library(strex)
  library(ggh4x)
  library(reshape2)
})



# breaks for plots

log_brks = c(1,10, 100, 1000, 10000, 100000) %>% log10()
log_labs = c('1', rep('10', 5))
log_labs = parse(text = paste(log_labs,
                              c("", "", "^2", "^3", "^4", "^5"),
                              sep = ""))

mito_brks = c(0.05, 0.1, 0.15, 0.3, 0.5) %>% qlogis()
mito_labs = c("5", "10", "15", "30", "50")




# this will append qc metrics to metadata of spatial seurat object
get_spot_qc = function(seu, mt_prefix = 'mt-'){


# get mito genes (library size and number of detected features should already be in Seurat object)
  mito_genes = rownames(seu)[grepl(mt_prefix, rownames(seu))]

# check if there are 13 mito genes
  assert("Mito genes were not extracted properly",
         length(mito_genes) == 13)

# get mito pct
  counts = GetAssayData(seu, slot = 'counts', assay = 'Spatial')

  mito_sum = Matrix::colSums(counts[mito_genes, ])
  all_sum = Matrix::colSums(counts)

  mito_pct = mito_sum/all_sum
  mito_pct_logit = qlogis((mito_sum + 1)/(all_sum + 2))

  mito_df = data.frame(mito_pct = mito_pct,
                       mito_pct_logit = mito_pct_logit)
  rownames(mito_df) = names(mito_pct)

  # add mito pct to metadata$
  seu_w_qc = AddMetaData(seu, mito_df)

  return(seu_w_qc)

}







qc_violin_plot = function(qc_meta, thr_vec = NULL){

  # transform qc meta

  violin_df = qc_meta %>%
    mutate(log_library_size = log10(nCount_Spatial),
           log_feature_number = log10(nFeature_Spatial)) %>%
    dplyr::select(-all_of(c('nCount_Spatial', 'nFeature_Spatial', 'mito_pct', 'batch')))

  violin_df_melt = reshape2::melt(violin_df, value.name = 'qc_val')

  levels(violin_df_melt$variable) = c('% mito', 'library size', '# of features')


  # get medians, quantiles..
  qc_meds = violin_df_melt %>%
    group_by(sample_id, variable) %>%
    summarise(q50 = median(qc_val, na.rm = TRUE),
              q10 = quantile(qc_val, 0.1, na.rm = TRUE),
              q90 = quantile(qc_val, 0.9, na.rm = TRUE),
              q025 = quantile(qc_val , 0.025, na.rm = TRUE),
              q975 = quantile(qc_val, 0,975, na.rm = TRUE))

  qc_meds$sample_id = factor(qc_meds$sample_id, levels = unique(qc_meds$sample_id))


  # make plot

  p = ggplot(violin_df_melt, aes(y = qc_val, x = sample_id)) +
   geom_violin(kernel = 'rectangular', adjust = 0.1,color = NA,  scale = 'width', fill = 'grey60') +
    geom_smooth(data = qc_meds, aes(y = q50, x = as.numeric(sample_id)),
                method = 'loess',  span = 0.50, color = 'black', alpha = 0.5,
                inherit.aes = FALSE) +
    facet_grid(.~ variable, scales = 'free', space = 'free_y') +
    scale_x_discrete(breaks = levels(qc_meds$sample_id), drop = FALSE) +
    theme_classic() +
    facetted_pos_scales(
      y = list(variable == 'library size' ~ scale_y_continuous(breaks = log_brks, labels = log_labs),
               variable == '# of features' ~ scale_y_continuous(breaks = log_brks, labels = log_labs),
               variable == '% mito' ~ scale_y_continuous(breaks = mito_brks, labels = mito_labs))

    ) +
    coord_flip() +
    labs(x = NULL, y = NULL)


  # extract thresholds if they are specified
  if(!is.null(thr_vec)){

    # check if names of thr_vec are ok
    assert("thr_vec names are incorrect",
           all(names(thr_vec) %in% c('nCount_Spatial', 'nFeature_Spatial', 'mito_pct')))

    lib_size_thr = thr_vec['nCount_Spatial'] %>% log10()
    feature_thr = thr_vec['nFeature_Spatial'] %>% log10()
    mito_thr = thr_vec['mito_pct'] %>% qlogis()

    thr_df = data.frame(variable = factor(c('library size', '# of features', '% mito')),
                        threshold = c(lib_size_thr, feature_thr, mito_thr))

    # add thresholds to plot

    p = p + geom_hline(data = thr_df,
                       aes(yintercept = threshold),
                       linetype = 'dashed',
                       color = 'red',
                       linewidth = 0.5)
  }


return(p)


}




