

library(corrplot)
library(ggplot2)
library(ggrepel)


# Construct the full path to each RDS file
nbh_ts_1_path <- file.path("..", "NBH_TS_1.rds")
nbh_ts_2_path <- file.path("..", "NBH_TS_2.rds")
rml6_ts_1_path <- file.path("..", "RML6_TS_1.rds")
rml6_ts_2_path <- file.path("..", "RML6_TS_2.rds")


# Load the RDS files
nbh_ts_1 <- readRDS(nbh_ts_1_path)
nbh_ts_2 <- readRDS(nbh_ts_2_path)
rml6_ts_1 <- readRDS(rml6_ts_1_path)
rml6_ts_2 <- readRDS(rml6_ts_2_path)

# Combine all datasets into a single list
all_data <- list(nbh_ts_1, nbh_ts_2, rml6_ts_1, rml6_ts_2)

# Extract the 'values' column from each data frame
all_values <- lapply(all_data, function(dataset) {
  do.call(cbind, lapply(dataset, function(x) x$values))
})

# Find the minimum number of rows among all matrices
min_rows <- min(sapply(all_values, nrow))

# Trim each matrix to the minimum number of rows
all_values_trimmed <- lapply(all_values, function(mat) mat[1:min_rows, , drop = FALSE])

# Combine the trimmed matrices into a single data frame
all_data_combined <- as.data.frame(do.call(cbind, all_values_trimmed))

# Calculate the correlation matrix
cor_matrix_all <- cor(all_data_combined, use = "pairwise.complete.obs")


# Customize the appearance of the correlation plot
corrplot(cor_matrix_all, method = "color", col = colorRampPalette(c("blue", "white", "red"))(20),
         addCoef.col = "black", tl.cex = 0.8, cl.cex = 0.8, number.cex = 0.7)

# Add a title
title("Correlation Matrix Heatmap")



# Plot scatterplot  of the 21 overlapping genes of the 2 Gpnmb deconvoluted cell-type________________________________
#____________________________________________________________________________________________________________________
gpnmb_values_rml6_ts_1 <- rml6_ts_1$Gpnmb
gpnmb_values_rml6_ts_2 <- rml6_ts_2$Gpnmb

# Extract Gpnmb values and genes from each dataset
# Remove columns 3 and 4 from gpnmb_values_rml6_ts_1
gpnmb_values_rml6_ts_1 <- gpnmb_values_rml6_ts_1[, -c(3, 4)]
gpnmb_values_rml6_ts_2 <- gpnmb_values_rml6_ts_2[, -c(3, 4)]

# Merge data frames based on 'genes' column
merged_df <- merge(gpnmb_values_rml6_ts_1, gpnmb_values_rml6_ts_2, by = "genes", suffixes = c("_rml6_ts_1", "_rml6_ts_2"))

# Remove 'genes' column, leaving only the 'values' columns
values_df <- merged_df[, c("values_rml6_ts_1", "values_rml6_ts_2")]

# Calculate correlation
correlation_value <- cor(values_df$values_rml6_ts_1, values_df$values_rml6_ts_2)

# Print the correlation value
cat("Correlation between the two data frames:", correlation_value, "\n")


# Find common genes
common_genes <- intersect(gpnmb_values_rml6_ts_1$genes, gpnmb_values_rml6_ts_2$genes)

# Filter data for common genes
common_data_rml6_ts_1 <- gpnmb_values_rml6_ts_1[gpnmb_values_rml6_ts_1$genes %in% common_genes, ]
common_data_rml6_ts_2 <- gpnmb_values_rml6_ts_2[gpnmb_values_rml6_ts_2$genes %in% common_genes, ]

# Create a data frame for ggplot
plot_data <- data.frame(
  values_rml6_ts_1 = common_data_rml6_ts_1$values,
  values_rml6_ts_2 = common_data_rml6_ts_2$values
)

# Identify genes that need special label positions
special_label_genes <- c("Serpina3n", "C1qa")

# Create a scatterplot with ggplot2 and ggrepel
svg(filename = "/Users/davidecaredio/spatial_prions/output/03st_deconvolution/Deconvolved_CellTypes/Cell_Type_Correlation/venn_diagram.svg")
ggplot(plot_data, aes(x = values_rml6_ts_1, y = values_rml6_ts_2)) +
  geom_point(size = 3, color = "blue") +
  geom_text_repel(data = plot_data[1:21, ], aes(label = common_genes), box.padding = 0.5, size = 8) +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  labs(
    title = "Scatterplot of Gpnmb values",
    x = "Gpnmb_Rep_1",
    y = "Gpnmb_Rep_2",
    caption = paste("Correlation =", round(cor(plot_data$values_rml6_ts_1, plot_data$values_rml6_ts_2), 2))
  ) +
  theme_minimal() +
  theme(
    legend.position = "topright",
    text = element_text(size = 26)
  )

dev.off()


# Plot ORA results of the 21 overlapping genes of the 2 Gpnmb deconvoluted cell-type_________________________________
#____________________________________________________________________________________________________________________
library(clusterProfiler)
library(biomaRt)
library(org.Mm.eg.db)
library(clusterProfiler)
perform_ORA <- function(gene_list, organism = "mouse", database = "GO_Biological_Process_2021", ontology = "BP") {
  # Check if the input gene_list is not empty
  if (length(gene_list) == 0) {
    warning("Input gene list is empty. Returning NULL.")
    return(NULL)
  }

  # Specify the Ensembl dataset for mouse
  ensembl_mart <- useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")

  # Function to get ENSEMBL IDs for a list of genes
  get_ensembl_ids <- function(gene_list, mart) {
    # Ensure the gene list is in upper case for better matching
    gene_list <- toupper(gene_list)

    # Use the getBM function to retrieve ENSEMBL IDs
    result <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                    filters = "external_gene_name",
                    values = gene_list,
                    mart = mart)

    return(result)
  }

  ensembl_ids <- tryCatch(
    expr = {
      print(paste("Input gene_list:", paste(gene_list, collapse = ", ")))
      ensembl_ids <- get_ensembl_ids(gene_list, ensembl_mart)
      print(ensembl_ids)
      ensembl_ids  # Return ensembl_ids
    },
    error = function(e) {
      warning("Error in retrieving ENSEMBL IDs. Returning NULL.")
      return(NULL)
    },
    timeout = 300  # Set timeout to 300 seconds (adjust as needed)
  )

  # Check if any genes can be mapped
  if (is.null(ensembl_ids) || nrow(ensembl_ids) == 0) {
    warning("No genes can be mapped. Returning NULL.")
    return(NULL)
  }

  # Perform ORA
  enr <- enrichGO(gene = ensembl_ids$ensembl_gene_id, OrgDb = org.Mm.eg.db, keyType = "ENSEMBL", ont = ontology)

  return(enr)
}

ora_results <- perform_ORA(common_genes)
custom_theme <- list(
  axis.text.x = element_text(angle = 90, hjust = 1),  # Angle x-axis labels
  axis.text.y = element_text(size = 10)  # Adjust y-axis text size
)
ora_data <- as.data.frame(ora_results)
ora_data <- ora_results@result
write_xlsx(ora_data, "Supplementary_Table_7.xlsx")


# Plot the results using dotplot
go_plot <- dotplot(result,showCategory = 15)

# Customize x-axis labels and y-axis text size
go_plot <- go_plot +
    theme(axis.text.x = custom_theme$axis.text.x, axis.text.y = custom_theme$axis.text.y)

ggsave(filename, go_plot, width = 6, height = 6)
print(go_plot)


# Order the terms based on Count in decreasing order
ordered_terms <- ora_results[order(ora_results$Count, decreasing = TRUE), ]

# Remove specific descriptions
terms_to_remove <- c("humoral immune response mediated by circulating immunoglobulin", "negative regulation of smooth muscle cell proliferation")
filtered_terms <- ordered_terms[!ordered_terms$Description %in% terms_to_remove, ]

# Select the top 15 terms
top_terms <- head(filtered_terms, 15)

# Create a dotplot using ggplot2
go_plot <- ggplot(top_terms, aes(x = reorder(GeneRatio, Count), y = reorder(Description, Count), size = Count, color = pvalue)) +
    geom_point() +
    scale_size_continuous(range = c(2, 10)) +
    scale_color_continuous(trans = 'log10') +  # Adjust color scale
    labs(
      title = "Enrichment Analysis Results",
      x = "GeneRatio",
      y = "Description"
    ) +
    theme_minimal() +
    theme(
      axis.text.y = element_text(size = 10),
      axis.text.x = element_text(size = 10),
      legend.position = "bottom"
    )

  filename <- "plot.svg"
ggsave(filename, plot = go_plot, width = 4, height = 6)

# Display the plot
print(go_plot)


# VennDiagram of the 2 Gpnmb deconvoluted cell-type__________________________________________________________________
#____________________________________________________________________________________________________________________
library(VennDiagram)

# Create data frames for each set
set1 <- data.frame(values = gpnmb_values_rml6_ts_1$values, genes = gpnmb_values_rml6_ts_1$genes)
set2 <- data.frame(values = gpnmb_values_rml6_ts_2$values, genes = gpnmb_values_rml6_ts_2$genes)

# Merge the two sets based on genes
merged_data <- merge(set1, set2, by = "genes", suffixes = c("_set1", "_set2"), all = TRUE)

# Create a contingency table
contingency_table <- table(!is.na(merged_data$values_set1), !is.na(merged_data$values_set2))

# Perform Fisher's Exact Test
fisher_test_result <- fisher.test(contingency_table)
region_colors <- c("Gpnmb_Rep_1" = "magenta", "Gpnmb_Rep_2" = "cyan")
# Create a Venn diagram with enhanced aesthetics

svg(filename = "/Users/davidecaredio/spatial_prions/output/03st_deconvolution/Deconvolved_CellTypes/Cell_Type_Correlation/venn_diagram.svg")
venn.diagram(
  x = list(set1 = set1$genes, set2 = set2$genes),
  filename = NULL, # Important: do NOT specify a filename here
  main = "Gene intersection of Gpnmb deconvolved transcriptional profile in the two replicates",
  main.cex = 2,
  col = region_colors,
  fill = c(alpha("magenta", 0.3), alpha('cyan', 0.3)),
  cex = 2,
  fontfamily = "sans",
  cat.cex = 2,
  cat.col = c("magenta", 'cyan')
)

# Close the SVG device
dev.off()


# Display the Venn diagram
grid.draw(venn.plot)

# Display statistics including p-value
intersect_count <- nrow(merged_data[complete.cases(merged_data), ])
set1_count <- nrow(set1)
set2_count <- nrow(set2)

cat("Intersection Count:", intersect_count, "\n")
cat("Set 1 Count:", set1_count, "\n")
cat("Set 2 Count:", set2_count, "\n")
cat("Fisher's Exact Test p-value:", fisher_test_result$p.value, "\n")

dev.off() # Close the SVG device to finalize the file









