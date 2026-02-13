library(ijtiff)
library(EBImage)
library(coloc)
library(stringr)
library(ijtiff)
library(EBImage)
library(readxl)
library(dplyr)
library(ggplot2)

# Function to read and split channels from a TIFF file
read_and_split_channels <- function(file_path) {
  img <- read_tif(file_path)
  dims <- dim(img)
  
  # Print dimensions for debugging
  print(paste("Dimensions of", file_path, ":", paste(dims, collapse = " x ")))
  
  # Check if the image has the expected 4 channels
  if (length(dims) == 4 && dims[3] == 4) { # Case where channels are in the third dimension
    channels <- list(
      DAPI = img[,,1,1],
      Gpnmb = img[,,2,1],
      Lgals3 = img[,,3,1],
      Iba1 = img[,,4,1]
    )
    return(channels)
  } else if (length(dims) == 4 && dims[4] == 4) { # Case where channels are in the fourth dimension
    channels <- list(
      DAPI = img[,,,1],
      Gpnmb = img[,,,2],
      Lgals3 = img[,,,3],
      Iba1 = img[,,,4]
    )
    return(channels)
  } else {
    stop(paste("Error: The file", file_path, "does not contain 4 channels."))
  }
}

# Function to calculate Pearson correlation coefficient
calculate_pearson <- function(channel1, channel2) {
  cor(as.numeric(channel1), as.numeric(channel2))
}

# Function to perform three-channel colocalization analysis
three_channel_colocalization <- function(channel1, channel2, channel3, threshold = 0.5) {
  coloc_mask <- (channel1 > threshold) & (channel2 > threshold) & (channel3 > threshold)
  colocalized_area <- sum(coloc_mask)
  total_area <- length(channel1)
  colocalization_fraction <- colocalized_area / total_area
  return(colocalization_fraction)
}

# Function to perform colocalization analysis on one image
process_image <- function(file_path) {
  channels <- read_and_split_channels(file_path)
  
  # Perform pairwise colocalization analysis
  gpnmb_lgals3_pearson <- calculate_pearson(channels$Gpnmb, channels$Lgals3)
  gpnmb_iba1_pearson <- calculate_pearson(channels$Gpnmb, channels$Iba1)
  lgals3_iba1_pearson <- calculate_pearson(channels$Lgals3, channels$Iba1)
  
  # Three-channel colocalization analysis
  colocalization_fraction <- three_channel_colocalization(
    channels$Gpnmb, channels$Lgals3, channels$Iba1, threshold = 0.5
  )
  
  return(list(
    Gpnmb_Lgals3_Pearson = gpnmb_lgals3_pearson,
    Gpnmb_Iba1_Pearson = gpnmb_iba1_pearson,
    Lgals3_Iba1_Pearson = lgals3_iba1_pearson,
    Three_Channel_Colocalization_Fraction = colocalization_fraction
  ))
}

# Directory containing your TIFF files
input_dir <- "~/Desktop/New_Spatial_Analyses/IF_Gpnmb_Lgals3_Iba1/TIFF/"
file_list <- list.files(input_dir, pattern = "\\.tif$", full.names = TRUE)

# Initialize a list to store results
results <- list()

# Loop through each file and process
for (file_path in file_list) {
  print(paste("Processing file:", file_path))
  result <- tryCatch({
    process_image(file_path)
  }, error = function(e) {
    print(paste("Error processing file:", file_path, ":", e$message))
    NULL
  })
  results[[basename(file_path)]] <- result
}

# Filter out NULL results
results <- Filter(Negate(is.null), results)

# Convert list to a data frame ensuring atomic vectors
results_df <- do.call(rbind, lapply(names(results), function(name) {
  data.frame(Image = name, 
             Gpnmb_Lgals3_Pearson = results[[name]]$Gpnmb_Lgals3_Pearson, 
             Gpnmb_Iba1_Pearson = results[[name]]$Gpnmb_Iba1_Pearson, 
             Lgals3_Iba1_Pearson = results[[name]]$Lgals3_Iba1_Pearson, 
             Three_Channel_Colocalization_Fraction = results[[name]]$Three_Channel_Colocalization_Fraction)
}))

# Print the results dataframe with full precision for debugging
print(results_df, digits = 22)

# Save results to a CSV file
write.csv(results_df, file = "~/Desktop/Gpnmb_Lgals3_Iba1/colocalization_results.csv", row.names = FALSE)


# Load the Excel file
file_path <- "colocalization_results.xlsx"
coloc_data <- read_excel(file_path)

# Print the first few rows to understand the structure
head(coloc_data)

# Function to threshold and count positive cells
count_positive_cells <- function(channel1, channel2, channel3, threshold = 20) {
  mask1 <- channel1 > threshold
  mask2 <- channel2 > threshold
  mask3 <- channel3 > threshold
  
  # Logical AND operation to find colocalized cells
  coloc_mask <- mask1 & mask2 & mask3
  
  # Count the number of positive cells
  positive_cell_count <- sum(coloc_mask)
  return(positive_cell_count)
}

# Function to process each image and count positive cells
process_image_and_count_cells <- function(file_path, threshold = 20) {
  img <- read_tif(file_path)
  
  # Assuming channels are in the third or fourth dimension as discussed previously
  dims <- dim(img)
  if (length(dims) == 4 && dims[3] == 4) {
    channels <- list(
      DAPI = img[,,1,1],
      Gpnmb = img[,,2,1],
      Lgals3 = img[,,3,1],
      Iba1 = img[,,4,1]
    )
  } else if (length(dims) == 4 && dims[4] == 4) {
    channels <- list(
      DAPI = img[,,,1],
      Gpnmb = img[,,,2],
      Lgals3 = img[,,,3],
      Iba1 = img[,,,4]
    )
  } else {
    stop(paste("Error: The file", file_path, "does not contain 4 channels."))
  }
  
  # Count positive cells
  positive_cells <- count_positive_cells(channels$Gpnmb, channels$Lgals3, channels$Iba1, threshold)
  
  return(positive_cells)
}

# Directory containing your TIFF files
input_dir <- "~/Desktop/New_Spatial_Analyses/IF_Gpnmb_Lgals3_Iba1/TIFF/"
file_list <- list.files(input_dir, pattern = "\\.tif$", full.names = TRUE)

# Initialize a list to store results
results <- data.frame(Image = character(), Positive_Cell_Count = integer())

# Loop through each file and process
for (file_path in file_list) {
  print(paste("Processing file:", file_path))
  positive_cells <- tryCatch({
    process_image_and_count_cells(file_path)
  }, error = function(e) {
    print(paste("Error processing file:", file_path, ":", e$message))
    NA
  })
  
  # Extract file name
  file_name <- basename(file_path)
  
  results <- rbind(results, data.frame(Image = file_name, Positive_Cell_Count = positive_cells))
}

# Print the results
print(results)

# Merge the results with the metadata from the Excel file
merged_results <- merge(results, coloc_data, by="Image")

# Print the merged results
print(merged_results)

# Check for the presence of both groups in each region
valid_regions <- merged_results %>%
  group_by(Region) %>%
  filter(n_distinct(Strain) == 2) %>%
  pull(Region) %>%
  unique()

# Filter the results to include only valid regions
filtered_results <- merged_results %>%
  filter(Region %in% valid_regions)

# Perform a t-test to compare prions and NBH conditions for each valid region
comparison_results <- filtered_results %>%
  group_by(Region) %>%
  summarize(
    Prions_Mean = mean(Positive_Cell_Count[Strain == "Prions"], na.rm = TRUE),
    NBH_Mean = mean(Positive_Cell_Count[Strain == "NBH"], na.rm = TRUE),
    p_value = t.test(Positive_Cell_Count ~ Strain)$p.value
  )

# Print the results
print(comparison_results)

# Create the enhanced boxplot with specific colors
ggplot(filtered_results, aes(x = Region, y = Positive_Cell_Count, fill = Strain)) +
  geom_boxplot(outlier.shape = NA, alpha = 1, position = position_dodge(width = 0.75)) +  # Boxplot with no outliers shown
  geom_jitter(aes(color = Strain), position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), alpha = 0.6, size = 1.5) +  # Jittered points
  theme_minimal(base_size = 15) +  # Minimal theme with larger base font size
  labs(       x = "Region",
       y = "Gpnmb-Lgals3-Iba1+ Pixel Count") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),  # Centered title
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14)
  ) +
  scale_fill_manual(values = c("NBH" = "lightgrey", "Prions" = "azure4")) +  # Custom box colors
  scale_color_manual(values = c("NBH" = "cadetblue", "Prions" = "red"))  # Custom dot colors

# Save the plot
ggsave("positive_pixel_count_comparison_boxplot_custom_colors.svg", width = 8, height = 3.5)
