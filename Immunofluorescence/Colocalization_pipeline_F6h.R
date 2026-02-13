library(tidyverse)
library(EBImage)
library(ggpubr)
library(tidyverse)
library(ggpubr)


data_dir <- getwd()
replicates <- 1:6
conditions <- c("sham", "stroke")

read_gray <- function(path) {
  img <- readImage(path)
  if (colorMode(img) != 0) img <- channel(img, "gray")
  img
}

# Morphological background subtraction:
# bg is estimated by a large opening (choose size bigger than objects)
bg_subtract <- function(img, brush_size = 35) {
  bg <- opening(img, makeBrush(brush_size, shape = "disc"))
  out <- img - bg
  out[out < 0] <- 0
  out
}

# Adaptive/local threshold -> binary
# w = window size (odd), offset tunes stringency
adaptive_bin <- function(img, w = 51, offset = 0.02) {
  # EBImage::thresh returns 0/1 image
  m <- thresh(img, w = w, h = w, offset = offset)
  m > 0
}

# Clean small speckles
clean_mask <- function(mask, min_size = 20) {
  mask <- opening(mask, makeBrush(3, "disc"))
  mask <- closing(mask, makeBrush(3, "disc"))
  # remove tiny components
  lab <- bwlabel(mask)
  tab <- table(lab)
  keep <- as.integer(names(tab)[tab >= min_size])
  mask & (lab %in% keep)
}

analyze_one <- function(rep, condition,
                        iba1_otsu = TRUE,
                        G_brush = 35, G_w = 51, G_offset = 0.02,
                        L_brush = 25, L_w = 51, L_offset = 0.02,
                        min_size = 20) {
  
  g_file <- file.path(data_dir, paste0("rep", rep, "_", condition, "_Gpnmb.tif"))
  i_file <- file.path(data_dir, paste0("rep", rep, "_", condition, "_Iba1.tif"))
  l_file <- file.path(data_dir, paste0("rep", rep, "_", condition, "_Lgals3.tif"))
  
  if (!all(file.exists(c(g_file, i_file, l_file)))) {
    message("Missing files for rep ", rep, " ", condition)
    return(NULL)
  }
  
  G <- read_gray(g_file)
  I <- read_gray(i_file)
  L <- read_gray(l_file)
  
  ## 1) Iba1 mask (reference compartment)
  if (iba1_otsu) {
    I_thr <- otsu(I)
    Ibin <- I > I_thr
  } else {
    Ibin <- adaptive_bin(I, w = 51, offset = 0.02)
    I_thr <- NA_real_
  }
  Ibin <- clean_mask(Ibin, min_size = min_size)
  
  ## 2) Gpnmb: background subtract + adaptive threshold
  Gbs <- bg_subtract(G, brush_size = G_brush)
  Gbin <- adaptive_bin(Gbs, w = G_w, offset = G_offset)
  Gbin <- clean_mask(Gbin, min_size = min_size)
  Gbin <- Gbin & Ibin
  
  ## 3) Lgals3: background subtract + adaptive threshold (works well for puncta)
  Lbs <- bg_subtract(L, brush_size = L_brush)
  Lbin <- adaptive_bin(Lbs, w = L_w, offset = L_offset)
  Lbin <- clean_mask(Lbin, min_size = min_size)
  Lbin <- Lbin & Ibin
  
  ## 4) Overlaps (inside Iba1)
  GL <- Gbin & Lbin
  triple <- Gbin & Lbin & Ibin  # same as GL since already gated to Ibin
  
  I_pix <- sum(Ibin)
  G_pix <- sum(Gbin)
  L_pix <- sum(Lbin)
  GL_pix <- sum(GL)
  
  tibble(
    Replicate = rep,
    Condition = condition,
    
    # record parameters (so you can report them)
    I_thr = I_thr,
    G_brush = G_brush, G_w = G_w, G_offset = G_offset,
    L_brush = L_brush, L_w = L_w, L_offset = L_offset,
    min_size = min_size,
    
    # counts within Iba1
    I_pixels = I_pix,
    G_in_I_pixels = G_pix,
    L_in_I_pixels = L_pix,
    GL_in_I_pixels = GL_pix,
    
    # the % overlaps (pick what you mean by "% overlap")
    pct_I_covered_by_GL = ifelse(I_pix > 0, 100 * GL_pix / I_pix, NA_real_),
    pct_G_covered_by_L  = ifelse(G_pix > 0, 100 * GL_pix / G_pix, NA_real_),
    pct_L_covered_by_G  = ifelse(L_pix > 0, 100 * GL_pix / L_pix, NA_real_)
  )
}

## Run batch
results <- expand_grid(Replicate = 1:6, Condition = c("sham","stroke")) %>%
  pmap_dfr(~ analyze_one(..1, ..2))

write_csv(results, file.path(data_dir, "overlap_Iba1_gated_adaptive.csv"))

## Plot (example: % of Iba1 area that is double-positive)
ggplot(results, aes(x = Condition, y = pct_I_covered_by_GL, fill = Condition)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.12, size = 2) +
  stat_compare_means(method = "wilcox.test", label = "p.format") +
  theme_classic() +
  ylab("% Iba1 area that is Gpnmb+ AND Lgals3+")




# ---- summarize mean ± SD + significance label ----
data_summary <- results %>%
  group_by(Condition) %>%
  summarise(
    Mean = mean(pct_I_covered_by_GL, na.rm = TRUE),
    StdDev = sd(pct_I_covered_by_GL, na.rm = TRUE),
    .groups = "drop"
  )

pval <- wilcox.test(pct_I_covered_by_GL ~ Condition, data = results)$p.value

sig_label <- dplyr::case_when(
  pval < 0.001 ~ "***",
  pval < 0.01  ~ "**",
  pval < 0.05  ~ "*",
  TRUE         ~ "ns"
)

data_summary <- data_summary %>%
  mutate(
    Significance = paste0("p = ", signif(pval, 2), " (", sig_label, ")"),
    y_text = Mean + StdDev + 0.05 * (max(results$pct_I_covered_by_GL, na.rm=TRUE) -
                                       min(results$pct_I_covered_by_GL, na.rm=TRUE))
  )

# ---- plot: grey bars (mean), points, errorbars, p-text ----
my_plot <- ggplot() +
  geom_bar(
    data = data_summary,
    aes(x = Condition, y = Mean, fill = Condition),
    stat = "identity",
    color = "black",
    width = 0.65
  ) +
  scale_fill_grey(start = 0.8, end = 0.2) +
  geom_errorbar(
    data = data_summary,
    aes(x = Condition, ymin = Mean - StdDev, ymax = Mean + StdDev),
    width = 0.18
  ) +
  geom_point(
    data = results,
    aes(x = Condition, y = pct_I_covered_by_GL, color = Condition),
    position = position_jitter(width = 0.12, height = 0),
    size = 2.2,
    show.legend = FALSE
  ) +
  geom_text(
    data = data_summary,
    aes(x = Condition, y = y_text, label = Significance),
    vjust = 0,
    size = 4
  ) +
  theme_minimal() +
  labs(
    x = "Condition",
    y = "% Iba1 area that is Gpnmb+ AND Lgals3+"
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 10, face = "bold"),
    panel.spacing = unit(1, "lines")
  )

my_plot



my_plot <- my_plot +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank()
  )

write_csv(results, file.path(data_dir, "overlap_Iba1_gated_adaptive.csv"))


# SVG (exact size: 54.5 x 18 mm)
ggsave(
  filename = "Gpnmb_Lgals3_overlap_boxplot.svg",
  plot = my_plot,
  width = 18,
  height = 54.5,
  units = "mm"
)


