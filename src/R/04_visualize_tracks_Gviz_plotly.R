
# File: src/R/04_visualize_tracks_Gviz_plotly.R
# Title: Visualize promoter TFBS density (ggplot + plotly) and optional genome tracks (Gviz)
# Description: Loads bin summary (Step 03) and produces static and interactive plots.
# Version: 1.0

# --- Dependencies ---
if (!requireNamespace("readr", quietly = TRUE))  install.packages("readr")
if (!requireNamespace("dplyr", quietly = TRUE))  install.packages("dplyr")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("plotly", quietly = TRUE))  install.packages("plotly")
if (!requireNamespace("htmlwidgets", quietly = TRUE)) install.packages("htmlwidgets")
# Gviz optional:
# if (!requireNamespace("Gviz", quietly = TRUE)) BiocManager::install("Gviz")

library(readr)
library(dplyr)
library(ggplot2)
library(plotly)
library(htmlwidgets)
# library(Gviz)  # optional

# --- Inputs / Outputs ---
bins_csv <- file.path("data", "processed", "bins_tfbs_summary.csv")
out_dir  <- file.path("data", "processed")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

out_png_density   <- file.path(out_dir, "tfbs_density_per_bin.png")
out_html_density  <- file.path(out_dir, "tfbs_density_per_bin_interactive.html")
out_png_weighted  <- file.path(out_dir, "tfbs_weighted_density_per_bin.png")
out_html_weighted <- file.path(out_dir, "tfbs_weighted_density_per_bin_interactive.html")

# --- Load data ---
if (!file.exists(bins_csv)) {
  stop("Bin summary not found: ", bins_csv, "
Run Step 03 to generate it.")
}
bins <- readr::read_csv(bins_csv, show_col_types = FALSE)

if (nrow(bins) == 0) {
  stop("Bin summary is empty. Check Step 03 output.")
}

# --- Ensure required columns ---
needed_cols <- c("bin", "tfs", "sum_affinity", "mean_affinity",
                 "sites_ge_threshold", "density_per_kb", "weighted_density_per_kb")
missing_cols <- setdiff(needed_cols, names(bins))
if (length(missing_cols) > 0) {
  stop("Missing columns in bin summary: ", paste(missing_cols, collapse = ", "))
}

# --- Order bins logically ---
bins <- bins %>%
  mutate(bin_start = as.numeric(sub("^(.*) to .*$", "\1", bin))) %>%
  arrange(bin_start)

# --- Plot 1: Density per kb (static) ---
p_density <- ggplot(bins, aes(x = bin, y = density_per_kb)) +
  geom_col(fill = "#2C7FB8") +
  coord_flip() +
  labs(title = "TFBS density per promoter bin",
       x = "Promoter interval (bp)",
       y = "Density per kb") +
  theme_minimal(base_size = 12)

ggsave(out_png_density, p_density, width = 8, height = 5, dpi = 300)
message("Saved static density plot: ", out_png_density)

# --- Plot 1 (interactive) ---
p_density_i <- ggplotly(p_density, tooltip = c("x", "y")) %>%
  style(hoverlabel = list(bgcolor = "white"))
saveWidget(as_widget(p_density_i), file = out_html_density, selfcontained = TRUE)
message("Saved interactive density plot: ", out_html_density)

# --- Plot 2: Weighted density per kb (static) ---
p_weighted <- ggplot(bins, aes(x = bin, y = weighted_density_per_kb)) +
  geom_col(fill = "#41B6C4") +
  coord_flip() +
  labs(title = "Weighted TFBS density per promoter bin",
       x = "Promoter interval (bp)",
       y = "Weighted density per kb") +
  theme_minimal(base_size = 12)

ggsave(out_png_weighted, p_weighted, width = 8, height = 5, dpi = 300)
message("Saved static weighted plot: ", out_png_weighted)

# --- Plot 2 (interactive) ---
p_weighted_i <- ggplotly(p_weighted, tooltip = c("x", "y")) %>%
  style(hoverlabel = list(bgcolor = "white"))
saveWidget(as_widget(p_weighted_i), file = out_html_weighted, selfcontained = TRUE)
message("Saved interactive weighted plot: ", out_html_weighted)

# --- Optional: Gviz genome tracks template ---
# Uncomment and configure if you have real genome positions.
# library(Gviz)
# genome <- "hg38"
# chr    <- "chr3"
# from   <- 10289734
# to     <- 10290734
# ideoTrack <- IdeogramTrack(genome = genome, chromosome = chr)
# axisTrack <- GenomeAxisTrack()
# message("Gviz block is provided as a template. Enable it after mapping bins to genome coordinates.")
