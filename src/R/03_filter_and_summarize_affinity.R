
# File: src/R/03_filter_and_summarize_affinity.R
# Title: Filter TFBS by normalized affinity ≥ 0.85 and summarize by promoter bins
# Description: Loads Step 02 results, re-applies the affinity cutoff if needed,
#              bins positions across the promoter, and computes per-bin summaries.
# Version: 1.0

# --- Dependencies ---
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("readr", quietly = TRUE)) install.packages("readr")
if (!requireNamespace("stringr", quietly = TRUE)) install.packages("stringr")

library(dplyr)
library(readr)
library(stringr)

# --- Parameters (edit as needed) ---
# Inputs produced in Step 02
in_rds <- file.path("data", "processed", "affinity_scan_results.rds")
in_csv <- file.path("data", "processed", "affinity_scan_results.csv")

# Final cutoff (as per your Methods)
normalized_affinity_cutoff <- 0.85

# Binning parameters
bin_size_bp <- 100      # bin width in base pairs
promoter_length <- 1000 # used to label bins; adjust if your sequence length differs

# Outputs
out_dir <- file.path("data", "processed")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_csv <- file.path(out_dir, "bins_tfbs_summary.csv")

# --- Load data (prefer RDS; fallback to CSV) ---
df <- NULL

if (file.exists(in_rds)) {
  res_list <- readRDS(in_rds)
  if (is.list(res_list) && length(res_list) > 0) {
    df <- dplyr::bind_rows(res_list, .id = "tf_name")
  }
}

if (is.null(df) && file.exists(in_csv)) {
  df <- readr::read_csv(in_csv, show_col_types = FALSE)
}

if (is.null(df) || nrow(df) == 0) {
  stop("No input data found. Please run Step 02 and ensure results exist in data/processed/.")
}

# --- Ensure required columns ---
needed_cols <- c("tf", "tf_id", "start", "end", "kmer", "score", "affinity_norm")
missing_cols <- setdiff(needed_cols, names(df))
if (length(missing_cols) > 0) {
  stop("Missing columns in input: ", paste(missing_cols, collapse = ", "))
}

# --- Re-apply final cutoff (defensive) ---
df_filt <- df %>% filter(affinity_norm >= normalized_affinity_cutoff)

if (nrow(df_filt) == 0) {
  warning("No TFBS entries meet the normalized affinity cutoff (≥ ", normalized_affinity_cutoff, ").")
  # Still write an empty summary with headers for reproducibility
  empty_summary <- tibble::tibble(
    bin = character(0),
    tfs = character(0),
    sum_affinity = numeric(0),
    mean_affinity = numeric(0),
    sites_ge_threshold = integer(0),
    density_per_kb = numeric(0),
    weighted_density_per_kb = numeric(0)
  )
  readr::write_csv(empty_summary, out_csv)
  message("Empty summary written to: ", out_csv)
  quit(save = "no")
}

# --- Create bins across the promoter ---
# If you used a simulated 1000 bp sequence in Step 02, positions start at 1..1000.
# Bin labeling: -1000 to -900, ..., -100 to 0 (relative scale; TSS assumed at 0 for labeling).
# For different promoter lengths, labels are scaled accordingly.

# Compute the bin index (0-based) from start coordinate
df_binned <- df_filt %>%
  mutate(
    bin_index = floor((start - 1) / bin_size_bp),
    bin_start_rel = (-promoter_length) + bin_index * bin_size_bp,
    bin_end_rel   = bin_start_rel + bin_size_bp,
    bin = paste0(bin_start_rel, " to ", bin_end_rel)
  )

# --- Summarize per bin ---
summary_bins <- df_binned %>%
  group_by(bin) %>%
  summarise(
    tfs = paste(unique(tf), collapse = ", "),
    sum_affinity = sum(affinity_norm),
    mean_affinity = mean(affinity_norm),
    sites_ge_threshold = dplyr::n(),
    density_per_kb = sites_ge_threshold / (bin_size_bp / 1000),
    weighted_density_per_kb = sum_affinity / (bin_size_bp / 1000),
    .groups = "drop"
  ) %>%
  arrange(bin)

# --- Save summary ---
readr::write_csv(summary_bins, out_csv)
message("Bin summary written to: ", out_csv)

# --- Optional: print a quick preview ---
print(utils::head(summary_bins, 10))
