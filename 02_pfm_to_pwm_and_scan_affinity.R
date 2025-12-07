
# File: src/R/02_pfm_to_pwm_and_scan_affinity.R
# Title: Convert PFMs to PWMs and scan promoter sequence to estimate binding affinity
# Description: Loads PFMs (from Step 01), converts to PWMs with log-odds (background 0.25 each base),
#              scans a promoter sequence, computes normalized affinity, and filters matches.
# Version: 1.0

# --- Dependencies ---
if (!requireNamespace("TFBSTools", quietly = TRUE)) install.packages("TFBSTools")
if (!requireNamespace("Biostrings", quietly = TRUE)) BiocManager::install("Biostrings")
if (!requireNamespace("GenomicRanges", quietly = TRUE)) BiocManager::install("GenomicRanges")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("readr", quietly = TRUE)) install.packages("readr")

library(TFBSTools)
library(Biostrings)
library(GenomicRanges)
library(dplyr)
library(readr)

# --- Parameters (edit as needed) ---
pfm_rds_file   <- file.path("data","external","PFMs_JASPAR_TRANSFAC_Hs.rds")
promoter_fasta <- file.path("data","raw","GHRL_promoter.fa")  # DO NOT upload real FASTA
use_simulated  <- TRUE   # Set to FALSE if you have a local FASTA (kept private)
set.seed(12345)

# Thresholds described in your methods
initial_score_threshold    <- 0.80  # initial score threshold
normalized_affinity_cutoff <- 0.85  # final normalized affinity cutoff

# Background probabilities (A,C,G,T)
background <- c(A=0.25, C=0.25, G=0.25, T=0.25)

# Output paths
out_dir <- file.path("data","processed")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_rds <- file.path(out_dir, "affinity_scan_results.rds")
out_csv <- file.path(out_dir, "affinity_scan_results.csv")

# --- Helper: convert PFM -> PWM (log-odds) ---
pfm_to_pwm_logodds <- function(pfm_matrix, bg = background, pseudocount = 0.8){
  freq <- sweep(pfm_matrix + pseudocount, 1, rowSums(pfm_matrix + pseudocount), FUN="/")
  pwm  <- log2(sweep(freq, 1, bg, FUN="/"))
  pwm
}

# --- Helper: score a k-mer with PWM ---
score_kmer <- function(kmer_chars, pwm){
  s <- 0
  for (pos in seq_len(ncol(pwm))) {
    base <- kmer_chars[pos]
    row_idx <- switch(base, A=1, C=2, G=3, T=4, NA)
    if (is.na(row_idx)) return(NA_real_)
    s <- s + pwm[row_idx, pos]
  }
  s
}

# --- Load PFMs ---
if (!file.exists(pfm_rds_file)) {
  stop("PFM RDS not found: ", pfm_rds_file, "
Run script 01 to generate it.")
}
pfm_list <- readRDS(pfm_rds_file)
if (length(pfm_list) == 0) stop("PFM list is empty. Check your PFM sources.")

# --- Load promoter sequence (simulate by default) ---
seq_obj <- NULL
if (!use_simulated && file.exists(promoter_fasta)) {
  dna <- readDNAStringSet(promoter_fasta)
  seq_obj <- dna[[1]]
  names(seq_obj) <- names(dna)[1]
} else {
  sim_seq <- DNAString(paste(sample(c("A","C","G","T"), size=1000, replace=TRUE), collapse=""))
  names(sim_seq) <- "GHRL_promoter_SIM"
  seq_obj <- sim_seq
}

# --- Scan sequence with each PWM, compute normalized affinity ---
results <- list()

for (i in seq_along(pfm_list)) {
  pfm <- pfm_list[[i]]
  mat <- as.matrix(Matrix(pfm))
  pwm <- pfm_to_pwm_logodds(mat, bg = background, pseudocount = 0.8)

  max_score <- sum(apply(pwm, 2, max))
  min_score <- sum(apply(pwm, 2, min))

  L <- ncol(pwm)
  seq_chars <- strsplit(as.character(seq_obj), "")[[1]]
  n_windows <- length(seq_chars) - L + 1
  if (n_windows <= 0) next

  scores <- numeric(n_windows)
  kmers  <- character(n_windows)

  for (start in seq_len(n_windows)) {
    kmer_chars <- seq_chars[start:(start+L-1)]
    kmers[start] <- paste0(kmer_chars, collapse="")
    scores[start] <- score_kmer(kmer_chars, pwm)
  }

  affinity_norm <- (scores - min_score) / (max_score - min_score)

  keep_idx <- which(affinity_norm >= initial_score_threshold)
  if (length(keep_idx) == 0) next

  df <- data.frame(
    tf            = TFBSTools::name(pfm),
    tf_id         = pfm@ID,
    start         = keep_idx,
    end           = keep_idx + L - 1,
    kmer          = kmers[keep_idx],
    score         = scores[keep_idx],
    affinity_norm = affinity_norm[keep_idx],
    stringsAsFactors = FALSE
  )

  df <- df[df$affinity_norm >= normalized_affinity_cutoff, , drop=FALSE]
  if (nrow(df) > 0) {
    results[[TFBSTools::name(pfm)]] <- df
  }
}

if (length(results) == 0) {
  message("No TFBS passed the normalized affinity cutoff (≥ ", normalized_affinity_cutoff, ").")
} else {
  combined <- bind_rows(results, .id = "tf_name")
  saveRDS(results, out_rds)
  write_csv(combined, out_csv)
  message("Saved RDS: ", out_rds)
  message("Saved CSV: ", out_csv)
}
