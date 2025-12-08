
#!/usr/bin/env Rscript
# Functional enrichment for predicted TFs in the human GHRL promoter
# - GO (BP, MF, CC), KEGG, Reactome
# - Benjamini-Hochberg correction; FDR < 0.05
# - Input: file with gene symbols (one per line or a column named 'gene_symbol')
# - Outputs: CSV tables + PNG dotplots (disabled with --dry-run)

suppressPackageStartupMessages({
  library(optparse)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ReactomePA)
  library(dplyr)
})

option_list <- list(
  make_option(c("--genes"), type="character", help="Path to file with TF gene symbols (txt or csv)."),
  make_option(c("--col"), type="character", default="gene_symbol", help="Column name with symbols (if CSV). Default: 'gene_symbol'."),
  make_option(c("--outdir"), type="character", default="results/enrichment", help="Output directory for CSV/PNG."),
  make_option(c("--qcut"), type="double", default=0.05, help="FDR cutoff (q-value). Default: 0.05"),
  make_option(c("--pAdj"), type="character", default="BH", help="Multiple testing correction. Default: BH"),
  make_option(c("--dry-run"), action="store_true", default=FALSE, help="Run without writing files.")
)
opt <- parse_args(OptionParser(option_list=option_list))

if (is.null(opt$genes)) {
  stop("Please provide --genes with the path to a file containing TF gene symbols.")
}

# -------- Load gene symbols --------
read_symbols <- function(path, colname="gene_symbol") {
  ext <- tools::file_ext(path)
  if (tolower(ext) == "csv") {
    df <- read.csv(path, stringsAsFactors = FALSE)
    if (!colname %in% colnames(df)) {
      stop(sprintf("Column '%s' not found in CSV. Available: %s", colname, paste(colnames(df), collapse=", ")))
    }
    syms <- df[[colname]]
  } else {
    syms <- readLines(path, warn = FALSE)
  }
  syms <- toupper(trimws(syms))
  syms <- syms[syms != "" & !is.na(syms)]
  unique(syms)
}

symbols <- read_symbols(opt$genes, opt$col)
if (length(symbols) == 0) stop("No gene symbols found.")

message(sprintf("Loaded %d unique gene symbols.", length(symbols)))

# -------- Map symbols -> Entrez IDs --------
# Use bitr for conversion; drop NAs/duplicates
map_df <- bitr(symbols, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
entrez <- unique(map_df$ENTREZID)
message(sprintf("Mapped %d/%d symbols to Entrez IDs.", length(entrez), length(symbols)))

if (length(entrez) < 3) {
  stop("Fewer than 3 Entrez IDs after mapping; enrichment may be unreliable.")
}

# -------- Helper: save result table --------
save_enrich <- function(ego, outfile) {
  if (is.null(ego) || nrow(as.data.frame(ego)) == 0) return(FALSE)
  df <- as.data.frame(ego)
  df <- df %>%
    arrange(p.adjust) %>%
    mutate(qvalue = ifelse("qvalue" %in% names(df), qvalue, NA_real_))
  write.csv(df, outfile, row.names = FALSE)
  TRUE
}

# -------- GO enrichment (BP/MF/CC) --------
ego_bp <- suppressWarnings(
  enrichGO(gene = entrez,
           OrgDb = org.Hs.eg.db,
           keyType = "ENTREZID",
           ont = "BP",
           pAdjustMethod = opt$pAdj,
           qvalueCutoff = opt$qcut,
           readable = TRUE)
)

ego_mf <- suppressWarnings(
  enrichGO(gene = entrez,
           OrgDb = org.Hs.eg.db,
           keyType = "ENTREZID",
           ont = "MF",
           pAdjustMethod = opt$pAdj,
           qvalueCutoff = opt$qcut,
           readable = TRUE)
)

ego_cc <- suppressWarnings(
  enrichGO(gene = entrez,
           OrgDb = org.Hs.eg.db,
           keyType = "ENTREZID",
           ont = "CC",
           pAdjustMethod = opt$pAdj,
           qvalueCutoff = opt$qcut,
           readable = TRUE)
)

# -------- KEGG enrichment --------
# clusterProfiler::enrichKEGG (org = 'hsa' for human)
ekegg <- suppressWarnings(
  enrichKEGG(gene = entrez,
             organism = "hsa",
             keyType = "ncbi-geneid",
             pAdjustMethod = opt$pAdj,
             qvalueCutoff = opt$qcut)
)

# -------- Reactome enrichment --------
ereact <- suppressWarnings(
  enrichPathway(gene = entrez,
                organism = "human",
                pAdjustMethod = opt$pAdj,
                qvalueCutoff = opt$qcut,
                readable = TRUE)
)

# -------- Create output dir --------
if (!opt$dry_run) dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

# -------- Save CSVs --------
if (!opt$dry_run) {
  saved_bp <- save_enrich(ego_bp, file.path(opt$outdir, "GO_BP_enrichment.csv"))
  saved_mf <- save_enrich(ego_mf, file.path(opt$outdir, "GO_MF_enrichment.csv"))
  saved_cc <- save_enrich(ego_cc, file.path(opt$outdir, "GO_CC_enrichment.csv"))
  saved_kg <- save_enrich(ekegg,  file.path(opt$outdir, "KEGG_enrichment.csv"))
  saved_rc <- save_enrich(ereact, file.path(opt$outdir, "Reactome_enrichment.csv"))
  message(sprintf("Saved tables -> %s", opt$outdir))
} else {
  message("Dry-run: results were NOT written.")
}

# -------- Plots (dotplot) --------
plot_and_save <- function(ego, outfile_png, title_txt) {
  if (is.null(ego) || nrow(as.data.frame(ego)) == 0) return(FALSE)
  suppressPackageStartupMessages(library(enrichplot))
  p <- dotplot(ego, showCategory = 20) + ggplot2::ggtitle(title_txt)
  ggplot2::ggsave(outfile_png, p, width = 8, height = 5, dpi = 200)
  TRUE
}

if (!opt$dry_run) {
  ok_bp <- plot_and_save(ego_bp, file.path(opt$outdir, "GO_BP_dotplot.png"), "GO Biological Process")
  ok_mf <- plot_and_save(ego_mf, file.path(opt$outdir, "GO_MF_dotplot.png"), "GO Molecular Function")
  ok_cc <- plot_and_save(ego_cc, file.path(opt$outdir, "GO_CC_dotplot.png"), "GO Cellular Component")
  ok_kg <- plot_and_save(ekegg,  file.path(opt$outdir, "KEGG_dotplot.png"),  "KEGG Pathways")
  ok_rc <- plot_and_save(ereact, file.path(opt$outdir, "Reactome_dotplot.png"), "Reactome Pathways")
  message(sprintf("Saved dotplots -> %s", opt$outdir))
}

# -------- Console summary --------
summ <- function(ego, name) {
  df <- as.data.frame(ego)
  k <- ifelse(is.null(df) || nrow(df) == 0, 0, nrow(df))
  message(sprintf("%s: %d enriched terms (q < %.3f, method=%s)", name, k, opt$qcut, opt$pAdj))
}

summ(ego_bp, "GO BP")
summ(ego_mf, "GO MF")
summ(ego_cc, "GO CC")
summ(ekegg,  "KEGG")
summ(ereact, "Reactome")

message("Done.")
