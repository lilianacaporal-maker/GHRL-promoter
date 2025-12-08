
#!/usr/bin/env Rscript
# Functional enrichment with explicit gene universe (GO BP/MF/CC, KEGG, Reactome)
# - BH correction; q (FDR) < 0.05
# - Inputs: TF list (genes of interest) + universe list (background)
# - Supports TXT (one ID per line) or CSV (specify column names)
# - Outputs: CSV tables + PNG dotplots (disabled with --dry-run)

suppressPackageStartupMessages({
  library(optparse)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ReactomePA)
  library(dplyr)
})

option_list <- list(
  make_option(c("--genes"), type="character", help="Path to TF gene list (txt/csv)."),
  make_option(c("--universe"), type="character", help="Path to universe gene list (txt/csv)."),
  make_option(c("--genes-col"), type="character", default="gene_symbol", help="Column name for TFs if CSV (default: gene_symbol)."),
  make_option(c("--univ-col"),  type="character", default="gene_symbol", help="Column name for universe if CSV (default: gene_symbol)."),
  make_option(c("--keytype"),   type="character", default="SYMBOL",
              help="Key type of input IDs (e.g., SYMBOL, ENTREZID, ENSEMBL). Default: SYMBOL."),
  make_option(c("--qcut"),      type="double",    default=0.05, help="FDR (q-value) cutoff. Default: 0.05."),
  make_option(c("--pAdj"),      type="character", default="BH", help="Multiple testing correction. Default: BH."),
  make_option(c("--minGS"),     type="integer",   default=10,   help="Minimum gene set size."),
  make_option(c("--maxGS"),     type="integer",   default=500,  help="Maximum gene set size."),
  make_option(c("--outdir"),    type="character", default="results/enrichment", help="Output directory."),
  make_option(c("--dry-run"),   action="store_true", default=FALSE, help="Run without writing files.")
)
opt <- parse_args(OptionParser(option_list=option_list))

if (is.null(opt$genes) || is.null(opt$universe)) {
  stop("Please provide --genes and --universe paths.")
}

read_ids <- function(path, colname=NULL) {
  ext <- tolower(tools::file_ext(path))
  if (ext == "csv") {
    df <- read.csv(path, stringsAsFactors = FALSE)
    if (is.null(colname) || !(colname %in% colnames(df))) {
      stop(sprintf("Column '%s' not found in CSV %s", colname, path))
    }
    ids <- df[[colname]]
  } else {
    ids <- readLines(path, warn = FALSE)
  }
  ids <- toupper(trimws(ids))
  ids <- ids[ids != "" & !is.na(ids)]
  unique(ids)
}

gene_ids <- read_ids(opt$genes, opt$genes_col)
univ_ids <- read_ids(opt$universe, opt$univ_col)

if (length(gene_ids) < 3) stop("Too few genes in --genes for enrichment.")
if (length(univ_ids) < 50) warning("Universe seems small; enrichment may be unstable.")

# --- Map to Entrez using the declared keytype ---
map_to_entrez <- function(ids, fromType) {
  suppressWarnings(
    bitr(ids, fromType=fromType, toType="ENTREZID", OrgDb=org.Hs.eg.db) %>%
      distinct(ENTREZID) %>%
      pull(ENTREZID)
  )
}

gene_entrez <- map_to_entrez(gene_ids, opt$keytype)
univ_entrez <- map_to_entrez(univ_ids,  opt$keytype)

if (length(gene_entrez) < 3) stop("Mapping to Entrez produced <3 IDs for genes.")
if (length(univ_entrez) < 50) warning("Universe mapped to <50 Entrez IDs.")

# ReactomePA expects character Entrez IDs (not numeric) for universe/gene
gene_entrez_chr <- as.character(gene_entrez)
univ_entrez_chr <- as.character(univ_entrez)

message(sprintf("Genes mapped to Entrez: %d | Universe: %d", length(gene_entrez), length(univ_entrez)))

# --- Helpers to save results and plots ---
save_enrich <- function(eres, outfile) {
  if (is.null(eres) || nrow(as.data.frame(eres)) == 0) return(FALSE)
  df <- as.data.frame(eres) %>% arrange(p.adjust) %>%
    mutate(qvalue = ifelse("qvalue" %in% names(.), qvalue, NA_real_))
  write.csv(df, outfile, row.names = FALSE)
  TRUE
}

plot_and_save <- function(eres, outfile_png, title_txt) {
  if (is.null(eres) || nrow(as.data.frame(eres)) == 0) return(FALSE)
  suppressPackageStartupMessages(library(enrichplot))
  p <- dotplot(eres, showCategory = 20) + ggplot2::ggtitle(title_txt)
  ggplot2::ggsave(outfile_png, p, width = 8, height = 5, dpi = 200)
  TRUE
}

if (!opt$dry_run) dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

# --- GO enrichment (BP/MF/CC) with explicit universe ---
ego_bp <- suppressWarnings(
  enrichGO(gene = gene_entrez,
           OrgDb = org.Hs.eg.db,
           keyType = "ENTREZID",
           ont = "BP",
           universe = univ_entrez,         # explicit background
           pAdjustMethod = opt$pAdj,
           qvalueCutoff = opt$qcut,
           pvalueCutoff = 1.0,             # let qvalueCutoff drive significance
           minGSSize = opt$minGS,
           maxGSSize = opt$maxGS,
           readable = TRUE)
)
ego_mf <- suppressWarnings(
  enrichGO(gene = gene_entrez, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
           ont = "MF", universe = univ_entrez,
           pAdjustMethod = opt$pAdj, qvalueCutoff = opt$qcut,
           pvalueCutoff = 1.0, minGSSize = opt$minGS, maxGSSize = opt$maxGS, readable = TRUE)
)
ego_cc <- suppressWarnings(
  enrichGO(gene = gene_entrez, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
           ont = "CC", universe = univ_entrez,
           pAdjustMethod = opt$pAdj, qvalueCutoff = opt$qcut,
           pvalueCutoff = 1.0, minGSSize = opt$minGS, maxGSSize = opt$maxGS, readable = TRUE)
)

# --- KEGG enrichment with explicit universe ---
ekegg <- suppressWarnings(
  enrichKEGG(gene = gene_entrez_chr,     # as character
             organism = "hsa",
             keyType = "ncbi-geneid",
             universe = univ_entrez_chr, # explicit background
             pAdjustMethod = opt$pAdj,
             qvalueCutoff = opt$qcut,
             pvalueCutoff = 1.0,
             minGSSize = opt$minGS,
             maxGSSize = opt$maxGS)
)

# --- Reactome enrichment with explicit universe (character Entrez IDs) ---
ereact <- suppressWarnings(
  enrichPathway(gene = gene_entrez_chr,
                organism = "human",
                universe = univ_entrez_chr,  # explicit background
                pAdjustMethod = opt$pAdj,
                qvalueCutoff = opt$qcut,
                pvalueCutoff = 1.0,
                minGSSize = opt$minGS,
                maxGSSize = opt$maxGS,
                readable = TRUE)
)

# --- Save tables & dotplots unless dry-run ---
if (!opt$dry_run) {
  saved_bp <- save_enrich(ego_bp, file.path(opt$outdir, "GO_BP_enrichment.csv"))
  saved_mf <- save_enrich(ego_mf, file.path(opt$outdir, "GO_MF_enrichment.csv"))
  saved_cc <- save_enrich(ego_cc, file.path(opt$outdir, "GO_CC_enrichment.csv"))
  saved_kg <- save_enrich(ekegg,  file.path(opt$outdir, "KEGG_enrichment.csv"))
  saved_rc <- save_enrich(ereact, file.path(opt$outdir, "Reactome_enrichment.csv"))

  ok_bp <- plot_and_save(ego_bp, file.path(opt$outdir, "GO_BP_dotplot.png"), "GO Biological Process")
  ok_mf <- plot_and_save(ego_mf, file.path(opt$outdir, "GO_MF_dotplot.png"), "GO Molecular Function")
  ok_cc <- plot_and_save(ego_cc, file.path(opt$outdir, "GO_CC_dotplot.png"), "GO Cellular Component")
  ok_kg <- plot_and_save(ekegg,  file.path(opt$outdir, "KEGG_dotplot.png"),  "KEGG Pathways")
  ok_rc <- plot_and_save(ereact, file.path(opt$outdir, "Reactome_dotplot.png"), "Reactome Pathways")

  message(sprintf("Saved tables & dotplots -> %s", opt$outdir))
} else {
  message("Dry-run: results were NOT written.")
}

# --- Console summary ---
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
