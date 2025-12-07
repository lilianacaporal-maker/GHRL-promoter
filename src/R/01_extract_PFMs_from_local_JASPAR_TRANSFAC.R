
# Title: Load local PFMs from JASPAR (2024) and TRANSFAC, save as R objects
# Description: Reads manually downloaded PFMs for human TFs from local files and saves them as a list of PFMatrix objects in RDS.
# Note: For TRANSFAC®, use PFMs from your institutional license and DO NOT upload them to the public repository.
# Version: 1.1

# --- Dependencies ---
# Install missing packages if needed
if (!requireNamespace("TFBSTools", quietly = TRUE)) install.packages("TFBSTools")
if (!requireNamespace("Biostrings", quietly = TRUE)) BiocManager::install("Biostrings")

library(TFBSTools)
library(Biostrings)

# --- Parameters (edit these paths to match your local setup) ---
jaspar_dir <- file.path("data", "external", "JASPAR_PFMs_local")
transfac_dir <- file.path("data", "external", "TRANSFAC_PFMs_local")

output_dir <- file.path("data", "external")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

output_rds <- file.path(output_dir, "PFMs_JASPAR_TRANSFAC_Hs.rds")

# --- Helper functions ---
read_jaspar_pfms_from_dir <- function(dir_path) {
  if (!dir.exists(dir_path)) {
    warning("JASPAR directory does not exist: ", dir_path)
    return(list())
  }

  files <- list.files(dir_path, full.names = TRUE)
  if (length(files) == 0) {
    warning("No files found in JASPAR directory: ", dir_path)
    return(list())
  }

  pfm_list <- list()
  for (f in files) {
    ext <- tolower(tools::file_ext(f))
    if (ext %in% c("pfm", "jaspar", "txt")) {
      m <- tryCatch(
        readJASPARMatrix(f, matrixClass = "PFMatrix"),
        error = function(e) {
          tryCatch(
            importMatrix(file = f, format = "JASPAR", matrixClass = "PFMatrix"),
            error = function(e2) {
              message("Failed to read as JASPAR text: ", basename(f), " | Error: ", e2$message)
              NULL
            }
          )
        }
      )
      if (!is.null(m)) {
        if (inherits(m, "PFMatrix")) {
          pfm_list[[TFBSTools::name(m) %||% basename(f)]] <- m
        } else if (is.list(m) && all(vapply(m, inherits, logical(1), "PFMatrix"))) {
          for (mm in m) {
            pfm_list[[TFBSTools::name(mm) %||% paste0(basename(f), "_", mm@ID)]] <- mm
          }
        }
      }
    } else if (ext %in% c("meme")) {
      m <- tryCatch(
        importMatrix(file = f, format = "MEME", matrixClass = "PFMatrix"),
        error = function(e) {
          message("Failed to read MEME file: ", basename(f), " | Error: ", e$message)
          NULL
        }
      )
      if (!is.null(m)) {
        if (inherits(m, "PFMatrix")) {
          pfm_list[[TFBSTools::name(m) %||% basename(f)]] <- m
        } else if (is.list(m) && all(vapply(m, inherits, logical(1), "PFMatrix"))) {
          for (mm in m) {
            pfm_list[[TFBSTools::name(mm) %||% paste0(basename(f), "_", mm@ID)]] <- mm
          }
        }
      }
    } else {
      message("Skipping unsupported file: ", basename(f))
    }
  }
  pfm_list
}

read_transfac_pfms_from_dir <- function(dir_path) {
  if (!dir.exists(dir_path)) {
    message("TRANSFAC directory does not exist: ", dir_path)
    return(list())
  }
  files <- list.files(dir_path, full.names = TRUE)
  if (length(files) == 0) {
    message("No files found in TRANSFAC directory: ", dir_path)
    return(list())
  }

  pfm_list <- list()
  for (f in files) {
    m <- tryCatch(
      importMatrix(file = f, format = "TRANSFAC", matrixClass = "PFMatrix"),
      error = function(e) {
        message("Failed to read TRANSFAC file: ", basename(f), " | Error: ", e$message)
        NULL
      }
    )
    if (!is.null(m)) {
      if (inherits(m, "PFMatrix")) {
        pfm_list[[TFBSTools::name(m) %||% basename(f)]] <- m
      } else if (is.list(m) && all(vapply(m, inherits, logical(1), "PFMatrix"))) {
        for (mm in m) {
          pfm_list[[TFBSTools::name(mm) %||% paste0(basename(f), "_", mm@ID)]] <- mm
        }
      }
    }
  }
  pfm_list
}

`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0 && !is.na(a) && a != "") a else b

jaspar_pfms <- read_jaspar_pfms_from_dir(jaspar_dir)
transfac_pfms <- read_transfac_pfms_from_dir(transfac_dir)

pfm_list <- c(jaspar_pfms, transfac_pfms)

if (length(pfm_list) == 0) {
  stop("No PFMs loaded. Please check your directories and file formats.")
}

saveRDS(pfm_list, output_rds)
message("PFMs saved to: ", output_rds)
message("Total PFMs loaded: ", length(pfm_list))
preview_names <- utils::head(names(pfm_list), 10)
message("Example TF names/IDs: ", paste(preview_names, collapse = ", "))
