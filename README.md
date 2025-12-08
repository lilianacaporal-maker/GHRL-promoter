
# GHRL Promoter Analysis

## Project Description
This repository contains the bioinformatic pipeline for the in-silico analysis of the human **GHRL** gene promoter. The workflow integrates promoter sequence retrieval, transcription factor binding site (TFBS) prediction, evolutionary conservation analysis, epigenetic profiling, and functional annotation.

---

## Initial TFBS Prediction
Transcription factor binding sites (TFBS) were initially identified using **Alibaba 2.1** ([www.gene-regulation.com](http://www.gene-regulation.com)), which relies on TRANSFAC® matrices.  
**Note:** Due to licensing restrictions, TRANSFAC® matrices and results are not included in this repository.  
To reproduce the analysis with open-access tools, we recommend using **RSAT** ([https://rsat.france-bioinformatique.fr](https://rsat.france-bioinformatique.fr)) and **JASPAR 2024** matrices.

---

## Methods Overview
- **Promoter sequence retrieval**: EPD database and UCSC Genome Browser.
- **TFBS prediction**: Alibaba 2.1 (TRANSFAC®), RSAT, and validation with JASPAR 2024 PFMs.
- **Affinity calculation**: Performed in R using Bioconductor packages (`TFBSTools`, `Biostrings`, `GenomicRanges`, `Gviz`).
- **Evolutionary conservation**: phastCons and phyloP scores from UCSC BigWig files.
- **Epigenetic profiling**: DNase-seq, ATAC-seq, H3K4me3, H3K27ac from ENCODE.

---

## Installation

### R Environment
```r
install.packages(c("BiocManager", "dplyr", "readxl", "plotly"))
BiocManager::install(c("TFBSTools", "Biostrings", "GenomicRanges", "Gviz"))
```

### Python Environment (optional)
```bash
python -m venv venv
source venv/bin/activate   # Linux/Mac
# venv\Scripts\activate   # Windows
pip install numpy pandas matplotlib seaborn pybigwig biopython tqdm click requests
```

---

## How to Reproduce

### R (TFBS prediction and affinity analysis)
```bash
Rscript scripts_R/tfbs_pipeline.R
# Outputs:
# results/tfbs_predictions.csv
# results/tfbs_bins_100bp.csv
```

### Python (conservation and epigenetic analysis)
```bash
# Download promoter sequence
python scripts/download_epd_sequence.py --tss chr3:10290734 --strand - --up 1000 --down 100 --out data/GHRL_promoter.fa

# Conservation analysis
python scripts/conservation_analysis.py --tss chr3:10290734 --up 1000 --down 100 --phastcons_bw data/bw/hg38_phastCons100way.bw --phylop_bw data/bw/hg38_phyloP100way.bw --out_prefix results/conservation

# Epigenetic profiling
python scripts/epigenetic_pipeline.py --tss chr3:10290734 --up 1000 --down 100 --dnase_bw <URL> --atac_bw <URL> --h3k4me3_bw <URL> --h3k27ac_bw <URL> --ctcf_bw <URL> --out_prefix results/epigenetics
```

---

## Data Availability
- [EPD Database](https://epd.epfl.ch)
- [UCSC Genome Browser](https://genome.ucsc.edu)
- [ENCODE Project](https://www.encodeproject.org)
- [JASPAR 2024](https://jaspar.elixir.no)

---

## Sensitive Data Notice
This repository does **not** include clinical data, unpublished figures, or proprietary TRANSFAC® matrices. For reproducibility, use open-access resources listed above.

---

## License
- **Code**: MIT License (see `LICENSE`)
- **Figures, text, and data**: Creative Commons Attribution 4.0 International (CC-BY 4.0)

---

## Citation
Please cite this work as described in `CITATION.cff`.


## 📂 Repository Structure

# GHRL-promoter

This repository contains resources to study the promoter region of the human **GHRL** (ghrelin) gene, including:

- **Transcription factor binding site (TFBS) prediction** using JASPAR and TRANSFAC matrices.
- **Evolutionary conservation analysis** with phastCons and phyloP scores.
- **Epigenetic profiling** (H3K4me3, H3K27ac, DNase, ATAC-seq) using ENCODE data.
- **Functional enrichment analysis** for predicted transcription factors.
- **Visualization scripts** for motifs, conservation tracks, and epigenetic marks.

---


## Data & Results Policy
This repository **does not** include raw data, external references (FASTA/BigWig), or generated results (CSV/figures), because the associated manuscript is under preparation and not yet published.

- **Do not commit** any files under `data/raw/`, `data/external/`, or `results/`.
- A `.gitignore` and a `pre-commit` hook are configured to prevent accidental publishing of sensitive or unpublished outputs.
- To reproduce locally, download public tracks (EPD/Ensembl, UCSC, ENCODE) and run the scripts with `--dry-run` for validation, or without it to generate outputs in your local machine. Outputs are **not tracked** by Git.
---


## 📂 Repository Structure
ghrelin-promoter-analysis/
├── README.md                # Project documentation
├── LICENSE                  # MIT for code + CC-BY for data
├── CITATION.cff             # Citation metadata
├── .gitignore               # Ignore heavy/sensitive files
├── data/                    # Promoter FASTA and links to open data
│   └── GHRL_promoter.fa
├── scripts_R/               # R scripts for TFBS and motif analysis
│   ├── tfbs_pipeline.R
│   └── find_inr.R
├── scripts/                 # Optional Python scripts for conservation/epigenetics
│   ├── download_epd_sequence.py
│   ├── conservation_analysis.py
│   ├── epigenetic_pipeline.py
│   └── summarize_tfbs_bins.py
├── results/                 # CSV outputs (no unpublished figures)
│   ├── tfbs_predictions.csv
│   └── tfbs_bins_100bp.csv
