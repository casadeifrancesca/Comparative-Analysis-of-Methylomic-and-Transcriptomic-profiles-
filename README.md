# Comparative Analysis of Methylation and Gene Expression Profiles in PDAC Patients

## Project Overview
This repository contains a multi-omic bioinformatics pipeline designed to investigate the molecular mechanisms underlying phenotypic differences in patients with **Pancreatic Ductal Adenocarcinoma (PDAC)** compared to healthy controls. 

The project integrates transcriptomic and epigenomic data to identify potential biomarkers through:
* **Differential Analysis**: Identification of DEGs (Expression) and DMGs (Methylation).
* **Multi-omic Integration**: Correlation between DNA methylation (Beta values) and gene expression (FPKM).
* **Epigenomics**: Mapping of differentially methylated probes within Transcription Factor Binding Sites (TFBS).
* **Functional Enrichment**: Gene Ontology Analysis.
* **Survival Analysis**: Evaluation of the prognostic impact of key genes.

---

## Requirements

### R Environment
**R version:** >= 4.0  
**Bioconductor packages:**
```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("minfi", "IlluminaHumanMethylation450kanno.ilmn12.hg19", 
                       "clusterProfiler", "org.Hs.eg.db", "AnnotationDbi", 
                       "GenomicRanges", "rtracklayer"))
```

**CRAN packages:** `data.table`, `dplyr`, `tibble`, `here`, `ggplot2`, `gridExtra`, `tidyr`, `gprofiler2`, `tidyverse`

### Python Environment (For plotting)
**Python version:** >= 3.8  
**Required libraries:** `pandas`, `matplotlib`, `seaborn`, `numpy`

---

## 📂 Project Structure
The pipeline consists of the following scripts, which should be executed in numerical order:

* **`01_differential_analysis.R`**: Data pre-processing, IQR filtering, and identification of DEGs and DMGs. Generates Volcano Plots.
* **`02_correlation_analysis.R`**: Spearman correlation analysis between gene expression and DNA methylation levels.
* **`03_functional_enrichment.R`**: Automated Gene Ontology (BP, MF, CC) enrichment analysis using the g:Profiler API.
* **`04_functional_enrichment_bubble_plot.py`**: Python-based visualization to generate high-resolution Bubble Plots from enrichment results.
* **`05_tfbs_analysis.R`**: Mapping of methylation probes to Transcription Factor Binding Sites (TFBS) and probe-level differential analysis.
* **`data/`**: Directory for input datasets (see Data Setup).
* **`results/`**: Directory containing the specific outputs of my analysis for the TCGA-PAAD dataset (plots, tables, and reports).

> **Note on Results:** The `results/` folder currently reflects my specific findings. If you run this pipeline with a different dataset, the scripts will automatically update or overwrite these files.

---

## Survival Analysis (KM Plotter)
The prognostic potential of the identified key genes was evaluated using the **Kaplan-Meier Plotter** online platform ([kmplot.com](http://www.kmplot.com)). 

* **Dataset:** TCGA-PAAD (Pancreatic Adenocarcinoma).
* **Stratification Strategy:** Patients were stratified into "High" and "Low" expression groups.
* **Cutoff Method:** **"Auto cutoff"** algorithm. This method identifies the optimal threshold between the 25th and 75th percentiles that maximizes the statistical significance of the group separation.
* **Statistical Metrics:** Log-rank P-value and Hazard Ratio (HR) with 95% Confidence Interval (CI).
* **Outputs:** Resulting Kaplan-Meier plots are stored in `results/survival_analysis/`.

---

## Data Setup & Usage

### Data Source
Analyses were performed using data retrieved from the **openGDC** database:
* **Gene Expression**: FPKM-normalized counts.
* **DNA Methylation**: Beta values (Illumina 450k platform).

### Setup Instructions
1. **Prepare Data**: Download the PDAC expression and methylation files and place them in the `data/` folder.
2. **Metadata**: Ensure `unified_tfbs.bed` and `cpg2gene__corrected.tsv` are present in the `data/` directory for TFBS mapping.
3. **Execution**: 
   * Open the `Analisi.Rproj` file in RStudio to set the working directory automatically via `here`.
   * Run the scripts sequentially from `01` to `05`.

