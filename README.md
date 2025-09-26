Agilent One-Color Microarray Analysis (E-MTAB-6646)

This repository contains an R pipeline for analyzing Agilent one-color microarray data, demonstrated on the E-MTAB-6646 dataset.
The workflow includes data import, background correction, normalization, probe filtering, differential expression analysis, and visualization.

Features:

1. Import and preprocess Agilent raw .txt files
2. Background correction (normexp) and quantile normalization
3. Removal of control and low-expression probes
4. Collapsing duplicate probes
5. Differential expression analysis using limma
6. Volcano plots with top gene labels
7. Export results to .csv, .pdf, and .png

Requirements:

1. R (≥ 4.0)
2. CRAN packages: ggplot2, dplyr
3. Bioconductor packages: limma, ggrepel
Packages are automatically checked/installed by the script.

Usage:

1. Clone the repository:
git clone https://github.com/shahabshahgaldi/microarray-analysis-E-MTAB-6646-
2. Open agilent_microarray_analysis.R in R or RStudio.
3. Update the setwd() line to point to your local directory containing raw Agilent .txt files.
4. Run the script.
5 Results will be written as CSV and volcano plots as PDF/PNG.

Output:

1. DE.results.preadipocytes.IAV.vs.mock.csv
2. DE.results.adipocytes.IAV.vs.mock.csv
3. Volcano plots (.pdf and .png)
4. QC plots (QC_plots.pdf)

License:

MIT License (feel free to adapt if you prefer GPL or another).
