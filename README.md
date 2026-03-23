# APPLE.R
APPLE is a comprehensive R package designed to analyze poly(A) tail lengths and alternative polyadenylation (APA) from sequencing data.
- [1. Introduction](#1-Introduction)
- [2. Installation of APPLE package](#2-Installation-of-APPLE-package)
- [3. Workflow of APPLE](#3.-Workflow-of-APPLE)
  - [3.1 Minimap2 alignment and filtering on the reads.(optional，Primer chopped)](#31-Minimap2-alignment-and-filtering-（optional）)
  - [3.2 Extracting polyAsite and detect polyA tail](#31-Extracting-polyAsite-and-detect-polyA-tail)
  - [3.3 Identifying clusters of alternative polyadenylation events and tails](#32-Identifying-clusters-of-alternative-polyadenylation-events-and-tails)
    - [3.3.1 Load raw polyA data](#331-Load-raw-polyA-data)
    - [3.3.2 Weighted density peak clustering](#332-Weighted-density-peak-clustering)
    - [3.3.3 Feature annotation and APA quantification](#333-Feature-annotation-and-APA-quantification)
    - [3.3.4 Filter low-confidence PolyA Clusters](#334-Filter-low-confidence-PolyA-Clusters)
    - [3.3.5 Map tail lengths to PolyA Clusters](#335-Map-tail-lengths-to-PolyA-Clusters)
  - [3.4 Statistical analysis of polyA tail length changes](#34-Statistical-analysis-of-polyA-tail-length-changes)
    - [3.4.1 Sample-level tests: t-test, Wilcoxon, and linear mixed models](#341-Sample-level-tests-t-test,-Wilcoxon-and-linear-mixed-models)
    - [3.4.2 Principal Component Analysis (PCA) on tail length matrices](#342-Principal-Component-Analysis-(PCA)-on-tail-length-matrices)
  - [3.5 Differential Expression Analysis of Poly(A) Sites](#35-Differential-Expression-Analysis-of-Poly(A)-Sites)
  - [3.6 Dynamic Analysis of APA at Gene Level](#36-Dynamic-Analysis-of-APA-at-Gene-Level)
- [4. Application of APPLE](#4-application-of-scdapa2-in-an-arabidopsis-dataset)


# 1. Introduction
Poly(A) tail length and alternative polyadenylation are key regulatory factors governing mRNA stability, translation, and localization. With advancements in 3' end sequencing technologies—such as PAT-seq, PAL-seq, and FLAM-seq—it is now possible to capture site-level information regarding poly(A) tail lengths and polyadenylation site usage. APPLE provides an integrated pipeline for processing such data; it identifies polyadenylation site clusters (PACs) and calculates their usage rates using a density peak clustering algorithm, and performs tail-length analyses on individual reads using statistical tests including the t-test, Wilcoxon test test, and linear mixed models.

# 2. Installation of APPLE package
APPLE relies on certain tools which are exclusive to the Linux environment. Therefore, it is advisable to install and utilize scDAPA2 within a Linux setting.
### [1]. Install dependencies
**Install standalone tools:** \
samtools (>=1.17), bedtools, minimap2(optional)
**Install R packages:** \
bedr, stringr, dplyr, tidyr,matrixStats, pbmcapply, FactoMineR, factoextra, ggplot2, uwot, lme4, lmerTest, GenomicRanges, GenomicFeatures, rtracklayer,Rsamtools, DESeq2, ggbio, readr, BiocGenerics, GenomeInfoDb,IRanges,  S4Vectors, SummarizedExperiment, methods, outliers, tidyselect, data.table, parallel
**If R >= 3.5.0:**
```         
if (!require("tools")) install.packages("tools")
if (!require("bedr")) install.packages("bedr")
if (!require("stringr")) install.packages("stringr")
if (!require("outliers")) install.packages("outliers")
if (!require("dplyr")) install.packages("dplyr")
if (!require("tidyr")) install.packages("tidyr")
if (!require("matrixStats")) install.packages("matrixStats")
if (!require("pbmcapply")) install.packages("pbmcapply")
if (!require("FactoMineR")) install.packages("FactoMineR")
if (!require("factoextra")) install.packages("factoextra")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("uwot")) install.packages("uwot")
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!require("plyranges")) BiocManager::install("plyranges")
if (!require("GenomicRanges")) BiocManager::install("GenomicRanges")
if (!require("GenomicFeatures")) BiocManager::install("GenomicFeatures")
if (!require("rtracklayer")) BiocManager::install("rtracklayer")
if (!require("Rsamtools")) BiocManager::install("Rsamtools")
if (!require("DESeq2")) BiocManager::install("DESeq2")
if (!require("ggbio")) BiocManager::install("ggbio")
if (!require("readr")) BiocManager::install("readr")
if (!require("stringr")) BiocManager::install("stringr")
```
### [2]. Install scDAPA2
```         
install.packages('devtools')
devtools::install_github("wangziiiiiii/APPLE")
```
# 3. Workflow of scDAPA2
The essential functions of APPLE include: (1) Minimap2 alignment and filtering on the reads.(optional) (2) Extracting polyAsite and detect polyA tail (3) Identifying clusters of alternative polyadenylation events and tails (4)Statistical analysis of polyA tail length changes (5)Differential Expression Analysis of Poly(A) Sites
## 3.1  Minimap2 alignment and filtering on the reads.(optional)
APPLE provides a wrapper function minimap2() that runs minimap2 on FASTQ files, filters reads by SAM flags, and sorts the output.
Note: It is recommended that the FASTQ files have already undergone primer chopped.
#### Usage
Process FASTQ files in a specified directory:
```
sam_files <- minimap2(reference = "path/to/genome.fa",
                      work_dir = "path/to/fastq/",
                      threads = 8,
                      filter_flags = 2308)
```
#### Arguments
```
 reference                     Reference genome file path
 work_dir                      Working directory path containing fq files
 threads                       Number of threads for minimap2, default is 4
 filter_flags                  SAM flags to filter out, default is 2308
```

## 3.2  Extracting polyAsite and detect polyA tail
After alignment, Extract_polyAsite() scans the SAM files, identifies reads containing poly(A) tails (via pt:i: tags), and writes BED‑like files with coordinates, strand, and tail lengths.
#### Usage
# Extract poly(A) tails from SAM files
```
results <- Extract_polyAsite(work_dir = "path/to/sam/",
                             intron_max = 50000,
                             min_tail_length = 6,
                             bedtools_path = "bedtools",
                             remove_temp_files = FALSE)
```
#### Arguments
```
 work_dir                      Working directory containing SAM files
 intron_max                    Maximum intron size, default is 50000
 min_tail_length               Minimum polyA tail length to consider, default is 6
 sample_size                   Number of lines to sample for pt tag detection, default is 1000
 bedtools_path                 Path to bedtools executable, default is "bedtools"
 remove_temp_files             Whether to remove temporary files, default is FALSE
```

## 3.3  Identifying clusters of alternative polyadenylation events and tails
### 3.3.1 Load raw polyA data
Load the .bed files produced by Extract_polyAsite() into a QuantifyPolyA object, which holds all raw poly(A) site data.
#### Usage
```
QpolyA <- Load.PolyA(files,dir = "path/to/bed/")
```
#### Arguments
```
files:                    A character vector specifying the names of BED files.
dir:                      A string specifying the directory containing BED files.
```
#### Output
```
A QuantifyPolyA object containing raw poly(A) site information in the @pre.polyA slot and tail length information in @tail_lengths.
```

### 3.3.2 Weighted density peak clustering
Cluster.PolyA() applies a weighted density peak clustering algorithm to group adjacent poly(A) sites into Poly(A) Clusters (PACs). The parameter max.gapwidth controls the maximum allowed gap between sites within a cluster. Clusters wider than max.gapwidth are further refined by a second clustering step.
#### Usage
```
QpolyA <- Cluster.PolyA(QpolyA, max.gapwidth = 24, mc.cores = 4)
```
#### Arguments
```
QpolyA:                   A QuantifyPolyA object containing clean poly(A) sites.
max.gapwidth:             Maximum distance between two adjacent sites in a PAC, default 24.
mc.cores:                 Number of cores for parallel clustering, default 4.
```
#### Output
```
An updated QuantifyPolyA object with PAC information stored in @polyA. Clusters that were split are recorded in @split.clusters.
```

### 3.3.3 Feature annotation and APA quantification
Annotate.PolyA() uses a genome annotation file (GTF/GFF) to assign each PAC to a gene and classify its location (e.g., 3’UTR, intron, intergenic). This is essential for downstream biological interpretation.
#### Usage
```
QpolyA <- Annotate.PolyA(QpolyA, gff = "path/to/annotation.gtf")
```
#### Arguments
```
QpolyA:                   A QuantifyPolyA object with PACs.
gff:                      A genome annotation file in GFF or GTF format (GTF recommended).
```
#### Output
```
An updated QuantifyPolyA object where the @polyA data frame includes additional columns: gene_id, distance, and type.
```


### 3.3.4 Filter low-confidence PolyA Clusters
Remove PACs with low read counts across samples using Filter.PolyA(). Only PACs with at least min_count reads in at least min_sample samples are retained.
#### Usage
```
QpolyA <- Filter.PolyA(QpolyA, min_count = 10, min_sample = 1)
```
#### Arguments
```
QpolyA:                   A QuantifyPolyA object with annotated PACs.
min_count:                Minimum read count in a PAC, default 10.
min_sample:               Minimum number of samples with `min_count` reads, default 1.
```
#### Output
```
A filtered QuantifyPolyA object. PACs not meeting criteria are removed from @polyA.
```

### 3.3.5 Map tail lengths to PolyA Clusters
After Filter, individual tail lengths must be assigned to the PACs they belong to. mapTail() performs this mapping, creating a sample‑wise table linking PAC IDs to concatenated tail lengths.
#### Usage
```
QpolyA <- mapTail(QpolyA, delimiter = ";")
```
#### Arguments
```
QpolyA:                   A QuantifyPolyA object with PACs defined.
delimiter:                Delimiter used to separate tail lengths in the output, default ";".
```
#### Output
```
A QuantifyPolyA object with tail length information mapped to clusters in the @cluster_tail_lengths slot. Each entry is a data frame containing columns: cluster_id, seqnames, start, end, strand, and all_tail_lengths.
```

## 3.4  Statistical analysis of polyA tail length changes
### 3.4.1 Sample-level tests: t-test, Wilcoxon, and linear mixed models
polyAlength() performs statistical comparisons for each PAC between a control group and multiple treatment groups. It uses all individual mRNA tail lengths as independent observations. The function can run t‑tests, Wilcoxon rank‑sum tests, and linear mixed models (LMM) with batch effects (via lib_id).
#### Usage
```
results <- polyAlength(QpolyA,
                       sample_info = sample_metadata,
                       test_methods = c("t_test", "wilcoxon", "lmm"),
                       min_mRNA_per_condition = 10,
                       logscale = TRUE,
                       control_group = "NC",
                       mc.cores = 4)
```
#### Arguments
```
QpolyA:                   A QuantifyPolyA object with cluster tail lengths.
sample_info:              A data.frame with columns `sample`, `condition`, and optionally `lib_id`.
test_methods:             Vector of tests to perform: "t_test", "wilcoxon", "lmm".
min_mRNA_per_condition:   Minimum number of mRNA molecules per condition for a PAC to be tested.
logscale:                 Whether to log2‑transform tail lengths before testing.
control_group:            Name of the control condition (must match values in `sample_info$condition`).
mc.cores:                 Number of cores for parallel processing.
```
#### Output
```
A data frame (wide format) with one row per PAC, containing:

(1) Summary statistics for control and each treatment (mean, median, sd, n).

(2) Test statistics, p‑values, and q‑values (FDR) for each test method.
```

### 3.4.2 Principal Component Analysis (PCA) on tail length matrices
Principal Component Analysis can be used to explore global patterns in poly(A) tail length variation across samples. APPLE provides a streamlined pipeline tail_pca() that aggregates tail lengths per PAC, handles missing values, performs PCA, and returns the coordinates for visualization.
#### Usage
```
pca_data <- tail_pca(QpolyA, sample_info,
                     aggregation_method = "mean",
                     max_missing = 0.2,
                     impute_method = "mean",
                     scale = TRUE,
                     center = TRUE)
```
#### Arguments
```
QpolyA:                   A QuantifyPolyA object that has been processed through clustering and tail length mapping.
sample_info:              A data.frame with sample metadata. Must contain columns `sample` and `condition` 
aggregation_method:       Method to aggregate tail lengths per PAC: either "mean" or "median". Default is "mean".
max_missing:              Maximum allowed proportion of missing values per PAC. PACs with more missing values are removed. Default is 0.2.
impute_method:            Method to handle remaining missing values: "mean" (impute with column mean), "knn" (k‑nearest neighbours, requires `impute` package), or "remove" (remove samples with any missing). Default is "mean".
scale:                    Logical; whether to scale variables to unit variance before PCA. Default is TRUE.
center:                   Logical; whether to center variables to zero mean before PCA. Default is TRUE.
```
#### Output
```
A data frame containing PCA coordinates for each sample (columns PC1, PC2, …), merged with the provided sample_info metadata. The returned object also has attributes variance_explained and cumulative_variance storing the proportion of variance explained by each principal component.
```


## 3.5 Differential Expression Analysis of Poly(A) Sites
APPLE integrates the DESeq2 tool to perform Variance Stabilizing Transformation (VST) and generate PCA and UMAP plots based on PAC counts. This function enables the identification of differential poly(A) site usage at the individual site level.
#### Usage
```
colData <- data.frame(
  condition = c("Control", "Control", "Treatment", "Treatment"),
  rownames = QpolyA@sample_names,
  type= as.factor('single')
)
results <- DESeq2.PolyA(QpolyA, colData)
```
#### Arguments
```
QpolyA:                  A QuantifyPolyA object containing PAC counts in the `@polyA` slot. PACs must have been filtered and annotated.
colData:                 A data.frame with sample metadata. Row names must match the sample names in `QpolyA@sample_names`, and must include a column named `condition` specifying the experimental groups.
```
#### Output
```
A list containing three elements:

(1) DESeq2.Result: The DESeq2 DESeqDataSet object after running DESeq().

(2) PCA.Plot: A PCA plot generated by factoextra::fviz_pca_ind(), colored by experimental condition.

(3) UMAP.Plot: A UMAP plot generated by uwot::umap(), colored by experimental condition.
```
## 3.6 Dynamic Analysis of APA at Gene Level
For focused analysis, APPLE provides several functions to quantify changes in the usage (relative abundance) of PACs within genes:

(1) Quantify.SplitAPA() – only for PACs that were split during clustering.

(2) Quantify.CanonicalAPA() – only for PACs annotated as 3’UTR or extended 3’UTR.

(3) Quantify.CNCAPA() – compares canonical (3’UTR) vs. non‑canonical PACs.

(4) Quantify.GeneAPA() – all PACs of a gene, regardless of annotation.

These functions compute a proportional difference (pd) and a correlation‑based statistic (r) for each gene, along with a p‑value (from chi‑squared test on usage counts).
#### Usage
```
apa_results <- Quantify.GeneAPA(QpolyA,
                                colData = sample_metadata,
                                contrast = c("condition", "Control", "Treatment"))
```
#### Arguments
```
QpolyA:                   A QuantifyPolyA object with annotated PACs.
colData:                  A data.frame with sample metadata, must include a `condition` column.
contrast:                 A three‑element vector: c(column, control_group, treatment_group).
```
#### Output
```
A tibble with columns: gene_id, pd, r, p.value.
```

# 4. Application of APPLE
