# APPLE.R

**Version 1.0**

*Analysis of poly(A) tail lengths and alternative polyadenylation in R.*

APPLE connects read processing, poly(A) site clustering, genomic annotation, tail-length comparisons, and gene-level APA analysis in one workflow.

| Analysis | Main outputs |
| --- | --- |
| Poly(A) sites | Annotated poly(A) site clusters (PACs) and sample counts |
| Poly(A) tails | PAC-level tail-length summaries, statistical comparisons, and PCA |
| APA dynamics | Changes in PAC usage, relative poly(A) position, and directional shifts |

**Contents**

- [1. Introduction](#1-introduction)
- [2. Installation](#2-installation)
- [3. Workflow](#3-workflow)
  - [3.1 Read alignment](#31-read-alignment)
  - [3.2 Site and tail extraction](#32-site-and-tail-extraction)
  - [3.3 PAC identification and annotation](#33-pac-identification-and-annotation)
  - [3.4 Tail-length analysis](#34-tail-length-analysis)
  - [3.5 Differential PAC counts](#35-differential-pac-counts)
  - [3.6 Gene-level APA analysis](#36-gene-level-apa-analysis)
- [4. Worked example](#4-worked-example)

## 1. Introduction

APPLE analyzes poly(A) site usage and tail-length measurements from sequencing data. It groups nearby sites into poly(A) site clusters (PACs), annotates their genomic context, and builds sample-level PAC count data. Downstream functions compare individual mRNA tail lengths, analyze PAC abundance, and summarize shifts in within-gene PAC usage.

These analyses address different questions: changes in tail length, changes in site abundance, and redistribution among sites within the same gene. Input reads and tail-length annotations must be prepared for the extraction workflow described below.

## 2. Installation

Use Linux for the complete workflow, which invokes command-line tools and uses multicore processing. Install external tools separately and make sure they are available on your PATH.

### Install dependencies

**Install standalone tools:** \
samtools (>=1.17), bedtools, and minimap2 (needed only for the alignment step)
**Install R packages:** \
bedr, stringr, dplyr, tidyr,matrixStats, pbmcapply, FactoMineR, factoextra, ggplot2, uwot, lme4, lmerTest, GenomicRanges, GenomicFeatures, rtracklayer,Rsamtools, DESeq2, ggbio, readr, BiocGenerics, GenomeInfoDb,IRanges,  S4Vectors, SummarizedExperiment, methods, outliers, tidyselect, data.table, parallel
**Install R dependencies:**

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

### Install APPLE

```
install.packages('devtools')
devtools::install_github("wangziiiiiii/APPLE.R")
```

The package declares R dependencies in DESCRIPTION; external executables must be installed separately. While this repository is private, GitHub authentication with access to APPLE.R is required for installation.

## 3. Workflow

Follow steps 3.1–3.3 to prepare annotated PACs, then use tail-length analysis, differential PAC count analysis, or gene-level APA analysis as appropriate for your question.

### 3.1 Read alignment

The minimap2() wrapper aligns FASTQ files to a reference genome, filters SAM records, and sorts the output. Skip this step if suitable SAM alignments are already available.
Remove primers before alignment. The wrapper accepts .fq, .fastq, and their gzip-compressed equivalents.

#### Usage

Process FASTQ files in a specified directory:

```
sam_files <- minimap2(reference = "/path/to/genome.fa",
                      work_dir = "/path/to/fastq/",
                      threads = 8,
                      filter_flags = 2308)
```

#### Arguments
```
reference                 Reference genome file path
work_dir                  Working directory path containing fq files
threads                   Number of threads for minimap2, default is 4
filter_flags              SAM flags to filter out, default is 2308
```

### 3.2 Site and tail extraction

Extract_polyAsite() scans SAM alignments and writes BED-like poly(A) site files with tail-length information. It samples records to detect pt:i: tags and uses those values when available; otherwise, it attempts sequence-based tail detection with regular expressions. Tag-based and sequence-based estimates depend on the upstream data preparation and should not be assumed interchangeable.

Use an absolute work_dir path: this function changes the R working directory. Set remove_temp_files = TRUE when loading all resulting .bed files from the directory, so intermediate files are excluded.

#### Usage

Extract poly(A) sites and tail lengths from SAM files:

```
results <- Extract_polyAsite(work_dir = "/path/to/sam/",
                             intron_max = 50000,
                             min_tail_length = 6,
                             bedtools_path = "bedtools",
                             remove_temp_files = FALSE)
```

#### Arguments
```
work_dir                  Working directory containing SAM files
intron_max                Maximum intron size, default is 50000
min_tail_length           Minimum polyA tail length to consider, default is 6
sample_size               Number of lines to sample for pt tag detection, default is 1000
bedtools_path             Path to bedtools executable, default is "bedtools"
remove_temp_files         Whether to remove temporary files, default is FALSE
```

### 3.3 PAC identification and annotation

Load site data, define PACs, assign genomic features, filter low-count clusters, and map tail lengths to the retained clusters.

#### 3.3.1 Load raw polyA data

Load the .bed files produced by Extract_polyAsite() into a QuantifyPolyA object, which holds all raw poly(A) site data.

#### Usage

```
QpolyA <- Load.PolyA(dir = "/path/to/bed/")
```

#### Arguments
```
files                     Optional vector of existing file paths; if supplied, dir is not prepended.
dir                       Directory to scan for .bed files when files is omitted.
```

#### Output
```
A QuantifyPolyA object containing raw poly(A) site information in the @pre.polyA slot and tail length information in @tail_lengths.
```

#### 3.3.2 Weighted density peak clustering

Cluster.PolyA() applies a weighted density peak clustering algorithm to group adjacent poly(A) sites into Poly(A) Clusters (PACs). The parameter max.gapwidth controls the maximum allowed gap between sites within a cluster. Clusters wider than max.gapwidth are further refined by a second clustering step.

#### Usage

```
QpolyA <- Cluster.PolyA(QpolyA, max.gapwidth = 24, mc.cores = 4)
```

#### Arguments
```
QpolyA                    A QuantifyPolyA object containing clean poly(A) sites.
max.gapwidth              Maximum distance between two adjacent sites in a PAC, default 24.
mc.cores                  Number of cores for parallel clustering, default 4.
```

#### Output
```
An updated QuantifyPolyA object with PAC information stored in @polyA. Clusters that were split are recorded in @split.clusters.
```

#### 3.3.3 Feature annotation and APA quantification

Annotate.PolyA() uses a genome annotation file (GTF/GFF) to assign each PAC to a gene and classify its location (e.g., 3’UTR, intron, intergenic). This is essential for downstream biological interpretation.

#### Usage

```
QpolyA <- Annotate.PolyA(QpolyA, gff = "/path/to/annotation.gtf")
```

#### Arguments
```
QpolyA                    A QuantifyPolyA object with PACs.
gff                       A genome annotation file in GFF or GTF format (GTF recommended).
```

#### Output
```
An updated QuantifyPolyA object where the @polyA data frame includes additional columns: gene_id, distance, and type.
```

#### 3.3.4 Filter low-confidence PolyA Clusters

Remove PACs with low read counts across samples using Filter.PolyA(). Only PACs with at least min_count reads in at least min_sample samples are retained.

#### Usage

```
QpolyA <- Filter.PolyA(QpolyA, min_count = 10, min_sample = 1)
```

#### Arguments
```
QpolyA                    A QuantifyPolyA object with annotated PACs.
min_count                 Minimum read count in a PAC, default 10.
min_sample                Minimum number of samples with `min_count` reads, default 1.
```

#### Output
```
A filtered QuantifyPolyA object. PACs not meeting criteria are removed from @polyA.
```

#### 3.3.5 Map tail lengths to PolyA Clusters

After filtering, map individual tail lengths to their PACs. mapTail() performs this mapping, creating a sample‑wise table linking PAC IDs to concatenated tail lengths.

#### Usage

```
QpolyA <- mapTail(QpolyA, delimiter = ";")
```

#### Arguments
```
QpolyA                    A QuantifyPolyA object with PACs defined.
delimiter                 Delimiter used to separate tail lengths in the output, default ";".
```

#### Output
```
A QuantifyPolyA object with tail length information mapped to clusters in the @cluster_tail_lengths slot. Each entry is a data frame containing columns: cluster_id, seqnames, start, end, strand, and all_tail_lengths.
```

### 3.4 Tail-length analysis

#### 3.4.1 Comparing tail lengths at the mRNA level

polyAlength() compares individual mRNA tail lengths within each PAC between a control condition and one or more treatment conditions. The t-test and Wilcoxon options pool mRNA observations within each condition; they do not test sample-level averages. The LMM option uses a random intercept for lib_id when multiple library IDs are available, and otherwise fits a linear model. Supply library IDs that reflect the experimental design.

With logscale = TRUE, the test helpers use log2-transformed tail lengths. Summary means and medians are calculated from the original tail lengths. The function applies Benjamini–Hochberg adjustment across tested PACs separately for each treatment and test method.

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
QpolyA                    A QuantifyPolyA object with cluster tail lengths.
sample_info               A data.frame with columns `sample`, `condition`, and optionally `lib_id`.
test_methods              Vector of tests to perform: "t_test", "wilcoxon", "lmm".
min_mRNA_per_condition    Minimum number of mRNA molecules per condition for a PAC to be tested.
logscale                  Whether to log2‑transform tail lengths before testing.
control_group             Name of the control condition (must match values in `sample_info$condition`).
mc.cores                  Number of cores for parallel processing.
```

#### Output
```
A data frame (wide format) with one row per PAC, containing:

(1) Summary statistics for control and each treatment (mean, median, sd, n).

(2) Test statistics, p‑values, and q‑values (FDR) for each test method.
```

#### 3.4.2 Principal Component Analysis (PCA) on tail length matrices

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
QpolyA                    A QuantifyPolyA object that has been processed through clustering and tail length mapping.
sample_info               A data.frame with sample metadata. Must contain columns `sample` and `condition`
aggregation_method        Method to aggregate tail lengths per PAC: either "mean" or "median". Default is "mean".
max_missing               Maximum allowed proportion of missing values per PAC. PACs with more missing values are removed. Default is 0.2.
impute_method             Method to handle remaining missing values: "mean" (impute with column mean), "knn" (k‑nearest neighbours, requires `impute` package), or "remove" (remove samples with any missing). Default is "mean".
scale                     Logical; whether to scale variables to unit variance before PCA. Default is TRUE.
center                    Logical; whether to center variables to zero mean before PCA. Default is TRUE.
```

#### Output
```
A data frame containing PCA coordinates for each sample (columns PC1, PC2, …), merged with the provided sample_info metadata. The returned object also has attributes variance_explained and cumulative_variance storing the percentage of variance explained by each principal component.
```

### 3.5 Differential PAC counts

DESeq2.PolyA() fits a DESeq2 model to PAC counts using design ~ condition. It also applies a variance-stabilizing transformation and generates PCA and UMAP plots. This analysis measures changes in PAC abundance; it does not directly model each PAC's proportion within its gene. Use section 3.6 to examine within-gene APA usage.

#### Usage

```
colData <- data.frame(
  condition = c("Control", "Control", "Treatment", "Treatment"),
  row.names = QpolyA@sample_names,
  type= as.factor('single')
)
results <- DESeq2.PolyA(QpolyA, colData)
```

#### Arguments
```
QpolyA                    A QuantifyPolyA object containing PAC counts in the `@polyA` slot. PACs must have been filtered and annotated.
colData                   A data.frame with sample metadata. Row names must match the sample names in `QpolyA@sample_names`, and must include a column named `condition` specifying the experimental groups.
```

#### Output
```
A list containing three elements:

(1) DESeq2.Result: The DESeq2 DESeqDataSet object after running DESeq().

(2) PCA.Plot: A PCA plot generated by factoextra::fviz_pca_ind(), colored by experimental condition.

(3) UMAP.Plot: A UMAP plot generated by uwot::umap(), colored by experimental condition.
```

### 3.6 Gene-level APA analysis

Quantify.GeneAPA() compares the distribution of PAC counts within each gene. It excludes intergenic PACs and retains genes with at least two remaining PACs. The contrast order is **column, control, treatment**.

```
apa_results <- Quantify.GeneAPA(
  QpolyA,
  colData = sample_metadata,
  contrast = c("condition", "Control", "Treatment")
)

gene_RPP <- compute_gene_RPP(
  polyA = QpolyA@polyA,
  sample_names = rownames(sample_metadata)
)

delta_RPP <- compute_delta_RPP(
  polyA_rank = gene_RPP,
  colData = sample_metadata,
  control_cond = "Control",
  treat_cond = "Treatment"
)

apa_results <- dplyr::left_join(apa_results, delta_RPP, by = "gene_id")

```

#### Arguments
```
QpolyA                    A QuantifyPolyA object containing annotated PACs.
colData                   Sample metadata; row names must match sample count columns.
                          Must contain a condition column.
contrast                  c(column, control_group, treatment_group).
polyA                     PAC data frame, typically QpolyA@polyA.
sample_names              Names of the sample count columns in polyA.
type_col                  PAC annotation column, default "type".
gene_id_col               Gene identifier column, default "gene_id".
strand_col                Strand column, default "strand".
center_col                PAC center coordinate column, default "center".
polyA_rank                Gene-level RPP table returned by compute_gene_RPP().
control_cond              Control condition name in colData$condition.
treat_cond                Treatment condition name in colData$condition.
```

#### Output
```
Quantify.GeneAPA()
  One row per gene with columns: gene_id, pd, r, p.value.

compute_gene_RPP()
  One row per gene with gene_id and one RPP column per sample.

compute_delta_RPP()
  One row per gene with columns: gene_id, delta_RPP.

After joining by gene_id
  Columns: gene_id, pd, r, p.value, delta_RPP.
```

compute_gene_RPP() excludes intergenic PACs and genes with fewer than two retained PACs.

#### Poly(A) site usage (PSU)

For gene $g$, PAC $k$, and sample $s$, define usage as:

$$
P_{g,s,k} = \mathrm{PSU}_{g,s,k}
= \frac{n_{g,s,k}}{\sum_{\ell=1}^{K_g} n_{g,s,\ell}}
$$

Here, $n$ is the PAC read count and $K_g$ is the number of retained PACs in the gene. With a positive gene total, usage proportions sum to 1 within each sample.

#### Proportion difference (PD)

PD measures the magnitude of the change in PAC usage. For each control–treatment replicate pair, the code sums absolute usage differences across PACs and divides by two. It then averages over all cross-condition replicate pairs:

$$
\mathrm{PD}_g =
\frac{1}{c_C c_T}
\sum_{i=1}^{c_C}\sum_{j=1}^{c_T}
\frac{1}{2}\sum_{k=1}^{K_g}
\left|P_{g,C_i,k}-P_{g,T_j,k}\right|
$$

$c_C$ and $c_T$ are the numbers of control and treatment replicates. For valid, nonzero sample totals, PD ranges from 0 to 1: 0 indicates identical usage distributions, and larger values indicate greater redistribution. **PD has no sign and does not describe shortening or lengthening.** The implementation uses the sum-based expression above, not the maximum difference at a single PAC.

#### RPP

Relative poly(A) position (RPP) summarizes proximal versus distal PAC usage. PACs are ordered along the direction of transcription: increasing center coordinates on the positive strand and decreasing coordinates on the negative strand. For distinct centers, the rank weight is:

$$
w_{g,k}=\frac{k-1}{K_g-1},
\qquad
\mathrm{RPP}_{g,s}=\sum_{k=1}^{K_g}w_{g,k}P_{g,s,k}
$$

The implementation uses percent_rank(), so tied centers share a rank. With positive sample totals, RPP ranges from 0 to 1. Larger values indicate more distal usage; smaller values indicate more proximal usage. RPP is a **rank-weighted position score**, not a physical length in nucleotides.

#### Delta RPP

The change in relative poly(A) position is reported as delta_RPP. compute_delta_RPP() subtracts the control-group mean RPP from the treatment-group mean:

$$
\Delta\mathrm{RPP}_g =
\frac{1}{c_T}\sum_{j=1}^{c_T}\mathrm{RPP}_{g,T_j}
-
\frac{1}{c_C}\sum_{i=1}^{c_C}\mathrm{RPP}_{g,C_i}
$$

| Result | Interpretation |
| --- | --- |
| delta_RPP > 0 | Shift toward distal PACs in treatment |
| delta_RPP < 0 | Shift toward proximal PACs in treatment |
| delta_RPP = 0 | No net change in the mean rank-weighted position; usage can still change |

The default analysis includes genic PACs outside the 3′UTR. Interpret a shift as **3′UTR lengthening or shortening only when the selected PAC set and annotation support that interpretation**. The difference retains its sign; it is not an absolute difference.

#### Correlation score and P values

The r column is a count-weighted correlation between condition identity (control = 1, treatment = 2) and PAC center coordinate, corrected for strand and averaged across replicate pairs. Positive values indicate a distal shift in treatment; negative values indicate a proximal shift. Unlike RPP, this score uses genomic coordinates rather than position ranks.

For each replicate pair, dynamicsDetect() runs a chi-squared test on a PAC-by-condition count table after removing rows with zero combined counts. If only one row remains, it assigns P = 1. The returned p.value is the **arithmetic mean of these pairwise P values**. It is not a P value from a joint generalized linear model, a formal combined P value, or an FDR-adjusted value.

The example screening rule **pd > 0.1 and p.value < 0.05** can be applied explicitly by the user. It is not applied automatically by the function and should not be treated as a universal significance criterion.

**Zero-count handling:** usage is undefined when a gene has zero total counts in a sample. The current PD routine skips such pairs but averages without removing missing values, which can yield NA. The RPP routine sums with na.rm = TRUE and can return 0 for zero-total samples; that value must not be interpreted as evidence of proximal usage. Check coverage in both conditions before interpreting results.

## 4. Worked example

This example uses the six FASTQ files in the repository's [test directory](test): two NC samples, two Fip1 samples, and two Fip2 samples. The supplied test reads have undergone primer removal. Small test subsets are intended to illustrate the workflow; filtering, statistical tests, and visualizations depend on the retained coverage.

#### Prepare files and paths

Download the test FASTQ files into one directory. Obtain the GRCh38 primary-assembly FASTA and the matching release 113 GTF from the [Ensembl human FASTA directory](https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/dna/) and [GTF directory](https://ftp.ensembl.org/pub/release-113/gtf/homo_sapiens/). Prepare an uncompressed FASTA and GTF, and use matching chromosome names.

Replace the paths below before running the example. Use absolute paths because extraction changes the working directory.

```
library(APPLE)

work_dir <- normalizePath("/path/to/test", mustWork = TRUE)
reference <- normalizePath(
  "/path/to/Homo_sapiens.GRCh38.dna.primary_assembly.fa",
  mustWork = TRUE
)
annotation <- normalizePath(
  "/path/to/Homo_sapiens.GRCh38.113.gtf",
  mustWork = TRUE
)

```

#### Align reads and extract sites

```
sam_files <- minimap2(
  reference = reference,
  work_dir = work_dir,
  threads = 4,
  filter_flags = 2308
)

extraction_stats <- Extract_polyAsite(
  work_dir = work_dir,
  intron_max = 50000,
  min_tail_length = 6,
  sample_size = 10000,
  bedtools_path = "bedtools",
  remove_temp_files = TRUE
)

```

#### Load, cluster, annotate, and filter PACs

Select the final BED files corresponding to the alignments, rather than unrelated BED files in the directory.

```
bed_files <- sub("\\.sam$", ".bed", unname(sam_files))
stopifnot(all(file.exists(bed_files)))

QpolyA <- Load.PolyA(files = bed_files)
QpolyA <- Cluster.PolyA(QpolyA, max.gapwidth = 24, mc.cores = 4)
QpolyA <- Annotate.PolyA(QpolyA, gff = annotation)
QpolyA <- Filter.PolyA(QpolyA, min_count = 10, min_sample = 2)
QpolyA <- mapTail(QpolyA)

```

#### Define sample metadata

Load.PolyA() derives sample names from BED file basenames. Preserve those exact names in the metadata. The condition labels below are derived from the NC, Fip1, and Fip2 filename prefixes in this test dataset.

```
sample_names <- QpolyA@sample_names
conditions <- sub("-.*$", "", sample_names)
stopifnot(
  length(sample_names) == 6L,
  all(conditions %in% c("NC", "Fip1", "Fip2")),
  all(table(factor(conditions, levels = c("NC", "Fip1", "Fip2"))) == 2L)
)

sample_info <- data.frame(
  sample = sample_names,
  condition = conditions,
  lib_id = sample_names,
  stringsAsFactors = FALSE
)

colData <- data.frame(
  condition = factor(conditions, levels = c("NC", "Fip1", "Fip2")),
  row.names = sample_names
)

```

Here, each sample is treated as a separate library. For other datasets, assign lib_id according to the actual experimental design.

#### Compare tail lengths

```
tail_results <- polyAlength(
  QpolyA = QpolyA,
  sample_info = sample_info,
  test_methods = c("t_test", "wilcoxon", "lmm"),
  min_mRNA_per_condition = 10,
  control_group = "NC",
  logscale = TRUE,
  mc.cores = 4
)

head(tail_results)

```

#### Explore tail-length variation

```
tail_coordinates <- tail_pca(
  QpolyA = QpolyA,
  sample_info = sample_info,
  aggregation_method = "mean",
  max_missing = 0.2,
  impute_method = "mean",
  scale = TRUE,
  center = TRUE
)

head(tail_coordinates)

```

PCA requires at least three samples and two retained PAC variables in the current implementation. Scaling also requires nonzero variance in the retained variables.

#### Analyze PAC counts

```
pac_results <- DESeq2.PolyA(QpolyA, colData)

```

The returned list contains DESeq2.Result, PCA.Plot, and UMAP.Plot. The wrapper runs both visualizations; UMAP and VST can require more samples or adequate counts than a small demonstration subset provides.

#### Quantify APA changes and join directional scores

```
gene_RPP <- compute_gene_RPP(
  polyA = QpolyA@polyA,
  sample_names = rownames(colData)
)

apa_by_condition <- lapply(c("Fip1", "Fip2"), function(treatment) {
  apa <- Quantify.GeneAPA(
    QpolyA,
    colData,
    contrast = c("condition", "NC", treatment)
  )

  delta <- compute_delta_RPP(
    polyA_rank = gene_RPP,
    colData = colData,
    control_cond = "NC",
    treat_cond = treatment
  )

  result <- dplyr::left_join(apa, delta, by = "gene_id")
  result$contrast <- paste(treatment, "vs NC")
  result
})

DEAPA_gene <- dplyr::bind_rows(apa_by_condition)

# Illustrative screening rule; see the P-value interpretation in section 3.6.
screened_APA <- dplyr::filter(
  DEAPA_gene,
  !is.na(pd),
  !is.na(p.value),
  pd > 0.1,
  p.value < 0.05
)

head(DEAPA_gene)
```

The joined table contains gene_id, pd, r, p.value, delta_RPP, and contrast. Check sample coverage before interpreting the direction or strength of a change.
