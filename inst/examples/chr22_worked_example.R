# APPLE 1.0.0 worked example using the bundled chromosome 22 BED files.
#
# Edit the reference and annotation paths before running. Results are written
# to APPLE-example-output in the current working directory.

library(APPLE)
library(dplyr)
library(tidyr)
library(ggplot2)

optional_packages <- c("DEXSeq", "BiocParallel")
missing_packages <- optional_packages[
  !vapply(optional_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop(
    "Install the optional Bioconductor packages first: ",
    paste(missing_packages, collapse = ", ")
  )
}

reference <- normalizePath(
  "/path/to/Homo_sapiens.GRCh38.dna.primary_assembly.fa",
  mustWork = TRUE
)
annotation <- normalizePath(
  "/path/to/Homo_sapiens.GRCh38.113.gtf",
  mustWork = TRUE
)

bed_dir <- system.file("extdata", "chr22", package = "APPLE")
if (!nzchar(bed_dir)) {
  stop("Bundled chromosome 22 BED files were not found. Reinstall APPLE.")
}

output_dir <- file.path(getwd(), "APPLE-example-output")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
cores <- if (.Platform$OS.type == "windows") 1L else 4L

control_group <- "NC"
treatment_groups <- c("Fip1", "Fip2")

# Load, clean, cluster, annotate, filter, and map tail lengths.
bed_files <- list.files(bed_dir, pattern = "\\.bed$", full.names = TRUE)
stopifnot(length(bed_files) == 6L)

QpolyA <- Load.PolyA(files = bed_files)
QpolyA <- Remove.IP(QpolyA, fasta = reference)
QpolyA <- Cluster.PolyA(QpolyA, max.gapwidth = 24, mc.cores = cores)
QpolyA <- Annotate.PolyA(QpolyA, gff = annotation)
QpolyA <- Filter.PolyA(QpolyA, min_count = 10, min_sample = 2)
QpolyA <- Map.Tail(QpolyA)

sample_names <- QpolyA@sample_names
conditions <- sub("-[12]$", "", sample_names)
stopifnot(
  length(sample_names) == 6L,
  all(conditions %in% c(control_group, treatment_groups))
)

sample_info <- data.frame(
  sample = sample_names,
  condition = conditions,
  lib_id = sample_names,
  stringsAsFactors = FALSE
)
colData <- data.frame(
  condition = factor(conditions, levels = c(control_group, treatment_groups)),
  type = factor("single"),
  row.names = sample_names
)

# Add gene names for readable result tables.
polyA <- QpolyA@polyA
gene_annotation <- as.data.frame(rtracklayer::import(annotation)) |>
  filter(type == "gene") |>
  select(gene_id, gene_name, gene_biotype) |>
  distinct()
polyA <- left_join(polyA, gene_annotation, by = "gene_id")
polyA$cluster_id <- rownames(QpolyA@polyA)

# Gene-level APA, PAS-level usage, and directional RPP scores.
APA_Gene <- polyA |>
  filter(type != "intergenic", !grepl("^ERCC", seqnames)) |>
  group_by(gene_id) |>
  filter(n() >= 2) |>
  ungroup()

polyA_PSU <- APA_Gene |>
  group_by(gene_id) |>
  mutate(across(all_of(sample_names), ~ .x / sum(.x))) |>
  ungroup()

gene_RPP <- compute_gene_RPP(
  polyA = QpolyA@polyA,
  sample_names = sample_names
)
samples_by_condition <- split(sample_names, conditions)

run_apa_contrast <- function(treatment) {
  control_samples <- samples_by_condition[[control_group]]
  treatment_samples <- samples_by_condition[[treatment]]
  selected_samples <- c(control_samples, treatment_samples)
  contrast_label <- paste0(treatment, " v.s. ", control_group)

  gene_result <- Quantify.GeneAPA(
    QpolyA,
    colData,
    contrast = c("condition", control_group, treatment)
  )

  dxd <- DEXSeq::DEXSeqDataSet(
    countData = APA_Gene[, selected_samples, drop = FALSE],
    sampleData = colData[selected_samples, , drop = FALSE],
    design = ~ sample + exon + condition:exon,
    featureID = APA_Gene$cluster_id,
    groupID = APA_Gene$gene_id,
    featureRanges = NULL,
    transcripts = NULL,
    alternativeCountData = NULL
  )
  bp <- if (.Platform$OS.type == "windows") {
    BiocParallel::SerialParam()
  } else {
    BiocParallel::MulticoreParam(workers = cores)
  }
  dxd <- DEXSeq::DEXSeq(dxd, BPPARAM = bp)

  gene_q <- DEXSeq::perGeneQValue(dxd)
  gene_q <- data.frame(
    gene_id = names(gene_q),
    qvalue = unname(gene_q),
    stringsAsFactors = FALSE
  )

  delta_rpp <- compute_delta_RPP(
    polyA_rank = gene_RPP,
    colData = colData,
    control_cond = control_group,
    treat_cond = treatment
  )

  gene_result <- gene_result |>
    left_join(gene_q, by = "gene_id") |>
    left_join(delta_rpp, by = "gene_id") |>
    mutate(contrast = contrast_label)

  pas_result <- as.data.frame(dxd) |>
    select(
      featureID,
      gene_id = groupID,
      exonBaseMean,
      dispersion,
      stat,
      pvalue,
      padj
    ) |>
    left_join(
      polyA_PSU |>
        mutate(
          delta_PSU =
            rowMeans(across(all_of(treatment_samples))) -
            rowMeans(across(all_of(control_samples)))
        ) |>
        select(featureID = cluster_id, delta_PSU),
      by = "featureID"
    ) |>
    mutate(contrast = contrast_label)

  list(gene = gene_result, pas = pas_result)
}

apa_results <- lapply(treatment_groups, run_apa_contrast)
DEAPA_gene <- bind_rows(lapply(apa_results, function(x) x[["gene"]])) |>
  left_join(gene_annotation, by = "gene_id")
DEAPA_PAS <- bind_rows(lapply(apa_results, function(x) x[["pas"]])) |>
  left_join(gene_annotation, by = "gene_id")

write.csv(
  DEAPA_gene,
  file.path(output_dir, "DEAPA_gene.csv"),
  quote = FALSE,
  row.names = FALSE
)
write.csv(
  DEAPA_PAS,
  file.path(output_dir, "DEAPA_PAS.csv"),
  quote = FALSE,
  row.names = FALSE
)

# Signed PD: positive values indicate distal shifts and negative values
# indicate proximal shifts when the annotation supports that interpretation.
GENE_sig <- DEAPA_gene |>
  mutate(
    pd_signed = sign(r) * pd,
    significance = case_when(
      abs(pd) >= 0.1 & qvalue < 0.05 & pd_signed > 0 ~ "Distal",
      abs(pd) >= 0.1 & qvalue < 0.05 & pd_signed < 0 ~ "Proximal",
      TRUE ~ "NS"
    )
  )
positive_q <- GENE_sig$qvalue[is.finite(GENE_sig$qvalue) & GENE_sig$qvalue > 0]
q_floor <- if (length(positive_q)) min(positive_q) else .Machine$double.xmin

p_deapa <- GENE_sig |>
  mutate(neg_log10_q = -log10(pmax(qvalue, q_floor))) |>
  ggplot(aes(x = pd_signed, y = neg_log10_q, color = significance)) +
  geom_point(size = 0.8, alpha = 0.5) +
  geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  facet_wrap(~ contrast, nrow = 1) +
  scale_color_manual(
    values = c(Distal = "#D55E00", Proximal = "#0072B2", NS = "#999999")
  ) +
  labs(
    title = "Differential polyA usage (DEAPA)",
    x = "Signed PD",
    y = expression(-log[10](italic(q)~value)),
    color = "Regulation"
  ) +
  theme_classic(base_size = 13)

ggsave(
  file.path(output_dir, "deapa-volcano.png"),
  p_deapa,
  width = 8,
  height = 5,
  dpi = 300
)

# Tail-length PCA and per-contrast tests.
tail_pca <- Tail.PCA(
  QpolyA,
  sample_info,
  summary_stat = "median",
  min_count_per_sample = 10,
  cores = cores,
  show_progress = TRUE
)
variance_explained <- tail_pca$pca$sdev^2 /
  sum(tail_pca$pca$sdev^2) * 100
pca_data <- as.data.frame(tail_pca$pca$x)
pca_data$sample <- rownames(pca_data)
pca_data <- left_join(pca_data, sample_info, by = "sample")

p_tail_pca <- ggplot(pca_data, aes(PC1, PC2, color = condition)) +
  geom_point(size = 4) +
  labs(
    title = "Poly(A) Length - PCA",
    x = sprintf("PC1 (%.1f%%)", variance_explained[1]),
    y = sprintf("PC2 (%.1f%%)", variance_explained[2])
  ) +
  theme_minimal(base_size = 13)

ggsave(
  file.path(output_dir, "tail-pca.png"),
  p_tail_pca,
  width = 7,
  height = 6,
  dpi = 300
)

run_tail_contrast <- function(treatment) {
  result <- Tail.DiffPair(
    QpolyA = QpolyA,
    sample_info = sample_info[
      sample_info$condition %in% c(control_group, treatment),
      ,
      drop = FALSE
    ],
    control_group = control_group,
    treatment_group = treatment,
    test_method = "t_test",
    min_mRNA_per_condition = 10,
    logscale = FALSE,
    mc.cores = cores
  ) |>
    mutate(
      sig = case_when(
        q_value < 0.05 & mean_diff > 15 ~ "Lengthening",
        q_value < 0.05 & mean_diff < -15 ~ "Shortening",
        TRUE ~ "NS"
      )
    ) |>
    left_join(
      polyA |>
        select(
          cluster_id,
          seqnames,
          gene_id,
          gene_name,
          gene_biotype,
          type
        ),
      by = "cluster_id"
    )

  output_name <- paste0("polyA_tail_diff.NC_", treatment, ".tsv")
  write.table(
    result,
    file.path(output_dir, output_name),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  positive_tail_q <- result$q_value[
    is.finite(result$q_value) & result$q_value > 0
  ]
  tail_q_floor <- if (length(positive_tail_q)) {
    min(positive_tail_q)
  } else {
    .Machine$double.xmin
  }
  plot_data <- result |>
    mutate(neg_log10_q = -log10(pmax(q_value, tail_q_floor)))

  p <- ggplot(plot_data, aes(mean_diff, neg_log10_q, color = sig)) +
    geom_point(size = 1.5, alpha = 0.6) +
    geom_vline(xintercept = c(-15, 15), linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    scale_color_manual(
      values = c(
        Lengthening = "#E41A1C",
        Shortening = "#377EB8",
        NS = "grey70"
      )
    ) +
    labs(
      title = paste(treatment, "vs", control_group),
      x = "Mean tail-length difference (nt)",
      y = expression(-log[10](italic(q)~value)),
      color = NULL
    ) +
    theme_bw(base_size = 13)

  ggsave(
    file.path(output_dir, paste0("tail-diff-", tolower(treatment), ".png")),
    p,
    width = 7,
    height = 6,
    dpi = 300
  )
  result
}

tail_results <- setNames(
  lapply(treatment_groups, run_tail_contrast),
  treatment_groups
)

message("Example completed. Results are in: ", normalizePath(output_dir))

