#' Differential alternative polyadenylation analysis with DEXSeq
#'
#' Perform pairwise differential APA analyses between one control condition
#' and one or more treatment conditions. Gene-level APPLE metrics are combined
#' with DEXSeq gene-level q-values, delta RPP, PAS-level DEXSeq statistics, and
#' delta PSU.
#'
#' @param QpolyA A `QuantifyPolyA` object generated after poly(A)-site
#'   clustering, annotation, and filtering.
#' @param colData A data frame containing sample information. Row names must
#'   match sample names in `QpolyA`, and a column named `condition` is required.
#' @param control A character string specifying the control condition.
#' @param treatment A character vector specifying one or more treatment
#'   conditions. If `NULL`, all conditions except `control` are analysed.
#' @param workers Number of parallel workers passed to
#'   [BiocParallel::MulticoreParam()].
#'
#' @return A list with three elements:
#'   \describe{
#'     \item{DEAPA_gene}{Gene-level results containing APPLE APA metrics,
#'       DEXSeq q-values, delta RPP, and the contrast label.}
#'     \item{DEAPA_PAS}{PAS-level DEXSeq statistics, delta PSU, and the
#'       contrast label.}
#'     \item{DEXSeq.Result}{A named list containing the fitted DEXSeq object
#'       for every treatment-versus-control comparison.}
#'   }
#'
#' @examples
#' \dontrun{
#' result <- DEAPA(
#'   QpolyA = QpolyA,
#'   colData = colData,
#'   control = "NC",
#'   treatment = c("EX1", "EX2"),
#'   workers = 10
#' )
#'
#' DEAPA_gene <- result$DEAPA_gene
#' DEAPA_PAS <- result$DEAPA_PAS
#' }
#'
#' @export
DEAPA <- function(QpolyA,
                  colData,
                  control,
                  treatment = NULL,
                  workers = 10) {
  if (!requireNamespace("DEXSeq", quietly = TRUE)) {
    stop(
      "Package 'DEXSeq' is required. Install it with ",
      "BiocManager::install('DEXSeq').",
      call. = FALSE
    )
  }
  if (!requireNamespace("BiocParallel", quietly = TRUE)) {
    stop(
      "Package 'BiocParallel' is required. Install it with ",
      "BiocManager::install('BiocParallel').",
      call. = FALSE
    )
  }

  if (!methods::is(QpolyA, "QuantifyPolyA")) {
    stop("QpolyA must be a QuantifyPolyA object.", call. = FALSE)
  }
  if (!is.data.frame(colData)) {
    stop("colData must be a data frame.", call. = FALSE)
  }
  if (!"condition" %in% colnames(colData)) {
    stop("colData must contain a 'condition' column.", call. = FALSE)
  }
  if (is.null(rownames(colData)) || any(!nzchar(rownames(colData)))) {
    stop("The row names of colData must contain sample names.", call. = FALSE)
  }
  if (anyDuplicated(rownames(colData))) {
    stop("The row names of colData must be unique.", call. = FALSE)
  }
  if (nrow(QpolyA@polyA) == 0L) {
    stop(
      "QpolyA@polyA is empty. Run clustering, annotation, and filtering first.",
      call. = FALSE
    )
  }

  sample_names <- rownames(colData)
  missing_samples <- setdiff(sample_names, QpolyA@sample_names)
  if (length(missing_samples)) {
    stop(
      "Samples in colData that are absent from QpolyA: ",
      paste(missing_samples, collapse = ", "),
      call. = FALSE
    )
  }
  missing_count_columns <- setdiff(sample_names, colnames(QpolyA@polyA))
  if (length(missing_count_columns)) {
    stop(
      "Sample columns absent from QpolyA@polyA: ",
      paste(missing_count_columns, collapse = ", "),
      call. = FALSE
    )
  }

  conditions <- as.character(colData$condition)
  if (anyNA(conditions) || any(!nzchar(conditions))) {
    stop("colData$condition cannot contain missing or empty values.", call. = FALSE)
  }
  available_conditions <- unique(conditions)

  if (length(control) != 1L || is.na(control) || !nzchar(control)) {
    stop("control must be one non-empty condition name.", call. = FALSE)
  }
  if (!control %in% available_conditions) {
    stop(
      "Control condition '", control, "' was not found in colData$condition.",
      call. = FALSE
    )
  }

  if (is.null(treatment)) {
    treatment <- setdiff(available_conditions, control)
  }
  treatment <- unique(as.character(treatment))
  if (!length(treatment) || anyNA(treatment) || any(!nzchar(treatment))) {
    stop("At least one valid treatment condition is required.", call. = FALSE)
  }
  if (control %in% treatment) {
    stop("The control condition cannot also be a treatment.", call. = FALSE)
  }
  unknown_treatments <- setdiff(treatment, available_conditions)
  if (length(unknown_treatments)) {
    stop(
      "Treatment conditions not found in colData$condition: ",
      paste(unknown_treatments, collapse = ", "),
      call. = FALSE
    )
  }

  if (length(workers) != 1L || !is.numeric(workers) ||
      is.na(workers) || !is.finite(workers) || workers < 1 ||
      workers != as.integer(workers)) {
    stop("workers must be a positive integer.", call. = FALSE)
  }
  workers <- as.integer(workers)

  polyA <- QpolyA@polyA
  required_columns <- c("type", "gene_id", "seqnames", "strand", "center")
  missing_columns <- setdiff(required_columns, colnames(polyA))
  if (length(missing_columns)) {
    stop(
      "QpolyA@polyA is missing required columns: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  cluster_id <- rownames(polyA)
  if (is.null(cluster_id) || any(!nzchar(cluster_id)) || anyDuplicated(cluster_id)) {
    stop("rownames(QpolyA@polyA) must be unique PAC identifiers.", call. = FALSE)
  }
  polyA$cluster_id <- cluster_id
  polyA$gene_id <- as.character(polyA$gene_id)

  APA_Gene <- polyA |>
    dplyr::filter(
      .data[["type"]] != "intergenic",
      !grepl("^ERCC", .data[["seqnames"]])
    ) |>
    dplyr::group_by(.data[["gene_id"]]) |>
    dplyr::filter(dplyr::n() >= 2L) |>
    dplyr::ungroup()

  if (!nrow(APA_Gene)) {
    stop("No genes with at least two eligible PACs were found.", call. = FALSE)
  }

  count_values <- as.matrix(APA_Gene[, sample_names, drop = FALSE])
  if (!is.numeric(count_values) || anyNA(count_values) ||
      any(!is.finite(count_values)) || any(count_values < 0) ||
      any(abs(count_values - round(count_values)) > sqrt(.Machine$double.eps))) {
    stop("PAC counts must be finite, non-negative integers.", call. = FALSE)
  }

  polyA_PSU <- APA_Gene |>
    dplyr::group_by(.data[["gene_id"]]) |>
    dplyr::mutate(
      dplyr::across(
        dplyr::all_of(sample_names),
        function(x) x / sum(x)
      )
    ) |>
    dplyr::ungroup()

  gene_RPP <- compute_gene_RPP(
    polyA = QpolyA@polyA,
    sample_names = sample_names
  )
  gene_RPP$gene_id <- as.character(gene_RPP$gene_id)

  samples_by_condition <- split(sample_names, conditions)
  BPPARAM <- BiocParallel::MulticoreParam(workers = workers)

  analyse_one_treatment <- function(treatment_name) {
    control_samples <- samples_by_condition[[control]]
    treatment_samples <- samples_by_condition[[treatment_name]]
    selected_samples <- c(control_samples, treatment_samples)
    contrast_label <- paste0(treatment_name, " v.s. ", control)

    gene_result <- Quantify.GeneAPA(
      QpolyA = QpolyA,
      colData = colData,
      contrast = c("condition", control, treatment_name)
    )
    gene_result <- dplyr::ungroup(gene_result)
    gene_result$gene_id <- as.character(gene_result$gene_id)

    count_matrix <- as.matrix(
      APA_Gene[, selected_samples, drop = FALSE]
    )
    storage.mode(count_matrix) <- "integer"

    dxd <- DEXSeq::DEXSeqDataSet(
      countData = count_matrix,
      sampleData = colData[selected_samples, , drop = FALSE],
      design = ~ sample + exon + condition:exon,
      featureID = APA_Gene$cluster_id,
      groupID = APA_Gene$gene_id,
      featureRanges = NULL,
      transcripts = NULL,
      alternativeCountData = NULL
    )
    dxd <- DEXSeq::DEXSeq(dxd, BPPARAM = BPPARAM)

    gene_qvalue <- DEXSeq::perGeneQValue(dxd)
    gene_qvalue <- data.frame(
      gene_id = as.character(names(gene_qvalue)),
      qvalue = unname(gene_qvalue),
      stringsAsFactors = FALSE
    )

    delta_rpp <- compute_delta_RPP(
      polyA_rank = gene_RPP,
      colData = colData,
      control_cond = control,
      treat_cond = treatment_name
    )
    delta_rpp$gene_id <- as.character(delta_rpp$gene_id)

    gene_result <- gene_result |>
      dplyr::left_join(gene_qvalue, by = "gene_id") |>
      dplyr::left_join(delta_rpp, by = "gene_id") |>
      dplyr::mutate(contrast = contrast_label)

    treatment_psu <- rowMeans(
      as.data.frame(polyA_PSU[, treatment_samples, drop = FALSE])
    )
    control_psu <- rowMeans(
      as.data.frame(polyA_PSU[, control_samples, drop = FALSE])
    )
    delta_psu <- data.frame(
      featureID = as.character(polyA_PSU$cluster_id),
      delta_PSU = treatment_psu - control_psu,
      stringsAsFactors = FALSE
    )

    dxd_table <- as.data.frame(dxd)
    pas_result <- data.frame(
      featureID = as.character(dxd_table$featureID),
      gene_id = as.character(dxd_table$groupID),
      exonBaseMean = dxd_table$exonBaseMean,
      dispersion = dxd_table$dispersion,
      stat = dxd_table$stat,
      pvalue = dxd_table$pvalue,
      padj = dxd_table$padj,
      stringsAsFactors = FALSE
    ) |>
      dplyr::left_join(delta_psu, by = "featureID") |>
      dplyr::mutate(contrast = contrast_label)

    list(gene = gene_result, pas = pas_result, dxd = dxd)
  }

  comparison_results <- lapply(treatment, analyse_one_treatment)
  names(comparison_results) <- paste0(treatment, "_vs_", control)

  list(
    DEAPA_gene = dplyr::bind_rows(
      lapply(comparison_results, function(x) x$gene)
    ),
    DEAPA_PAS = dplyr::bind_rows(
      lapply(comparison_results, function(x) x$pas)
    ),
    DEXSeq.Result = lapply(comparison_results, function(x) x$dxd)
  )
}

