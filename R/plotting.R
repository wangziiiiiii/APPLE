#' Plot poly(A) tail-length PCA
#'
#' Reproduces the PCA plot used by the APPLE example workflow.
#'
#' @param x A result returned by [Tail.PCA()].
#' @param colors Named colors for sample conditions.
#' @param point_size Point size.
#' @return A `ggplot` object.
#' @export
Plot.TailPCA <- function(
    x,
    colors = c("NC" = "#E41A1C", "EX1" = "#4DAF4A", "EX2" = "#377EB8"),
    point_size = 4.5
) {
  if (!is.list(x) || is.null(x$pca) || is.null(x$sample_info)) {
    stop("x must be a result returned by Tail.PCA().", call. = FALSE)
  }

  variance_explained <- x$pca$sdev^2 / sum(x$pca$sdev^2) * 100
  plot_data <- as.data.frame(x$pca$x)
  plot_data$sample <- rownames(plot_data)
  plot_data <- merge(plot_data, x$sample_info, by = "sample", all.x = TRUE)

  conditions <- unique(as.character(plot_data$condition))
  missing_colors <- setdiff(conditions, names(colors))
  if (length(missing_colors) > 0) {
    colors <- c(
      colors,
      stats::setNames(
        grDevices::hcl.colors(length(missing_colors), "Dark 3"),
        missing_colors
      )
    )
  }

  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = .data[["PC1"]], y = .data[["PC2"]])
  ) +
    ggplot2::geom_point(
      ggplot2::aes(color = .data[["condition"]]),
      size = point_size,
      stroke = 1.2,
      alpha = 0.9
    ) +
    ggplot2::labs(
      title = "Poly(A) Length - PCA",
      x = paste0("PC", 1, " (", round(variance_explained[1], 1), "%)"),
      y = paste0("PC", 2, " (", round(variance_explained[2], 1), "%)")
    ) +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, size = 13, face = "bold"),
      legend.position = "right",
      axis.title = ggplot2::element_text(size = 12),
      legend.text = ggplot2::element_text(size = 12, colour = "black"),
      legend.key = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(color = "gray90", linewidth = 0.5),
      panel.grid.minor = ggplot2::element_line(color = "gray95", linewidth = 0.2),
      panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.5),
      plot.margin = ggplot2::margin(1, 1, 1, 1, "cm")
    )
}

#' Plot poly(A) tail-length density
#'
#' @param x A result returned by [Tail.PCA()].
#' @param colors Named colors for sample conditions.
#' @return A `ggplot` object.
#' @export
Plot.TailDensity <- function(
    x,
    colors = c("NC" = "#E64B35", "EX1" = "#4DBBD5", "EX2" = "#00A087")
) {
  if (!is.list(x) || is.null(x$filtered_data)) {
    stop("x must be a result returned by Tail.PCA().", call. = FALSE)
  }

  required <- c("cluster_id", "sample", "summary_value")
  missing <- setdiff(required, colnames(x$filtered_data))
  if (length(missing) > 0) {
    stop("x$filtered_data is missing: ", paste(missing, collapse = ", "), call. = FALSE)
  }

  tail_table <- data.table::dcast(
    data.table::as.data.table(x$filtered_data),
    cluster_id ~ sample,
    value.var = "summary_value",
    fun.aggregate = mean,
    fill = NA_real_
  )

  tail_long <- tidyr::pivot_longer(
    as.data.frame(tail_table),
    cols = -dplyr::all_of("cluster_id"),
    names_to = "sample",
    values_to = "tail_length"
  )
  tail_long <- tail_long[!is.na(tail_long$tail_length), , drop = FALSE]
  tail_long$condition <- stringr::str_extract(tail_long$sample, "^[^-]+")

  tail_med <- dplyr::summarise(
    dplyr::group_by(tail_long, .data[["cluster_id"]], .data[["condition"]]),
    median_tail = stats::median(.data[["tail_length"]], na.rm = TRUE),
    .groups = "drop"
  )

  conditions <- unique(as.character(tail_med$condition))
  missing_colors <- setdiff(conditions, names(colors))
  if (length(missing_colors) > 0) {
    colors <- c(
      colors,
      stats::setNames(
        grDevices::hcl.colors(length(missing_colors), "Dark 3"),
        missing_colors
      )
    )
  }

  ggplot2::ggplot(
    tail_med,
    ggplot2::aes(x = .data[["median_tail"]], color = .data[["condition"]])
  ) +
    ggplot2::geom_density(alpha = 0.25, linewidth = 1) +
    ggplot2::scale_color_manual(values = colors, name = "Condition") +
    ggplot2::labs(
      title = "Poly(A) tail length distribution",
      x = "Median poly(A) tail length (nt)",
      y = "Density"
    ) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(
      legend.position = c(0.95, 0.95),
      legend.justification = c(1, 1),
      legend.box = "horizontal",
      legend.background = ggplot2::element_rect(
        fill = scales::alpha("white", 0.8),
        color = NA
      )
    )
}

#' Plot differential poly(A) tail lengths
#'
#' @param result A result returned by [Tail.DiffPair()].
#' @param title Plot title.
#' @param control_group Control-group label.
#' @param q_cutoff Adjusted P-value cutoff.
#' @param mean_diff_cutoff Absolute mean-difference cutoff in nucleotides.
#' @return A `ggplot` object.
#' @export
Plot.TailVolcano <- function(
    result,
    title,
    control_group = "NC",
    q_cutoff = 0.05,
    mean_diff_cutoff = 15
) {
  required <- c("mean_diff", "q_value")
  missing <- setdiff(required, colnames(result))
  if (length(missing) > 0) {
    stop("result is missing: ", paste(missing, collapse = ", "), call. = FALSE)
  }

  result <- as.data.frame(result)
  result$sig <- dplyr::case_when(
    result$q_value < q_cutoff & result$mean_diff > mean_diff_cutoff ~ "Lengthening",
    result$q_value < q_cutoff & result$mean_diff < -mean_diff_cutoff ~ "Shortening",
    TRUE ~ "NS"
  )

  sig_counts <- dplyr::count(result, .data[["sig"]])
  count_for <- function(label) {
    value <- sig_counts$n[sig_counts$sig == label]
    if (length(value) == 0) 0 else value
  }
  lengthening_count <- count_for("Lengthening")
  shortening_count <- count_for("Shortening")
  ns_count <- count_for("NS")
  max_abs_d <- max(abs(result$mean_diff), na.rm = TRUE)

  ggplot2::ggplot(
    result,
    ggplot2::aes(
      x = .data[["mean_diff"]],
      y = -log10(.data[["q_value"]]),
      color = .data[["sig"]]
    )
  ) +
    ggplot2::geom_point(size = 1.8, alpha = 0.6) +
    ggplot2::geom_hline(
      yintercept = -log10(q_cutoff),
      linetype = "dashed",
      color = "gray50",
      linewidth = 0.6
    ) +
    ggplot2::geom_vline(
      xintercept = c(-mean_diff_cutoff, mean_diff_cutoff),
      linetype = "dashed",
      color = "gray50",
      linewidth = 0.6
    ) +
    ggplot2::annotate(
      "text", x = -max_abs_d * 0.95, y = 97,
      label = shortening_count, hjust = 0, size = 5,
      color = "#377EB8", fontface = "bold"
    ) +
    ggplot2::annotate(
      "text", x = max_abs_d * 0.95, y = 97,
      label = lengthening_count, hjust = 1, size = 5,
      color = "#E41A1C", fontface = "bold"
    ) +
    ggplot2::annotate(
      "text", x = 0, y = 97, label = ns_count,
      hjust = 0.5, size = 4.5, color = "gray50", fontface = "italic"
    ) +
    ggplot2::scale_x_continuous(limits = c(-max_abs_d, max_abs_d)) +
    ggplot2::scale_y_continuous(limits = c(0, 100)) +
    ggplot2::scale_color_manual(
      values = c("Lengthening" = "#E41A1C", "Shortening" = "#377EB8", "NS" = "gray70")
    ) +
    ggplot2::labs(
      title = title,
      subtitle = paste0(control_group, " as control group"),
      x = "Mean Diff",
      y = expression(-Log[10](q - value))
    ) +
    ggplot2::guides(color = "none") +
    ggplot2::theme_bw(base_family = "Arial", base_size = 14) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, size = 16, face = "bold"),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 12, color = "gray40"),
      axis.title = ggplot2::element_text(size = 14),
      axis.text = ggplot2::element_text(size = 12)
    )
}

#' Plot PAS-expression PCA
#'
#' @param x A result returned by [DESeq2.PolyA()] or a DESeqDataSet.
#' @return A `ggplot` object.
#' @export
Plot.PASPCA <- function(x) {
  dds <- if (is.list(x) && !is.null(x$DESeq2.Result)) x$DESeq2.Result else x
  vsd <- DESeq2::vst(dds, blind = FALSE)
  plot_data <- DESeq2::plotPCA(vsd, returnData = TRUE, ntop = nrow(vsd))

  ggplot2::ggplot(plot_data, ggplot2::aes(x = .data[["PC1"]], y = .data[["PC2"]])) +
    ggplot2::geom_point(
      size = 3.5,
      ggplot2::aes(color = .data[["condition"]])
    ) +
    ggplot2::scale_shape_manual(values = c(15, 16, 17, 18, 19)) +
    cols4all::scale_color_discrete_c4a_cat("cols4all.friendly13") +
    ggplot2::xlab("PC1") +
    ggplot2::ylab("PC2") +
    ggplot2::ggtitle("PAS Expression - PCA") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, size = 16, face = "bold"),
      axis.title = ggplot2::element_text(size = 14, colour = "black", hjust = 0.5),
      axis.text = ggplot2::element_text(size = 12, colour = "black"),
      legend.text = ggplot2::element_text(size = 12, colour = "black"),
      legend.title = ggplot2::element_blank(),
      title = ggplot2::element_text(size = 12, colour = "black"),
      plot.margin = ggplot2::margin(1, 1, 1, 1, "cm")
    )
}

#' Plot differential PAS expression
#'
#' @param result A DESeq2 result table, typically returned by `lfcShrink()`.
#' @param treatment_group Treatment-group label.
#' @param control_group Control-group label.
#' @param padj_cutoff Adjusted P-value cutoff.
#' @param lfc_cutoff Absolute log2 fold-change cutoff.
#' @return A `ggplot` object.
#' @export
Plot.PASVolcano <- function(
    result,
    treatment_group,
    control_group = "NC",
    padj_cutoff = 0.05,
    lfc_cutoff = 1
) {
  result <- as.data.frame(result)
  required <- c("log2FoldChange", "padj")
  missing <- setdiff(required, colnames(result))
  if (length(missing) > 0) {
    stop("result is missing: ", paste(missing, collapse = ", "), call. = FALSE)
  }

  result <- result[!is.na(result$padj), , drop = FALSE]
  result$significance <- dplyr::case_when(
    result$padj < padj_cutoff & result$log2FoldChange > lfc_cutoff ~ "Up",
    result$padj < padj_cutoff & result$log2FoldChange < -lfc_cutoff ~ "Down",
    TRUE ~ "NS"
  )
  result$neg_log10_padj <- -log10(result$padj)

  n_down <- sum(result$significance == "Down")
  n_up <- sum(result$significance == "Up")
  n_ns <- sum(result$significance == "NS")
  max_abs_d <- max(abs(result$log2FoldChange), na.rm = TRUE)

  ggplot2::ggplot(
    result,
    ggplot2::aes(
      x = .data[["log2FoldChange"]],
      y = .data[["neg_log10_padj"]],
      color = .data[["significance"]]
    )
  ) +
    ggplot2::geom_point(size = 1.8, alpha = 0.6) +
    ggplot2::geom_hline(
      yintercept = -log10(padj_cutoff),
      linetype = "dashed",
      color = "gray50",
      linewidth = 0.6
    ) +
    ggplot2::geom_vline(
      xintercept = c(-lfc_cutoff, lfc_cutoff),
      linetype = "dashed",
      color = "gray50",
      linewidth = 0.6
    ) +
    ggplot2::annotate(
      "text", x = -max_abs_d * 0.95, y = 97,
      label = paste0("Down: ", n_down), hjust = 0,
      size = 5, color = "#1B9E77", fontface = "bold"
    ) +
    ggplot2::annotate(
      "text", x = 0, y = 97,
      label = paste0("NS: ", n_ns), hjust = 0.5,
      size = 4.5, color = "gray50", fontface = "italic"
    ) +
    ggplot2::annotate(
      "text", x = max_abs_d * 0.95, y = 97,
      label = paste0("Up: ", n_up), hjust = 1,
      size = 5, color = "#D95F02", fontface = "bold"
    ) +
    ggplot2::scale_x_continuous(limits = c(-max_abs_d, max_abs_d)) +
    ggplot2::scale_y_continuous(limits = c(0, 100)) +
    ggplot2::scale_color_manual(
      values = c("Down" = "#1B9E77", "NS" = "gray70", "Up" = "#D95F02")
    ) +
    ggplot2::labs(
      title = "PAS Expression Change",
      subtitle = paste0(control_group, " as control group"),
      x = bquote(log[2] ~ Fold ~ Change ~ (.(treatment_group) / .(control_group))),
      y = expression(-Log[10](P - adj))
    ) +
    ggplot2::guides(color = "none") +
    ggplot2::theme_bw(base_family = "Arial", base_size = 14) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, size = 16, face = "bold"),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 12, color = "gray40"),
      axis.title = ggplot2::element_text(size = 14),
      axis.text = ggplot2::element_text(size = 12)
    )
}

#' Plot gene-level differential APA
#'
#' @param result Gene-level APA result table.
#' @param pd_cutoff Absolute PD cutoff.
#' @param q_cutoff Adjusted P-value cutoff.
#' @param facet_columns Number of facet columns. By default, one per contrast.
#' @return A `ggplot` object.
#' @export
Plot.DEAPAVolcano <- function(
    result,
    pd_cutoff = 0.1,
    q_cutoff = 0.05,
    facet_columns = NULL
) {
  required <- c("pd", "r", "qvalue", "contrast")
  missing <- setdiff(required, colnames(result))
  if (length(missing) > 0) {
    stop("result is missing: ", paste(missing, collapse = ", "), call. = FALSE)
  }

  gene_sig <- dplyr::mutate(
    dplyr::ungroup(result),
    pd_signed = sign(.data[["r"]]) * .data[["pd"]],
    significance = dplyr::case_when(
      abs(.data[["pd"]]) >= pd_cutoff & .data[["qvalue"]] < q_cutoff & .data[["pd_signed"]] > 0 ~ "Distal",
      abs(.data[["pd"]]) >= pd_cutoff & .data[["qvalue"]] < q_cutoff & .data[["pd_signed"]] < 0 ~ "Proximal",
      TRUE ~ "NS"
    )
  )

  gene_counts <- tidyr::pivot_wider(
    dplyr::count(
      dplyr::filter(gene_sig, .data[["significance"]] != "NS"),
      .data[["contrast"]],
      .data[["significance"]]
    ),
    names_from = "significance",
    values_from = "n",
    values_fill = 0
  )
  for (nm in c("Distal", "Proximal")) {
    if (!nm %in% colnames(gene_counts)) gene_counts[[nm]] <- 0L
  }
  gene_counts$label <- paste0(
    "Distal: ", gene_counts$Distal,
    "\nProximal: ", gene_counts$Proximal
  )

  positive_q <- gene_sig$qvalue[gene_sig$qvalue > 0]
  min_nonzero_logq <- if (length(positive_q) > 0) {
    -log10(min(positive_q, na.rm = TRUE))
  } else {
    -log10(.Machine$double.xmin)
  }
  gene_sig$neg_log10_q <- ifelse(
    gene_sig$qvalue == 0,
    min_nonzero_logq,
    -log10(gene_sig$qvalue)
  )

  if (is.null(facet_columns)) {
    facet_columns <- length(unique(gene_sig$contrast))
  }

  ggplot2::ggplot(
    gene_sig,
    ggplot2::aes(x = .data[["pd_signed"]], y = .data[["neg_log10_q"]])
  ) +
    ggplot2::geom_point(
      ggplot2::aes(color = .data[["significance"]]),
      size = 0.8,
      alpha = 0.5
    ) +
    ggplot2::scale_color_manual(
      values = c("Distal" = "#D55E00", "Proximal" = "#0072B2", "NS" = "#999999"),
      name = "Regulation"
    ) +
    ggplot2::geom_vline(
      xintercept = c(-pd_cutoff, pd_cutoff),
      linetype = "dashed",
      color = "grey40",
      linewidth = 0.5
    ) +
    ggplot2::geom_hline(
      yintercept = -log10(q_cutoff),
      linetype = "dashed",
      color = "grey40",
      linewidth = 0.5
    ) +
    ggplot2::geom_text(
      data = gene_counts,
      ggplot2::aes(x = Inf, y = Inf, label = .data[["label"]]),
      inherit.aes = FALSE,
      hjust = -0.1,
      vjust = 1.5,
      size = 3.8,
      color = "grey20"
    ) +
    ggplot2::facet_wrap(ggplot2::vars(.data[["contrast"]]), ncol = facet_columns) +
    ggplot2::labs(
      x = expression("Signed PD" ~ (sign(r) %*% PD)),
      y = expression(-log[10](italic(q) ~ value)),
      title = "Differential polyA usage (DEAPA) volcano plot"
    ) +
    ggplot2::theme_classic(base_size = 16) +
    ggplot2::theme(
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(size = 13, face = "bold"),
      axis.title = ggplot2::element_text(size = 14),
      axis.text = ggplot2::element_text(size = 11, color = "grey30"),
      legend.position = "bottom",
      legend.title = ggplot2::element_text(size = 13),
      legend.text = ggplot2::element_text(size = 12),
      plot.title = ggplot2::element_text(size = 16, face = "bold", hjust = 0.5),
      panel.border = ggplot2::element_rect(color = "grey70", fill = NA, linewidth = 0.5),
      panel.spacing = grid::unit(1.5, "lines")
    )
}

#' Plot counts of differential APA events
#'
#' @param result Gene- or PAS-level APA result table.
#' @param level Either `"gene"` for Distal/Proximal genes or `"PAS"` for Up/Down PASs.
#' @param effect_cutoff Absolute PD or delta-PSU cutoff.
#' @param adjusted_p_cutoff Adjusted P-value cutoff.
#' @return A `ggplot` object.
#' @export
Plot.DEAPACounts <- function(
    result,
    level = c("gene", "PAS"),
    effect_cutoff = 0.1,
    adjusted_p_cutoff = 0.05
) {
  level <- match.arg(level)
  result <- dplyr::ungroup(result)

  if (level == "gene") {
    required <- c("pd", "r", "qvalue", "contrast")
    missing <- setdiff(required, colnames(result))
    if (length(missing) > 0) {
      stop("result is missing: ", paste(missing, collapse = ", "), call. = FALSE)
    }

    result <- dplyr::mutate(
      result,
      pd_signed = sign(.data[["r"]]) * .data[["pd"]],
      significance = dplyr::case_when(
        abs(.data[["pd"]]) >= effect_cutoff & .data[["qvalue"]] < adjusted_p_cutoff & .data[["pd_signed"]] > 0 ~ "Distal",
        abs(.data[["pd"]]) >= effect_cutoff & .data[["qvalue"]] < adjusted_p_cutoff & .data[["pd_signed"]] < 0 ~ "Proximal",
        TRUE ~ "NS"
      )
    )
    classes <- c("Distal", "Proximal")
    colors <- c("Distal" = "#D55E00", "Proximal" = "#0072B2")
    y_label <- "Number of significant genes"
    plot_title <- "Number of DEAPA genes per contrast"
    legend_position <- c(0.55, 0.65)
  } else {
    required <- c("delta_PSU", "padj", "contrast")
    missing <- setdiff(required, colnames(result))
    if (length(missing) > 0) {
      stop("result is missing: ", paste(missing, collapse = ", "), call. = FALSE)
    }

    result <- dplyr::mutate(
      result,
      significance = dplyr::case_when(
        .data[["delta_PSU"]] >= effect_cutoff & .data[["padj"]] < adjusted_p_cutoff ~ "Up",
        .data[["delta_PSU"]] <= -effect_cutoff & .data[["padj"]] < adjusted_p_cutoff ~ "Down",
        TRUE ~ "NS"
      )
    )
    classes <- c("Up", "Down")
    colors <- c("Up" = "#D55E00", "Down" = "#0072B2")
    y_label <- "Number of significant PASs"
    plot_title <- "Number of differential PAS usage per contrast"
    legend_position <- "top"
  }

  counts <- tidyr::pivot_wider(
    dplyr::count(
      dplyr::filter(result, .data[["significance"]] != "NS"),
      .data[["contrast"]],
      .data[["significance"]]
    ),
    names_from = "significance",
    values_from = "n",
    values_fill = 0
  )
  for (nm in classes) {
    if (!nm %in% colnames(counts)) counts[[nm]] <- 0L
  }

  counts_long <- tidyr::pivot_longer(
    counts,
    cols = dplyr::all_of(classes),
    names_to = "significance",
    values_to = "count"
  )

  ggplot2::ggplot(
    counts_long,
    ggplot2::aes(
      x = .data[["contrast"]],
      y = .data[["count"]],
      fill = .data[["significance"]]
    )
  ) +
    ggplot2::geom_col(
      position = ggplot2::position_dodge(width = 0.7),
      width = 0.6,
      linewidth = 0.3,
      color = "grey30"
    ) +
    ggplot2::scale_fill_manual(values = colors, name = "Regulation") +
    ggplot2::geom_text(
      ggplot2::aes(label = .data[["count"]]),
      position = ggplot2::position_dodge(width = 0.7),
      vjust = -0.4,
      size = 4
    ) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.15))) +
    ggplot2::labs(x = NULL, y = y_label, title = plot_title) +
    ggplot2::theme_classic(base_size = 14) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 12, color = "black"),
      axis.text.y = ggplot2::element_text(size = 11, color = "black"),
      axis.title.y = ggplot2::element_text(size = 14),
      plot.title = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5),
      legend.position = legend_position,
      legend.justification = c(0, 0),
      legend.background = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = 11),
      legend.title = ggplot2::element_text(size = 12),
      legend.key.size = grid::unit(0.7, "cm"),
      plot.margin = ggplot2::margin(t = 15, r = 10, b = 10, l = 10, unit = "pt")
    )
}
