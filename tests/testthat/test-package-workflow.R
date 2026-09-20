make_tail_fixture <- function() {
  samples <- c("NC-1", "NC-2", "TR-1", "TR-2")
  tails <- lapply(seq_along(samples), function(i) {
    data.frame(
      cluster_id = c("PAC1", "PAC2", "constant"),
      all_tail_lengths = c(
        paste(seq_len(12) + 10 + i * 3, collapse = ";"),
        paste(seq_len(12) + 40 + i^2, collapse = ";"),
        paste(seq_len(12) + 20, collapse = ";")
      )
    )
  })
  names(tails) <- samples
  list(
    object = methods::new("QuantifyPolyA", sample_names = samples,
                          cluster_tail_lengths = tails),
    metadata = data.frame(sample = samples,
                          condition = rep(c("NC", "TR"), each = 2),
                          lib_id = samples)
  )
}

test_that("installed namespace exposes the new public API", {
  exports <- getNamespaceExports("APPLE")
  expect_true(all(c(
    "Remove.IP", "Map.Tail", "Tail.PCA", "Tail.DiffPair", "DEAPA",
    "DEXSeq.PolyA",
    "Plot.TailPCA", "Plot.TailDensity", "Plot.TailVolcano",
    "Plot.PASPCA", "Plot.PASVolcano", "Plot.DEAPAVolcano",
    "Plot.DEAPACounts"
  ) %in% exports))
  expect_false(any(c("mapTail", "polyAlength", "tail_pca") %in% exports))
})

test_that("BED loading and tail mapping retain sample identifiers and values", {
  path <- tempfile(fileext = ".bed")
  on.exit(unlink(path))
  write.table(data.frame("chr1", "+", 100, 2, 80, "12,24"),
              path, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  object <- Load.PolyA(files = path)
  expect_identical(object@sample_names, tools::file_path_sans_ext(basename(path)))
  object@polyA <- data.frame(seqnames = "chr1", start = 95, end = 105,
                             strand = "+", row.names = "PAC1")
  mapped <- Map.Tail(object)
  expect_equal(mapped@cluster_tail_lengths[[1]]$all_tail_lengths, "12;24")
  expect_equal(mapped@cluster_tail_lengths[[1]]$cluster_id, "PAC1")
})

test_that("PD magnitude and strand-aware RPP direction match a known shift", {
  sites <- data.frame(type = rep("three_prime_UTR", 2), gene_id = "g1",
                      strand = "+", center = c(100, 200),
                      NC = c(80, 20), TR = c(20, 80))
  rpp <- compute_gene_RPP(sites, c("NC", "TR"))
  expect_equal(rpp$NC, 0.2)
  expect_equal(rpp$TR, 0.8)
  metadata <- data.frame(condition = c("NC", "TR"), row.names = c("NC", "TR"))
  delta <- compute_delta_RPP(rpp, metadata, "NC", "TR")
  expect_equal(delta$delta_RPP, 0.6)
  sites$strand <- "-"
  reverse <- compute_delta_RPP(compute_gene_RPP(sites, c("NC", "TR")), metadata, "NC", "TR")
  expect_equal(reverse$delta_RPP, -0.6)
  warn <- getOption("warn")
  on.exit(options(warn = warn))
  result <- dynamicsDetect(sites["NC"], sites["TR"], c("+", "+"), sites$center)
  expect_equal(result$pd, 0.6)
  expect_gt(result$r, 0)
  expect_equal(result$p.value, stats::chisq.test(as.matrix(sites[c("NC", "TR")]))$p.value)
})

test_that("tail PCA returns a scaled PCA and excludes constant clusters", {
  fixture <- make_tail_fixture()
  result <- Tail.PCA(fixture$object, fixture$metadata, cores = 1, show_progress = FALSE)
  expect_s3_class(result$pca, "prcomp")
  expect_equal(dim(result$matrix), c(2L, 4L))
  expect_setequal(rownames(result$pca$x), fixture$metadata$sample)
  expect_false("constant" %in% rownames(result$matrix))
})

test_that("pairwise tests retain raw summaries and apply BH correction", {
  fixture <- make_tail_fixture()
  for (method in c("t_test", "wilcoxon", "lmm")) {
    metadata <- fixture$metadata
    if (method == "lmm") metadata$lib_id <- NULL # Exercise documented lm fallback.
    result <- Tail.DiffPair(fixture$object, metadata, "NC", "TR",
                            test_method = method, mc.cores = 1)
    expect_equal(nrow(result), 3L)
    expect_equal(result$q_value, p.adjust(result$p_value, "BH"))
    expect_equal(result$fold_change, result$mean_treatment / result$mean_control)
    expect_equal(result$mean_diff, result$mean_treatment - result$mean_control)
    expect_gt(result$mean_diff[result$cluster_id == "PAC1"], 0)
  }
})

