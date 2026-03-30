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

##############################################
#           Data pre-processing            #
##############################################
#' Batch align fq files using minimap2 and filter/sort SAM files
#' @name minimap2
#' @param reference Reference genome file path
#' @param work_dir Working directory path containing fq files
#' @param threads Number of threads for minimap2, default is 4
#' @param filter_flags SAM flags to filter out, default is 2308
#' @return Vector of generated filtered SAM file paths
#' @export
minimap2 <- function(reference, work_dir, threads = 4, filter_flags = 2308) {

  # Check if reference genome exists
  if (!file.exists(reference)) {
    stop("Reference genome file does not exist: ", reference)
  }

  # Find fq files
  fq_files <- list.files(work_dir, pattern = "\\.(fq|fastq)(\\.gz)?$",
                         full.names = TRUE, ignore.case = TRUE)

  if (length(fq_files) == 0) {
    stop("No fq files found")
  }

  message("Found ", length(fq_files), " fq files")

  # Process each file
  filtered_sam_files <- sapply(fq_files, function(fq_file) {
    base_name <- sub("\\.(fq|fastq)(\\.gz)?$", "", fq_file)
    filtered_sam_file <- paste0(base_name, ".chop.aln.sam")

    message("Processing: ", basename(fq_file))

    # Step 1: Run minimap2 and generate SAM
    message("  Running minimap2...")
    minimap2_status <- system2("minimap2",
                               args = c("-y", "-ax", "splice", "-uf", "-t", threads,
                                        reference, fq_file),
                               stdout = filtered_sam_file)

    # Check if minimap2 failed
    if (minimap2_status != 0) {
      stop("MINIMAP2 ALIGNMENT FAILED for file: ", basename(fq_file))
    }

    # Step 2: Filter SAM file
    message("  Filtering SAM (removing flags ", filter_flags, ")...")

    # Create temporary file for filtered content
    temp_sam <- tempfile(pattern = "filtered_", fileext = ".sam")

    # Use samtools to filter and maintain SAM format
    filter_status <- system2("samtools",
                             args = c("view", "-h", "-F", filter_flags, filtered_sam_file),
                             stdout = temp_sam)

    # Check if samtools filter failed
    if (filter_status != 0) {
      if (file.exists(temp_sam)) file.remove(temp_sam)
      stop("SAMTOOLS FILTERING FAILED for file: ", basename(fq_file))
    }

    # Step 3: Sort the filtered SAM file
    message("  Sorting filtered SAM...")

    # Use samtools to sort and maintain SAM format
    sort_status <- system2("samtools",
                           args = c("sort", "--threads", threads, "-O", "SAM", temp_sam),
                           stdout = filtered_sam_file)

    # Check if samtools sort failed
    if (sort_status != 0) {
      if (file.exists(temp_sam)) file.remove(temp_sam)
      stop("SAMTOOLS SORTING FAILED for file: ", basename(fq_file))
    }

    # Clean up temporary file
    file.remove(temp_sam)

    message("  Successfully generated filtered SAM: ", basename(filtered_sam_file))

    return(filtered_sam_file)
  })

  message("Completed! Generated ", length(filtered_sam_files), " filtered SAM files")
  return(invisible(filtered_sam_files))
}

##############################################
#           Extract_polyAsite           #
##############################################
#' Process SAM files to detect polyA tails and generate processed BED files
#' @name Extract_polyAsite
#' @param work_dir Working directory containing SAM files
#' @param intron_max Maximum intron size, default is 50000
#' @param min_tail_length Minimum polyA tail length to consider, default is 6
#' @param sample_size Number of lines to sample for pt tag detection, default is 1000
#' @param bedtools_path Path to bedtools executable, default is "bedtools"
#' @param remove_temp_files Whether to remove temporary files, default is FALSE
#' @return List containing processing statistics for each file
#' @export
Extract_polyAsite <- function(work_dir, intron_max = 50000, min_tail_length = 6,
                              sample_size = 1000, bedtools_path = "bedtools",
                              remove_temp_files = FALSE) {

  # Load required libraries
  if (!requireNamespace("readr", quietly = TRUE)) {
    stop("Please install the 'readr' package")
  }
  if (!requireNamespace("stringr", quietly = TRUE)) {
    stop("Please install the 'stringr' package")
  }

  library(readr)
  library(stringr)

  # Set working directory
  setwd(work_dir)

  # Find all SAM files in the directory
  sam_files <- list.files(work_dir, pattern = "\\.sam$", full.names = TRUE, ignore.case = TRUE)

  if (length(sam_files) == 0) {
    stop("No SAM files found in directory: ", work_dir)
  }

  message("Found ", length(sam_files), " SAM files to process")

  # Initialize results list to store statistics for each file
  all_results <- list()

  # Process each SAM file
  for (sam_file_path in sam_files) {
    message("\nProcessing: ", basename(sam_file_path))

    # Generate output BED file name
    base_name <- tools::file_path_sans_ext(sam_file_path)
    initial_bed_file <- paste0(base_name, ".init.bed")   # 初始文件，标记为临时
    sorted_bed_file <- paste0(base_name, ".sort.bed")    # 排序文件，临时
    final_bed_file <- paste0(base_name, ".bed")          # 最终结果文件

    # Process single SAM file to generate initial BED
    file_results <- process_single_sam(
      sam_file_path = sam_file_path,
      bed_file_path = initial_bed_file,   # 注意这里使用 initial_bed_file
      intron_max = intron_max,
      min_tail_length = min_tail_length,
      sample_size = sample_size
    )
    # Step 4: Sort and process BED file
    message("  Sorting and processing BED file...")

    # Check if initial BED file was created and has content
    if (!file.exists(initial_bed_file) || file.info(initial_bed_file)$size == 0) {
      warning("Initial BED file is empty or missing for: ", basename(sam_file_path))
      next
    }

    # Sort BED file
    sort_command <- paste("sort -k 1,1 -k6,6r -k2,2n", initial_bed_file, ">", sorted_bed_file)
    sort_status <- system(sort_command)

    if (sort_status != 0) {
      stop("BED FILE SORTING FAILED for file: ", basename(sam_file_path))
    }

    # Check if sorted file was created
    if (!file.exists(sorted_bed_file) || file.info(sorted_bed_file)$size == 0) {
      stop("SORTED BED FILE ERROR: No sorted BED file created for: ", basename(sam_file_path))
    }
    message("  Running bedtools groupby...")
    groupby_status <- system2(bedtools_path,
                              args = c("groupby", "-i", sorted_bed_file,
                                       "-g", "1,6,2", "-c", "2,3,7", "-o", "count,mode,collapse"),
                              stdout = final_bed_file)

    if (groupby_status != 0) {
      stop("BEDTOOLS GROUPBY FAILED for file: ", basename(sam_file_path))
    }


    # Check if final BED file was created
    if (!file.exists(final_bed_file) || file.info(final_bed_file)$size == 0) {
      stop("FINAL BED FILE ERROR: No final BED file created for: ", basename(sam_file_path))
    }

    # Remove temporary files if requested
    if (remove_temp_files) {
      if (file.exists(initial_bed_file)) file.remove(initial_bed_file)
      if (file.exists(sorted_bed_file)) file.remove(sorted_bed_file)
    }

    # Add BED processing status to results
    file_results$bed_processing_success <- TRUE
    file_results$final_bed_file <- final_bed_file

    message("  Successfully generated processed BED file: ", basename(final_bed_file))

    # Store results
    all_results[[basename(sam_file_path)]] <- file_results
  }

  # Print summary for all files
  print_summary(all_results)

  return(invisible(all_results))
}

#' Process a single SAM file and generate BED output
#'
#' @param sam_file_path Path to input SAM file
#' @param bed_file_path Path to output BED file
#' @param intron_max Maximum intron size
#' @param min_tail_length Minimum polyA tail length
#' @param sample_size Number of lines to sample for pt tag detection
#'
#' @return List containing processing statistics
process_single_sam <- function(sam_file_path, bed_file_path, intron_max, min_tail_length, sample_size) {

  # First, check if the file contains pt tags by sampling lines
  message("  Checking for pt tags in file...")
  file_has_pt_tags <- check_pt_tags(sam_file_path, sample_size)

  if (file_has_pt_tags) {
    message("  File contains pt tags - using pt tag values for tail length")
  } else {
    message("  No pt tags found - using regex detection for tail length")
  }

  # Open files
  sam_file <- file(sam_file_path, "r")
  bed_file <- file(bed_file_path, "w")

  # Initialize counters
  stats <- list(
    total_mapping_records = 0,
    reads_flag0 = 0,
    reads_flag16 = 0,
    reads_flag0_pass = 0,
    reads_flag16_pass = 0,
    abnormal_as = 0,
    passed_reads = 0,
    failed_reads = 0,
    pt_tag_used = 0,
    regex_used = 0,
    file_has_pt_tags = file_has_pt_tags
  )

  start_time <- Sys.time()

  # Process SAM file line by line
  while (TRUE) {
    line = readLines(sam_file, n = 1)
    if (length(line) == 0) {
      break
    }

    # Skip header lines
    if (substr(line, start = 1, stop = 1) == '@') {
      next
    } else {
      stats$total_mapping_records <- stats$total_mapping_records + 1

      fields = str_split(line, '\t', simplify = TRUE)
      cigar_num = str_extract_all(fields[6], '\\d+', simplify = TRUE)
      cigar_num = as.numeric(cigar_num)
      cigar_mode = str_extract_all(fields[6], '\\D', simplify = TRUE)

      # Check for abnormal splicing
      idx = which(cigar_mode == 'N')
      if (length(idx)) {
        if (any(cigar_num[idx] > intron_max)) {
          stats$abnormal_as <- stats$abnormal_as + 1
          next
        }
        # Additional filtering conditions
        if (cigar_num[idx[1]] > 5000) {
          left_match = sum(cigar_num[cigar_mode[1:(idx[1]-1)] == 'M'])
          if (left_match < 30) {
            stats$abnormal_as <- stats$abnormal_as + 1
            next
          }
        }
        if (cigar_num[idx[length(idx)]] > 5000) {
          sub_idx = (idx[length(idx)]+1):length(cigar_mode)
          right_match = sum(cigar_num[sub_idx][cigar_mode[sub_idx] == 'M'])
          if (right_match < 30) {
            stats$abnormal_as <- stats$abnormal_as + 1
            next
          }
        }
      }

      # Calculate genomic coordinates
      idx = which(cigar_mode %in% c('M', 'D', 'N'))
      five_prime_end = as.numeric(fields[4]) - 1
      three_prime_end = as.numeric(fields[4]) + sum(cigar_num[idx])

      # Initialize variables
      tail_length <- 0
      has_polyA_tail <- FALSE
      method_used <- "none"

      # Use different detection methods based on file content
      if (file_has_pt_tags) {
        # File has pt tags - use pt tag method
        pt_tag <- fields[grepl("^pt:i:", fields)]
        if (length(pt_tag) > 0) {
          tail_length_str <- str_extract(pt_tag[1], "pt:i:(\\d+)")
          if (!is.na(tail_length_str)) {
            tail_length <- as.numeric(str_replace(tail_length_str, "pt:i:", ""))
            has_polyA_tail <- TRUE
            method_used <- "pt_tag"
            stats$pt_tag_used <- stats$pt_tag_used + 1

            # Update strand counts for pt tag
            if (fields[2] == '0') {
              stats$reads_flag0 <- stats$reads_flag0 + 1
              stats$reads_flag0_pass <- stats$reads_flag0_pass + 1
            } else if (fields[2] == '16') {
              stats$reads_flag16 <- stats$reads_flag16 + 1
              stats$reads_flag16_pass <- stats$reads_flag16_pass + 1
            }
          }
        }
      } else {
        # File has no pt tags - use regex method
        if (fields[2] == '0') {
          stats$reads_flag0 <- stats$reads_flag0 + 1
          if (cigar_mode[length(cigar_mode)] == 'S') {
            softclip_seq = str_sub(fields[10], str_length(fields[10]) - cigar_num[length(cigar_mode)] + 1)
            if (str_detect(softclip_seq, paste0('A{', min_tail_length, ',}'))) {
              a_tail <- str_extract(softclip_seq, 'A+')
              tail_length <- str_length(a_tail)
              has_polyA_tail <- TRUE
              method_used <- "regex_forward"
              stats$reads_flag0_pass <- stats$reads_flag0_pass + 1
              stats$regex_used <- stats$regex_used + 1
            }
          }
        } else if (fields[2] == '16') {
          stats$reads_flag16 <- stats$reads_flag16 + 1
          if (cigar_mode[1] == 'S') {
            softclip_seq = str_sub(fields[10], 1, cigar_num[1])
            if (str_detect(softclip_seq, paste0('T{', min_tail_length, ',}'))) {
              t_tail <- str_extract(softclip_seq, 'T+')
              tail_length <- str_length(t_tail)
              has_polyA_tail <- TRUE
              method_used <- "regex_reverse"
              stats$reads_flag16_pass <- stats$reads_flag16_pass + 1
              stats$regex_used <- stats$regex_used + 1
            }
          }
        }
      }

      # Output to BED file if polyA tail detected
      if (has_polyA_tail) {
        stats$passed_reads <- stats$passed_reads + 1
        if (fields[2] == '0') {
          writeLines(paste(c(fields[3], three_prime_end, five_prime_end, fields[1], fields[5], '+', tail_length), collapse = "\t"), bed_file, sep = '\n')
        } else if (fields[2] == '16') {
          writeLines(paste(c(fields[3], five_prime_end, three_prime_end, fields[1], fields[5], '-', tail_length), collapse = "\t"), bed_file, sep = '\n')
        }
      } else {
        stats$failed_reads <- stats$failed_reads + 1
      }
    }
  }

  end_time <- Sys.time()

  # Close files
  close(sam_file)
  close(bed_file)

  # Calculate processing time
  stats$processing_time <- round(as.numeric(difftime(end_time, start_time, units = "secs")), 2)
  stats$pass_rate <- round(stats$passed_reads / stats$total_mapping_records * 100, 2)

  return(stats)
}

#' Check if SAM file contains pt tags by sampling lines
#'
#' @param sam_file_path Path to SAM file
#' @param sample_size Number of lines to sample
#'
#' @return TRUE if pt tags found, FALSE otherwise
check_pt_tags <- function(sam_file_path, sample_size = 1000) {
  sam_file <- file(sam_file_path, "r")
  lines_checked <- 0
  pt_tag_found <- FALSE

  while (lines_checked < sample_size) {
    line <- readLines(sam_file, n = 1)
    if (length(line) == 0) {
      break
    }

    # Skip header lines
    if (substr(line, start = 1, stop = 1) == '@') {
      next
    }

    # Check for pt tag in data lines
    fields <- str_split(line, '\t', simplify = TRUE)
    pt_tag <- fields[grepl("^pt:i:", fields)]

    if (length(pt_tag) > 0) {
      pt_tag_found <- TRUE
      break
    }

    lines_checked <- lines_checked + 1
  }

  close(sam_file)
  return(pt_tag_found)
}

#' Print summary of processing results
#'
#' @param results_list List containing processing statistics for all files
print_summary <- function(results_list) {
  cat("\n=== OVERALL PROCESSING SUMMARY ===\n")

  for (file_name in names(results_list)) {
    stats <- results_list[[file_name]]

    cat("\nFile: ", file_name, "\n")
    cat("  File has pt tags:", stats$file_has_pt_tags, "\n")
    cat("  Total mapping records:", stats$total_mapping_records, "\n")
    cat("  Abnormal splicing skipped:", stats$abnormal_as, "\n")
    cat("  Passed reads:", stats$passed_reads, "\n")
    cat("  Failed reads:", stats$failed_reads, "\n")
    cat("  Pass rate:", stats$pass_rate, "%\n")
    cat("  PT tag used:", stats$pt_tag_used, "\n")
    cat("  Regex used:", stats$regex_used, "\n")
    cat("  Processing time:", stats$processing_time, "seconds\n")

    # Add BED processing status
    if (!is.null(stats$bed_processing_success)) {
      cat("  BED processing:", ifelse(stats$bed_processing_success, "SUCCESS", "FAILED"), "\n")
      if (stats$bed_processing_success) {
        cat("  Final BED file:", basename(stats$final_bed_file), "\n")
      }
    }

    cat("  Detailed strand statistics:\n")
    cat("    Forward strand - total:", stats$reads_flag0, ", passed:", stats$reads_flag0_pass, "\n")
    cat("    Reverse strand - total:", stats$reads_flag16, ", passed:", stats$reads_flag16_pass, "\n")
  }
}

#' Class QuantifyPolyA.
#'
#' Class \code{QuantifyPolyA.} defines a poly(A) dataset.
#'
#' @name QuantifyPolyA-class
#' @rdname QuantifyPolyA-class
#' @slot sample_names A character vector storing names of samples.
#' @slot pre.polyA A list storing un-clustered poly(A) sites.
#' @slot polyA A data.frame storing PACs (Poly(A) Site Clusters) generated from the weighted density peak clustering.
#' @slot simple.clusters A data.frame storing PACs generated from the k-nt iteratively clustering.
#' @slot split.clusters A data.frame storing only the split PACs.
#' @slot tail_lengths A list storing poly(A) tail length information for each sample.
#' @slot cluster_tail_lengths A list storing poly(A) tail length information mapped to clusters for each sample.
#' @exportClass QuantifyPolyA
#'
setClass("QuantifyPolyA",slots=list(sample_names="character",
                                    pre.polyA="list",
                                    polyA='data.frame',
                                    simple.clusters='data.frame',
                                    split.clusters ='data.frame',
                                    tail_lengths='list',
                                    cluster_tail_lengths='list'))
#' Method show.
#'
#' @rdname show-methods
#' @aliases show,ANY-method
#' @param object A QuantifyPolyA object.
#' @importFrom utils head
setMethod("show",
          "QuantifyPolyA",
          function(object) {
            cat('Number of samples:',length(object@sample_names),'\n')
            cat('Number of PAS in each sample:\n')
            print(sapply(object@pre.polyA, nrow))

            if (nrow(object@polyA)>0){
              print(head(object@polyA[, 1:min(11, ncol(object@polyA)), drop=F],2))
              cat(sprintf('%s...[%d x %d]%s','@polyA',nrow(object@polyA), ncol(object@polyA),'\n'))
            }

            if (nrow(object@simple.clusters)>0){
              print(head(object@simple.clusters[, 1:8, drop=F],2))
              cat(sprintf('%s...[%d x %d]%s','@simple.clusters',nrow(object@simple.clusters), ncol(object@simple.clusters),'\n'))
            }

            if (nrow(object@split.clusters)>0){
              print(head(object@split.clusters[ , , drop=F],2))
              cat(sprintf('%s...[%d x %d]%s','@split.clusters',nrow(object@split.clusters), ncol(object@split.clusters),'\n'))
            }

            if (length(object@tail_lengths)>0){
              cat('Poly(A) tail length information stored for', length(object@tail_lengths), 'samples\n')
            }

            if (length(object@cluster_tail_lengths)>0){
              cat('Poly(A) tail length information mapped to clusters for', length(object@cluster_tail_lengths), 'samples\n')
            }
          }
)


##############################################
#           Load raw poly(A) data            #
##############################################
#' @title Load raw poly(A) site data
#' @description Load the raw poly(A) site information extracted from 3' end sequencing data.
#' @name Load.PolyA
#' @usage Load.PolyA(files,dir)
#' @param files A character vector specifying the names of files containing poly(A) site information, each file contains 6 columns: 'seqnames','strand','coord','score','five_prime_end','tail_lengths'.
#' @param dir A string specifying the directory of files containing poly(A) site information.
#' @return A QuantifyPolyA object containing all raw poly(A) site information.
#' @importFrom tools file_path_sans_ext
#' @importFrom utils read.table
#' @importFrom methods new
#' @export
#'
Load.PolyA <- function(files,dir){
  # Check parameters.
  if (missing(files)){
    if (missing(dir)){
      stop("Parameter 'files' or 'dir' should be provided!")
    } else{
      if (!dir.exists(dir)) stop(paste('Directory',dir,'does not exist!'))
      files = list.files(dir,'\\.bed$',full.names = T)
      if (length(files) == 0) stop(paste('No bed files was found in directory',dir,'!'))
    }
  } else{
    for (file in files) {
      if (!file.exists(file)) stop(paste('Input file',file,'does not exist!'))
    }
  }

  # Load poly(A) file.
  alt_names = file_path_sans_ext(basename(files))
  pre.polyA = list()
  tail_lengths = list()

  for (i in 1:length(files)) {
    data = read.table(file = files[i],sep = '\t',stringsAsFactors = FALSE)

    # Check if the file has 6 columns (new format with tail lengths)
    if (ncol(data) == 6) {
      colnames(data) = c('seqnames','strand','coord','score','five_prime_end','tail_lengths')

      # Extract poly(A) tail length information from the 6th column
      # The 6th column contains comma-separated tail lengths
      tail_info <- str_split(data$tail_lengths, ",", simplify = FALSE)

      # Store tail length information
      tail_lengths[[alt_names[i]]] <- data.frame(
        seqnames = data$seqnames,
        strand = data$strand,
        coord = data$coord,
        score = data$score,
        five_prime_end = data$five_prime_end,
        tail_lengths = data$tail_lengths,
        stringsAsFactors = FALSE
      )

      # For APA analysis, we only use the first 5 columns and ignore tail lengths
      pre.polyA[[alt_names[i]]] = data[,c('seqnames','strand','coord','score','five_prime_end')]
    } else {
      # Old format (4 or 5 columns)
      pre.polyA[[alt_names[i]]] = data
      if (ncol(data) == 4) {
        colnames(pre.polyA[[alt_names[i]]]) = c('seqnames','strand','coord','score')
        pre.polyA[[alt_names[i]]]$five_prime_end = NA
      } else if (ncol(data) == 5) {
        colnames(pre.polyA[[alt_names[i]]]) = c('seqnames','strand','coord','score','five_prime_end')
      }
    }
  }
  print('Load polyA files finished!')
  QpolyA = new("QuantifyPolyA",
               sample_names = alt_names,
               pre.polyA = pre.polyA,
               tail_lengths = tail_lengths,
               cluster_tail_lengths = list())
  return(QpolyA)
}

##############################################
#          Cluster poly(A) sites             #
##############################################
#' @title Generate PACs
#' @description Cluster poly(A) sites into clusters using a weighted density peak clustering algorithm.
#' @name Cluster.PolyA
#' @usage Cluster.PolyA(QpolyA, max.gapwidth = 24, mc.cores = 4)
#' @param QpolyA A QuantifyPolyA object containing clean poly(A) sites.
#' @param max.gapwidth A cutoff limiting the max distance between two adjacent poly(A) sites in a poly(A) site cluster (PAC), default value is 24.
#' @param mc.cores An integer indicating the number of processors/cores to perform the clustering, default value is 4.
#' @importFrom pbmcapply pbmcmapply
#' @importFrom dplyr bind_rows
#' @export
#'
Cluster.PolyA <- function(QpolyA, max.gapwidth = 24, mc.cores = 4){
  # Check parameters.
  if (class(QpolyA) != "QuantifyPolyA") stop(paste('QpolyA should be a QuantifyPolyA object!'))
  if (!is.numeric(max.gapwidth)) stop("'max.gapwidth' is not a number!")
  if (max.gapwidth<=0) stop("'max.gapwidth' should be larger than 0, the default value is 24!")
  if (!is.numeric(mc.cores)) stop("'mc.cores' is not a number!")
  if (mc.cores<=0) stop("'mc.cores' should be larger than 0, the default value is 4!")

  # Declare
  seqnames = strand = coord = score = five_prime_end = sum.wts = NULL

  # Combine data of different samples.
  polyA = bind_rows(QpolyA@pre.polyA) %>%
    group_by(seqnames,strand,coord) %>%
    summarise(score = sum(score),five_prime_end = round(median(five_prime_end)),.groups='drop')

  # Simple clustering by distance.
  print('Simple clustering by distance!')
  points.gr = buildGenomicRanges(seqname = polyA$seqnames, position = polyA$coord,
                                 score = polyA$score, strand = polyA$strand, five_prime_end = polyA$five_prime_end)

  simple.clusters = simpleCluster(points.gr, max.gapwidth = max.gapwidth)
  simple.clusters.df = as.data.frame(simple.clusters)
  simple.clusters.df$split_label = NA

  # Re-clustering of polyA sites with large width using weighted density peak calling algorithm.
  print('Re-clustering by weighted density peak clustering algorithm!')
  idx = which(simple.clusters.df$width > max.gapwidth)
  pos = extractList(points.gr@ranges@start,simple.clusters$revmap[idx])
  wts = extractList(points.gr$score,simple.clusters$revmap[idx])
  five_prime_end = extractList(points.gr$five_prime_end,simple.clusters$revmap[idx])

  split.clusters = pbmcmapply(findPeaks,pos,wts,five_prime_end,SIMPLIFY = F, mc.cores = mc.cores)

  lens = sapply(split.clusters,nrow)
  split.clusters.df = bind_rows(split.clusters[lens!=1], .id = "split_label")

  split.clusters.df$seqnames = rep(simple.clusters.df$seqnames[idx[lens!=1]],times = lens[lens!=1])
  split.clusters.df$strand = rep(simple.clusters.df$strand[idx[lens!=1]],times = lens[lens!=1])
  split.clusters.df$width = split.clusters.df$end - split.clusters.df$start + 1
  split.clusters.df = dplyr::rename(split.clusters.df,score=sum.wts)

  # Generate final clusters.
  polyA = simple.clusters.df[-idx[lens!=1],c("seqnames","start","end","width","strand","score","center","five_prime_end","split_label")]
  polyA = rbind(polyA,split.clusters.df[,c("seqnames","start","end","width","strand","score","center","five_prime_end","split_label")])

  QpolyA@simple.clusters = simple.clusters.df
  QpolyA@split.clusters = split.clusters.df
  QpolyA@polyA = polyA


  return(QpolyA)
}

##############################################
#    Map tail lengths to clusters            #
##############################################
#' @title Map tail lengths to clusters
#' @description Map poly(A) tail length information to PAS clusters
#' @name mapTailLengthsToClusters
#' @param QpolyA A QuantifyPolyA object containing clean poly(A) sites and PACs.
#' @param delimiter The delimiter to use for tail lengths (default: ";")
#' @return A QuantifyPolyA object with tail length information mapped to clusters.
#' @importFrom GenomicRanges GRanges findOverlaps
#' @importFrom IRanges IRanges
#' @importFrom dplyr group_by summarise
#' @export
#'
mapTail <- function(QpolyA, delimiter = ";") {
  cluster_tail_lengths <- list()

  clusters.gr <- GRanges(
    seqnames = QpolyA@polyA$seqnames,
    ranges = IRanges(start = QpolyA@polyA$start, end = QpolyA@polyA$end),
    strand = QpolyA@polyA$strand,
    cluster_id = rownames(QpolyA@polyA)
  )

  for (sample_name in QpolyA@sample_names) {
    if (sample_name %in% names(QpolyA@tail_lengths)) {
      tail_data <- QpolyA@tail_lengths[[sample_name]]

      tail.gr <- GRanges(
        seqnames = tail_data$seqnames,
        ranges = IRanges(start = tail_data$coord, end = tail_data$coord),
        strand = tail_data$strand,
        tail_lengths = tail_data$tail_lengths
      )

      overlaps <- findOverlaps(tail.gr, clusters.gr, ignore.strand = FALSE)

      if (length(overlaps) > 0) {
        mapping_df <- data.frame(
          cluster_id = clusters.gr$cluster_id[subjectHits(overlaps)],
          tail_lengths = tail_data$tail_lengths[queryHits(overlaps)],
          stringsAsFactors = FALSE
        )

        cluster_tails <- mapping_df %>%
          group_by(cluster_id) %>%
          summarise(
            all_tail_lengths = paste(
              sapply(tail_lengths, function(x) {
                if (grepl(",", x)) {

                  gsub(",", delimiter, x)
                } else {
                  x
                }
              }),
              collapse = delimiter
            ),
            .groups = 'drop'
          )

        cluster_tails$seqnames <- QpolyA@polyA[cluster_tails$cluster_id, "seqnames"]
        cluster_tails$start <- QpolyA@polyA[cluster_tails$cluster_id, "start"]
        cluster_tails$end <- QpolyA@polyA[cluster_tails$cluster_id, "end"]
        cluster_tails$strand <- QpolyA@polyA[cluster_tails$cluster_id, "strand"]
        cluster_tails <- cluster_tails[, c("cluster_id", "seqnames", "start", "end", "strand", "all_tail_lengths")]
        cluster_tail_lengths[[sample_name]] <- cluster_tails
      }
    }
  }

  QpolyA@cluster_tail_lengths <- cluster_tail_lengths
  return(QpolyA)
}
##############################################
#    Annotate and quantify polyA sites       #
##############################################
#' @title Annotate the PACs
#' @description Annotate the PACs based on existed genome annotations, and generate the PAC table of all experimental samples.
#' @name Annotate.PolyA
#' @usage Annotate.PolyA(QpolyA,gff,seq.levels=NA)
#' @param QpolyA A QuantifyPolyA object containing clean poly(A) sites and PACs.
#' @param gff A genome annotation file in GFF or GTF format, GTF format is recommended.
#' @param seq.levels A character vector of the names of Chromosomes/Contigs. This parameter is useful when the Chromosomes/Contigs names in genome file (used for short reads mapping) and genome annotation file are not consistent.
#' @importFrom S4Vectors countQueryHits aggregate
#' @export
#'
Annotate.PolyA <- function(QpolyA,gff,seq.levels=NA){
  # Check parameters.
  if (class(QpolyA) != "QuantifyPolyA") stop(paste('QpolyA should be a QuantifyPolyA object!'))
  if (!file.exists(gff)) stop(paste('Annotation file',gff,'does not exist!'))

  # Annotate polyA sites
  print('Annotate polyA sites!')
  polyA = polyAsite.annotation(QpolyA@polyA,gff,seq.levels)

  # Generate full table of polyA sites
  print('Generate full table of polyA sites!')
  polyA.gr = GRanges(seqnames = polyA$seqnames,ranges = IRanges(start = polyA$start, end = polyA$end),strand = polyA$strand)

  # Declare
  score = NULL

  for (name in QpolyA@sample_names) {
    points.gr = GRanges(seqnames = QpolyA@pre.polyA[[name]]$seqnames,
                        ranges = IRanges(start = QpolyA@pre.polyA[[name]]$coord, end = QpolyA@pre.polyA[[name]]$coord),
                        strand = QpolyA@pre.polyA[[name]]$strand,
                        score = QpolyA@pre.polyA[[name]]$score)

    hits = findOverlaps(polyA.gr,points.gr,ignore.strand=FALSE)
    agg = S4Vectors::aggregate(points.gr, hits, score=sum(score))
    polyA[,name] = 0
    polyA[countQueryHits(hits) > 0L,name] = agg$score
  }
  QpolyA@polyA = polyA
  return(QpolyA)
}

##############################################
#       Filter low count polyA sites         #
##############################################
#' @title Filter out PACs
#' @description Filter out PACs with low count of supporting reads across multiple samples.
#' @name Filter.PolyA
#' @usage Filter.PolyA(QpolyA,min_count = 10, min_sample = 1)
#' @param QpolyA A QuantifyPolyA object containing clean poly(A) sites and PAC table.
#' @param min_count A fixed value specifying the minimum read count in a PAC.
#' @param min_sample A fixed value specifying the minimum number of samples with min_count reads in a PAC.
#' @export
#'
Filter.PolyA <- function(QpolyA,min_count = 10,min_sample = 1){
  # Check parameters.
  if (class(QpolyA) != "QuantifyPolyA") stop(paste('QpolyA should be a QuantifyPolyA object!'))
  if (!is.numeric(min_count)) stop("'min_count' is not a number!")
  if (min_count<=0) stop("'min_count' should be larger than 0, the default value is 10!")
  if (!is.numeric(min_sample)) stop("'min_sample' is not a number!")
  if (min_sample<=0) stop("'min_sample' should be larger than 0, the default value is 1!")

  QpolyA@polyA = subset(QpolyA@polyA, rowSums(QpolyA@polyA[,QpolyA@sample_names]>=min_count) >= min_sample)
  rownames(QpolyA@polyA) = paste('PA',1:nrow(QpolyA@polyA),sep = '')
  return(QpolyA)
}

##############################################
#          construct genome ranges           #
##############################################
#' Construct genomic ranges
#' @name buildGenomicRanges
#' @usage buildGenomicRanges(seqname,position,score,strand)
#' @param seqname A vector containing the sequence (Chromosomes/Contigs) names of poly(A) sites.
#' @param position A numeric vector containing the genomic positions of poly(A) sites.
#' @param score A numeric vector containing the numbers of reads supporting each poly(A) site.
#' @param strand A character vector containing the strand information of poly(A) sites.
#' @return A Granges object.
#' @export
#'
buildGenomicRanges <- function(seqname,position,score,strand = '*',five_prime_end){
  points.gr = GRanges(seqnames = seqname,ranges = IRanges(start = position,width = 1),
                      strand = strand,score = score,five_prime_end = five_prime_end)
  points.gr = sort(points.gr)
}


##############################################
#       simple clustering by distance        #
##############################################
#' Cluster poly(A) sites into groups by a certain distance iterative.
#' @name simpleCluster
#' @param points.gr A GRanges object constructed from poly(A) sites.
#' @param max.gapwidth A cutoff to limit the maximum distance between two adjacent sites in a PAC.
#' @return A reduced GRanges object.
#' @importFrom GenomicRanges reduce
#' @importFrom IRanges extractList
#' @importFrom BiocGenerics start
#' @export
#'
simpleCluster <- function(points.gr,max.gapwidth=24){

  # Cluster points by distance
  range.gr = reduce(points.gr,min.gapwidth=max.gapwidth,with.revmap=T,ignore.strand=FALSE)

  # Sum the score
  range.gr$score = sum(extractList(points.gr$score,range.gr$revmap))
  range.gr$five_prime_end = round(vapply(
    extractList(points.gr$five_prime_end, range.gr$revmap),
    median,
    numeric(1)
  ))
  # Select the center
  idx = BiocGenerics::which.max(extractList(points.gr$score,range.gr$revmap))
  range.gr$center = start(points.gr)[idx+c(0,cumsum(lengths(range.gr$revmap))[1:(length(range.gr$revmap)-1)])]

  # Return result
  return(range.gr)
}


##############################################
#             findPeaks function             #
##############################################
#' Cluster the one-dimensional poly(A) sites into groups based on local density
#' @name findPeaks
#' @usage findPeaks(sub_pos,sub_wts,min_delta=24)
#' @param sub_pos The genomic coordinates of poly(A) sites.
#' @param sub_wts The numbers of reads support each poly(A) site in 'sub_pos'.
#' @param min_delta A cutoff to limit the minimum distance between two centers.
#' @return A tibble object with five columns 'group', 'start', 'end', 'sum.wts', 'center'.
#' @importFrom dplyr tibble %>% group_by summarise
#' @importFrom outliers scores
#' @export
#'
findPeaks <- function(sub_pos,sub_wts,sub_five_prime_end,min_delta=24){

  # Declare
  group = pos = wts = NULL

  ND <- length(sub_pos)

  if(ND<=2){
    res <- tibble(group=1,start=min(sub_pos),end=max(sub_pos),sum.wts=sum(sub_wts),center=sub_pos[which.max(sub_wts)])
    return(res)
  }

  dc <- 0.5
  #paste('Computing Rho with gaussian kernel of radius: ',dc,collapse = '')

  rho = rep(0,ND)
  #
  # Gaussian kernel
  #
  for(i in 1:(ND-1)){
    for(j in (i+1):ND){
      tmp_dist = abs(sub_pos[i]-sub_pos[j])
      tmp_wts = sub_wts[i]*sub_wts[j]
      rho[i]=rho[i]+exp(-(tmp_dist/dc)^2)*sub_wts[j];
      rho[j]=rho[j]+exp(-(tmp_dist/dc)^2)*sub_wts[i];
    }
  }
  for(i in 1:ND){
    rho[i]=rho[i]+sub_wts[i];
  }

  maxd <- abs(sub_pos[1] - sub_pos[ND])
  rho_sorted = sort(rho,decreasing = T,index.return=T)
  ordrho <- rho_sorted$ix
  rho_sorted <- rho_sorted$x

  delta <- rep(-1.,ND)
  nneigh <- rep(0,ND)

  for(ii in 2:ND){
    delta[ordrho[ii]] = maxd
    for(jj in 1:(ii-1)){
      tmp_dist = abs(sub_pos[ordrho[ii]]-sub_pos[ordrho[jj]])
      if(tmp_dist<=delta[ordrho[ii]]){
        delta[ordrho[ii]] = tmp_dist
        nneigh[ordrho[ii]] = ordrho[jj]
      }
    }
  }
  delta[ordrho[1]]=max(delta)

  decision.data <- data.frame(rho=rho,delta=delta,rho.delta=rho*delta)
  outlier <- rho>max(rho)/2 & delta>max(delta)/2 & delta>min_delta

  tmp <- decision.data$rho.delta
  tmp[which(outlier)] <- 0
  outlier <- (outlier | scores(tmp, type="chisq", prob=0.99))& delta>min_delta

  if(!sum(outlier)){
    res <- tibble(group=1,start=min(sub_pos),end=max(sub_pos),sum.wts=sum(sub_wts),center=sub_pos[which.max(sub_wts)],five_prime_end=round(median(sub_five_prime_end)))
    return(res)
  }

  NCLUST = 0
  cl = rep(-1,ND)
  icl = c()
  for(i in 1:ND){
    if(outlier[i]){
      NCLUST = NCLUST + 1
      cl[i] = NCLUST
      icl[NCLUST] = i
    }
  }
  #paste('NUMBER OF CLUSTERS: ',NCLUST,collapse = '')

  # Assignation
  for(i in 1:ND){
    if (cl[ordrho[i]]==-1){
      cl[ordrho[i]] = cl[nneigh[ordrho[i]]];
    }
  }

  center <- rep(1,ND)
  center[icl] <- 2
  res <- data.frame(pos= sub_pos,wts = sub_wts,group = cl,center=center,five_prime_end = sub_five_prime_end)
  res <- res %>% group_by(group) %>% summarise(start=min(pos),end=max(pos),sum.wts = sum(wts),center = pos[center==2],five_prime_end = round(median(five_prime_end)),.groups = 'drop')
}


##############################################
#            annotation function             #
##############################################
#' Annotate the poly(A) sites based on the genomic structures
#' @name polyAsite.annotation
#' @usage polyAsite.annotation(polyA, gff, seq.levels=NA)
#' @param polyA A data.frame containing the info of clustered poly(A) sites - PAC.
#' @param gff A genome annotation file in GFF or GTF format, GTF format is recommended.
#' @param seq.levels A character vector of the names of Chromosomes/Contigs. This parameter is useful when the Chromosomes/Contigs names in genome file (used for short reads mapping) and genome annotation file are not consistent.
#' @return The polyA data.frame will be return with three additional columns 'gene_id', 'distance, 'type'.
#' @importFrom GenomicFeatures makeTxDbFromGFF cds exons genes threeUTRsByTranscript fiveUTRsByTranscript
#' @importFrom GenomicRanges GRanges findOverlaps follow distance
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors countQueryHits aggregate
#' @importFrom BiocGenerics width
#' @importFrom GenomeInfoDb seqlevels seqlevels<-
#' @export
#'
polyAsite.annotation <- function(polyA, gff, seq.levels=NA){

  # build granges object for polyA
  polyA.gr = GRanges(seqnames = polyA$seqnames,ranges = IRanges(start = polyA$start, end = polyA$end),strand = polyA$strand)

  # import annotation info from GFF/GTF file
  txdb = makeTxDbFromGFF(gff)
  if(!sum(is.na(seq.levels))){
    seqlevels(txdb) = seq.levels
  }

  # extract specific regions
  genes.gr = genes(txdb) # gene
  threeUTRs.gr = unlist(threeUTRsByTranscript(txdb)) # 3' UTR
  fiveUTRs.gr = unlist(fiveUTRsByTranscript(txdb)) # 5' UTR
  cds.gr = cds(txdb) # CDS
  exons.gr = exons(txdb) # EXON

  threeUTR.mean.len = mean(width(threeUTRs.gr),na.rm = TRUE) # mean length of 3' UTR

  # annotation by gene
  # Declare
  gene_id = NULL
  hits.gene = findOverlaps(polyA.gr,genes.gr,ignore.strand=FALSE)
  agg <- S4Vectors::aggregate(genes.gr, hits.gene, gene_id=BiocGenerics::paste(gene_id,collapse=';'))
  polyA$gene_id[countQueryHits(hits.gene) > 0L] <- agg$gene_id

  # calculate the distance between polyA site and its most upstream gene (distance of polyA site with gene annotation in previous step was set to 0)
  hits.upstream = follow(polyA.gr,genes.gr,ignore.strand=FALSE)
  idx = !is.na(hits.upstream)
  polyA$distance[idx] = distance(polyA.gr[idx],genes.gr[hits.upstream[idx]],ignore.strand=FALSE) + 1
  polyA$distance[!is.na(polyA$gene_id)] = 0

  # set gene_id for polyA site without gene annotation in previous step
  idx = idx & is.na(polyA$gene_id) # choose elements without gene annotation and having a precedent gene
  polyA$gene_id[idx] = genes.gr$gene_id[hits.upstream[idx]]

  # set type to 'ext_3UTR' if polyA site locating in a intergenic region but within 'threeUTR.mean.len' bp of the downstream of a specific gene
  idx = polyA$distance <= threeUTR.mean.len*2 & polyA$distance > 0
  polyA$type[idx] = 'ext_3UTR'


  # start annotation
  # annotate 3' UTR
  idx = which(polyA$distance == 0)
  if(length(idx)){
    hits.threeUTR = findOverlaps(polyA.gr[idx],threeUTRs.gr,ignore.strand=FALSE)
    polyA$type[idx[unique(hits.threeUTR@from)]] = '3UTR'
  }

  # annotate 5' UTR
  idx = which(polyA$distance == 0 & is.na(polyA$type))
  if(length(idx)){
    hits.fiveUTR = findOverlaps(polyA.gr[idx],fiveUTRs.gr,ignore.strand=FALSE)
    polyA$type[idx[unique(hits.fiveUTR@from)]] = '5UTR'
  }

  # annotate CDS
  idx = which(polyA$distance == 0 & is.na(polyA$type))
  if(length(idx)){
    hits.cds = findOverlaps(polyA.gr[idx],cds.gr,ignore.strand=FALSE)
    polyA$type[idx[unique(hits.cds@from)]] = 'CDS'
  }

  # annotate exon
  idx = which(polyA$distance == 0 & is.na(polyA$type))
  if(length(idx)){
    hits.exon = findOverlaps(polyA.gr[idx],exons.gr,ignore.strand=FALSE)
    polyA$type[idx[unique(hits.exon@from)]] = 'exon'
  }

  # annotate intron
  idx = which(polyA$distance == 0 & is.na(polyA$type))
  if(length(idx)){
    polyA$type[idx] = 'intron'
  }

  # annotate intergenic
  idx = which(is.na(polyA$type))
  if(length(idx)){
    polyA$type[idx] = 'intergenic'
  }


  ####################
  # refine annotation
  # construct reads range
  reads.gr = GRanges(seqnames = polyA$seqnames,
                     ranges = IRanges(start = rowMins(cbind(polyA$start,polyA$five_prime_end)) ,
                                      end = rowMaxs(cbind(polyA$end,polyA$five_prime_end))),
                     strand = polyA$strand,
                     pac_id = rownames(polyA))
  # find overlap between reads and genes
  intersect_regions = join_overlap_intersect_directed(reads.gr,genes.gr)
  intersect_regions.df = as.data.frame(intersect_regions)
  genes.gr.df = as.data.frame(genes.gr)
  genes.gr.df = dplyr::rename(genes.gr.df,gene_width = width)
  intersect_regions.df = left_join(intersect_regions.df,genes.gr.df[,c('gene_id','gene_width')],by=c('gene_id'))

  # assign to the best match gene
  assign_res = intersect_regions.df %>% group_by(pac_id) %>% summarise(gene_id_2 = select_best_match(.data))

  polyA$pac_id = rownames(polyA)
  polyA = left_join(polyA,assign_res,by=c('pac_id'))

  idx = which(!(((polyA$gene_id==polyA$gene_id_2)&(polyA$type!='intergenic'))|(is.na(polyA$gene_id_2)&polyA$type=='intergenic')))

  for (i in idx) {
    if(polyA$gene_id[i]==polyA$gene_id_2[i]){
      polyA$type[i] = 'ext_3UTR'
    }else{
      polyA$gene_id[i] = polyA$gene_id_2[i]
      tmp_gr = as_granges(polyA[i,c('seqnames','start','end','strand')])
      sub_gene = genes.gr[polyA$gene_id_2[i]]
      if(isEmpty(find_overlaps(tmp_gr,sub_gene)) ){
        polyA$type[i] = 'ext_3UTR'
      }
    }
  }

  # return result
  return(polyA[,!(names(polyA) %in% c('pac_id','gene_id_2','five_prime_end'))])
}

select_best_match <- function(object){
  # only 1 match
  if(length(object$gene_id)==1){
    return(object$gene_id[1])
  }

  # >=2 matches
  # 1st consider the overall coverage difference
  coverage_diff = max(object$width) - object$width
  idx = which(coverage_diff<=1000)
  if(length(idx)>1){
    # 2nd consider the 3' end distance
    if(object$strand[1] == '+'){
      distance_diff = max(object$end) - object$end
    }else{
      distance_diff = object$start - min(object$start)
    }
    idx = which(distance_diff<=500)
    if(length(idx)>1){
      # 3rd consider gene body length and overlap length
      idx = which.min(abs(object$gene_width - object$width))
      return(object$gene_id[idx])
    }else{
      return(object$gene_id[idx])
    }
  }else{
    return(object$gene_id[idx])
  }
}


##############################################
#    Save tail length information            #
##############################################
#' @title Save tail length information
#' @description Save poly(A) tail length information for each sample separately
#' @name Save.TailLengths
#' @usage Save.TailLengths(QpolyA, output_dir)
#' @param QpolyA A QuantifyPolyA object containing tail length information.
#' @param output_dir A string specifying the directory to save tail length information.
#' @export
#'
Save.TailLengths <- function(QpolyA, output_dir) {
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Save raw tail length information
  for (sample_name in QpolyA@sample_names) {
    if (sample_name %in% names(QpolyA@tail_lengths)) {
      tail_data <- QpolyA@tail_lengths[[sample_name]]
      output_file <- file.path(output_dir, paste0(sample_name, "_raw_tail_lengths.txt"))
      write.table(tail_data, output_file, sep = "\t", row.names = FALSE, quote = FALSE)
      cat("Saved raw tail length information for", sample_name, "to", output_file, "\n")
    }
  }

  # Save cluster tail length information
  for (sample_name in QpolyA@sample_names) {
    if (sample_name %in% names(QpolyA@cluster_tail_lengths)) {
      cluster_tail_data <- QpolyA@cluster_tail_lengths[[sample_name]]
      output_file <- file.path(output_dir, paste0(sample_name, "_cluster_tail_lengths.txt"))
      write.table(cluster_tail_data, output_file, sep = "\t", row.names = FALSE, quote = FALSE)
      cat("Saved cluster tail length information for", sample_name, "to", output_file, "\n")
    }
  }
}

##############################################
#          Motif search function             #
##############################################
#' @title Motif search
#' @description Search canonical poly(A) signals at the upstream of poly(A) sites.
#' @name Motif.Search
#' @param Sites A data.frame containing the information of poly(A) sites.
#' @param fasta A string specifies the location and name of genome file in FASTA format.
#' @return A list containing the types and genomic locations of signals.
#' @importFrom stringr str_extract_all str_locate_all
#' @export
#'
Motif.Search <- function(Sites,fasta){
  # Check parameters.
  if (!file.exists(fasta)) stop(paste('Fasta file',fasta,'does not exist!'))

  # Construct a data.frame to build a bed file for extracting sequences
  seq.bed = data.frame(
    chr = as.character(Sites$seqnames),
    start = Sites$start,
    end = Sites$end,
    strand = as.character(Sites$strand),
    stringsAsFactors = FALSE
  )

  # Check index file of fasta file
  if (!file.exists(paste0(fasta,'.fai'))) indexFa(fasta)

  # Sort
  idx = order(seq.bed$chr,seq.bed$start,decreasing = FALSE)
  seq.bed.sort = seq.bed[idx,]

  # Extract sequences
  seq.df = get.fasta(seq.bed.sort,fasta = fasta,verbose = FALSE,check.chr = FALSE)

  Sites$seq[idx] = seq.df$sequence
  Sites$index[idx] = seq.df$index

  signal.str = c()
  signal.pos = c()
  signal.str[Sites$strand=='+'] = str_extract_all(Sites$seq[Sites$strand=='+'],'.ATAAA|A.TAAA|AA.AAA|AAT.AA|AATA.A|AATAA.', simplify = F)
  signal.pos[Sites$strand=='+'] = str_locate_all(Sites$seq[Sites$strand=='+'],'.ATAAA|A.TAAA|AA.AAA|AAT.AA|AATA.A|AATAA.')
  signal.str[Sites$strand=='-'] = str_extract_all(Sites$seq[Sites$strand=='-'],'.TTATT|T.TATT|TT.ATT|TTT.TT|TTTA.T|TTTAT.', simplify = F)
  signal.pos[Sites$strand=='-'] = str_locate_all(Sites$seq[Sites$strand=='-'],'.TTATT|T.TATT|TT.ATT|TTT.TT|TTTA.T|TTTAT.')

  # num = sum(sapply(plus.signal, length)>0) +  sum(sapply(minus.signal, length)>0)
  #
  # minus.signal = DNAStringSet(unlist(minus.signal))
  # minus.signal = as.character(reverseComplement(minus.signal))
  #
  # logo = ggplot() + geom_logo( c(unlist(plus.signal),minus.signal)) + theme_logo()
  # annotate('rect', xmin = 0, xmax = 7, ymin = -0.05, ymax = 1.3, alpha = .1, col='black', fill='yellow')

  return(list(strings = signal.str, positions = signal.pos))
}


##############################################
#          Perform DESeq2 analysis           #
##############################################
#' Perform basic differential expression analysis using DESeq2 package
#' @name DESeq2.PolyA
#' @usage DESeq2.PolyA(QpolyA,colData)
#' @param QpolyA A QuantifyPolyA object.
#' @param colData A data.frame object containing the information of experimental design.
#' @importFrom DESeq2 DESeq counts vst DESeqDataSetFromMatrix
#' @importFrom grDevices colorRampPalette
#' @importFrom factoextra fviz_pca_ind
#' @importFrom ggplot2 coord_fixed xlab ylab
#' @importFrom uwot umap
#' @importFrom FactoMineR PCA
#' @importFrom SummarizedExperiment assay
#' @importFrom stats dist
#' @importFrom S4Vectors isEmpty
#' @export
#'
DESeq2.PolyA <- function(QpolyA,colData){
  # Check parameters.
  if (class(QpolyA) != "QuantifyPolyA") stop(paste('QpolyA should be a QuantifyPolyA object!'))
  if (!is.data.frame(colData)) stop("'colData' is not a data.frame!")
  if (isEmpty(colData$condition)) stop("'colData' does not contain a 'condition' column!")
  if (!isEmpty(setdiff(rownames(colData),QpolyA@sample_names))) stop("Inconsistency between rownames of 'colData' and sample names!")

  # Declare
  condition = X1 = X2 = NULL

  # Build DESeqDataSet
  dds <- DESeqDataSetFromMatrix(countData = QpolyA@polyA[,rownames(colData)],
                                colData = colData,
                                design = ~ condition)
  # Run DESeq
  dds <- DESeq(dds)

  # Normalized table
  polyA.normailized = cbind(QpolyA@polyA[,c(1:11)],counts(dds,normalized=TRUE))

  ### Data transformations and visualization
  ## Count data transformations
  vsd <- vst(dds, blind=FALSE)

  ## Data quality assessment by sample clustering and visualization
  # Heatmap of the sample-to-sample distances
  # library(pheatmap)
  # sampleDists <- dist(t(assay(vsd)))
  # sampleDistMatrix <- as.matrix(sampleDists)
  # rownames(sampleDistMatrix) <- vsd$condition
  # colnames(sampleDistMatrix) <- NULL
  # colors <- colorRampPalette(rev(brewer.pal(9, "Blues")) )(255)
  # sample.heatmap = pheatmap(sampleDistMatrix,
  #                           clustering_distance_rows=sampleDists,
  #                           clustering_distance_cols=sampleDists,
  #                           col=colors)

  # Principal component plot of the samples
  # library(FactoMineR)
  # library(factoextra)
  pcaData <- PCA(t(assay(vsd)), graph = FALSE)
  sample.pca = fviz_pca_ind(pcaData,
                            geom.ind = "point", # show points only (nbut not "text")
                            col.ind = factor(colData$condition), # color by groups
                            addEllipses = FALSE, # Concentration ellipses
                            legend.title = "Groups"
  )

  # tSNE plot of the samples
  # library(Rtsne) # Load package
  # set.seed(100) # Sets seed for reproducibility
  # tsne_out <- Rtsne(as.matrix(t(assay(vsd))),perplexity = 3) # Run TSNE
  #
  # sample.tsne = ggplot(data.frame(tsne_out$Y,condition=colData$condition), aes(X1, X2, color=condition)) +
  #   geom_point(size=3) +
  #   xlab(paste0("tSNE1")) +
  #   ylab(paste0("tSNE2")) +
  #   coord_fixed()

  # umap plot of the samples
  # library(uwot)
  # run UMAP algorithm
  umap_out <- umap(as.matrix(t(assay(vsd))), n_neighbors = 2, init = "spca")

  sample.umap = ggplot(data.frame(umap_out,condition=colData$condition), aes(X1, X2, color=condition)) +
    geom_point(size=3) +
    xlab(paste0("Dim1")) +
    ylab(paste0("Dim2")) +
    coord_fixed()

  return(list(DESeq2.Result = dds,
              PCA.Plot = sample.pca,
              UMAP.Plot = sample.umap))
}


#########################################################
#    Quantify APA dynamics among split polyA sites      #
#########################################################
#' @title Quantify APA dynamics among split PACs
#' @description Quantify APA dynamics among split PACs identified by the weighted density peak clustering.
#' @name Quantify.SplitAPA
#' @usage Quantify.SplitAPA(QpolyA,colData,contrast)
#' @param QpolyA A QuantifyPolyA object generated after poly(A) site clustering and annotation.
#' @param colData A data.frame object containing the information of experimental design.
#' @param contrast A three element vector specifying the two samples to be compared, e.g. contrast=c("condition","Brain","UHR").
#' @return A tibble object containing APA dynamics metrics of all split PACs.
#' @importFrom dplyr filter n do
#' @importFrom S4Vectors isEmpty
#' @export
#'
Quantify.SplitAPA <- function(QpolyA,colData,contrast){
  # Check parameters.
  if (class(QpolyA) != "QuantifyPolyA") stop(paste('QpolyA should be a QuantifyPolyA object!'))
  if (!is.data.frame(colData)) stop("'colData' is not a data.frame!")
  if (isEmpty(colData$condition)) stop("'colData' does not contain a 'condition' column!")
  if (!isEmpty(setdiff(rownames(colData),QpolyA@sample_names))) stop("Inconsistency between rownames of 'colData' and sample names!")

  control = rownames(colData)[colData[,contrast[1]]==contrast[2]]
  if (isEmpty(control)) stop(paste('Sample name',contrast[2],'is not valid!'))
  treat = rownames(colData)[colData[,contrast[1]]==contrast[3]]
  if (isEmpty(treat)) stop(paste('Sample name',contrast[3],'is not valid!'))

  # Declare
  split_label = . = NULL

  # Select out split PACs
  split.Sites = QpolyA@polyA %>% filter(!is.na(split_label)) %>%
    group_by(split_label) %>%
    filter(n()>=2)

  res = split.Sites %>% group_by(split_label) %>% do(dynamicsDetect(.[,control],.[,treat],.$strand,.$center))
}


############################################################
#    Quantify APA dynamics among canonical polyA sites     #
############################################################
#' @title Quantify APA dynamics among canonical PACs
#' @description Quantify APA dynamics among canonical PACs of a specific gene, canonical PACs represent PACs located at 3' UTR.
#' @name Quantify.CanonicalAPA
#' @usage Quantify.CanonicalAPA(QpolyA,colData,contrast)
#' @param QpolyA A QuantifyPolyA object generated after poly(A) site clustering and annotation.
#' @param colData A data.frame object containing the information of experimental design.
#' @param contrast A three element vector specifying the two samples to be compared, e.g. contrast=c("condition","Brain","UHR").
#' @return A tibble object containing APA dynamics metrics of all APA genes.
#' @importFrom S4Vectors isEmpty
#' @export
#'
Quantify.CanonicalAPA <- function(QpolyA,colData,contrast){
  # Check parameters.
  if (class(QpolyA) != "QuantifyPolyA") stop(paste('QpolyA should be a QuantifyPolyA object!'))
  if (!is.data.frame(colData)) stop("'colData' is not a data.frame!")
  if (isEmpty(colData$condition)) stop("'colData' does not contain a 'condition' column!")
  if (!isEmpty(setdiff(rownames(colData),QpolyA@sample_names))) stop("Inconsistency between rownames of 'colData' and sample names!")

  control = rownames(colData)[colData[,contrast[1]]==contrast[2]]
  if (isEmpty(control)) stop(paste('Sample name',contrast[2],'is not valid!'))
  treat = rownames(colData)[colData[,contrast[1]]==contrast[3]]
  if (isEmpty(treat)) stop(paste('Sample name',contrast[3],'is not valid!'))

  # Declare
  gene_id = type = . = NULL

  canonical.Sites = QpolyA@polyA %>% filter(type %in% c('3UTR','ext_3UTR')) %>%
    group_by(gene_id) %>%
    filter(n()>=2)

  res = canonical.Sites %>% group_by(gene_id) %>% do(dynamicsDetect(.[,control],.[,treat],.$strand,.$center))
}


################################################################################
#    Quantify APA dynamics between canonical and non-canonical polyA sites     #
################################################################################
#' Quantify the APA dynamics between canonical and non-canonical PACs in a gene.
#' @name Quantify.CNCAPA
#' @usage Quantify.CNCAPA(QpolyA,colData,contrast)
#' @param QpolyA A QuantifyPolyA object generated after poly(A) site clustering and annotation.
#' @param colData A data.frame object containing the information of experimental design.
#' @param contrast A two element vector specifying the two samples to be compared, e.g. contrast=c("condition","Brain","UHR").
#' @return A tibble object containing APA dynamics metrics of all APA genes.
#' @importFrom dplyr vars full_join mutate summarise_at
#' @importFrom tidyselect all_of
#' @importFrom S4Vectors isEmpty
#' @export
#'
Quantify.CNCAPA <- function(QpolyA,colData,contrast){
  # Check parameters.
  if (class(QpolyA) != "QuantifyPolyA") stop(paste('QpolyA should be a QuantifyPolyA object!'))
  if (!is.data.frame(colData)) stop("'colData' is not a data.frame!")
  if (isEmpty(colData$condition)) stop("'colData' does not contain a 'condition' column!")
  if (!isEmpty(setdiff(rownames(colData),QpolyA@sample_names))) stop("Inconsistency between rownames of 'colData' and sample names!")

  control = rownames(colData)[colData[,contrast[1]]==contrast[2]]
  if (isEmpty(control)) stop(paste('Sample name',contrast[2],'is not valid!'))
  treat = rownames(colData)[colData[,contrast[1]]==contrast[3]]
  if (isEmpty(treat)) stop(paste('Sample name',contrast[3],'is not valid!'))

  # Declare
  type = gene_id = canonical.label = strand = center = . = NULL

  # Group PACs
  gene.Sites = QpolyA@polyA %>% filter(type != 'intergenic') %>%
    mutate(canonical.label = type %in% c('3UTR','ext_3UTR')) %>%
    group_by(gene_id,canonical.label)

  gene.Sites.Count = gene.Sites %>% summarise_at(vars(all_of(c(control,treat))),sum) %>%
    group_by(gene_id) %>%
    filter(n()>=2)

  gene.Sites.Info = gene.Sites %>% summarise(strand = strand[1],center = mean(center),.groups = 'drop') %>%
    group_by(gene_id) %>%
    filter(n()>=2)

  gene.Sites = full_join(gene.Sites.Info,gene.Sites.Count,by=c('gene_id','canonical.label'))

  res = gene.Sites %>% group_by(gene_id) %>% do(dynamicsDetect(.[,control],.[,treat],.$strand,.$center))
}


############################################################
#        Quantify APA dynamics across a whole gene         #
############################################################
#' Quantify the APA dynamics across a whole gene
#' @name Quantify.GeneAPA
#' @usage Quantify.GeneAPA(QpolyA,colData,contrast)
#' @param QpolyA A QuantifyPolyA object generated after poly(A) site clustering and annotation.
#' @param colData A data.frame object containing the information of experimental design.
#' @param contrast A two element vector specifying the two samples to be compared, e.g. contrast=c("condition","Brain","UHR").
#' @return A tibble object containing APA dynamics metrics of all APA genes.
#' @importFrom S4Vectors isEmpty
#' @export
#'
Quantify.GeneAPA <- function(QpolyA,colData,contrast){
  # Check parameters.
  if (class(QpolyA) != "QuantifyPolyA") stop(paste('QpolyA should be a QuantifyPolyA object!'))
  if (!is.data.frame(colData)) stop("'colData' is not a data.frame!")
  if (isEmpty(colData$condition)) stop("'colData' does not contain a 'condition' column!")
  if (!isEmpty(setdiff(rownames(colData),QpolyA@sample_names))) stop("Inconsistency between rownames of 'colData' and sample names!")

  control = rownames(colData)[colData[,contrast[1]]==contrast[2]]
  if (isEmpty(control)) stop(paste('Sample name',contrast[2],'is not valid!'))
  treat = rownames(colData)[colData[,contrast[1]]==contrast[3]]
  if (isEmpty(treat)) stop(paste('Sample name',contrast[3],'is not valid!'))

  # Declare
  type = gene_id = . = NULL

  gene.Sites = QpolyA@polyA %>% filter(type != 'intergenic') %>%
    group_by(gene_id) %>%
    filter(n()>=2)

  res = gene.Sites %>% group_by(gene_id) %>% do(dynamicsDetect(.[,control],.[,treat],.$strand,.$center))
}

################################################
    #       compute_gene_RPP       #
################################################
#' Compute gene-level relative polyadenylation proportion (RPP) scores
#'
#' This function processes APA (alternative polyadenylation) site data to
#' compute gene-level RPP scores. It filters out intergenic sites, retains
#' only genes with at least two APA sites, calculates a rank-based weight
#' (`ps_rank`) based on site position and strand, normalizes counts per
#' sample to proportions, weights each site by its rank, and finally sums
#' the weighted proportions per gene.
#'
#' @param polyA A data frame containing APA site information. Must include:
#'   - Columns named as in `sample_names` (counts per sample).
#'   - Columns specified by `type_col`, `gene_id_col`, `strand_col`, `center_col`.
#'   - Row names should be unique site identifiers (used as `cluster_id`).
#' @param sample_names Character vector of sample names (must be columns in `polyA`).
#' @param type_col Name of the column containing site type (default: "type").
#' @param gene_id_col Name of the column containing gene identifier (default: "gene_id").
#' @param strand_col Name of the column containing strand information (default: "strand").
#' @param center_col Name of the column containing genomic center position (default: "center").
#'
#' @return A data frame with one row per gene and columns:
#'   - `gene_id`: gene identifier (from `gene_id_col`).
#'   - One column per sample (named as in `sample_names`) containing the RPP score.
#'
#' @importFrom dplyr %>% filter group_by mutate if_else percent_rank summarise_at vars all_of
#' @importFrom rlang !! sym
#'
#' @examples
#' \dontrun{
#' # Assuming `polyA` is a data frame with APA site counts and metadata,
#' # and `sample_names` is a vector of sample column names.
#' gene_rpp <- compute_gene_RPP(polyA, sample_names)
#' }
#'
#' @export
compute_gene_RPP <- function(polyA,
                             sample_names,
                             type_col = "type",
                             gene_id_col = "gene_id",
                             strand_col = "strand",
                             center_col = "center") {

  # Input validation
  if (!is.data.frame(polyA)) stop("polyA must be a data frame")
  if (!all(sample_names %in% colnames(polyA))) {
    missing <- setdiff(sample_names, colnames(polyA))
    stop("Sample columns missing from polyA: ", paste(missing, collapse = ", "))
  }
  required_cols <- c(type_col, gene_id_col, strand_col, center_col)
  if (!all(required_cols %in% colnames(polyA))) {
    missing <- setdiff(required_cols, colnames(polyA))
    stop("Required columns missing from polyA: ", paste(missing, collapse = ", "))
  }

  # Make a copy to avoid modifying original
  rpp_data <- polyA

  # Add cluster_id from row names (if row names are unique site IDs)
  if (is.null(rownames(rpp_data))) {
    rpp_data$cluster_id <- seq_len(nrow(rpp_data))
  } else {
    rpp_data$cluster_id <- rownames(rpp_data)
  }

  # Filter and compute site-level weighted proportions
  rpp_data <- rpp_data %>%
    dplyr::filter(!!sym(type_col) != "intergenic") %>%
    dplyr::group_by(!!sym(gene_id_col)) %>%
    dplyr::filter(dplyr::n() > 1) %>%
    dplyr::mutate(ps_rank = dplyr::if_else(!!sym(strand_col) == "+",
                                           dplyr::percent_rank(!!sym(center_col)),
                                           dplyr::percent_rank(-!!sym(center_col)))) %>%
    dplyr::mutate_at(dplyr::vars(dplyr::all_of(sample_names)), ~ . / sum(.)) %>%
    dplyr::mutate_at(dplyr::vars(dplyr::all_of(sample_names)), ~ . * ps_rank)

  # Aggregate to gene level
  gene_rpp <- rpp_data %>%
    dplyr::group_by(!!sym(gene_id_col)) %>%
    dplyr::summarise_at(dplyr::vars(dplyr::all_of(sample_names)), ~ sum(., na.rm = TRUE))

  # Rename the grouping column to 'gene_id' for consistency
  colnames(gene_rpp)[1] <- "gene_id"

  return(gene_rpp)
}

#' Compute delta RPP from precomputed gene-level RPP scores
#'
#' This function calculates the difference in mean relative polyadenylation
#' proportion (RPP) between two conditions (e.g., treatment vs. control)
#' using a precomputed gene-level RPP matrix. The delta RPP is defined as
#' `mean(treat group RPP) - mean(control group RPP)`.
#'
#' @param polyA_rank A data frame containing gene-level RPP scores.
#'   Must include a column `gene_id` and one column per sample, named
#'   exactly as the sample identifiers in `colData`.
#' @param colData A data frame with sample metadata. Row names must be
#'   sample names (matching column names in `polyA_rank`). It must contain
#'   a column named `condition` that specifies group membership.
#' @param control_cond Character string, the name of the control condition
#'   as it appears in `colData$condition`.
#' @param treat_cond Character string, the name of the treatment condition
#'   as it appears in `colData$condition`.
#'
#' @return A data frame with two columns:
#'   \item{gene_id}{Gene identifier (from `polyA_rank$gene_id`).}
#'   \item{delta_RPP}{Difference in mean RPP (treat minus control).}
#'
#' @examples
#' \dontrun{
#' # Assuming polyA_rank and colData are available
#' delta <- compute_delta_RPP(polyA_rank, colData, "NC", "Fip1")
#' head(delta)
#' }
#'
#' @export
compute_delta_RPP <- function(polyA_rank, colData, control_cond, treat_cond) {
  # Input validation
  if (!is.data.frame(polyA_rank)) stop("polyA_rank must be a data frame")
  if (!is.data.frame(colData)) stop("colData must be a data frame")
  if (!"gene_id" %in% colnames(polyA_rank)) stop("polyA_rank must contain a 'gene_id' column")
  if (!"condition" %in% colnames(colData)) stop("colData must contain a 'condition' column")

  # Extract sample names for each condition
  control_samples <- rownames(colData)[colData$condition == control_cond]
  treat_samples   <- rownames(colData)[colData$condition == treat_cond]

  # Check that sample columns exist in polyA_rank
  if (!all(control_samples %in% colnames(polyA_rank))) {
    missing <- setdiff(control_samples, colnames(polyA_rank))
    stop("Control samples missing from polyA_rank: ", paste(missing, collapse = ", "))
  }
  if (!all(treat_samples %in% colnames(polyA_rank))) {
    missing <- setdiff(treat_samples, colnames(polyA_rank))
    stop("Treatment samples missing from polyA_rank: ", paste(missing, collapse = ", "))
  }

  # Compute row means and delta
  delta <- rowMeans(polyA_rank[, treat_samples, drop = FALSE]) -
    rowMeans(polyA_rank[, control_samples, drop = FALSE])

  # Return result as a data frame
  data.frame(gene_id = polyA_rank$gene_id, delta_RPP = delta, stringsAsFactors = FALSE)
}


################################################
#       Detect dynamics of poly(A) sites       #
################################################
#' Detect APA dynamics between two samples
#' @name dynamicsDetect
#' @param sample_x A data.frame object with multiple rows and columns, each row represents a PAC, each column represents a replicate.
#' @param sample_y A data.frame object with multiple rows and columns, each row represents a PAC, each column represents a replicate.
#' @param sample_strand A vector containing the strand information of PACs.
#' @param sample_center A vector containing the coordinates (centers) information of PACs.
#' @return A data.frame object containing three columns 'pd', 'r', 'p.value'.
#' @importFrom stats chisq.test
#' @importFrom dplyr pull
#' @export
#'
dynamicsDetect <- function(sample_x,sample_y,sample_strand,sample_center){
  options(warn = -1)

  status = 1
  if(sample_strand[1]=='-') status = -1

  n_x = ncol(sample_x)
  n_y = ncol(sample_y)
  name_x = colnames(sample_x)
  name_y = colnames(sample_y)
  diff_x_y = data.frame(x = rep(name_x,each = n_y),
                        y = rep(name_y,times = n_x),
                        pd = rep(NA,n_x*n_y),
                        r = rep(NA,n_x*n_y),
                        p.value = rep(NA,n_x*n_y), stringsAsFactors = FALSE)

  for (x in 1:n_x) {
    ind.sum_x = sum(sample_x[,x])
    if(ind.sum_x==0) next
    ind.pet_x = sample_x[,x]/ind.sum_x

    for (y in 1:n_y) {
      ind.sum_y = sum(sample_y[,y])
      if(ind.sum_y==0) next
      ind.pet_y = sample_y[,y]/ind.sum_y

      diff_x_y[(x-1)*n_y+y,'pd'] = sum(abs(ind.pet_x - ind.pet_y))/2
      diff_x_y[(x-1)*n_y+y,'r'] = FU(pull(sample_x,x),pull(sample_y,y),sample_center)*status

      test.data = cbind(sample_x[,x],sample_y[,y])
      test.data = test.data[which(rowSums(test.data)>0),]
      if(dim(test.data)[1]==1) {diff_x_y[(x-1)*n_y+y,'p.value'] = 1}
      else {diff_x_y[(x-1)*n_y+y,'p.value'] = chisq.test(test.data)$p.value}
    }
  }
  #diff_x_y
  data.frame(pd = mean(diff_x_y$pd), r = mean(diff_x_y$r), p.value = mean(diff_x_y$p.value))
}


################################################
#                FU function                   #
################################################
#' Calculate the Pearson correlation coefficient of two PAC vectors.
#' @name FU
#' @param s1 A numeric vector of PAC usages of sample X.
#' @param s2 A numeric vector of PAC usages of sample Y.
#' @param score A numeric vector containing the genomic coordinates of PACs.
#' @return A numeric value in a range of [-1,1].
#' @export
#'
FU <- function(s1,s2,score) {
  p1=s1/sum(s1+s2)
  p2=s2/sum(s1+s2)
  p=rbind(p1,p2)
  u=c(1,2)
  v=score
  ubar=sum(p1)*u[1]+sum(p2)*u[2]
  vbar=sum((p1+p2)*v)

  fz=0
  for (i in 1:2) {
    for (j in 1:length(v)) {
      fz=fz+(u[i]-ubar)*(v[j]-vbar)*p[i,j]
    }
  }

  fm1=0
  for (i in 1:2) {
    fm1=fm1+(u[i]-ubar)^2*sum(p[i,])
  }

  fm2=0
  for (j in 1:length(v)) {
    fm2=fm2+(v[j]-vbar)^2*sum(p[,j])
  }

  r=fz/sqrt(fm1*fm2)
  return(r)
}

#' Prepare tail length matrix for PCA
#' @name prepare_tail_length_matrix
#' @description Create sample x PAS cluster matrix with aggregated tail lengths
#' @param QpolyA QuantifyPolyA object
#' @param sample_info Sample metadata
#' @param aggregation_method Method to aggregate tail lengths: "mean", "median"
#' @return A matrix for PCA analysis
#'
prepare_tail_length_matrix <- function(QpolyA, sample_info, aggregation_method = "mean") {

  analysis_data <- prepare_cluster_tail_data_simple(QpolyA, sample_info, cores = 1)

  if (nrow(analysis_data) == 0) {
    stop("No data available for PCA")
  }

  cluster_means <- analysis_data %>%
    group_by(sample, cluster_id) %>%
    summarise(
      aggregated_tail_length = case_when(
        aggregation_method == "mean" ~ mean(tail_length),
        aggregation_method == "median" ~ median(tail_length),
        TRUE ~ mean(tail_length)
      ),
      .groups = 'drop'
    )

  tail_matrix <- cluster_means %>%
    pivot_wider(
      names_from = cluster_id,
      values_from = aggregated_tail_length,
      values_fill = NA
    ) %>%
    as.data.frame()

  rownames(tail_matrix) <- tail_matrix$sample
  tail_matrix$sample <- NULL

  tail_matrix <- as.matrix(tail_matrix)

  message("Created tail length matrix: ", nrow(tail_matrix), " samples x ",
          ncol(tail_matrix), " PAS clusters")

  return(tail_matrix)
}

perform_tail_length_pca <- function(tail_matrix, scale = TRUE, center = TRUE) {

  if (ncol(tail_matrix) < 2) {
    stop("Need at least 2 PAS clusters for PCA")
  }

  if (nrow(tail_matrix) < 3) {
    stop("Need at least 3 samples for meaningful PCA")
  }

  message("Performing PCA on ", nrow(tail_matrix), " samples and ",
          ncol(tail_matrix), " PAS clusters")

  pca_result <- prcomp(tail_matrix, scale = scale, center = center)


  variance_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2) * 100

  message("PCA completed. Top 5 PCs explain:")
  for (i in 1:min(5, length(variance_explained))) {
    message("  PC", i, ": ", round(variance_explained[i], 2), "%")
  }


  pca_result$variance_explained <- variance_explained
  pca_result$cumulative_variance <- cumsum(variance_explained)

  return(pca_result)
}



#' Handle missing values in tail length matrix
#' @name handle_missing_values
#' @description Impute or remove missing values for PCA
#' @param tail_matrix Tail length matrix
#' @param max_missing Maximum allowed missing proportion per cluster (default: 0.2)
#' @param impute_method Imputation method: "mean", "knn", "remove"
#' @return Cleaned matrix ready for PCA
#'
handle_missing_values <- function(tail_matrix, max_missing = 0.2, impute_method = "mean") {

  missing_prop <- apply(tail_matrix, 2, function(x) sum(is.na(x)) / length(x))

  clusters_to_keep <- missing_prop <= max_missing
  tail_matrix_clean <- tail_matrix[, clusters_to_keep, drop = FALSE]

  message("Removed ", sum(!clusters_to_keep), " PAS clusters with >",
          max_missing*100, "% missing values")
  message("Retained ", ncol(tail_matrix_clean), " PAS clusters")

  if (any(is.na(tail_matrix_clean))) {
    if (impute_method == "mean") {
      for (j in 1:ncol(tail_matrix_clean)) {
        col_means <- mean(tail_matrix_clean[, j], na.rm = TRUE)
        tail_matrix_clean[is.na(tail_matrix_clean[, j]), j] <- col_means
      }
      message("Imputed missing values with column means")

    } else if (impute_method == "knn" && requireNamespace("impute", quietly = TRUE)) {
      tail_matrix_clean <- impute::impute.knn(tail_matrix_clean)$data
      message("Imputed missing values using KNN")

    } else if (impute_method == "remove") {
      complete_cases <- complete.cases(tail_matrix_clean)
      tail_matrix_clean <- tail_matrix_clean[complete_cases, , drop = FALSE]
      message("Removed ", sum(!complete_cases), " samples with missing values")
    }
  }

  return(tail_matrix_clean)
}

#' Extract PCA results for visualization
#' @name extract_pca_results
#' @description Extract coordinates and metadata for PCA plotting
#' @param pca_result PCA result from prcomp
#' @param sample_info Sample metadata
#' @param n_pcs Number of principal components to extract (default: 5)
#' @return Data frame with PCA coordinates and sample metadata
#'
extract_pca_results <- function(pca_result, sample_info, n_pcs = 5) {

  pca_coords <- as.data.frame(pca_result$x)

  pcs_to_keep <- paste0("PC", 1:min(n_pcs, ncol(pca_coords)))
  pca_coords <- pca_coords[, pcs_to_keep, drop = FALSE]

  pca_coords$sample <- rownames(pca_coords)
  pca_data <- merge(pca_coords, sample_info, by = "sample", all.x = TRUE)

  attr(pca_data, "variance_explained") <- pca_result$variance_explained[1:n_pcs]
  attr(pca_data, "cumulative_variance") <- pca_result$cumulative_variance[1:n_pcs]

  return(pca_data)
}
#' Simple Poly(A) Tail Length PCA Pipeline
#'
#' A simple wrapper function that connects all steps of tail length PCA analysis.
#'
#' @param QpolyA A QuantifyPolyA object
#' @param sample_info A data.frame with sample metadata
#' @param aggregation_method Aggregation method for tail lengths (default: "mean")
#' @param max_missing Maximum missing proportion (default: 0.2)
#' @param impute_method Imputation method (default: "mean")
#' @param scale Whether to scale PCA (default: TRUE)
#' @param center Whether to center PCA (default: TRUE)
#' @param color_by Column to color by (default: "condition")
#' @param shape_by Column to shape by (default: NULL)
#' @param pc_x Which PC to plot on x-axis (default: 1)
#' @param pc_y Which PC to plot on y-axis (default: 2)
#'
#' @return A ggplot object with PCA plot
#' @export
#'
tail_pca <- function(QpolyA, sample_info,
                     aggregation_method = "mean",
                     max_missing = 0.2,
                     impute_method = "mean",
                     scale = TRUE,
                     center = TRUE
                     ) {

  # Step 1: Prepare tail length matrix
  tail_matrix <- prepare_tail_length_matrix(QpolyA, sample_info, aggregation_method)

  # Step 2: Handle missing values
  tail_matrix_clean <- handle_missing_values(tail_matrix, max_missing, impute_method)

  # Step 3: Perform PCA
  pca_result <- perform_tail_length_pca(tail_matrix_clean, scale, center)

  # Step 4: Extract results for plotting
  pca_data <- extract_pca_results(pca_result, sample_info)

  return(pca_data)
}

if (!require("lme4")) install.packages("lme4")
if (!require("lmerTest")) install.packages("lmerTest")
if (!require("dplyr")) install.packages("dplyr")
if (!require("pbmcapply")) install.packages("pbmcapply")

library(lme4)
library(lmerTest)
library(dplyr)
library(pbmcapply)

#' Optimized data preparation - simplified columns
#' @name prepare_cluster_tail_data_simple
#' @description Fast parallel processing with only essential columns
#' @param QpolyA A QuantifyPolyA object
#' @param sample_info Data frame with sample metadata
#' @param cores Number of CPU cores to use
#' @return A data frame with individual mRNA tail lengths (cluster_id, sample, condition, tail_length, lib_id)
#'
prepare_cluster_tail_data_simple <- function(QpolyA, sample_info,
                                             cores = 4,
                                             show_progress = TRUE) {

  if (!require("parallel")) install.packages("parallel")
  library(parallel)

  sample_names <- names(QpolyA@cluster_tail_lengths)

  message("Using ", cores, " CPU cores")
  message("Processing ", length(sample_names), " samples")

  if (is.null(QpolyA@cluster_tail_lengths)) {
    stop("QpolyA@cluster_tail_lengths is NULL")
  }

  sample_info_map <- sample_info
  if (!"lib_id" %in% colnames(sample_info_map)) {
    sample_info_map$lib_id <- NA_character_
  }

  process_sample_simple <- function(sample_name) {
    tryCatch({
      cluster_data <- QpolyA@cluster_tail_lengths[[sample_name]]

      if (is.null(cluster_data) || nrow(cluster_data) == 0) {
        return(NULL)
      }

      sample_meta <- sample_info_map[sample_info_map$sample == sample_name, ]
      if (nrow(sample_meta) == 0) {
        return(NULL)
      }

      condition_val <- sample_meta$condition[1]
      lib_id_val <- sample_meta$lib_id[1]

      result_list <- lapply(seq_len(nrow(cluster_data)), function(i) {
        tail_str <- cluster_data$all_tail_lengths[i]
        if (is.na(tail_str) || tail_str == "") return(NULL)

        tails <- as.numeric(strsplit(tail_str, ";", fixed = TRUE)[[1]])
        valid_tails <- tails[!is.na(tails) & tails > 0]
        if (length(valid_tails) == 0) return(NULL)

        data.frame(
          cluster_id = rep(cluster_data$cluster_id[i], length(valid_tails)),
          sample = rep(sample_name, length(valid_tails)),
          condition = rep(condition_val, length(valid_tails)),
          tail_length = valid_tails,
          lib_id = rep(lib_id_val, length(valid_tails)),
          stringsAsFactors = FALSE
        )
      })

      valid_items <- result_list[!sapply(result_list, is.null)]
      if (length(valid_items) == 0) return(NULL)

      do.call(rbind, valid_items)

    }, error = function(e) {
      message("Error in sample ", sample_name, ": ", e$message)
      return(NULL)
    })
  }


  message("Starting parallel processing...")

  if (show_progress && requireNamespace("pbmcapply", quietly = TRUE)) {
    results <- pbmcapply::pbmclapply(
      sample_names,
      process_sample_simple,
      mc.cores = cores,
      ignore.interactive = !interactive()
    )
  } else {
    results <- mclapply(sample_names, process_sample_simple, mc.cores = cores)
  }

  valid_results <- results[!sapply(results, is.null)]

  if (length(valid_results) == 0) {
    warning("No valid data obtained from any sample")
    return(data.frame())
  }

  message("Merging results from ", length(valid_results), " samples")

  if (requireNamespace("data.table", quietly = TRUE)) {
    final_data <- data.table::rbindlist(valid_results, use.names = TRUE, fill = TRUE)
    final_data <- as.data.frame(final_data)
  } else {
    final_data <- dplyr::bind_rows(valid_results)
  }

  message("Processing completed! Total rows: ", nrow(final_data))
  return(final_data)
}

##############################################
#           Statistical test function    #
##############################################

#' Perform t-test using all individual mRNA tail lengths for multiple treatments vs NC
#' @name perform_t_test_all_mRNA_multi
#' @param cluster_data Data for a single cluster
#' @param logscale Whether to log2 transform tail lengths
#' @param control_group Name of control group (default: "NC")
#' @return A list of results for each treatment vs control comparison
perform_t_test_all_mRNA_multi <- function(cluster_data, logscale = TRUE, control_group = "NC") {

  if (logscale) {
    cluster_data$tail_length <- log2(cluster_data$tail_length)
  }

  # 获取所有条件
  all_conditions <- unique(cluster_data$condition)

  # 检查控制组是否存在
  if (!control_group %in% all_conditions) {
    stop("Control group '", control_group, "' not found in data. Available conditions: ",
         paste(all_conditions, collapse = ", "))
  }

  # 获取处理组（所有非控制组的条件）
  treatment_groups <- setdiff(all_conditions, control_group)

  if (length(treatment_groups) == 0) {
    stop("No treatment groups found. Only found control group: ", control_group)
  }

  result_list <- list()

  # 对每个处理组与NC进行比较
  for (treatment in treatment_groups) {
    control_data <- cluster_data$tail_length[cluster_data$condition == control_group]
    treatment_data <- cluster_data$tail_length[cluster_data$condition == treatment]

    # 确保数据存在
    if (length(control_data) == 0) {
      warning("No data found for control group '", control_group, "' in this cluster")
      next
    }
    if (length(treatment_data) == 0) {
      warning("No data found for treatment group '", treatment, "' in this cluster")
      next
    }

    # 执行t检验
    t_result <- t.test(control_data, treatment_data)

    n1 <- length(control_data)
    n2 <- length(treatment_data)
    pooled_sd <- sqrt(((n1-1)*var(control_data) + (n2-1)*var(treatment_data)) / (n1+n2-2))
    cohens_d <- (mean(control_data) - mean(treatment_data)) / pooled_sd

    if (logscale) {
      fold_change <- 2^(mean(treatment_data) - mean(control_data))
      log2_fc <- mean(treatment_data) - mean(control_data)
    } else {
      fold_change <- mean(treatment_data) / mean(control_data)
      log2_fc <- log2(fold_change)
    }

    result_list[[treatment]] <- list(
      statistic = t_result$statistic,
      p_value = t_result$p.value,
      estimate = fold_change,
      log2_fc = log2_fc,
      conf_int = t_result$conf.int,
      cohens_d = cohens_d,
      n_control = n1,
      n_treatment = n2,
      method = paste("Student's t-test: ", control_group, " vs ", treatment)
    )
  }

  return(result_list)
}

#' Perform Wilcoxon test for multiple treatments vs NC
#' @name perform_wilcoxon_all_mRNA_multi
#' @param cluster_data Data for a single cluster
#' @param logscale Whether to log2 transform tail lengths
#' @param control_group Name of control group (default: "NC")
perform_wilcoxon_all_mRNA_multi <- function(cluster_data, logscale = TRUE, control_group = "NC") {

  if (logscale) {
    cluster_data$tail_length <- log2(cluster_data$tail_length)
  }

  # 获取所有条件
  all_conditions <- unique(cluster_data$condition)

  # 检查控制组是否存在
  if (!control_group %in% all_conditions) {
    stop("Control group '", control_group, "' not found in data.")
  }

  # 获取处理组
  treatment_groups <- setdiff(all_conditions, control_group)

  if (length(treatment_groups) == 0) {
    stop("No treatment groups found.")
  }

  result_list <- list()

  # 对每个处理组与NC进行比较
  for (treatment in treatment_groups) {
    control_data <- cluster_data$tail_length[cluster_data$condition == control_group]
    treatment_data <- cluster_data$tail_length[cluster_data$condition == treatment]

    if (length(control_data) == 0 || length(treatment_data) == 0) {
      next
    }

    wilcox_result <- wilcox.test(control_data, treatment_data, exact = FALSE)

    n_total <- length(control_data) + length(treatment_data)
    z_value <- qnorm(wilcox_result$p.value / 2)
    effect_size_r <- abs(z_value) / sqrt(n_total)

    median_diff <- median(treatment_data) - median(control_data)
    if (logscale) {
      fold_change <- 2^median_diff
      log2_fc <- median_diff
    } else {
      fold_change <- median(treatment_data) / median(control_data)
      log2_fc <- log2(fold_change)
    }

    result_list[[treatment]] <- list(
      statistic = wilcox_result$statistic,
      p_value = wilcox_result$p.value,
      estimate = fold_change,
      log2_fc = log2_fc,
      effect_size_r = effect_size_r,
      n_control = length(control_data),
      n_treatment = length(treatment_data),
      method = paste("Wilcoxon Rank Sum Test: ", control_group, " vs ", treatment)
    )
  }

  return(result_list)
}

#' Perform LMM for multiple treatments vs NC
#' @name perform_lmm_multi
#' @param cluster_data Data for a single cluster
#' @param logscale Whether to log2 transform tail lengths
#' @param control_group Name of control group (default: "NC")
#' @description LMM implementation for multiple treatment groups
perform_lmm_multi <- function(cluster_data, logscale = TRUE, control_group = "NC") {

  if (logscale) {
    cluster_data$tail_length <- log2(cluster_data$tail_length)
  }

  # 获取所有条件
  all_conditions <- unique(cluster_data$condition)

  # 检查控制组是否存在
  if (!control_group %in% all_conditions) {
    stop("Control group '", control_group, "' not found in data.")
  }

  # 获取处理组
  treatment_groups <- setdiff(all_conditions, control_group)

  if (length(treatment_groups) == 0) {
    stop("No treatment groups found.")
  }

  result_list <- list()

  # 对每个处理组分别运行LMM
  for (treatment in treatment_groups) {
    tryCatch({
      # 只取当前处理组和控制组的数据
      subset_data <- cluster_data[cluster_data$condition %in% c(control_group, treatment), ]

      # 设置因子水平，以控制组为参考
      subset_data$condition <- factor(subset_data$condition,
                                      levels = c(control_group, treatment))

      # 构建模型公式
      if ("lib_id" %in% colnames(subset_data) &&
          length(unique(subset_data$lib_id)) > 1) {
        model_formula <- as.formula("tail_length ~ condition + (1 | lib_id)")
      } else {
        model_formula <- as.formula("tail_length ~ condition")
      }

      # 拟合模型
      if ("lib_id" %in% colnames(subset_data) &&
          length(unique(subset_data$lib_id)) > 1) {
        res <- lme4::lmer(model_formula, data = subset_data)
        coefficients <- summary(res)$coefficients
      } else {
        res <- lm(model_formula, data = subset_data)
        coefficients <- summary(res)$coefficients
      }

      # 检查是否有条件效应
      if (nrow(coefficients) < 2) {
        result_list[[treatment]] <- list(
          estimate = NA,
          std_error = NA,
          t_value = NA,
          p_value = NA,
          method = paste("LMM: ", control_group, " vs ", treatment, " - NO_CONDITION_EFFECT")
        )
        next
      }

      # 提取条件效应的估计值
      if (nrow(coefficients) > 1) {
        estimate <- coefficients[2, 1]  # Estimate for treatment
        std_error <- coefficients[2, 2]
        t_value <- coefficients[2, 3]

        # 计算p值
        if ("lib_id" %in% colnames(subset_data) &&
            length(unique(subset_data$lib_id)) > 1) {
          # 对于混合模型，使用lmerTest获取p值
          res_test <- lmerTest::lmer(model_formula, data = subset_data)
          p_value <- summary(res_test)$coefficients[2, 5]
        } else {
          # 对于线性模型
          df <- nrow(subset_data) - 2
          p_value <- 2 * pt(abs(t_value), df = df, lower.tail = FALSE)
        }

        result_list[[treatment]] <- list(
          estimate = estimate,
          std_error = std_error,
          t_value = t_value,
          p_value = p_value,
          n_control = sum(subset_data$condition == control_group),
          n_treatment = sum(subset_data$condition == treatment),
          method = paste("Linear Mixed Model: ", control_group, " vs ", treatment)
        )
      }

    }, error = function(e) {
      result_list[[treatment]] <- list(
        estimate = NA,
        std_error = NA,
        t_value = NA,
        p_value = NA,
        method = paste("LMM - ERROR for ", control_group, " vs ", treatment, ": ", e$message)
      )
    })
  }

  return(result_list)
}

##############################################
#         Main Analytic Functions
##############################################

#' Analyze PAS clusters using all mRNA tail lengths as independent observations
#' @name analyze_pas_clusters_mRNA_level_multi
#' @description Statistical analysis for multiple treatment groups vs control
#' @param QpolyA A QuantifyPolyA object
#' @param sample_info Data frame with sample metadata
#' @param test_methods Statistical methods to use: "t_test", "wilcoxon", "lmm"
#' @param min_mRNA_per_condition Minimum mRNA molecules per condition (default: 10)
#' @param logscale Whether to log2 transform tail lengths
#' @param control_group Name of control group (default: "NC")
#' @param mc.cores Number of cores for parallel processing
#' @return Data frame with statistical results for each PAS cluster and each treatment
#' @export
#'
polyAlength <- function(QpolyA, sample_info,
                        test_methods = c("t_test", "wilcoxon", "lmm"),
                        min_mRNA_per_condition = 10,
                        logscale = TRUE,
                        control_group = "NC",
                        mc.cores = 4) {

  if (class(QpolyA) != "QuantifyPolyA") {
    stop("QpolyA should be a QuantifyPolyA object!")
  }

  if (length(QpolyA@cluster_tail_lengths) == 0) {
    stop("No cluster tail length information found. Run mapTailLengthsToClusters first.")
  }

  required_cols <- c("sample", "condition")
  if (!all(required_cols %in% colnames(sample_info))) {
    stop("sample_info must contain columns: ", paste(required_cols, collapse = ", "))
  }

  if (!control_group %in% sample_info$condition) {
    warning("Control group '", control_group, "' not found in sample_info conditions. Found: ",
            paste(unique(sample_info$condition), collapse = ", "))
  }

  analysis_data <- prepare_cluster_tail_data_simple(QpolyA, sample_info, cores = mc.cores)

  if (nrow(analysis_data) == 0) {
    warning("No data available for statistical analysis")
    return(data.frame())
  }

  cluster_ids <- unique(analysis_data$cluster_id)

  message("Performing mRNA-level statistical tests on ", length(cluster_ids), " PAS clusters...")
  message("Total mRNA observations: ", nrow(analysis_data))
  message("Average mRNA per cluster: ", round(nrow(analysis_data) / length(cluster_ids), 1))
  message("Control group: ", control_group)

  all_treatments <- setdiff(unique(analysis_data$condition), control_group)
  message("Treatment groups: ", paste(all_treatments, collapse = ", "))

  results <- pbmcapply::pbmclapply(cluster_ids, function(cluster_id) {

    cluster_data <- analysis_data[analysis_data$cluster_id == cluster_id, ]

    if (!control_group %in% cluster_data$condition) {
      return(NULL)
    }

    treatment_groups <- setdiff(unique(cluster_data$condition), control_group)

    if (length(treatment_groups) == 0) {
      return(NULL)
    }

    mRNA_counts <- table(cluster_data$condition)

    if (mRNA_counts[control_group] < min_mRNA_per_condition) {
      return(NULL)
    }

    valid_treatments <- treatment_groups[treatment_groups %in% names(mRNA_counts) &
                                           mRNA_counts[treatment_groups] >= min_mRNA_per_condition]

    if (length(valid_treatments) == 0) {
      return(NULL)
    }

    result_row <- data.frame(
      cluster_id = cluster_id,
      n_samples = length(unique(cluster_data$sample)),
      n_total_mRNA = nrow(cluster_data),
      stringsAsFactors = FALSE
    )

    control_data <- cluster_data$tail_length[cluster_data$condition == control_group]
    result_row$n_control <- as.numeric(mRNA_counts[control_group])
    result_row$mean_control <- mean(control_data, na.rm = TRUE)
    result_row$median_control <- median(control_data, na.rm = TRUE)
    result_row$sd_control <- sd(control_data, na.rm = TRUE)

    for (treatment in all_treatments) {

      if (treatment %in% valid_treatments) {


        treatment_data <- cluster_data$tail_length[cluster_data$condition == treatment]

        result_row[[paste0(treatment, "_n")]] <- as.numeric(mRNA_counts[treatment])
        result_row[[paste0(treatment, "_mean")]] <- mean(treatment_data, na.rm = TRUE)
        result_row[[paste0(treatment, "_median")]] <- median(treatment_data, na.rm = TRUE)
        result_row[[paste0(treatment, "_sd")]] <- sd(treatment_data, na.rm = TRUE)

        if (logscale) {
          log2_fc <- mean(treatment_data, na.rm = TRUE) - mean(control_data, na.rm = TRUE)
          fc <- 2^log2_fc
        } else {
          fc <- mean(treatment_data, na.rm = TRUE) / mean(control_data, na.rm = TRUE)
          log2_fc <- log2(fc)
        }

        result_row[[paste0(treatment, "_fold_change")]] <- fc
        result_row[[paste0(treatment, "_log2_fc")]] <- log2_fc

        median_diff <- median(treatment_data, na.rm = TRUE) - median(control_data, na.rm = TRUE)
        result_row[[paste0(treatment, "_median_diff")]] <- median_diff

        for (test_method in test_methods) {
          tryCatch({
            if (test_method == "t_test") {
              t_results <- perform_t_test_all_mRNA_multi(cluster_data, logscale, control_group)
              if (treatment %in% names(t_results)) {
                t_result <- t_results[[treatment]]
                result_row[[paste0(treatment, "_t_statistic")]] <- t_result$statistic
                result_row[[paste0(treatment, "_t_p_value")]] <- t_result$p_value
                result_row[[paste0(treatment, "_t_cohens_d")]] <- t_result$cohens_d
              }
            } else if (test_method == "wilcoxon") {
              w_results <- perform_wilcoxon_all_mRNA_multi(cluster_data, logscale, control_group)
              if (treatment %in% names(w_results)) {
                w_result <- w_results[[treatment]]
                result_row[[paste0(treatment, "_w_statistic")]] <- w_result$statistic
                result_row[[paste0(treatment, "_w_p_value")]] <- w_result$p_value
                result_row[[paste0(treatment, "_w_effect_size")]] <- w_result$effect_size_r
              }
            } else if (test_method == "lmm") {
              lmm_results <- perform_lmm_multi(cluster_data, logscale, control_group)
              if (treatment %in% names(lmm_results)) {
                lmm_result <- lmm_results[[treatment]]
                result_row[[paste0(treatment, "_lmm_estimate")]] <- lmm_result$estimate
                result_row[[paste0(treatment, "_lmm_std_error")]] <- lmm_result$std_error
                result_row[[paste0(treatment, "_lmm_t_value")]] <- lmm_result$t_value
                result_row[[paste0(treatment, "_lmm_p_value")]] <- lmm_result$p_value
              }
            }
          }, error = function(e) {
            NULL
          })
        }

      } else {

        result_row[[paste0(treatment, "_n")]] <- NA_real_
        result_row[[paste0(treatment, "_mean")]] <- NA_real_
        result_row[[paste0(treatment, "_median")]] <- NA_real_
        result_row[[paste0(treatment, "_sd")]] <- NA_real_
        result_row[[paste0(treatment, "_fold_change")]] <- NA_real_
        result_row[[paste0(treatment, "_log2_fc")]] <- NA_real_
        result_row[[paste0(treatment, "_median_diff")]] <- NA_real_  # 新增


        if ("t_test" %in% test_methods) {
          result_row[[paste0(treatment, "_t_statistic")]] <- NA_real_
          result_row[[paste0(treatment, "_t_p_value")]] <- NA_real_
          result_row[[paste0(treatment, "_t_cohens_d")]] <- NA_real_
        }
        if ("wilcoxon" %in% test_methods) {
          result_row[[paste0(treatment, "_w_statistic")]] <- NA_real_
          result_row[[paste0(treatment, "_w_p_value")]] <- NA_real_
          result_row[[paste0(treatment, "_w_effect_size")]] <- NA_real_
        }
        if ("lmm" %in% test_methods) {
          result_row[[paste0(treatment, "_lmm_estimate")]] <- NA_real_
          result_row[[paste0(treatment, "_lmm_std_error")]] <- NA_real_
          result_row[[paste0(treatment, "_lmm_t_value")]] <- NA_real_
          result_row[[paste0(treatment, "_lmm_p_value")]] <- NA_real_
        }
      }
    }

    return(result_row)

  }, mc.cores = mc.cores)


  results <- results[!sapply(results, is.null)]

  if (length(results) == 0) {
    warning("No valid results obtained from statistical tests")
    return(data.frame())
  }


  final_results <- dplyr::bind_rows(results)


  rownames(final_results) <- NULL


  for (treatment in all_treatments) {
    t_p_col <- paste0(treatment, "_t_p_value")
    if (t_p_col %in% colnames(final_results)) {
      q_col <- paste0(treatment, "_t_q_value")
      final_results[[q_col]] <- p.adjust(final_results[[t_p_col]], method = "BH")
    }

    w_p_col <- paste0(treatment, "_w_p_value")
    if (w_p_col %in% colnames(final_results)) {
      q_col <- paste0(treatment, "_w_q_value")
      final_results[[q_col]] <- p.adjust(final_results[[w_p_col]], method = "BH")
    }

    lmm_p_col <- paste0(treatment, "_lmm_p_value")
    if (lmm_p_col %in% colnames(final_results)) {
      q_col <- paste0(treatment, "_lmm_q_value")
      final_results[[q_col]] <- p.adjust(final_results[[lmm_p_col]], method = "BH")
    }
  }

  final_results <- final_results[order(final_results$cluster_id), ]

  message("mRNA-level analysis completed. Results for ", nrow(final_results),
          " clusters in wide format.")

  return(final_results)
}
##############################################
#           Results summary function    #
##############################################

#' Generate summary of mRNA-level statistical results for multiple treatments
#' @name summarize_mRNA_level_results_multi
summarize_mRNA_level_results_multi <- function(results, alpha = 0.05, control_group = "NC") {

  if (nrow(results) == 0) {
    return("No results to summarize")
  }


  mean_cols <- grep("_mean$", colnames(results), value = TRUE)
  treatments <- unique(gsub("_mean$", "", mean_cols))

  treatments <- setdiff(treatments, control_group)

  if (length(treatments) == 0) {
    message("No treatment groups found in results")
    return(invisible(NULL))
  }

  summary_list <- list()

  for (treatment in treatments) {
    treatment_summary <- data.frame(
      treatment = treatment,
      n_clusters = nrow(results),
      stringsAsFactors = FALSE
    )

    t_p_col <- paste0(treatment, "_t_p_value")
    if (t_p_col %in% colnames(results)) {
      n_sig_t <- sum(results[[t_p_col]] < alpha, na.rm = TRUE)
      t_q_col <- paste0(treatment, "_t_q_value")
      n_sig_t_q <- if (t_q_col %in% colnames(results))
        sum(results[[t_q_col]] < alpha, na.rm = TRUE) else NA

      treatment_summary$t_test_total <- sum(!is.na(results[[t_p_col]]))
      treatment_summary$t_test_sig_p <- n_sig_t
      treatment_summary$t_test_sig_q <- n_sig_t_q


      t_cohens_col <- paste0(treatment, "_t_cohens_d")
      if (t_cohens_col %in% colnames(results)) {
        treatment_summary$mean_cohens_d <- mean(results[[t_cohens_col]], na.rm = TRUE)
      }
    }


    w_p_col <- paste0(treatment, "_w_p_value")
    if (w_p_col %in% colnames(results)) {
      n_sig_w <- sum(results[[w_p_col]] < alpha, na.rm = TRUE)
      w_q_col <- paste0(treatment, "_w_q_value")
      n_sig_w_q <- if (w_q_col %in% colnames(results))
        sum(results[[w_q_col]] < alpha, na.rm = TRUE) else NA

      treatment_summary$wilcoxon_total <- sum(!is.na(results[[w_p_col]]))
      treatment_summary$wilcoxon_sig_p <- n_sig_w
      treatment_summary$wilcoxon_sig_q <- n_sig_w_q
    }


    lmm_p_col <- paste0(treatment, "_lmm_p_value")
    if (lmm_p_col %in% colnames(results)) {
      n_sig_lmm <- sum(results[[lmm_p_col]] < alpha, na.rm = TRUE)
      lmm_q_col <- paste0(treatment, "_lmm_q_value")
      n_sig_lmm_q <- if (lmm_q_col %in% colnames(results))
        sum(results[[lmm_q_col]] < alpha, na.rm = TRUE) else NA

      treatment_summary$lmm_total <- sum(!is.na(results[[lmm_p_col]]))
      treatment_summary$lmm_sig_p <- n_sig_lmm
      treatment_summary$lmm_sig_q <- n_sig_lmm_q
    }


    log2_fc_col <- paste0(treatment, "_log2_fc")
    if (log2_fc_col %in% colnames(results)) {
      treatment_summary$mean_log2_fc <- mean(results[[log2_fc_col]], na.rm = TRUE)
      treatment_summary$median_log2_fc <- median(results[[log2_fc_col]], na.rm = TRUE)
      treatment_summary$sd_log2_fc <- sd(results[[log2_fc_col]], na.rm = TRUE)

      up_regulated <- sum(results[[log2_fc_col]] > 0, na.rm = TRUE)
      down_regulated <- sum(results[[log2_fc_col]] < 0, na.rm = TRUE)
      treatment_summary$up_regulated <- up_regulated
      treatment_summary$down_regulated <- down_regulated
      treatment_summary$up_percent <- round(up_regulated / nrow(results) * 100, 1)
      treatment_summary$down_percent <- round(down_regulated / nrow(results) * 100, 1)
    }
    median_diff_col <- paste0(treatment, "_median_diff")
    if (median_diff_col %in% colnames(results)) {
      treatment_summary$mean_median_diff <- mean(results[[median_diff_col]], na.rm = TRUE)
      treatment_summary$median_median_diff <- median(results[[median_diff_col]], na.rm = TRUE)
      treatment_summary$sd_median_diff <- sd(results[[median_diff_col]], na.rm = TRUE)
    }

    summary_list[[treatment]] <- treatment_summary
  }

  summary_df <- dplyr::bind_rows(summary_list)

  cat("=== mRNA-LEVEL STATISTICAL TEST SUMMARY (", control_group, " vs treatments) ===\n", sep = "")
  print(summary_df)

  return(invisible(summary_df))
}

#' APPLE: Analysis of Poly(A) Lengths and Expression
#'
#' @description
#' The 'APPLE' package provides a comprehensive workflow for analyzing poly(A) tail
#' lengths and alternative polyadenylation (APA) from sequencing data. It supports
#' the entire analysis pipeline: from raw SAM/BED file processing, poly(A) site
#' clustering, tail length extraction, to statistical testing for differential
#' polyadenylation between multiple experimental conditions.
#'
#' @details
#' The main features of APPLE include:
#' \itemize{
#'   \item **Data import and preprocessing**: Functions to read SAM files, detect
#'         poly(A) tails (using pt tags or regex), and generate BED files with
#'         tail length information (\code{\link{Extract_polyAsite}}, \code{\link{Load.PolyA}}).
#'   \item **Poly(A) site clustering**: Weighted density peak clustering to define
#'         Poly(A) Clusters (PACs) (\code{\link{Cluster.PolyA}}).
#'   \item **Tail length mapping**: Map individual mRNA tail lengths to PACs
#'         (\code{\link{mapTail}}).
#'   \item **Annotation and filtering**: Annotate PACs with genomic features
#'         (\code{\link{Annotate.PolyA}}) and filter low-count PACs (\code{\link{Filter.PolyA}}).
#'   \item **Statistical analysis**: Perform differential polyadenylation analysis
#'         using various approaches:
#'         \itemize{
#'           \item mRNA-level t-test, Wilcoxon, or linear mixed models for multiple
#'                 treatment groups vs control (\code{\link{polyAlength}}).
#'           \item Quantify APA dynamics among split, canonical, or whole-gene PACs
#'                 (\code{\link{Quantify.SplitAPA}}, \code{\link{Quantify.CanonicalAPA}},
#'                 \code{\link{Quantify.CNCAPA}}, \code{\link{Quantify.GeneAPA}}).
#'           \item Principal Component Analysis on tail length matrices
#'                 (\code{\link{tail_pca}}).
#'         }
#'   \item **Visualization**: PCA plots, UMAP, and DESeq2-based sample quality
#'         assessment (\code{\link{DESeq2.PolyA}}).
#'   \item **Utility functions**: Save tail length tables (\code{\link{Save.TailLengths}}),
#'         search for poly(A) signals (\code{\link{Motif.Search}}), and generate
#'         analysis summaries (\code{\link{summarize_mRNA_level_results_multi}}).
#' }
#'
#' The package defines an S4 class \code{\linkS4class{QuantifyPolyA}} that holds
#' all raw and processed data throughout the workflow.
#'
#' @section Package options:
#' No specific options are currently defined.
#'
#' @section Dependencies:
#' APPLE relies on several CRAN and Bioconductor packages, including:
#' \code{lme4}, \code{lmerTest}, \code{ggplot2}, \code{dplyr}, \code{tidyr},
#' \code{matrixStats}, \code{data.table}, \code{parallel}, \code{pbmcapply},
#' \code{GenomicRanges}, \code{Rsamtools}, \code{rtracklayer}, \code{plyranges},
#' \code{DESeq2}, \code{ggbio}, \code{FactoMineR}, \code{factoextra}, \code{uwot},
#' and others. Most of these are automatically installed when installing the package
#' via \code{BiocManager::install("APPLE")} (once submitted to Bioconductor) or
#' manually from CRAN/Bioconductor.
#'
#' @note
#' This package is designed for users familiar with poly(A) tail sequencing assays
#' (e.g., PAT-seq, PAL-seq, FLAM-seq) and R/Bioconductor.
#'
#' @author
  #' * Contributor 1
  #' * Contributor 2
#'
#' @docType package
#' @name APPLE-package
#' @aliases APPLE
#' @importFrom methods new
#' @importFrom stats as.formula complete.cases lm median p.adjust prcomp pt qnorm sd t.test var wilcox.test
#' @importFrom utils install.packages write.table
#' @importFrom dplyr %>% bind_rows group_by summarise filter n mutate full_join left_join rename case_when
#' @importFrom tidyr pivot_wider
#' @importFrom ggplot2 ggplot aes geom_point xlab ylab coord_fixed
#' @importFrom GenomicRanges GRanges findOverlaps
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors queryHits subjectHits
#' @importFrom BiocGenerics start end width
#' @importFrom tools file_path_sans_ext
#' @importFrom stringr str_split str_extract str_replace str_sub str_length str_detect
#' @importFrom parallel mclapply
#' @importFrom pbmcapply pbmclapply
#' @importFrom matrixStats rowMins rowMaxs
#' @importFrom DESeq2 DESeqDataSetFromMatrix DESeq counts vst
#' @importFrom factoextra fviz_pca_ind
#' @importFrom uwot umap
#' @importFrom FactoMineR PCA
#' @importFrom SummarizedExperiment assay
#' @importFrom S4Vectors isEmpty
#' @importFrom BiocGenerics which.max
#' @importFrom outliers scores
#' @importFrom plyranges join_overlap_intersect_directed as_granges find_overlaps
#' @importFrom rtracklayer import
#' @importFrom Rsamtools indexFa
#' @importFrom bedr bedr
"_PACKAGE"
