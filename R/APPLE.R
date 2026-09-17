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


  # 记录当前工作目录，函数退出时恢复
  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)

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

    # Chromosome and strand identifiers are labels, even when a BED subset
    # contains only numeric chromosome names (for example, chromosome 22).
    data[[1]] = as.character(data[[1]])
    data[[2]] = as.character(data[[2]])

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
#       Remove internal priming events       #
##############################################
#' @title Remove IP
#' @description Remove poly(A) sites which potentially result from internal priming events.
#' @name Remove.IP
#' @param QpolyA A QuantifyPolyA object containing all raw poly(A) sites.
#' @param fasta A string specifying the location and name of the corresponding genome file in FASTA format.
#' @param flank_len An integer specifying the neighborhood to search for template polyA.
#' @param win_size An integer specifying the size of the sliding window.
#' @param min_A An integer specifying the minimum number of base A in a window.
#' @return A QuantifyPolyA object containing poly(A) sites with IP removed.
#' @export
#'
Remove.IP <- function(QpolyA,fasta,flank_len=15,win_size=10,min_A=8){
  # Check parameters.
  if (!is(QpolyA, "QuantifyPolyA")) stop(paste('QpolyA should be a QuantifyPolyA object!'))
  if (!file.exists(fasta)) stop(paste('Fasta file',fasta,'does not exist!'))
  if (!is.numeric(flank_len)) stop("'flank_len' is not a number!")
  if (flank_len<=0) stop("'flank_len' should be larger than 0, the default value is 15!")
  if (!is.numeric(win_size)) stop("'win_size' is not a number!")
  if (win_size<=0) stop("'win_size' should be larger than 0, the default value is 10!")
  if (!is.numeric(min_A)) stop("'min_A' is not a number!")
  if (min_A<=0) stop("'min_A' should be larger than 0, the default value is 8!")
  
  # Remove internal priming events.
  for (alt_name in QpolyA@sample_names) {
    is.IP = is.internal.priming(QpolyA@pre.polyA[[alt_name]],fasta = fasta,flank_len,win_size,min_A)
    print(paste(sum(is.IP),'internal priming events were found in sample', alt_name,'!'))
    QpolyA@pre.polyA[[alt_name]] = QpolyA@pre.polyA[[alt_name]][!is.IP,]

    if (alt_name %in% names(QpolyA@tail_lengths)) {
      QpolyA@tail_lengths[[alt_name]] = QpolyA@tail_lengths[[alt_name]][!is.IP,]
    }
  }
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
  if (!is(QpolyA, "QuantifyPolyA")) stop(paste('QpolyA should be a QuantifyPolyA object!'))
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
#' @name Map.Tail
#' @param QpolyA A QuantifyPolyA object containing clean poly(A) sites and PACs.
#' @param delimiter The delimiter to use for tail lengths (default: ";")
#' @return A QuantifyPolyA object with tail length information mapped to clusters.
#' @importFrom GenomicRanges GRanges findOverlaps
#' @importFrom IRanges IRanges
#' @importFrom dplyr group_by summarise
#' @export
#'
Map.Tail <- function(QpolyA, delimiter = ";") {
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
  if (!is(QpolyA, "QuantifyPolyA")) stop(paste('QpolyA should be a QuantifyPolyA object!'))
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
  if (!is(QpolyA, "QuantifyPolyA")) stop(paste('QpolyA should be a QuantifyPolyA object!'))
  if (!is.numeric(min_count)) stop("'min_count' is not a number!")
  if (min_count<=0) stop("'min_count' should be larger than 0, the default value is 10!")
  if (!is.numeric(min_sample)) stop("'min_sample' is not a number!")
  if (min_sample<=0) stop("'min_sample' should be larger than 0, the default value is 1!")

  QpolyA@polyA = subset(QpolyA@polyA, rowSums(QpolyA@polyA[,QpolyA@sample_names]>=min_count) >= min_sample)
  rownames(QpolyA@polyA) = paste('PA',1:nrow(QpolyA@polyA),sep = '')
  return(QpolyA)
}


##############################################
#       Check internal priming events        #
##############################################
#' Check whether the poly(A) sites are resulted from internal priming artifacts or not.
#' @name is.internal.priming
#' @param polyA A data.frame containing the information of raw poly(A) sites.
#' @param fasta A string specifies the location and name of genome file in FASTA format.
#' @param flank_len An integer specifying the neighborhood to search for template polyA.
#' @param win_size An integer specifying the size of the sliding window.
#' @param min_A An integer specifying the minimum number of base A in a window.
#' @return A logic vector indicates whether the poly(A) sites are resulted from internal priming artifacts or not.
#' @importFrom dplyr left_join
#' @importFrom bedr get.fasta
#' @importFrom Rsamtools indexFa
#' @importFrom stringr str_detect str_split
#' @importFrom matrixStats rowCumsums
#' @export
#'
is.internal.priming <- function(polyA,fasta,flank_len=15,win_size=10,min_A=8){
  # Construct a data.frame to build a bed file for extracting sequences
  seq.bed = data.frame(
    chr = polyA$seqnames,
    start = polyA$coord-flank_len,
    end = polyA$coord+flank_len,
    strand = polyA$strand,
    stringsAsFactors = FALSE
  )
  
  # Check index file of fasta file
  if (!file.exists(paste0(fasta,'.fai'))) indexFa(fasta)
  
  # Get info of fasta file
  fai = read.table(file = paste0(fasta,'.fai'))
  colnames(fai) = c('chr','len','offset','linebase','linewidth')

  # Keep join keys type-stable for BED subsets containing numeric-only
  # chromosome names and FASTA indexes that also contain X, Y, or MT.
  seq.bed$chr = as.character(seq.bed$chr)
  fai$chr = as.character(fai$chr)
  
  # Check for valid regions i.e. start < 0 or end > seq.len
  seq.bed = left_join(seq.bed,fai[,c(1,2)], by='chr')
  seq.bed$start[seq.bed$start<=0] = 1
  idx = seq.bed$end>seq.bed$len
  seq.bed$end[idx] = seq.bed$len[idx]
  
  # Sort
  idx = order(seq.bed$chr,seq.bed$start,decreasing = FALSE)
  seq.bed = seq.bed[idx,]
  
  # Extract sequences
  seq.df = get.fasta(seq.bed,fasta = fasta,verbose = FALSE,check.chr = FALSE,check.sort = FALSE,check.valid = FALSE)
  
  # Check internal priming events
  is.IP = rep(TRUE,nrow(seq.bed))
  plus.idx = which(seq.bed$strand=='+')
  minus.idx = which(seq.bed$strand=='-')
  
  is.IP[plus.idx] = str_detect(seq.df$sequence[plus.idx],'A{6,}')
  is.IP[minus.idx] = str_detect(seq.df$sequence[minus.idx],'T{6,}')
  
  plus.chars = str_split(seq.df$sequence[plus.idx],'',simplify = T)=='A'
  plus.chars = rowCumsums(plus.chars,na.rm=TRUE)
  plus.chars = rowSums(cbind(plus.chars[,win_size],plus.chars[,(win_size+1):(2*flank_len)] - plus.chars[,1:(2*flank_len-win_size)])>=min_A,na.rm = TRUE)
  
  is.IP[plus.idx] = is.IP[plus.idx]|plus.chars
  
  minus.chars = str_split(seq.df$sequence[minus.idx],'',simplify = T)=='T'
  minus.chars = rowCumsums(minus.chars,na.rm=TRUE)
  minus.chars = rowSums(cbind(minus.chars[,win_size],minus.chars[,(win_size+1):(2*flank_len)] - minus.chars[,1:(2*flank_len-win_size)])>=min_A,na.rm = TRUE)
  
  is.IP[minus.idx] = is.IP[minus.idx]|minus.chars
  
  # return result
  is.IP[idx] = is.IP
  return(is.IP)
}


##############################################
#          construct genome ranges           #
##############################################
#' Construct genomic ranges
#' @name buildGenomicRanges
#' @param seqname A vector containing the sequence (Chromosomes/Contigs) names of poly(A) sites.
#' @param position A numeric vector containing the genomic positions of poly(A) sites.
#' @param score A numeric vector containing the numbers of reads supporting each poly(A) site.
#' @param strand A character vector containing the strand information of poly(A) sites.
#' @param five_prime_end A numeric vector of read 5-prime coordinates.
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
  range.gr = GenomicRanges::reduce(points.gr,min.gapwidth=max.gapwidth,with.revmap=T,ignore.strand=FALSE)

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
#' @param sub_pos The genomic coordinates of poly(A) sites.
#' @param sub_wts The numbers of reads support each poly(A) site in 'sub_pos'.
#' @param sub_five_prime_end Read 5-prime coordinates corresponding to sub_pos.
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
    res <- tibble(group=1,start=min(sub_pos),end=max(sub_pos),sum.wts=sum(sub_wts),center=sub_pos[which.max(sub_wts)],five_prime_end = round(median(sub_five_prime_end)))
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
#' @importFrom txdbmaker makeTxDbFromGFF
#' @importFrom GenomicFeatures cds exons genes threeUTRsByTranscript fiveUTRsByTranscript
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
  txdb =makeTxDbFromGFF(gff)
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

  print('start refinement!')
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
  if (!is(QpolyA, "QuantifyPolyA")) stop(paste('QpolyA should be a QuantifyPolyA object!'))
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
  polyA.normailized = cbind(QpolyA@polyA[,c(1:11)],DESeq2::counts(dds,normalized=TRUE))

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
  if (!is(QpolyA, "QuantifyPolyA")) stop(paste('QpolyA should be a QuantifyPolyA object!'))
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
  if (!is(QpolyA, "QuantifyPolyA")) stop(paste('QpolyA should be a QuantifyPolyA object!'))
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
  if (!is(QpolyA, "QuantifyPolyA")) stop(paste('QpolyA should be a QuantifyPolyA object!'))
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
  if (!is(QpolyA, "QuantifyPolyA")) stop(paste('QpolyA should be a QuantifyPolyA object!'))
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
#' @importFrom rlang sym
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
#' delta <- compute_delta_RPP(polyA_rank, colData, "NC", "EX1")
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
#' @return A numeric value between -1 and 1.
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


#' Principal component analysis of poly(A) tail lengths
#'
#' Summarize positive tail lengths per cluster and sample, filter low-count,
#' incomplete and constant clusters, and perform centered, scaled PCA.
#' @param QpolyA A QuantifyPolyA object processed by Map.Tail().
#' @param sample_info A data frame containing sample and condition columns.
#'   An optional lib_id column identifies libraries.
#' @param summary_stat Summary statistic, either "mean" or "median".
#' @param min_count_per_sample Minimum positive tail observations per cluster
#'   and sample.
#' @param cores Number of cores. Use 1 on Windows.
#' @param show_progress Whether to show progress during sample processing.
#' @return A list with matrix (clusters by samples), sample_info, pca (a
#'   prcomp object), filtered_data, and removed_clusters.
#' @export
Tail.PCA <- function(QpolyA,
                     sample_info,
                     summary_stat = c("mean", "median"),
                     min_count_per_sample = 10,
                     cores = 4,
                     show_progress = TRUE) {

  if (!requireNamespace("parallel", quietly = TRUE))
    stop("parallel is required!")
  if (!requireNamespace("data.table", quietly = TRUE))
    stop("data.table is required!")

  if (is.null(QpolyA@cluster_tail_lengths))
    stop("QpolyA@cluster_tail_lengths is NULL!")

  summary_stat <- match.arg(summary_stat)

  sample_info_dt <- data.table::as.data.table(sample_info)

  # 校验
  if (!"sample" %in% colnames(sample_info_dt))
    stop("sample_info must contain a column named 'sample'")
  if (!"condition" %in% colnames(sample_info_dt))
    stop("sample_info must contain a column named 'condition'")

  missing_samples <- setdiff(sample_info_dt$sample, names(QpolyA@cluster_tail_lengths))
  if (length(missing_samples) > 0)
    stop("Samples in sample_info not found in QpolyA: ",
         paste(missing_samples, collapse = ", "))

  sample_names <- sample_info_dt$sample

  if (!"lib_id" %in% colnames(sample_info_dt))
    sample_info_dt[, lib_id := NA_character_]

  # 预提取数据，避免并行时传递整个 QpolyA
  sample_data_list <- lapply(sample_names, function(snm) {
    dt <- QpolyA@cluster_tail_lengths[[snm]]
    if (is.null(dt) || nrow(dt) == 0) return(NULL)
    meta <- sample_info_dt[sample == snm]
    if (nrow(meta) == 0) return(NULL)
    list(
      cluster_id = dt$cluster_id,
      all_tail_lengths = dt$all_tail_lengths,
      sample = snm,
      condition = meta$condition[1],
      lib_id = meta$lib_id[1]
    )
  })
  names(sample_data_list) <- sample_names
  sample_data_list <- sample_data_list[!sapply(sample_data_list, is.null)]

  if (length(sample_data_list) == 0)
    stop("No samples have data in cluster_tail_lengths.")

  # process_sample 直接接收单个样本的小数据块
  process_sample <- function(sample_data) {
    tail_list <- strsplit(sample_data$all_tail_lengths, ";", fixed = TRUE)
    cluster_ids <- sample_data$cluster_id

    res_list <- vector("list", length(cluster_ids))
    for (i in seq_along(cluster_ids)) {
      tails <- as.numeric(tail_list[[i]])
      tails <- tails[!is.na(tails) & tails > 0]
      n <- length(tails)
      if (n < min_count_per_sample) next

      val <- if (summary_stat == "mean") mean(tails) else median(tails)
      res_list[[i]] <- data.table::data.table(
        cluster_id = cluster_ids[i],
        sample = sample_data$sample,
        condition = sample_data$condition,
        lib_id = sample_data$lib_id[1],
        summary_value = val,
        count = n
      )
    }
    out <- data.table::rbindlist(res_list, use.names = TRUE, fill = FALSE)
    if (nrow(out) == 0) {
      out <- data.table::data.table(
        cluster_id = character(),
        sample = character(),
        condition = character(),
        lib_id = character(),
        summary_value = numeric(),
        count = integer()
      )
    }
    out
  }

  message("Using ", cores, " CPU cores")
  message("Processing ", length(sample_data_list), " samples")
  message("Per-sample minimum count threshold = ", min_count_per_sample)

  if (show_progress && requireNamespace("pbmcapply", quietly = TRUE)) {
    result_list <- pbmcapply::pbmclapply(sample_data_list, process_sample,
                                         mc.cores = cores,
                                         ignore.interactive = !interactive())
  } else {
    result_list <- parallel::mclapply(sample_data_list, process_sample,
                                      mc.cores = cores)
  }

  # 移除 NULL 和无效结果
  result_list <- result_list[!sapply(result_list, is.null)]
  is_valid <- sapply(result_list, function(x) inherits(x, "data.table") || is.data.frame(x))
  if (any(!is_valid)) {
    warning(sum(!is_valid), " sample(s) returned invalid data and were skipped.")
    result_list <- result_list[is_valid]
  }

  if (length(result_list) == 0)
    stop("No data after per-sample filtering.")

  all_dt <- data.table::rbindlist(result_list, use.names = TRUE, fill = TRUE)
  message("Total (cluster-sample) pairs after filtering: ", nrow(all_dt))

  if (nrow(all_dt) == 0)
    stop("All clusters filtered out. Lower min_count_per_sample?")

  # 构建宽矩阵 (clusters x samples)
  wide_mat <- data.table::dcast(all_dt,
                                cluster_id ~ sample,
                                value.var = "summary_value",
                                fun.aggregate = mean,
                                fill = NA)

  cluster_ids <- wide_mat$cluster_id
  mat <- as.matrix(wide_mat[, -1])
  rownames(mat) <- cluster_ids
  colnames(mat) <- names(wide_mat)[-1]

  # 移除有 NA 的 cluster（至少一个样本不满足 min_count 阈值）
  complete_clusters <- apply(mat, 1, function(row) !any(is.na(row)))
  mat_complete <- mat[complete_clusters, , drop = FALSE]
  removed_clusters <- sum(!complete_clusters)
  if (removed_clusters > 0) {
    message(removed_clusters, " clusters removed due to insufficient count in at least one sample.")
    message(sum(complete_clusters), " clusters kept (present in all samples).")
  } else {
    message("All clusters are present in all samples.")
  }

  if (nrow(mat_complete) == 0)
    stop("No clusters remain after removing incomplete ones. Lower min_count_per_sample?")

  # 移除零方差的 cluster
  col_var <- apply(mat_complete, 1, var, na.rm = TRUE)
  const_cols <- which(col_var == 0 | is.na(col_var))
  if (length(const_cols) > 0) {
    message(length(const_cols), " constant clusters (zero variance) removed before PCA.")
    mat_complete <- mat_complete[-const_cols, , drop = FALSE]
  }

  if (nrow(mat_complete) == 0)
    stop("No variable clusters remain. Cannot perform PCA.")

  # PCA (samples x clusters)
  mat_for_pca <- t(mat_complete)
  pca_res <- prcomp(mat_for_pca, center = TRUE, scale. = TRUE)

  # 样本元数据（按 pca 行顺序）
  sample_metadata <- unique(all_dt[, .(sample, condition, lib_id)])
  sample_metadata <- sample_metadata[sample %in% rownames(pca_res$x), ]

  list(
    matrix = mat_complete,
    sample_info = sample_metadata,
    pca = pca_res,
    filtered_data = all_dt,
    removed_clusters = names(which(!complete_clusters))
  )
}


#' Convert the raw tail length data into a long table
#' @name prepare_cluster_tail_data_fast
#' @description Fast parallel processing with only essential columns
#' @param QpolyA A QuantifyPolyA object
#' @param sample_info Data frame with sample metadata
#' @param cores Number of CPU cores to use
#' @param show_progress Whether to show processing progress.
#' @return A data table with individual mRNA tail lengths (cluster_id, sample, condition, tail_length, lib_id)
#'
# ====================== 优化的数据准备函数 ======================
prepare_cluster_tail_data_fast <- function(QpolyA, sample_info, cores = 4, show_progress = TRUE) {

  if (!requireNamespace("parallel", quietly = TRUE))
    stop("parallel is required for this function!")
  if (!requireNamespace("data.table", quietly = TRUE))
    stop("data.table is required for this function!")
  if (is.null(QpolyA@cluster_tail_lengths))
    stop("QpolyA@cluster_tail_lengths is NULL!")

  sample_info_dt <- data.table::as.data.table(sample_info)

  if (!"sample" %in% colnames(sample_info_dt))
    stop("sample_info must contain a column named 'sample'")
  if (!"condition" %in% colnames(sample_info_dt))
    stop("sample_info must contain a column named 'condition'")

  missing_samples <- setdiff(sample_info_dt$sample, names(QpolyA@cluster_tail_lengths))
  if (length(missing_samples) > 0)
    stop("Samples in sample_info not found in QpolyA: ",
         paste(missing_samples, collapse = ", "))

  sample_names <- sample_info_dt$sample

  if (!"lib_id" %in% colnames(sample_info_dt))
    sample_info_dt[, lib_id := NA_character_]

  # 预提取每个样本的数据，打包成列表元素
  sample_data_list <- lapply(sample_names, function(snm) {
    dt <- QpolyA@cluster_tail_lengths[[snm]]
    if (is.null(dt) || nrow(dt) == 0) return(NULL)
    meta <- sample_info_dt[sample == snm]
    if (nrow(meta) == 0) return(NULL)
    list(
      cluster_id = dt$cluster_id,
      all_tail_lengths = dt$all_tail_lengths,
      sample = snm,
      condition = meta$condition[1],
      lib_id = meta$lib_id[1]
    )
  })
  names(sample_data_list) <- sample_names
  sample_data_list <- sample_data_list[!sapply(sample_data_list, is.null)]

  if (length(sample_data_list) == 0)
    stop("No samples have data in cluster_tail_lengths.")

  # 核心优化：直接对 sample_data_list 的每个元素做 mclapply
  # 这样每个子进程只拿到自己需要的那一个元素，而不是整个列表
  process_sample <- function(sample_data) {
    tail_list <- strsplit(sample_data$all_tail_lengths, ";", fixed = TRUE)
    cluster_id_rep <- rep(sample_data$cluster_id, lengths(tail_list))
    tail_vec <- as.numeric(unlist(tail_list, use.names = FALSE))
    valid <- !is.na(tail_vec) & tail_vec > 0
    if (!any(valid)) return(NULL)

    data.table::data.table(
      cluster_id = cluster_id_rep[valid],
      sample = sample_data$sample,
      condition = sample_data$condition,
      tail_length = tail_vec[valid],
      lib_id = sample_data$lib_id
    )
  }

  message("Using ", cores, " CPU cores")
  message("Processing ", length(sample_data_list), " samples")

  if (show_progress && requireNamespace("pbmcapply", quietly = TRUE)) {
    result_list <- pbmcapply::pbmclapply(
      sample_data_list,   # <--- 直接传列表元素，不是 names
      process_sample,
      mc.cores = cores, ignore.interactive = !interactive()
    )
  } else {
    result_list <- parallel::mclapply(
      sample_data_list,   # <--- 直接传列表元素，不是 names
      process_sample,
      mc.cores = cores
    )
  }

  result_list <- result_list[!sapply(result_list, is.null)]
  if (length(result_list) == 0) return(data.table::data.table())

  final_dt <- data.table::rbindlist(result_list)
  message("Processing completed! Total rows: ", nrow(final_dt))
  final_dt[]
}


##############################################
#           Statistical test function    #
##############################################
#' Pairwise comparison between treatment and control for each PAS cluster
#' @name Tail.DiffPair
#' @description Statistical analysis for tail length
#' @param QpolyA QuantifyPolyA object
#' @param sample_info data.frame with columns: sample, condition, and optionally lib_id
#' @param control_group character, name of control condition
#' @param treatment_group character, name of treatment condition
#' @param test_method character, one of "t_test", "wilcoxon", "lmm"
#' @param logscale logical, if TRUE apply log2 transformation before testing
#' @param min_mRNA_per_condition integer, minimum number of mRNA molecules per group to include cluster
#' @param mc.cores integer, number of cores for parallel processing
#' @return data.frame with one row per cluster, containing descriptive statistics (raw scale) and test results
#' @export
Tail.DiffPair <- function(QpolyA, sample_info,
                          control_group, treatment_group,
                          test_method = c("t_test", "wilcoxon", "lmm"),
                          logscale = TRUE,
                          min_mRNA_per_condition = 10,
                          mc.cores = 4) {
  
  test_method <- match.arg(test_method)
  
  # 参数检查
  if (!is(QpolyA, "QuantifyPolyA"))
    stop("QpolyA must be a QuantifyPolyA object!")
  required_cols <- c("sample", "condition")
  if (!all(required_cols %in% colnames(sample_info)))
    stop("sample_info must contain columns: sample, condition!")
  if (!control_group %in% sample_info$condition)
    stop("control_group '", control_group, "' not found in sample_info$condition!")
  if (!treatment_group %in% sample_info$condition)
    stop("treatment_group '", treatment_group, "' not found in sample_info$condition!")
  if (test_method == "lmm") {
    if (!requireNamespace("lme4", quietly = TRUE))
      stop("Package 'lme4' is required for LMM. Please install it.")
    if (!"lib_id" %in% colnames(sample_info))
      message("Note: 'lib_id' column not found in sample_info. LMM will fall back to linear model.")
  }
  
  # 数据准备
  message("Preparing data...")
  dt_all <- prepare_cluster_tail_data_fast(QpolyA, sample_info, cores = mc.cores)
  if (nrow(dt_all) == 0) stop("No data available after preparation.")
  data.table::setDT(dt_all)
  
  # 只保留所需两组
  dt_all <- dt_all[condition %in% c(control_group, treatment_group)]
  if (nrow(dt_all) == 0) stop("No data for the specified groups.")
  
  cluster_ids <- unique(dt_all$cluster_id)
  message("Processing ", length(cluster_ids), " clusters (", control_group, " vs ", treatment_group, ")")
  
  # 并行处理每个 cluster
  dt_all <- split(dt_all, by = "cluster_id")

  results_list <- pbmcapply::pbmclapply(dt_all, function(cluster_dt) {
    cid <- cluster_dt$cluster_id[1] 
    
    # 样本量检查
    n_ctrl <- sum(cluster_dt$condition == control_group)
    n_trt  <- sum(cluster_dt$condition == treatment_group)
    if (n_ctrl < min_mRNA_per_condition || n_trt < min_mRNA_per_condition)
      return(NULL)
    
    # ---------- 描述统计（基于原始尾长）----------
    raw_vals <- split(cluster_dt$tail_length, cluster_dt$condition)
    control_raw <- raw_vals[[control_group]]
    treat_raw   <- raw_vals[[treatment_group]]
    
    mean_ctrl <- mean(control_raw);   mean_trt <- mean(treat_raw)
    median_ctrl <- median(control_raw); median_trt <- median(treat_raw)
    sd_ctrl <- sd(control_raw);       sd_trt <- sd(treat_raw)
    fold_change_raw <- mean_trt / mean_ctrl
    mean_diff_raw <- mean_trt - mean_ctrl
    median_diff_raw <- median_trt - median_ctrl
    
    # ---------- 检验数据（可能 log2 转换）----------
    if (logscale) {
      cluster_dt[, test_length := log2(tail_length)]
    } else {
      cluster_dt[, test_length := tail_length]
    }
    test_vals <- split(cluster_dt$test_length, cluster_dt$condition)
    control_test <- test_vals[[control_group]]
    treat_test   <- test_vals[[treatment_group]]
    
    # ---------- 执行检验 ----------
    if (test_method == "t_test") {
      t_res <- t.test(treat_test, control_test)
      n1 <- length(control_test); n2 <- length(treat_test)
      pooled_sd <- sqrt(((n1-1)*var(control_test) + (n2-1)*var(treat_test)) / (n1+n2-2))
      cohens_d <- (mean(treat_test) - mean(control_test)) / pooled_sd

      mean_diff_test <- mean(treat_test) - mean(control_test)
      log2_fc <- if(logscale) mean_diff_test else log2(fold_change_raw)
      
      res <- data.table::data.table(
        cluster_id = cid,
        n_control = n_ctrl, n_treatment = n_trt,
        mean_control = mean_ctrl, mean_treatment = mean_trt,
        median_control = median_ctrl, median_treatment = median_trt,
        sd_control = sd_ctrl, sd_treatment = sd_trt,
        fold_change = fold_change_raw,
        mean_diff = mean_diff_raw,
        median_diff = median_diff_raw,
        log2_fc = log2_fc,
        statistic = t_res$statistic,
        p_value = t_res$p.value,
        cohens_d = cohens_d,
        method = paste("t-test", if(logscale) "on log2 scale" else "on raw scale")
      )
      
    } else if (test_method == "wilcoxon") {
      w_res <- wilcox.test(treat_test, control_test, exact = FALSE)
      n1 <- length(control_test); n2 <- length(treat_test)
      W <- w_res$statistic

      # 期望和标准差（带 ties 校正）
      all_vals <- c(control_test, treat_test)
      ties_tab <- table(all_vals)
      ties_corr <- sum(ties_tab^3 - ties_tab) / (12 * (n1+n2) * (n1+n2-1))
      mu <- n1 * n2 / 2
      sigma <- sqrt(n1 * n2 * ((n1 + n2 + 1) - ties_corr) / 12)

      # 连续性修正
      z_val <- (W - mu - 0.5 * sign(W - mu)) / sigma

      # 带方向的效应量
      r <- z_val / sqrt(n1 + n2)

      median_diff_test <- median(treat_test) - median(control_test)
      log2_fc <- if(logscale) median_diff_test else log2(fold_change_raw)
      
      res <- data.table::data.table(
        cluster_id = cid,
        n_control = n_ctrl, n_treatment = n_trt,
        mean_control = mean_ctrl, mean_treatment = mean_trt,
        median_control = median_ctrl, median_treatment = median_trt,
        sd_control = sd_ctrl, sd_treatment = sd_trt,
        fold_change = fold_change_raw,
        mean_diff = mean_diff_raw,
        median_diff = median_diff_raw,
        log2_fc = log2_fc,
        statistic = w_res$statistic,
        p_value = w_res$p.value,
        effect_size_r = r,
        method = paste("Wilcoxon", if(logscale) "on log2 scale" else "on raw scale")
      )
      
    } else { # lmm
      # 确保 condition 为因子，对照组为参考水平
      cluster_dt[, condition := factor(condition, levels = c(control_group, treatment_group))]
      
      # 判断是否使用混合模型
      if ("lib_id" %in% colnames(cluster_dt) && uniqueN(cluster_dt$lib_id) > 1) {
        # 混合模型
        model <- tryCatch(
          lme4::lmer(test_length ~ condition + (1 | lib_id), data = cluster_dt),
          error = function(e) NULL
        )
        if (!is.null(model)) {
          coefs <- summary(model)$coefficients
          if (nrow(coefs) > 1) {
            estimate <- coefs[2, 1]
            std_error <- coefs[2, 2]
            t_value <- coefs[2, 3]
            if (requireNamespace("lmerTest", quietly = TRUE)) {
              model_test <- lmerTest::lmer(test_length ~ condition + (1 | lib_id), data = cluster_dt)
              p_value <- summary(model_test)$coefficients[2, 5]
            } else {
              # 近似自由度
              df <- nrow(cluster_dt) - length(unique(cluster_dt$lib_id)) - 1
              p_value <- 2 * pt(abs(t_value), df = df, lower.tail = FALSE)
            }
            log2_fc <- estimate
          } else {
            estimate <- std_error <- t_value <- p_value <- log2_fc <- NA_real_
          }
        } else {
          estimate <- std_error <- t_value <- p_value <- log2_fc <- NA_real_
        }
      } else {
        # 普通线性模型
        model <- lm(test_length ~ condition, data = cluster_dt)
        coefs <- summary(model)$coefficients
        if (nrow(coefs) > 1) {
          estimate <- coefs[2, 1]
          std_error <- coefs[2, 2]
          t_value <- coefs[2, 3]
          p_value <- coefs[2, 4]
          log2_fc <- estimate
        } else {
          estimate <- std_error <- t_value <- p_value <- log2_fc <- NA_real_
        }
      }
      
      res <- data.table::data.table(
        cluster_id = cid,
        n_control = n_ctrl, n_treatment = n_trt,
        mean_control = mean_ctrl, mean_treatment = mean_trt,
        median_control = median_ctrl, median_treatment = median_trt,
        sd_control = sd_ctrl, sd_treatment = sd_trt,
        fold_change = fold_change_raw,
        mean_diff = mean_diff_raw,
        median_diff = median_diff_raw,
        log2_fc = if(logscale) log2_fc else log2(fold_change_raw),
        estimate = estimate,
        std_error = std_error,
        t_value = t_value,
        p_value = p_value,
        method = paste("LMM", if(logscale) "on log2 scale" else "on raw scale")
      )
    }
    
    return(res)
  }, mc.cores = mc.cores)
  
  # 合并结果
  results_list <- results_list[!sapply(results_list, is.null)]
  if (length(results_list) == 0) {
    warning("No clusters passed the minimum mRNA filter.")
    return(data.frame())
  }
  
  final_dt <- data.table::rbindlist(results_list, fill = TRUE)
  
  # 多重假设校正 (Benjamini-Hochberg)
  if ("p_value" %in% names(final_dt)) {
    final_dt[, q_value := p.adjust(p_value, method = "BH")]
  }
  
  final_dt <- final_dt[order(cluster_id)]
  message("Completed. Results for ", nrow(final_dt), " clusters.")
  as.data.frame(final_dt)
}



##############################################
#           Results summary function    #
##############################################

#' Generate summary of mRNA-level statistical results for multiple treatments
#' @name summarize_mRNA_level_results_multi
#' @param results Legacy wide-format results with treatment-specific column names.
#' @param alpha Significance threshold for adjusted P values.
#' @param control_group Name of the control condition.
#' @return A summary data frame, a message for empty input, or invisible NULL
#'   when no treatment columns are present. This internal helper expects the
#'   legacy wide format, not the output of Tail.DiffPair().
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

#' APPLE: poly(A) sites and tail-length analysis
#'
#' Align reads, extract poly(A) sites and tails, remove internal priming,
#' cluster and annotate sites, and compare tail lengths and APA usage.
#' @seealso [minimap2()], [Extract_polyAsite()], [Load.PolyA()],
#'   [Remove.IP()], [Cluster.PolyA()], [Annotate.PolyA()], [Filter.PolyA()],
#'   [Map.Tail()], [Tail.PCA()], [Tail.DiffPair()], [Quantify.GeneAPA()]
#' @importFrom methods new is setClass setMethod
#' @importFrom stats as.formula complete.cases lm median p.adjust prcomp pt qnorm sd t.test var wilcox.test
#' @importFrom utils write.table
#' @importFrom dplyr %>% bind_rows group_by summarise filter n mutate full_join left_join rename case_when
#' @importFrom tidyr pivot_wider
#' @importFrom ggplot2 ggplot aes geom_point xlab ylab coord_fixed
#' @importFrom GenomicRanges GRanges findOverlaps
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors queryHits subjectHits
#' @importFrom BiocGenerics start end width which.max
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
#' @importFrom outliers scores
#' @importFrom plyranges join_overlap_intersect_directed as_granges find_overlaps
#' @importFrom rtracklayer import
#' @importFrom Rsamtools indexFa
#' @importFrom bedr bedr
#' @importFrom data.table := uniqueN
#' @importFrom rlang .data
"_PACKAGE"

############################################################
#       Differential alternative polyadenylation          #
############################################################

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

############################################################
#                    Plotting functions                    #
############################################################

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

