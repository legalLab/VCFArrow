#' @title read_vcf
#'
#' @description
#' Read VCF file and store its content within a VCFArrow object
#'
#' @author Tomas Hrbek April 2026
#'
#' @param vcf_file -> VCF file
#' @param chunk_size -> number of variants to read in at a time, default 50000 (integer)
#'
#' @return VCFArrow object
#'
#' @details
#' This function read a VCF file into an S4 class object in chunks,
#' returning a VCFArrow object.
#' It accepts both uncompressed and gz compressed files.
#' The GT field is stored as an Apache Arrow in Long format.
#' Various metrics are precalculated for fast and easy filtering.
#' GT slot content stored in a TEMP directory for lazy loading.
#'
#' @examples
#' vcf_filter_rank(vcf_file = my_vcf, chunk_size = 100000)
#' vcf_filter_rank(my_vcf, 100000)
#' vcf_filter_rank(my_vcf)
#'

read_vcf <- function(vcf_file, chunk_size = 100000) {

  # helper function
  parse_format_header <- function(header) {

    fmt_lines <- grep("^##FORMAT=", header, value = TRUE)

    fmt_ids <- sub(".*ID=([^,]+).*", "\\1", fmt_lines)
    fmt_types <- sub(".*Type=([^,]+).*", "\\1", fmt_lines)

    return(data.frame(ID = fmt_ids, Type = fmt_types, stringsAsFactors = FALSE))
  }

  # generate a unique path name within the session tmp folder
  tmp_dir <- tempfile(pattern = "arrow_vcf_")

  # physically create the directory at that path
  dir.create(tmp_dir)

  # set up a conditional file connection
  con <- if (grepl("\\.gz$", vcf_file)) gzfile(vcf_file, "rt") else file(vcf_file, "rt")
  on.exit(close(con))

  # read header
  header <- character()
  repeat {
    line <- readLines(con, n = 1)

    if (length(line) == 0) {
      cli::cli_abort("Unexpected EOF before header")
    }
    if (startsWith(line, "##")) {
      header <- c(header, line)
    } else if (startsWith(line, "#CHROM")) {
      header <- c(header, line)
      break
    } else {
      cli::cli_abort("Malformed VCF: missing #CHROM line")
    }
  }

  # extract samples
  header_fields <- strsplit(header[length(header)], "\t", fixed = TRUE)[[1]]
  samples <- header_fields[-(1:9)]

  # placeholder for group information
  groups <- rep.int(NA_character_, length(samples))

  # storage
  variant_buffer <- list()
  info_buffer <- list()
  format_buffer <- list()
  chunk_id <- 1

  # chunk size message
  cli::cli_alert_info("VCF is being read in chunks of {chunk_size} variants")

  # set up progress bar
  cli::cli_progress_bar("Reading in VCF chunk", total = NA)

  repeat {
    lines <- readLines(con, n = chunk_size)
    if (length(lines) == 0) break

    # progress update
    cli::cli_progress_update()

    n <- length(lines)

    parsed <- parse_vcf_cpp(lines, length(samples))

    # build long-format table
    n <- nrow(parsed$a1)
    n_samples <- length(samples)

    # repeat row ids for each sample
    row_ids <- as.integer(rep(seq_len(n) + (chunk_id - 1) * chunk_size, each = n_samples))

    # repeat sample names for each variant
    sample_ids <- rep(samples, times = n)

    # flatten matrices (column-major → need transpose first)
    a1_vec <- as.vector(t(parsed$a1))
    a2_vec <- as.vector(t(parsed$a2))
    phased_vec <- as.vector(t(parsed$phased))
    fmt_vec <- as.vector(t(parsed$fmt))
    DP_vec <- as.vector(t(parsed$DP))
    GQ_vec <- as.vector(t(parsed$GQ))
    ADR_vec <- as.vector(t(parsed$ADR))

    gt_long <- data.frame(
      .row_id = row_ids,
      sample = sample_ids,
      a1 = a1_vec,
      a2 = a2_vec,
      phased = phased_vec,
      fmt = fmt_vec,
      DP = DP_vec,
      GQ = GQ_vec,
      ADR = ADR_vec,
      stringsAsFactors = FALSE
    )

    arrow::write_feather(gt_long, file.path(tmp_dir, paste0("chunk_", chunk_id, ".arrow")))
    # option to store as Apache parquet
    #arrow::write_parquet(gt_long, file.path(tmp_dir, paste0("chunk_", chunk_id, ".parquet")))

    # FORMAT is per-variant → for memory efficiency keep in a separate FORMAT lookup
    format_df <- data.frame(
      FORMAT = parsed$format,
      .row_id = as.integer(seq_len(n) + (chunk_id - 1) * chunk_size),
      stringsAsFactors = FALSE
    )
    format_buffer[[chunk_id]] <- format_df

    variant_buffer[[chunk_id]] <- as.matrix(parsed$variants)
    info_buffer[[chunk_id]] <- as.character(parsed$info)

    chunk_id <- chunk_id + 1
  }

  # end of progress bar
  cli::cli_progress_done()

  # combine FORMAT field
  format_df_all <- do.call(rbind, format_buffer)

  # combine metadata
  variants_df <- as.data.frame(do.call(rbind, variant_buffer), stringsAsFactors = FALSE)
  colnames(variants_df) <- c("CHROM","POS","ID","REF","ALT","QUAL","FILTER")
  variants_df$POS <- suppressWarnings(as.integer(variants_df$POS))
  info_vec <- unlist(info_buffer, use.names = FALSE)

  # detect if Rk exists anywhere and make a slot
  has_rk <- any(stringr::str_detect(info_vec, "(^|;)Rk="))

  if (has_rk) {
    variants_df$Rk <- suppressWarnings(
      as.numeric(stringr::str_remove(
        stringr::str_extract(info_vec, "Rk=[^;]+"), "Rk="))
    )
  } else {
    variants_df$Rk <- NA_real_
  }

  # extract variant information for filtering
  variants_df <- variants_df |>
    dplyr::mutate(
      n_alt = stringr::str_count(ALT, ",") + 1,
      is_biallelic = n_alt == 1,
      is_indel = (nchar(REF) != 1 |
                    vapply(strsplit(ALT, ","), function(x) any(nchar(x) != 1), logical(1)))
    )

  variants_df$.row_id <- as.integer(seq_len(nrow(variants_df)))

  # open Arrow Dataset (lazy, unified view)
  gt_arrow <- arrow::open_dataset(tmp_dir, format = "feather")
  # option to store as Apache parquet
  #gt_arrow <- arrow::open_dataset(tmp_dir, format = "parquet")

  vcf_arrow <- .new_vcfarrow(
    header,
    info_vec,
    format_df_all,
    variants_df,
    gt_arrow,
    samples,
    groups,
    tmp_dir
  )

  return(vcf_arrow)
}
