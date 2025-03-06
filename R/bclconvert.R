#' BclconvertReports R6 Class
#'
#' @description
#' Reads and writes tidy versions of files within the `Reports` directory output
#' from BCLConvert v4.2.7. See the DRAGEN v4.2 documentation at
#' https://support-docs.illumina.com/SW/dragen_v42/Content/SW/DRAGEN/OutputFiles.htm.
#'
#' @examples
#' \dontrun{
#' p1 <- "240816_A01052_0220_AHM7VHDSXC/202408195d4f2fc4/Reports"
#' b <- here::here("nogit/bcl_convert", p1) |>
#'   BclconvertReports$new()
#' b$path
#' b$contents
#' d <- b$read()
#' b$write(d, out_dir = tempdir(), prefix = "sampleA", out_format = "tsv")
#' }
#'
#' @export
BclconvertReports <- R6::R6Class(
  "BclconvertReports",
  public = list(
    #' @field path Path to the `Reports` directory.
    #' @field contents Tibble with file path, basename, and size.
    path = NULL,
    contents = NULL,
    #' @description Create a new BclconvertReports object.
    #' @param path Path to the `Reports` directory.
    initialize = function(path = NULL) {
      stopifnot(is.character(path), length(path) == 1)
      self$path <- normalizePath(path)
      self$contents <- fs::dir_info(path) |>
        dplyr::mutate(
          bname = basename(.data$path),
          size = as.character(trimws(.data$size))
        ) |>
        dplyr::select("path", "bname", "size")
    },
    #' @description Print details about the BclconvertReports directory.
    #' @param ... (ignored).
    print = function(...) {
      bnames <- self$contents |>
        dplyr::mutate(
          low = tolower(.data$bname),
        ) |>
        dplyr::arrange(.data$low) |>
        dplyr::mutate(
          n = dplyr::row_number(),
          bn = glue("{.data$n}. {.data$bname} ({.data$size})")
        ) |>
        dplyr::pull("bn")
      cat("#--- BclconvertReports ---#\n")
      cat(glue("Path: {self$path}"), "\n")
      cat("Contents:\n")
      cat(bnames, sep = "\n")
      invisible(self)
    },

    #' @description Read Adapter_Cycle_Metrics.csv file.
    #'
    #' - lane: lane number.
    #' - sampleid: sample ID from sample sheet.
    #' - indexes: index/index2 from sample sheet for this sample.
    #' - read: read number.
    #' - cycle: cycle number.
    #' - cluster_n: number of clusters where the adapter was detected
    #'   to begin precisely at this cycle.
    #' - cluster_pct: percentage of all clusters where the adapter was detected
    #'   to begin precisely at this cycle.
    #' @param x (`character(1)`)\cr
    #'   Path to Adapter_Cycle_Metrics.csv file.
    read_adaptercyclemetrics = function(x) {
      cnames <- list(
        old = c(
          "Lane",
          "Sample_ID",
          "index",
          "index2",
          "ReadNumber",
          "Cycle",
          "NumClustersWithAdapterAtCycle",
          "% At Cycle"
        ),
        new = c(
          "lane",
          "sampleid",
          "indexes",
          "read",
          "cycle",
          "cluster_n",
          "cluster_pct"
        )
      )
      ctypes <- list(
        old = "ccccccdd",
        new = "cccccdd"
      )
      if (!file.exists(x)) {
        return(empty_tbl(cnames$new, ctypes$new))
      }
      d <- readr::read_csv(x, col_types = ctypes$old)
      assertthat::assert_that(all(colnames(d) == cnames$old))
      d |>
        dplyr::mutate(indexes = paste0(.data$index, "-", .data$index2)) |>
        dplyr::select(-c("index", "index2")) |>
        dplyr::relocate("indexes", .after = "Sample_ID") |>
        rlang::set_names(cnames$new)
    },

    #' @description Read Adapter_Metrics.csv file.
    #'
    #' - lane: lane number.
    #' - sampleid: sample ID from sample sheet.
    #' - indexes: index/index2 from sample sheet for this sample.
    #' - readnum: read number.
    #' - adapter_bases: total number of bases trimmed as adapter from the read.
    #' - sample_bases: total number of bases not trimmed from the read.
    #' - adapter_bases_pct: percentage of bases trimmed as adapter from the read.
    #' @param x (`character(1)`)\cr
    #'   Path to Adapter_Metrics.csv file.
    read_adaptermetrics = function(x) {
      cnames <- list(
        old = c(
          "Lane",
          "Sample_ID",
          "index",
          "index2",
          "ReadNumber",
          "AdapterBases",
          "SampleBases",
          "% Adapter Bases"
        ),
        new = c(
          "lane",
          "sampleid",
          "indexes",
          "readnum",
          "adapter_bases",
          "sample_bases",
          "adapter_bases_pct"
        )
      )
      ctypes <- list(
        old = "cccccddd",
        new = "ccccddd"
      )
      if (!file.exists(x)) {
        return(empty_tbl(cnames$new, ctypes$new))
      }
      d <- readr::read_csv(x, col_types = ctypes$old)
      assertthat::assert_that(all(colnames(d) == cnames$old))
      d |>
        dplyr::mutate(
          indexes = ifelse(
            is.na(.data$index),
            NA_character_,
            paste0(.data$index, "-", .data$index2)
          )
        ) |>
        dplyr::select(-c("index", "index2")) |>
        dplyr::relocate("indexes", .after = "Sample_ID") |>
        rlang::set_names(cnames$new)
    },

    #' @description Read Demultiplex_Stats.csv file.
    #'
    #' - lane: lane number.
    #' - sampleid: sample ID from sample sheet.
    #' - indexes: index/index2 from sample sheet for this sample.
    #' - reads_n: total number of pass-filter reads mapping to this sample for the lane.
    #' - perfect_idxreads_n: number of mapped reads with barcodes matching the indexes exactly.
    #' - one_mismatch_idxreads_n: number of mapped reads with barcodes matched with one base mismatched.
    #' - two_mismatch_idxreads_n: number of mapped reads with barcodes matched with two bases mismatched.
    #' - reads_pct: percentage of pass-filter reads mapping to this sample for the lane.
    #' - perfect_idxreads_pct: percentage of mapped reads with barcodes matching the indexess exactly.
    #' - one_mismatch_idxreads_pct: percentage of mapped reads with one mismatch to the indexes.
    #' - two_mismatch_idxreads_pct: percentage of mapped reads with two mismatches to the indexes.
    #' @param x (`character(1)`)\cr
    #'   Path to Demultiplex_Stats.csv file.
    read_demultiplexstats = function(x) {
      cnames <- list(
        old = c(
          "Lane",
          "SampleID",
          "Index",
          "# Reads",
          "# Perfect Index Reads",
          "# One Mismatch Index Reads",
          "# Two Mismatch Index Reads",
          "% Reads",
          "% Perfect Index Reads",
          "% One Mismatch Index Reads",
          "% Two Mismatch Index Reads"
        ),
        new = c(
          "lane",
          "sampleid",
          "indexes",
          "reads_n",
          "perfect_idxreads_n",
          "one_mismatch_idxreads_n",
          "two_mismatch_idxreads_n",
          "reads_pct",
          "perfect_idxreads_pct",
          "one_mismatch_idxreads_pct",
          "two_mismatch_idxreads_pct"
        )
      )
      ctypes <- list(
        old = "cccdddddddd",
        new = "cccdddddddd"
      )
      if (!file.exists(x)) {
        return(empty_tbl(cnames$new, ctypes$new))
      }
      d <- readr::read_csv(x, col_types = ctypes$old)
      assertthat::assert_that(all(colnames(d) == cnames$old))
      d |>
        rlang::set_names(cnames$new)
    },

    #' @description Read Demultiplex_Tile_Stats.csv file.
    #'
    #' - lane: lane number.
    #' - sampleid: sample ID from sample sheet.
    #' - indexes: index/index2 from sample sheet for this sample.
    #' - tile: tile number.
    #' - reads_n: total number of pass-filter reads mapping to this sample for the lane.
    #' - reads_pct: percentage of pass-filter reads mapping to this sample for the lane.
    #' @param x (`character(1)`)\cr
    #'   Path to Demultiplex_Tile_Stats.csv file.
    read_demultiplextilestats = function(x) {
      cnames <- list(
        old = c("Lane", "SampleID", "Index", "Tile", "# Reads", "% Reads"),
        new = c("lane", "sampleid", "indexes", "tile", "reads_n", "reads_pct")
      )
      ctypes <- list(
        old = "ccccdd",
        new = "ccccdd"
      )
      if (!file.exists(x)) {
        return(empty_tbl(cnames$new, ctypes$new))
      }
      d <- readr::read_csv(x, col_types = ctypes$old)
      assertthat::assert_that(all(colnames(d) == cnames$old))
      d |>
        rlang::set_names(cnames$new)
    },

    #' @description Read Quality_Metrics.csv file.
    #'
    #' - lane: lane number.
    #' - sampleid: sample ID from sample sheet.
    #' - indexes: index/index2 from sample sheet for this sample.
    #' - readnum: read number (1 or 2).
    #' - yield: number of bases mapping.
    #' - yieldq30: number of bases with quality score >= 30 mapping.
    #' - qscore_sum: sum of quality scores of bases mapping.
    #' - qscore_mean_pf: mean quality score of bases mapping.
    #' - q30_pct: percentage of bases with quality score >= 30 mapping.
    #' @param x (`character(1)`)\cr
    #'   Path to Quality_Metrics.csv file.
    read_qualitymetrics = function(x) {
      cnames <- list(
        old = c(
          "Lane",
          "SampleID",
          "index",
          "index2",
          "ReadNumber",
          "Yield",
          "YieldQ30",
          "QualityScoreSum",
          "Mean Quality Score (PF)",
          "% Q30"
        ),
        new = c(
          "lane",
          "sampleid",
          "indexes",
          "readnum",
          "yield",
          "yieldq30",
          "qscore_sum",
          "qscore_mean_pf",
          "q30_pct"
        )
      )
      ctypes <- list(
        old = "cccccddddd",
        new = "ccccddddd"
      )
      if (!file.exists(x)) {
        return(empty_tbl(cnames$new, ctypes$new))
      }
      d <- readr::read_csv(x, col_types = ctypes$old)
      assertthat::assert_that(all(colnames(d) == cnames$old))
      d |>
        dplyr::mutate(indexes = paste0(.data$index, "-", .data$index2)) |>
        dplyr::select(-c("index", "index2")) |>
        dplyr::relocate("indexes", .after = "SampleID") |>
        rlang::set_names(cnames$new)
    },

    #' @description Read Quality_Tile_Metrics.csv file.
    #'
    #' - lane: lane number.
    #' - sampleid: sample ID from sample sheet.
    #' - indexes: index/index2 from sample sheet for this sample.
    #' - readnum: read number (1 or 2).
    #' - tile: tile number.
    #' - yield: number of bases mapping.
    #' - yieldq30: number of bases with quality score >= 30 mapping.
    #' - qscore_sum: sum of quality scores of bases mapping.
    #' - qscore_mean_pf: mean quality score of bases mapping.
    #' - q30_pct: percentage of bases with quality score >= 30 mapping.
    #' @param x (`character(1)`)\cr
    #'   Path to Quality_Tile_Metrics.csv file.
    read_qualitytilemetrics = function(x) {
      cnames <- list(
        old = c(
          "Lane",
          "SampleID",
          "index",
          "index2",
          "ReadNumber",
          "Tile",
          "Yield",
          "YieldQ30",
          "QualityScoreSum",
          "Mean Quality Score (PF)",
          "% Q30"
        ),
        new = c(
          "lane",
          "sampleid",
          "indexes",
          "readnum",
          "tile",
          "yield",
          "yieldq30",
          "qscore_sum",
          "qscore_mean_pf",
          "q30_pct"
        )
      )
      ctypes <- list(
        old = "ccccccddddd",
        new = "cccccddddd"
      )
      if (!file.exists(x)) {
        return(empty_tbl(cnames$new, ctypes$new))
      }
      d <- readr::read_csv(x, col_types = ctypes$old)
      assertthat::assert_that(all(colnames(d) == cnames$old))
      d |>
        dplyr::mutate(indexes = paste0(.data$index, "-", .data$index2)) |>
        dplyr::select(-c("index", "index2")) |>
        dplyr::relocate("indexes", .after = "SampleID") |>
        rlang::set_names(cnames$new)
    },

    #' @description Read Index_Hopping_Counts.csv file.
    #'
    #' - lane: lane number.
    #' - sampleid: sample ID from sample sheet.
    #' - indexes: index/index2 from sample sheet for this sample.
    #' - reads_n: total number of pass-filter reads mapping to the indexes.
    #' - reads_hopped_pct: percentage of hopped pass-filter reads mapping to the indexes.
    #' - reads_pct: percentage of all pass-filter reads mapping to the indexes.
    #' @param x (`character(1)`)\cr
    #'   Path to Index_Hopping_Counts.csv file.
    read_indexhoppingcounts = function(x) {
      cnames <- list(
        old = c(
          "Lane",
          "SampleID",
          "index",
          "index2",
          "# Reads",
          "% of Hopped Reads",
          "% of All Reads"
        ),
        new = c(
          "lane",
          "sampleid",
          "indexes",
          "reads_n",
          "reads_hopped_pct",
          "reads_pct"
        )
      )
      ctypes <- list(
        old = "ccccd",
        new = "cccddd"
      )
      if (!file.exists(x)) {
        return(empty_tbl(cnames$new, ctypes$new))
      }
      d <- readr::read_csv(x, col_types = ctypes$old)
      assertthat::assert_that(all(colnames(d) == cnames$old))
      d |>
        dplyr::mutate(indexes = paste0(.data$index, "-", .data$index2)) |>
        dplyr::select(-c("index", "index2")) |>
        dplyr::relocate("indexes", .after = "SampleID") |>
        rlang::set_names(cnames$new)
    },

    #' @description Read Top_Unknown_Barcodes.csv file.
    #'
    #' - lane: lane number.
    #' - indexes: index/index2 of this unlisted sequence.
    #' - reads_n: total number of pass-filter reads mapping to the indexes.
    #' - unknownbcodes_pct: percentage of unknown pass-filter reads mapping to the indexes.
    #' @param x (`character(1)`)\cr
    #'   Path to Top_Unknown_Barcodes.csv file.
    read_topunknownbarcodes = function(x) {
      cnames <- list(
        old = c(
          "Lane",
          "index",
          "index2",
          "# Reads",
          "% of Unknown Barcodes",
          "% of All Reads"
        ),
        new = c("lane", "indexes", "reads_n", "unknownbcodes_pct", "reads_pct")
      )
      ctypes <- list(
        old = "cccddd",
        new = "ccddd"
      )
      if (!file.exists(x)) {
        return(empty_tbl(cnames$new, ctypes$new))
      }
      d <- readr::read_csv(x, col_types = ctypes$old)
      assertthat::assert_that(all(colnames(d) == cnames$old))
      d |>
        dplyr::mutate(indexes = paste0(.data$index, "-", .data$index2)) |>
        dplyr::select(-c("index", "index2")) |>
        dplyr::relocate("indexes", .after = "Lane") |>
        rlang::set_names(cnames$new)
    },

    #' @description Read fastq_list.csv file.
    #'
    #' - rgid: read group.
    #' - rgsm: sample ID.
    #' - rglb: library.
    #' - lane: flow cell lane.
    #' - readnum: read number (1 or 2).
    #' - filepath: path to the FASTQ file.
    #' @param x (`character(1)`)\cr
    #'   Path to fastq_list.csv file.
    read_fastqlist = function(x) {
      cnames <- list(
        old = c("RGID", "RGSM", "RGLB", "Lane", "Read1File", "Read2File"),
        new = c("rgid", "rgsm", "rglb", "lane", "readnum", "filepath")
      )
      ctypes <- list(
        old = c("cccccc"),
        new = c("cccccc")
      )
      if (!file.exists(x)) {
        return(empty_tbl(cnames$new, ctypes$new))
      }
      d <- readr::read_csv(x, col_types = ctypes$old)
      assertthat::assert_that(all(colnames(d) == cnames$old))
      d |>
        tidyr::pivot_longer(
          c("Read1File", "Read2File"),
          names_to = "readnum",
          values_to = "filepath"
        ) |>
        dplyr::mutate(readnum = sub("Read(.)File", "\\1", .data$readnum)) |>
        rlang::set_names(cnames$new)
    },

    #' @description
    #' Reads contents of `Reports` directory output by BCLConvert.
    #'
    #' @return A list of tibbles.
    #' @export
    read = function() {
      # now return all as list elements
      p <- self$path
      list(
        adapter_cycle_metrics = self$read_adaptercyclemetrics(file.path(
          p,
          "Adapter_Cycle_Metrics.csv"
        )),
        adapter_metrics = self$read_adaptermetrics(file.path(
          p,
          "Adapter_Metrics.csv"
        )),
        demultiplex_stats = self$read_demultiplexstats(file.path(
          p,
          "Demultiplex_Stats.csv"
        )),
        demultiplex_tile_stats = self$read_demultiplextilestats(file.path(
          p,
          "Demultiplex_Tile_Stats.csv"
        )),
        quality_metrics = self$read_qualitymetrics(file.path(
          p,
          "Quality_Metrics.csv"
        )),
        quality_tile_metrics = self$read_qualitytilemetrics(file.path(
          p,
          "Quality_Tile_Metrics.csv"
        )),
        index_hopping_counts = self$read_indexhoppingcounts(file.path(
          p,
          "Index_Hopping_Counts.csv"
        )),
        top_unknown_barcodes = self$read_topunknownbarcodes(file.path(
          p,
          "Top_Unknown_Barcodes.csv"
        )),
        fastq_list = self$read_fastqlist(file.path(p, "fastq_list.csv"))
      )
    },

    #' @description
    #' Writes tidied contents of `Reports` directory output by BCLConvert.
    #'
    #' @param d Parsed object from `self$read()`.
    #' @param prefix Prefix of output file(s).
    #' @param out_dir Output directory.
    #' @param out_format Format of output file(s).
    #' @param drid dracarys ID to use for the dataset (e.g. `wfrid.123`, `prid.456`).
    write = function(
      d,
      out_dir = NULL,
      prefix,
      out_format = "tsv",
      drid = NULL
    ) {
      if (!is.null(out_dir)) {
        prefix <- file.path(out_dir, prefix)
      }
      d_write <- d |>
        tibble::enframe(name = "section") |>
        dplyr::rowwise() |>
        dplyr::mutate(
          section_low = tolower(.data$section),
          p = glue("{prefix}_{.data$section_low}"),
          out = list(write_dracarys(
            obj = .data$value,
            prefix = .data$p,
            out_format = out_format,
            drid = drid
          ))
        ) |>
        dplyr::ungroup() |>
        dplyr::select("section", "value") |>
        tibble::deframe()
      invisible(d_write)
    }
  )
)

#' BclconvertReports375 R6 Class
#'
#' @description
#' Reads and writes tidy versions of files within the `Reports` directory output
#' from BCLConvert v4.2.7. See the DRAGEN v4.2 documentation at
#' https://support-docs.illumina.com/SW/dragen_v42/Content/SW/DRAGEN/OutputFiles.htm.
#'
#' @examples
#' \dontrun{
#' p1 <- "nogit/bcl_convert/WGS_TsqNano/Reports"
#' b <- here::here(p1) |>
#'   BclconvertReports375$new()
#' b$path
#' b$contents
#' d <- b$read()
#' b$write(d, out_dir = tempdir(), prefix = "sampleA", out_format = "tsv")
#' }
#'
#' @export
BclconvertReports375 <- R6::R6Class(
  "BclconvertReports375",
  public = list(
    #' @field path Path to the `Reports` directory.
    #' @field contents Tibble with file path, basename, and size.
    path = NULL,
    contents = NULL,
    #' @description Create a new BclconvertReports375 object.
    #' @param path Path to the `Reports` directory.
    initialize = function(path = NULL) {
      stopifnot(is.character(path), length(path) == 1)
      self$path <- normalizePath(path)
      self$contents <- fs::dir_info(path) |>
        dplyr::mutate(
          bname = basename(.data$path),
          size = as.character(trimws(.data$size))
        ) |>
        dplyr::select("path", "bname", "size")
    },
    #' @description Print details about the BclconvertReports375 directory.
    #' @param ... (ignored).
    print = function(...) {
      bnames <- self$contents |>
        dplyr::mutate(
          low = tolower(.data$bname),
        ) |>
        dplyr::arrange(.data$low) |>
        dplyr::mutate(
          n = dplyr::row_number(),
          bn = glue("{.data$n}. {.data$bname} ({.data$size})")
        ) |>
        dplyr::pull("bn")
      cat("#--- BclconvertReports ---#\n")
      cat(glue("Path: {self$path}"), "\n")
      cat("Contents:\n")
      cat(bnames, sep = "\n")
      invisible(self)
    },

    #' @description Read Adapter_Metrics.csv file.
    #'
    #' @param x (`character(1)`)\cr
    #'   Path to Adapter_Metrics.csv file.
    read_adaptermetrics = function(x) {
      cnames <- list(
        old = c(
          "Lane",
          "Sample_ID",
          "index",
          "index2",
          "R1_AdapterBases",
          "R1_SampleBases",
          "R2_AdapterBases",
          "R2_SampleBases",
          "# Reads"
        ),
        new = c(
          "lane",
          "sampleid",
          "indexes",
          "adapter_bases_r1",
          "sample_bases_r1",
          "adapter_bases_r2",
          "sample_bases_r2",
          "reads_n"
        )
      )
      ctypes <- list(
        old = "ccccddddd",
        new = "cccddddd"
      )
      if (!file.exists(x)) {
        return(empty_tbl(cnames$new, ctypes$new))
      }
      d <- readr::read_csv(x, col_types = ctypes$old)
      assertthat::assert_that(all(colnames(d) == cnames$old))
      d |>
        dplyr::mutate(
          indexes = ifelse(
            is.na(.data$index),
            NA_character_,
            paste0(.data$index, "-", .data$index2)
          )
        ) |>
        dplyr::select(-c("index", "index2")) |>
        dplyr::relocate("indexes", .after = "Sample_ID") |>
        rlang::set_names(cnames$new)
    },

    #' @description Read Demultiplex_Stats.csv file.
    #'
    #' @param x (`character(1)`)\cr
    #'   Path to Demultiplex_Stats.csv file.
    read_demultiplexstats = function(x) {
      cnames <- list(
        old = c(
          "Lane",
          "SampleID",
          "Index",
          "# Reads",
          "# Perfect Index Reads",
          "# One Mismatch Index Reads",
          "# of >= Q30 Bases (PF)",
          "Mean Quality Score (PF)"
        ),
        new = c(
          "lane",
          "sampleid",
          "indexes",
          "reads_n",
          "perfect_idxreads_n",
          "one_mismatch_idxreads_n",
          "q30_bases_n",
          "qscore_mean_pf"
        )
      )
      ctypes <- list(
        old = "cccddddd",
        new = "cccddddd"
      )
      if (!file.exists(x)) {
        return(empty_tbl(cnames$new, ctypes$new))
      }
      d <- readr::read_csv(x, col_types = ctypes$old)
      assertthat::assert_that(all(colnames(d) == cnames$old))
      d |>
        rlang::set_names(cnames$new)
    },

    #' @description Read Index_Hopping_Counts.csv file.
    #'
    #' @param x (`character(1)`)\cr
    #'   Path to Index_Hopping_Counts.csv file.
    read_indexhoppingcounts = function(x) {
      cnames <- list(
        old = c("Lane", "SampleID", "index", "index2", "# Reads"),
        new = c("lane", "sampleid", "indexes", "reads_n")
      )
      ctypes <- list(
        old = "ccccd",
        new = "cccd"
      )
      if (!file.exists(x)) {
        return(empty_tbl(cnames$new, ctypes$new))
      }
      d <- readr::read_csv(x, col_types = ctypes$old)
      assertthat::assert_that(all(colnames(d) == cnames$old))
      d |>
        dplyr::mutate(indexes = paste0(.data$index, "-", .data$index2)) |>
        dplyr::select(-c("index", "index2")) |>
        dplyr::relocate("indexes", .after = "SampleID") |>
        rlang::set_names(cnames$new)
    },

    #' @description Read Top_Unknown_Barcodes.csv file.
    #'
    #' @param x (`character(1)`)\cr
    #'   Path to Top_Unknown_Barcodes.csv file.
    read_topunknownbarcodes = function(x) {
      cnames <- list(
        old = c("Lane", "index", "index2", "# Reads"),
        new = c("lane", "indexes", "reads_n")
      )
      ctypes <- list(
        old = "cccd",
        new = "ccd"
      )
      if (!file.exists(x)) {
        return(empty_tbl(cnames$new, ctypes$new))
      }
      d <- readr::read_csv(x, col_types = ctypes$old)
      assertthat::assert_that(all(colnames(d) == cnames$old))
      d |>
        dplyr::mutate(indexes = paste0(.data$index, "-", .data$index2)) |>
        dplyr::select(-c("index", "index2")) |>
        dplyr::relocate("indexes", .after = "Lane") |>
        rlang::set_names(cnames$new)
    },

    #' @description Read fastq_list.csv file.
    #'
    #' - rgid: read group.
    #' - rgsm: sample ID.
    #' - rglb: library.
    #' - lane: flow cell lane.
    #' - readnum: read number (1 or 2).
    #' - filepath: path to the FASTQ file.
    #' @param x (`character(1)`)\cr
    #'   Path to fastq_list.csv file.
    read_fastqlist = function(x) {
      cnames <- list(
        old = c("RGID", "RGSM", "RGLB", "Lane", "Read1File", "Read2File"),
        new = c("rgid", "rgsm", "rglb", "lane", "readnum", "filepath")
      )
      ctypes <- list(
        old = c("cccccc"),
        new = c("cccccc")
      )
      if (!file.exists(x)) {
        return(empty_tbl(cnames$new, ctypes$new))
      }
      d <- readr::read_csv(x, col_types = ctypes$old)
      assertthat::assert_that(all(colnames(d) == cnames$old))
      d |>
        tidyr::pivot_longer(
          c("Read1File", "Read2File"),
          names_to = "readnum",
          values_to = "filepath"
        ) |>
        dplyr::mutate(readnum = sub("Read(.)File", "\\1", .data$readnum)) |>
        rlang::set_names(cnames$new)
    },

    #' @description
    #' Reads contents of `Reports` directory output by BCLConvert.
    #'
    #' @return A list of tibbles.
    #' @export
    read = function() {
      # now return all as list elements
      p <- self$path
      list(
        adapter_metrics = self$read_adaptermetrics(file.path(
          p,
          "Adapter_Metrics.csv"
        )),
        demultiplex_stats = self$read_demultiplexstats(file.path(
          p,
          "Demultiplex_Stats.csv"
        )),
        index_hopping_counts = self$read_indexhoppingcounts(file.path(
          p,
          "Index_Hopping_Counts.csv"
        )),
        top_unknown_barcodes = self$read_topunknownbarcodes(file.path(
          p,
          "Top_Unknown_Barcodes.csv"
        )),
        fastq_list = self$read_fastqlist(file.path(p, "fastq_list.csv"))
      )
    },

    #' @description
    #' Writes tidied contents of `Reports` directory output by BCLConvert.
    #'
    #' @param d Parsed object from `self$read()`.
    #' @param prefix Prefix of output file(s).
    #' @param out_dir Output directory.
    #' @param out_format Format of output file(s).
    #' @param drid dracarys ID to use for the dataset (e.g. `wfrid.123`, `prid.456`).
    write = function(
      d,
      out_dir = NULL,
      prefix,
      out_format = "tsv",
      drid = NULL
    ) {
      if (!is.null(out_dir)) {
        prefix <- file.path(out_dir, prefix)
      }
      d_write <- d |>
        tibble::enframe(name = "section") |>
        dplyr::rowwise() |>
        dplyr::mutate(
          section_low = tolower(.data$section),
          p = glue("{prefix}_{.data$section_low}"),
          out = list(write_dracarys(
            obj = .data$value,
            prefix = .data$p,
            out_format = out_format,
            drid = drid
          ))
        ) |>
        dplyr::ungroup() |>
        dplyr::select("section", "value") |>
        tibble::deframe()
      invisible(d_write)
    }
  )
)
