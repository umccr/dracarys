#' Wf bcl_convert R6 Class
#'
#' @description
#' Contains methods for reading and processing files output from the UMCCR
#' `bcl_convert` workflow.
#'
#' @examples
#' \dontrun{
#' indir <- file.path()
#' sample_id <- "PTC_ctTSO240429"
#' library_id <- "L2400482"
#' d <- TsoCombinedVariantOutputFile$new(x)
#' d$read()
#' }
#' @export
Wf_bcl_convert <- R6::R6Class(
  "Wf_bcl_convert",
  public = list(
    #' @field indir Input directory containing Reports per assay type
    #' (e.g. /primary_data/240607_A01052_0209_BHLHFTDSXC/2024061140802544/).
    indir = NULL,

    #' @description Create a new Wf_bcl_convert object.
    #' @param indir Input directory containing Reports per assay type.
    initialize = function(indir) {
      self$indir <- indir
    },
    #' @description Print details about the Workflow
    #' @param ... (ignored).
    print = function(...) {

    }
  )
)

#' BclconvertReports R6 Class
#'
#' @description
#' Reads and writes tidy versions of files within the `Reports` directory output
#' from BCLConvert.
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

    #' @description
    #' Reads contents of `Reports` directory output by BCLConvert.
    #'
    #' @return A list of tibbles.
    #' @export
    read = function() {
      p <- self$path
      read_adaptercyclemetrics <- function(x) {
        cnames <- list(
          old = c(
            "Lane", "Sample_ID", "index", "index2", "ReadNumber", "Cycle",
            "NumClustersWithAdapterAtCycle", "% At Cycle"
          ),
          new = c(
            "lane", "sampleid", "barcode", "readnum", "cycle",
            "clustadapt_n", "cycle_pct"
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
          dplyr::mutate(barcode = paste0(.data$index, "-", .data$index2)) |>
          dplyr::select(-c("index", "index2")) |>
          dplyr::relocate("barcode", .after = "Sample_ID") |>
          rlang::set_names(cnames$new)
      }
      read_topunknownbarcodes <- function(x) {
        cnames <- list(
          old = c("Lane", "index", "index2", "# Reads", "% of Unknown Barcodes", "% of All Reads"),
          new = c("lane", "barcode", "reads_n", "unknownbcodes_pct", "reads_pct")
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
          dplyr::mutate(barcode = paste0(.data$index, "-", .data$index2)) |>
          dplyr::select(-c("index", "index2")) |>
          dplyr::relocate("barcode", .after = "Lane") |>
          rlang::set_names(cnames$new)
      }
      read_adaptermetrics <- function(x) {
        cnames <- list(
          old = c(
            "Lane", "Sample_ID", "index", "index2", "ReadNumber",
            "AdapterBases", "SampleBases", "% Adapter Bases"
          ),
          new = c(
            "lane", "sampleid", "barcode", "readnum", "adapter_bases",
            "sample_bases", "adapter_bases_pct"
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
          dplyr::mutate(barcode = ifelse(
            is.na(.data$index), NA_character_, paste0(.data$index, "-", .data$index2)
          )) |>
          dplyr::select(-c("index", "index2")) |>
          dplyr::relocate("barcode", .after = "Sample_ID") |>
          rlang::set_names(cnames$new)
      }
      read_indexhoppingcounts <- function(x) {
        cnames <- list(
          old = c(
            "Lane", "SampleID", "index", "index2", "# Reads",
            "% of Hopped Reads", "% of All Reads"
          ),
          new = c(
            "lane", "sampleid", "barcode",
            "reads_n", "reads_hopped_pct", "reads_pct"
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
          dplyr::mutate(barcode = paste0(.data$index, "-", .data$index2)) |>
          dplyr::select(-c("index", "index2")) |>
          dplyr::relocate("barcode", .after = "SampleID") |>
          rlang::set_names(cnames$new)
      }
      read_demultiplexstats <- function(x) {
        cnames <- list(
          old = c(
            "Lane", "SampleID", "Index", "# Reads", "# Perfect Index Reads",
            "# One Mismatch Index Reads", "# Two Mismatch Index Reads",
            "% Reads", "% Perfect Index Reads", "% One Mismatch Index Reads",
            "% Two Mismatch Index Reads"
          ),
          new = c(
            "lane", "sampleid", "barcode", "reads_n", "perfect_idxreads_n",
            "one_mismatch_idxreads_n", "two_mismatch_idxreads_n",
            "reads_pct", "perfect_idxreads_pct",
            "one_mismatch_idxreads_pct", "two_mismatch_idxreads_pct"
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
      }
      read_fastqlist <- function(x) {
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
          tidyr::pivot_longer(c("Read1File", "Read2File"), names_to = "readnum", values_to = "filepath") |>
          dplyr::mutate(readnum = sub("Read(.)File", "\\1", .data$readnum)) |>
          rlang::set_names(cnames$new)
      }
      # now return all as list elements
      ac <- read_adaptercyclemetrics(file.path(p, "Adapter_Cycle_Metrics.csv"))
      am <- read_adaptermetrics(file.path(p, "Adapter_Metrics.csv"))
      ds <- read_demultiplexstats(file.path(p, "Demultiplex_Stats.csv"))
      ih <- read_indexhoppingcounts(file.path(p, "Index_Hopping_Counts.csv"))
      ub <- read_topunknownbarcodes(file.path(p, "Top_Unknown_Barcodes.csv"))
      fq <- read_fastqlist(file.path(p, "fastq_list.csv"))
      list(
        adapter_cycle_metrics = ac,
        adapter_metrics = am,
        demultiplex_stats = ds,
        index_hopping_counts = ih,
        top_unknown_barcodes = ub,
        fastq_list = fq
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
    write = function(d, out_dir = NULL, prefix, out_format = "tsv", drid = NULL) {
      if (!is.null(out_dir)) {
        prefix <- file.path(out_dir, prefix)
      }
      d_write <- d |>
        tibble::enframe(name = "section") |>
        dplyr::rowwise() |>
        dplyr::mutate(
          section_low = tolower(.data$section),
          p = glue("{prefix}_{.data$section_low}"),
          out = list(write_dracarys(obj = .data$value, prefix = .data$p, out_format = out_format, drid = drid))
        ) |>
        dplyr::ungroup() |>
        dplyr::select("section", "value") |>
        tibble::deframe()
      invisible(d_write)
    }
  )
)
