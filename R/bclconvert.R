#' bcl_convert Wf R6 Class
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
#' b <- BclconvertReports$new(here::here("nogit/bcl_convert/WGS_TsqNano/Reports"))
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
      req_fnames <- c(
        "Adapter_Metrics.csv", "Demultiplex_Stats.csv",
        "Index_Hopping_Counts.csv", "Top_Unknown_Barcodes.csv"
      )
      assertthat::assert_that(
        all(req_fnames %in% self$contents[["bname"]])
      )
      .read_topunknownbarcodes <- function(x) {
        d <- readr::read_csv(x, col_types = "cccd")
        assertthat::assert_that(all(colnames(d) == c("Lane", "index", "index2", "# Reads")))
        d |>
          rlang::set_names(c("lane", "index1", "index2", "n_reads")) |>
          dplyr::mutate(barcode = glue("{.data$index1}-{.data$index2}") |> as.character()) |>
          dplyr::select("lane", "barcode", "n_reads")
      }
      .read_adaptermetrics <- function(x) {
        d <- readr::read_csv(x, col_types = "ccccddddd")
        old_nms <- c(
          "Lane", "Sample_ID", "index", "index2", "R1_AdapterBases",
          "R1_SampleBases", "R2_AdapterBases", "R2_SampleBases", "# Reads"
        )
        assertthat::assert_that(all(colnames(d) == old_nms))
        d |>
          dplyr::rename(
            index1 = "index", n_reads = "# Reads", SampleID = "Sample_ID", lane = "Lane"
          ) |>
          dplyr::mutate(barcode = ifelse(
            is.na(.data$index1), NA_character_, glue("{.data$index1}-{.data$index2}")
          )) |>
          dplyr::select(
            "lane", "SampleID", "barcode", "n_reads",
            "R1_AdapterBases", "R2_AdapterBases",
            "R1_SampleBases", "R2_SampleBases"
          )
      }
      .read_indexhoppingcounts <- function(x) {
        d <- readr::read_csv(x, col_types = "ccccd")
        old_nms <- c("Lane", "SampleID", "index", "index2", "# Reads")
        assertthat::assert_that(all(colnames(d) == old_nms))
        d |>
          dplyr::rename(index1 = "index", n_reads = "# Reads", lane = "Lane") |>
          dplyr::mutate(barcode = glue("{.data$index1}-{.data$index2}")) |>
          dplyr::select("lane", "SampleID", "barcode", "n_reads")
      }
      .read_demultiplexstats <- function(x) {
        nms <- tibble::tribble(
          ~new_nm, ~old_nm, ~class,
          "lane", "Lane", "c",
          "SampleID", "SampleID", "c",
          "barcode", "Index", "c",
          "n_reads", "# Reads", "d",
          "n_perfect_idxreads", "# Perfect Index Reads", "d",
          "n_one_mismatch_idxreads", "# One Mismatch Index Reads", "d",
          "n_q30_bases", "# of >= Q30 Bases (PF)", "d",
          "mean_quality_score", "Mean Quality Score (PF)", "d"
        )
        lookup <- tibble::deframe(nms[c("new_nm", "old_nm")])
        d <- readr::read_csv(x, col_types = nms[["class"]])
        assertthat::assert_that(all(colnames(d) == nms[["old_nm"]]))
        d |>
          dplyr::rename(dplyr::all_of(lookup))
      }
      .read_fastqlist <- function(x) {
        nms <- tibble::tribble(
          ~new_nm, ~old_nm, ~class,
          "rgid", "RGID", "c",
          "rgsm", "RGSM", "c",
          "rglb", "RGLB", "c",
          "lane", "Lane", "c",
          "1", "Read1File", "c",
          "2", "Read2File", "c"
        )
        lookup <- tibble::deframe(nms[c("new_nm", "old_nm")])
        d <- readr::read_csv(x, col_types = readr::cols(.default = "c"))
        assertthat::assert_that(all(colnames(d) == nms[["old_nm"]]))
        d |>
          dplyr::rename(dplyr::all_of(lookup)) |>
          tidyr::pivot_longer(c("1", "2"), names_to = "read", values_to = "path")
      }

      am <- .read_adaptermetrics(file.path(p, "Adapter_Metrics.csv"))
      ds <- .read_demultiplexstats(file.path(p, "Demultiplex_Stats.csv"))
      ih <- .read_indexhoppingcounts(file.path(p, "Index_Hopping_Counts.csv"))
      ub <- .read_topunknownbarcodes(file.path(p, "Top_Unknown_Barcodes.csv"))
      fq <- .read_fastqlist(file.path(p, "fastq_list.csv"))
      list(
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
