#' BclconvertReports R6 Class
#'
#' @description
#' Reads the files within the `Reports` directory output from BCLConvert.
#'
#' @examples
#' \dontrun{
#' b <- BclconvertReports$new()
#' }
#'
#' @export
BclconvertReports <- R6::R6Class(
  "BclconvertReports",
  inherit = File,
  public = list(
    #' @description
    #' Reads contents of `Reports` directory output by BCLConvert.
    #'
    #' @return A list of tibbles.
    #' @export
    read = function() {
      x <- self$path
    },

    #' @description
    #' Writes a tidy version of the `wgs_contig_mean_cov_<phenotype>.csv` file output
    #' from DRAGEN.
    #'
    #' @param d Parsed object from `self$read()`.
    #' @param prefix Prefix of output file(s).
    #' @param out_dir Output directory.
    #' @param out_format Format of output file(s).
    #' @param drid dracarys ID to use for the dataset (e.g. `wfrid.123`, `prid.456`).
    write = function(d, out_dir = NULL, prefix, out_format = "tsv", drid = NULL) {
    }
  )
)

#' Read BCLConvert Top Unknown Barcodes
#'
#' Reads the `Top_Unknown_Barcodes.csv` file in the `Reports` directory
#' output by BCLConvert.
#'
#' @param x Path to `Top_Unknown_Barcodes.csv` file.
#'
#' @return Tibble
#'
#' @examples
#' \dontrun{
#' x <- here::here("nogit/bcl_convert/WGS_TsqNano/Reports/Top_Unknown_Barcodes.csv")
#' bclconvert_read_topunknownbarcodes(x)
#' }
#' @export
bclconvert_read_topunknownbarcodes <- function(x) {
  d <- readr::read_csv(x, col_types = "cccd")
  assertthat::assert_that(all(colnames(d) == c("Lane", "index", "index2", "# Reads")))
  d |>
    rlang::set_names(c("lane", "index1", "index2", "n_reads")) |>
    dplyr::mutate(barcode = glue("{.data$index1}-{.data$index2}")) |>
    dplyr::select("lane", "barcode", "n_reads")
}

#' Read BCLConvert Adapter Metrics
#'
#' Reads the `Adapter_Metrics.csv` file in the `Reports` directory
#' output by BCLConvert.
#'
#' @param x Path to `Adapter_Metrics.csv` file.
#'
#' @return Tibble
#'
#' @examples
#' \dontrun{
#' x <- here::here("nogit/bcl_convert/WGS_TsqNano/Reports/Adapter_Metrics.csv")
#' bclconvert_read_adaptermetrics(x)
#' }
#' @export
bclconvert_read_adaptermetrics <- function(x) {
  d <- readr::read_csv(x, col_types = "ccccddddd")
  old_nms <- c(
    "Lane", "Sample_ID", "index", "index2", "R1_AdapterBases",
    "R1_SampleBases", "R2_AdapterBases", "R2_SampleBases", "# Reads"
  )
  assertthat::assert_that(all(colnames(d) == old_nms))
  d |>
    dplyr::rename(index1 = "index", n_reads = "# Reads") |>
    dplyr::mutate(barcode = ifelse(
      is.na(.data$index1), NA_character_, glue("{.data$index1}-{.data$index2}")
    )) |>
    dplyr::select(
      "Lane", "Sample_ID", "barcode", "n_reads",
      "R1_AdapterBases", "R2_AdapterBases",
      "R1_SampleBases", "R2_SampleBases"
    )
}
