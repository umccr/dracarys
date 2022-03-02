#' VCMetricsFile R6 Class
#'
#' @description
#' Contains methods for reading and displaying contents of
#' the `vc_metrics.csv` file output from DRAGEN, which contains variant calling metrics.
#' reported for each sample in multi sample VCF and gVCF files.
#' Based on the run case, metrics are reported either as standard VARIANT
#' CALLER or JOINT CALLER. Metrics are reported both for the raw
#' (PREFILTER) and hard filtered (POSTFILTER) VCFs.
#' PON (Panel of Normals) and COSMIC filtered variants are counted as
#' though they are PASS variants in the POSTFILTER VCFs metrics,
#' which may result in higher than expected variant counts in the
#' POSTFILTER VCF metrics.
#'
#' @examples
#' x <- system.file("extdata/wgs/SEQC-II.vc_metrics.csv.gz", package = "dracarys")
#' vm <- VCMetricsFile$new(x)
#' vm$read() # or read(vm)
#'
#' @export
VCMetricsFile <- R6::R6Class("VCMetricsFile", inherit = File, public = list(
  #' @description
  #' Reads the `vc_metrics.csv` file output from DRAGEN.
  #'
  #' @return tibble with the following columns:
  #'   - category
  #'   - sample
  #'   - var: variable name
  #'   - count: count value
  #'   - pct: percent value
  read = function() {
    x <- self$path
    d <- readr::read_lines(x)
    assertthat::assert_that(grepl("VARIANT CALLER", d[1]))

    d |>
      tibble::enframe(name = "name", value = "value") |>
      tidyr::separate(.data$value,
        into = c("category", "sample", "extra"), sep = ",", extra = "merge"
      ) |>
      tidyr::separate(.data$extra,
        into = c("var", "value"), sep = ",", extra = "merge"
      ) |>
      tidyr::separate(.data$value,
        into = c("count", "pct"), sep = ",", fill = "right", convert = TRUE
      ) |>
      dplyr::mutate(
        category = dplyr::case_when(
          grepl("SUMMARY", .data$category) ~ "summary",
          grepl("PREFILTER", .data$category) ~ "prefilter",
          grepl("POSTFILTER", .data$category) ~ "postfilter",
          TRUE ~ "unknown"
        )
      ) |>
      dplyr::select(.data$category, .data$sample, .data$var, .data$count, .data$pct)
  }
))
