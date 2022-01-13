#' PloidyEstimationMetricsFile R6 Class
#'
#' @description
#' Contains methods for reading contents of
#' the `ploidy_estimation_metrics.csv` file output from DRAGEN.
#'
#' @examples
#' x <- system.file("extdata/SEQC-II.ploidy_estimation_metrics.csv", package = "dracarys")
#' pem <- PloidyEstimationMetricsFile$new(x)
#' pem$read() # or read(pem)
#'
#' @export
PloidyEstimationMetricsFile <- R6::R6Class("PloidyEstimationMetricsFile", inherit = File, public = list(
  #' @description
  #' Reads the `ploidy_estimation_metrics.csv` file output from DRAGEN.
  #'
  #' @return tibble with the following columns:
  #'     - label: sample label (inferred from file name)
  #'     - var: variable of interest (e.g. X median coverage)
  #'     - value: value for specific variable (e.g. X median coverage
  #'       variable with a value of  50)
  read = function() {
    x <- self$path
    d <- readr::read_lines(x)
    assertthat::assert_that(grepl("PLOIDY ESTIMATION", d[1]))

    b <- self$bname()
    label <- sub("(.*)\\.ploidy_estimation_metrics.*", "\\1", b)

    d |>
      tibble::as_tibble_col(column_name = "value") |>
      tidyr::separate(.data$value, into = c("dummy1", "dummy2", "var", "value"), sep = ",", convert = FALSE) |>
      dplyr::mutate(label = tidyselect::all_of(label)) |>
      dplyr::select(.data$label, .data$var, .data$value)
  }
))
