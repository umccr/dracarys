#' TimeMetricsFile R6 Class
#'
#' @description
#' Contains methods for reading contents of
#' the `time_metrics.csv` file output from DRAGEN, which contains
#' a breakdown of the run duration for each DRAGEN process.
#'
#' @examples
#' x <- system.file("extdata/wgs/SEQC-II.time_metrics.csv.gz", package = "dracarys")
#' tm <- TimeMetricsFile$new(x)
#' tm$read() # or read(tm)
#' @export
TimeMetricsFile <- R6::R6Class("TimeMetricsFile", inherit = File, public = list(
  #' @description Reads the `time_metrics.csv` file.
  #' @return tibble with the following columns:
  #'   - Step: DRAGEN step
  #'   - Time: time in HH:MM
  read = function() {
    x <- self$path
    cn <- c("dummy1", "dummy2", "Step", "time_hrs", "time_sec")
    ct <- readr::cols(.default = "c", time_hrs = readr::col_time(format = "%T"), time_sec = "d")
    d <- readr::read_csv(x, col_names = cn, col_types = ct)
    assertthat::assert_that(d$dummy1[1] == "RUN TIME", is.na(d$dummy2[1]))
    assertthat::assert_that(inherits(d$time_hrs, "hms"))

    label <- sub(".time_metrics.csv.*", "", basename(x))

    d |>
      dplyr::mutate(
        Step = tools::toTitleCase(sub("Time ", "", .data$Step)),
        Label = label,
        Time = substr(.data$time_hrs, 1, 5)
      ) |>
      dplyr::select(.data$Label, .data$Step, .data$Time)
  }
))

#' Process Multiple TimeMetricsFile Objects
#'
#' Processes multiple TimeMetricsFile objects.
#'
#' @param x Atomic vector with one or more TimeMetricsFile objects.
#' @param id ID for each input, which is used to disambiguate files
#' generated from same samples. Default: index from 1 to length of `x`.
#' @return tibble with the following columns:
#'   - Step: DRAGEN step
#'   - Time: time in HH:MM
#'
#' @examples
#' x <- system.file("extdata/wgs/SEQC-II.time_metrics.csv.gz", package = "dracarys")
#' x <- TimeMetricsFile$new(x)
#' (tm <- time_metrics_process(c(x, x), id = c("run1", "run2")))
#'
#' @testexamples
#' expect_equal(nrow(tm), 2)
#'
#' @export
time_metrics_process <- function(x, id = seq_len(length(x))) {
  assertthat::assert_that(all(purrr::map_lgl(x, ~ inherits(.x, "TimeMetricsFile"))))
  x |>
    purrr::map(read) |>
    purrr::set_names(id) |>
    dplyr::bind_rows(.id = "ID") |>
    tidyr::pivot_wider(id_cols = c("ID", "Label"), names_from = "Step", values_from = "Time") |>
    dplyr::relocate(.data$`Total Runtime`, .after = .data$Label)
}
