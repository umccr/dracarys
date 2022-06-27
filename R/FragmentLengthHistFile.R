#' FragmentLengthHistFile R6 Class
#'
#' @description
#' Contains methods for reading and plotting contents of
#' the `fragment_length_hist.csv` file output from DRAGEN.
#'
#' @examples
#' x <- system.file("extdata/wgs/SEQC-II.fragment_length_hist.csv.gz", package = "dracarys")
#' fl <- FragmentLengthHistFile$new(x)
#' fl$read() # or read(fl)
#' fl$plot() # or plot(fl)
#' @export
FragmentLengthHistFile <- R6::R6Class("FragmentLengthHistFile", inherit = File, public = list(
  #' @description Reads the `fragment_length_hist.csv` file, which contains the
  #' fragment length distribution for each sample.
  #' @return A tibble with the following columns:
  #' - fragmentLength: estimated fragment length
  #' - count: number of reads with estimated fragment length
  #' - sample: name of sample
  read = function() {
    x <- self$path
    d <- readr::read_lines(x)
    assertthat::assert_that(grepl("#Sample", d[1]))

    d |>
      tibble::enframe() |>
      dplyr::mutate(
        sample = dplyr::if_else(
          grepl("#Sample", .data$value),
          sub("#Sample: (.*)", "\\1", .data$value),
          NA_character_
        )
      ) |>
      tidyr::fill(.data$sample, .direction = "down") |>
      dplyr::filter(!grepl("#Sample: |FragmentLength,Count", .data$value)) |>
      tidyr::separate(.data$value, c("fragmentLength", "count"), convert = TRUE) |>
      dplyr::select(-.data$name)
  },


  #' @description Plots the fragment length distributions as given in the
  #' `fragment_length_hist.csv` file.
  #'
  #' @param min_count Minimum read count to be plotted (Default: 10).
  #' @return A ggplot2 plot containing fragment lengths on X axis and read counts
  #'   on Y axis for each sample.
  plot = function(min_count = 10) {
    assertthat::assert_that(is.numeric(min_count), min_count >= 0)
    d <- self$read() |>
      dplyr::bind_rows() |>
      dplyr::filter(.data$count >= min_count)

    d |>
      ggplot2::ggplot(ggplot2::aes(x = .data$fragmentLength, y = .data$count, colour = sample)) +
      ggplot2::geom_line() +
      ggplot2::labs(title = "Fragment Length Distribution") +
      ggplot2::xlab("Fragment Length (bp)") +
      ggplot2::ylab(glue::glue("Read Count (min: {min_count})")) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        legend.position = c(0.9, 0.9),
        legend.justification = c(1, 1),
        panel.grid.minor = ggplot2::element_blank(),
        panel.grid.major = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(colour = "#2c3e50", size = 14, face = "bold")
      )
  }
))
