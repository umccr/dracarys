#' FineHistFile R6 Class
#'
#' @description
#' Contains methods for reading and displaying contents of
#' the `wgs_fine_hist_<phenotype>.csv` file output from DRAGEN.
#' This file contains two columns: Depth and Overall.
#' The value in the Depth column ranges from 0 to 1000+ and the Overall
#' column indicates the number of loci covered at the corresponding depth.
#'
#' @examples
#' x1 <- system.file("extdata/wgs/SEQC-II.wgs_fine_hist_normal.csv.gz", package = "dracarys")
#' x2 <- system.file("extdata/wgs/SEQC-II.wgs_fine_hist_tumor.csv.gz", package = "dracarys")
#' ch1 <- FineHistFile$new(x1)
#' ch2 <- FineHistFile$new(x2)
#' read(ch1)
#' read(ch2)
#' plot(ch1)
#' plot(ch2)
#'
#' @export
FineHistFile <- R6::R6Class("FineHistFile", inherit = File, public = list(
  #' @description
  #' Reads the `wgs_fine_hist_<phenotype>.csv` file output from DRAGEN.
  #' @return tibble with three columns:
  #'   - label
  #'   - depth
  #'   - number of loci with given depth
  read = function() {
    x <- self$path
    b <- self$bname()
    d <- readr::read_csv(x, col_types = "cd")
    assertthat::assert_that(all(colnames(d) == c("Depth", "Overall")))
    suffix <- dplyr::if_else(
      grepl("_normal\\.csv", b), "_N",
      dplyr::if_else(grepl("_tumor\\.csv", b), "_T", "")
    )
    nm <- sub("(.*)\\.wgs_fine_hist.*", "\\1", b)
    label <- paste0(nm, suffix)

    # there's a max Depth of 2000+, so convert to numeric for easier plotting
    d |>
      dplyr::mutate(
        label = tidyselect::all_of(label),
        Depth = ifelse(grepl("+", .data$Depth), sub("(\\d*)\\+", "\\1", .data$Depth), .data$Depth),
        Depth = as.integer(.data$Depth)
      ) |>
      dplyr::select(.data$label, depth = .data$Depth, n_loci = .data$Overall)
  },

  #' @description Plots the `wgs_fine_hist_<phenotype>.csv` files.
  #' @param x_lim X axis range to plot.
  #' @return A ggplot2 object with depth of coverage on X axis,
  #' and number of loci with that depth on Y axis.
  plot = function(x_lim = c(0, 300)) {
    assertthat::assert_that(length(x_lim) == 2)

    cov <- self$read()

    cov |>
      ggplot2::ggplot(ggplot2::aes(
        x = .data$depth, y = .data$n_loci,
        colour = .data$label, group = .data$label
      )) +
      ggplot2::geom_line() +
      ggplot2::coord_cartesian(xlim = x_lim) +
      ggplot2::scale_y_continuous(labels = scales::comma) +
      ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(n = 8)) +
      ggplot2::theme_minimal() +
      ggplot2::labs(title = "Coverage Distribution", colour = "Label") +
      ggplot2::xlab("Depth of Coverage") +
      ggplot2::ylab("Number of Loci with Given Coverage") +
      ggplot2::theme(
        legend.position = c(0.9, 0.9),
        legend.justification = c(1, 1),
        panel.grid.minor = ggplot2::element_blank(),
        panel.grid.major.y = ggplot2::element_blank(),
        # axis.text.x = ggplot2::element_text(angle = 0, vjust = 1, hjust = 1),
        plot.title = ggplot2::element_text(colour = "#2c3e50", size = 14, face = "bold")
      )
  }
))
