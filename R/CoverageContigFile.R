#' CoverageContigFile R6 Class
#'
#' @description
#' Contains methods for reading and displaying contents of
#' the `wgs_contig_mean_cov_<phenotype>.csv` file output from DRAGEN.
#' This file contains the estimated coverage for all contigs, and an autosomal
#' estimated coverage.
#'
#' @examples
#' x1 <- system.file("extdata/SEQC-II.wgs_contig_mean_cov_normal.csv", package = "dracarys")
#' x2 <- system.file("extdata/SEQC-II.wgs_contig_mean_cov_tumor.csv", package = "dracarys")
#' cc1 <- CoverageContigFile$new(x1)
#' cc2 <- CoverageContigFile$new(x2)
#' read(cc1)
#' read(cc2, keep_alt = TRUE)
#' plot(cc1)
#'
#' @export
CoverageContigFile <- R6::R6Class("CoverageContigFile", inherit = File, public = list(
  #' @description
  #' Reads the `wgs_contig_mean_cov_<phenotype>.csv` file output from DRAGEN.
  #'
  #' @param keep_alt Keep the ALT + Mito chromosomes?
  #' @return tibble with the following columns:
  #'   - label: file label.
  #'   - chrom: contig name.
  #'   - n_bases: number of bases aligned to contig (excludes bases from
  #'   duplicate marked reads, reads with MAPQ=0, and clipped bases).
  #'   - coverage: col2 / contig length
  read = function(keep_alt = FALSE) {
    x <- self$path
    b <- self$bname()
    suffix <- dplyr::if_else(
      grepl("_normal\\.csv", b), "_N",
      dplyr::if_else(grepl("_tumor\\.csv", b), "_T", "")
    )
    nm <- sub("(.*)\\.wgs_contig_mean_cov.*", "\\1", b)
    label <- paste0(nm, suffix)

    readr::read_csv(x, col_names = c("chrom", "n_bases", "coverage"), col_types = "cdd") |>
      dplyr::filter(
        if (!keep_alt) {
          !grepl("chrM|MT|_|Autosomal|HLA-", .data$chrom)
        } else {
          TRUE
        }
      ) |>
      dplyr::mutate(label = tidyselect::all_of(label)) |>
      dplyr::select(.data$label, .data$chrom, .data$n_bases, .data$coverage)
  },

  #' @description Plots the `wgs_contig_mean_cov_<phenotype>.csv` files.
  #' @param top_alt_n Number of top covered alt contigs to plot per phenotype.
  #' @return A ggplot2 object with chromosomes on X axis, and coverage on Y axis.
  plot = function(top_alt_n = 15) {
    assertthat::assert_that(length(top_alt_n) == 1, top_alt_n >= 0, is.numeric(top_alt_n))
    cov_contig <- self$read(keep_alt = TRUE)

    # Display chr1-22, X, Y at top (M goes to bottom).
    # Display top 20 of the rest, plus rest as 'other', at bottom
    main_chrom1 <- c(1:22, "X", "Y")
    main_chrom2 <- c(paste0("chr", main_chrom1))
    main_chrom <- c(main_chrom1, main_chrom2, "Autosomal regions")

    cov_contig <- cov_contig |>
      dplyr::mutate(
        panel = dplyr::if_else(.data$chrom %in% main_chrom, "main", "alt"),
        panel = factor(.data$panel, levels = c("main", "alt"))
      )

    main_panel <- cov_contig |>
      dplyr::filter(.data$panel == "main") |>
      dplyr::select(.data$label, .data$chrom, .data$coverage, .data$panel)
    alt_panel <- cov_contig |>
      dplyr::filter(.data$panel == "alt") |>
      dplyr::select(.data$label, .data$chrom, .data$coverage, .data$panel)

    top_alt <- alt_panel |>
      dplyr::group_by(.data$label) |>
      dplyr::top_n(top_alt_n, wt = .data$coverage) |>
      dplyr::arrange(dplyr::desc(.data$coverage)) |>
      dplyr::pull(.data$chrom) |>
      unique()

    alt_panel2 <- alt_panel |>
      dplyr::mutate(alt_group = dplyr::if_else(.data$chrom %in% top_alt, "top", "bottom"))

    alt_panel_final <- alt_panel2 |>
      dplyr::group_by(.data$alt_group, .data$label) |>
      dplyr::summarise(mean_cov = mean(.data$coverage)) |>
      dplyr::inner_join(alt_panel2, by = c("alt_group", "label")) |>
      dplyr::mutate(
        chrom = dplyr::if_else(.data$alt_group == "bottom", "OTHER", .data$chrom),
        coverage = dplyr::if_else(.data$alt_group == "bottom", .data$mean_cov, .data$coverage)
      ) |>
      dplyr::distinct() |>
      dplyr::ungroup() |>
      dplyr::select(.data$label, .data$chrom, .data$coverage, .data$panel)

    chrom_fac_levels <- c(main_chrom, "chrM", "MT", top_alt[!top_alt %in% c("chrM", "MT")], "OTHER")
    d <- dplyr::bind_rows(main_panel, alt_panel_final) |>
      dplyr::mutate(chrom = factor(.data$chrom, levels = chrom_fac_levels))

    d |>
      ggplot2::ggplot(
        ggplot2::aes(
          x = .data$chrom, y = .data$coverage,
          colour = .data$label, group = .data$label
        )
      ) +
      ggplot2::geom_point() +
      ggplot2::geom_line() +
      ggplot2::scale_y_continuous(
        limits = c(0, NA), expand = c(0, 0), labels = scales::comma,
        breaks = scales::pretty_breaks(n = 8)
      ) +
      ggplot2::theme_minimal() +
      ggplot2::labs(title = "Mean Coverage Per Chromosome", colour = "Label") +
      ggplot2::xlab("Chromosome") +
      ggplot2::ylab("Coverage") +
      ggplot2::theme(
        legend.position = "top",
        panel.grid.minor = ggplot2::element_blank(),
        panel.grid.major.y = ggplot2::element_blank(),
        strip.background = ggplot2::element_blank(),
        strip.text.x = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1, size = 6),
        plot.title = ggplot2::element_text(colour = "#2c3e50", size = 14, face = "bold"),
        panel.spacing = ggplot2::unit(2, "lines")
      ) +
      ggplot2::facet_wrap(ggplot2::vars(.data$panel), nrow = 2, scales = "free")
  }
))
