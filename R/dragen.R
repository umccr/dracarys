#' ContigMeanCovFile R6 Class
#'
#' @description
#' Contains methods for reading and displaying contents of
#' the `wgs_contig_mean_cov_<phenotype>.csv` file output from DRAGEN.
#' This file contains the estimated coverage for all contigs, and an autosomal
#' estimated coverage.
#'
#' @examples
#' x1 <- system.file("extdata/wgs/SEQC-II.wgs_contig_mean_cov_normal.csv.gz", package = "dracarys")
#' x2 <- system.file("extdata/wgs/SEQC-II.wgs_contig_mean_cov_tumor.csv.gz", package = "dracarys")
#' cc1 <- ContigMeanCovFile$new(x1)
#' cc2 <- ContigMeanCovFile$new(x2)
#' d1 <- cc1$read()
#' d2 <- cc2$read()
#' cc1$write(d1, out_dir = tempdir(), prefix = "seqc_n", out_format = "tsv")
#' cc2$write(d2, out_dir = tempdir(), prefix = "seqc_t", out_format = "tsv")
#'
#' cc1$plot(d1)
#' cc2$plot(d2)
#'
#' @export
ContigMeanCovFile <- R6::R6Class(
  "ContigMeanCovFile",
  inherit = File,
  public = list(
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
    read = function(keep_alt = TRUE) {
      x <- self$path
      readr::read_csv(x, col_names = c("chrom", "n_bases", "coverage"), col_types = "cdd") |>
        dplyr::filter(
          if (!keep_alt) {
            !grepl("chrM|MT|_|Autosomal|HLA-|EBV", .data$chrom)
          } else {
            TRUE
          }
        )
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
      if (!is.null(out_dir)) {
        prefix <- file.path(out_dir, prefix)
      }
      write_dracarys(obj = d, prefix = prefix, out_format = out_format, drid = drid)
    },


    #' @description Plots the `wgs_contig_mean_cov_<phenotype>.csv` files.
    #' @param d Parsed object from `self$read()`.
    #' @param top_alt_n Number of top covered alt contigs to plot per phenotype.
    #' @return A ggplot2 object with chromosomes on X axis, and coverage on Y axis.
    plot = function(d, top_alt_n = 15) {
      assertthat::assert_that(top_alt_n >= 0)

      # Display chr1-22, X, Y at top (M goes to bottom).
      # Display top 20 of the rest, plus rest as 'other', at bottom
      main_chrom1 <- c(1:22, "X", "Y")
      main_chrom2 <- c(paste0("chr", main_chrom1))
      main_chrom <- c(main_chrom1, main_chrom2, "Autosomal regions")

      d <- d |>
        dplyr::mutate(
          panel = dplyr::if_else(.data$chrom %in% main_chrom, "main", "alt"),
          panel = factor(.data$panel, levels = c("main", "alt"))
        ) |>
        dplyr::select("chrom", "coverage", "panel")

      main_panel <- d |> dplyr::filter(.data$panel == "main")
      alt_panel <- d |> dplyr::filter(.data$panel == "alt")
      top_alt <- alt_panel |>
        dplyr::top_n(top_alt_n, wt = .data$coverage) |>
        dplyr::arrange(dplyr::desc(.data$coverage)) |>
        dplyr::pull(.data$chrom) |>
        unique()

      alt_panel2 <- alt_panel |>
        dplyr::mutate(alt_group = dplyr::if_else(.data$chrom %in% top_alt, "top", "bottom"))

      alt_panel_final <- alt_panel2 |>
        dplyr::group_by(.data$alt_group) |>
        dplyr::summarise(mean_cov = mean(.data$coverage)) |>
        dplyr::inner_join(alt_panel2, by = c("alt_group")) |>
        dplyr::mutate(
          chrom = dplyr::if_else(.data$alt_group == "bottom", "OTHER", .data$chrom),
          coverage = dplyr::if_else(.data$alt_group == "bottom", .data$mean_cov, .data$coverage)
        ) |>
        dplyr::distinct() |>
        dplyr::ungroup() |>
        dplyr::select("chrom", "coverage", "panel")

      chrom_fac_levels <- c(main_chrom, "chrM", "MT", top_alt[!top_alt %in% c("chrM", "MT")], "OTHER")
      d <- dplyr::bind_rows(main_panel, alt_panel_final) |>
        dplyr::mutate(chrom = factor(.data$chrom, levels = chrom_fac_levels))

      d |>
        dplyr::mutate(label = "sampleA") |>
        ggplot2::ggplot(
          ggplot2::aes(
            x = .data$chrom, y = .data$coverage, group = .data$label,
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
  )
)

#' CoverageMetricsFile R6 Class
#'
#' @description
#' Contains methods for reading and displaying contents of
#' the `wgs_coverage_metrics_<phenotype>.csv` file output from DRAGEN.
#' This file contains read depth of coverage metrics.
#'
#' @examples
#' x1 <- system.file("extdata/wgs/SEQC-II.wgs_coverage_metrics_normal.csv.gz", package = "dracarys")
#' x2 <- system.file("extdata/wgs/SEQC-II.wgs_coverage_metrics_tumor.csv.gz", package = "dracarys")
#' cm1 <- CoverageMetricsFile$new(x1)
#' cm2 <- CoverageMetricsFile$new(x2)
#' d1 <- read(cm1)
#' d2 <- read(cm2)
#' cm1$write(d1, out_dir = tempdir(), prefix = "seqc_n", out_format = "tsv")
#' cm2$write(d2, out_dir = tempdir(), prefix = "seqc_t", out_format = "tsv")
#'
#' @export
CoverageMetricsFile <- R6::R6Class(
  "CoverageMetricsFile",
  inherit = File,
  public = list(
    #' @description
    #' Reads the `wgs_coverage_metrics_<phenotype>.csv` file output from DRAGEN.
    #'
    #' @return tibble with the following columns:
    #'   - label: file label.
    #'   - var: variable name.
    #'   - var_abbrev: variable abbreviation.
    #'   - pct: percentage value.
    #'   - count: count value.
    read = function() {
      abbrev_nm <- c(
        "Aligned bases"                                       = "bases_aligned_dragen",
        "Aligned bases in genome"                             = "bases_aligned_in_genome_dragen",
        "Average alignment coverage over genome"              = "cov_alignment_avg_over_genome_dragen",
        "Uniformity of coverage (PCT > 0.2*mean) over genome" = "cov_uniformity_over_genome_pct_gt02mean_dragen",
        "PCT of genome with coverage [100x: inf)"             = "cov_genome_pct_100x_inf_dragen",
        "PCT of genome with coverage [ 50x: inf)"             = "cov_genome_pct_50x_inf_dragen",
        "PCT of genome with coverage [ 20x: inf)"             = "cov_genome_pct_20x_inf_dragen",
        "PCT of genome with coverage [ 15x: inf)"             = "cov_genome_pct_15x_inf_dragen",
        "PCT of genome with coverage [ 10x: inf)"             = "cov_genome_pct_10x_inf_dragen",
        "PCT of genome with coverage [  3x: inf)"             = "cov_genome_pct_3x_inf_dragen",
        "PCT of genome with coverage [  1x: inf)"             = "cov_genome_pct_1x_inf_dragen",
        "PCT of genome with coverage [  0x: inf)"             = "cov_genome_pct_0x_inf_dragen",
        "PCT of genome with coverage [ 50x:100x)"             = "cov_genome_pct_50x_100x_dragen",
        "PCT of genome with coverage [ 20x: 50x)"             = "cov_genome_pct_20x_50x_dragen",
        "PCT of genome with coverage [ 15x: 20x)"             = "cov_genome_pct_15x_20x_dragen",
        "PCT of genome with coverage [ 10x: 15x)"             = "cov_genome_pct_10x_15x_dragen",
        "PCT of genome with coverage [  3x: 10x)"             = "cov_genome_pct_3x_10x_dragen",
        "PCT of genome with coverage [  1x:  3x)"             = "cov_genome_pct_1x_3x_dragen",
        "PCT of genome with coverage [  0x:  1x)"             = "cov_genome_pct_0x_1x_dragen",
        "Average chr X coverage over genome"                  = "cov_avg_x_over_genome_dragen",
        "Average chr Y coverage over genome"                  = "cov_avg_y_over_genome_dragen",
        "Average mitochondrial coverage over genome"          = "cov_avg_mt_over_genome_dragen",
        "Average autosomal coverage over genome"              = "cov_avg_auto_over_genome_dragen",
        "Median autosomal coverage over genome"               = "cov_median_auto_over_genome_dragen",
        "Mean/Median autosomal coverage ratio over genome"    = "cov_mean_median_auto_ratio_over_genome_dragen",
        "Aligned reads"                                       = "reads_aligned_dragen",
        "Aligned reads in genome"                             = "reads_aligned_in_genome_dragen"
      )

      x <- self$path
      raw <- readr::read_lines(x)
      assertthat::assert_that(grepl("COVERAGE SUMMARY", d[1]))

      raw |>
        tibble::as_tibble_col(column_name = "value") |>
        tidyr::separate_wider_delim(
          "value",
          delim = ",", too_few = "align_start",
          names = c("category", "dummy1", "var", "value", "pct")
        ) |>
        dplyr::mutate(
          var = dplyr::recode(.data$var, !!!abbrev_nm),
          value = as.numeric(.data$value)
        ) |>
        # pct just shows 100% for a couple rows
        dplyr::select("var", "value") |>
        tidyr::pivot_wider(names_from = "var", values_from = "value")
    },
    #' @description
    #' Writes a tidy version of the `wgs_coverage_metrics_<phenotype>.csv` file output
    #' from DRAGEN
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
      write_dracarys(obj = d, prefix = prefix, out_format = out_format, drid = drid)
    }
  )
)

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
#' d1 <- read(ch1)
#' d2 <- read(ch2)
#' ch1$plot(d1)
#' ch2$plot(d2)
#' ch1$write(d1, out_dir = tempdir(), prefix = "seqc_n", out_format = "tsv")
#' ch2$write(d2, out_dir = tempdir(), prefix = "seqc_t", out_format = "tsv")
#' @export
FineHistFile <- R6::R6Class(
  "FineHistFile",
  inherit = File,
  public = list(
    #' @description
    #' Reads the `wgs_fine_hist_<phenotype>.csv` file output from DRAGEN.
    #' @return tibble with three columns:
    #'   - label
    #'   - depth
    #'   - number of loci with given depth
    read = function() {
      x <- self$path
      d <- readr::read_csv(x, col_types = "cd")
      assertthat::assert_that(all(colnames(d) == c("Depth", "Overall")))

      # there's a max Depth of 2000+, so convert to numeric for easier plotting
      d |>
        dplyr::mutate(
          Depth = ifelse(grepl("+", .data$Depth), sub("(\\d*)\\+", "\\1", .data$Depth), .data$Depth),
          Depth = as.integer(.data$Depth)
        ) |>
        dplyr::select(depth = "Depth", n_loci = "Overall")
    },
    #' @description
    #' Writes a tidy version of the `wgs_fine_hist_<phenotype>.csv` file output
    #' from DRAGEN
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
      write_dracarys(obj = d, prefix = prefix, out_format = out_format, drid = drid)
    },

    #' @description Plots the `wgs_fine_hist_<phenotype>.csv` files.
    #' @param d Parsed object from `self$read()`.
    #' @param x_lim X axis range to plot.
    #' @return A ggplot2 object with depth of coverage on X axis,
    #' and number of loci with that depth on Y axis.
    plot = function(d, x_lim = c(0, 300)) {
      assertthat::assert_that(length(x_lim) == 2)
      d |>
        ggplot2::ggplot(ggplot2::aes(x = .data$depth, y = .data$n_loci)) +
        ggplot2::geom_line() +
        ggplot2::coord_cartesian(xlim = x_lim) +
        ggplot2::scale_y_continuous(labels = scales::label_comma()) +
        ggplot2::scale_x_continuous(n.breaks = 8) +
        ggplot2::labs(title = "Coverage Distribution", colour = "Label") +
        ggplot2::xlab("Depth of Coverage") +
        ggplot2::ylab("Number of Loci with Given Coverage") +
        ggplot2::theme_minimal() +
        ggplot2::theme(
          legend.position = c(0.9, 0.9),
          legend.justification = c(1, 1),
          panel.grid.minor = ggplot2::element_blank(),
          panel.grid.major.y = ggplot2::element_blank(),
          plot.title = ggplot2::element_text(colour = "#2c3e50", size = 14, face = "bold")
        )
    }
  )
)

#' FragmentLengthHistFile R6 Class
#'
#' @description
#' Contains methods for reading and plotting contents of
#' the `fragment_length_hist.csv` file output from DRAGEN.
#'
#' @examples
#' x <- system.file("extdata/wgs/SEQC-II.fragment_length_hist.csv.gz", package = "dracarys")
#' fl <- FragmentLengthHistFile$new(x)
#' d <- fl$read() # or read(fl)
#' fl$plot(d) # or plot(fl)
#' fl$write(d |> dplyr::filter(count > 10), out_dir = tempdir(), prefix = "seqc_fl")
#' @export
FragmentLengthHistFile <- R6::R6Class(
  "FragmentLengthHistFile",
  inherit = File,
  public = list(
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
        tidyr::fill("sample", .direction = "down") |>
        dplyr::filter(!grepl("#Sample: |FragmentLength,Count", .data$value)) |>
        tidyr::separate_wider_delim(cols = "value", names = c("fragmentLength", "count"), delim = ",") |>
        dplyr::mutate(
          count = as.numeric(.data$count),
          fragmentLength = as.numeric(.data$fragmentLength)
        ) |>
        dplyr::select("sample", "fragmentLength", "count")
    },
    #' @description
    #' Writes a tidy version of the `fragment_length_hist.csv` file output
    #' from DRAGEN.
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
      write_dracarys(obj = d, prefix = prefix, out_format = out_format, drid = drid)
    },


    #' @description Plots the fragment length distributions as given in the
    #' `fragment_length_hist.csv` file.
    #'
    #' @param d Parsed object from `self$read()`.
    #' @param min_count Minimum read count to be plotted (Default: 10).
    #' @return A ggplot2 plot containing fragment lengths on X axis and read counts
    #'   on Y axis for each sample.
    plot = function(d, min_count = 10) {
      assertthat::assert_that(min_count >= 0)
      d <- d |>
        dplyr::filter(.data$count >= min_count)

      d |>
        ggplot2::ggplot(ggplot2::aes(x = .data$fragmentLength, y = .data$count, colour = sample)) +
        ggplot2::geom_line() +
        ggplot2::labs(title = "Fragment Length Distribution") +
        ggplot2::xlab("Fragment Length (bp)") +
        ggplot2::ylab(glue("Read Count (min: {min_count})")) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
          legend.position = c(0.9, 0.9),
          legend.justification = c(1, 1),
          panel.grid.minor = ggplot2::element_blank(),
          panel.grid.major = ggplot2::element_blank(),
          plot.title = ggplot2::element_text(colour = "#2c3e50", size = 14, face = "bold")
        )
    }
  )
)

#' MappingMetricsFile R6 Class
#'
#' @description
#' Contains methods for reading and displaying contents of
#' the `mapping_metrics.csv` file output from DRAGEN.
#' This file contains mapping and aligning metrics, like the metrics computed by
#' the Samtools Flagstat command. These metrics are available on an aggregate
#' level (over all input data), and on a per read group level.
#' Unless explicitly stated, the metrics units are in reads (i.e., not in
#' terms of pairs or alignments).
#'
#' @examples
#' x <- system.file("extdata/wgs/SEQC-II.mapping_metrics.csv.gz", package = "dracarys")
#' mm <- MappingMetricsFile$new(x)
#' d <- mm$read() # or read(mm)
#' mm$write(d, out_dir = tempdir(), prefix = "seqc_mm", out_format = "tsv")
#'
#' @export
MappingMetricsFile <- R6::R6Class(
  "MappingMetricsFile",
  inherit = File,
  public = list(
    #' @description
    #' Reads the `mapping_metrics.csv` file output from DRAGEN.
    #'
    #' @return tibble with the following columns:
    #'     - category: summary or read group
    #'     - Phenotype: e.g. tumor, normal
    #'     - RG: read group
    #'     - var: metric variable
    #'     - var_abbrev: metric variable abbreviation
    #'     - count: count of reads
    #'     - pct: percentage of reads
    read = function() {
      abbrev_nm <- c(
        "Total Reads per RG" = "Tot",
        "Number of duplicate marked reads" = "Dup",
        "Number of duplicate marked and mate reads removed" = "Dup Rem",
        "Number of unique reads (excl. duplicate marked reads)" = "Unique",
        "Reads with mate sequenced" = "Mated",
        "Reads without mate sequenced" = "noMated",
        "QC-failed reads" = "Failed",
        "Mapped reads" = "Mapped",
        "Mapped reads R1" = "R1map",
        "Mapped reads R2" = "R2map",
        "Number of unique & mapped reads (excl. duplicate marked reads)" = "UniqueMap",
        "Unmapped reads" = "Unmapped",
        "Singleton reads (itself mapped; mate unmapped)" = "Singleton",
        "Paired reads (itself & mate mapped)" = "PairedMap",
        "Properly paired reads" = "PairedProper",
        "Not properly paired reads (discordant)" = "PairedDisc",
        "Paired reads mapped to different chromosomes" = "DiffChrom",
        "Paired reads mapped to different chromosomes (MAPQ>=10)" = "DiffChrom MQ10",
        "Reads with MAPQ [40:inf)" = "MQ 40+",
        "Reads with MAPQ [30:40)" = "MQ 30-40",
        "Reads with MAPQ [20:30)" = "MQ 20-30",
        "Reads with MAPQ [10:20)" = "MQ 10-20",
        "Reads with MAPQ [ 0:10)" = "MQ 0-10",
        "Reads with MAPQ NA (Unmapped reads)" = "MQ NA",
        "Reads with indel R1" = "R1 indel",
        "Reads with indel R2" = "R2 indel",
        "Total bases" = "TotalBases",
        "Total bases R1" = "TotBasesR1",
        "Total bases R2" = "TotBasesR2",
        "Mapped bases R1" = "MappedBasesR1",
        "Mapped bases R2" = "MappedBasesR2",
        "Soft-clipped bases R1" = "SoftClipR1",
        "Soft-clipped bases R2" = "SoftClipR2",
        "Mismatched bases R1" = "MismatchR1",
        "Mismatched bases R2" = "MismatchR2",
        "Mismatched bases R1 (excl. indels)" = "MismatchR1 NI",
        "Mismatched bases R2 (excl. indels)" = "MismatchR2 NI",
        "Q30 bases" = "Q30",
        "Q30 bases R1" = "R1Q30",
        "Q30 bases R2" = "R2Q30",
        "Q30 bases (excl. dups & clipped bases)" = "Q30 nondup",
        "Total alignments" = "TotAlign",
        "Secondary alignments" = "SecAlign",
        "Supplementary (chimeric) alignments" = "ChimericAlign",
        "Estimated read length" = "Read Length",
        "Bases in reference genome" = "Genome Bases",
        "Bases in target bed [% of genome]" = "BedBases %Genome",
        "Average sequenced coverage over genome" = "Coverage Avg",
        "Insert length: mean" = "InsertLength Mean",
        "Insert length: median" = "InsertLength Median",
        "Insert length: standard deviation" = "InsertLength StdDev",
        "Provided sex chromosome ploidy" = "Ploidy SexChrom",
        "DRAGEN mapping rate [mil. reads/second]" = "Map Rate"
      )

      x <- self$path
      d <- readr::read_lines(x)
      assertthat::assert_that(grepl("MAPPING/ALIGNING", d[1]))

      d |>
        tibble::as_tibble_col(column_name = "value") |>
        tidyr::separate_wider_delim("value", names = c("category", "RG", "extra"), delim = ",", too_many = "merge") |>
        tidyr::separate_wider_delim("extra", names = c("var", "value"), delim = ",", too_many = "merge") |>
        tidyr::separate_wider_delim("value", names = c("count", "pct"), delim = ",", too_few = "align_start") |>
        dplyr::mutate(
          count = dplyr::na_if(.data$count, "NA"),
          count = as.numeric(.data$count),
          pct = as.numeric(.data$pct),
          Phenotype = dplyr::case_when(
            grepl("TUMOR", .data$category) ~ "tumor",
            grepl("NORMAL", .data$category) ~ "normal",
            TRUE ~ "unknown"
          ),
          category = dplyr::case_when(
            grepl("ALIGNING SUMMARY", .data$category) ~ "summary",
            grepl("ALIGNING PER RG", .data$category) ~ "readgroup",
            TRUE ~ "unknown"
          ),
          RG = ifelse(.data$RG == "", "TOTAL", .data$RG),
          var = ifelse(grepl("Total.*reads", .data$var), "Total Reads per RG", .data$var),
          var_abbrev = dplyr::recode(.data$var, !!!abbrev_nm)
        ) |>
        dplyr::select(
          "category", "Phenotype", "RG",
          "var", "var_abbrev", "count", "pct"
        )
    },
    #' @description
    #' Writes a tidy version of the `mapping_metrics.csv` file output
    #' from DRAGEN.
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
      write_dracarys(obj = d, prefix = prefix, out_format = out_format, drid = drid)
    }
  )
)

#' PloidyEstimationMetricsFile R6 Class
#'
#' @description
#' Contains methods for reading contents of
#' the `ploidy_estimation_metrics.csv` file output from DRAGEN.
#'
#' @examples
#' x <- system.file("extdata/wgs/SEQC-II.ploidy_estimation_metrics.csv.gz", package = "dracarys")
#' pem <- PloidyEstimationMetricsFile$new(x)
#' d <- pem$read() # or read(pem)
#' pem$write(d, out_dir = tempdir(), prefix = "seqc_ploidy", out_format = "tsv")
#'
#' @export
PloidyEstimationMetricsFile <- R6::R6Class(
  "PloidyEstimationMetricsFile",
  inherit = File,
  public = list(
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
      d <- d |>
        tibble::as_tibble_col(column_name = "value") |>
        tidyr::separate_wider_delim("value", names = c("dummy1", "dummy2", "var", "value"), delim = ",") |>
        dplyr::select("var", "value") |>
        tidyr::pivot_wider(names_from = "var", values_from = "value")
      # now convert all except 'Ploidy estimation' to numeric
      cols1 <- colnames(d)[colnames(d) != "Ploidy estimation"]
      d |>
        dplyr::mutate(dplyr::across(dplyr::all_of(cols1), as.numeric))
    },
    #' @description
    #' Writes a tidy version of the `ploidy_estimation_metrics.csv` file output
    #' from DRAGEN.
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
      write_dracarys(obj = d, prefix = prefix, out_format = out_format, drid = drid)
    }
  )
)

#' ReplayFile R6 Class
#'
#' @description
#' Contains methods for reading contents of
#' the `replay.json` file output from DRAGEN, which contains the DRAGEN command
#' line, parameters and version for the specific run.
#'
#' @examples
#' x <- system.file("extdata/wgs/SEQC-II-replay.json.gz", package = "dracarys")
#' r <- ReplayFile$new(x)
#' d <- r$read() # or read(r)
#' r$write(d, out_dir = tempdir(), prefix = "seqc_replay", out_format = "tsv")
#' @export
ReplayFile <- R6::R6Class(
  "ReplayFile",
  inherit = File,
  public = list(
    #' @description Reads the `replay.json` file.
    #' @return A list with the following elements:
    #'   - `command_line`: character of DRAGEN command line used.
    #'   - `dragen_config`: tibble of parameters used for the DRAGEN run.
    #'   - `system`: tibble with dragen_version, nodename, and kernel_release.
    #'   - `label`: character of sample label (inferred from file name)
    #'   - `hash_table_build`: tibble with details about the DRAGEN hash table build.
    read = function() {
      x <- self$path
      res <- x |>
        jsonlite::read_json(simplifyVector = TRUE) |>
        purrr::map_if(is.data.frame, tibble::as_tibble)

      req_elements <- c("command_line", "hash_table_build", "dragen_config", "system")
      assertthat::assert_that(all(names(res) %in% req_elements))
      res[["system"]] <- res[["system"]] |>
        tibble::as_tibble_row()
      res[["hash_table_build"]] <- res[["hash_table_build"]] |>
        tibble::as_tibble_row()
      # we don't care if the columns are characters, no analysis likely to be done on dragen options
      # (though never say never!)
      res[["dragen_config"]] <- res[["dragen_config"]] |>
        tidyr::pivot_wider(names_from = "name", values_from = "value")

      dplyr::bind_cols(res)
    },
    #' @description
    #' Writes a tidy version of the `replay.json` file output
    #' from DRAGEN.
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
      write_dracarys(obj = d, prefix = prefix, out_format = out_format, drid = drid)
    }
  )
)

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
#' d <- tm$read() # or read(tm)
#' tm$write(d, out_dir = tempdir(), prefix = "seqc_time", out_format = "tsv")
#' @export
TimeMetricsFile <- R6::R6Class(
  "TimeMetricsFile",
  inherit = File,
  public = list(
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
      d |>
        dplyr::mutate(
          Step = tools::toTitleCase(sub("Time ", "", .data$Step)),
          Time = substr(.data$time_hrs, 1, 5)
        ) |>
        dplyr::select("Step", "Time") |>
        tidyr::pivot_wider(names_from = "Step", values_from = "Time") |>
        dplyr::relocate("Total Runtime")
    },
    #' @description
    #' Writes a tidy version of the `time_metrics.csv` file output
    #' from DRAGEN.
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
      write_dracarys(obj = d, prefix = prefix, out_format = out_format, drid = drid)
    }
  )
)

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
#' p <- system.file("extdata/wgs/SEQC-II.time_metrics.csv.gz", package = "dracarys")
#' x <- TimeMetricsFile$new(p)
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
    dplyr::bind_rows(.id = "ID")
}

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
#' d <- vm$read() # or read(vm)
#' vm$write(d, out_dir = tempdir(), prefix = "seqc_vc", out_format = "tsv")
#'
#' @export
VCMetricsFile <- R6::R6Class(
  "VCMetricsFile",
  inherit = File,
  public = list(
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
        tibble::as_tibble_col(column_name = "value") |>
        tidyr::separate_wider_delim("value", names = c("category", "sample", "extra"), delim = ",", too_many = "merge") |>
        tidyr::separate_wider_delim("extra", names = c("var", "value"), delim = ",", too_many = "merge") |>
        tidyr::separate_wider_delim("value", names = c("count", "pct"), delim = ",", too_few = "align_start") |>
        dplyr::mutate(
          count = dplyr::na_if(.data$count, "NA"),
          count = as.numeric(.data$count),
          pct = as.numeric(.data$pct),
          category = dplyr::case_when(
            grepl("SUMMARY", .data$category) ~ "summary",
            grepl("PREFILTER", .data$category) ~ "prefilter",
            grepl("POSTFILTER", .data$category) ~ "postfilter",
            TRUE ~ "unknown"
          )
        ) |>
        dplyr::select("category", "sample", "var", "count", "pct")
    },
    #' @description
    #' Writes a tidy version of the `vc_metrics.csv` file output
    #' from DRAGEN.
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
      write_dracarys(obj = d, prefix = prefix, out_format = out_format, drid = drid)
    }
  )
)

#' TrimmerMetricsFile R6 Class
#'
#' @description
#' Contains methods for reading and displaying contents of
#' the `trimmer_metrics.csv` file output from DRAGEN
#'
#' @examples
#' x <- system.file("extdata/wgs/SEQC-II.trimmer_metrics.csv.gz", package = "dracarys")
#' tm <- TrimmerMetricsFile$new(x)
#' d <- tm$read()
#' tm$write(d, out_dir = tempdir(), prefix = "seqc_tm", out_format = "tsv")
#'
#' @export
TrimmerMetricsFile <- R6::R6Class(
  "TrimmerMetricsFile",
  inherit = File,
  public = list(
    #' @description
    #' Reads the `trimmer_metrics.csv` file output from DRAGEN.
    #'
    #' @return tibble (TODO)
    read = function() {
      x <- self$path
      d <- readr::read_lines(x)
      assertthat::assert_that(grepl("TRIMMER STATISTICS", d[1]))
      trimmer_abbrev_nm <- c(
        "Total input reads"                              = "reads_tot_input_dragen",
        "Total input bases"                              = "bases_tot_dragen",
        "Total input bases R1"                           = "bases_r1_dragen",
        "Total input bases R2"                           = "bases_r2_dragen",
        "Average input read length"                      = "read_len_avg_dragen",
        "Total trimmed reads"                            = "reads_trimmed_tot_dragen",
        "Total trimmed bases"                            = "bases_trimmed_tot_dragen",
        "Average bases trimmed per read"                 = "bases_trimmed_avg_per_read_dragen",
        "Average bases trimmed per trimmed read"         = "bases_trimmed_avg_per_trimmedread_dragen",
        "Remaining poly-G K-mers R1 3prime"              = "polygkmers3r1_remaining_dragen",
        "Remaining poly-G K-mers R2 3prime"              = "polygkmers3r2_remaining_dragen",
        "Poly-G soft trimmed reads unfiltered R1 3prime" = "polyg_soft_trimmed_reads_unfilt_3r1_dragen",
        "Poly-G soft trimmed reads unfiltered R2 3prime" = "polyg_soft_trimmed_reads_unfilt_3r2_dragen",
        "Poly-G soft trimmed reads filtered R1 3prime"   = "polyg_soft_trimmed_reads_filt_3r1_dragen",
        "Poly-G soft trimmed reads filtered R2 3prime"   = "polyg_soft_trimmed_reads_filt_3r2_dragen",
        "Poly-G soft trimmed bases unfiltered R1 3prime" = "polyg_soft_trimmed_bases_unfilt_3r1_dragen",
        "Poly-G soft trimmed bases unfiltered R2 3prime" = "polyg_soft_trimmed_bases_unfilt_3r2_dragen",
        "Poly-G soft trimmed bases filtered R1 3prime"   = "polyg_soft_trimmed_bases_filt_3r1_dragen",
        "Poly-G soft trimmed bases filtered R2 3prime"   = "polyg_soft_trimmed_bases_filt_3r2_dragen",
        "Total filtered reads"                           = "reads_tot_filt_dragen",
        "Reads filtered for minimum read length R1"      = "reads_filt_minreadlenr1_dragen",
        "Reads filtered for minimum read length R2"      = "reads_filt_minreadlenr2_dragen"
      )

      d |>
        tibble::as_tibble_col(column_name = "value") |>
        tidyr::separate_wider_delim("value", names = c("category", "extra", "var", "count", "pct"), delim = ",", too_few = "align_start") |>
        dplyr::mutate(
          count = as.numeric(.data$count),
          pct = round(as.numeric(.data$pct), 2),
          var = dplyr::recode(.data$var, !!!trimmer_abbrev_nm)
        ) |>
        dplyr::select("var", "count", "pct") |>
        tidyr::pivot_longer(c("count", "pct")) |>
        dplyr::filter(!is.na(.data$value)) |>
        dplyr::mutate(
          name = dplyr::if_else(.data$name == "count", "", "_pct"),
          var = glue("{.data$var}{.data$name}")
        ) |>
        dplyr::select("var", "value") |>
        tidyr::pivot_wider(names_from = "var", values_from = "value")
    },
    #' @description
    #' Writes a tidy version of the `trimmer_metrics.csv` file output
    #' from DRAGEN.
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
      write_dracarys(obj = d, prefix = prefix, out_format = out_format, drid = drid)
    }
  )
)

# TODO:
# - wgs_hist.csv
#
# wts
# - quant_metrics.csv
# - quant.transcript_fragment_lengths.txt
# - quant.transcript_coverage.txt
# - saturation.txt
