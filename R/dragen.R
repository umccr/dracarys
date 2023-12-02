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
#' cc1$plot(d1)
#' cc2$plot(d2)
#'
#' @export
ContigMeanCovFile <- R6::R6Class("ContigMeanCovFile", inherit = File, public = list(
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
    b <- self$bname()
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
))

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
#' read(cm1)
#' read(cm2)
#'
#' @export
CoverageMetricsFile <- R6::R6Class("CoverageMetricsFile", inherit = File, public = list(
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
      "Aligned bases" = "Aln bases",
      "Aligned bases in genome" = "Aln Bases Genome",
      "Average alignment coverage over genome" = "Avg Cov Genome",
      "Average chr X coverage over genome" = "Avg Cov chrX",
      "Average chr Y coverage over genome" = "Avg Cov chrY",
      "Average mitochondrial coverage over genome" = "Avg Cov chrM",
      "Average autosomal coverage over genome" = "Avg Cov Autos",
      "Median autosomal coverage over genome" = "Med Cov Autos",
      "Mean/Median autosomal coverage ratio over genome" = "Cov Ratio",
      "Aligned reads" = "Aln Reads",
      "Aligned reads in genome" = "Aln Reads Genome"
    )

    x <- self$path
    b <- self$bname()
    suffix <- dplyr::if_else(
      grepl("_normal\\.csv", b), "_N",
      dplyr::if_else(grepl("_tumor\\.csv", b), "_T", "")
    )
    nm <- sub("(.*)\\.wgs_coverage_metrics.*", "\\1", b)
    label <- paste0(nm, suffix)

    d <- readr::read_lines(x)
    assertthat::assert_that(grepl("COVERAGE SUMMARY", d[1]))

    d <- d |>
      tibble::enframe(name = "name", value = "value") |>
      tidyr::separate(.data$value, into = c("category", "dummy1", "extra"), sep = ",", extra = "merge") |>
      tidyr::separate(.data$extra, into = c("var", "value"), sep = ",", extra = "merge") |>
      dplyr::mutate(label = tidyselect::all_of(label))

    pct <- d |>
      dplyr::filter(grepl("PCT", .data$var)) |>
      dplyr::mutate(
        value = as.numeric(.data$value),
        var_abbrev = dplyr::case_when(
          grepl("PCT of genome", .data$var) ~ sub("PCT of genome with coverage", "%genome", var),
          grepl("Uniformity", .data$var) ~ "uniformity (% > 0.2*mean)",
          TRUE ~ "FOO"
        )
      ) |>
      dplyr::select("label", "var", "var_abbrev", pct = "value")

    cnt <- d |>
      dplyr::filter(!grepl("PCT", .data$var)) |>
      tidyr::separate(.data$value, into = c("count", "pct"), sep = ",", fill = "right", convert = TRUE) |>
      dplyr::mutate(var_abbrev = dplyr::recode(.data$var, !!!abbrev_nm)) |>
      dplyr::select("label", "var", "var_abbrev", "count", "pct")

    dplyr::bind_rows(pct, cnt)
  }
))

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
      dplyr::select("label", depth = "Depth", n_loci = "Overall")
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
      dplyr::select(-"name")
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
))

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
#' mm$read() # or read(mm)
#'
#' @export
MappingMetricsFile <- R6::R6Class("MappingMetricsFile", inherit = File, public = list(
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
      tibble::enframe(name = "name", value = "value") |>
      tidyr::separate(.data$value, into = c("category", "RG", "extra"), sep = ",", extra = "merge") |>
      tidyr::separate(.data$extra, into = c("var", "value"), sep = ",", extra = "merge") |>
      tidyr::separate(.data$value, into = c("count", "pct"), sep = ",", fill = "right", convert = TRUE) |>
      dplyr::mutate(
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
  }
))

#' PloidyEstimationMetricsFile R6 Class
#'
#' @description
#' Contains methods for reading contents of
#' the `ploidy_estimation_metrics.csv` file output from DRAGEN.
#'
#' @examples
#' x <- system.file("extdata/wgs/SEQC-II.ploidy_estimation_metrics.csv.gz", package = "dracarys")
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
      dplyr::select("label", "var", "value")
  }
))

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
#' r$read() # or read(r)
#' @export
ReplayFile <- R6::R6Class("ReplayFile", inherit = File, public = list(
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
      tibble::as_tibble() |>
      tidyr::pivot_longer(dplyr::everything())
    res[["hash_table_build"]] <- res[["hash_table_build"]] |>
      tibble::as_tibble() |>
      tidyr::pivot_longer(dplyr::everything())
    res[["label"]] <- sub("-replay.json.*", "", basename(x))

    res
  }
))

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
      dplyr::select("Label", "Step", "Time")
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
    dplyr::relocate("Total Runtime", .after = "Label")
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
      dplyr::select("category", "sample", "var", "count", "pct")
  }
))
