#' Tidy TSO ctDNA Results
#'
#' Tidies TSO500 ctDNA results into a list of tibbles and writes individual tibbles to
#' TSV and/or Parquet format.
#'
#' @param in_dir Directory path to TSO500 ctDNA workflow results (can be GDS or local).
#' @param prefix Prefix of output file(s).
#' @param out_dir Output directory.
#' @param gds_local_dir If `indir` is a GDS directory, 'recognisable' files
#' will be first downloaded to this directory.
#' @param dryrun Just list the files that will be downloaded (def: FALSE).
#' @param token ICA access token (by default uses $ICA_ACCESS_TOKEN env var).
#' @param out_format Format of output (tsv, parquet, both) (def: tsv).
#'
#' @return Tibble with path to input file and the resultant tidy object.
#' @examples
#' \dontrun{
#' in_dir <- paste0(
#'   "gds://production/analysis_data/SBJ02858/tso_ctdna_tumor_only/",
#'   "20221104b7ad0b38/L2201560/Results/PRJ222206_L2201560/"
#' )
#' in_dir <- here::here(glue("nogit/tso/2022-12-13/SBJ02858/dracarys_gds_sync"))
#' out_dir <- file.path(in_dir, "../out")
#' gds_local_dir <- NULL
#' prefix <- "SBJ02858"
#' dryrun <- F
#' tso_tidy(in_dir = in_dir, out_dir = out_dir, prefix = prefix)
#' }
#' @export
tso_tidy <- function(in_dir, out_dir, prefix, gds_local_dir = NULL, out_format = "tsv",
                     dryrun = FALSE, token = Sys.getenv("ICA_ACCESS_TOKEN")) {
  output_format_valid(out_format)
  pat <- "tso__"
  e <- emojifont::emoji

  if (grepl("^gds://", in_dir)) {
    if (is.null(gds_local_dir)) {
      gds_local_dir <- file.path(out_dir, "dracarys_gds_sync")
    }
    dr_gds_download(
      gdsdir = in_dir, outdir = gds_local_dir, token = token,
      pattern = pat, dryrun = dryrun
    )
    # Use the downloaded results
    in_dir <- gds_local_dir
  } else {
    # in_dir is not gds
    if (!is.null(gds_local_dir)) {
      cli::cli_warn(glue(
        "You have specified 'gds_local_dir' to download GDS results,\n",
        "but your input directory is local - ignoring {gds_local_dir}"
      ))
    }
  }
  if (dryrun) {
    cli::cli_inform("You have specified 'dryrun' - terminating {e('ghost')}!")
    return(NULL)
  } else {
    d <- fs::dir_ls(in_dir) |>
      tibble::as_tibble_col(column_name = "path") |>
      dplyr::mutate(
        bname = basename(.data$path),
        type = purrr::map_chr(.data$bname, match_regex)
      ) |>
      dplyr::filter(!is.na(.data$type), grepl(pat, .data$type))

    if (nrow(d) == 0) {
      regex <- FILE_REGEX |>
        dplyr::filter(grepl("tso__", .data$name)) |>
        dplyr::pull("regex")
      msg <- paste(
        "No TSO files for dracarys were found in {.file {in_dir}}.",
        "See current supported regexes for TSO: {regex}."
      )
      cli::cli_abort(msg)
    }

    cli::cli_alert_info("{date_log()} {e('dragon')} {.emph {prefix}}: Start tidying TSO dir: {.file {in_dir}}")
    res <- d |>
      dplyr::select("type", "path", "bname") |>
      dplyr::rowwise() |>
      dplyr::mutate(
        env = list(tso_funcall(.data$type)),
        obj = list(.data$env$new(.data$path)),
        has_plot = "plot" %in% names(.data$env[["public_methods"]]),
        obj_parsed = list(.data$obj$read()),
        obj_parsed2 = list(.data$obj$write(.data$obj_parsed, out_dir = out_dir, prefix = prefix, out_format = out_format)),
        plot = ifelse(.data$has_plot, list(.data$obj$plot(.data$obj_parsed)), list(NULL))
      ) |>
      dplyr::select("type", "path", "obj", dat = "obj_parsed", "plot")
    cli::cli_alert_success("{date_log()} {e('tada')} {.emph {prefix}}: TSO tidy results at: {.file {out_dir}}")
    return(invisible(res))
  }
}

tso_funcall <- function(x) {
  l <- list(
    "tso__align_collapse_fusion_caller_metrics" = TsoAlignCollapseFusionCallerMetricsFile,
    "tso__target_region_coverage" = TsoTargetRegionCoverageFile,
    "tso__fragment_length_hist" = TsoFragmentLengthHistFile,
    "tso__msi" = TsoMsiFile,
    "tso__tmb" = TsoTmbFile,
    "tso__tmb_trace_tsv" = TsoTmbTraceTsvFile,
    "tso__fusions_csv" = TsoFusionsCsvFile,
    "tso__sample_analysis_results" = TsoSampleAnalysisResultsFile
  )
  if (!x %in% names(l)) {
    return(NULL)
  }
  l[[x]]
}

#' TsoTmbTraceTsvFile R6 Class
#'
#' @description
#' Contains methods for reading and displaying contents of the
#' `TMB_Trace.tsv` file output from TSO.
#'
#' @examples
#' x <- system.file("extdata/tso/sample705_TMB_Trace.tsv", package = "dracarys")
#' d <- TsoTmbTraceTsvFile$new(x)
#' d_parsed <- d$read() # or read(d)
#' d$write(d_parsed, out_dir = tempdir(), prefix = "sample705", out_format = "both")
#' @export
TsoTmbTraceTsvFile <- R6::R6Class(
  "TsoTmbTraceTsvFile",
  inherit = File,
  public = list(
    #' @description
    #' Reads the `TMB_Trace.tsv` file output from TSO.
    #'
    #' @return tibble with the following columns:
    #' * FragmentLength
    #' * Count
    read = function() {
      x <- self$path
      ct <- readr::cols(
        Chromosome = "c", Position = "i", RefCall = "c", AltCall = "c",
        VAF = "d", Depth = "d", CytoBand = "c", GeneName = "c",
        VariantType = "c", CosmicIDs = "c", MaxCosmicCount = "d",
        AlleleCountsGnomadExome = "d", AlleleCountsGnomadGenome = "d",
        AlleleCounts1000Genomes = "d", MaxDatabaseAlleleCounts = "d",
        GermlineFilterDatabase = "l", GermlineFilterProxi = "l",
        CodingVariant = "l", Nonsynonymous = "l", IncludedInTMBNumerator = "l"
      )
      readr::read_tsv(x, col_types = ct)
    },

    #' @description
    #' Writes a tidy version of the `TMB_Trace.tsv` file output from TSO.
    #'
    #' @param d Parsed object from `self$read()`.
    #' @param prefix Prefix of output file(s).
    #' @param out_dir Output directory.
    #' @param out_format Format of output file(s) (one of 'tsv' (def.),
    #' 'parquet', 'both').
    write = function(d, out_dir, prefix, out_format = "tsv") {
      prefix <- file.path(out_dir, prefix)
      prefix2 <- glue("{prefix}_TMB_Trace")
      write_dracarys(obj = d, prefix = prefix2, out_format = out_format)
    }
  )
)

#' TsoFragmentLengthHistFile R6 Class
#'
#' @description
#' Contains methods for reading and displaying contents of the
#' `fragment_length_hist.json.gz` file output from TSO.
#'
#' @examples
#' x <- system.file("extdata/tso/sample705.fragment_length_hist.json.gz", package = "dracarys")
#' fl <- TsoFragmentLengthHistFile$new(x)
#' d_parsed <- fl$read() # or read(fl)
#' fl$plot(d_parsed, 5)
#' fl$write(d_parsed, out_dir = tempdir(), prefix = "sample705", out_format = "both")
#' @export
TsoFragmentLengthHistFile <- R6::R6Class(
  "TsoFragmentLengthHistFile",
  inherit = File,
  public = list(
    #' @description
    #' Reads the `fragment_length_hist.json.gz` file output from TSO.
    #'
    #' @return tibble with the following columns:
    #' * FragmentLength
    #' * Count
    read = function() {
      x <- self$path
      j <- jsonlite::read_json(x)
      assertthat::assert_that(
        all(names(j[[1]] %in% c("FragmentLength", "Count")))
      )
      j |>
        purrr::map(tibble::as_tibble) |>
        dplyr::bind_rows()
    },

    #' @description
    #' Writes a tidy version of the `fragment_length_hist.json.gz` file output
    #' from TSO.
    #'
    #' @param d Parsed object from `self$read()`.
    #' @param prefix Prefix of output file(s).
    #' @param out_dir Output directory.
    #' @param out_format Format of output file(s) (one of 'tsv' (def.),
    #' 'parquet', 'both').
    #'
    write = function(d, out_dir, prefix, out_format = "tsv") {
      prefix <- file.path(out_dir, prefix)
      prefix2 <- glue("{prefix}_fragment_length_hist")
      write_dracarys(obj = d, prefix = prefix2, out_format = out_format)
    },


    #' @description Plots the fragment length distributions as given in the
    #' `fragment_length_hist.json.gz` file.
    #'
    #' @param d Parsed object from `self$read()`.
    #' @param min_count Minimum read count to be plotted (Default: 10).
    #' @return A ggplot2 plot containing fragment lengths on X axis and read counts
    #'   on Y axis for each sample.
    plot = function(d, min_count = 10) {
      assertthat::assert_that(is.numeric(min_count), min_count >= 0)
      d |>
        dplyr::filter(.data$Count >= min_count) |>
        ggplot2::ggplot(ggplot2::aes(x = .data$FragmentLength, y = .data$Count)) +
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

#' TsoTargetRegionCoverageFile R6 Class
#'
#' @description
#' Contains methods for reading and displaying contents of the
#' `TargetRegionCoverage.json.gz` file output from TSO.
#'
#' @examples
#' x <- system.file("extdata/tso/sample705.TargetRegionCoverage.json.gz", package = "dracarys")
#' trc <- TsoTargetRegionCoverageFile$new(x)
#' d_parsed <- trc$read() # or read(trc)
#' trc$plot(d_parsed, 90) # or plot(trc, d_parsed, 90)
#' trc$write(d_parsed, out_dir = tempdir(), prefix = "sample705", out_format = "both")
#' @export
TsoTargetRegionCoverageFile <- R6::R6Class(
  "TsoTargetRegionCoverageFile",
  inherit = File,
  public = list(
    #' @description
    #' Reads the `TargetRegionCoverage.json.gz` file output from TSO.
    #'
    #' @return tibble with the following columns:
    #' * ConsensusReadDepth
    #' * BasePair
    #' * Percentage
    read = function() {
      x <- self$path
      l2tib <- function(l) {
        tibble::as_tibble(l) |>
          dplyr::mutate(dplyr::across(dplyr::everything(), ~ as.character(.)))
      }
      j <- jsonlite::read_json(x)
      assertthat::assert_that(
        all(names(j[[1]] %in% c("ConsensusReadDepth", "BasePair", "Percentage")))
      )
      j |>
        purrr::map(l2tib) |>
        dplyr::bind_rows() |>
        dplyr::mutate(Percentage = as.numeric(sub("%", "", .data$Percentage)))
    },

    #' @description
    #' Writes a tidy version of the `TargetRegionCoverage.json.gz` file output
    #' from TSO.
    #'
    #' @param d Parsed object from `self$read()`.
    #' @param prefix Prefix of output file(s).
    #' @param out_dir Output directory.
    #' @param out_format Format of output file(s) (one of 'tsv' (def.),
    #' 'parquet', 'both').
    #'
    write = function(d, out_dir, prefix, out_format = "tsv") {
      prefix <- file.path(out_dir, prefix)
      prefix2 <- glue("{prefix}_TargetRegionCoverage")
      write_dracarys(obj = d, prefix = prefix2, out_format = out_format)
    },

    #' @description Plots the `TargetRegionCoverage.json.gz` file.
    #'
    #' @param d Parsed object from `self$read()`.
    #' @param min_pct Minimum percentage to be plotted (Default: 2).
    #' @importFrom ggplot2 ggplot
    #' @return A ggplot2 plot containing read depth on X axis and percentage
    #'   covered on Y axis.
    plot = function(d, min_pct = 2) {
      assertthat::assert_that(is.numeric(min_pct), min_pct >= 0)
      d <- d |>
        dplyr::filter(
          !.data$ConsensusReadDepth == "TargetRegion",
          .data$Percentage >= min_pct
        ) |>
        dplyr::select(dp = "ConsensusReadDepth", pct = "Percentage") |>
        dplyr::mutate(dp = as.numeric(sub("X", "", .data$dp)))
      d |>
        ggplot2::ggplot(ggplot2::aes(x = dp, y = pct, label = dp)) +
        ggplot2::geom_point() +
        ggplot2::geom_line() +
        ggrepel::geom_text_repel() +
        ggplot2::labs(title = "Percentage of Target Region Covered by Given Read Depth") +
        ggplot2::xlab("Depth") +
        ggplot2::ylab("Percentage") +
        ggplot2::theme_minimal() +
        ggplot2::theme(
          legend.position = c(0.9, 0.9),
          legend.justification = c(1, 1),
          panel.grid.minor = ggplot2::element_blank(),
          plot.title = ggplot2::element_text(colour = "#2c3e50", size = 14, face = "bold")
        )
    }
  )
)

#' TsoAlignCollapseFusionCallerMetricsFile R6 Class
#'
#' @description
#' Contains methods for reading and displaying contents of the
#' `AlignCollapseFusionCaller_metrics.json.gz` file output from TSO.
#'
#' @examples
#' x <- system.file("extdata/tso/sample705.AlignCollapseFusionCaller_metrics.json.gz",
#'   package = "dracarys"
#' )
#' m <- TsoAlignCollapseFusionCallerMetricsFile$new(x)
#' d_parsed <- m$read() # or read(m)
#' m$write(d_parsed, out_dir = tempdir(), prefix = "sample705", out_format = "both")
#' @export
TsoAlignCollapseFusionCallerMetricsFile <- R6::R6Class(
  "TsoAlignCollapseFusionCallerMetricsFile",
  inherit = File,
  public = list(
    #' @description
    #' Reads the `AlignCollapseFusionCaller_metrics.json.gz` file output from TSO.
    #'
    #' @return tibble with the following columns:
    #' * section: name of original JSON element
    #' * name: name of metric
    #' * value: value of metric
    #' * percent: percentage
    read = function() {
      x <- self$path
      l2tib <- function(el) {
        # list to tibble and turn cols to char
        fun1 <- function(l) {
          # handle silly NULLs..
          if (l[["value"]] |> is.null()) {
            l[["value"]] <- NA
          }
          tibble::as_tibble(l) |>
            dplyr::mutate(dplyr::across(dplyr::everything(), ~ as.character(.)))
        }
        el |>
          purrr::map(fun1) |>
          dplyr::bind_rows()
      }
      # j <- jsonlite::read_json(x) # cannot handle Infinity
      j <- RJSONIO::fromJSON(x, simplify = FALSE) # turns Infinity to NULL
      j |>
        purrr::map(l2tib) |>
        dplyr::bind_rows(.id = "section")
    },

    #' @description
    #' Writes a tidy version of the `AlignCollapseFusionCaller_metrics.json.gz` file
    #' output from TSO.
    #'
    #' Histo is the majority from UmiStatistics section, write out separately.
    #' Histo of num supporting fragments: Num of families with 0/1/2/3... raw reads.
    #' Histo of unique UMIs per fragment pos: Num of pos with 0/1/2/3... UMI seqs.
    #'
    #' @param d Parsed object from `self$read()`.
    #' @param prefix Prefix of output file(s).
    #' @param out_dir Output directory.
    #' @param out_format Format of output file(s) (one of 'tsv' (def.),
    #' 'parquet', 'both').
    #'
    write = function(d, out_dir, prefix, out_format = "tsv") {
      prefix <- file.path(out_dir, prefix)
      dhist <- self$histoprep(d)
      dmain <- d |>
        dplyr::filter(!grepl("Hist", .data$name))

      p <- glue("{prefix}_AlignCollapseFusionCaller_metrics_")
      p_hist <- glue("{p}_hist")
      p_main <- glue("{p}_main")
      write_dracarys(obj = dhist, prefix = p_hist, out_format = out_format)
      write_dracarys(obj = dmain, prefix = p_main, out_format = out_format)
    },
    #' @description
    #' Prepares the UmiStatistics histogram data from the
    #' `AlignCollapseFusionCaller_metrics.json.gz` file output from TSO.
    #'
    #' - Histo is the majority from UmiStatistics section, deal with it separately.
    #' - Histo of num supporting fragments: Num of families with 0/1/2/3... raw reads.
    #' - Histo of unique UMIs per fragment pos: Num of pos with 0/1/2/3... UMI seqs.
    #'
    #' @param d Parsed object from `self$read()`.
    histoprep = function(d) {
      dhist <- d |>
        dplyr::filter(grepl("Hist", .data$name))
      assertthat::assert_that(all(dhist[["section"]] == "UmiStatistics"))
      assertthat::assert_that(all(dhist[["percent"]] %in% NA))
      dhist |>
        dplyr::mutate(
          name = sub("Histogram of ", "", .data$name),
          name = gsub(" ", "_", .data$name),
          value = as.numeric(.data$value)
        ) |>
        dplyr::group_by(name) |>
        dplyr::mutate(num = dplyr::row_number()) |>
        dplyr::ungroup() |>
        dplyr::select(c("name", "num", "value"))
    },
    #' @description
    #' Generates the UmiStatistics Histogram plots from the
    #' `AlignCollapseFusionCaller_metrics.json.gz` file output from TSO.
    #'
    #' Histo is the majority from UmiStatistics section, deal with it separately.
    #' Histo of num supporting fragments: Num of families with 0/1/2/3... raw reads.
    #' Histo of unique UMIs per fragment pos: Num of pos with 0/1/2/3... UMI seqs.
    #' @param d Parsed object from `self$read()`.
    #' @param max_num Maximum number to display in both plots.
    #' @return Both histogram plot objects.
    plot = function(d, max_num = 15) {
      h <- self$histoprep(d)
      # 15 seems like a good cutoff for both plots
      p1 <- h |>
        dplyr::filter(
          .data$name == "num_supporting_fragments",
          .data$num <= max_num
        ) |>
        ggplot2::ggplot(ggplot2::aes(x = num, y = value)) +
        ggplot2::geom_line() +
        ggplot2::theme_bw() +
        ggplot2::labs(title = "Number of families with 0/1/2/3... raw reads.") +
        ggplot2::xlab("Families") +
        ggplot2::ylab("Reads")
      p2 <- h |>
        dplyr::filter(
          name == "unique_UMIs_per_fragment_position",
          .data$num <= max_num
        ) |>
        ggplot2::ggplot(ggplot2::aes(x = num, y = value)) +
        ggplot2::geom_line() +
        ggplot2::theme_bw() +
        ggplot2::labs(title = "Number of positions with 0/1/2/3... UMI sequences.") +
        ggplot2::xlab("Positions") +
        ggplot2::ylab("UMI Sequences")
      list(
        p_num_supporting_fragments = p1,
        p_unique_umis_per_frag_pos = p2
      )
    }
  )
)

#' TsoTmbFile R6 Class
#'
#' @description
#' Contains methods for reading and displaying contents of the
#' `tmb.json.gz` file output from TSO.
#'
#' @examples
#' x <- system.file("extdata/tso/sample705.tmb.json.gz", package = "dracarys")
#' tmb <- TsoTmbFile$new(x)
#' d_parsed <- tmb$read() # or read(tmb)
#' tmb$write(d_parsed, out_dir = tempdir(), prefix = "sample705", out_format = "both")
#' @export
TsoTmbFile <- R6::R6Class(
  "TsoTmbFile",
  inherit = File,
  public = list(
    #' @description
    #' Reads the `tmb.json.gz` file output from TSO.
    #'
    #' @return tibble with the following columns:
    #' * TmbPerMb
    #' * AdjustedTmbPerMb
    #' * NonsynonymousTmbPerMb
    #' * AdjustedNonsynonymousTmbPerMb
    #' * SomaticCodingVariantsCount
    #' * NonsynonymousSomaticCodingVariantsCount
    #' * TotalRegionSizeMb
    #' * CodingRegionSizeMb
    read = function() {
      x <- self$path
      # j <- jsonlite::read_json(x) # cannot handle NaN
      j <- RJSONIO::fromJSON(x) # turns NaN to NULL
      # not interested in Settings element
      j[["Settings"]] <- NULL
      # handle silly NULLs
      j <- lapply(j, function(x) ifelse(is.null(x), NA, x))
      tibble::as_tibble_row(j) |>
        dplyr::mutate(dplyr::across(dplyr::everything(), as.numeric))
    },

    #' @description
    #' Writes a tidy version of the `tmb.json.gz` file output from TSO.
    #'
    #' @param d Parsed object from `self$read()`.
    #' @param prefix Prefix of output file(s).
    #' @param out_dir Output directory.
    #' @param out_format Format of output file(s) (one of 'tsv' (def.),
    #' 'parquet', 'both').
    #'
    write = function(d, out_dir, prefix, out_format = "tsv") {
      prefix <- file.path(out_dir, prefix)
      prefix2 <- glue("{prefix}_tmb")
      write_dracarys(obj = d, prefix = prefix2, out_format = out_format)
    }
  )
)

#' TsoFusionsCsvFile R6 Class
#'
#' @description
#' Contains methods for reading and displaying contents of the
#' `Fusions.csv` file output from TSO.
#'
#' @examples
#' x <- system.file("extdata/tso/sample705_Fusions.csv", package = "dracarys")
#' fus <- TsoFusionsCsvFile$new(x)
#' d_parsed <- fus$read() # or read(fus)
#' fus$write(d_parsed, out_dir = tempdir(), prefix = "sample705", out_format = "both")
#' @export
TsoFusionsCsvFile <- R6::R6Class(
  "TsoFusionsCsvFile",
  inherit = File,
  public = list(
    #' @description
    #' Reads the `Fusions.csv` file output from TSO.
    #'
    #' @return tibble with several columns.
    read = function() {
      x <- self$path
      ct <- readr::cols(
        Sample = "c", Name = "c", Chr1 = "c", Pos1 = "d", Chr2 = "c",
        Pos2 = "d", Direction = "c", Alt_Depth = "d", BP1_Depth = "d",
        BP2_Depth = "d", Total_Depth = "d", VAF = "d", Gene1 = "c", Gene2 = "c",
        Contig = "c", Filter = "c", Is_Cosmic_GenePair = "l"
      )
      readr::read_csv(x, col_types = ct, comment = "#")
    },

    #' @description
    #' Writes a tidy version of the `Fusions.csv` file output from TSO.
    #'
    #' @param d Parsed object from `self$read()`.
    #' @param prefix Prefix of output file(s).
    #' @param out_dir Output directory.
    #' @param out_format Format of output file(s) (one of 'tsv' (def.),
    #' 'parquet', 'both').
    write = function(d, out_dir, prefix, out_format = "tsv") {
      prefix <- file.path(out_dir, prefix)
      prefix2 <- glue("{prefix}_Fusions")
      write_dracarys(obj = d, prefix = prefix2, out_format = out_format)
    }
  )
)

#' TsoMsiFile R6 Class
#'
#' @description
#' Contains methods for reading and displaying contents of the
#' `msi.json.gz` file output from TSO.
#'
#' @examples
#' x <- system.file("extdata/tso/sample705.msi.json.gz", package = "dracarys")
#' msi <- TsoMsiFile$new(x)
#' d_parsed <- msi$read() # or read(msi)
#' msi$write(d_parsed, out_dir = tempdir(), prefix = "sample705", out_format = "both")
#' @export
TsoMsiFile <- R6::R6Class(
  "TsoMsiFile",
  inherit = File,
  public = list(
    #' @description
    #' Reads the `msi.json.gz` file output from TSO.
    #'
    #' @return tibble with the following columns:
    #'   - label:
    read = function() {
      x <- self$path
      j <- jsonlite::read_json(x)
      # not interested in Settings element
      j[["Settings"]] <- NULL
      if (is.null(j[["ResultMessage"]])) {
        j[["ResultMessage"]] <- NA_character_
      }
      if (j[["PercentageUnstableSites"]] == "NaN") {
        j[["PercentageUnstableSites"]] <- NA_real_
      }
      tibble::as_tibble_row(j) |>
        dplyr::mutate(ResultIsValid = as.character(.data$ResultIsValid))
    },

    #' @description
    #' Writes a tidy version of the `msi.json.gz` file output from TSO.
    #'
    #' @param d Parsed object from `self$read()`.
    #' @param prefix Prefix of output file(s).
    #' @param out_dir Output directory.
    #' @param out_format Format of output file(s) (one of 'tsv' (def.),
    #' 'parquet', 'both').
    #'
    write = function(d, out_dir, prefix, out_format = "tsv") {
      prefix <- file.path(out_dir, prefix)
      prefix2 <- glue("{prefix}_msi")
      write_dracarys(obj = d, prefix = prefix2, out_format = out_format)
    }
  )
)

#' TsoSampleAnalysisResultsFile R6 Class
#'
#' @description
#' Contains methods for reading and displaying contents of the
#' `SampleAnalysisResults.json.gz` file output from TSO.
#'
#' @examples
#' x <- system.file("extdata/tso/sample705_SampleAnalysisResults.json.gz", package = "dracarys")
#' res <- TsoSampleAnalysisResultsFile$new(x)
#' d_parsed <- res$read() # or read(res)
#' res$write(d_parsed, tempfile(), "both")
#' @export
TsoSampleAnalysisResultsFile <- R6::R6Class(
  "TsoSampleAnalysisResultsFile",
  inherit = File,
  public = list(
    #' @description
    #' Reads the `SampleAnalysisResults.json.gz` file output from TSO.
    #'
    #' @return list of tibbles
    read = function() {
      x <- self$path
      j <- jsonlite::read_json(x)
      dat <- j[["data"]]
      ## sampleInformation
      sample_info <- dat[["sampleInformation"]] |>
        tibble::as_tibble_row()
      ## softwareConfiguration
      sw_conf <- dat[["softwareConfiguration"]]
      sw_nl_data_sources <- sw_conf[["nirvanaVersionList"]][[1]][["dataSources"]] |>
        purrr::map(tibble::as_tibble_row) |>
        dplyr::bind_rows()
      # get rid of it to grab the remaining elements
      sw_conf[["nirvanaVersionList"]][[1]][["dataSources"]] <- NULL
      sw_nl_rest <- sw_conf[["nirvanaVersionList"]][[1]] |>
        tibble::as_tibble()
      sw_conf[["nirvanaVersionList"]] <- NULL
      sw_rest <- tibble::as_tibble(sw_conf)
      sw_all <- dplyr::bind_cols(sw_rest, sw_nl_rest)
      sw <- list(
        data_sources = sw_nl_data_sources,
        other = sw_all
      )

      ## biomarkers
      biom <- dat[["biomarkers"]]
      biom_list <- list()
      if ("microsatelliteInstability" %in% names(biom)) {
        msi <- biom[["microsatelliteInstability"]]
        assertthat::assert_that(
          purrr::is_list(msi, n = 3),
          all(c("msiPercentUnstableSites", "additionalMetrics") %in% names(msi)),
          msi[["additionalMetrics"]][[1]][["name"]] == "SumJsd"
        )
        biom_list[["msi_pct_unstable_sites"]] <- msi[["msiPercentUnstableSites"]]
        biom_list[["msi_SumJsd"]] <- msi[["additionalMetrics"]][[1]][["value"]]
      }
      if ("tumorMutationalBurden" %in% names(biom)) {
        tmb <- biom[["tumorMutationalBurden"]]
        amet <- tmb[["additionalMetrics"]] |>
          purrr::map(tibble::as_tibble_row) |>
          dplyr::bind_rows() |>
          dplyr::select("name", "value") |>
          tidyr::pivot_wider(names_from = "name", values_from = "value") |>
          as.list()
        assertthat::assert_that(
          all(c("CodingRegionSizeMb", "SomaticCodingVariantsCount") %in% names(amet))
        )
        biom_list[["tmb_per_mb"]] <- tmb[["tumorMutationalBurdenPerMegabase"]]
        biom_list[["tmb_coding_region_sizemb"]] <- amet[["CodingRegionSizeMb"]]
        biom_list[["tmb_somatic_coding_variants_count"]] <- amet[["SomaticCodingVariantsCount"]]
      }
      biom_tbl <- tibble::as_tibble_row(biom_list)
      empty_tbl <- function(cnames) {
        cnames |>
          purrr::map_dfc(setNames, object = list(logical()))
      }
      if (length(biom_list) == 0) {
        biom_tbl <- c(
          "msi_pct_unstable_sites", "msi_SumJsd", "tmb_per_mb",
          "tmb_coding_region_sizemb", "tmb_somatic_coding_variants_count"
        ) |>
          empty_tbl()
      }

      ## sampleMetrics
      qc2tib <- function(el) {
        el[["metrics"]] |>
          purrr::map(tibble::as_tibble) |>
          dplyr::bind_rows()
      }
      smet <- dat[["sampleMetrics"]]
      smet_em <- smet[["expandedMetrics"]][[1]][["metrics"]] |>
        purrr::map(tibble::as_tibble_row) |>
        dplyr::bind_rows() |>
        dplyr::select("name", "value") |>
        tidyr::pivot_wider(names_from = "name", values_from = "value")
      smet_qc <- smet[["qualityControlMetrics"]]
      smet_nms <- purrr::map_chr(smet_qc, "name")
      smet_qc <- smet_qc |>
        purrr::map(qc2tib) |>
        purrr::set_names(smet_nms) |>
        dplyr::bind_rows(.id = "metric") |>
        # try to tidy up a bit
        dplyr::mutate(
          metric = dplyr::case_when(
            .data$metric == "DNA Library QC Metrics" ~ "dna_qc",
            .data$metric == "DNA Library QC Metrics for Small Variant Calling and TMB" ~ "dna_qc_smallv_tmb",
            .data$metric == "DNA Library QC Metrics for Copy Number Variant Calling" ~ "dna_qc_cnv",
            .data$metric == "DNA Library QC Metrics for MSI and Fusions" ~ "dna_qc_msi_fus",
            TRUE ~ .data$metric
          ),
          name = glue("{.data$metric}_{.data$name}")
        ) |>
        dplyr::select("name", "value", "LSL", "USL") |>
        tidyr::pivot_longer(c("value", "LSL", "USL"), names_to = "variable") |>
        tidyr::pivot_wider(names_from = c("name", "variable"), values_from = "value")

      qc <- dplyr::bind_cols(smet_qc, smet_em)
      snvs <- tso_snv(dat[["variants"]][["smallVariants"]])
      if (nrow(snvs) == 0) {
        snvs <- c(
          "chrom", "pos", "ref", "alt", "af", "qual", "dp_tot", "dp_alt",
          "transcript", "source", "bioType", "aminoAcids", "cdnaPos", "codons",
          "cdsPos", "exons", "geneId", "hgnc", "hgvsc", "hgvsp", "isCanonical",
          "polyPhenScore", "polyPhenPrediction", "proteinId", "proteinPos",
          "siftScore", "siftPrediction", "consequence", "introns"
        ) |>
          empty_tbl()
      }

      cnvs <- tso_cnv(dat[["variants"]][["copyNumberVariants"]])
      if (nrow(cnvs) == 0) {
        cnvs <- c(
          "foldChange", "qual", "copyNumberType", "gene", "chromosome",
          "startPosition", "endPosition"
        ) |>
          empty_tbl()
      }

      res <- list(
        sample_info = sample_info,
        software_config = sw,
        biomarkers = biom_tbl,
        qc = qc,
        snvs = snvs,
        cnvs = cnvs
      )
      res
    },
    #' @description
    #' Writes a tidy version of the `SampleAnalysisResults.json.gz` file output
    #' from TSO.
    #'
    #' @param d Parsed object from `self$read()`.
    #' @param prefix Prefix of output file(s).
    #' @param out_dir Output directory.
    #' @param out_format Format of output file(s) (one of 'tsv' (def.),
    #' 'parquet', 'both').
    #'
    write = function(d, out_dir, prefix, out_format = "tsv") {
      prefix <- file.path(out_dir, prefix)
      p <- glue("{prefix}_SampleAnalysisResults")
      l <- list(
        sample_info = list(
          obj = d[["sample_info"]],
          pref = glue("{p}_sample_info")
        ),
        qc = list(
          obj = d[["qc"]],
          pref = glue("{p}_qc")
        ),
        biomarkers = list(
          obj = d[["biomarkers"]],
          pref = glue("{p}_biomarkers")
        ),
        sw_conf_datasources = list(
          obj = d[["software_config"]][["data_sources"]],
          pref = glue("{p}_sw_conf_datasources")
        ),
        sw_conf_other = list(
          obj = d[["software_config"]][["other"]],
          pref = glue("{p}_sw_conf_other")
        ),
        snv = list(
          obj = d[["snvs"]],
          pref = glue("{p}_smallv")
        ),
        cnv = list(
          obj = d[["cnvs"]],
          pref = glue("{p}_cnv")
        )
      )
      purrr::map(l, function(k) {
        write_dracarys(obj = k[["obj"]], prefix = k[["pref"]], out_format = out_format)
      })
    }
  )
)

tso_snv <- function(snvs) {
  # snvs is an array of snv elements
  snv_info <- function(snv) {
    main <- tibble::tibble(
      chrom = snv[["vcfChromosome"]],
      pos = snv[["vcfPosition"]],
      ref = snv[["vcfRefAllele"]],
      alt = snv[["vcfAltAllele"]],
      af = snv[["vcfVariantFrequency"]],
      qual = snv[["quality"]],
      dp_tot = snv[["totalDepth"]],
      dp_alt = snv[["altAlleleDepth"]]
    )
    # each snv has a single-element array nirvana
    nirv <- snv[["nirvana"]][[1]]
    # nirvana has a transcript array
    txs <- nirv[["transcripts"]]
    get_tx_info <- function(tx) {
      cons <- unlist(tx$consequence) |> paste(collapse = ",")
      tx$consequence <- NULL
      tibble::as_tibble_row(tx) |>
        dplyr::mutate(consequence = cons)
    }
    if (length(txs) > 0) {
      tx_info <- txs |>
        purrr::map_dfr(get_tx_info) |>
        dplyr::bind_rows()
    } else {
      tx_info <- NULL
    }
    dplyr::bind_cols(main, tx_info)
  }
  purrr::map_dfr(snvs, snv_info)
}

tso_cnv <- function(cnvs) {
  purrr::map_dfr(cnvs, tibble::as_tibble)
}
