#' Dracarys Tso Tidy
#'
#' Generate tidier representations of TSO500 ctDNA workflow outputs.
#'
#' @param indir Path to TSO500 ctDNA workflow results.
#' @param outdir Path to output dracarys results.
#' @param dryrun Just list the files that will be downloaded?
#' @param out_format Format of output (tsv, parquet, both) (def: tsv).
#' @param token ICA access token (by default uses $ICA_ACCESS_TOKEN env var).
#'
#' @return Generates tidy TSV and/or Parquet representations of the workflow results.
#' @examples
#' \dontrun{
#' indir <- "gds://production/analysis_data/SBJ02858/tso_ctdna_tumor_only/20221104b7ad0b38/L2201560/Results/PRJ222206_L2201560/"
#' indir <- "gds://production/analysis_data/SBJ00476/tso_ctdna_tumor_only/20211217c0a5d1f8/L2100345/Results/PRJ200603_L2100345/"
#' indir <- "gds://production/analysis_data/SBJ00006/tso_ctdna_tumor_only/202210111b001d7f/L2201499/Results/NTC_ctTSO221004_L2201499/"
#' outdir <- here::here(glue("nogit/tso/SBJ00476"))
#' outdir <- here::here(glue("nogit/tso/SBJ00006"))
#' dryrun <- TRUE
#' token <- Sys.getenv("ICA_ACCESS_TOKEN_PROD")
#' dracarys_tso(indir, outdir, dryrun = F, token = token)
#' }
#' @export
dracarys_tso <- function(indir, outdir, out_format = "tsv", dryrun = FALSE,
                         token = Sys.getenv("ICA_ACCESS_TOKEN")) {
  output_format_valid(out_format)
  e <- emojifont::emoji
  # cli::cli_div(theme = list(
  #   span.file = list(color = "blue"),
  #   span.emph = list(color = "orange")
  # ))
  cli::cli_alert_info("{date_log()} {e('dragon')} Start tidying {.file {indir}} {e('fire')}")
  gds_local_dir <- if (grepl("^gds://", indir)) file.path(outdir, "dracarys_gds_sync") else NULL
  # main dracarys function
  d <- tso_tidy(indir = indir, gds_local_dir = gds_local_dir, dryrun = dryrun, token = token)
  if (!dryrun) {
    fs::dir_create(outdir)
    res <- file.path(outdir, "res.rds")
    readr::write_rds(d, res)

    # if (out_format %in% c("tsv", "both")) {
    #   tsv_out <- file.path(outdir, glue("{prefix}.tsv"))
    #   readr::write_tsv(d1, tsv_out)
    # }
    # if (out_format %in% c("parquet", "both")) {
    #   parquet_out <- file.path(outdir, glue("{prefix}.parquet"))
    #   arrow::write_parquet(d1, parquet_out)
    # }

    cli::cli_alert_success("{date_log()} {e('rocket')} End tidying {.file {indir}} {e('comet')}!")
    cli::cli_alert_info("{date_log()} {e('tada')} Path to output directory with results for {.emph {indir}}: {.file {outdir}}")
  }
}

tso_write <- function(d, outdir, prefix, format = "tsv") {
  tso_fnames <- FILE_REGEX |>
    dplyr::filter(grepl("tso__", .data$name)) |>
    dplyr::pull("name")
  assertthat::assert_that(purrr::is_list(d, n = length(tso_fnames)))
  output_format_valid(format)
}

#' Tidy TSO ctDNA Results
#'
#' Tidies TSO ctDNA results into a list of tibbles.
#'
#' @param indir Directory path to TSO500 ctDNA workflow results (must end
#' with `/`).
#' @param gds_local_dir If `indir` is a GDS directory, 'recognisable' files
#' will be first downloaded to this directory.
#' @param dryrun Just list the files that will be downloaded?
#' @param token ICA access token (by default uses $ICA_ACCESS_TOKEN env var).
#'
#' @return List of tibbles.
#' @examples
#' \dontrun{
#' indir <- "gds://production/analysis_data/SBJ02858/tso_ctdna_tumor_only/20221104b7ad0b38/L2201560/Results/PRJ222206_L2201560/"
#' gds_local_dir <- here::here(glue("nogit/tso/SBJ02858"))
#' dryrun <- F
#' tso_tidy(indir, gds_local_dir, dryrun)
#' }
#' @export
tso_tidy <- function(indir, gds_local_dir = NULL, dryrun = FALSE,
                     token = Sys.getenv("ICA_ACCESS_TOKEN")) {
  # - List indir files
  # - If GDS, download locally to gdslocaldir
  #   - Grab only the 'recognisable' ones
  # - Apply the tidy functions to each
  # - Export as list of tibbles
  pat <- "tso__"
  e <- emojifont::emoji
  if (grepl("^gds://", indir)) {
    assertthat::assert_that(
      !is.null(gds_local_dir),
      msg = "You need to specify a local directory to download the GDS files."
    )
    print(dr_gds_download(
      gdsdir = indir, outdir = gds_local_dir, token = token,
      pattern = pat, dryrun = dryrun
    ))
    # Use the downloaded results
    indir <- gds_local_dir
  } else {
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
    d <- fs::dir_ls(indir) |>
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
        "No TSO files for dracarys were found in {.file {indir}}.",
        "See current supported regexes for TSO: {regex}."
      )
      cli::cli_abort(msg)
    }

    d <- d |>
      dplyr::select("type", "path", "bname") |>
      dplyr::rowwise() |>
      dplyr::mutate(
        res = list(tso_funcall(.data$type)$new(.data$path)$read()) |>
          purrr::set_names(.data$type)
      )

    d[["res"]]
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
#' d$read() # or read(d)
#' d$write(prefix = tempfile(), out_format = "both")
#' @export
TsoTmbTraceTsvFile <- R6::R6Class("TsoTmbTraceTsvFile",
  inherit = File, public = list(
    #' @description
    #' Reads the `TMB_Trace.tsv` file output from TSO.
    #'
    #' @return tibble with the following columns:
    #' * FragmentLength
    #' * Count
    read = function() {
      x <- self$path
      ct <- readr::cols(
        Chromosome = "c",
        Position = "i",
        RefCall = "c",
        AltCall = "c",
        VAF = "d",
        Depth = "d",
        CytoBand = "c",
        GeneName = "c",
        VariantType = "c",
        CosmicIDs = "c",
        MaxCosmicCount = "d",
        AlleleCountsGnomadExome = "d",
        AlleleCountsGnomadGenome = "d",
        AlleleCounts1000Genomes = "d",
        MaxDatabaseAlleleCounts = "d",
        GermlineFilterDatabase = "l",
        GermlineFilterProxi = "l",
        CodingVariant = "l",
        Nonsynonymous = "l",
        IncludedInTMBNumerator = "l"
      )
      readr::read_tsv(x, col_types = ct)
    },

    #' @description
    #' Writes a tidy version of the `TMB_Trace.tsv` file output from TSO.
    #'
    #' @param prefix Prefix path of output file(s).
    #' @param out_format Format of output file(s) (one of 'tsv' (def.),
    #' 'parquet', 'both').
    write = function(prefix, out_format = "tsv") {
      d <- self$read()
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
#' fl$read() # or read(fl)
#' fl$plot(5)
#' fl$write(tempfile(), out_format = "both")
#' @export
TsoFragmentLengthHistFile <- R6::R6Class("TsoFragmentLengthHistFile",
  inherit = File, public = list(
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
    #' @param prefix Prefix path of output file(s).
    #' @param out_format Format of output file(s) (one of 'tsv' (def.),
    #' 'parquet', 'both').
    #'
    write = function(prefix, out_format = "tsv") {
      d <- self$read()
      prefix2 <- glue("{prefix}_fragment_length_hist")
      write_dracarys(obj = d, prefix = prefix2, out_format = out_format)
    },


    #' @description Plots the fragment length distributions as given in the
    #' `fragment_length_hist.json.gz` file.
    #'
    #' @param min_count Minimum read count to be plotted (Default: 10).
    #' @return A ggplot2 plot containing fragment lengths on X axis and read counts
    #'   on Y axis for each sample.
    plot = function(min_count = 10) {
      assertthat::assert_that(is.numeric(min_count), min_count >= 0)
      d <- self$read() |>
        dplyr::filter(.data$Count >= min_count)

      d |>
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
#' trc$read() # or read(trc)
#' trc$plot(90) # or plot(trc, 90)
#' trc$write(prefix = tempfile(), out_format = "both")
#' @export
TsoTargetRegionCoverageFile <- R6::R6Class("TsoTargetRegionCoverageFile",
  inherit = File, public = list(
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
    #' @param prefix Prefix path of output file(s).
    #' @param out_format Format of output file(s) (one of 'tsv' (def.),
    #' 'parquet', 'both').
    #'
    write = function(prefix, out_format = "tsv") {
      d <- self$read()
      prefix2 <- glue("{prefix}_TargetRegionCoverage")
      write_dracarys(obj = d, prefix = prefix2, out_format = out_format)
    },

    #' @description Plots the `TargetRegionCoverage.json.gz` file.
    #'
    #' @param min_pct Minimum percentage to be plotted (Default: 2).
    #' @return A ggplot2 plot containing read depth on X axis and percentage
    #'   covered on Y axis.
    plot = function(min_pct = 2) {
      assertthat::assert_that(is.numeric(min_pct), min_pct >= 0)
      d <- self$read() |>
        dplyr::filter(
          !.data$ConsensusReadDepth == "TargetRegion",
          .data$Percentage >= min_pct
        ) |>
        dplyr::select(dp = "ConsensusReadDepth", pct = "Percentage") |>
        dplyr::mutate(dp = as.numeric(sub("X", "", .data$dp)))
      d |>
        ggplot2::ggplot(aes(x = dp, y = pct, label = dp)) +
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
#' m$read() # or read(m)
#' m$write(prefix = tempfile(), out_format = "both")
#' @export
TsoAlignCollapseFusionCallerMetricsFile <- R6::R6Class("TsoAlignCollapseFusionCallerMetricsFile",
  inherit = File, public = list(
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
    #' @param prefix Prefix path of output file(s).
    #' @param out_format Format of output file(s) (one of 'tsv' (def.),
    #' 'parquet', 'both').
    #'
    write = function(prefix, out_format = "tsv") {
      d <- self$read()
      dhist <- d |>
        dplyr::filter(grepl("Hist", .data$name))
      assertthat::assert_that(all(dhist[["section"]] == "UmiStatistics"))
      assertthat::assert_that(all(dhist[["percent"]] %in% NA))
      dhist <- dhist |>
        dplyr::mutate(
          name = sub("Histogram of ", "", .data$name),
          name = gsub(" ", "_", .data$name),
          value = as.numeric(.data$value)
        ) |>
        dplyr::group_by(name) |>
        dplyr::mutate(num = dplyr::row_number()) |>
        dplyr::ungroup() |>
        dplyr::select(c("name", "num", "value"))
      dmain <- d |>
        dplyr::filter(!grepl("Hist", .data$name))

      p <- glue("{prefix}_AlignCollapseFusionCaller_metrics_")
      p_hist <- glue("{p}_hist")
      p_main <- glue("{p}_main")
      write_dracarys(obj = dhist, prefix = p_hist, out_format = out_format)
      write_dracarys(obj = dmain, prefix = p_main, out_format = out_format)
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
#' tmb$read() # or read(tmb)
#' tmb$write(tempfile(), "both")
#' @export
TsoTmbFile <- R6::R6Class("TsoTmbFile", inherit = File, public = list(
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
  #' @param prefix Prefix path of output file(s).
  #' @param out_format Format of output file(s) (one of 'tsv' (def.),
  #' 'parquet', 'both').
  #'
  write = function(prefix, out_format = "tsv") {
    d <- self$read()
    prefix2 <- glue("{prefix}_tmb")
    write_dracarys(obj = d, prefix = prefix2, out_format = out_format)
  }
))

#' TsoFusionsCsvFile R6 Class
#'
#' @description
#' Contains methods for reading and displaying contents of the
#' `Fusions.csv` file output from TSO.
#'
#' @examples
#' x <- system.file("extdata/tso/sample705_Fusions.csv", package = "dracarys")
#' fus <- TsoFusionsCsvFile$new(x)
#' fus$read() # or read(fus)
#' fus$write(tempfile(), "both")
#' @export
TsoFusionsCsvFile <- R6::R6Class("TsoFusionsCsvFile", inherit = File, public = list(
  #' @description
  #' Reads the `Fusions.csv` file output from TSO.
  #'
  #' @return tibble with the following columns:
  #'   - label:
  read = function() {
    x <- self$path
    readr::read_csv(x, col_types = readr::cols(.default = "c"), comment = "#")
  },

  #' @description
  #' Writes a tidy version of the `Fusions.csv` file output from TSO.
  #'
  #' @param prefix Prefix path of output file(s).
  #' @param out_format Format of output file(s) (one of 'tsv' (def.),
  #' 'parquet', 'both').
  #'
  write = function(prefix, out_format = "tsv") {
    d <- self$read()
    prefix2 <- glue("{prefix}_Fusions")
    write_dracarys(obj = d, prefix = prefix2, out_format = out_format)
  }
))

#' TsoMsiFile R6 Class
#'
#' @description
#' Contains methods for reading and displaying contents of the
#' `msi.json.gz` file output from TSO.
#'
#' @examples
#' x <- system.file("extdata/tso/sample705.msi.json.gz", package = "dracarys")
#' msi <- TsoMsiFile$new(x)
#' msi$read() # or read(msi)
#' msi$write(tempfile(), "both")
#' @export
TsoMsiFile <- R6::R6Class("TsoMsiFile", inherit = File, public = list(
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
  #' @param prefix Prefix path of output file(s).
  #' @param out_format Format of output file(s) (one of 'tsv' (def.),
  #' 'parquet', 'both').
  #'
  write = function(prefix, out_format = "tsv") {
    d <- self$read()
    prefix2 <- glue("{prefix}_msi")
    write_dracarys(obj = d, prefix = prefix2, out_format = out_format)
  }
))

#' TsoSampleAnalysisResultsFile R6 Class
#'
#' @description
#' Contains methods for reading and displaying contents of the
#' `SampleAnalysisResults.json.gz` file output from TSO.
#'
#' @examples
#' x <- system.file("extdata/tso/sample705_SampleAnalysisResults.json.gz", package = "dracarys")
#' res <- TsoSampleAnalysisResultsFile$new(x)
#' res$read() # or read(res)
#' @export
TsoSampleAnalysisResultsFile <- R6::R6Class("TsoSampleAnalysisResultsFile", inherit = File, public = list(
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

    ## sampleMetrics
    smet <- dat[["sampleMetrics"]]
    smet_em <- smet[["expandedMetrics"]][[1]][["metrics"]] |>
      purrr::map(tibble::as_tibble_row) |>
      dplyr::bind_rows() |>
      dplyr::select("name", "value") |>
      tidyr::pivot_wider(names_from = "name", values_from = "value")

    qc2tib <- function(el) {
      el[["metrics"]] |>
        purrr::map(tibble::as_tibble) |>
        dplyr::bind_rows()
    }

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

    v <- dat[["variants"]]
    vars <- list(
      # fusions are more comprehensive in the Fusions.csv file
      snvs = tso_snv(v[["smallVariants"]]),
      cnvs = tso_cnv(v[["copyNumberVariants"]])
    )

    res <- list(
      sample_info = sample_info,
      software_config = sw,
      biomarkers = biom_tbl,
      sample_metrics_qc = smet_qc,
      sample_metrics_expanded = smet_em,
      vars = vars
    )
    res
  }
))

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
