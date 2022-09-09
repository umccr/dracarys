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
        dplyr::bind_rows() |>
        dplyr::mutate(FragmentLength = as.character(.data$FragmentLength))
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
        dplyr::bind_rows()
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
          tibble::as_tibble(l) |>
            dplyr::mutate(dplyr::across(dplyr::everything(), ~ as.character(.)))
        }
        el |>
          purrr::map(fun1) |>
          dplyr::bind_rows()
      }
      j <- jsonlite::read_json(x)
      j |>
        purrr::map(l2tib) |>
        dplyr::bind_rows(.id = "section")
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
    j <- jsonlite::read_json(x)
    # not interested in Settings element
    j[["Settings"]] <- NULL
    tibble::as_tibble_row(j) |>
      dplyr::mutate(dplyr::across(dplyr::everything(), as.numeric))
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
    tibble::as_tibble_row(j) |>
      dplyr::mutate(ResultIsValid = as.character(.data$ResultIsValid))
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
  #' @return tibble with the following columns:
  #' * `TOTAL_PF_READS`
  #' * `MEAN_FAMILY_SIZE`
  #' * `MEDIAN_TARGET_COVERAGE`
  #' * `PCT_CHIMERIC_READS`
  #' * `PCT_EXON_500X`
  #' * `PCT_EXON_1500X`
  #' * `PCT_READ_ENRICHMENT`
  #' * `PCT_USABLE_UMI_READS`
  #' * `MEAN_TARGET_COVERAGE`
  #' * `PCT_ALIGNED_READS`
  #' * `PCT_CONTAMINATION_EST`
  #' * `PCT_TARGET_0.4X_MEAN`
  #' * `PCT_TARGET_500X`
  #' * `PCT_TARGET_1000X`
  #' * `PCT_TARGET_1500X`
  #' * `PCT_DUPLEXFAMILIES`
  #' * `MEDIAN_INSERT_SIZE`
  #' * `MAX_SOMATIC_AF`
  read = function() {
    x <- self$path
    # just read sampleMetrics for now
    d <- jsonlite::read_json(x)[["data"]][["sampleMetrics"]]
    dplyr::bind_rows(d[["expandedMetrics"]][[1]][["metrics"]]) |>
      dplyr::select(.data$name, .data$value) |>
      tidyr::pivot_wider(names_from = .data$name, values_from = .data$value)
  }
))
