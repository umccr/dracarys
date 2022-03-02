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
  #'   - TmbPerMb
  #'   - AdjustedTmbPerMb
  #'   - NonsynonymousTmbPerMb
  #'   - AdjustedNonsynonymousTmbPerMb
  #'   - SomaticCodingVariantsCount
  #'   - NonsynonymousSomaticCodingVariantsCount
  #'   - TotalRegionSizeMb
  #'   - CodingRegionSizeMb
  read = function() {
    x <- self$path
    j <- jsonlite::read_json(x)
    # not interested in Settings element
    j[["Settings"]] <- NULL
    tibble::as_tibble_row(j) |>
      dplyr::mutate(dplyr::across(.cols = dplyr::everything(), as.numeric))
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

# d <-
#   tibble(
#     x = list.files(here::here("nogit/tso/tmb"),
#                    pattern = "tmb\\.json\\.gz$",
#                    recursive = TRUE, full.names = TRUE)) |>
#   rowwise() |>
#   mutate(sample = sub(".tmb.json.gz", "", basename(x)),
#          y = list(tso_read_tmb(x))) |>
#   unnest(y) |>
#   select(-x) |>
#   pivot_longer(TmbPerMb:CodingRegionSizeMb)
#
# theme_set(theme_bw())
# d |>
#   ggplot(aes(x = "", y = value, label = sample, colour = sample)) +
#   geom_violin(fill = "transparent", colour = "grey80", alpha = 0.4) +
#   # geom_point() +
#   ggforce::geom_sina(aes(
#     group = name,
#     colour = sample,
#   ), seed = 42) +
#   facet_wrap(~name, scales = "free")
#
# d |>
#   ggplot(aes(x = value)) +
#   geom_histogram(fill = "purple", colour = "black") +
#   facet_wrap(~name, scales = "free")

# arrow::read_parquet(here::here("nogit/tso/tmb/56.TMB.parquet"))

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
TsoSampleAnalysisResultsFile <- R6::R6Class("TsoMsiFile", inherit = File, public = list(
  #' @description
  #' Reads the `msi.json.gz` file output from TSO.
  #'
  #' @return tibble with the following columns:
  #'         - TOTAL_PF_READS
  #'         - MEAN_FAMILY_SIZE
  #'         - MEDIAN_TARGET_COVERAGE,
  #'         - PCT_CHIMERIC_READS
  #'         - PCT_EXON_500X
  #'         - PCT_EXON_1500X
  #'         - PCT_READ_ENRICHMENT,
  #'         - PCT_USABLE_UMI_READS
  #'         - MEAN_TARGET_COVERAGE
  #'         - PCT_ALIGNED_READS,
  #'         - PCT_CONTAMINATION_EST
  #'         - PCT_TARGET_0.4X_MEAN
  #'         - PCT_TARGET_500X,
  #'         - PCT_TARGET_1000X
  #'         - PCT_TARGET_1500X
  #'         - PCT_DUPLEXFAMILIES,
  #'         - MEDIAN_INSERT_SIZE
  #'         - MAX_SOMATIC_AF
  read = function() {
    x <- self$path
    # just read sampleMetrics for now
    d <- jsonlite::read_json(x)[["data"]][["sampleMetrics"]]
    dplyr::bind_rows(d[["expandedMetrics"]][[1]][["metrics"]]) |>
      dplyr::select(.data$name, .data$value) |>
      tidyr::pivot_wider(.data$name)
  }
))
