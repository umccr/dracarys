#' TsoCombinedVariantOutputFile R6 Class
#'
#' @description
#' Contains methods for reading and displaying contents of the
#' `CombinedVariantOutput.tsv` file output from TSO.
#'
#' @examples
#' x <- system.file("extdata/tso/sample705_CombinedVariantOutput.tsv", package = "dracarys")
#' d <- TsoCombinedVariantOutputFile$new(x)
#' d$read() # or read(d)
#' @export
TsoCombinedVariantOutputFile <- R6::R6Class("TsoCombinedVariantOutputFile",
  inherit = File, public = list(
    #' @description
    #' Reads the `CombinedVariantOutput.tsv` file output from TSO.
    #'
    #' @return list of following tibbles:
    #' * FragmentLength
    #' * Count
    read = function() {
      x <- self$path
      .read_section <- function(section, start, end, lines) {
        snms <- c(
          "Analysis Details", "Sequencing Run Details", "TMB", "MSI",
          "Copy Number Variants", "DNA Fusions", "Small Variants"
        )
        assertthat::assert_that(section %in% snms)
        chunk <- lines[start:end]
        s <- paste(chunk, collapse = "\n")
        d <- tibble::tibble() # empty tibble by default
        if (section %in% c("Analysis Details", "TMB", "MSI")) {
          if (section == "MSI") {
            # single row (as far as I've seen), so the collapse above doesn't touch it
            s <- paste0(s, "\n")
          }
          d <- s |>
            readr::read_tsv(
              col_types = readr::cols(.default = "c"),
              col_names = c("variable", "value", "empty")
            ) |>
            dplyr::select(.data$variable, .data$value)
        } else if (section %in% c("Copy Number Variants", "DNA Fusions", "Small Variants")) {
          if (length(chunk) > 2) {
            d <- s |>
              readr::read_tsv(
                col_types = readr::cols(.default = "c"),
                col_names = TRUE
              )
          }
        }
        return(d)
      } # read_section end
      regx <- list(blank = "^\\s*$", na = "NA\t\t", header = "^\\[(.*)\\]\t\t$")
      lines <- readr::read_lines(x)
      empty_lines <- grepl(regx$blank, lines)
      lines <- lines[!empty_lines]
      headers <- base::which(grepl(regx$header, lines))
      headers_names <- sub(regx$header, "\\1", lines[headers])
      tsv_chunks <- tibble::tibble(
        start = c(headers, length(lines) + 1),
        nm = c(headers_names, "BLANK")
      ) |>
        dplyr::mutate(end = dplyr::lead(.data$start) - 1) |>
        dplyr::filter(!is.na(.data$end)) |>
        dplyr::mutate(start = .data$start + 1) |>
        dplyr::rowwise() |>
        dplyr::mutate(d = list(.read_section(.data$nm, .data$start, .data$end, lines))) |>
        dplyr::select(nm, d)

      tsv_chunks[["d"]] |>
        purrr::set_names(tsv_chunks[["nm"]])
    }
  )
)

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
#' fl$plot(5)
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
      amet <- msi[["additionalMetrics"]] |>
        purrr::map(tibble::as_tibble_row) |>
        dplyr::bind_rows()
      biom_list[["msi"]] <- list(
        msi_pct_unstable_sites = msi[["msiPercentUnstableSites"]],
        additional_metrics = amet
      )
    }
    if ("tumorMutationalBurden" %in% names(biom)) {
      tmb <- biom[["tumorMutationalBurden"]]
      amet <- tmb[["additionalMetrics"]] |>
        purrr::map(tibble::as_tibble_row) |>
        dplyr::bind_rows()
      biom_list[["tmb"]] <- list(
        tmb_per_mb = tmb[["tumorMutationalBurdenPerMegabase"]],
        additional_metrics = amet
      )
    }

    ## sampleMetrics
    smet <- dat[["sampleMetrics"]]
    smet_em <- smet[["expandedMetrics"]][[1]][["metrics"]] |>
      purrr::map(tibble::as_tibble_row) |>
      dplyr::bind_rows() |>
      dplyr::select(.data$name, .data$value) |>
      tidyr::pivot_wider(names_from = .data$name, values_from = .data$value)

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
      dplyr::bind_rows(.id = "metric")

    v <- dat[["variants"]]
    vars <- list(
      # fusions are more comprehensive in the Fusions.csv file
      snvs = tso_snv(v[["smallVariants"]]),
      cnvs = tso_cnv(v[["copyNumberVariants"]])
    )

    res <- list(
      sample_info = sample_info,
      software_config = sw,
      biomarkers = biom_list,
      sample_metrics_qc = smet_qc,
      sample_metrics_expanded = smet_em,
      vars = vars
    )
    res
  }
))
