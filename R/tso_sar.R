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
#' res$write(d_parsed, tempfile(), "tsv")
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
      j <- read_jsongz_jsonlite(x)
      dat <- j[["data"]]
      ## sampleInformation
      sampleinfo <- dat[["sampleInformation"]] |>
        tibble::as_tibble_row()
      ## softwareConfiguration
      sw_conf <- dat[["softwareConfiguration"]]
      sw_nl_data_sources <- sw_conf[["nirvanaVersionList"]][[1]][["dataSources"]] |>
        purrr::map(tibble::as_tibble_row) |>
        dplyr::bind_rows()
      # get rid of it to grab the remaining elements
      sw_conf[["nirvanaVersionList"]][[1]][["dataSources"]] <- NULL
      sw_nl_rest <- sw_conf[["nirvanaVersionList"]][[1]] |>
        tibble::as_tibble_row()
      sw_conf[["nirvanaVersionList"]] <- NULL
      sw_rest <- tibble::as_tibble_row(sw_conf)
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
        biom_list[["msi_sum_jsd"]] <- msi[["additionalMetrics"]][[1]][["value"]]
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
      empty_tbl2 <- function(cnames) {
        cnames |>
          purrr::map_dfc(setNames, object = list(logical()))
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
        dplyr::mutate(name = tolower(.data$name)) |>
        dplyr::select("name", "value")
      smet_qc <- smet[["qualityControlMetrics"]]
      smet_nms <- purrr::map_chr(smet_qc, "name")
      smet_qc <- smet_qc |>
        purrr::map(qc2tib) |>
        dplyr::bind_rows() |>
        dplyr::distinct() |>
        dplyr::mutate(name = tolower(.data$name)) |>
        dplyr::select("name", "value")
      qc <- dplyr::bind_rows(smet_qc, smet_em) |>
        tidyr::pivot_wider(names_from = "name", values_from = "value")
      # also bind qc with biom_tbl if non-empty
      if (length(biom_list) != 0) {
        qc <- dplyr::bind_cols(qc, biom_tbl)
      }

      snvs <- tso_snv(dat[["variants"]][["smallVariants"]])
      if (nrow(snvs) == 0) {
        snvs <- c(
          "chrom", "pos", "ref", "alt", "af", "qual", "dp_tot", "dp_alt",
          "transcript", "source", "bioType", "aminoAcids", "cdnaPos", "codons",
          "cdsPos", "exons", "geneId", "hgnc", "hgvsc", "hgvsp", "isCanonical",
          "polyPhenScore", "polyPhenPrediction", "proteinId", "proteinPos",
          "siftScore", "siftPrediction", "consequence", "introns"
        ) |>
          empty_tbl2()
      }

      cnvs <- tso_cnv(dat[["variants"]][["copyNumberVariants"]])
      if (nrow(cnvs) == 0) {
        cnvs <- c(
          "foldChange", "qual", "copyNumberType", "gene", "chromosome",
          "startPosition", "endPosition"
        ) |>
          empty_tbl2()
      }

      res <- list(
        sampleinfo = sampleinfo,
        qc = qc,
        swconfds = sw[["data_sources"]],
        swconfother = sw[["other"]],
        snv = snvs,
        cnv = cnvs
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
    #' @param out_format Format of output file(s).
    #' @param drid dracarys ID to use for the dataset (e.g. `wfrid.123`, `prid.456`).
    write = function(d, out_dir = NULL, prefix, out_format = "tsv", drid = NULL) {
      if (!is.null(out_dir)) {
        prefix <- file.path(out_dir, prefix)
      }
      d_write <- d |>
        tibble::enframe(name = "section") |>
        dplyr::rowwise() |>
        dplyr::mutate(
          section_low = tolower(.data$section),
          p = glue("{prefix}_{.data$section_low}"),
          out = list(write_dracarys(obj = .data$value, prefix = .data$p, out_format = out_format, drid = drid))
        ) |>
        dplyr::ungroup() |>
        dplyr::select("section", "value") |>
        tibble::deframe()
      invisible(d_write)
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
