#' TsoCopyNumberVariantsVcfFile R6 Class
#'
#' @description
#' Contains methods for reading and displaying contents of the
#' `CopyNumberVariants.vcf.gz` file output from TSO.
#'
#' @examples
#' \dontrun{
#' x <- normalizePath("CopyNumberVariants.vcf.gz")
#' d <- TsoCopyNumberVariantsVcfFile$new(x)
#' d_parsed <- d$read() # or read(d)
#' d$write(d_parsed, out_dir = tempdir(), prefix = "sample705", out_format = "tsv")
#' d$write(d_parsed, prefix = "FOO", out_format = "delta", drid = "wfr.123")
#' }
#' @export
TsoCopyNumberVariantsVcfFile <- R6::R6Class(
  "TsoCopyNumberVariantsVcfFile",
  inherit = File,
  public = list(
    #' @description
    #' Reads the `CopyNumberVariants.vcf.gz` file output from TSO.
    #' @param only_pass Only include PASS variants (def: TRUE).
    #' @param alias Substitute sample names with S1/S2/... alias (def: TRUE).
    #'
    #' @return tibble with variants.
    read = function(only_pass = TRUE, alias = TRUE) {
      x <- self$path
      bcftools_parse_vcf(x, only_pass = only_pass, alias = alias)
    },
    #' @description
    #' Writes a tidy version of the `CopyNumberVariants.vcf.gz` file output from TSO.
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
      # prefix2 <- glue("{prefix}copynumber_variants")
      write_dracarys(obj = d, prefix = prefix, out_format = out_format, drid = drid)
    }
  )
)

#' TsoMergedSmallVariantsVcfFile R6 Class
#'
#' @description
#' Contains methods for reading and displaying contents of the
#' `MergedSmallVariants.vcf.gz` file output from TSO.
#'
#' @examples
#' \dontrun{
#' x <- "MergedSmallVariants.vcf.gz"
#' d <- TsoMergedSmallVariantsVcfFile$new(x)
#' d_parsed <- d$read() # or read(d)
#' d$write(d_parsed, out_dir = tempdir(), prefix = "sample705", out_format = "tsv")
#' d$write(d_parsed, prefix = "FOO", out_format = "delta", drid = "wfr.123")
#' }
#' @export
TsoMergedSmallVariantsVcfFile <- R6::R6Class(
  "TsoMergedSmallVariantsVcfFile",
  inherit = File,
  public = list(
    #' @description
    #' Reads the `MergedSmallVariants.vcf.gz` file output from TSO.
    #' @param only_pass Only include PASS variants (def: TRUE).
    #' @param alias Substitute sample names with S1/S2/... alias (def: TRUE).
    #'
    #' @return tibble with variants.
    read = function(only_pass = TRUE, alias = TRUE) {
      x <- self$path
      bcftools_parse_vcf(x, only_pass = only_pass, alias = alias)
    },
    #' @description
    #' Writes a tidy version of the `MergedSmallVariants.vcf.gz` file output from TSO.
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
      # prefix2 <- glue("{prefix}merged_small_variants")
      write_dracarys(obj = d, prefix = prefix, out_format = out_format, drid = drid)
    }
  )
)

#' TsoMergedSmallVariantsGenomeVcfFile R6 Class
#'
#' @description
#' Contains methods for reading and displaying contents of the
#' `MergedSmallVariants.genome.vcf.gz` file output from TSO.
#'
#' @examples
#' \dontrun{
#' x <- "MergedSmallVariants.genome.vcf.gz"
#' d <- TsoMergedSmallVariantsGenomeVcfFile$new(x)
#' d_parsed <- d$read() # or read(d)
#' d$write(d_parsed, out_dir = tempdir(), prefix = "sample705", out_format = "tsv")
#' }
#' @export
TsoMergedSmallVariantsGenomeVcfFile <- R6::R6Class(
  "TsoMergedSmallVariantsGenomeVcfFile",
  inherit = File,
  public = list(
    #' @description
    #' Reads the `MergedSmallVariants.genome.vcf.gz` file output from TSO.
    #' @param only_pass Only include PASS variants (def: TRUE).
    #' @param alias Substitute sample names with S1/S2/... alias (def: TRUE).
    #'
    #' @return tibble with variants.
    read = function(only_pass = TRUE, alias = TRUE) {
      x <- self$path
      bcftools_parse_vcf(x, only_pass = only_pass, alias = alias)
    },
    #' @description
    #' Writes a tidy version of the `MergedSmallVariants.genome.vcf.gz` file output from TSO.
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
      # prefix2 <- glue("{prefix}merged_small_variants")
      write_dracarys(obj = d, prefix = prefix, out_format = out_format, drid = drid)
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
#' d_parsed <- d$read() # or read(d)
#' d$write(d_parsed, out_dir = tempdir(), prefix = "sample705", out_format = "tsv")
#' @export
TsoTmbTraceTsvFile <- R6::R6Class(
  "TsoTmbTraceTsvFile",
  inherit = File,
  public = list(
    #' @description
    #' Reads the `TMB_Trace.tsv` file output from TSO.
    #'
    #' @return tibble with variants.
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
    #' @param out_format Format of output file(s).
    #' @param drid dracarys ID to use for the dataset (e.g. `wfrid.123`, `prid.456`).
    write = function(d, out_dir = NULL, prefix, out_format = "tsv", drid = NULL) {
      if (!is.null(out_dir)) {
        prefix <- file.path(out_dir, prefix)
      }
      # prefix2 <- glue("{prefix}tmb_trace")
      write_dracarys(obj = d, prefix = prefix, out_format = out_format, drid = drid)
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
#' fl$write(d_parsed, out_dir = tempdir(), prefix = "sample705", out_format = "tsv")
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
      j <- read_jsongz_jsonlite(x)
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
    #' @param out_format Format of output file(s).
    #' @param drid dracarys ID to use for the dataset (e.g. `wfrid.123`, `prid.456`).
    write = function(d, out_dir = NULL, prefix, out_format = "tsv", drid = NULL) {
      if (!is.null(out_dir)) {
        prefix <- file.path(out_dir, prefix)
      }
      # prefix2 <- glue("{prefix}fragment_length_hist")
      write_dracarys(obj = d, prefix = prefix, out_format = out_format, drid = drid)
    },

    #' @description Plots the fragment length distributions as given in the
    #' `fragment_length_hist.json.gz` file.
    #'
    #' @param d Parsed object from `self$read()`.
    #' @param min_count Minimum read count to be plotted (def: 10).
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
#' trc$plot(d_parsed, 0) # or plot(trc, d_parsed, 90)
#' trc$write(d_parsed, out_dir = tempdir(), prefix = "sample705", out_format = "tsv")
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
      j <- read_jsongz_jsonlite(x)
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
    #' @param out_format Format of output file(s).
    #' @param drid dracarys ID to use for the dataset (e.g. `wfrid.123`, `prid.456`).
    write = function(d, out_dir = NULL, prefix, out_format = "tsv", drid = NULL) {
      if (!is.null(out_dir)) {
        prefix <- file.path(out_dir, prefix)
      }
      # prefix2 <- glue("{prefix}target_region_coverage")
      write_dracarys(obj = d, prefix = prefix, out_format = out_format, drid = drid)
    },

    #' @description Plots the `TargetRegionCoverage.json.gz` file.
    #'
    #' @param d Parsed object from `self$read()`.
    #' @param min_pct Minimum percentage to be plotted (def: 2).
    #' @importFrom ggplot2 ggplot
    #' @importFrom ggrepel geom_text_repel
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
#' tmb$write(d_parsed, out_dir = tempdir(), prefix = "sample705", out_format = "tsv")
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
      j <- read_jsongz_rjsonio(x) # turns NaN to NULL
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
    #' @param out_format Format of output file(s).
    #' @param drid dracarys ID to use for the dataset (e.g. `wfrid.123`, `prid.456`).
    write = function(d, out_dir = NULL, prefix, out_format = "tsv", drid = NULL) {
      if (!is.null(out_dir)) {
        prefix <- file.path(out_dir, prefix)
      }
      # prefix2 <- glue("{prefix}tmb")
      write_dracarys(obj = d, prefix = prefix, out_format = out_format, drid = drid)
    }
  )
)

#' TsoFusionsCsvFile R6 Class
#'
#' @description
#' Contains methods for reading and displaying contents of the
#' `Fusions.csv` file output from TSO.
#' Returns a tibble where the columns are (based on the input file header metadata):
#'
#' - Sample: input sample ID.
#' - Name: Fusion name as reported by manta.
#' - Chr1/Chr2: The chromosome of the 1st/2nd breakend.
#' - Pos1/Pos2: The position of the 1st/2nd breakend.
#' - Direction: The direction of how the breakends are joined together.
#' - Alt_Depth: The number of read-pairs supporting the fusion call.
#' - BP1_Depth/BP2_Depth: Number of read-pairs aligned to the 1st/2nd breakend.
#' - Total_Depth: Max number of read-pairs aligned to a fusion breakend.
#' - VAF: Variant allele frequency.
#' - Gene1/Gene2: Genes that overlap the 1st/2nd breakend.
#' - Contig: The fusion contig.
#' - Filter: Indicates whether the fusion has passed all of the fusion filters.
#' - Is_Cosmic_GenePair: Indicates whether the gene pair has been reported by Cosmic(True/False).
#'
#' @examples
#' x <- system.file("extdata/tso/sample705_Fusions.csv", package = "dracarys")
#' fus <- TsoFusionsCsvFile$new(x)
#' d_parsed <- fus$read() # or read(fus)
#' fus$write(d_parsed, out_dir = tempdir(), prefix = "sample705", out_format = "tsv")
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
      res <- readr::read_csv(x, col_types = ct, comment = "#")
      if (nrow(res) == 0) {
        return(empty_tbl(cnames = names(ct)))
      }
      return(res)
    },

    #' @description
    #' Writes a tidy version of the `Fusions.csv` file output from TSO.
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
      # prefix2 <- glue("{prefix}fusions")
      write_dracarys(obj = d, prefix = prefix, out_format = out_format, drid = drid)
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
#' msi$write(d_parsed, out_dir = tempdir(), prefix = "sample705", out_format = "tsv")
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
      j <- read_jsongz_jsonlite(x)
      # not interested in Settings element
      j[["Settings"]] <- NULL
      j[["ResultMessage"]] <- j[["ResultMessage"]] %||% NA_character_
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
    #' @param out_format Format of output file(s).
    #' @param drid dracarys ID to use for the dataset (e.g. `wfrid.123`, `prid.456`).
    write = function(d, out_dir = NULL, prefix, out_format = "tsv", drid = NULL) {
      if (!is.null(out_dir)) {
        prefix <- file.path(out_dir, prefix)
      }
      # prefix2 <- glue("{prefix}msi")
      write_dracarys(obj = d, prefix = prefix, out_format = out_format, drid = drid)
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
