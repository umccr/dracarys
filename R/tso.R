#' Wf_tso_ctdna_tumor_only R6 Class
#'
#' @description
#' Reads and writes tidy versions of files from the `tso_ctdna_tumor_only` workflow.
#'
#' @examples
#' \dontrun{
#'
#' #---- Local ----#
#' p <- file.path(
#'   "~/icav1/g/production/analysis_data/SBJ04651/tso_ctdna_tumor_only",
#'   "20240223d1951163/L2400183/Results"
#' )
#' SampleID <- "PRJ230876"
#' LibraryID <- "L2400183"
#' prefix <- glue("{SampleID}__{LibraryID}")
#' t1 <- Wf_tso_ctdna_tumor_only$new(path = p, SampleID = SampleID, LibraryID = LibraryID)
#' t1$list_files(max_files = 20)
#' t1$list_files_filter_relevant(max_files = 300)
#' d <- t1$download_files(max_files = 100, dryrun = F)
#' d_tidy <- t1$tidy_files(d)
#' d_write <- t1$write(
#'   d_tidy,
#'   outdir = file.path(p, "dracarys_tidy"),
#'   prefix = prefix,
#'   format = "tsv"
#' )
#'
#' #---- GDS ----#
#' p <- file.path(
#'   "gds://production/analysis_data/SBJ05563/tso_ctdna_tumor_only",
#'   "20240914d41300cd/L2401388/Results"
#' )
#' SampleID <- "PRJ241446"
#' LibraryID <- "L2401388"
#' prefix <- glue("{SampleID}__{LibraryID}")
#' outdir <- file.path(sub("gds:/", "~/icav1/g", p))
#' token <- Sys.getenv("ICA_ACCESS_TOKEN")
#' t2 <- Wf_tso_ctdna_tumor_only$new(path = p, SampleID = SampleID, LibraryID = LibraryID)
#' t2$list_files(max_files = 100)
#' t2$list_files_filter_relevant(max_files = 100)
#' d <- t2$download_files(
#'   outdir = outdir, ica_token = token,
#'   max_files = 100, dryrun = F
#' )
#' d_tidy <- t2$tidy_files(d)
#' d_write <- t2$write(
#'   d_tidy,
#'   outdir = file.path(outdir, "dracarys_tidy"),
#'   prefix = prefix,
#'   format = "tsv"
#' )
#' }
#' @export
Wf_tso_ctdna_tumor_only <- R6::R6Class(
  "Wf_tso_ctdna_tumor_only",
  inherit = Wf,
  public = list(
    #' @field SampleID The SampleID of the tumor sample (needed for path lookup).
    #' @field LibraryID The LibraryID of the tumor sample (needed for path lookup).
    SampleID = NULL,
    LibraryID = NULL,
    #' @description Create a new Wf_tso_ctdna_tumor_only object.
    #' @param path Path to directory with raw workflow results (from GDS, S3, or
    #' local filesystem).
    #' @param SampleID The SampleID of the tumor sample (needed for path lookup).
    #' @param LibraryID The LibraryID of the sample (needed for path lookup).
    initialize = function(path = NULL, SampleID = NULL, LibraryID = NULL) {
      wname <- "tso_ctdna_tumor_only"
      pref <- glue("{SampleID}_{LibraryID}")
      regexes <- tibble::tribble(
        ~regex, ~fun,
        glue("{pref}/{pref}.SampleAnalysisResults\\.json\\.gz$"), "sar",
        glue("{pref}/{pref}_TMB_Trace\\.tsv$"), "tmbt",
        glue("{pref}/{pref}.AlignCollapseFusionCaller_metrics\\.json\\.gz$"), "acfc",
        glue("{pref}/{pref}_MergedSmallVariants\\.vcf\\.gz$"), "msv",
        glue("{pref}/{pref}_MergedSmallVariants\\.vcf\\.gz\\.tbi$"), "DOWNLOAD_ONLY",
        # glue("{pref}/{pref}_MergedSmallVariants\\.genome\\.vcf\\.gz$"), "DOWNLOAD_ONLY",
        # glue("{pref}/{pref}_MergedSmallVariants\\.genome\\.vcf\\.gz\\.tbi$"), "DOWNLOAD_ONLY",
        glue("{pref}/{pref}_CombinedVariantOutput\\.tsv$"), "cvo",
        glue("{pref}/{pref}_CopyNumberVariants\\.vcf\\.gz$"), "cnv",
        glue("{pref}/{pref}_CopyNumberVariants\\.vcf\\.gz\\.tbi$"), "DOWNLOAD_ONLY",
        glue("{pref}/{pref}.fragment_length_hist\\.json\\.gz$"), "flh",
        glue("{pref}/{pref}.TargetRegionCoverage\\.json\\.gz$"), "trc",
        glue("{pref}/{pref}.tmb\\.json\\.gz$"), "tmb",
        glue("{pref}/{pref}.msi\\.json\\.gz$"), "msi",
        glue("{pref}/{pref}_Fusions\\.csv$"), "fus"
      ) |>
        dplyr::mutate(
          fun = paste0("read_", .data$fun),
          fun = ifelse(.data$fun == "read_DOWNLOAD_ONLY", "DOWNLOAD_ONLY", .data$fun)
        )

      super$initialize(path = path, wname = wname, regexes = regexes)
      self$SampleID <- SampleID
      self$LibraryID <- LibraryID
    },
    #' @description Print details about the Workflow.
    #' @param ... (ignored).
    print = function(...) {
      res <- tibble::tribble(
        ~var, ~value,
        "path", self$path,
        "wname", self$wname,
        "filesystem", self$filesystem,
        "SampleID", self$SampleID,
        "LibraryID", self$LibraryID
      )
      print(res)
      invisible(self)
    },
    #' @description Read `SampleAnalysisResults.json.gz` file.
    #' @param x Path to file.
    read_sar = function(x) {
      TsoSampleAnalysisResultsFile$new(x)$read()
    },
    #' @description Read `TMB_Trace.tsv` file.
    #' @param x Path to file.
    read_tmbt = function(x) {
      dat <- TsoTmbTraceTsvFile$new(x)$read()
      tibble::tibble(name = "tmbtrace", data = list(dat))
    },
    #' @description Read `AlignCollapseFusionCaller_metrics.json.gz` file.
    #' @param x Path to file.
    read_acfc = function(x) {
      TsoAlignCollapseFusionCallerMetricsFile$new(x)$read()
    },
    #' @description Read `MergedSmallVariants.vcf.gz` file.
    #' @param x Path to file.
    read_msv = function(x) {
      dat <- TsoMergedSmallVariantsVcfFile$new(x)$read(only_pass = FALSE, alias = FALSE)
      tibble::tibble(name = "mergedsmallv", data = list(dat))
    },
    #' @description Read `MergedSmallVariants.genome.vcf.gz` file.
    #' @param x Path to file.
    read_msvg = function(x) {
      dat <- TsoMergedSmallVariantsGenomeVcfFile$new(x)$read(only_pass = FALSE, alias = FALSE)
      tibble::tibble(name = "mergedsmallvg", data = list(dat))
    },
    #' @description Read `CombinedVariantOutput.tsv` file.
    #' @param x Path to file.
    read_cvo = function(x) {
      dat <- TsoCombinedVariantOutputFile$new(x)$read()
      tibble::tibble(name = "combinedvaro", data = list(dat))
    },
    #' @description Read `CopyNumberVariants.vcf.gz` file.
    #' @param x Path to file.
    read_cnv = function(x) {
      dat <- TsoCopyNumberVariantsVcfFile$new(x)$read(only_pass = FALSE, alias = FALSE)
      tibble::tibble(name = "cnv", data = list(dat))
    },
    #' @description Read `fragment_length_hist.json.gz` file.
    #' @param x Path to file.
    read_flh = function(x) {
      j <- read_jsongz_jsonlite(x)
      cnames <- c("FragmentLength", "Count")
      # handle SBJ00006...
      if (length(j) == 0) {
        return(empty_tbl(cnames = cnames))
      }
      assertthat::assert_that(
        all(names(j[[1]] %in% cnames))
      )
      dat <- j |>
        purrr::map(tibble::as_tibble) |>
        dplyr::bind_rows()
      tibble::tibble(name = "fraglenhist", data = list(dat))
    },
    #' @description Read `TargetRegionCoverage.json.gz` file.
    #' @param x Path to file.
    read_trc = function(x) {
      dat <- TsoTargetRegionCoverageFile$new(x)$read()
      tibble::tibble(name = "targetcvg", data = list(dat))
    },
    #' @description Read `tmb.json.gz` file.
    #' @param x Path to file.
    read_tmb = function(x) {
      j <- read_jsongz_rjsonio(x) # turns NaN to NULL
      # not interested in Settings element
      j[["Settings"]] <- NULL
      # handle silly NULLs
      j <- lapply(j, function(x) ifelse(is.null(x), NA_real_, x))
      dat <- tibble::as_tibble_row(j) |>
        dplyr::mutate(dplyr::across(dplyr::everything(), as.numeric))
      tibble::tibble(name = "tmb", data = list(dat))
    },
    #' @description Read `msi.json.gz` file.
    #' @param x Path to file.
    read_msi = function(x) {
      dat <- TsoMsiFile$new(x)$read()
      tibble::tibble(name = "msi", data = list(dat))
    },
    #' @description Read `Fusions.csv` file.
    #' @param x Path to file.
    read_fus = function(x) {
      dat <- tso_fusions_read(x)
      tibble::tibble(name = "fusions", data = list(dat))
    }
  ) # end public
)

#' Read TSO CombinedVariantOutput File
#'
#' Reads `CombinedVariantOutput.tsv` output from the TSO500 workflow and extracts
#' only the Small Variants section (due to inconsistencies with other sections).
#'
#' @param x Path to file.
#'
#' @export
tso_combinedvaro_smallv_read <- function(x) {
  nm_map <- c(
    gene = "Gene", chrom = "Chromosome", pos = "Genomic Position",
    ref = "Reference Call", alt = "Alternative Call",
    vaf = "Allele Frequency", dp = "Depth",
    pdot = "P-Dot Notation", cdot = "C-Dot Notation",
    csq = "Consequence(s)", exons = "Affected Exon(s)"
  )
  # read full file
  ln <- readr::read_lines(x, skip_empty_rows = TRUE)
  # now construct a tibble of small variants
  smallv <- grep("\\[Small Variants\\]", ln)
  # handle 0 variants
  if (length(smallv) == 0 || ln[(smallv + 2)] == "NA\t\t") {
    return(empty_tbl(names(nm_map)))
  }
  d <- ln[(smallv + 1):length(ln)] |>
    I() |> # read parsed data as-is
    readr::read_tsv(
      col_names = TRUE, col_types = readr::cols(
        .default = "c",
        "Genomic Position" = "i",
        "Allele Frequency" = "d",
        "Depth" = "d"
      )
    ) |>
    dplyr::rename(dplyr::any_of(nm_map))
  d[]
}

#' Read TSO TMB_Trace File
#'
#' Reads the `TMB_Trace.tsv` file output from the TSO500 workflow.
#'
#' @examples
#' x <- system.file("extdata/tso/sample705_TMB_Trace.tsv", package = "dracarys")
#' tso_tmbt_read(x)
#' @export
tso_tmbt_read <- function(x) {
  ct <- list(
    Chromosome = "c", Position = "d", RefCall = "c", AltCall = "c",
    VAF = "d", Depth = "d", CytoBand = "c", GeneName = "c",
    VariantType = "c", CosmicIDs = "c", MaxCosmicCount = "d",
    AlleleCountsGnomadExome = "d", AlleleCountsGnomadGenome = "d",
    AlleleCounts1000Genomes = "d", MaxDatabaseAlleleCounts = "d",
    GermlineFilterDatabase = "l", GermlineFilterProxi = "l",
    CodingVariant = "l", Nonsynonymous = "l", IncludedInTMBNumerator = "l"
  )
  # cttso v2 has introduced a few extra columns
  ct2 <- list(
    Chromosome = "c", Position = "d", RefCall = "c", AltCall = "c",
    VAF = "d", Depth = "d", CytoBand = "c", GeneName = "c",
    VariantType = "c", CosmicIDs = "c", MaxCosmicCount = "d",
    ClinVarIDs = "c", ClinVarSignificance = "c",
    AlleleCountsGnomadExome = "d", AlleleCountsGnomadGenome = "d",
    AlleleCounts1000Genomes = "d", MaxDatabaseAlleleCounts = "d",
    GermlineFilterDatabase = "c", GermlineFilterProxi = "c",
    Nonsynonymous = "c", withinValidTmbRegion = "c",
    IncludedInTMBNumerator = "c", Status = "c", ProteinChange = "c",
    CDSChange = "c", Exons = "c", Consequence = "c"
  )
  hdr <- readr::read_tsv(x, n_max = 0, show_col_types = FALSE)
  if (all(names(ct2) %in% colnames(hdr))) {
    ct <- ct2
  }
  d <- readr::read_tsv(x, col_types = ct)
  d[]
}

#' Plot Fragment Length Hist
#'
#' Plots the fragment length distributions as given in the
#' `fragment_length_hist` file.
#'
#' @param d Parsed tibble.
#' @param min_count Minimum read count to be plotted (def: 10).
#'
#' @return A ggplot2 plot containing fragment lengths on X axis and read counts
#' on Y axis for each sample.
tso_fraglenhist_plot <- function(d, min_count = 10) {
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

#' Read TSO TargetRegionCoverage File
#'
#' Reads the `TargetRegionCoverage.json.gz` file output from the TSO500 workflow.
#'
#' @return tibble with the following columns:
#'
#' - ConsensusReadDepth
#' - BasePair
#' - Percentage
#'
#' @examples
#' x <- system.file("extdata/tso/sample705.TargetRegionCoverage.json.gz", package = "dracarys")
#' d <- tso_targetregcvg_read(x)
#' tso_targetregcvg_plot(d, min_pct = 0)
#' @export
tso_targetregcvg_read <- function(x) {
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
}

#' Plot TargetRegionCoverage
#'
#' Plots stuff from the `TargetRegionCoverage.json.gz` file output from the
#' TSO500 workflow.
#'
#' @param d Parsed tibble.
#' @param min_pct Minimum percentage to be plotted (def: 2).
#' @return A ggplot2 plot containing read depth on X axis and percentage
#'   covered on Y axis.
tso_targetregcvg_plot <- function(d, min_pct = 2) {
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
#' tso_fusions_read(x)
#' @export
tso_fusions_read <- function(x) {
  # extra column (at the end) in v2
  ct <- list(
    Sample = "c", Name = "c", Chr1 = "c", Pos1 = "d", Chr2 = "c",
    Pos2 = "d", Direction = "c", Alt_Depth = "d", BP1_Depth = "d",
    BP2_Depth = "d", Total_Depth = "d", VAF = "d", Gene1 = "c", Gene2 = "c",
    Contig = "c", Filter = "c", Is_Cosmic_GenePair = "l"
  )
  ct2 <- list("Fusion Directionality Known" = "c")
  hdr <- readr::read_csv(x, n_max = 0, comment = "#", show_col_types = FALSE)
  if (all(names(ct2) %in% colnames(hdr))) {
    ct <- c(ct, ct2)
  }
  res <- readr::read_csv(x, col_types = ct, comment = "#")
  if (nrow(res) == 0) {
    return(empty_tbl(cnames = names(ct), ctypes = ct))
  }
  return(res[])
}

#' Read TSO MSI JSON File
#'
#' Reads `msi.json.gz` file output from the TSO500 workflow.
#'
#' @examples
#' x <- system.file("extdata/tso/sample705.msi.json.gz", package = "dracarys")
#' tso_msi_read(x)
#' @export
tso_msi_read <- function(x) {
  j <- read_jsongz_jsonlite(x)
  # not interested in Settings element
  j[["Settings"]] <- NULL
  j[["ResultMessage"]] <- j[["ResultMessage"]] %||% NA_character_
  if (j[["PercentageUnstableSites"]] == "NaN") {
    j[["PercentageUnstableSites"]] <- NA_real_
  }
  tibble::as_tibble_row(j) |>
    dplyr::mutate(ResultIsValid = as.character(.data$ResultIsValid))
}
