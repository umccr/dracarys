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
#'
#' #---- GDS ----#
#' p <- file.path(
#'   "gds://production/analysis_data/SBJ04651/tso_ctdna_tumor_only",
#'   "20240223d1951163/L2400183/Results"
#' )
#'
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
#'   outdir = file.path(p, "dracarys_tidy"),
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
    #' @description Tidy given files.
    #' @param x Tibble with `localpath` to file and the function `type` to parse it.
    tidy_files = function(x) {
      assertthat::assert_that(is.data.frame(x))
      assertthat::assert_that(all(c("type", "localpath") %in% colnames(x)))
      d1 <- x |>
        dplyr::filter(.data$type != "DOWNLOAD_ONLY") |>
        dplyr::rowwise() |>
        dplyr::mutate(
          data = list(dr_func_eval(f = .data$type, v = .data$type, envir = self)(.data$localpath))
        ) |>
        dplyr::ungroup()
      d1 |>
        dplyr::select("data") |>
        tidyr::unnest("data")
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
      dat <- TsoMergedSmallVariantsVcfFile$new(x)$read()
      tibble::tibble(name = "mergedsmallv", data = list(dat))
    },
    #' @description Read `MergedSmallVariants.genome.vcf.gz` file.
    #' @param x Path to file.
    read_msvg = function(x) {
      dat <- TsoMergedSmallVariantsGenomeVcfFile$new(x)$read()
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
      dat <- TsoCopyNumberVariantsVcfFile$new(x)$read()
      tibble::tibble(name = "cnv", data = list(dat))
    },
    #' @description Read `fragment_length_hist.json.gz` file.
    #' @param x Path to file.
    read_flh = function(x) {
      dat <- TsoFragmentLengthHistFile$new(x)$read()
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
      dat <- TsoTmbFile$new(x)$read()
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
      dat <- TsoFusionsCsvFile$new(x)$read()
      tibble::tibble(name = "fusions", data = list(dat))
    }
  ) # end public
)

#' TsoCombinedVariantOutputFile R6 Class
#'
#' @description
#' Contains methods for reading and displaying contents of the
#' `CombinedVariantOutput.tsv` file output from TSO.
#'
#' @examples
#' \dontrun{
#' x <- "/path/to/CombinedVariantOutput.tsv"
#' d <- TsoCombinedVariantOutputFile$new(x)
#' d$read()
#' }
#' @export
TsoCombinedVariantOutputFile <- R6::R6Class(
  "TsoCombinedVariantOutputFile",
  inherit = File,
  public = list(
    #' @description
    #' Reads the `CombinedVariantOutput.tsv` file output from TSO and extracts
    #' only the Small Variants section (due to inconsistencies with other sections).
    #' @return tibble with variants.
    read = function() {
      x <- self$path
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
    },
    #' @description
    #' Writes a tidy version of the `CombinedVariantOutput.tsv` (only Small Variants)
    #' file output from TSO.
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
      d <- readr::read_tsv(x, col_types = ct)
      d[]
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
      cnames <- c("FragmentLength", "Count")
      # handle SBJ00006...
      if (length(j) == 0) {
        return(empty_tbl(cnames = cnames))
      }
      assertthat::assert_that(
        all(names(j[[1]] %in% cnames))
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
      return(res[])
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
