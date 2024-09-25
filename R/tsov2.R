#' Wf_tso_ctdna_tumor_only_v2 R6 Class
#'
#' @description
#' Reads and writes tidy versions of files from the `tso_ctdna_tumor_only_v2` workflow.
#'
#' @examples
#' \dontrun{
#'
#' #---- Local ----#
#' p <- file.path(
#'   "~/s3/pipeline-prod-cache-503977275616-ap-southeast-2/byob-icav2/production",
#'   "analysis/cttsov2/20240915ff0295ed",
#'   "Results"
#' )
#' LibraryID <- "L2401290"
#' t1 <- Wf_tso_ctdna_tumor_only_v2$new(path = p, LibraryID = LibraryID)
#' t1$list_files(max_files = 20)
#' t1$list_files_filter_relevant(max_files = 300)
#' d <- t1$download_files(max_files = 100, dryrun = F)
#' d_tidy <- t1$tidy_files(d)
#' d_write <- t1$write(
#'   d_tidy,
#'   outdir = file.path(p, "dracarys_tidy"),
#'   prefix = LibraryID,
#'   format = "tsv"
#' )
#'
#' #---- S3 ----#
#' p <- file.path(
#'   "s3://pipeline-prod-cache-503977275616-ap-southeast-2/byob-icav2/production",
#'   "analysis/cttsov2/20240915ff0295ed"
#' )
#' LibraryID <- "L2401290"
#' outdir <- sub("s3:/", "~/s3", p)
#' t2 <- Wf_tso_ctdna_tumor_only_v2$new(path = p, LibraryID = LibraryID)
#' t2$list_files(max_files = 500)
#' t2$list_files_filter_relevant(max_files = 500)
#' d <- t2$download_files(
#'   outdir = outdir,
#'   max_files = 500,
#'   dryrun = F
#' )
#' d_tidy <- t2$tidy_files(d)
#' d_write <- t2$write(
#'   d_tidy,
#'   outdir = file.path(outdir, "dracarys_tidy"),
#'   prefix = LibraryID,
#'   format = "tsv"
#' )
#' }
#' @export
Wf_tso_ctdna_tumor_only_v2 <- R6::R6Class(
  "Wf_tso_ctdna_tumor_only_v2",
  inherit = Wf,
  public = list(
    #' @field LibraryID The LibraryID of the tumor sample (needed for path lookup).
    LibraryID = NULL,
    #' @description Create a new Wf_tso_ctdna_tumor_only_v2 object.
    #' @param path Path to directory with raw workflow results (from S3 or
    #' local filesystem).
    #' @param LibraryID The LibraryID of the sample (needed for path lookup).
    initialize = function(path = NULL, LibraryID = NULL) {
      wname <- "tso_ctdna_tumor_only_v2"
      pref <- LibraryID
      res <- glue("Results/{pref}")
      li <- "Logs_Intermediates"
      dc <- glue("{li}/DragenCaller/{pref}")
      # Results
      reg1 <- tibble::tribble(
        ~regex, ~fun,
        glue("{res}/{pref}\\.cnv\\.vcf\\.gz$"), "cnv",
        glue("{res}/{pref}\\.cnv\\.vcf\\.gz\\.tbi$"), "DOWNLOAD_ONLY",
        glue("{res}/{pref}\\.exon_cov_report\\.tsv$"), "cvgrepe",
        glue("{res}/{pref}\\.gene_cov_report\\.tsv$"), "cvgrepg",
        glue("{res}/{pref}\\.hard-filtered\\.vcf\\.gz$"), "hardfilt",
        glue("{res}/{pref}\\.hard-filtered\\.vcf\\.gz\\.tbi$"), "DOWNLOAD_ONLY",
        glue("{res}/{pref}\\.microsat_output\\.json$"), "msi",
        glue("{res}/{pref}\\.tmb.trace\\.tsv$"), "tmbt",
        glue("{res}/{pref}_CombinedVariantOutput\\.tsv$"), "cvo",
        glue("{res}/{pref}_Fusions\\.csv$"), "fus",
        glue("{res}/{pref}_MetricsOutput\\.tsv$"), "DOWNLOAD_ONLY",
        glue("{res}/{pref}_SmallVariants_Annotated\\.json\\.gz$"), "DOWNLOAD_ONLY",
        glue("{li}/SampleAnalysisResults/{pref}_SampleAnalysisResults\\.json$"), "sar"
      )
      # DragenCaller
      reg2 <- tibble::tribble(
        ~regex, ~fun,
        glue("{dc}/{pref}\\-replay\\.json$"), "DOWNLOAD_ONLY",
        glue("{dc}/{pref}\\.cnv_metrics.csv$"), "DOWNLOAD_ONLY",
        glue("{dc}/{pref}\\.exon_contig_mean_cov\\.csv$"), "DOWNLOAD_ONLY",
        glue("{dc}/{pref}\\.exon_coverage_metrics\\.csv$"), "DOWNLOAD_ONLY",
        glue("{dc}/{pref}\\.exon_fine_hist\\.csv$"), "DOWNLOAD_ONLY",
        glue("{dc}/{pref}\\.exon_hist\\.csv$"), "DOWNLOAD_ONLY",
        glue("{dc}/{pref}\\.exon_overall_mean_cov\\.csv$"), "DOWNLOAD_ONLY",
        glue("{dc}/{pref}\\.fastqc_metrics\\.csv$"), "DOWNLOAD_ONLY",
        glue("{dc}/{pref}\\.fragment_length_hist\\.csv$"), "DOWNLOAD_ONLY",
        glue("{dc}/{pref}\\.gc_metrics\\.csv$"), "DOWNLOAD_ONLY",
        glue("{dc}/{pref}\\.gvcf_metrics\\.csv$"), "DOWNLOAD_ONLY",
        glue("{dc}/{pref}\\.mapping_metrics\\.csv$"), "DOWNLOAD_ONLY",
        glue("{dc}/{pref}\\.microsat_diffs\\.txt$"), "DOWNLOAD_ONLY",
        glue("{dc}/{pref}\\.microsat_output\\.json$"), "DOWNLOAD_ONLY",
        glue("{dc}/{pref}\\.sv_metrics\\.csv$"), "DOWNLOAD_ONLY",
        glue("{dc}/{pref}\\.target_bed_contig_mean_cov\\.csv$"), "DOWNLOAD_ONLY",
        glue("{dc}/{pref}\\.target_bed_coverage_metrics\\.csv$"), "DOWNLOAD_ONLY",
        glue("{dc}/{pref}\\.target_bed_fine_hist\\.csv$"), "DOWNLOAD_ONLY",
        glue("{dc}/{pref}\\.target_bed_hist\\.csv$"), "DOWNLOAD_ONLY",
        glue("{dc}/{pref}\\.target_bed_overall_mean_cov\\.csv$"), "DOWNLOAD_ONLY",
        glue("{dc}/{pref}\\.time_metrics\\.csv$"), "DOWNLOAD_ONLY",
        glue("{dc}/{pref}\\.tmb_contig_mean_cov\\.csv$"), "DOWNLOAD_ONLY",
        glue("{dc}/{pref}\\.tmb_coverage_metrics\\.csv$"), "DOWNLOAD_ONLY",
        glue("{dc}/{pref}\\.tmb_fine_hist\\.csv$"), "DOWNLOAD_ONLY",
        glue("{dc}/{pref}\\.tmb_hist\\.csv$"), "DOWNLOAD_ONLY",
        glue("{dc}/{pref}\\.tmb_overall_mean_cov\\.csv$"), "DOWNLOAD_ONLY",
        glue("{dc}/{pref}\\.trimmer_metrics\\.csv$"), "DOWNLOAD_ONLY",
        glue("{dc}/{pref}\\.umi_metrics\\.csv$"), "DOWNLOAD_ONLY",
        glue("{dc}/{pref}\\.vc_metrics\\.csv$"), "DOWNLOAD_ONLY",
        glue("{dc}/{pref}\\.wgs_contig_mean_cov\\.csv$"), "DOWNLOAD_ONLY",
        glue("{dc}/{pref}\\.wgs_coverage_metrics\\.csv$"), "DOWNLOAD_ONLY",
        glue("{dc}/{pref}\\.wgs_fine_hist\\.csv$"), "DOWNLOAD_ONLY",
        glue("{dc}/{pref}\\.wgs_hist\\.csv$"), "DOWNLOAD_ONLY",
        glue("{dc}/{pref}\\.wgs_overall_mean_cov\\.csv$"), "DOWNLOAD_ONLY"
      )
      regexes <- dplyr::bind_rows(reg1, reg2) |>
        dplyr::mutate(
          fun = paste0("read_", .data$fun),
          fun = ifelse(.data$fun == "read_DOWNLOAD_ONLY", "DOWNLOAD_ONLY", .data$fun)
        )

      super$initialize(path = path, wname = wname, regexes = regexes)
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
    #' @description Read `cnv.vcf` file.
    #' @param x Path to file.
    read_cnv = function(x) {
      dat <- TsoCopyNumberVariantsVcfFile$new(x)$read(only_pass = FALSE, alias = FALSE)
      tibble::tibble(name = "cnv", data = list(dat))
    },
    #' @description Read `exon_cov_report.tsv` file.
    #' @param x Path to file.
    read_cvgrepe = function(x) {
      ctypes <- list(
        chr = "c", start = "d", end = "d", gene = "c",
        mean_coverage = "d", median_coverage = "d",
        min_coverage = "d", max_coverage = "d"
      )
      dat <- x |>
        readr::read_tsv(col_types = ctypes)
      tibble::tibble(name = "cvgrepe", data = list(dat[]))
    },
    #' @description Read `gene_cov_report.tsv` file.
    #' @param x Path to file.
    read_cvgrepg = function(x) {
      ctypes <- list(
        chr = "c", start = "d", end = "d", gene = "c",
        mean_coverage = "d", min_coverage = "d", max_coverage = "d"
      )
      dat <- x |>
        readr::read_tsv(col_types = ctypes)
      tibble::tibble(name = "cvgrepg", data = list(dat[]))
    },
    #' @description Read `hard-filtered.vcf` file.
    #' @param x Path to file.
    read_hardfilt = function(x) {
      dat <- bcftools_parse_vcf(x, only_pass = FALSE, alias = FALSE)
      tibble::tibble(name = "hardfilt", data = list(dat))
    },
    #' @description Read `microsat_output.json` file.
    #' @param x Path to file.
    read_msi = function(x) {
      dat <- TsoMsiFile$new(x)$read()
      tibble::tibble(name = "msi", data = list(dat))
    },
    #' @description Read `TMB_Trace.tsv` file.
    #' @param x Path to file.
    read_tmbt = function(x) {
      ctypes <- list(
        Chromosome = "c", Position = "d", RefCall = "c", AltCall = "c",
        VAF = "d", Depth = "d", CytoBand = "c", GeneName = "c",
        VariantType = "c", CosmicIDs = "c", MaxCosmicCount = "d",
        ClinVarIDs = "c", ClinVarSignificance = "c", AlleleCountsGnomadExome = "d",
        AlleleCountsGnomadGenome = "d", AlleleCounts1000Genomes = "d",
        MaxDatabaseAlleleCounts = "d", GermlineFilterDatabase = "c",
        GermlineFilterProxi = "c", Nonsynonymous = "c", withinValidTmbRegion = "c",
        IncludedInTMBNumerator = "c", Status = "c", ProteinChange = "c",
        CDSChange = "c", Exons = "c", Consequence = "c"
      )
      dat <- readr::read_tsv(x, col_types = ctypes)[]
      tibble::tibble(name = "tmbtrace", data = list(dat))
    },
    #' @description Read `CombinedVariantOutput.tsv` file.
    #' @param x Path to file.
    read_cvo = function(x) {
      dat <- TsoCombinedVariantOutputFile$new(x)$read()
      tibble::tibble(name = "combinedvaro", data = list(dat))
    },
    #' @description Read `Fusions.csv` file.
    #' @param x Path to file.
    read_fus = function(x) {
      ct <- readr::cols(
        Sample = "c", Name = "c", Chr1 = "c", Pos1 = "d", Chr2 = "c",
        Pos2 = "d", Direction = "c", Alt_Depth = "d", BP1_Depth = "d",
        BP2_Depth = "d", Total_Depth = "d", VAF = "d", Gene1 = "c", Gene2 = "c",
        Contig = "c", Filter = "c", Is_Cosmic_GenePair = "l",
        "Fusion Directionality Known" = "c"
      )
      dat <- readr::read_csv(x, col_types = ct, comment = "#")
      tibble::tibble(name = "fusions", data = list(dat[]))
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
    }
  ) # end public
)
