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
#'   "analysis/cttsov2/20240922ced7560e",
#'   "Results"
#' )
#' LibraryID <- "L2401321"
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
#'   "analysis/cttsov2/20240922ced7560e",
#'   "Results"
#' )
#' LibraryID <- "L2401321"
#' outdir <- sub("s3:/", "~/s3", p)
#' t2 <- Wf_tso_ctdna_tumor_only_v2$new(path = p, LibraryID = LibraryID)
#' t2$list_files(max_files = 100)
#' t2$list_files_filter_relevant(max_files = 100)
#' d <- t2$download_files(
#'   outdir = outdir,
#'   max_files = 100,
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
      regexes <- tibble::tribble(
        ~regex, ~fun,
        glue("{pref}\\.cnv\\.vcf$"), "cnv",
        glue("{pref}\\.exon_cov_report\\.tsv$"), "cvgrepe",
        glue("{pref}\\.gene_cov_report\\.tsv$"), "cvgrepg",
        glue("{pref}\\.hard-filtered\\.vcf$"), "hardfilt",
        glue("{pref}\\.microsat_output\\.json$"), "msi",
        glue("{pref}\\.tmb.trace\\.tsv$"), "tmbt",
        glue("{pref}_CombinedVariantOutput\\.tsv$"), "cvo",
        glue("{pref}_Fusions\\.csv$"), "fus"
      ) |>
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
      ct <- readr::cols(
        Sample = "c", Name = "c", Chr1 = "c", Pos1 = "d", Chr2 = "c",
        Pos2 = "d", Direction = "c", Alt_Depth = "d", BP1_Depth = "d",
        BP2_Depth = "d", Total_Depth = "d", VAF = "d", Gene1 = "c", Gene2 = "c",
        Contig = "c", Filter = "c", Is_Cosmic_GenePair = "l",
        "Fusion Directionality Known" = "c"
      )
      dat <- readr::read_csv(x, col_types = ct, comment = "#")
      tibble::tibble(name = "fusions", data = list(dat))
    },
    #' @description Read `hard-filtered.vcf` file.
    #' @param x Path to file.
    read_hardfilt = function(x) {
      dat <- bcftools_parse_vcf(x, only_pass = FALSE, alias = FALSE)
      tibble::tibble(name = "hardfilt", data = list(dat))
    }
  ) # end public
)
