#' Wf_cttsov2 R6 Class
#'
#' @description
#' Reads and writes tidy versions of files from the `cttsov2` workflow.
#'
#' @examples
#' \dontrun{
#'
#' #---- Local ----#
#' p <- file.path(
#'   "~/s3/pipeline-prod-cache-503977275616-ap-southeast-2/byob-icav2/production",
#'   "analysis/cttsov2/20250117c5b9baa8"
#' )
#' prefix <- "L2500039"
#' t1 <- Wf_cttsov2$new(path = p, prefix = prefix)
#' t1$list_files(max_files = 100)
#' t1$dragenObj$list_files(max_files = 100)
#' t1$list_files_filter_relevant(max_files = 300)
#' t1$dragenObj$list_files_filter_relevant(max_files = 300)
#' d1 <- t1$download_files(max_files = 100, dryrun = F)
#' d2 <- t1$dragenObj$download_files(max_files = 100, dryrun = F)
#' d1_tidy <- t1$tidy_files(d1)
#' d2_tidy <- t1$dragenObj$tidy_files(d2)
#' d_write1 <- t1$write(
#'   d1_tidy,
#'   outdir = file.path(p, "dracarys_tidy"),
#'   prefix = prefix,
#'   format = "tsv"
#' )
#' d_write2 <- t1$dragenObj$write(
#'   d2_tidy,
#'   outdir = file.path(p, "dracarys_tidy"),
#'   prefix = prefix,
#'   format = "tsv"
#' )
#'
#' #---- S3 ----#
#' p <- file.path(
#'   "s3://pipeline-prod-cache-503977275616-ap-southeast-2/byob-icav2/production",
#'   "analysis/cttsov2/20250117c5b9baa8"
#' )
#' prefix <- "L2500039"
#' outdir <- sub("s3:/", "~/s3", p)
#' t2 <- Wf_cttsov2$new(path = p, prefix = prefix)
#' t2$list_files(max_files = 500)
#' t2$list_files_filter_relevant(max_files = 500)
#' t2$dragenObj$list_files_filter_relevant(max_files = 300)
#' d1 <- t2$download_files(
#'   outdir = outdir,
#'   max_files = 500,
#'   dryrun = FALSE
#' )
#' d2 <- t2$dragenObj$download_files(
#'   outdir = outdir,
#'   max_files = 500,
#'   dryrun = FALSE
#' )
#' d_tidy1 <- t2$tidy_files(d1)
#' d_tidy2 <- t2$dragenObj$tidy_files(d2)
#' d_write <- t2$write(
#'   d_tidy,
#'   outdir = file.path(outdir, "dracarys_tidy"),
#'   prefix = prefix,
#'   format = "tsv"
#' )
#' }
#' @export
Wf_cttsov2 <- R6::R6Class(
  "Wf_cttsov2",
  inherit = Wf,
  public = list(
    #' @field prefix The LibraryID prefix of the tumor sample (needed for path lookup).
    prefix = NULL,
    #' @field dragenObj dragen object.
    dragenObj = NULL,
    #' @description Create a new Wf_cttsov2 object.
    #' @param path Path to directory with raw workflow results (from S3 or
    #' local filesystem).
    #' @param prefix The LibraryID prefix of the tumor sample (needed for path lookup).
    initialize = function(path = NULL, prefix = NULL) {
      wname <- "cttsov2"
      pref <- prefix
      res <- glue("Results/{pref}")
      li <- "Logs_Intermediates"
      dc <- glue("{li}/DragenCaller/{pref}")
      path <- sub("/$", "", path) # remove potential trailing slash
      self$dragenObj <- Wf_dragen$new(
        path = file.path(path, dc),
        prefix = glue("{dc}/{prefix}")
      )
      # Results
      # fmt: skip
      regexes <- tibble::tribble(
        ~regex, ~fun,
        glue("{res}/{pref}\\.cnv\\.vcf\\.gz$"), "read_cnv",
        glue("{res}/{pref}\\.cnv\\.vcf\\.gz\\.tbi$"), "DOWNLOAD_ONLY-cnvvcfi",
        glue("{res}/{pref}\\.exon_cov_report\\.tsv$"), "read_cvgrepe",
        glue("{res}/{pref}\\.gene_cov_report\\.tsv$"), "read_cvgrepg",
        glue("{res}/{pref}\\.hard-filtered\\.vcf\\.gz$"), "read_hardfilt",
        glue("{res}/{pref}\\.hard-filtered\\.vcf\\.gz\\.tbi$"), "DOWNLOAD_ONLY-hardfiltvcfi",
        # glue("{res}/{pref}\\.microsat_output\\.json$"), "msi", # in DragenCaller
        glue("{res}/{pref}\\.tmb.trace\\.tsv$"), "read_tmbt",
        glue("{res}/{pref}_CombinedVariantOutput\\.tsv$"), "read_cvo",
        glue("{res}/{pref}_Fusions\\.csv$"), "read_fus",
        glue("{res}/{pref}_MetricsOutput\\.tsv$"), "DOWNLOAD_ONLY-metricsoutput",
        # glue("{res}/{pref}_SmallVariants_Annotated\\.json\\.gz$"), "DOWNLOAD_ONLY-smallvannjson",
        glue("{li}/SampleAnalysisResults/{pref}_SampleAnalysisResults\\.json$"), "read_sar"
      )
      super$initialize(path = path, wname = wname, regexes = regexes)
      self$prefix <- prefix
    },
    #' @description Print details about the Workflow.
    #' @param ... (ignored).
    print = function(...) {
      # fmt: skip
      res <- tibble::tribble(
        ~var, ~value,
        "path", private$.path,
        "wname", private$.wname,
        "filesystem", private$.filesystem,
        "prefix", self$prefix
      )
      print(res)
      invisible(self)
    },
    #' @description Read `SampleAnalysisResults.json.gz` file.
    #' @param x Path to file.
    read_sar = function(x) {
      tso_sar_read(x)
    },
    #' @description Read `TMB_Trace.tsv` file.
    #' @param x Path to file.
    read_tmbt = function(x) {
      dat <- tso_tmbt_read(x)
      tibble::tibble(name = "tmbtrace", data = list(dat))
    },
    #' @description Read `CombinedVariantOutput.tsv` file.
    #' @param x Path to file.
    read_cvo = function(x) {
      dat <- tso_combinedvaro_smallv_read(x)
      tibble::tibble(name = "combinedvaro", data = list(dat))
    },
    #' @description Read `cnv.vcf` file.
    #' @param x Path to file.
    read_cnv = function(x) {
      dat <- bcftools_parse_vcf(x, only_pass = FALSE, alias = TRUE)
      tibble::tibble(name = "cnv", data = list(dat))
    },
    #' @description Read `exon_cov_report.tsv` file.
    #' @param x Path to file.
    read_cvgrepe = function(x) {
      ctypes <- list(
        chr = "c",
        start = "d",
        end = "d",
        gene = "c",
        mean_coverage = "d",
        median_coverage = "d",
        min_coverage = "d",
        max_coverage = "d"
      )
      dat <- x |>
        readr::read_tsv(col_types = ctypes)
      tibble::tibble(name = "cvgrepe", data = list(dat[]))
    },
    #' @description Read `gene_cov_report.tsv` file.
    #' @param x Path to file.
    read_cvgrepg = function(x) {
      ctypes <- list(
        chr = "c",
        start = "d",
        end = "d",
        gene = "c",
        mean_coverage = "d",
        min_coverage = "d",
        max_coverage = "d"
      )
      dat <- x |>
        readr::read_tsv(col_types = ctypes)
      tibble::tibble(name = "cvgrepg", data = list(dat[]))
    },
    #' @description Read `hard-filtered.vcf` file.
    #' @param x Path to file.
    read_hardfilt = function(x) {
      dat <- bcftools_parse_vcf(x, only_pass = FALSE, alias = TRUE)
      tibble::tibble(name = "hardfilt", data = list(dat))
    },
    #' @description Read `microsat_output.json` file.
    #' @param x Path to file.
    read_msi = function(x) {
      dat <- tso_msi_read(x)
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

#' Wf_cttsov2 Download Tidy and Write
#'
#' Downloads files from the `cttsov2` workflow and writes them in a tidy format.
#'
#' @param path Path to directory with raw workflow results (S3 or
#' local filesystem).
#' @param prefix The LibraryID prefix of the sample (needed for path lookup).
#' @param outdir Path to output directory with raw files.
#' @param outdir_tidy Path to output directory with tidy files.
#' @param format Format of output files.
#' @param max_files Max number of files to list.
#' @param dryrun If TRUE, just list the files that will be downloaded (don't
#' download them).
#' @return Tibble of tidy tibbles.
#'
#' @examples
#' \dontrun{
#' #---- Local ----#
#' p <- file.path(
#'   "~/s3/pipeline-prod-cache-503977275616-ap-southeast-2/byob-icav2/production",
#'   "analysis/cttsov2/20240915ff0295ed"
#' )
#' prefix <- "L2401290"
#' outdir <- file.path(p, "dracarys_tidy")
#' d <- dtw_Wf_cttsov2(
#'   path = p, prefix = prefix, outdir = outdir,
#'   format = "tsv",
#'   dryrun = F
#' )
#'
#' p <- file.path(
#'   "s3://pipeline-prod-cache-503977275616-ap-southeast-2/byob-icav2/production",
#'   "analysis/cttsov2/20240915ff0295ed"
#' )
#' prefix <- "L2401290"
#' outdir <- sub("s3:/", "~/s3", p)
#' d <- dtw_Wf_cttsov2(
#'   path = p, prefix = prefix, outdir = outdir,
#'   format = "tsv",
#'   dryrun = F
#' )
#' }
#' @export
dtw_Wf_cttsov2 <- function(
  path,
  prefix,
  outdir,
  outdir_tidy = file.path(outdir, "dracarys_tidy"),
  format = "rds",
  max_files = 1000,
  dryrun = FALSE
) {
  obj <- Wf_cttsov2$new(path = path, prefix = prefix)
  d_dl1 <- obj$download_files(
    outdir = outdir,
    max_files = max_files,
    dryrun = dryrun
  )
  d_dl2 <- obj$dragenObj$download_files(
    outdir = outdir,
    max_files = max_files,
    dryrun = dryrun
  )
  if (!dryrun) {
    d_tidy1 <- obj$tidy_files(d_dl1)
    d_tidy2 <- obj$dragenObj$tidy_files(d_dl2)
    d_write1 <- obj$write(
      d_tidy1,
      outdir = outdir_tidy,
      prefix = prefix,
      format = format
    )
    d_write2 <- obj$dragenObj$write(
      d_tidy2,
      outdir = outdir_tidy,
      prefix = prefix,
      format = format
    )
    return(dplyr::bind_rows(d_write1, d_write2))
  }
  return(dplyr::bind_rows(d_dl1, d_dl2))
}
