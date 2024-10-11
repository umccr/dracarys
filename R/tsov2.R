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
#'   "analysis/cttsov2/20240915ff0295ed"
#' )
#' prefix <- "L2401290"
#' t1 <- Wf_tso_ctdna_tumor_only_v2$new(path = p, prefix = prefix)
#' t1$list_files(max_files = 100)
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
#' #---- S3 ----#
#' p <- file.path(
#'   "s3://pipeline-prod-cache-503977275616-ap-southeast-2/byob-icav2/production",
#'   "analysis/cttsov2/20240915ff0295ed"
#' )
#' prefix <- "L2401290"
#' outdir <- sub("s3:/", "~/s3", p)
#' t2 <- Wf_tso_ctdna_tumor_only_v2$new(path = p, prefix = prefix)
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
#'   prefix = prefix,
#'   format = "tsv"
#' )
#' }
#' @export
Wf_tso_ctdna_tumor_only_v2 <- R6::R6Class(
  "Wf_tso_ctdna_tumor_only_v2",
  inherit = Wf,
  public = list(
    #' @field prefix The LibraryID prefix of the tumor sample (needed for path lookup).
    prefix = NULL,
    #' @description Create a new Wf_tso_ctdna_tumor_only_v2 object.
    #' @param path Path to directory with raw workflow results (from S3 or
    #' local filesystem).
    #' @param prefix The LibraryID prefix of the tumor sample (needed for path lookup).
    initialize = function(path = NULL, prefix = NULL) {
      wname <- "tso_ctdna_tumor_only_v2"
      pref <- prefix
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
        # glue("{res}/{pref}_SmallVariants_Annotated\\.json\\.gz$"), "DOWNLOAD_ONLY",
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
      self$prefix <- prefix
    },
    #' @description Print details about the Workflow.
    #' @param ... (ignored).
    print = function(...) {
      res <- tibble::tribble(
        ~var, ~value,
        "path", self$path,
        "wname", self$wname,
        "filesystem", self$filesystem,
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

#' Wf_tso_ctdna_tumor_only_v2 Download Tidy and Write
#'
#' Downloads files from the `tso_ctdna_tumor_only_v2` workflow and writes them in a tidy format.
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
#' p <- file.path(
#'   "s3://pipeline-prod-cache-503977275616-ap-southeast-2/byob-icav2/production",
#'   "analysis/cttsov2/20240915ff0295ed"
#' )
#' prefix <- "L2401290"
#' outdir <- sub("s3:/", "~/s3", p)
#' d <- dtw_Wf_tso_ctdna_tumor_only_v2(
#'   path = p, prefix = prefix, outdir = outdir,
#'   format = "tsv",
#'   dryrun = F
#' )
#' }
#' @export
dtw_Wf_tso_ctdna_tumor_only_v2 <- function(path, prefix, outdir,
                                           outdir_tidy = file.path(outdir, "dracarys_tidy"),
                                           format = "rds",
                                           max_files = 1000,
                                           dryrun = FALSE) {
  obj <- Wf_tso_ctdna_tumor_only_v2$new(path = path, prefix = prefix)
  d_dl <- obj$download_files(
    outdir = outdir, max_files = max_files, dryrun = dryrun
  )
  if (!dryrun) {
    d_tidy <- obj$tidy_files(d_dl)
    d_write <- obj$write(
      d_tidy,
      outdir = outdir_tidy,
      prefix = prefix,
      format = format
    )
    return(d_write)
  }
  return(d_dl)
}
