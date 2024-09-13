#' Wf_sash R6 Class
#'
#' @description
#' Reads and writes tidy versions of files from the `sash` workflow
#'
#' @examples
#' \dontrun{
#'
#' #---- Local ----#
#' p1 <- "~/s3/org.umccr.data.oncoanalyser/analysis_data/SBJ05571/sash"
#' p2 <- "202408270b93455e/L2401308_L2401307"
#' p <- normalizePath(file.path(p1, p2))
#' SubjectID <- "SBJ05571"
#' SampleID_tumor <- "MDX240307"
#' prefix <- glue("{SubjectID}__{SampleID_tumor}")
#' s1 <- Wf_sash$new(path = p, SubjectID = SubjectID, SampleID_tumor = SampleID_tumor)
#' s1$list_files(max_files = 20)
#' s1$list_files_filter_relevant(max_files = 300)
#' d <- s1$download_files(max_files = 1000, dryrun = F)
#' d_tidy <- s1$tidy_files(d)
#' d_write <- s1$write(
#'   d_tidy,
#'   outdir = file.path(p, "dracarys_tidy"),
#'   prefix = glue("{SubjectID}_{SampleID_tumor}"),
#'   format = "tsv"
#' )
#'
#' #---- S3 ----#
#' p1 <- "s3://org.umccr.data.oncoanalyser/analysis_data/SBJ05571/sash"
#' p2 <- "202408270b93455e/L2401308_L2401307"
#' p <- file.path(p1, p2)
#' SubjectID <- "SBJ05571"
#' SampleID_tumor <- "MDX240307"
#' prefix <- glue("{SubjectID}__{SampleID_tumor}")
#' s1 <- Wf_sash$new(path = p, SubjectID = SubjectID, SampleID_tumor = SampleID_tumor)
#' s1$list_files(max_files = 20)
#' s1$list_files_filter_relevant()
#' outdir <- sub("s3:/", "~/s3", p)
#' d <- s1$download_files(outdir = outdir, max_files = 1000, dryrun = F)
#' d_tidy <- s1$tidy_files(d)
#' d_write <- s1$write(
#'   d_tidy,
#'   outdir = file.path(p, "dracarys_tidy"),
#'   prefix = glue("{SubjectID}__{SampleID_tumor}"),
#'   format = "tsv"
#' )
#' }
#'
#' @export
Wf_sash <- R6::R6Class(
  "Wf_sash",
  inherit = Wf,
  public = list(
    #' @field SubjectID The SubjectID of the sample (needed for path lookup).
    #' @field SampleID_tumor The SampleID of the tumor sample (needed for path lookup).
    SubjectID = NULL,
    SampleID_tumor = NULL,
    #' @description Create a new Wf_sash object.
    #' @param path Path to directory with raw workflow results (from GDS, S3, or
    #' local filesystem).
    #' @param SubjectID The SubjectID of the sample (needed for path lookup).
    #' @param SampleID_tumor The SampleID of the tumor sample (needed for path lookup).
    initialize = function(path = NULL, SubjectID = NULL, SampleID_tumor = NULL) {
      wname <- "sash"
      regexes <- tibble::tribble(
        ~regex, ~fun,
        "cancer_report/cancer_report_tables/hrd/.*-chord\\.tsv\\.gz$", "hrdchordtsv",
        "cancer_report/cancer_report_tables/hrd/.*-hrdetect\\.tsv\\.gz$", "hrdetecttsv",
        "cancer_report/cancer_report_tables/hrd/.*-dragen\\.tsv\\.gz$", "hrddragentsv",
        "cancer_report/cancer_report_tables/sigs/.*-snv_2015\\.tsv\\.gz$", "sigssnv2015tsv",
        "cancer_report/cancer_report_tables/sigs/.*-snv_2020\\.tsv\\.gz$", "sigssnv2020tsv",
        "cancer_report/cancer_report_tables/sigs/.*-dbs\\.tsv\\.gz$", "sigsdbstsv",
        "cancer_report/cancer_report_tables/sigs/.*-indel\\.tsv\\.gz$", "sigsindeltsv",
        "cancer_report/cancer_report_tables/.*-qc_summary\\.tsv\\.gz$", "qcsummarytsv",
        "smlv_somatic/report/pcgr/.*\\.pcgr_acmg\\.grch38\\.json\\.gz$", "pcgrjson"
      ) |>
        dplyr::mutate(fun = paste0("read_", .data$fun))

      super$initialize(path = path, wname = wname, regexes = regexes)
      self$SubjectID <- SubjectID
      self$SampleID_tumor <- SampleID_tumor
    },
    #' @description Print details about the Workflow.
    #' @param ... (ignored).
    print = function(...) {
      res <- tibble::tribble(
        ~var, ~value,
        "path", self$path,
        "wname", self$wname,
        "filesystem", self$filesystem,
        "SubjectID", self$SubjectID,
        "SampleID_tumor", self$SampleID_tumor
      )
      print(res)
      invisible(self)
    },
    #' @description List dracarys files under given path
    #' @param max_files Max number of files to list (for gds/s3 only).
    #' @param ica_token ICA access token (def: $ICA_ACCESS_TOKEN env var).
    #' @param ... Passed on to the `gds_list_files_filter_relevant` or
    #' the `s3_list_files_filter_relevant` function.
    list_files_filter_relevant = function(max_files = 1000,
                                          ica_token = Sys.getenv("ICA_ACCESS_TOKEN"), ...) {
      path <- self$path
      dir1 <- file.path(path, glue("{self$SubjectID}_{self$SampleID_tumor}"))
      f1 <- super$list_files_filter_relevant(path = dir1, max_files = 500)
      return(f1)
    },
    #' @description Download files from GDS/S3 to local filesystem.
    #' @param outdir Path to output directory.
    #' @param ica_token ICA access token (def: $ICA_ACCESS_TOKEN env var).
    #' @param max_files Max number of files to list.
    #' @param dryrun If TRUE, just list the files that will be downloaded (don't
    #' download them).
    #' @param recursive Should files be returned recursively _in and under_ the specified
    #' GDS directory, or _only directly in_ the specified GDS directory (def: TRUE via ICA API).
    download_files = function(outdir, ica_token = Sys.getenv("ICA_ACCESS_TOKEN"),
                              max_files = 1000, dryrun = FALSE, recursive = NULL) {
      super$download_files(
        outdir = outdir, ica_token = ica_token, max_files = max_files,
        dryrun = dryrun, recursive = recursive,
        list_filter_fun = self$list_files_filter_relevant
      )
    },
    #' @description Read `pcgr.json.gz` file.
    #' @param x Path to file.
    read_pcgrjson = function(x) {
      j <- read_jsongz_jsonlite(x)
      tmb <-
        j[["content"]][["tmb"]][["variant_statistic"]] %||%
        j[["content"]][["tmb"]][["v_stat"]] %||%
        list(tmb_estimate = NA, n_tmb = NA)
      tmb <- purrr::flatten(tmb) |>
        tibble::as_tibble_row() |>
        dplyr::select("tmb_estimate", "n_tmb")
      msi <- j[["content"]][["msi"]][["prediction"]][["msi_stats"]]
      # handle nulls
      msi <- msi %||% list(fracIndels = NA, predicted_class = NA)
      msi <- purrr::flatten(msi) |>
        tibble::as_tibble_row() |>
        dplyr::select("fracIndels", "predicted_class")
      metrics <- dplyr::bind_cols(msi, tmb)
      return(metrics)
    },
    #' @description Read `dragen.tsv.gz` cancer report hrd file.
    #' @param x Path to file.
    read_hrddragentsv = function(x) {
      ct <- readr::cols(.default = "d", Sample = "c")
      read_tsvgz(x, col_types = ct)
    },
    #' @description Read `chord.tsv.gz` cancer report hrd file.
    #' @param x Path to file.
    read_hrdchordtsv = function(x) {
      ct <- readr::cols_only(
        p_hrd = "d",
        hr_status = "c",
        hrd_type = "c",
        p_BRCA1 = "d",
        p_BRCA2 = "d"
      )
      read_tsvgz(x, col_types = ct)
    },
    #' @description Read `hrdetect.tsv.gz` cancer report hrd file.
    #' @param x Path to file.
    read_hrdetecttsv = function(x) {
      ct <- readr::cols(
        .default = "d",
        sample = "c"
      )
      read_tsvgz(x, col_types = ct) |>
        dplyr::select(-c("sample"))
    },
    #' @description Read signature cancer report file.
    #' @param x Path to file.
    read_sigstsv = function(x) {
      ct <- readr::cols(
        .default = "d",
        Signature = "c"
      )
      read_tsvgz(x, col_types = ct)
    },
    #' @description Read `snv_2015.tsv.gz` sigs cancer report file.
    #' @param x Path to file.
    read_sigssnv2015tsv = function(x) {
      self$read_sigstsv(x)
    },
    #' @description Read `snv_2020.tsv.gz` sigs cancer report file.
    #' @param x Path to file.
    read_sigssnv2020tsv = function(x) {
      self$read_sigstsv(x)
    },
    #' @description Read `dbs.tsv.gz` sigs cancer report file.
    #' @param x Path to file.
    read_sigsdbstsv = function(x) {
      self$read_sigstsv(x)
    },
    #' @description Read `indel.tsv.gz` sigs cancer report file.
    #' @param x Path to file.
    read_sigsindeltsv = function(x) {
      self$read_sigstsv(x)
    },
    #' @description Read `qc_summary.tsv.gz` cancer report file.
    #' @param x Path to file.
    read_qcsummarytsv = function(x) {
      d <- read_tsvgz(x, col_types = readr::cols(.default = "c"))
      d |>
        dplyr::select("variable", "value") |>
        tidyr::pivot_wider(names_from = "variable", values_from = "value") |>
        dplyr::rename(MSI_mb_tmp = "MSI (indels/Mb)") |>
        dplyr::mutate(
          purity_hmf = sub("(.*) \\(.*\\)", "\\1", .data$Purity) |> as.numeric(),
          ploidy_hmf = sub("(.*) \\(.*\\)", "\\1", .data$Ploidy) |> as.numeric(),
          msi_mb_hmf = sub(".* \\((.*)\\)", "\\1", .data$MSI_mb_tmp) |> as.numeric(),
          contamination_hmf = as.numeric(.data$Contamination),
          deleted_genes_hmf = as.numeric(.data$DeletedGenes),
          msi_hmf = sub("(.*) \\(.*\\)", "\\1", .data$MSI_mb_tmp),
          tmb_hmf = sub("(.*) \\(.*\\)", "\\1", .data$TMB) |> as.numeric(),
          tml_hmf = sub("(.*) \\(.*\\)", "\\1", .data$TML) |> as.numeric(),
          hypermutated = ifelse("Hypermutated" %in% d$variable, .data[["Hypermutated"]], NA) |> as.character()
        ) |>
        dplyr::select(
          qc_status_hmf = "QC_Status",
          sex_hmf = "Gender",
          "purity_hmf", "ploidy_hmf", "msi_hmf", "msi_mb_hmf",
          "contamination_hmf",
          "deleted_genes_hmf", "tmb_hmf", "tml_hmf",
          wgd_hmf = "WGD",
          "hypermutated"
        )
    }
  ) # end public
)

#' sash Download Tidy and Write
#'
#' Downloads files from the `sash` workflow and writes them in a tidy format.
#'
#' @param path Path to directory with raw workflow results (from GDS, S3, or
#' local filesystem).
#' @param SubjectID The SubjectID of the sample (needed for path lookup).
#' @param SampleID_tumor The SampleID of the tumor sample (needed for path lookup).
#' @param outdir Path to output directory.
#' @param format Format of output files.
#' @param max_files Max number of files to list.
#' @param ica_token ICA access token (def: $ICA_ACCESS_TOKEN env var).
#' @param dryrun If TRUE, just list the files that will be downloaded (don't
#' download them).
#' @return List where each element is a tidy tibble of a sash file.
#'
#' @examples
#' \dontrun{
#' SubjectID <- "SBJ03043"
#' SampleID_tumor <- "PRJ230004"
#' p1_gds <- glue("gds://production/analysis_data/{SubjectID}/umccrise")
#' p <- file.path(p1_gds, "20240830ec648f40/L2300064__L2300063")
#' outdir <- file.path(sub("gds:/", "~/icav1/g", p))
#' token <- Sys.getenv("ICA_ACCESS_TOKEN")
#' d <- Wf_sash_download_tidy_write(
#'   path = p, SubjectID = SubjectID, SampleID_tumor = SampleID_tumor,
#'   outdir = outdir,
#'   dryrun = F
#' )
#' }
#' @export
Wf_sash_download_tidy_write <- function(path, SubjectID, SampleID_tumor,
                                        outdir, format = "rds", max_files = 1000,
                                        ica_token = Sys.getenv("ICA_ACCESS_TOKEN"),
                                        dryrun = FALSE) {
  s <- Wf_sash$new(
    path = path, SubjectID = SubjectID, SampleID_tumor = SampleID_tumor
  )
  d_dl <- s$download_files(
    outdir = outdir, ica_token = ica_token,
    max_files = max_files, dryrun = dryrun
  )
  if (!dryrun) {
    d_tidy <- s$tidy_files(d_dl)
    d_write <- s$write(
      d_tidy,
      outdir = file.path(outdir, "dracarys_tidy"),
      prefix = glue("{SubjectID}__{SampleID_tumor}"),
      format = format
    )
    return(d_write)
  }
  return(d_dl)
}
