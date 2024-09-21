#' @title Workflow
#'
#' @description Workflow is a base R6 class representing a bioinformatic
#' workflow run from a UMCCR workflow manager.
#'
#' A workflow has:
#'
#' - a directory path with all the raw output files (either on GDS, S3 or
#' local filesystem)
#' - a subset of files that are of interest for ingestion
#'   - tibble with full path and basename columns
#' - a set of parsers that can parse and tidy those files
#'   - each parser takes a path and returns a tidy tibble
#' - a list of tidy tibbles (or a tibble with nested tibbles)
#'
#' @examples
#' \dontrun{
#' regexes <- tibble::tribble(
#'   ~regex, ~fun,
#'   "-chord\\.tsv\\.gz$", "UmChordTsvFile",
#'   "-hrdetect\\.tsv\\.gz$", "UmHrdetectTsvFile",
#'   "-snv_2015\\.tsv\\.gz$", "UmSigsSnvFile",
#'   "-snv_2020\\.tsv\\.gz$", "UmSigsSnvFile",
#'   "-dbs\\.tsv\\.gz$", "UmSigsDbsFile",
#'   "-indel\\.tsv\\.gz$", "UmSigsIndelFile",
#'   "-qc_summary\\.tsv\\.gz$", "UmQcSumFile",
#' )
#'
#' #---- LOCAL ----#
#' p1_local <- "~/icav1/g/production/analysis_data"
#' p <- file.path(p1_local, "SBJ01155/umccrise/202408300c218043/L2101566__L2101565")
#' um1 <- Wf$new(path = p, wname = "umccrise", regexes = regexes)
#' um1$list_files(max_files = 10)
#' um1$list_files_filter_relevant(max_files = 10)
#'
#' #---- GDS ----#
#' p1_gds <- "gds://production/analysis_data"
#' p <- file.path(p1_gds, "SBJ03043/umccrise/20240830ec648f40/L2300064__L2300063")
#' outdir <- file.path(sub("gds:/", "~/icav1/g", p))
#' token <- Sys.getenv("ICA_ACCESS_TOKEN")
#' um2 <- Wf$new(path = p, wname = "umccrise", regexes = regexes)
#' um2$list_files(max_files = 10)
#' um2$list_files_filter_relevant(ica_token = token, max_files = 500)
#' d <- um2$download_files(
#'   outdir = outdir, ica_token = token,
#'   max_files = 1000, dryrun = T
#' )
#' d_tidy <- um2$tidy_files(d)
#'
#' #---- S3 ----#
#' p1_s3 <- "s3://org.umccr.data.oncoanalyser/analysis_data/SBJ05570/sash/202408275fce06c3"
#' p2_s3 <- "L2401304_L2401303/SBJ05570_MDX240299/cancer_report/cancer_report_tables"
#' p <- file.path(p1_s3, p2_s3)
#' outdir <- sub("s3:/", "~/s3", p)
#' um3 <- Wf$new(path = p, wname = "sash", regexes = regexes)
#' um3$list_files(max_files = 10)
#' um3$list_files_filter_relevant(max_files = 50)
#' d <- um3$download_files(outdir = outdir, regexes = regexes, max_files = 50, dryrun = F)
#' }
#'
#' @export
Wf <- R6::R6Class(
  "Wf",
  public = list(
    #' @field path Path to directory with raw workflow results (from GDS, S3, or
    #' local filesystem).
    #' @field wname Name of workflow (e.g. umccrise, sash).
    #' @field filesystem  Filesystem of `path` (gds/s3/local).
    #' @field regexes Tibble with file `regex` and `fun`ction to parse it.
    path = NULL,
    wname = NULL,
    filesystem = NULL,
    regexes = NULL,
    #' @description Create a new Workflow object.
    #' @param path Path to directory with raw workflow results.
    #' @param wname Name of workflow.
    #' @param regexes Tibble with file `regex` and `fun`ction to parse it.
    initialize = function(path = NULL, wname = NULL, regexes = NULL) {
      wnames <- c(
        "bcl_convert",
        "tso_ctdna_tumor_only",
        "wgs_alignment_qc",
        "wts_alignment_qc",
        "wts_tumor_only",
        "wgs_tumor_normal",
        "umccrise",
        "rnasum",
        "star_alignment",
        "oncoanalyser_wts",
        "oncoanalyser_wgs",
        "oncoanalyser_wgts_existing_both",
        "sash"
      )
      assertthat::assert_that(wname %in% wnames)
      self$path <- path
      self$wname <- wname
      self$filesystem <- dplyr::case_when(
        grepl("^gds://", path) ~ "gds",
        grepl("^s3://", path) ~ "s3",
        .default = "local"
      )
      self$regexes <- regexes
    },
    #' @description Print details about the Workflow.
    #' @param ... (ignored).
    print = function(...) {
      res <- tibble::tribble(
        ~var, ~value,
        "path", self$path,
        "wname", self$wname,
        "filesystem", self$filesystem
      )
      print(res)
      invisible(self)
    },
    #' @description List all files under given path.
    #' @param path Path with raw results.
    #' @param max_files Max number of files to list (for gds/s3 only).
    #' @param ica_token ICA access token (def: $ICA_ACCESS_TOKEN env var).
    #' @param ... Passed on to `gds_list_files_dir` function.
    list_files = function(path = self$path, max_files = 1000,
                          ica_token = Sys.getenv("ICA_ACCESS_TOKEN"), ...) {
      if (self$filesystem == "gds") {
        d <- gds_list_files_dir(
          gdsdir = path, token = ica_token, page_size = max_files, ...
        )
      } else if (self$filesystem == "s3") {
        d <- s3_list_files_dir(s3dir = path, max_objects = max_files)
      } else {
        d <- local_list_files_dir(localdir = path, max_files = max_files)
      }
      return(d)
    },
    #' @description List dracarys files under given path
    #' @param path Path with raw results.
    #' @param max_files Max number of files to list (for gds/s3 only).
    #' @param ica_token ICA access token (def: $ICA_ACCESS_TOKEN env var).
    #' @param ... Passed on to the `gds_list_files_filter_relevant` or
    #' the `s3_list_files_filter_relevant` function.
    list_files_filter_relevant = function(path = self$path, max_files = 1000,
                                          ica_token = Sys.getenv("ICA_ACCESS_TOKEN"), ...) {
      regexes <- self$regexes
      assertthat::assert_that(!is.null(regexes))
      if (self$filesystem == "gds") {
        d <- gds_list_files_filter_relevant(
          gdsdir = path, regexes = regexes, token = ica_token, page_size = max_files, ...
        )
      } else if (self$filesystem == "s3") {
        d <- s3_list_files_filter_relevant(
          s3dir = path, regexes = regexes, max_objects = max_files, ...
        )
      } else {
        d <- local_list_files_filter_relevant(
          localdir = path, regexes = regexes, max_files = max_files
        )
      }
      d
    },
    #' @description Download files from GDS/S3 to local filesystem.
    #' @param path Path with raw results.
    #' @param outdir Path to output directory.
    #' @param ica_token ICA access token (def: $ICA_ACCESS_TOKEN env var).
    #' @param max_files Max number of files to list.
    #' @param dryrun If TRUE, just list the files that will be downloaded (don't
    #' download them).
    #' @param recursive Should files be returned recursively _in and under_ the specified
    #' GDS directory, or _only directly in_ the specified GDS directory (def: TRUE via ICA API).
    download_files = function(path = self$path, outdir, ica_token = Sys.getenv("ICA_ACCESS_TOKEN"),
                              max_files = 1000, dryrun = FALSE, recursive = NULL) {
      # TODO: add envvar checker
      regexes <- self$regexes
      assertthat::assert_that(!is.null(regexes))
      if (self$filesystem == "gds") {
        d <- dr_gds_download(
          gdsdir = path, outdir = outdir, regexes = regexes, token = ica_token,
          page_size = max_files, dryrun = dryrun, recursive = recursive
        )
        if (!dryrun) {
          self$filesystem <- "local"
          self$path <- outdir
        }
      } else if (self$filesystem == "s3") {
        d <- dr_s3_download(
          s3dir = path, outdir = outdir, regexes = regexes,
          max_objects = max_files, dryrun = dryrun
        )
        if (!dryrun) {
          self$filesystem <- "local"
          self$path <- outdir
        }
      } else {
        d <- self$list_files_filter_relevant(regexes = regexes, max_files = max_files)
      }
      return(d)
    },
    #' @description Tidy given files.
    #' @param x Tibble with `localpath` to file and the function `type` to parse it.
    tidy_files = function(x) {
      # awesomeness
      tidy_files(x, envir = self)
    },
    #' @description Write tidy data.
    #' @param x Tibble with tidy `data` list-column.
    #' @param outdir Directory path to output tidy files.
    #' @param prefix Prefix of output files.
    #' @param format Format of output files.
    #' @param drid dracarys ID to use for the dataset (e.g. `wfrid.123`, `prid.456`).
    write = function(x, outdir = NULL, prefix = NULL, format = "tsv", drid = NULL) {
      assertthat::assert_that(!is.null(prefix))
      if (!is.null(outdir)) {
        prefix <- file.path(outdir, prefix)
      }
      d_write <- x |>
        dplyr::rowwise() |>
        dplyr::mutate(
          p = glue("{prefix}_{.data$name}"),
          out = list(write_dracarys(obj = .data$data, prefix = .data$p, out_format = format, drid = drid))
        ) |>
        dplyr::ungroup() |>
        dplyr::select("name", "data") |>
        tibble::deframe()
      invisible(d_write)
    }
  ) # end public
)
