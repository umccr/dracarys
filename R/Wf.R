#' @title Workflow
#'
#' @description Workflow is a base R6 class representing a bioinformatic
#' workflow run from a UMCCR workflow manager.
#'
#' A workflow has:
#'
#' - an output directory path with all the result output files (either on GDS, S3 or
#' local filesystem)
#' - a subset of files that are of interest for ingestion
#'   - tibble with full path and basename columns
#' - a set of parsers that can parse and tidy those files
#'   - each parser takes a path and returns a tidy tibble
#' - a list of tidy tibbles (or a tibble with nested tibbles)
#'
#' @examples
#' \dontrun{
#' p1 <- "~/icav1/g/production/analysis_data"
#' p <- file.path(p1, "SBJ01155/umccrise/202408300c218043/L2101566__L2101565")
#' um <- Wf$new(p, "umccrise")
#' }
#'
#' @export
Wf <- R6::R6Class(
  "Wf",
  public = list(
    #' @field path (`character(1)`)\cr
    #' Path to directory with raw workflow results (from GDS, S3, or local filesystem).
    #' @field wname (`character(1)`)\cr
    #' Name of workflow (e.g. umccrise, sash).
    #' @field filesystem (`character(1)`)\cr
    #' Filesystem of `path`.
    path = NULL,
    wname = NULL,
    filesystem = NULL,
    #' @description Create a new Workflow object.
    #' @param path Output directory path with results.
    #' @param wname Name of workflow.
    initialize = function(path = NULL, wname = NULL) {
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
    #' @param max_files Maximum number of files to list.
    #' @param ica_token ICA access token (def: $ICA_ACCESS_TOKEN env var).
    #' @param ... Passed on to `gds_list_files_dir` function.
    list_files = function(max_files = 1000, ica_token = Sys.getenv("ICA_ACCESS_TOKEN"), ...) {
      path <- self$path
      if (self$filesystem == "gds") {
        d <- gds_list_files_dir(
          gdsdir = path, token = ica_token, page_size = max_files, ...
        )
      } else if (self$filesystem == "s3") {
        d <- s3_list_files_dir(s3dir = path, max_objects = max_files)
      } else {
        d <- local_list_files_dir(localdir = path)
      }
      return(d)
    },
    #' @description List dracarys files under given path
    #' @param regexes Tibble with `regex` and `fun`ction name.
    #' @param max_files Maximum number of files to list.
    #' @param ica_token ICA access token (def: $ICA_ACCESS_TOKEN env var).
    #' @param ... Passed on to the `gds_list_files_filter_relevant` or
    #' the `s3_list_files_filter_relevant` function.
    list_files_filter_relevant = function(regexes = NULL,
                                          max_files = 1000,
                                          ica_token = Sys.getenv("ICA_ACCESS_TOKEN"), ...) {
      assertthat::assert_that(!is.null(regexes))
      path <- self$path
      if (self$filesystem == "gds") {
        d <- gds_list_files_filter_relevant(
          gdsdir = path, regexes = regexes, token = ica_token, page_size = max_files, ...
        )
      } else if (self$filesystem == "s3") {
        d <- s3_list_files_filter_relevant(
          s3dir = path, regexes = regexes, max_objects = max_files, ...
        )
      } else {
        d <- local_list_files_filter_relevant(localdir = path, regexes = regexes)
      }
      d
    },
    #' @description Download files from GDS/S3 to local filesystem.
    #' @param outdir Path to output directory.
    #' @param regexes Tibble with `regex` and `fun`ction name.
    #' @param ica_token ICA access token (def: $ICA_ACCESS_TOKEN env var).
    #' @param max_files Maximum number of files to list.
    #' @param dryrun If TRUE, just list the files that will be downloaded (don't
    #' download them).
    #' @param recursive Should files be returned recursively _in and under_ the specified
    #' GDS directory, or _only directly in_ the specified GDS directory (def: TRUE via ICA API).
    download_files = function(outdir, regexes = NULL,
                              ica_token = Sys.getenv("ICA_ACCESS_TOKEN"),
                              max_files = 1000, dryrun = FALSE, recursive = NULL) {
      # TODO: add envvar checker
      path <- self$path
      assertthat::assert_that(!is.null(regexes))
      if (self$filesystem == "gds") {
        d <- dr_gds_download(
          gdsdir = path, outdir = outdir, regexes = regexes, token = ica_token,
          page_size = max_files, dryrun = dryrun, recursive = recursive
        )
        self$filesystem <- "local"
        self$path <- outdir
      } else if (self$filesystem == "s3") {
        d <- dr_s3_download(
          s3dir = path, outdir = outdir, regexes = regexes,
          max_objects = max_files, dryrun = dryrun
        )
        self$filesystem <- "local"
        self$path <- outdir
      } else {
        d <- self$list_files_filter_relevant(regexes = regexes)
      }
      return(d)
    }
  ) # end public
)
