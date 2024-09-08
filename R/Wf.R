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
    #' Output directory path with results.
    #' @field type (`character(1)`)\cr
    #' Type of workflow (e.g. umccrise, sash).
    #' @field filesystem (`character(1)`)\cr
    #' Filesystem of `path`.
    path = NULL,
    type = NULL,
    filesystem = NULL,
    #' @description Create a new Workflow object.
    #' @param path Output directory path with results.
    #' @param type Type of workflow.
    initialize = function(path = NULL, type = NULL) {
      types <- c(
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
      assertthat::assert_that(type %in% types)
      self$path <- path
      self$type <- type
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
        "type", self$type,
        "filesystem", self$filesystem
      )
      print(res)
      invisible(self)
    },
    #' @description List all files under given path.
    #' @param max_objects Maximum number of objects to list.
    list_files = function(max_objects = 1000, ica_token = Sys.getenv("ICA_ACCESS_TOKEN")) {
      path <- self$path
      if (self$filesystem == "gds") {
        d <- gds_list_files_dir(gdsdir = path, page_size = max_objects, token = ica_token)
      } else if (self$filesystem == "s3") {
        d <- s3_list_files_dir(s3dir = path, max_objects = max_objects)
      } else {
        d <- local_list_files_dir(localdir = path)
      }
      return(d)
    },
    #' @description List dracarys files under given path
    #' @param page_size Page size
    #' @param regexes Tibble with `regex` and `fun`ction name.
    list_files_filter_relevant = function(regexes = NULL, max_objects = 1000, ica_token = Sys.getenv("ICA_ACCESS_TOKEN"), ...) {
      assertthat::assert_that(!is.null(regexes))
      path <- self$path
      if (self$filesystem == "gds") {
        d <- gds_list_files_filter_relevant(gdsdir = path, token = ica_token, page_size = max_objects, regexes = regexes, ...)
      } else if (self$filesystem == "s3") {
        # "type", "bname", "size", "date_utc", "path"
        d <- s3_list_files_filter_relevant(s3dir = path, max_objects = max_objects, regexes = regexes, ...)
      } else {
        # "type", "bname", "path"
        d <- local_list_files_filter_relevant(path = path, regexes = regexes)
      }
      d
    },
    download_files = function(ica_token = Sys.getenv("ICA_ACCESS_TOKEN")) {
      # TODO: add envvar checker
      path <- self$path
      if (self$filesystem == "gds") {
        d <- dr_gds_download(gdsdir = path, token = ica_token)
      } else if (self$filesystem == "s3") {
        d <- s3_download_files(s3dir = path)
      } else {
        d <- local_download_files(localdir = path)
      }
    }
  ) # end public
)
