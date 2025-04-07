#' @title Workflow
#'
#' @description Workflow is a base R6 class representing a bioinformatic
#' workflow run from a UMCCR workflow manager.
#'
#' A workflow has:
#'
#' - a directory path with all the raw output files (either on S3 or
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
  private = list(
    .path = NULL,
    .wname = NULL,
    .regexes = NULL,
    .filesystem = NULL
  ),
  active = list(
    #' @field regexes Get/Set regexes. Tibble with file `regex` and `fun`ction
    #' to parse it.
    regexes = function(value) {
      if (missing(value)) {
        private$.regexes
      } else {
        assertthat::assert_that(
          tibble::is_tibble(value),
          all(c("regex", "fun") %in% colnames(value))
        )
        private$.regexes <- value
      }
    }
  ),
  public = list(
    #' @description Create a new Workflow object.
    #' @param path Path to directory with raw workflow results.
    #' @param wname Name of workflow.
    #' @param regexes Tibble with file `regex` and `fun`ction to parse it.
    initialize = function(path = NULL, wname = NULL, regexes = NULL) {
      wnames <- c(
        "bcl_convert",
        "tso_ctdna_tumor_only",
        "cttsov2",
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
      subwnames <- c("dragen")
      assertthat::assert_that(wname %in% c(wnames, subwnames))
      path <- sub("/$", "", path) # remove potential trailing slash
      private$.path <- path
      private$.wname <- wname
      private$.filesystem <- dplyr::case_when(
        grepl("^s3://", path) ~ "s3",
        .default = "local"
      )
      assertthat::assert_that(
        tibble::is_tibble(regexes),
        all(c("regex", "fun") %in% colnames(regexes))
      )
      private$.regexes <- regexes
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
        "nregexes", as.character(nrow(private$.regexes))
      )
      print(res)
      invisible(self)
    },
    #' @description List all files under given path.
    #' @param path Path with raw results.
    #' @param max_files Max number of files to list.
    list_files = function(path = private$.path, max_files = 1000) {
      if (private$.filesystem == "s3") {
        d <- s3_list_files_dir(s3dir = path, max_objects = max_files)
      } else {
        d <- local_list_files_dir(localdir = path, max_files = max_files)
      }
      return(d)
    },
    #' @description List dracarys files under given path
    #' @param path Path with raw results.
    #' @param max_files Max number of files to list.
    #' @param ... Passed on to `s3_list_files_filter_relevant`.
    list_files_filter_relevant = function(
      path = private$.path,
      max_files = 1000,
      ...
    ) {
      regexes <- private$.regexes
      assertthat::assert_that(!is.null(regexes))
      if (private$.filesystem == "s3") {
        d <- s3_list_files_filter_relevant(
          s3dir = path,
          regexes = regexes,
          max_objects = max_files,
          ...
        )
      } else {
        d <- local_list_files_filter_relevant(
          localdir = path,
          regexes = regexes,
          max_files = max_files
        )
      }
      d
    },
    #' @description For DOWNLOAD_ONLY files, just return the input path.
    #' @param x Path with raw results.
    #' @param suffix Suffix.
    DOWNLOAD_ONLY = function(x, suffix = "") {
      tibble::tibble(
        name = glue("DOWNLOAD_ONLY{suffix}"),
        data = list(tibble::tibble(input_path = x))
      )
    },
    #' @description Download files from S3 to local filesystem.
    #' @param path Path with raw results.
    #' @param outdir Path to output directory.
    #' @param max_files Max number of files to list.
    #' @param dryrun If TRUE, just list the files that will be downloaded (don't
    #' download them).
    download_files = function(
      path = private$.path,
      outdir,
      max_files = 1000,
      dryrun = FALSE
    ) {
      regexes <- private$.regexes
      assertthat::assert_that(!is.null(regexes))
      if (private$.filesystem == "s3") {
        d <- dr_s3_download(
          s3dir = path,
          outdir = outdir,
          regexes = regexes,
          max_objects = max_files,
          dryrun = dryrun
        )
        if (!dryrun) {
          private$.filesystem <- "local"
          private$.path <- outdir
        }
      } else {
        d <- self$list_files_filter_relevant(
          regexes = regexes,
          max_files = max_files
        )
      }
      return(d)
    },
    #' @description Tidy given files.
    #' @param x Tibble with `localpath` to file and the function `type` to parse it.
    tidy_files = function(x) {
      # awesomeness
      tidy_files(x, envir = self) |>
        dplyr::arrange(.data$name)
    },
    #' @description Write tidy data.
    #' @param x Tibble with tidy `data` list-column.
    #' @param outdir Directory path to output tidy files.
    #' @param prefix Prefix of output files.
    #' @param format Format of output files.
    #' @param drid dracarys ID to use for the dataset (e.g. `wfrid.123`, `prid.456`).
    #' @param dbconn Database connection object.
    write = function(
      x,
      outdir = NULL,
      prefix = NULL,
      format = "tsv",
      drid = NULL,
      dbconn = NULL
    ) {
      assertthat::assert_that(!is.null(prefix))
      assertthat::assert_that(all(c("name", "data") %in% colnames(x)))
      if (!is.null(outdir)) {
        prefix <- file.path(outdir, prefix)
      }
      d_write <- x |>
        dplyr::rowwise() |>
        # for db we want the tibble name
        dplyr::mutate(
          p = ifelse(
            grepl("DOWNLOAD_ONLY", .data$name),
            as.character(.data$data |> dplyr::pull("input_path")),
            ifelse(
              format == "db",
              as.character(.data$name),
              as.character(glue("{prefix}_{.data$name}"))
            )
          ),
          out = ifelse(
            !grepl("DOWNLOAD_ONLY", .data$name),
            list(write_dracarys(
              obj = .data$data,
              prefix = .data$p,
              out_format = format,
              drid = drid,
              dbconn = dbconn
            )),
            list(.data$data)
          )
        ) |>
        dplyr::ungroup() |>
        dplyr::select("name", "data", prefix = "p")
      invisible(d_write)
    }
  ) # end public
)
