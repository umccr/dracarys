#' List Files in Local Directory
#'
#' Lists files in a local directory.
#'
#' @param localdir Path to local directory.
#'
#' @return A tibble with file basename, size, last modification timestamp
#' and full path.
#' @examples
#' localdir <- system.file("R", package = "dracarys")
#' x <- local_list_files_dir(localdir)
#' @testexamples
#' expect_equal(names(x), c("bname", "size", "lastmodified", "path"))
#' @export
local_list_files_dir <- function(localdir) {
  fs::dir_info(path = localdir, recurse = TRUE, type = "file") |>
    dplyr::mutate(
      bname = basename(.data$path),
      lastmodified = .data$modification_time
    ) |>
    dplyr::select("bname", "size", "lastmodified", "path")
}

#' List Relevant Files In Local Directory
#'
#' Lists relevant files in a local directory.
#'
#' @param path Path to local directory.
#' @param regexes Tibble with `regex` and `fun`ction name (see example).
#' @return A tibble with file type, basename, size, last modified timestamp, and
#' path.
#'
#' @examples
#' path <- system.file("extdata/tso", package = "dracarys")
#' regexes <- tibble::tibble(regex = "multiqc_data\\.json$", fun = "MultiqcFile")
#' x <- local_list_files_filter_relevant(path, regexes)
#' @testexamples
#' expect_equal(nrow(x), 1)
#' @export
local_list_files_filter_relevant <- function(path, regexes = DR_FILE_REGEX) {
  local_list_files_dir(localdir = path) |>
    dplyr::mutate(
      type = purrr::map_chr(.data$bname, \(x) match_regex(x, regexes = regexes))
    ) |>
    dplyr::filter(!is.na(.data$type)) |>
    dplyr::select("type", "bname", "size", "lastmodified", "path")
}

#' Print current timestamp for logging
#'
#' @return Current timestamp as character.
#' @export
date_log <- function() {
  as.character(glue('[{format(Sys.time(), "%Y-%m-%dT%H:%M:%S%Z")}]'))
}

#' Session Information Kable
#'
#' Session information kables for vignettes.
#'
#' @param pkgs Vector of R packages to display in the vignette. By default returns all.
#'
#' @return A list with two kables containing information about the platform and
#' the specified packages.
#' @export
session_info_kable <- function(pkgs = NULL) {
  si <- session_info_tbls(pkgs)
  si_pl <- si$si_pl
  si_pkg <- si$si_pkg
  list(
    si_pl = knitr::kable(si_pl, caption = "Platform information."),
    si_pkg = knitr::kable(si_pkg, caption = "Main packages used.")
  )
}

session_info_tbls <- function(pkgs = NULL) {
  si <- sessioninfo::session_info(include_base = TRUE)
  assertthat::assert_that(all(c("platform", "packages") %in% names(si)))
  si_pl <- unclass(si[["platform"]]) |>
    unlist() |>
    tibble::enframe(name = "name", value = "value")
  si_pkg <- unclass(si[["packages"]]) |>
    dplyr::as_tibble() |>
    dplyr::select(
      "package",
      version = "ondiskversion", datestamp = "date", "source"
    )
  if (!is.null(pkgs)) {
    si_pkg <- si_pkg |>
      dplyr::filter(.data$package %in% pkgs)
  }
  list(
    si_pl = si_pl,
    si_pkg = si_pkg
  )
}

#' Output Format is Valid
#'
#' Checks that the specified output format is valid.
#' @param x Output format.
#' @export
dr_output_format_valid <- function(x) {
  format_choices <- c("tsv", "parquet", "delta", "rds")
  assertthat::assert_that(
    is.null(x) | all(x %in% format_choices),
    msg = paste0(
      "Output format should be one or more of ",
      paste(format_choices, collapse = ", "), ", or _only_ NULL."
    )
  )
}

#' Write Local R Dataframe to Spark-backed Table
#'
#' First converts local R dataframe to Spark DataFrame using SparkR,
#' then _appends_ it to the specified table.
#'
#' @param rdf Local R dataframe
#' @param outpath Path to output table.
#' @param drid dracarys ID to use for the dataset (e.g. `wfrid.123`, `prid.456`)
#' @param ... Additional arguments for `SparkR::write.df`.
#' @export
#' @examples
#' \dontrun{
#' rdf <- mtcars
#' rdf2tab(rdf, "dev.wf1.mtcars", drid = "wfr.123")
#' }
rdf2tab <- function(rdf, outpath, drid = NULL, ...) {
  assertthat::assert_that(!is.null(drid))
  if (nrow(rdf) == 0) {
    # don't write empty dataframes
    return(NULL)
  }
  rdf <- rdf |>
    dplyr::mutate(dr_id = drid) |>
    dplyr::select("dr_id", dplyr::everything())
  sdf <- SparkR::createDataFrame(rdf)
  SparkR::write.df(sdf, path = outpath, mode = "append", ...)
}

write_dracarys <- function(obj, prefix, out_format, drid = NULL) {
  prefix <- as.character(prefix) # glue is error prone in Spark
  dr_output_format_valid(out_format)
  if ("delta" %in% out_format) {
    rdf2tab(rdf = obj, outpath = prefix, drid = drid)
    return(invisible(obj)) # skip other outputs
  }
  if ("tsv" %in% out_format) {
    fs::dir_create(dirname(prefix))
    tsv_out <- glue("{prefix}.tsv.gz")
    readr::write_tsv(obj, tsv_out)
  }
  if ("parquet" %in% out_format) {
    fs::dir_create(dirname(prefix))
    parquet_out <- glue("{prefix}.parquet")
    arrow::write_parquet(obj, parquet_out)
  }
  if ("rds" %in% out_format) {
    fs::dir_create(dirname(prefix))
    rds_out <- glue("{prefix}.rds")
    readr::write_rds(obj, rds_out)
  }
  # also gets returned in case of NULL out_format
  return(invisible(obj))
}

#' Write List of Tidy Tibbles
#'
#' @param list_of_tbls List of tidy tibbles.
#' @param out_dir Output directory.
#' @param prefix Prefix of output file(s).
#' @param out_format Format of output file(s).
#' @param drid dracarys ID to use for the dataset (e.g. `wfrid.123`, `prid.456`).
#'
#' @return Tibble with nested objects that have been written to the output directory.
#' @export
write_dracarys_list_of_tbls <- function(list_of_tbls, out_dir = NULL, prefix = NULL, out_format = "tsv", drid = NULL) {
  assertthat::assert_that(!is.null(prefix))
  if (!is.null(out_dir)) {
    prefix <- file.path(out_dir, prefix)
  }
  d_write <- list_of_tbls |>
    tibble::enframe(name = "section") |>
    dplyr::rowwise() |>
    dplyr::mutate(
      section_low = tolower(.data$section),
      p = glue("{prefix}_{.data$section_low}"),
      out = list(write_dracarys(obj = .data$value, prefix = .data$p, out_format = out_format, drid = drid))
    ) |>
    dplyr::ungroup() |>
    dplyr::select("section", "value") |>
    tibble::deframe()
  invisible(d_write)
}

#' Create Empty Tibble Given Column Names
#'
#' From https://stackoverflow.com/a/62535671/2169986. Useful for handling
#' edge cases with empty data.
#'
#' @param ctypes Character vector of column types corresponding to `cnames`.
#' @param cnames Character vector of column names to use.
#'
#' @return A tibble with 0 rows and the given column names.
#' @export
empty_tbl <- function(cnames, ctypes = readr::cols(.default = "c")) {
  readr::read_csv("\n", col_names = cnames, col_types = ctypes)
}

read_tsvgz <- function(x, ...) {
  if (is_url(x)) {
    res <- base::url(x) |>
      base::gzcon() |>
      readr::read_tsv(...)
    return(res)
  }
  readr::read_tsv(x, ...)
}

read_jsongz_jsonlite <- function(x, ...) {
  if (is_url(x)) {
    # https://github.com/jeroen/jsonlite/issues/414
    res <- base::url(x) |>
      base::gzcon() |>
      jsonlite::parse_json(...)
    return(res)
  }
  jsonlite::read_json(x, ...)
}

read_jsongz_rjsonio <- function(x, ...) {
  if (is_url(x)) {
    # https://github.com/umccr/dracarys/issues/74
    res <- base::url(x) |>
      base::gzcon() |>
      readr::read_lines() |>
      RJSONIO::fromJSON(...)
    return(res)
  }
  RJSONIO::fromJSON(x, ...)
}

#' Grep File Pattern
#'
#' @param path Path to look for file.
#' @param regexp A regular expression (e.g. [.]csv$) passed on to `grep()` to filter paths.
#'
#' @return The path to the file or an empty string if no match is found.
#' @export
grep_file <- function(path = ".", regexp) {
  x <- fs::dir_ls(path, recurse = TRUE, type = "file", regexp = regexp)
  if (length(x) > 1) {
    fnames <- paste(x, collapse = ", ")
    cli::cli_abort("More than 1 match found for {regexp} ({fnames}). Aborting.")
  }
  if (length(x) == 0) {
    return("") # file.exists("") returns FALSE
  }
  return(x)
}

#' @noRd
dummy1 <- function() {
  # Solves R CMD check: Namespaces in Imports field not imported from
  scales::pretty_breaks
  argparse::ArgumentParser
  here::here
}
