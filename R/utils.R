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

#' Read Google LIMS
#'
#' Reads UMCCR's Google LIMS spreadsheet.
#'
#' @return Tibble with all columns and rows from the Google LIMS spreadsheet.
#' @export
glims_read <- function() {
  lims_key <- googledrive::drive_find("^Google LIMS$", shared_drive = "LIMS")$id
  lims <- lims_key |>
    googlesheets4::read_sheet("Sheet1", na = c(".", "", "-"), col_types = "c")
  lims |> readr::type_convert(col_types = readr::cols(.default = "c", Timestamp = "T"))
}


#' @noRd
dummy1 <- function() {
  # Solves R CMD check: Namespaces in Imports field not imported from
  scales::pretty_breaks
  argparse::ArgumentParser
  here::here
}
