#' List Files in Local Directory
#'
#' Lists files in a local directory.
#'
#' @param localdir Path to local directory.
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
#' @inheritParams local_list_files_dir
#' @param regexes Tibble with `regex` and `fun`ction name (see example).
#' @return A tibble with file type, basename, size, last modified timestamp, and
#' path.
#'
#' @examples
#' localdir <- system.file("extdata/tso", package = "dracarys")
#' regexes <- tibble::tibble(regex = "multiqc_data\\.json$", fun = "MultiqcFile")
#' x <- local_list_files_filter_relevant(localdir, regexes)
#' @testexamples
#' expect_equal(nrow(x), 1)
#' @export
local_list_files_filter_relevant <- function(localdir, regexes = DR_FILE_REGEX) {
  local_list_files_dir(localdir = localdir) |>
    dplyr::mutate(
      type = purrr::map_chr(.data$bname, \(x) match_regex(x, regexes = regexes))
    ) |>
    dplyr::filter(!is.na(.data$type)) |>
    dplyr::select("type", "bname", "size", "lastmodified", localpath = "path")
}
