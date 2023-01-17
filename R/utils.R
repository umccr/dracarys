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
#' @param pkgs Vector of R packages to display in the vignette.
#'
#' @return A list with two kables containing information about the platform and
#' the specified packages (`pkgs`).
#' @export
session_info_kable <- function(pkgs) {
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
    ) |>
    dplyr::filter(.data$package %in% pkgs)

  list(
    si_pl = knitr::kable(si_pl, caption = "Platform information."),
    si_pkg = knitr::kable(si_pkg, caption = "Main packages used.")
  )
}

output_format_valid <- function(x) {
  format_choices <- c("tsv", "parquet", "both")
  assertthat::assert_that(
    length(x) == 1,
    x %in% format_choices
  )
}

write_dracarys <- function(obj, prefix, out_format) {
  output_format_valid(out_format)
  fs::dir_create(dirname(prefix))
  if (out_format %in% c("tsv", "both")) {
    tsv_out <- glue("{prefix}.tsv.gz")
    readr::write_tsv(obj, tsv_out)
  }
  if (out_format %in% c("parquet", "both")) {
    parquet_out <- glue("{prefix}.parquet")
    arrow::write_parquet(obj, parquet_out)
  }
  obj
}
