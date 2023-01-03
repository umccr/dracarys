#' Generate TSO Report
#'
#' Generates a HTML report with results from the UMCCR TSO500 ctDNA workflow.
#'
#' @param prefix Prefix for input files generated from UMCCR TSO500 ctDNA workflow (i.e. '/path/to/sbjx'.ext, where 'sbjx'
#' is the common prefix of all relevant files in directory '/path/to/').
#' @param out_file Name of output HTML file (needs '.html' suffix) (def: `tso_report.html`).
#' @param quiet Suppress warning and other messages.
#'
#' @return Path to rendered HTML report.
#' @export
tso_rmd <- function(prefix = NULL, out_file = "tso_report.html", quiet = FALSE) {
  assertthat::assert_that(
    quiet %in% c(FALSE, TRUE),
    !is.null(prefix),
    is.character(out_file),
    tools::file_ext(out_file) == "html"
  )
  tmp_dir <- tempdir()
  rmd_dir <- system.file("rmd/tso", package = "dracarys")
  rmd_dir_b <- basename(rmd_dir)
  fs::dir_copy(rmd_dir, tmp_dir, overwrite = TRUE)
  rmd_file <- file.path(tmp_dir, "tso.Rmd")
  out_dir <- dirname(out_file)
  fs::dir_create(out_dir)
  # suppress DT large size warning
  options(DT.warn.size = FALSE)
  pars <- list(prefix = prefix)
  rmarkdown::render(
    input = rmd_file, output_dir = out_dir, output_file = I(out_file),
    params = pars, quiet = quiet
  )
}
