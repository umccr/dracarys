#' Generate TSO Report
#'
#' @param prefix Prefix for input files i.e. '/path/to/file1'.ext, where 'file1'
#' is the common prefix of all files.
#' @param out_file Name of output HTML file (needs '.html' suffix) (def: `{tumor_name}_cancer_report.html`).
#' @param quiet Suppress warning and other messages.
#'
#' @return Path to rendered HTML report.
#' @export
tso_rmd <- function(prefix = NULL, out_file = NULL, quiet = FALSE) {
  assertthat::assert_that(
    quiet %in% c(FALSE, TRUE),
    !is.null(prefix)
  )
  if (!is.null(out_file)) {
    assertthat::assert_that(
      is.character(out_file),
      tools::file_ext(out_file) == "html"
    )
  } else {
    out_file <- glue::glue("tso_report.html")
  }
  tmp_dir <- tempdir()
  # R's file.copy('foo', 'bar/baz') copies 'foo' to 'bar/baz/foo'
  rmd_dir <- system.file("rmd/tso", package = "dracarys")
  rmd_dir_b <- basename(rmd_dir)
  fs::dir_copy(rmd_dir, tmp_dir) # /path/to/rmd/tso -> /tmp/tso
  rmd_file <- file.path(tmp_dir, "tso", "tso.Rmd")
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
