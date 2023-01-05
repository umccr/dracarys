#' Generate TSO Report
#'
#' Generates a HTML report with results from the UMCCR TSO500 ctDNA workflow.
#' @param tidy_rds Path to RDS object output from `tso_tidy`.
#' @param out_dir Output directory.
#' @param prefix Prefix of output file(s).
#' @param quiet Suppress warning and other messages.
#'
#' @return Path to rendered HTML report.
#' @examples
#' \dontrun{
#' tidy_rds <- "~/tmp/tso/SBJ02858.rds"
#' tso_tidy(...) |> saveRDS(tidy_rds)
#' quiet <- TRUE
#' out_dir <- file.path("~/tmp/tso/out2")
#' prefix <- "SBJ02858"
#' tso_rmd(tidy_rds = tidy_rds, out_dir = out_dir, prefix = prefix, quiet = quiet)
#' }
#' @export
tso_rmd <- function(tidy_rds, out_dir, prefix, quiet = FALSE) {
  e <- emojifont::emoji
  assertthat::assert_that(
    rlang::is_logical(quiet), file.exists(tidy_rds), fs::path_ext(tidy_rds) == "rds"
  )
  tmp_dir <- tempdir()
  rmd_dir <- system.file("rmd/tso", package = "dracarys")
  fs::dir_copy(rmd_dir, tmp_dir, overwrite = TRUE)
  rmd_file <- file.path(tmp_dir, "tso.Rmd")
  fs::dir_create(out_dir)
  out_html <- glue("{prefix}_tso_dracarys_report.html")
  # suppress DT large size warning
  options(DT.warn.size = FALSE)
  rmarkdown::render(
    input = rmd_file, output_dir = out_dir, output_file = I(out_html),
    params = list(rds = tidy_rds), quiet = quiet
  )
  html_path <- file.path(out_dir, out_html)
  cli::cli_alert_success("{date_log()} {e('sparkles')} {.emph {prefix}}: TSO HTML reportr at: {.file {html_path}}")
}
