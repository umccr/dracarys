#' Generate TSO Report
#'
#' Generates a HTML report with results from the UMCCR TSO500 ctDNA workflow.
#'
#' @inheritParams tso_tidy
#' @param quiet Suppress warning and other messages.
#' @param token ICA access token (by default uses $ICA_ACCESS_TOKEN env var).
#'
#' @return Path to rendered HTML report.
#' @examples
#' \dontrun{
#' in_dir <- here::here(glue("nogit/tso/2022-12-13/SBJ02858/dracarys_gds_sync"))
#' gds_local_dir <- NULL
#' out_format <- "tsv"
#' quiet <- TRUE
#' out_dir <- file.path(in_dir, "../out")
#' prefix <- "SBJ02858"
#' tso_rmd(in_dir = in_dir, out_dir = out_dir, prefix = prefix, quiet = quiet)
#' }
#' @export
tso_rmd <- function(in_dir, out_dir, prefix, gds_local_dir = NULL, out_format = "tsv",
                    token = Sys.getenv("ICA_ACCESS_TOKEN"), quiet = FALSE) {
  assertthat::assert_that(quiet %in% c(FALSE, TRUE))
  tmp_dir <- tempdir()
  rmd_dir <- system.file("rmd/tso", package = "dracarys")
  fs::dir_copy(rmd_dir, tmp_dir, overwrite = TRUE)
  rmd_file <- file.path(tmp_dir, "tso.Rmd")
  out_dir_rds <- file.path(out_dir, "rds")
  out_dir_html <- file.path(out_dir, "html")
  out_dir_dr <- file.path(out_dir, "dracarys")
  c(out_dir_dr, out_dir_rds, out_dir_html) |> fs::dir_create()
  dat <- tso_tidy(
    in_dir = in_dir, out_dir = out_dir_dr, prefix = prefix, gds_local_dir = gds_local_dir,
    out_format = out_format, dryrun = FALSE, token = token
  )
  dat_rds <- file.path(out_dir_rds, glue("{prefix}_tso_dracarys_data.rds"))
  readr::write_rds(dat, dat_rds)
  out_html <- glue("{prefix}_tso_dracarys_report.html")
  # suppress DT large size warning
  options(DT.warn.size = FALSE)
  pars <- list(rds = dat_rds)
  rmarkdown::render(
    input = rmd_file, output_dir = out_dir_html, output_file = I(out_html),
    params = pars, quiet = quiet
  )
}
