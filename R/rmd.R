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
#' indir <- here::here(glue("nogit/tso/2022-12-13/SBJ02858/dracarys_gds_sync"))
#' gds_local_dir <- NULL
#' out_format <- "tsv"
#' quiet <- TRUE
#' outprefix <- file.path(indir, "out/SBJ02858")
#' tso_rmd(indir = indir, outprefix = outprefix, quiet = quiet)
#' }
#' @export
tso_rmd <- function(indir, outprefix, gds_local_dir = NULL, out_format = "tsv",
                    token = Sys.getenv("ICA_ACCESS_TOKEN"), quiet = FALSE) {
  assertthat::assert_that(quiet %in% c(FALSE, TRUE))
  tmp_dir <- tempdir()
  rmd_dir <- system.file("rmd/tso", package = "dracarys")
  fs::dir_copy(rmd_dir, tmp_dir, overwrite = TRUE)
  rmd_file <- file.path(tmp_dir, "tso.Rmd")
  out_dir <- dirname(outprefix) # /path/to/prefix/sampleA -> /path/to/prefix
  out_bname <- basename(outprefix) # /path/to/prefix/sampleA -> sampleA
  # write the RDS and HTML at same level as out_dir
  out_dir_rds <- file.path(dirname(out_dir), "rds")
  out_dir_html <- file.path(dirname(out_dir), "html")
  c(out_dir, out_dir_rds, out_dir_html) |> fs::dir_create()
  dat <- tso_tidy(
    indir = indir, outprefix = outprefix, gds_local_dir = gds_local_dir,
    out_format = out_format, dryrun = FALSE, token = token
  )
  dat_rds <- file.path(out_dir_rds, glue("{out_bname}_dracarys.rds"))
  readr::write_rds(dat, dat_rds)
  out_html <- glue("{out_bname}_dracarys.html")
  # suppress DT large size warning
  options(DT.warn.size = FALSE)
  pars <- list(rds = dat_rds)
  rmarkdown::render(
    input = rmd_file, output_dir = out_dir_html, output_file = I(out_html),
    params = pars, quiet = quiet
  )
}
