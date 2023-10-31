# portal workflow meta subset
require(dracarys)
require(here)
require(dplyr)
require(purrr)
require(readr)


wfs <- c(
  "bcl_convert", "rnasum", "tso_ctdna_tumor_only",
  "umccrise", "wgs_alignment_qc", "wgs_tumor_normal", "wts_tumor_only",
  "wts_alignment_qc",
  "oncoanalyser_wgs", "oncoanalyser_wgts_existing_both",
  "oncoanalyser_wts", "sash", "star_alignment"
)

account <- "stg"
get_top_succeeded <- function(wf, num_row = 10, num_top = 4) {
  dracarys::portal_meta_read(params = glue::glue("&type_name={wf}"), account = account, rows = num_row) |>
    dplyr::filter(.data$end_status == "Succeeded") |>
    dplyr::slice_head(n = num_top)
}

# get top 10 rows, then get top 4 successful runs
d <- wfs |>
  purrr::map(\(x) get_top_succeeded(x, 10, 4)) |>
  dplyr::bind_rows()
# leave dates as character
d |>
  readr::write_csv(here::here("inst/extdata/portal_meta_top4.csv"))
# date_fmt <- "%Y-%m-%dT%H:%M:%S"
