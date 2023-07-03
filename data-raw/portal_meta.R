# portal workflow meta subset
require(dracarys)
require(here)
require(dplyr)
require(readr)

pmeta <- here::here("nogit/data_portal/workflows/2023-07-02.csv")
wfs <- c(
  "bcl_convert", "rnasum", "tso_ctdna_tumor_only",
  "umccrise", "wgs_alignment_qc", "wgs_tumor_normal", "wts_tumor_only"
)
readr::read_csv(pmeta, col_types = readr::cols(.default = "c")) |>
  dplyr::filter(.data$type_name %in% wfs, .data$end_status == "Succeeded") |>
  dplyr::group_by(.data$type_name) |>
  dplyr::slice_head(n = 4) |>
  dplyr::ungroup() |>
  readr::write_csv(here::here("inst/extdata/portal_meta_top4.csv"))
