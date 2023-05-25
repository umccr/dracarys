# portal workflow meta subset
require(dracarys)
require(here)
require(dplyr)
require(readr)

pmeta <- here::here("nogit/data_portal/2023-05-21_workflows.csv")
wfs <- c(
  "bcl_convert", "rnasum", "tso_ctdna_tumor_only",
  "umccrise", "wgs_alignment_qc", "wgs_tumor_normal", "wts_tumor_only"
)
dracarys:::portal_meta_read(pmeta) |>
  dplyr::filter(.data$type_name %in% wfs, .data$end_status == "Succeeded") |>
  dplyr::group_by(.data$type_name) |>
  dplyr::slice_head(n = 4) |>
  readr::write_csv(here::here("inst/extdata/portal_meta_top4.csv"))
