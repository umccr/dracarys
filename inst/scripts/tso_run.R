require(tidyverse)
require(jsonlite)
require(here)
require(glue)

pmeta <- here::here("nogit/data_portal/2023-05-07_workflows_5f43da12-0b6a-41ce-86c9-9ee62df8792e.csv") |>
  dracarys::cttso_metadata()

x <- pmeta |>
  dplyr::mutate(
    year1 = lubridate::year(start),
    month1 = sprintf("%02d", lubridate::month(start)),
    outdir = here::here(glue("nogit/warehouse/cttso/{year1}/{month1}/{wfr_id}")),
    local_indir = file.path(outdir, "dracarys_gds_sync")
  ) |>
  dplyr::select(
    -c(year1, month1)
  )

token <- Sys.getenv("ICA_ACCESS_TOKEN_PRO")
dryrun <- TRUE
dryrun <- FALSE

for (i in seq_len(nrow(x))) {
  print(i)
  print(x$gds_indir[i])
  dracarys::umccr_tidy(in_dir = x$gds_indir[i], out_dir = x$outdir[i], prefix = x$LibraryID_w_rerun[i], dryrun = dryrun, token = token)
  # print(x$local_indir[i])
  # dracarys::umccr_tidy(in_dir = x$local_indir[i], out_dir = x$outdir[i], prefix = x$sbj[i], dryrun = dryrun, token = token)
}
