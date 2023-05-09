require(tidyverse)
require(jsonlite)
require(here)
require(glue)

pmeta <- "nogit/data_portal/2023-05-07_workflows_5f43da12-0b6a-41ce-86c9-9ee62df8792e.csv"

#---- cttso (succeeded)----#

lims <- read_tsv("~/Downloads/Google LIMS - Sheet1.tsv")
table(wf1$LibraryID_w_rerun %in% lims$LibraryID)
lims |> filter(LibraryID %in% wf1$LibraryID_w_rerun)

x <- wf1 |>
  mutate(
    year1 = lubridate::year(start),
    month1 = sprintf("%02d", lubridate::month(start)),
    outdir = here(glue("nogit/warehouse/cttso/{year1}/{month1}/{wfr_id}")),
    local_indir = file.path(outdir, "dracarys_gds_sync")
  ) |>
  select(
    -c(year1, month1)
  )

token <- Sys.getenv("ICA_ACCESS_TOKEN_PRO")
dryrun <- TRUE
dryrun <- FALSE

# for (i in seq_len(nrow(x))) {
for (i in 11:20) {
  print(i)
  print(x$gds_indir[i])
  dracarys::umccr_tidy(in_dir = x$gds_indir[i], out_dir = x$outdir[i], prefix = x$LibraryID_w_rerun[i], dryrun = dryrun, token = token)
  # print(x$local_indir[i])
  # dracarys::umccr_tidy(in_dir = x$local_indir[i], out_dir = x$outdir[i], prefix = x$sbj[i], dryrun = dryrun, token = token)
}
