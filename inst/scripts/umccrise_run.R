require(dracarys)
require(here)
require(dplyr)
require(readr)

# SQL
#   select * from data_portal.data_portal_gdsfile where regexp_like(path, 'cancer_report_tables/.*qc_summary.tsv.gz') order by time_created desc;

d <- here("nogit/umccrise/sql/7c348d27-adcb-427f-a8e8-549bcc0c8490_2023-01-18.csv") |>
  read_csv(col_names = TRUE)


x <- d |>
  mutate(
    sbj = sub("/analysis_data/(SBJ.*?)/.*", "\\1", path),
    dir = dirname(path),
    gds_indir = glue("gds://{volume_name}{dir}/")
  ) |>
  group_by(sbj) |>
  # TODO: check if results already exist for same sbj
  mutate(
    n_samp = n(),
    sbj2 = if_else(n_samp > 1, glue("{sbj}_{dplyr::row_number()}"), sbj)
  ) |>
  ungroup() |>
  select(sbj, sbj2, gds_indir, date = time_created) |>
  arrange(sbj2) |>
  mutate(
    outdir = here(glue("nogit/umccrise/2023-01-18/{sbj2}")),
    local_indir = file.path(outdir, "dracarys_gds_sync")
  ) |>
  select(sbj2, gds_indir, outdir, local_indir, date)


token <- Sys.getenv("ICA_ACCESS_TOKEN_PROD")
dryrun <- F
for (i in 101:276) {
  print(i)
  # print(x$gds_indir[i])
  # print(x$local_indir[i])
  umccr_tidy(in_dir = x$gds_indir[i], out_dir = x$outdir[i], prefix = x$sbj2[i], dryrun = dryrun, token = token)
  # umccr_tidy(in_dir = x$local_indir[i], out_dir = x$outdir[i], prefix = x$sbj2[i], dryrun = FALSE, token = token)
}
