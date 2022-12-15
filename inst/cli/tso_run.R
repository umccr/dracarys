require(here)
require(glue)
require(dplyr)
require(readr)
require(dracarys)

# SQL
# select * from data_portal.data_portal_gdsfile where regexp_like(path, 'SampleAnalysisResults.json.gz') order by time_created desc;

d <- here("nogit/tso/sql/81b58b53-5d3e-4ff3-bb37-e02e61d91838.csv") |>
  read_csv(col_names = TRUE)

x <- d |>
  mutate(
    sbj = sub("/analysis_data/(SBJ.*?)/.*", "\\1", path),
    dir = dirname(path),
    gds_indir = glue("gds://{volume_name}{dir}/")
  ) |>
  group_by(sbj) |>
  mutate(
    n_samp = n(),
    sbj2 = if_else(n_samp > 1, glue("{sbj}_{dplyr::row_number()}"), sbj)
  ) |>
  ungroup() |>
  select(sbj, sbj2, gds_indir, date = time_created) |>
  arrange(sbj2) |>
  mutate(outdir = here(glue("nogit/tso/2022-12-13/{sbj2}"))) |>
  select(sbj2, gds_indir, outdir)

token <- Sys.getenv("ICA_ACCESS_TOKEN_PROD")
for (i in 1:20) {
  dracarys::dracarys_tso(indir = x$gds_indir[i], outdir = x$outdir[i], dryrun = FALSE, token = token)
}
