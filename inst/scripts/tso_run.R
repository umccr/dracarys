require(here)
require(glue)
require(dplyr)
require(readr)
require(dracarys)

# SQL
# select * from data_portal.data_portal_gdsfile where regexp_like(path, 'SampleAnalysisResults.json.gz') order by time_created desc;

# d <- here("nogit/tso/sql/81b58b53-5d3e-4ff3-bb37-e02e61d91838.csv") |>
#   read_csv(col_names = TRUE)
# d1 <- here("nogit/tso/sql/7991c298-b26e-481c-b8f5-df530303b3d5.csv") |>
#   read_csv(col_names = TRUE)
# just keep the new tso samples
# d <- d1 |>
#   anti_join(d)
d <- here("nogit/tso/sql/7991c298-b26e-481c-b8f5-df530303b3d5.csv") |>
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
    outdir = here(glue("nogit/tso/2022-12-13/{sbj2}")),
    local_indir = file.path(outdir, "dracarys_gds_sync")
  ) |>
  select(sbj2, gds_indir, outdir, local_indir)


token <- Sys.getenv("ICA_ACCESS_TOKEN_PROD")
for (i in 319:320) {
  print(i)
  # print(x$gds_indir[i])
  print(x$local_indir[i])
  dracarys::dracarys_tso(indir = x$gds_indir[i], outdir = x$outdir[i], dryrun = FALSE, token = token)
  # dracarys::dracarys_tso(indir = x$local_indir[i], outdir = x$outdir[i], dryrun = FALSE, token = token)
  # dracarys::tso_rmd(prefix = x$pref[i], out_file = x$out_html)
}
