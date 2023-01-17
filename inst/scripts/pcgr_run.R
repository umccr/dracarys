require(dracarys)
require(here)
require(dplyr)
require(readr)

# SQL
#  select * from data_portal.data_portal_gdsfile where (regexp_like(path, '-somatic.pcgr.json.gz') AND NOT regexp_like(path, 'pcgr_run')) order by time_created desc;

d <- here("nogit/pcgr/sql/279893c1-79e4-44d7-808e-e34635bd9e50_2023-01-16.csv") |>
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
    outdir = here(glue("nogit/pcgr/2023-01-16/{sbj2}")),
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

# x |>
#   mutate(
#     cmd = glue("./dracarys.R tso -i {gds_indir} -o {outdir} -r {outdir}/report_dir -p {sbj2} --rds_dir {outdir}/rds_dir --quiet_rmd")
#   ) |>
#   select(cmd) |>
#   write_tsv(here("inst/cli/run.sh"), col_names = FALSE)

res <- read_rds(here("nogit/pcgr/rds/res_2023-01-17.rds"))

dplyr::left_join(x, res, by = c("sbj2" = "sbj")) |>
  dplyr::select(sbj = sbj2, date, fracIndels, predicted_class, tmb_estimate, n_tmb, gds_indir) |>
  readr::write_tsv(here("nogit/pcgr/res_2023-01-17.tsv"))
