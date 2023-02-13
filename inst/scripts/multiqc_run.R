require(dracarys)
require(here)
require(dplyr)
require(readr)

# SQL
# select * from data_portal.data_portal_gdsfile where regexp_like(path, 'multiqc_data.json') order by time_created desc;

d <- here("nogit/multiqc/sql/e7fd7995-d280-4457-a1ac-0cf46bf723dd_gds_multiqcjson_query_2023-02-03.csv") |>
  read_csv(col_names = TRUE)

date1 <- "2023-02-07"

x <- d |>
  filter(!grepl("bclconvert|interop", path)) |>
  mutate(
    sbj = sub("/analysis_data/(SBJ.*?)/.*", "\\1", path),
    workflow = sub("/analysis_data/SBJ.*/(.*?)/.*", "\\1", path),
    dir = dirname(path),
    gds_indir = glue("gds://{volume_name}{dir}/"),
    unique_hash = substr(unique_hash, 1, 6),
    time_created = as.Date(time_created)
  ) |>
  select(sbj, workflow, gds_indir, time_created, unique_hash) |>
  filter(workflow == "umccrise") |>
  mutate(
    outdir = here(glue("nogit/multiqc/{date1}/{sbj}/{time_created}_{unique_hash}")),
    local_indir = file.path(outdir, "dracarys_gds_sync")
  ) |>
  arrange(sbj, time_created) |>
  select(sbj, gds_indir, outdir, local_indir)


token <- Sys.getenv("ICA_ACCESS_TOKEN_PROD")
dryrun <- TRUE
dryrun <- FALSE
for (i in 201:288) {
  print(i)
  print(x$gds_indir[i])
  # print(x$local_indir[i])
  dracarys::umccr_tidy(in_dir = x$gds_indir[i], out_dir = x$outdir[i], prefix = x$sbj[i], dryrun = dryrun, token = token)
  # dracarys::tso_tidy(in_dir = x$local_indir[i], out_dir = x$outdir[i], prefix = x$sbj2[i], dryrun = FALSE, token = token)
}

x |>
  mutate(
    cmd = glue("./dracarys.R tso -i {gds_indir} -o {outdir} -r {outdir}/report_dir -p {sbj2} --rds_dir {outdir}/rds_dir --quiet_rmd")
  ) |>
  select(cmd) |>
  write_tsv(here("inst/cli/run.sh"), col_names = FALSE)
