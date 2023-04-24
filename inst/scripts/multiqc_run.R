require(dracarys)
require(here)
require(glue)
require(dplyr)
require(readr)

# SQL
# select * from data_portal.data_portal_gdsfile where regexp_like(path, 'multiqc_data.json') order by time_created desc;
d <- glue("nogit/multiqc/sql/2c527c86-7dee-4377-b18c-5ef01fca6375_gds_multiqcjson_query_2023-04-24.csv") |>
  here() |>
  read_csv(col_names = TRUE)

wf <- c("umccrise", "wgs_alignment_qc", "wgs_tumor_normal")

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
  filter(time_created >= "2023-04-23") |>
  filter(workflow %in% wf) |>
  mutate(
    outdir = here(glue("nogit/warehouse/{workflow}/{sbj}/{time_created}_{unique_hash}")),
    local_indir = file.path(outdir, "dracarys_gds_sync")
  ) |>
  arrange(sbj, time_created) |>
  select(sbj, gds_indir, outdir, local_indir, time_created)


token <- Sys.getenv("ICA_ACCESS_TOKEN_PRO")
dryrun <- TRUE
dryrun <- FALSE

for (i in seq_len(nrow(x))) {
  print(i)
  print(x$gds_indir[i])
  dracarys::umccr_tidy(in_dir = x$gds_indir[i], out_dir = x$outdir[i], prefix = x$sbj[i], dryrun = dryrun, token = token)
  # print(x$local_indir[i])
  # dracarys::umccr_tidy(in_dir = x$local_indir[i], out_dir = x$outdir[i], prefix = x$sbj[i], dryrun = dryrun, token = token)
}
