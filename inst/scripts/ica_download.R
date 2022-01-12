# Download DRAGEN qc metrics files from ICA
require(tidyverse)
require(glue)
require(here)
require(fs)
require(jsonlite)

j <- here("nogit/seqcii/seqc_dragen_outputs.json")
l <- jsonlite::read_json(j)[["items"]]
p <- c("replay.json", "fragment_length_hist.csv",
       "mapping_metrics.csv", "ploidy_estimation_metrics.csv",
       "time_metrics.csv", "vc_metrics.csv",
       "wgs_contig_mean_cov_normal.csv", "wgs_contig_mean_cov_tumor.csv",
       "wgs_coverage_metrics_normal.csv", "wgs_coverage_metrics_tumor.csv",
       "wgs_fine_hist_normal.csv", "wgs_fine_hist_tumor.csv") |>
  paste(collapse = "|")
bind_rows(l) |>
  mutate(gds = paste0("gds://", volumeName, path),
         size = fs::as_fs_bytes(sizeInBytes)) |>
  select(gds, name, size) |>
  filter(grepl(p, gds)) |>
  mutate(
    dname = basename(dirname(gds)),
    dname = sub("SEQC-II_dragen_", "", dname),
    target = file.path(here("nogit/seqcii"), dname, name)) |>
  select(gds, target) |>
  rowwise() |>
  mutate(cmd = system(glue("ica files download {gds} {target}")))
