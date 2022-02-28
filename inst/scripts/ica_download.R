# Download DRAGEN qc metrics files from ICA
require(fs)
require(glue)
require(here)
require(httr)
require(jsonlite)
require(tidyverse)

# dragen wgs
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

# tso
tso_download_ica <- function(sname) {
  base_url <- "https://aps2.platform.illumina.com/v1"
  path <- glue::glue("/analysis_data/{sname}/tso_ctdna_tumor_only/*")
  token <- Sys.getenv("ICA_ACCESS_TOKEN")
  volname <- "production"
  r1 <- httr::GET(glue("{base_url}/files?volume.name={volname}&path={path}&pageSize=100"),
                  httr::add_headers(Authorization = glue("Bearer {token}")),
                  httr::accept_json())
  parsed <- jsonlite::fromJSON(httr::content(r1, "text"), simplifyVector = FALSE)
  ftypes <- c("AlignCollapseFusionCaller_metrics", "TargetRegionCoverage", "MSI",
              "FragmentLengthHist",  "CombinedVariantOutput",
              "TMB",  "TMBtrace_json", "TMBtrace_tsv",
              "FailedExonCoverageQC_json", "FailedExonCoverageQC_txt",
              "Fusions_csv", "Fusions_json", "SampleAnalysisResults")

  map_df(parsed[["items"]], function(x) c(path = x[["path"]], size = x[["sizeInBytes"]])) |>
    mutate(
      size = fs::as_fs_bytes(size),
      bname = basename(path),
      type = tso_guess_file_type(bname),
      path = glue("gds://{volname}{path}"),
      sample = sname) |>
    filter(type %in% ftypes) |>
    mutate(dname = basename(dirname(path))) |>
    select(sample, dname, type, size, path, bname)
}

samples <- c(
  "SBJ00998", "SBJ00999", "SBJ01001", "SBJ01003", "SBJ01032", "SBJ01040",
  "SBJ01043", "SBJ01046", "SBJ01059", "SBJ01131", "SBJ01132", "SBJ01133",
  "SBJ01134", "SBJ01135", "SBJ01136", "SBJ01137", "SBJ01138", "SBJ01139",
  "SBJ01140", "SBJ01141", "SBJ01142", "SBJ01143", "SBJ01144", "SBJ01146",
  "SBJ01150", "SBJ01156", "SBJ01157", "SBJ01158", "SBJ01159", "SBJ01160",
  "SBJ01162", "SBJ01163", "SBJ01164", "SBJ01195", "SBJ01196", "SBJ01197",
  "SBJ01211", "SBJ01212", "SBJ01213", "SBJ01214", "SBJ01215", "SBJ01216",
  "SBJ01217", "SBJ01218", "SBJ01219", "SBJ01220", "SBJ01221", "SBJ01222",
  "SBJ01223", "SBJ01224"
)

# 13 files * 50 samples = 650 files
d <- map_df(samples, tso_download_ica)
outdir <- here::here("nogit/tso")
d |>
  mutate(out = file.path(outdir, sample, bname)) |>
  rowwise() |>
  mutate(
    cmd = system(glue("ica files download {path} {out}")))

# SBJ00705 will be used as a test example
samples <- "SBJ00705"
map_df(samples, tso_download_ica) |>
  filter(dname == "PTC_ctTSO220118_L2200054") |>
  mutate(out = file.path(outdir, sample, bname)) |>
  rowwise() |>
  mutate(
    cmd = system(glue("ica files download {path} {out}")))
