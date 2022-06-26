# Download files from ICA
require(fs)
require(glue)
require(here)
require(httr)
require(jsonlite)
require(tidyverse)

#---- dragen wgs ----#
ica_download_wgs <- function(sname) {
  base_url <- "https://aps2.platform.illumina.com/v1"
  path <- glue::glue("/validation_data/wgs/{sname}/analysis/somatic/*")
  token <- Sys.getenv("ICA_ACCESS_TOKEN")
  volname <- "development"
  res <- httr::GET(glue::glue("{base_url}/files?volume.name={volname}&path={path}&pageSize=100"),
                   httr::add_headers(Authorization = glue::glue("Bearer {token}")),
                   httr::accept_json())
  j <- jsonlite::fromJSON(httr::content(res, "text"), simplifyVector = FALSE)[["items"]]
  d <- purrr::map_df(j, function(x) c(path = x[["path"]], size = x[["sizeInBytes"]]))
  d |>
    dplyr::mutate(
      size = fs::as_fs_bytes(.data$size),
      bname = basename(.data$path),
      type = purrr::map_chr(bname, match_regex),
      path = glue::glue("gds://{volname}{.data$path}"),
      sample = sname,
      dname = basename(dirname(.data$path))
      ) |>
    dplyr::select(.data$sample, .data$dname, .data$type, .data$size, .data$path, .data$bname) |>
    dplyr::filter(!is.na(.data$type))
}

nm <- "SEQC50"
d <- ica_download_wgs(nm)
outdir <- here::here("nogit/wgs")
d |>
  dplyr::mutate(out = file.path(outdir, sample, bname)) |>
  dplyr::rowwise() |>
  dplyr::mutate(
    cmd = system(glue::glue("ica files download {path} {out}")))

#---- tso ----#
ica_download_tso <- function(sname) {
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

#---- dragen wts ----# (---TODO: CHECK THAT THIS WORKS AGAIN---)
ica_download_wgs <- function(sname) {
  base_url <- "https://aps2.platform.illumina.com/v1"
  path <- glue::glue("/analysis_data/{sname}/dragen_wts/*")
  token <- Sys.getenv("ICA_ACCESS_TOKEN")
  volname <- "development"
  res <- httr::GET(glue::glue("{base_url}/files?volume.name={volname}&path={path}&pageSize=100"),
                   httr::add_headers(Authorization = glue::glue("Bearer {token}")),
                   httr::accept_json())
  j <- jsonlite::fromJSON(httr::content(res, "text"), simplifyVector = FALSE)[["items"]]
  d <- purrr::map_df(j, function(x) c(path = x[["path"]], size = x[["sizeInBytes"]]))
  d |>
    dplyr::mutate(
      size = fs::as_fs_bytes(.data$size),
      bname = basename(.data$path),
      type = purrr::map_chr(bname, match_regex),
      path = glue::glue("gds://{volname}{.data$path}"),
      sample = sname,
      dname = basename(dirname(.data$path))
      ) |>
    dplyr::select(.data$sample, .data$dname, .data$type, .data$size, .data$path, .data$bname) |>
    dplyr::filter(!is.na(.data$type))
}
