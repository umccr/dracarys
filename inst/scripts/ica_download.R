# List/Download files from ICA for testing
require(fs)
require(here)
require(httr)
require(jsonlite)
require(tidyverse)
require(dracarys)

#---- functions ----#

dr_download_multiqc <- function(gdsdir, outdir, token = Sys.getenv("ICA_TOKEN_PROD")) {
  fs::dir_create(outdir)
  d <- gds_files_list(gdsdir = gdsdir, token = token) |>
    dplyr::mutate(type = purrr::map_chr(bname, dracarys::match_regex)) |>
    dplyr::select(.data$dname, .data$type, .data$size, .data$path, .data$bname)

  # download dracarys files to outdir/{dname}.json
  d |>
    dplyr::filter(.data$type == "multiqc") |>
    dplyr::mutate(out = file.path(outdir, glue("{dname}.json"))) |>
    dplyr::rowwise() |>
    dplyr::mutate(cmd = gds_file_download(path, out, token))
}

#---- download ----#
#---- dragen wgs validation ----#
samples_wgs <- c("2016.249.18.WH.P025", "SBJ00303") # SEQC50
for (sname in samples_wgs) {
  dr_download(
    gdsdir = glue("gds://development/validation_data/wgs/{sname}/analysis/dragen_somatic/"),
    outdir = here(glue("nogit/wgs/dragen/{sname}")),
    token = Sys.getenv("ICA_TOKEN_PROD")
  )
}

#---- tso ----#
samples_tso <- c(
  # "SBJ00595", "SBJ00998", "SBJ00999" ,"SBJ01001", "SBJ01003", "SBJ01032",
  # "SBJ01040", "SBJ01043", "SBJ01046", "SBJ01059", "SBJ01131", "SBJ01132",
  # "SBJ01133", "SBJ01134", "SBJ01135", "SBJ01136", "SBJ01137", "SBJ01138",
  # "SBJ01139", "SBJ01140", "SBJ01141", "SBJ01142", "SBJ01143", "SBJ01144",
  # "SBJ01146", "SBJ01150", "SBJ01156", "SBJ01157", "SBJ01158", "SBJ01159",
  # "SBJ01160", "SBJ01162", "SBJ01163", "SBJ01164", "SBJ01195", "SBJ01196",
  # "SBJ01197", "SBJ01211", "SBJ01212", "SBJ01213", "SBJ01214", "SBJ01215",
  # "SBJ01216", "SBJ01217", "SBJ01218", "SBJ01219", "SBJ01220", "SBJ01221",
  # "SBJ01222", "SBJ01223", "SBJ01224", "SBJ02679", "SBJ02299", "SBJ02298",
  # "SBJ02096", "SBJ02067", "SBJ02559", "SBJ02558"
)
for (sname in samples_tso) {
  dr_download(
    gdsdir = glue("gds://production/analysis_data/{sname}/tso_ctdna_tumor_only/"),
    outdir = here(glue("nogit/tso/{sname}")),
    token = Sys.getenv("ICA_TOKEN_PROD")
  )
}

#---- dragen wts ----#
samples_wts <- c(
  "SBJ02298", "SBJ02299", "SBJ02300"
)
for (sname in samples_wts) {
  dr_download(
    gdsdir = glue("gds://production/analysis_data/{sname}/wts_tumor_only/"),
    outdir = here(glue("nogit/wts/dragen/{sname}")),
    token = Sys.getenv("DRACARYS_TOKEN_PROD")
  )
}

#---- multiqc ----#
samples <- c(
  "SBJ02402",
  "SBJ02403",
  "SBJ02404",
  "SBJ02405",
  "SBJ02406",
  "SBJ02407"
)

for (sname in samples) {
  workflow <- "wts_tumor_only"
  dr_download_multiqc(
    gdsdir = glue("gds://production/analysis_data/{sname}/{workflow}/"),
    outdir = here(glue("nogit/multiqc/dragen/{workflow}/{sname}")),
    token = Sys.getenv("DRACARYS_TOKEN_PROD")
  )
}

for (sname in samples) {
  workflow <- "wgs_alignment_qc"
  dr_download_multiqc(
    gdsdir = glue("gds://production/analysis_data/{sname}/{workflow}/"),
    outdir = here(glue("nogit/multiqc/dragen/{workflow}/{sname}")),
    token = Sys.getenv("DRACARYS_TOKEN_PROD")
  )
}

for (sname in samples) {
  workflow <- "wgs_tumor_normal"
  dr_download_multiqc(
    gdsdir = glue("gds://production/analysis_data/{sname}/{workflow}/"),
    outdir = here(glue("nogit/multiqc/dragen/{workflow}/{sname}")),
    token = Sys.getenv("DRACARYS_TOKEN_PROD")
  )
}

for (sname in samples) {
  workflow <- "umccrise"
  dr_download_multiqc(
    gdsdir = glue("gds://production/analysis_data/{sname}/{workflow}/"),
    outdir = here(glue("nogit/multiqc/dragen/{workflow}/{sname}")),
    token = Sys.getenv("DRACARYS_TOKEN_PROD")
  )
}
