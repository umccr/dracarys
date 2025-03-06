#!/usr/bin/env Rscript

{
  require(dplyr)
  require(dracarys, include.only = "umccr_tidy")
  require(glue, include.only = "glue")
  require(here, include.only = "here")
  require(rportal, include.only = c("orca_workflow_list"))
  require(stringr, include.only = "str_remove_all")
  require(tidyr, include.only = "unnest")
  require(fs, include.only = "dir_create")
}

# make sure you have logged into AWS
c("AWS_ACCESS_KEY_ID", "AWS_SECRET_ACCESS_KEY", "AWS_REGION") |>
  rportal::envvar_defined() |>
  stopifnot()
token <- rportal::orca_jwt() |>
  rportal::jwt_validate()
date_start <- "2024-12-21"
date_end <- "2025-01-15"
dates <- seq(from = as.Date(date_start), to = as.Date(date_end), by = "day") |>
  stringr::str_remove_all("-") |>
  paste(collapse = "|")
wf0 <- rportal::orca_workflow_list(
  wf_name = "wgts-qc",
  token = token,
  page_size = 50
)
# get pld
wf1 <- wf0 |>
  filter(grepl(dates, .data$portalRunId)) |>
  rowwise() |>
  mutate(
    pld = list(rportal::orca_wfrid2payload(
      wfrid = .data$orcabusId,
      token = token
    ))
  ) |>
  ungroup()
# tidy pld
wf2 <- wf1 |>
  rowwise() |>
  mutate(pld_tidy = list(rportal::pld_wgtsqc(.data$pld))) |>
  ungroup() |>
  select(
    workflowRunId = "orcabusId",
    portalRunId,
    currentStateTimestamp,
    pld_tidy
  ) |>
  tidyr::unnest(pld_tidy)

query_limsrow_libids <- function(libids) {
  stopifnot(!is.null(libids), all(grepl("^L", libids)))
  libids <- unique(libids) |>
    paste(collapse = "|")
  q1 <- glue("WHERE REGEXP_LIKE(\"library_id\", '{libids}');")
  rportal::portaldb_query_limsrow(q1)
}

lims0 <- query_limsrow_libids(wf2$libraryId)

lims1 <- lims0 |>
  tidyr::separate_wider_delim(
    library_id,
    delim = "_",
    names = c("library_id", "topup_or_rerun"),
    too_few = "align_start"
  ) |>
  select(
    individualId = "subject_id",
    libraryId = "library_id",
    sampleId = "sample_id",
    sampleName = "sample_name",
    subjectId = "external_subject_id",
    externalSampleId = "external_sample_id",
    projectName = "project_name",
    projectOwner = "project_owner",
    phenotype,
    type,
    source,
    assay,
    quality,
    workflow
  ) |>
  distinct()

wf_lims <- wf2 |>
  left_join(lims1, by = "libraryId") |>
  select(
    "libraryId",
    "individualId",
    "sampleId",
    "sampleName",
    "subjectId",
    "externalSampleId",
    "projectName",
    "projectOwner",
    lane = "input_lane",
    "phenotype",
    "sampleType",
    date = "currentStateTimestamp",
    "source",
    "assay",
    "quality",
    "workflow",
    "portalRunId",
    "output_dragenAlignmentOutputUri",
    "input_read1FileUri",
    "input_read2FileUri",
  ) |>
  mutate(rownum = row_number()) |>
  relocate("rownum")

# set up progress bar for the dtw function
nticks <- nrow(wf_lims)
bar_width <- 50
pb <- progress::progress_bar$new(
  format = "[:bar] :current/:total (:percent) elapsed :elapsedfull eta :eta",
  total = nticks,
  clear = FALSE,
  show_after = 0,
  width = bar_width
)
# wrapping the dtw function to use the progress bar
fun1 <- function(path, prefix, outdir) {
  pb$tick(0)
  res <- dracarys::dtw_Wf_dragen(
    path = path,
    prefix = prefix,
    outdir = outdir,
    format = "rds",
    max_files = 1000,
    dryrun = FALSE
  )
  pb$tick()
  return(res)
}

data_tidy <- wf_lims |>
  rowwise() |>
  mutate(
    indir = .data$output_dragenAlignmentOutputUri,
    outdir = file.path(sub("s3://", "", .data$indir)),
    outdir = fs::as_fs_path(file.path(normalizePath("~/s3"), .data$outdir))
    # indir = outdir # for when debugging locally
  ) |>
  mutate(
    data_tidy = list(
      fun1(
        path = .data$indir,
        prefix = .data$libraryId,
        outdir = .data$outdir
      )
    )
  ) |>
  ungroup()

outdir1 <- fs::dir_create("inst/reports/wgts-qc/nogit/tidy_data_rds")
date1 <- "2025-01-14"
data_tidy |>
  saveRDS(here(glue("{outdir1}/{date1}_wgts.rds")))

#---- for debugging/changing parsers ----#
data_tidy <- readRDS(here(glue("{outdir1}/{date1}_wgts.rds")))
data_tidy2 <- data_tidy |>
  select(-c(indir, outdir, data_tidy)) |>
  rowwise() |>
  mutate(
    indir = .data$output_dragenAlignmentOutputUri,
    outdir = file.path(sub("s3://", "", .data$indir)),
    outdir = fs::as_fs_path(file.path(normalizePath("~/s3"), .data$outdir)),
    indir = outdir, # for when debugging locally
  ) |>
  mutate(
    data_tidy = list(
      dracarys::dtw_Wf_dragen(
        path = .data$indir,
        prefix = .data$libraryId,
        outdir = .data$outdir,
        format = "rds",
        max_files = 1000,
        dryrun = FALSE
      )
    )
  ) |>
  ungroup()
data_tidy2 |>
  saveRDS(here(glue("{outdir1}/{date1}_wgts.rds")))
