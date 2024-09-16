#!/usr/bin/env Rscript

{
  require(dplyr)
  require(assertthat, include.only = "assert_that")
  require(dracarys, include.only = "Wf_sash_download_tidy_write")
  require(glue, include.only = "glue")
  require(here, include.only = "here")
  require(rportal, include.only = c("portaldb_query_workflow"))
  require(tidyr, include.only = "separate_wider_delim")
}

query_workflow_sash <- function(start_date, end_date) {
  q1 <- glue(
    "WHERE \"type_name\" = 'sash' ",
    "AND \"start\" >= date(\'{start_date}\') ",
    "AND \"end\" <= date(\'{end_date}\') ",
    "ORDER BY \"start\" DESC;"
  )
  rportal::portaldb_query_workflow(q1)
}

query_limsrow_libids <- function(libids) {
  assertthat::assert_that(!is.null(libids), all(grepl("^L", libids)))
  libids <- unique(libids) |>
    paste(collapse = "|")
  q1 <- glue("WHERE REGEXP_LIKE(\"library_id\", '{libids}');")
  rportal::portaldb_query_limsrow(q1)
}

# first read in the workflows table, extract metadata, then join with lims
start_date <- "2024-08-29"
end_date <- "2024-09-07"
meta_raw <- query_workflow_sash(start_date, end_date)
meta <- meta_raw |>
  rportal::meta_sash()
lims_raw <- query_limsrow_libids(meta$LibraryID_tumor)
lims <- lims_raw |>
  tidyr::separate_wider_delim(
    library_id,
    delim = "_", names = c("library_id", "topup_or_rerun"), too_few = "align_start"
  ) |>
  select(
    subject_id, library_id, sample_id, sample_name,
    external_subject_id, external_sample_id,
    project_name, project_owner,
    source, quality, workflow
  ) |>
  distinct()
table(lims$library_id %in% meta$LibraryID_tumor) # double-check

meta_lims <- meta |>
  left_join(lims, by = c("LibraryID_tumor" = "library_id")) |>
  mutate(rownum = row_number()) |>
  select(
    rownum, wfr_id, version, end_status, start, end, portal_run_id, SubjectID, LibraryID_tumor, LibraryID_normal,
    SampleID_tumor, SampleID_normal, s3_outdir_sash, external_subject_id, external_sample_id,
    project_owner, project_name, source, quality, workflow
  )
meta_lims |>
  saveRDS(here(glue("inst/rmd/umccr_workflows/sash/nogit/meta/{start_date}_{end_date}.rds")))

d <- meta_lims |>
  rowwise() |>
  mutate(
    # indir = .data$s3_outdir_sash,
    outdir = file.path(sub("s3://", "", .data$indir)),
    outdir = file.path(normalizePath("~/s3"), .data$outdir),
    indir = outdir, # for when debugging locally
    res = list(
      dracarys::Wf_sash_download_tidy_write(
        path = .data$indir, SubjectID = .data$SubjectID, SampleID_tumor = .data$SampleID_tumor,
        outdir = .data$outdir, max_files = 1000, dryrun = FALSE
      )
    )
  ) |>
  ungroup()

d |>
  saveRDS(here(glue("inst/rmd/umccr_workflows/sash/nogit/results_{start_date}_{end_date}.rds")))
