#!/usr/bin/env Rscript

{
  require(dplyr)
  require(dracarys, include.only = "umccr_tidy")
  require(glue, include.only = "glue")
  require(here, include.only = "here")
  require(rportal, include.only = c("awsvault_profile"))
}

# log into aws umccr prod account
rportal::awsvault_profile("upro")

query_workflow_alignqc <- function(start_date) {
  wfs <- c("wgs_alignment_qc", "wts_alignment_qc") |>
    shQuote() |>
    paste(collapse = ", ")
  q1 <- glue(
    "WHERE \"type_name\" IN ({wfs}) AND  \"start\" > date(\'{start_date}\') ",
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
start_date <- "2024-05-30"
p_raw <- query_workflow_alignqc(start_date)

wgs <- p_raw |>
  rportal::meta_wgs_alignment_qc(status = "Succeeded")
wts <- p_raw |>
  rportal::meta_wts_alignment_qc(status = "Succeeded")
p <- bind_rows(wgs, wts)
lims <- query_limsrow_libids(p$LibraryID)

d <- p |>
  left_join(lims, by = c("SubjectID" = "subject_id", "LibraryID" = "library_id")) |>
  select(
    "SubjectID", "LibraryID", "SampleID", "Lane", "phenotype", "type", "source",
    "quality", "assay", "external_subject_id", "project_name", "project_owner",
    "workflow", "start", "end", "portal_run_id", "gds_outdir_dragen"
  ) |>
  mutate(rownum = row_number())

tidy_script <- system.file("cli/dracarys.R", package = "dracarys")


meta <- d |>
  relocate(rownum) |>
  rowwise() |>
  mutate(
    indir = gds_outdir_dragen,
    outdir = file.path(sub("gds://", "", .data$indir)),
    outdir = file.path(normalizePath("~/icav1/g"), .data$outdir),
    # indir = file.path(outdir, "dracarys_gds_sync"), # for when debugging locally
    cmd = system(
      glue(
        "echo ---{.data$rownum}--- && ",
        "{tidy_script} tidy --in_dir {.data$indir} ",
        "--out_dir {.data$outdir} --prefix {.data$SampleID} ",
        "--format rds"
      )
    )
  ) |>
  ungroup()

meta |>
  saveRDS(here(glue("inst/rmd/umccr_workflows/alignment_qc/nogit/meta/{start_date}_wgts.rds")))
