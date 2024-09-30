#!/usr/bin/env Rscript

{
  require(dplyr)
  require(dracarys, include.only = "umccr_tidy")
  require(glue, include.only = "glue")
  require(here, include.only = "here")
  require(rportal, include.only = c("portaldb_query_workflow"))
}

# make sure you have logged into AWS and ICA
c("AWS_ACCESS_KEY_ID", "AWS_SECRET_ACCESS_KEY", "AWS_REGION") |>
  rportal::envvar_defined() |>
  stopifnot()
icav1_token <- Sys.getenv("ICA_ACCESS_TOKEN") |>
  dracarys::ica_token_validate()
# this helps keep annoying reticulate prompt away
Sys.setenv(RETICULATE_PYTHON = Sys.getenv("CONDA_PYTHON_EXE"))

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
start_date <- "2024-09-27"
p_raw <- query_workflow_alignqc(start_date)

wgs <- p_raw |>
  rportal::meta_wgs_alignment_qc(status = "Succeeded")
wts <- p_raw |>
  rportal::meta_wts_alignment_qc(status = "Succeeded")
p <- bind_rows(wgs, wts)
lims_raw <- query_limsrow_libids(p$LibraryID)

lims <- lims_raw |>
  tidyr::separate_wider_delim(
    library_id,
    delim = "_", names = c("library_id", "topup_or_rerun"), too_few = "align_start"
  ) |>
  select(
    subject_id, library_id, sample_id, sample_name,
    external_subject_id, external_sample_id,
    project_name, project_owner, phenotype, type,
    source, assay, quality, workflow
  ) |>
  distinct()

d <- p |>
  left_join(lims, by = c("SubjectID" = "subject_id", "LibraryID" = "library_id")) |>
  select(
    "SubjectID", "LibraryID", "SampleID", "lane", "phenotype", "type", "source",
    "assay", "workflow", "external_subject_id", "project_name", "project_owner",
    "start", "end", "portal_run_id", "gds_outdir_dragen", "fq1", "fq2"
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
        "--token {icav1_token} ",
        "--format rds"
      )
    )
  ) |>
  ungroup()

meta |>
  saveRDS(here(glue("inst/rmd/umccr_workflows/alignment_qc/nogit/meta/{start_date}_wgts.rds")))
