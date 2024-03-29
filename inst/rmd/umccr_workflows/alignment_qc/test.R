#!/usr/bin/env Rscript

{
  require(DBI)
  require(dplyr)
  require(dracarys)
  require(fs)
  require(glue)
  require(here)
  require(RAthena)
  require(readr)
}

# log into aws umccr prod account
dracarys:::awsvault_profile("upro")

dp_workflow_read <- function(start_date) {
  RAthena::RAthena_options(clear_s3_resource = FALSE)
  con <- DBI::dbConnect(
    RAthena::athena(),
    work_group = "data_portal",
    rstudio_conn_tab = FALSE
  )
  q1 <- glue(
    'SELECT * FROM "data_portal"."data_portal"."data_portal_workflow" ',
    'WHERE "start" > date(\'{start_date}\') AND REGEXP_LIKE("type_name", \'alignment_qc\') ',
    'ORDER BY "start" DESC;'
  )
  d <- RAthena::dbGetQuery(con, q1) |>
    tibble::as_tibble()
  DBI::dbDisconnect(con)
  d
}

# TODO: use dracarys::glims_portal_read
dp_lims_read <- function(libids) {
  assertthat::assert_that(!is.null(libids), all(grepl("^L", libids)))
  libids <- unique(libids)
  RAthena::RAthena_options(clear_s3_resource = FALSE)
  con <- DBI::dbConnect(
    RAthena::athena(),
    work_group = "data_portal",
    rstudio_conn_tab = FALSE
  )
  q_quote <- shQuote(paste(libids, collapse = "|"))
  q1 <- glue(
    'SELECT * FROM "data_portal"."data_portal"."data_portal_limsrow" where REGEXP_LIKE("library_id", {q_quote});'
  )
  d <- RAthena::dbGetQuery(con, q1) |>
    tibble::as_tibble()
  DBI::dbDisconnect(con)
  d
}

# first read in the workflows table, extract metadata, then join with lims
start_date <- "2024-03-23"
p_raw <- dp_workflow_read(start_date)

wgs <- p_raw |>
  filter(type_name == "wgs_alignment_qc") |>
  dracarys::meta_wgs_alignment_qc(status = "Succeeded")
wts <- p_raw |>
  filter(type_name == "wts_alignment_qc") |>
  dracarys::meta_wts_alignment_qc(status = "Succeeded")
p <- dplyr::bind_rows(wgs, wts)
lims <- dp_lims_read(p$LibraryID)

d <- p |>
  dplyr::left_join(lims, by = c("SubjectID" = "subject_id", "LibraryID" = "library_id")) |>
  dplyr::select(
    SubjectID, LibraryID, SampleID, Lane, phenotype, type, source,
    quality, assay, external_subject_id, project_name, project_owner,
    workflow, start, end, portal_run_id, gds_outdir_dragen
  )

tidy_script <- system.file("cli/dracarys.R", package = "dracarys")
meta <- d |>
  dplyr::rowwise() |>
  dplyr::mutate(
    indir = gds_outdir_dragen,
    outdir = file.path(sub("gds://", "", .data$indir)),
    outdir = file.path(normalizePath("~/icav1/g"), .data$outdir),
    # indir = file.path(outdir, "dracarys_gds_sync"), # for when debugging locally
    cmd = system(glue::glue("{tidy_script} tidy --in_dir {.data$indir} --out_dir {.data$outdir} --prefix {.data$SampleID} --format rds"))
  ) |>
  dplyr::ungroup()

meta |>
  saveRDS(here::here(glue::glue("inst/rmd/umccr_workflows/alignment_qc/nogit/meta/{start_date}_wgts.rds")))
