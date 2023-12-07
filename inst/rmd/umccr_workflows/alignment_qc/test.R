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

dp_workflow_read <- function() {
  # assertthat::assert_that(!is.null(sbj), all(grepl("^SBJ", sbj)))
  RAthena::RAthena_options(clear_s3_resource = FALSE)
  con <- DBI::dbConnect(
    RAthena::athena(),
    work_group = "data_portal",
    rstudio_conn_tab = FALSE
  )
  # q_quote <- shQuote(paste(glue("umccr__automated__umccrise__{sbj}"), collapse = "|"))
  q1 <- glue(
    #' SELECT * FROM "data_portal"."data_portal"."data_portal_workflow" where REGEXP_LIKE("wfr_name", {q_quote});'
    'SELECT * FROM "data_portal"."data_portal"."data_portal_workflow" WHERE "type_name" = \'wgs_alignment_qc\' ORDER BY "start" DESC LIMIT 10;'
  )
  d <- RAthena::dbGetQuery(con, q1) |>
    tibble::as_tibble()
  d |>
    dracarys::meta_wgs_alignment_qc(status = "Succeeded")
}

p <- dp_workflow_read()

d <- p |>
  dplyr::select(SubjectID, LibraryID, SampleID, Lane, start, end, portal_run_id, gds_outdir_dragen)

tidy_script <- system.file("cli/dracarys.R", package = "dracarys")
res <- d |>
  dplyr::rowwise() |>
  dplyr::mutate(
    indir = gds_outdir_dragen,
    outdir = file.path(sub("gds://", "", .data$indir)),
    outdir = file.path(normalizePath("~/icav1/g"), .data$outdir),
    # indir = file.path(outdir, "dracarys_gds_sync"), # for when debugging locally
    cmd = system(glue::glue("{tidy_script} tidy --in_dir {.data$indir} --out_dir {.data$outdir} --prefix {.data$SampleID} --format rds"))
  ) |>
  dplyr::ungroup()
