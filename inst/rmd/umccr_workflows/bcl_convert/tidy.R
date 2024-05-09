#!/usr/bin/env Rscript

{
  require(dplyr)
  require(dracarys, include.only = "MultiqcFile")
  require(glue, include.only = "glue")
  require(here, include.only = "here")
  require(rportal, include.only = c("awsvault_profile"))
}

# log into aws umccr prod account
rportal::awsvault_profile("upro")
query_workflow_bclconvert <- function(start_date) {
  q1 <- glue(
    'WHERE "start" > date(\'{start_date}\') AND REGEXP_LIKE("type_name", \'bcl_convert\') ',
    'ORDER BY "start" DESC;'
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
start_date <- "2024-04-14"

#-- workflows --#
# p_raw <- query_workflow_bclconvert(start_date)
p_rds <- here(glue("nogit/bcl_convert/rds/{start_date}_p.rds"))
# fs::dir_create(dirname(p_rds))
# saveRDS(p_raw, p_rds)
p_raw <- readRDS(p_rds)
p <- p_raw |>
  rportal::meta_bcl_convert(status = "Succeeded")
#-- lims --#
# lims_raw <- query_limsrow_libids(p$LibraryID)
lims_rds <- here(glue("nogit/bcl_convert/rds/{start_date}_lims.rds"))
# fs::dir_create(dirname(lims_rds))
# saveRDS(lims_raw, lims_rds)
lims_raw <- readRDS(lims_rds)
lims <- lims_raw |>
  select(
    SubjectID = "subject_id", LibraryID = "library_id", SampleID = "sample_id",
    "type", "phenotype", "source",
    ProjectOwner = "project_owner", ProjectName = "project_name", "workflow"
  ) |>
  distinct()

d <- p |>
  left_join(lims, by = c("LibraryID", "SampleID")) |>
  select(
    "SubjectID", "LibraryID", "SampleID", "phenotype",
    "runfolder_name", "batch_name",
    "gds_outdir_reports",
    "gds_indir_bcl", "gds_outdir_multiqc", "gds_outdir_multiqc_interop",
    "type", "source",
    "ProjectName", "ProjectOwner",
    "workflow", "start", "end", "portal_run_id",
  )

report_dirs <- d |>
  distinct(gds_outdir_reports)

x <- report_dirs |>
  rowwise() |>
  mutate(
    contents = list(dracarys::gds_files_list(gds_outdir_reports, token = Sys.getenv("ICA_ACCESS_TOKEN"), page_size = 20))
  ) |>
  ungroup()


# tidy_script <- system.file("cli/dracarys.R", package = "dracarys")
# meta <- d |>
#   rowwise() |>
#   mutate(
#     indir = gds_outdir_dragen,
#     outdir = file.path(sub("gds://", "", .data$indir)),
#     outdir = file.path(normalizePath("~/icav1/g"), .data$outdir),
#     # indir = file.path(outdir, "dracarys_gds_sync"), # for when debugging locally
#     cmd = system(glue("{tidy_script} tidy --in_dir {.data$indir} --out_dir {.data$outdir} --prefix {.data$SampleID} --format rds"))
#   ) |>
#   ungroup()
#
# meta |>
#   saveRDS(here(glue("inst/rmd/umccr_workflows/alignment_qc/nogit/meta/{start_date}_wgts.rds")))
