require(dracarys)
require(rportal, include.only = "portaldb_query_workflow")
require(here, include.only = "here")
require(glue, include.only = "glue")
require(readr, include.only = "read_rds")
require(dplyr)
require(tidyr, include.only = "separate_wider_delim")

start_date <- "2024-08-29"
query_workflow_umccrise <- function(start_date) {
  q1 <- glue(
    "WHERE \"type_name\" = 'umccrise' AND  \"start\" > date(\'{start_date}\') ",
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

# p_raw <- query_workflow_umccrise(start_date)
p_raw_rds <- here(glue("nogit/data_portal/workflows/{start_date}.rds"))
# saveRDS(p_raw, file = p_raw_rds)
p_raw <- readr::read_rds(p_raw_rds)

p <- p_raw |>
  rportal::meta_umccrise(status = "Succeeded")
# lims_raw <- query_limsrow_libids(p$LibraryID_tumor)
lims_raw_rds <- here(glue("nogit/data_portal/lims/{start_date}.rds"))
# saveRDS(lims_raw, file = lims_raw_rds)
# L2100192 is L2100192_rerun in the lims, 15 libs are rerun/topup/topup2
lims_raw <- readr::read_rds(lims_raw_rds)
lims <- lims_raw |>
  tidyr::separate_wider_delim(
    library_id,
    delim = "_", names = c("library_id", "topup_or_rerun"), too_few = "align_start"
  ) |>
  select(
    subject_id, library_id, sample_id, sample_name,
    external_subject_id, external_sample_id,
    project_name, project_owner,
    source, quality
  ) |>
  distinct()
table(lims$library_id %in% p$LibraryID_tumor) # double-check

d <- p |>
  left_join(lims, by = c("LibraryID_tumor" = "library_id")) |>
  mutate(gds_outdir_umccrise = glue("{.data$gds_outdir_umccrise}/{.data$SubjectID}__{.data$SampleID_tumor}")) |>
  select(
    wfr_id, version, end_status, start, end, portal_run_id, SubjectID, LibraryID_tumor, LibraryID_normal,
    SampleID_tumor, SampleID_normal, gds_outdir_umccrise, gds_indir_dragen_somatic, external_subject_id, external_sample_id,
    project_owner, project_name, source, quality
  )
d

saveRDS(d, file = here(glue("nogit/data_portal/workflows/umccrise_tidy_{start_date}.rds")))
