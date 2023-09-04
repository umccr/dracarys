require(dracarys)
require(here)
require(glue)
require(dplyr)
require(readr)

# read last 1000 umccrise runs from portal
# 475 from 2022-01-24 until 2023-09-03, of which 449 Succeeded
date1 <- "2023-09-04"
pmeta_raw_rds <- here(glue("nogit/umccrise/rds/portal_meta/{date1}_pmeta_raw.rds"))
# pmeta_raw <- dracarys::portal_meta_read(rows = 1000, params = "&type_name=umccrise")
# saveRDS(pmeta_raw, file = pmeta_raw_rds)
pmeta <- readr::read_rds(pmeta_raw_rds) |>
  dracarys::meta_umccrise(status = "Succeeded")
lims_raw_rds <- here(glue("nogit/umccrise/rds/lims/{date1}_lims_raw.rds"))
# lims_raw <- dracarys::glims_read()
# saveRDS(lims_raw, file = lims_raw_rds)
lims_raw <- readr::read_rds(lims_raw_rds)
lims <- lims_raw |>
  filter(Type == "WGS") |>
  filter(LibraryID %in% c(pmeta$LibraryID_normal, pmeta$LibraryID_tumor))
table(pmeta$LibraryID_tumor %in% lims$LibraryID)
table(pmeta$LibraryID_normal %in% lims$LibraryID)

# The final results sit under gds_outdir_umccrise/<SubjectID>__<SampleID_tumor>/
# We need to get the SampleID_tumor for runs before 2023-04-07. We can do that
# by using the LibraryID_tumor to match up with the glims.
missing_tumor_sampleid <- pmeta |>
  filter(end < "2023-04-07") |>
  pull(LibraryID_tumor)

table(missing_tumor_sampleid %in% lims$LibraryID)
libid2sampid <- lims |>
  filter(LibraryID %in% missing_tumor_sampleid) |>
  select(LibraryID_tumor = LibraryID, SampleID_tumor = SampleID)

d <- pmeta |>
  left_join(libid2sampid, by = "LibraryID_tumor") |>
  mutate(SampleID_tumor = if_else(is.na(SampleID_tumor.x), SampleID_tumor.y, SampleID_tumor.x)) |>
  select(-c(SampleID_tumor.x, SampleID_tumor.y)) |>
  relocate(SampleID_tumor, .before = SampleID_normal) |>
  mutate(gds_outdir_umccrise = glue("{.data$gds_outdir_umccrise}/{.data$SubjectID}__{.data$SampleID_tumor}"))
d

# final portal meta for umccrise runs
saveRDS(d, file = here(glue("nogit/umccrise/rds/portal_meta/{date1}_pmeta_final.rds")))
