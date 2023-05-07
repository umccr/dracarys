require(tidyverse)
require(jsonlite)

wf <- read_csv("nogit/data_portal/2023-05-07_workflows_5f43da12-0b6a-41ce-86c9-9ee62df8792e.csv")

#---- cttso (succeeded)----#
wf <- wf |>
  filter(type_name == "tso_ctdna_tumor_only") |>
  select(-c(sample_name, type_name, notified, partition_name)) |>
  filter(end_status == "Succeeded")

wf1 <- wf |>
  rowwise() |>
  mutate(
    i1 = list(jsonlite::fromJSON(input)),
    o1 = list(jsonlite::fromJSON(output))
  ) |>
  ungroup() |>
  mutate(
    tso500_samples = map(i1, "tso500_samples"),
    outdir = map(o1, list("output_results_dir", "location"))
  ) |>
  unnest(tso500_samples) |>
  mutate(
    libid1 = sub(".*_(L.*)", "\\1", .data$sample_id),
    rerun = grepl("rerun", libid1),
    subjectid = sub("umccr__automated__tso_ctdna_tumor_only__(SBJ.*)__L.*", "\\1", wfr_name),
    libid = sub("umccr__automated__tso_ctdna_tumor_only__SBJ.*__(L.*)__.*", "\\1", wfr_name) # equal to libid1 wo _rerun
  ) |>
  unnest(outdir) |>
  select(
    SubjectID = subjectid,
    LibraryID = libid,
    LibraryID_w_rerun = libid1,
    SampleID = sample_name,
    gds_outdir = outdir,
    wfr_name, wfr_id, version, sequence_run_id, batch_run_id, start, end, portal_run_id, rerun
  )

lims <- read_tsv("~/Downloads/Google LIMS - Sheet1.tsv")
table(wf1$LibraryID_w_rerun %in% lims$LibraryID)
lims |> filter(LibraryID %in% wf1$LibraryID_w_rerun)
