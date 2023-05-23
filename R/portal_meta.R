#' Metadata for wgs align-qc workflow
#'
#' @param status Workflow status to keep (default: Succeeded).
#' @param pmeta Portal workflows metadata table in TSV format.
#'
#' @return A tibble with metadata per wgs align-qc workflow run.
#' @examples
#' \dontrun{
#' pmeta <- here::here("nogit/data_portal/2023-05-21_workflows.csv")
#' meta_wgs_alignment_qc(pmeta)
#' }
#' @export
meta_wgs_alignment_qc <- function(pmeta, status = "Succeeded") {
  wf <- readr::read_csv(pmeta) |>
    dplyr::filter(.data$type_name == "wgs_alignment_qc") |>
    dplyr::select(-c("sample_name", "type_name", "notified", "partition_name")) |>
    dplyr::filter(.data$end_status %in% status)

  d <- wf |>
    dplyr::rowwise() |>
    dplyr::mutate(
      i1 = list(jsonlite::fromJSON(.data$input)),
      o1 = list(jsonlite::fromJSON(.data$output))
    ) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      libid = purrr::map_chr(.data$i1, list("fastq_list_rows", "rglb")),
      rgsm = purrr::map_chr(.data$i1, list("fastq_list_rows", "rgsm")),
      outdir_alignqc = purrr::map_chr(.data$o1, list("dragen_alignment_output_directory", "location")),
      outdir_multiqc = purrr::map_chr(.data$o1, list("multiqc_output_directory", "location")),
      subjectid = sub("umccr__automated__wgs_alignment_qc__(SBJ.*)__L.*", "\\1", .data$wfr_name),
      libid2 = sub("umccr__automated__wgs_alignment_qc__SBJ.*__(L.*)__.*", "\\1", .data$wfr_name),
      lane = sub(".*__(.*)", "\\1", .data$libid2),
    ) |>
    dplyr::select(
      SubjectID = "subjectid",
      LibraryID = "libid",
      Lane = "lane",
      SampleID = "rgsm",
      "wfr_name", "wfr_id", "version", "sequence_run_id", "batch_run_id",
      "start", "end", "portal_run_id", "outdir_alignqc", "outdir_multiqc"
    )
  d
}

#' Metadata for cttso workflow
#'
#' @param status Workflow status to keep (default: Succeeded).
#' @param pmeta Portal workflows metadata table in TSV format.
#'
#' @return A tibble with metadata per cttso workflow run.
#' @examples
#' \dontrun{
#' pmeta <- here::here("nogit/data_portal/2023-05-07_workflows_5f43da12-0b6a-41ce-86c9-9ee62df8792e.csv")
#' cttso_metadata(pmeta)
#' }
#' @export
cttso_metadata <- function(pmeta, status = c("Succeeded")) {
  # retrieve successful cttso workflow runs
  wf <- readr::read_csv(pmeta) |>
    dplyr::filter(.data$type_name == "tso_ctdna_tumor_only") |>
    dplyr::select(-c("sample_name", "type_name", "notified", "partition_name")) |>
    dplyr::filter(.data$end_status %in% dplyr::all_of(status))

  # grab libid/sampleid from the input meta, and outdir from the output meta
  d <- wf |>
    dplyr::rowwise() |>
    dplyr::mutate(
      i1 = list(jsonlite::fromJSON(.data$input)),
      o1 = list(jsonlite::fromJSON(.data$output))
    ) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      tso500_samples = purrr::map(.data$i1, "tso500_samples"),
      outdir = purrr::map(.data$o1, list("output_results_dir", "location"))
    ) |>
    tidyr::unnest(.data$tso500_samples) |>
    dplyr::mutate(
      libid1 = sub(".*_(L.*)", "\\1", .data$sample_id),
      rerun = grepl("rerun", .data$libid1),
      subjectid = sub("umccr__automated__tso_ctdna_tumor_only__(SBJ.*)__L.*", "\\1", .data$wfr_name),
      libid = sub("umccr__automated__tso_ctdna_tumor_only__SBJ.*__(L.*)__.*", "\\1", .data$wfr_name) # equal to libid1 wo _rerun
    ) |>
    tidyr::unnest(.data$outdir) |>
    dplyr::select(
      SubjectID = "subjectid",
      LibraryID = "libid",
      LibraryID_w_rerun = "libid1",
      SampleID = "sample_name",
      gds_indir = "outdir",
      "wfr_name", "wfr_id", "version", "sequence_run_id", "batch_run_id",
      "start", "end", "portal_run_id", "rerun"
    )
  d
}
