#' Metadata for wgs_alignment_qc workflow
#'
#' @param pmeta Path to portal workflows metadata table, or tibble with already parsed data.
#' @param status Workflow status to keep (default: Succeeded).
#'
#' @return A tibble with metadata per wgs align-qc workflow run.
#' @examples
#' \dontrun{
#' pmeta <- here::here("nogit/data_portal/2023-05-21_workflows.csv")
#' meta_wgs_alignment_qc(pmeta)
#' }
#' @export
meta_wgs_alignment_qc <- function(pmeta, status = "Succeeded") {
  # retrieve workflow runs with the given type and status
  wf <- portal_meta_read(pmeta) |>
    dplyr::filter(
      .data$type_name == "wgs_alignment_qc",
      .data$end_status %in% dplyr::all_of(status)
    )
  d <- wf |>
    meta_io_fromjson() |>
    dplyr::mutate(
      libid1 = purrr::map_chr(.data$input, list("fastq_list_rows", "rglb")),
      rgsm = purrr::map_chr(.data$input, list("fastq_list_rows", "rgsm")),
      gds_outdir_dragen = purrr::map_chr(.data$output, list("dragen_alignment_output_directory", "location")),
      gds_outdir_multiqc = purrr::map_chr(.data$output, list("multiqc_output_directory", "location")),
      SubjectID = sub("umccr__automated__wgs_alignment_qc__(SBJ.*)__L.*", "\\1", .data$wfr_name),
      libid2 = sub("umccr__automated__wgs_alignment_qc__SBJ.*__(L.*)__.*", "\\1", .data$wfr_name),
      Lane = sub(".*__(.*)", "\\1", .data$libid2),
    )
  d |>
    dplyr::select(
      tidyselect::all_of(meta_main_cols()),
      SubjectID,
      LibraryID = "libid1",
      SampleID = "rgsm",
      Lane,
      gds_outdir_dragen,
      gds_outdir_multiqc,
    )
}

#' Metadata for tso_ctdna_tumor_only workflow
#'
#' @param pmeta Path to portal workflows metadata table, or tibble with already parsed data.
#' @param status Workflow status to keep (default: Succeeded).
#'
#' @return A tibble with metadata per cttso workflow run.
#' @examples
#' \dontrun{
#' pmeta <- here::here("nogit/data_portal/2023-05-21_workflows.csv")
#' meta_tso_ctdna_tumor_only(pmeta)
#' }
#' @export
meta_tso_ctdna_tumor_only <- function(pmeta, status = c("Succeeded")) {
  # retrieve workflow runs with the given type and status
  wf <- portal_meta_read(pmeta) |>
    dplyr::filter(
      .data$type_name == "tso_ctdna_tumor_only",
      .data$end_status %in% dplyr::all_of(status)
    )
  # grab libid/sampleid from the input meta, and outdir from the output meta
  d <- wf |>
    meta_io_fromjson() |>
    dplyr::mutate(
      sample_id = purrr::map_chr(.data$input, list("tso500_samples", "sample_id")),
      sample_name2 = purrr::map_chr(.data$input, list("tso500_samples", "sample_name")),
      gds_outdir = purrr::map_chr(.data$output, list("output_results_dir", "location")),
      libid1 = sub(".*_(L.*)", "\\1", .data$sample_id),
      rerun = grepl("rerun", .data$libid1),
      subjectid = sub("umccr__automated__tso_ctdna_tumor_only__(SBJ.*)__L.*", "\\1", .data$wfr_name),
      libid = sub("umccr__automated__tso_ctdna_tumor_only__SBJ.*__(L.*)__.*", "\\1", .data$wfr_name) # equal to libid1 wo _rerun
    )
  d |>
    dplyr::select(
      tidyselect::all_of(meta_main_cols()),
      SubjectID = "subjectid",
      LibraryID = "libid",
      SampleID = "sample_name2",
      gds_outdir,
      cttso_rerun = "rerun"
    )
}

meta_io_fromjson <- function(pmeta) {
  pmeta <- portal_meta_read(pmeta)
  pmeta |>
    dplyr::rowwise() |>
    dplyr::mutate(
      input = list(jsonlite::fromJSON(.data$input)),
      output = list(jsonlite::fromJSON(.data$output))
    ) |>
    dplyr::ungroup()
}

portal_meta_read <- function(pmeta) {
  # if already parsed, just return it
  if (inherits(pmeta, "data.frame")) {
    assertthat::assert_that(
      all(c(
        "wfr_name", "type_name", "start", "end", "input", "output",
        "portal_run_id", "wfr_id", "wfl_id", "wfv_id", "version", "end_status"
      ) %in% colnames(pmeta))
    )
    return(pmeta)
  }
  assertthat::assert_that(file.exists(pmeta))
  # keep all character except start/end
  ctypes <- readr::cols(
    .default = "c",
    start = readr::col_datetime(format = ""),
    end = readr::col_datetime(format = ""),
  )
  readr::read_csv(pmeta, col_types = ctypes)
}

meta_main_cols <- function() {
  c(
    "id", "wfr_name", "wfr_id", "version", "sequence_run_id", "batch_run_id",
    "start", "end", "portal_run_id"
  )
}
