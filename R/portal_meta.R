#' Metadata for bcl_convert workflow
#'
#' @param pmeta Path to portal workflows metadata table, or tibble with already parsed data.
#' @param status Workflow status to keep (default: Succeeded).
#'
#' @return A tibble with metadata per workflow run.
#' @examples
#' pmeta <- system.file("extdata/portal_meta_top4.csv", package = "dracarys")
#' (m <- meta_bcl_convert(pmeta))
#' @testexamples
#' expect_equal(sum(!is.na(m$topup_or_rerun)), 2)
#' expect_equal(length(unique(m$portal_run_id)), 4)
#' @export
meta_bcl_convert <- function(pmeta, status = "Succeeded") {
  # retrieve workflow runs with the given type and status
  wf <- portal_meta_read(pmeta) |>
    dplyr::filter(
      .data$type_name == "bcl_convert",
      .data$end_status %in% status
    )
  d <- wf |>
    meta_io_fromjson() |>
    dplyr::mutate(
      batch_name = purrr::map(.data$input, list("settings_by_samples", "batch_name")),
      samples = purrr::map(.data$input, list("settings_by_samples", "samples")),
      runfolder_name = purrr::map_chr(.data$input, "runfolder_name"),
      gds_outdir_multiqc = purrr::map_chr(.data$output, list("bclconvert_multiqc_out", "location"), .default = NA),
      gds_outdir_multiqc_interop = purrr::map_chr(.data$output, list("interop_multiqc_out", "location"), .default = NA),
    ) |>
    dplyr::rowwise() |>
    dplyr::mutate(
      samples2 = list(purrr::set_names(.data$samples, .data$batch_name)),
      samples2 = list(.data$samples2 |> purrr::list_flatten() |> tibble::enframe(name = "batch_name", value = "sample"))
    ) |>
    dplyr::ungroup() |>
    dplyr::select(-c("samples", "batch_name")) |>
    tidyr::unnest("samples2") |>
    tidyr::unnest("sample")

  # need to take care of following patterns:
  # - sampleid_libid
  # - sampleid_libid_topup(2/3)
  # - sampleid_libid_rerun(2)
  # - sampleidA_sampleidB_libid_..
  # So we just say .*_(L.*) will be libid1, then split that based on topup/rerun
  d |>
    tidyr::separate_wider_regex("sample", c(sampleid = ".*", "_", libid1 = "L.*"), cols_remove = FALSE) |>
    tidyr::separate_wider_regex("libid1", c(libid2 = ".*", "_", topup_or_rerun = ".*"), cols_remove = FALSE, too_few = "align_start") |>
    dplyr::select(
      dplyr::all_of(meta_main_cols()),
      SampleID = "sampleid",
      LibraryID = "libid2",
      "topup_or_rerun",
      "batch_name",
      "runfolder_name",
      "gds_outdir_multiqc",
      "gds_outdir_multiqc_interop",
    )
}

#' Metadata for wts_tumor_only workflow
#'
#' @param pmeta Path to portal workflows metadata table, or tibble with already parsed data.
#' @param status Workflow status to keep (default: Succeeded).
#'
#' @return A tibble with metadata per workflow run.
#' @examples
#' pmeta <- system.file("extdata/portal_meta_top4.csv", package = "dracarys")
#' (m <- meta_wts_tumor_only(pmeta))
#' @testexamples
#' expect_equal(length(unique(m$portal_run_id)), 4)
#' expect_equal(length(unique(m$LibraryID)), 4)
#' @export
meta_wts_tumor_only <- function(pmeta, status = "Succeeded") {
  # retrieve workflow runs with the given type and status
  wf <- portal_meta_read(pmeta) |>
    dplyr::filter(
      .data$type_name == "wts_tumor_only",
      .data$end_status %in% status
    )
  d <- wf |>
    meta_io_fromjson() |>
    dplyr::mutate(
      rglb = purrr::map_chr(.data$input, \(x) unique(x[["fastq_list_rows"]][["rglb"]])),
      rgsm = purrr::map_chr(.data$input, \(x) unique(x[["fastq_list_rows"]][["rgsm"]])),
      lane = purrr::map_chr(.data$input, \(x) paste(x[["fastq_list_rows"]][["lane"]], collapse = ",")),
      lane = as.character(.data$lane),
      gds_outdir_dragen = purrr::map_chr(.data$output, list("dragen_transcriptome_output_directory", "location"), .default = NA),
      gds_outdir_multiqc = purrr::map_chr(.data$output, list("multiqc_output_directory", "location"), .default = NA),
      gds_outdir_arriba = purrr::map_chr(.data$output, list("arriba_output_directory", "location"), .default = NA),
      gds_outdir_qualimap = purrr::map_chr(.data$output, list("qualimap_output_directory", "location"), .default = NA),
      SubjectID = sub("umccr__automated__wts_tumor_only__(SBJ.*)__L.*", "\\1", .data$wfr_name)
    )
  d |>
    dplyr::select(
      dplyr::all_of(meta_main_cols()),
      "SubjectID",
      LibraryID = "rglb",
      SampleID = "rgsm",
      Lane = "lane",
      "gds_outdir_dragen",
      "gds_outdir_multiqc",
      "gds_outdir_arriba",
      "gds_outdir_qualimap"
    )
}

#' Metadata for wgs_alignment_qc workflow
#'
#' @param pmeta Path to portal workflows metadata table, or tibble with already parsed data.
#' @param status Workflow status to keep (default: Succeeded).
#'
#' @return A tibble with metadata per workflow run.
#' @examples
#' pmeta <- system.file("extdata/portal_meta_top4.csv", package = "dracarys")
#' (m <- meta_wgs_alignment_qc(pmeta))
#' @testexamples
#' expect_equal("Lane" %in% colnames(m), TRUE)
#' @export
meta_wgs_alignment_qc <- function(pmeta, status = "Succeeded") {
  # retrieve workflow runs with the given type and status
  wf <- portal_meta_read(pmeta) |>
    dplyr::filter(
      .data$type_name == "wgs_alignment_qc",
      .data$end_status %in% status
    )
  d <- wf |>
    meta_io_fromjson() |>
    dplyr::mutate(
      rglb = purrr::map_chr(.data$input, list("fastq_list_rows", "rglb")),
      rgsm = purrr::map_chr(.data$input, list("fastq_list_rows", "rgsm")),
      lane = purrr::map_int(.data$input, list("fastq_list_rows", "lane")),
      lane = as.character(.data$lane),
      gds_outdir_dragen = purrr::map_chr(.data$output, list("dragen_alignment_output_directory", "location")),
      gds_outdir_multiqc = purrr::map_chr(.data$output, list("multiqc_output_directory", "location")),
      SubjectID = sub("umccr__automated__wgs_alignment_qc__(SBJ.*)__L.*", "\\1", .data$wfr_name),
    )
  d |>
    dplyr::select(
      dplyr::all_of(meta_main_cols()),
      "SubjectID",
      LibraryID = "rglb",
      SampleID = "rgsm",
      Lane = "lane",
      "gds_outdir_dragen",
      "gds_outdir_multiqc",
    )
}

#' Metadata for tso_ctdna_tumor_only workflow
#'
#' @param pmeta Path to portal workflows metadata table, or tibble with already parsed data.
#' @param status Workflow status to keep (default: Succeeded).
#'
#' @return A tibble with metadata per workflow run.
#' @examples
#' pmeta <- system.file("extdata/portal_meta_top4.csv", package = "dracarys")
#' (m <- meta_tso_ctdna_tumor_only(pmeta))
#' @testexamples
#' expect_equal(length(unique(m$portal_run_id)), 4)
#' @export
meta_tso_ctdna_tumor_only <- function(pmeta, status = c("Succeeded")) {
  # retrieve workflow runs with the given type and status
  wf <- portal_meta_read(pmeta) |>
    dplyr::filter(
      .data$type_name == "tso_ctdna_tumor_only",
      .data$end_status %in% status
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
      dplyr::all_of(meta_main_cols()),
      SubjectID = "subjectid",
      LibraryID = "libid",
      SampleID = "sample_name2",
      "gds_outdir",
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
