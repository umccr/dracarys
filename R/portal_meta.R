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
#' expect_equal(sum(!is.na(m$topup_or_rerun)), 1)
#' expect_equal(length(unique(m$portal_run_id)), 4)
#' @export
meta_bcl_convert <- function(pmeta, status = "Succeeded") {
  # retrieve workflow runs with the given type and status
  type <- "bcl_convert"
  wf <- portal_meta_read(pmeta) |>
    dplyr::filter(
      .data$type_name == type,
      .data$end_status %in% status
    )
  if (nrow(wf) == 0) {
    return(wf)
  }
  d <- wf |>
    meta_io_fromjson() |>
    dplyr::mutate(
      # input
      batch_name = purrr::map(.data$input, list("settings_by_samples", "batch_name")),
      samples = purrr::map(.data$input, list("settings_by_samples", "samples")),
      runfolder_name = purrr::map_chr(.data$input, "runfolder_name"),
      # output
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
      -dplyr::any_of(c("batch_run")), # NA for bcl_convert
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
  type <- "wts_tumor_only"
  wf <- portal_meta_read(pmeta) |>
    dplyr::filter(
      .data$type_name == type,
      .data$end_status %in% status
    )
  if (nrow(wf) == 0) {
    return(wf)
  }
  d <- wf |>
    meta_io_fromjson() |>
    dplyr::mutate(
      # input
      rglb = purrr::map_chr(.data$input, \(x) unique(x[["fastq_list_rows"]][["rglb"]])),
      rgsm = purrr::map_chr(.data$input, \(x) unique(x[["fastq_list_rows"]][["rgsm"]])),
      lane = purrr::map_chr(.data$input, \(x) paste(x[["fastq_list_rows"]][["lane"]], collapse = ",")),
      lane = as.character(.data$lane),
      # output
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

#' Metadata for rnasum workflow
#'
#' @param pmeta Path to portal workflows metadata table, or tibble with already parsed data.
#' @param status Workflow status to keep (default: Succeeded).
#'
#' @return A tibble with metadata per workflow run.
#' @examples
#' pmeta <- system.file("extdata/portal_meta_top4.csv", package = "dracarys")
#' (m <- meta_rnasum(pmeta))
#' @testexamples
#' expect_equal(m$rnasum_dataset[1], "PANCAN")
#' expect_equal(basename(m$gds_outfile_rnasum_html[4]), "MDX230467.RNAseq_report.html")
#' @export
meta_rnasum <- function(pmeta, status = "Succeeded") {
  # retrieve workflow runs with the given type and status
  type <- "rnasum"
  wf <- portal_meta_read(pmeta) |>
    dplyr::filter(
      .data$type_name == type,
      .data$end_status %in% status
    )
  if (nrow(wf) == 0) {
    return(wf)
  }
  d <- wf |>
    meta_io_fromjson() |>
    dplyr::mutate(
      # input
      gds_indir_dragen = purrr::map_chr(.data$input, list("dragen_transcriptome_directory", "location"), .default = NA),
      gds_indir_umccrise = purrr::map_chr(.data$input, list("umccrise_directory", "location"), .default = NA),
      gds_indir_arriba = purrr::map_chr(.data$input, list("arriba_directory", "location"), .default = NA),
      rnasum_sample_name = purrr::map_chr(.data$input, "sample_name", .default = NA),
      rnasum_dataset = purrr::map_chr(.data$input, "dataset", .default = NA),
      rnasum_report_dir = purrr::map_chr(.data$input, "report_directory", .default = NA),
      sbjid1 = sub("(SBJ.*)__L.*", "\\1", .data$rnasum_report_dir),
      libid1 = sub("(SBJ.*)__(L.*)", "\\2", .data$rnasum_report_dir),
      # output
      gds_outfile_rnasum_html = purrr::map_chr(.data$output, list("rnasum_html", "location"), .default = NA),
      gds_outdir_rnasum = purrr::map_chr(.data$output, list("rnasum_output_directory", "location"), .default = NA),
    )
  d |>
    dplyr::select(
      dplyr::all_of(meta_main_cols()),
      -dplyr::any_of(c("sequence_run", "batch_run")), # NA for rnasum
      SubjectID = "sbjid1",
      LibraryID = "libid1",
      SampleID = "rnasum_sample_name",
      "rnasum_dataset",
      "gds_indir_dragen",
      "gds_indir_umccrise",
      "gds_indir_arriba",
      "gds_outdir_rnasum",
      "gds_outfile_rnasum_html"
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
  type <- "wgs_alignment_qc"
  wf <- portal_meta_read(pmeta) |>
    dplyr::filter(
      .data$type_name == type,
      .data$end_status %in% status
    )
  if (nrow(wf) == 0) {
    return(wf)
  }
  d <- wf |>
    meta_io_fromjson() |>
    dplyr::mutate(
      # input
      rglb = purrr::map_chr(.data$input, list("fastq_list_rows", "rglb")),
      rgsm = purrr::map_chr(.data$input, list("fastq_list_rows", "rgsm")),
      lane = purrr::map_int(.data$input, list("fastq_list_rows", "lane")),
      lane = as.character(.data$lane),
      # output
      gds_outdir_dragen = purrr::map_chr(.data$output, list("dragen_alignment_output_directory", "location"), .default = NA),
      gds_outdir_multiqc = purrr::map_chr(.data$output, list("multiqc_output_directory", "location"), .default = NA),
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

#' Metadata for wts_alignment_qc workflow
#'
#' @param pmeta Path to portal workflows metadata table, or tibble with already parsed data.
#' @param status Workflow status to keep (default: Succeeded).
#'
#' @return A tibble with metadata per workflow run.
#' @examples
#' pmeta <- system.file("extdata/portal_meta_top4.csv", package = "dracarys")
#' (m <- meta_wts_alignment_qc(pmeta))
#' @testexamples
#' expect_equal("Lane" %in% colnames(m), TRUE)
#' @export
meta_wts_alignment_qc <- function(pmeta, status = "Succeeded") {
  # retrieve workflow runs with the given type and status
  type <- "wts_alignment_qc"
  wf <- portal_meta_read(pmeta) |>
    dplyr::filter(
      .data$type_name == type,
      .data$end_status %in% status
    )
  if (nrow(wf) == 0) {
    return(wf)
  }
  d <- wf |>
    meta_io_fromjson() |>
    dplyr::mutate(
      # input
      rglb = purrr::map_chr(.data$input, list("fastq_list_rows", "rglb")),
      rgsm = purrr::map_chr(.data$input, list("fastq_list_rows", "rgsm")),
      lane = purrr::map_int(.data$input, list("fastq_list_rows", "lane")),
      lane = as.character(.data$lane),
      # output
      gds_outdir_dragen = purrr::map_chr(.data$output, list("dragen_alignment_output_directory", "location"), .default = NA),
      gds_outdir_multiqc = purrr::map_chr(.data$output, list("multiqc_output_directory", "location"), .default = NA),
      SubjectID = sub("umccr__automated__wts_alignment_qc__(SBJ.*)__L.*", "\\1", .data$wfr_name),
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

#' Metadata for wgs_tumor_normal workflow
#'
#' @param pmeta Path to portal workflows metadata table, or tibble with already parsed data.
#' @param status Workflow status to keep (default: Succeeded).
#'
#' @return A tibble with metadata per workflow run.
#' @examples
#' pmeta <- system.file("extdata/portal_meta_top4.csv", package = "dracarys")
#' (m <- meta_wgs_tumor_normal(pmeta))
#' @testexamples
#' expect_equal("SubjectID" %in% colnames(m), TRUE)
#' @export
meta_wgs_tumor_normal <- function(pmeta, status = "Succeeded") {
  # retrieve workflow runs with the given type and status
  type <- "wgs_tumor_normal"
  wf <- portal_meta_read(pmeta) |>
    dplyr::filter(
      .data$type_name == type,
      .data$end_status %in% status
    )
  extract_fqlist_el <- function(x, fq, y) {
    val <- x[[fq]][[y]]
    if (!is.null(val)) {
      return(paste(unique(val), collapse = ","))
    } else {
      return(NA)
    }
  }
  if (nrow(wf) == 0) {
    return(wf)
  }
  d <- wf |>
    meta_io_fromjson() |>
    dplyr::mutate(
      # input
      SampleID_normal = purrr::map_chr(.data$input, extract_fqlist_el, "fastq_list_rows", "rgsm"),
      SampleID_tumor = purrr::map_chr(.data$input, extract_fqlist_el, "tumor_fastq_list_rows", "rgsm"),
      LibraryID_normal = purrr::map_chr(.data$input, extract_fqlist_el, "fastq_list_rows", "rglb"),
      LibraryID_tumor = purrr::map_chr(.data$input, extract_fqlist_el, "tumor_fastq_list_rows", "rglb"),
      gds_infile_dragen_ref_tar = purrr::map_chr(.data$input, list("reference_tar", "location"), .default = NA),
      # output
      gds_outdir_dragen_germline = purrr::map_chr(.data$output, list("dragen_germline_output_directory", "location"), .default = NA), # NA in old runs
      gds_outdir_dragen_somatic = purrr::map_chr(.data$output, list("dragen_somatic_output_directory", "location"), .default = NA),
      gds_outdir_multiqc = purrr::map_chr(.data$output, list("multiqc_output_directory", "location"), .default = NA),
      gds_outfile_dragen_germline_snv_vcf = purrr::map_chr(.data$output, list("germline_snv_vcf_out", "location"), .default = NA), # NA in old runs
      gds_outfile_dragen_somatic_snv_vcf = purrr::map_chr(.data$output, list("somatic_snv_vcf_out", "location"), .default = NA),
      gds_outfile_dragen_somatic_snv_vcf_hardfilt = purrr::map_chr(.data$output, list("somatic_snv_vcf_hard_filtered_out", "location"), .default = NA),
      gds_outfile_dragen_somatic_sv_vcf = purrr::map_chr(.data$output, list("somatic_structural_vcf_out", "location"), .default = NA),
      SubjectID = sub("umccr__automated__wgs_tumor_normal__(SBJ.....)__L.*", "\\1", .data$wfr_name) # infer from wfr name
    )
  d |>
    dplyr::select(
      dplyr::all_of(meta_main_cols()),
      -dplyr::any_of(c("sequence_run", "batch_run")), # NA for wgs_tumor_normal
      "SubjectID",
      "LibraryID_tumor",
      "LibraryID_normal",
      "SampleID_tumor",
      "SampleID_normal",
      "gds_infile_dragen_ref_tar",
      "gds_outdir_dragen_somatic",
      "gds_outdir_dragen_germline",
      "gds_outdir_multiqc",
      "gds_outfile_dragen_germline_snv_vcf",
      "gds_outfile_dragen_somatic_snv_vcf",
      "gds_outfile_dragen_somatic_snv_vcf_hardfilt",
      "gds_outfile_dragen_somatic_sv_vcf"
    )
}

#' Metadata for umccrise workflow
#'
#' @param pmeta Path to portal workflows metadata table, or tibble with already parsed data.
#' @param status Workflow status to keep (default: Succeeded).
#'
#' @return A tibble with metadata per workflow run.
#' @examples
#' pmeta <- system.file("extdata/portal_meta_top4.csv", package = "dracarys")
#' (m <- meta_umccrise(pmeta))
#' @testexamples
#' expect_equal(all(c("LibraryID_normal", "LibraryID_tumor") %in% colnames(m)), TRUE)
#' @export
meta_umccrise <- function(pmeta, status = "Succeeded") {
  # retrieve workflow runs with the given type and status
  # The input/output json objects changed from 2023-04-07, so need to handle those too
  type <- "umccrise"
  wf <- portal_meta_read(pmeta) |>
    dplyr::filter(
      .data$type_name == type,
      .data$end_status %in% status
    )
  extract_fqlist_el <- function(x, y) {
    val <- x[["fastq_list_rows_germline"]][[y]]
    if (!is.null(val)) {
      return(paste(unique(val), collapse = ","))
    } else {
      return(NA)
    }
  }
  if (nrow(wf) == 0) {
    return(wf)
  }
  d <- wf |>
    meta_io_fromjson() |>
    dplyr::mutate(
      ## input first
      # new runs
      dragen_normal_sampleid2 = purrr::map_chr(.data$input, "dragen_normal_id", .default = NA),
      dragen_tumor_sampleid2 = purrr::map_chr(.data$input, "dragen_tumor_id", .default = NA),
      gds_indir_dragen_somatic = purrr::map_chr(.data$input, list("dragen_somatic_directory", "location"), .default = NA), # same for old
      gds_indir_dragen_germline2 = purrr::map_chr(.data$input, list("dragen_germline_directory", "location"), .default = NA),
      umccrise_outdir_name2 = purrr::map_chr(.data$input, "output_directory_name", .default = NA),
      umccrise_genomes_tar2 = purrr::map_chr(.data$input, list("genomes_tar", "location"), .default = NA),
      sbjid2 = purrr::map_chr(.data$input, "subject_identifier", .default = NA),
      dragen_tumor_libid2 = sub("(L.*)__(L.*)", "\\1", .data$umccrise_outdir_name2),
      dragen_normal_libid2 = sub("(L.*)__(L.*)", "\\2", .data$umccrise_outdir_name2),
      # old runs
      dragen_normal_sampleid1 = purrr::map_chr(.data$input, extract_fqlist_el, "rgsm"),
      dragen_normal_lane1 = purrr::map_chr(.data$input, extract_fqlist_el, "lane"),
      umccrise_outdir_name1 = purrr::map_chr(.data$input, "output_directory_umccrise", .default = NA),
      umccrise_genomes_tar1 = purrr::map_chr(.data$input, list("reference_tar_umccrise", "location"), .default = NA),
      sbjid1 = purrr::map_chr(.data$input, "subject_identifier_umccrise", .default = NA),
      dragen_tumor_libid1 = sub("(L.*)__(L.*)", "\\1", .data$umccrise_outdir_name1),
      dragen_normal_libid1 = sub("(L.*)__(L.*)", "\\2", .data$umccrise_outdir_name1),
      # merge new with old
      SubjectID = dplyr::if_else(is.na(.data$sbjid1), .data$sbjid2, .data$sbjid1),
      LibraryID_normal = dplyr::if_else(
        is.na(.data$dragen_normal_libid1), .data$dragen_normal_libid2, .data$dragen_normal_libid1
      ),
      LibraryID_tumor = dplyr::if_else(
        is.na(.data$dragen_tumor_libid1), .data$dragen_tumor_libid2, .data$dragen_tumor_libid1
      ),
      SampleID_normal = dplyr::if_else(
        is.na(.data$dragen_normal_sampleid1), .data$dragen_normal_sampleid2, .data$dragen_normal_sampleid1
      ),
      SampleID_tumor = .data$dragen_tumor_sampleid2, # NA for old runs
      Lane_normal = .data$dragen_normal_lane1, # NA for new runs
      gds_indir_dragen_germline = .data$gds_indir_dragen_germline2, # NA for old runs
      gds_infile_genomes_tar = dplyr::if_else(
        is.na(.data$umccrise_genomes_tar1), .data$umccrise_genomes_tar2, .data$umccrise_genomes_tar1
      ),
      ## output second
      gds_outdir_umccrise2 = purrr::map_chr(.data$output, list("output_directory", "location"), .default = NA),
      gds_outdir_umccrise1 = purrr::map_chr(.data$output, list("umccrise_output_directory", "location"), .default = NA),
      gds_outdir_umccrise = dplyr::if_else(
        is.na(.data$gds_outdir_umccrise1), .data$gds_outdir_umccrise2, .data$gds_outdir_umccrise1
      )
    )
  d |>
    dplyr::select(
      meta_main_cols(),
      -dplyr::any_of(c("sequence_run", "batch_run")), # NA for umccrise
      "SubjectID",
      "LibraryID_tumor",
      "LibraryID_normal",
      "SampleID_normal",
      "SampleID_tumor",
      "gds_outdir_umccrise",
      "gds_indir_dragen_somatic",
      "gds_indir_dragen_germline",
      "gds_infile_genomes_tar"
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
#' expect_equal(length(unique(m$portal_run_id)), 2)
#' @export
meta_tso_ctdna_tumor_only <- function(pmeta, status = c("Succeeded")) {
  # retrieve workflow runs with the given type and status
  type <- "tso_ctdna_tumor_only"
  wf <- portal_meta_read(pmeta) |>
    dplyr::filter(
      .data$type_name == type,
      .data$end_status %in% status
    )
  if (nrow(wf) == 0) {
    return(wf)
  }
  # grab libid/sampleid from the input meta, and outdir from the output meta
  d <- wf |>
    meta_io_fromjson() |>
    dplyr::mutate(
      # input
      sample_id = purrr::map_chr(.data$input, list("tso500_samples", "sample_id"), .default = NA),
      sample_name2 = purrr::map_chr(.data$input, list("tso500_samples", "sample_name"), .default = NA),
      # output
      gds_outdir = purrr::map_chr(.data$output, list("output_results_dir", "location"), .default = NA),
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

#' Metadata for star_alignment workflow
#'
#' @param pmeta Path to portal workflows metadata table, or tibble with already parsed data.
#' @param status Workflow status to keep (default: Succeeded).
#'
#' @return A tibble with metadata per workflow run.
#' @examples
#' pmeta <- system.file("extdata/portal_meta_top4.csv", package = "dracarys")
#' (m <- meta_star_alignment(pmeta))
#' @testexamples
#' expect_equal(all(c("s3_outdir_star", "LibraryID") %in% colnames(m)), TRUE)
#' @export
meta_star_alignment <- function(pmeta, status = "Succeeded") {
  # retrieve workflow runs with the given type and status
  type <- "star_alignment"
  wf <- portal_meta_read(pmeta) |>
    dplyr::filter(
      .data$type_name == type,
      .data$end_status %in% status
    )
  if (nrow(wf) == 0) {
    return(wf)
  }
  # grab libid/sampleid from the input meta, and outdir from the output meta
  d <- wf |>
    meta_io_fromjson() |>
    dplyr::mutate(
      # input
      SubjectID = purrr::map_chr(.data$input, "subject_id", .default = NA),
      SampleID = purrr::map_chr(.data$input, "sample_id", .default = NA),
      LibraryID = purrr::map_chr(.data$input, "library_id", .default = NA),
      gds_fq_fwd = purrr::map_chr(.data$input, "fastq_fwd", .default = NA),
      gds_fq_rev = purrr::map_chr(.data$input, "fastq_rev", .default = NA),
      # output
      s3_outdir_star = purrr::map_chr(.data$output, "output_directory", .default = NA)
    )
  d |>
    dplyr::select(
      meta_main_cols(),
      "SubjectID",
      "LibraryID",
      "SampleID",
      "s3_outdir_star",
      "gds_fq_fwd",
      "gds_fq_rev"
    )
}


#' Metadata for sash workflow
#'
#' @param pmeta Path to portal workflows metadata table, or tibble with already parsed data.
#' @param status Workflow status to keep (default: Succeeded).
#'
#' @return A tibble with metadata per workflow run.
#' @examples
#' pmeta <- system.file("extdata/portal_meta_top4.csv", package = "dracarys")
#' (m <- meta_sash(pmeta))
#' @testexamples
#' expect_equal(all(c("s3_indir_oncoanalyser", "LibraryID_tumor", "s3_outdir_sash") %in% colnames(m)), TRUE)
#' @export
meta_sash <- function(pmeta, status = "Succeeded") {
  # retrieve workflow runs with the given type and status
  type <- "sash"
  wf <- portal_meta_read(pmeta) |>
    dplyr::filter(
      .data$type_name == type,
      .data$end_status %in% status
    )
  if (nrow(wf) == 0) {
    return(wf)
  }
  # grab libid/sampleid from the input meta, and outdir from the output meta
  d <- wf |>
    meta_io_fromjson() |>
    dplyr::mutate(
      # input
      SubjectID = purrr::map_chr(.data$input, "subject_id", .default = NA),
      SampleID_tumor = purrr::map_chr(.data$input, "tumor_sample_id", .default = NA),
      SampleID_normal = purrr::map_chr(.data$input, "normal_sample_id", .default = NA),
      LibraryID_tumor = purrr::map_chr(.data$input, "tumor_library_id", .default = NA),
      LibraryID_normal = purrr::map_chr(.data$input, "normal_library_id", .default = NA),
      gds_indir_dragen_somatic = purrr::map_chr(.data$input, "dragen_somatic_dir", .default = NA),
      gds_indir_dragen_germline = purrr::map_chr(.data$input, "dragen_germline_dir", .default = NA),
      s3_indir_oncoanalyser = purrr::map_chr(.data$input, "oncoanalyser_dir", .default = NA),
      # output
      s3_outdir_sash = purrr::map_chr(.data$output, "output_directory", .default = NA)
    )
  d |>
    dplyr::select(
      meta_main_cols(),
      "SubjectID",
      "LibraryID_tumor",
      "LibraryID_normal",
      "SampleID_tumor",
      "SampleID_normal",
      "s3_outdir_sash",
      "s3_indir_oncoanalyser",
      "gds_indir_dragen_somatic",
      "gds_indir_dragen_germline"
    )
}

#' Metadata for oncoanalyser_wgs workflow
#'
#' @param pmeta Path to portal workflows metadata table, or tibble with already parsed data.
#' @param status Workflow status to keep (default: Succeeded).
#'
#' @return A tibble with metadata per workflow run.
#' @examples
#' pmeta <- system.file("extdata/portal_meta_top4.csv", package = "dracarys")
#' (m <- meta_oncoanalyser_wgs(pmeta))
#' @testexamples
#' expect_equal(all(c("s3_outdir_oncoanalyser", "LibraryID_tumor", "gds_bam_tumor") %in% colnames(m)), TRUE)
#' @export
meta_oncoanalyser_wgs <- function(pmeta, status = "Succeeded") {
  # retrieve workflow runs with the given type and status
  type <- "oncoanalyser_wgs"
  wf <- portal_meta_read(pmeta) |>
    dplyr::filter(
      .data$type_name == type,
      .data$end_status %in% status
    )
  if (nrow(wf) == 0) {
    return(wf)
  }
  # grab libid/sampleid from the input meta, and outdir from the output meta
  d <- wf |>
    meta_io_fromjson() |>
    dplyr::mutate(
      # input
      SubjectID = purrr::map_chr(.data$input, "subject_id", .default = NA),
      SampleID_tumor = purrr::map_chr(.data$input, "tumor_wgs_sample_id", .default = NA),
      SampleID_normal = purrr::map_chr(.data$input, "normal_wgs_sample_id", .default = NA),
      LibraryID_tumor = purrr::map_chr(.data$input, "tumor_wgs_library_id", .default = NA),
      LibraryID_normal = purrr::map_chr(.data$input, "normal_wgs_library_id", .default = NA),
      gds_bam_tumor = purrr::map_chr(.data$input, "tumor_wgs_bam", .default = NA),
      gds_bam_normal = purrr::map_chr(.data$input, "normal_wgs_bam", .default = NA),
      # output
      s3_outdir_oncoanalyser = purrr::map_chr(.data$output, "output_directory", .default = NA)
    )
  d |>
    dplyr::select(
      meta_main_cols(),
      "SubjectID",
      "LibraryID_tumor",
      "LibraryID_normal",
      "SampleID_tumor",
      "SampleID_normal",
      "s3_outdir_oncoanalyser",
      "gds_bam_tumor",
      "gds_bam_normal"
    )
}

#' Metadata for oncoanalyser_wgts_existing_both workflow
#'
#' @param pmeta Path to portal workflows metadata table, or tibble with already parsed data.
#' @param status Workflow status to keep (default: Succeeded).
#'
#' @return A tibble with metadata per workflow run.
#' @examples
#' pmeta <- system.file("extdata/portal_meta_top4.csv", package = "dracarys")
#' (m <- meta_oncoanalyser_wgts_existing_both(pmeta))
#' @testexamples
#' expect_equal(all(c("s3_outdir_oncoanalyser", "LibraryID_tumor_wts", "gds_bam_tumor_wgs") %in% colnames(m)), TRUE)
#' @export
meta_oncoanalyser_wgts_existing_both <- function(pmeta, status = "Succeeded") {
  # retrieve workflow runs with the given type and status
  type <- "oncoanalyser_wgts_existing_both"
  wf <- portal_meta_read(pmeta) |>
    dplyr::filter(
      .data$type_name == type,
      .data$end_status %in% status
    )
  if (nrow(wf) == 0) {
    return(wf)
  }
  # grab libid/sampleid from the input meta, and outdir from the output meta
  d <- wf |>
    meta_io_fromjson() |>
    dplyr::mutate(
      # input
      SubjectID = purrr::map_chr(.data$input, "subject_id", .default = NA),
      SampleID_tumor_wgs = purrr::map_chr(.data$input, "tumor_wgs_sample_id", .default = NA),
      SampleID_normal_wgs = purrr::map_chr(.data$input, "normal_wgs_sample_id", .default = NA),
      SampleID_tumor_wts = purrr::map_chr(.data$input, "tumor_wts_sample_id", .default = NA),
      LibraryID_tumor_wgs = purrr::map_chr(.data$input, "tumor_wgs_library_id", .default = NA),
      LibraryID_normal_wgs = purrr::map_chr(.data$input, "normal_wgs_library_id", .default = NA),
      LibraryID_tumor_wts = purrr::map_chr(.data$input, "tumor_wts_library_id", .default = NA),
      gds_bam_tumor_wgs = purrr::map_chr(.data$input, "tumor_wgs_bam", .default = NA),
      gds_bam_normal_wgs = purrr::map_chr(.data$input, "normal_wgs_bam", .default = NA),
      s3_bam_tumor_wts = purrr::map_chr(.data$input, "tumor_wts_bam", .default = NA),
      s3_indir_oncoanalyser_wgs = purrr::map_chr(.data$input, "existing_wgs_dir", .default = NA),
      s3_indir_oncoanalyser_wts = purrr::map_chr(.data$input, "existing_wts_dir", .default = NA),
      # output
      s3_outdir_oncoanalyser = purrr::map_chr(.data$output, "output_directory", .default = NA)
    )
  d |>
    dplyr::select(
      meta_main_cols(),
      "SubjectID",
      "LibraryID_tumor_wgs",
      "LibraryID_normal_wgs",
      "LibraryID_tumor_wts",
      "SampleID_tumor_wgs",
      "SampleID_normal_wgs",
      "SampleID_tumor_wts",
      "s3_outdir_oncoanalyser",
      "s3_indir_oncoanalyser_wgs",
      "s3_indir_oncoanalyser_wts",
      "gds_bam_tumor_wgs",
      "gds_bam_normal_wgs",
      "s3_bam_tumor_wts"
    )
}

#' Metadata for oncoanalyser_wts workflow
#'
#' @param pmeta Path to portal workflows metadata table, or tibble with already parsed data.
#' @param status Workflow status to keep (default: Succeeded).
#'
#' @return A tibble with metadata per workflow run.
#' @examples
#' pmeta <- system.file("extdata/portal_meta_top4.csv", package = "dracarys")
#' (m <- meta_oncoanalyser_wts(pmeta))
#' @testexamples
#' expect_equal(all(c("s3_outdir_oncoanalyser", "LibraryID", "s3_bam") %in% colnames(m)), TRUE)
#' @export
meta_oncoanalyser_wts <- function(pmeta, status = "Succeeded") {
  # retrieve workflow runs with the given type and status
  type <- "oncoanalyser_wts"
  wf <- portal_meta_read(pmeta) |>
    dplyr::filter(
      .data$type_name == type,
      .data$end_status %in% status
    )
  if (nrow(wf) == 0) {
    return(wf)
  }
  # grab libid/sampleid from the input meta, and outdir from the output meta
  d <- wf |>
    meta_io_fromjson() |>
    dplyr::mutate(
      # input
      mode = purrr::map_chr(.data$input, "mode", .default = NA),
      SubjectID = purrr::map_chr(.data$input, "subject_id", .default = NA),
      SampleID = purrr::map_chr(.data$input, "tumor_wts_sample_id", .default = NA),
      LibraryID = purrr::map_chr(.data$input, "tumor_wts_library_id", .default = NA),
      s3_bam = purrr::map_chr(.data$input, "tumor_wts_bam", .default = NA),
      # output
      s3_outdir_oncoanalyser = purrr::map_chr(.data$output, "output_directory", .default = NA)
    )
  d |>
    dplyr::select(
      meta_main_cols(),
      "SubjectID",
      "LibraryID",
      "SampleID",
      "s3_bam",
      "s3_outdir_oncoanalyser",
    )
}

meta_io_fromjson <- function(pmeta) {
  pmeta <- portal_meta_read(pmeta)
  pmeta |>
    dplyr::rowwise() |>
    dplyr::mutate(
      input = list(jsonlite::fromJSON(.data$input)),
      output = ifelse(!is.na(.data$output), list(jsonlite::fromJSON(.data$output)), NA)
    ) |>
    dplyr::ungroup()
}

meta_main_cols <- function() {
  c(
    "id", "wfr_name", "wfr_id", "version", "end_status", # "sequence_run", "batch_run",
    "start", "end", "portal_run_id"
  )
}

#' Read ICA Workflows Metadata via Portal API
#'
#' Reads ICA Workflows Metadata via Portal API using awscurl. See
#' https://github.com/okigan/awscurl for required `AWS_` environment variables.
#'
#' @param rows Number of rows to return.
#' @param params String containing additional params to pass to the `/workflows`
#' endpoint, e.g. `'&type_name=bclconvert'`.
#' @param pmeta Path to downloaded portal metadata file, or already parsed metadata tibble.
#' @param account UMCCR portal account (one of "prod", "dev", "stg").
#'
#' @return A tibble of the results from the given query.
#' @export
#'
#' @examples
#' \dontrun{
#' portal_meta_read(params = "&type_name=rnasum", rows = 4)
#' }
portal_meta_read <- function(pmeta = NULL, rows = 100, params = "", account = "prod") {
  assertthat::assert_that(account %in% c("prod", "dev", "stg"))
  au_tz <- "Australia/Melbourne"
  utc_tz <- "UTC"
  if (!is.null(pmeta)) {
    # if already parsed, just return it
    if (inherits(pmeta, "data.frame")) {
      assertthat::assert_that(
        all(c(
          "wfr_name", "type_name", "start", "end", "input", "output",
          "portal_run_id", "wfr_id", "wfl_id", "wfv_id", "version", "end_status"
        ) %in% colnames(pmeta))
      )
      return(pmeta)
    } else if (file.exists(pmeta)) {
      # local file downloaded via Portal
      # keep all character except start/end
      # and change to AEST timezone
      date_fmt_z <- "%Y-%m-%dT%H:%M:%SZ"
      ctypes <- readr::cols(
        .default = "c",
        start = readr::col_datetime(format = date_fmt_z),
        end = readr::col_datetime(format = date_fmt_z),
      )
      res <- pmeta |>
        readr::read_csv(col_types = ctypes) |>
        dplyr::mutate(
          start = lubridate::with_tz(.data$start, tz = au_tz),
          end = lubridate::with_tz(.data$end, tz = au_tz)
        )
      return(res)
    } else {
      stop("pmeta should be an already parsed dataframe or a local file.")
    }
  } # else pmeta is NULL, so read via portal API

  base_url <- glue("https://api.portal.{account}.umccr.org/iam")
  url1 <- utils::URLencode(glue("{base_url}/workflows?rowsPerPage={rows}{params}"))
  awscurl_cmd <- glue(
    "awscurl '{url1}' ",
    "--header 'Accept: application/json'"
  )
  message(glue("Running {awscurl_cmd}"))
  j <- system(awscurl_cmd, intern = TRUE)
  date_fmt <- "%Y-%m-%dT%H:%M:%S"
  d <- j |>
    jsonlite::fromJSON() |>
    purrr::pluck("results") |>
    tibble::as_tibble()
  d |>
    dplyr::mutate(
      start = as.POSIXct(.data$start, tz = utc_tz, format = date_fmt),
      end = as.POSIXct(.data$end, tz = utc_tz, format = date_fmt),
      start = lubridate::with_tz(.data$start, tz = au_tz),
      end = lubridate::with_tz(.data$end, tz = au_tz)
    )
}

#' Read ICA Workflows Metadata via Athena
#'
#' Reads the ICA Workflows Metadata for the given workflow run IDs.
#'
#' @param wfrids Character vector of wfr IDs to query.
#'
#' @return Tibble with a row per wfr ID.
#'
#' @examples
#' \dontrun{
#' wfrids <- c("wfr.1e764ca00e7a43a69e2424f250a34868")
#' portal_meta_read_athena(wfrids)
#' }
#' @export
portal_meta_read_athena <- function(wfrids = NULL) {
  assertthat::assert_that(!is.null(wfrids), all(grepl("^wfr\\.", wfrids)))
  RAthena::RAthena_options(clear_s3_resource = FALSE)
  con <- DBI::dbConnect(
    RAthena::athena(),
    work_group = "data_portal",
    rstudio_conn_tab = FALSE
  )
  wfrids_quote <- paste(shQuote(wfrids), collapse = ", ")
  q1 <- glue(
    'SELECT * FROM "data_portal"."data_portal"."data_portal_workflow" where wfr_id in ({wfrids_quote});'
  )
  RAthena::dbGetQuery(con, q1) |>
    tibble::as_tibble()
}


portal_meta_read_athena_tmp <- function(x = NULL) {
  assertthat::assert_that(!is.null(x))
  RAthena::RAthena_options(clear_s3_resource = FALSE)
  con <- DBI::dbConnect(
    RAthena::athena(),
    work_group = "data_portal",
    rstudio_conn_tab = FALSE
  )
  q1 <- glue(
    'SELECT * FROM "data_portal"."data_portal"."data_portal_workflow" {x}'
  )
  RAthena::dbGetQuery(con, q1) |>
    tibble::as_tibble()
}

#' Read Portal Limsrow Table
#'
#' Reads rows from the data portal's limsrow table, given a set of `LibraryID`s
#' to query.
#'
#' @param libids Character vector of LibraryID values to query for.
#'
#' @return Tibble with all rows from the data portal limsrow table where
#' there are hits with the `library_id` column.
#' @export
glims_portal_read <- function(libids) {
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

#' Read Google LIMS
#'
#' Reads UMCCR's Google LIMS spreadsheet.
#'
#' @return Tibble with all columns and rows from the Google LIMS spreadsheet.
#' @export
glims_excel_read <- function() {
  lims_key <- googledrive::drive_find("^Google LIMS$", shared_drive = "LIMS")$id
  lims <- lims_key |>
    googlesheets4::read_sheet("Sheet1", na = c(".", "", "-"), col_types = "c")
  lims |> readr::type_convert(col_types = readr::cols(.default = "c", Timestamp = "T"))
}
