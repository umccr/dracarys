#' Match File to Regex
#'
#' Matches a given file with the regexes found in FILE_REGEX and if there is
#' a match, it returns the 'name' of that match.
#'
#' @param x File to match.
#' @param regexes Tibble with regex and name.
#'
#' @return The 'name' of the matching regex from FILE_REGEX, or NA if there is
#' no match made.
#'
#' @examples
#' match_regex("foo.msi.json.gz")
#' match_regex("foo.fake.tsv")
#'
#' @export
match_regex <- function(x, regexes = FILE_REGEX) {
  assertthat::assert_that(
    inherits(regexes, "data.frame"),
    all(c("regex", "name") %in% colnames(regexes))
  )
  # d[grepl(FILE_REGEX[["regex"]][i], d[["path"]]), "type"] <- FILE_REGEX[["name"]][i]
  for (i in seq_len(nrow(regexes))) {
    if (grepl(regexes[["regex"]][i], x)) {
      return(regexes[["name"]][i])
    }
  }
  return(NA_character_)
}

FILE_REGEX <- tibble::tribble(
  ~regex, ~name, ~fun,
  "AlignCollapseFusionCaller_metrics\\.json\\.gz$", "tso__align_collapse_fusion_caller_metrics", "TsoAlignCollapseFusionCallerMetricsFile",
  "TargetRegionCoverage\\.json\\.gz$", "tso__target_region_coverage", "TsoTargetRegionCoverageFile",
  "fragment_length_hist\\.json\\.gz$", "tso__fragment_length_hist", "TsoFragmentLengthHistFile",
  "msi\\.json\\.gz$", "tso__msi", "TsoMsiFile",
  "tmb\\.json\\.gz$", "tso__tmb", "TsoTmbFile",
  "TMB_Trace\\.tsv$", "tso__tmb_trace_tsv", "TsoTmbTraceTsvFile",
  "_Fusions\\.csv$", "tso__fusions_csv", "TsoFusionsCsvFile",
  "SampleAnalysisResults\\.json\\.gz$", "tso__sample_analysis_results", "TsoSampleAnalysisResultsFile",
  "MergedSmallVariants\\.vcf\\.gz$", "tso__mergedsmallvariants_vcf", "TsoMergedSmallVariantsVcfFile",
  "fastqc_metrics\\.csv$", "dragen/fastqc_metrics", NULL,
  "insert-stats\\.tab$", "dragen/insert_stats", NULL,
  "sv_metrics\\.csv$", "dragen/sv_metrics", NULL,
  "trimmer_metrics\\.csv$", "dragen/trimmer_metrics", NULL,
  "wgs_contig_mean_cov(_tumor|_normal|)\\.csv$", "dragen/contig_mean_cov", NULL,
  "wgs_fine_hist(_tumor|_normal|)\\.csv$", "dragen/fine_hist", NULL,
  "wgs_hist(_tumor|_normal|)\\.csv$", "dragen/hist", NULL,
  "wgs_overall_mean_cov(_tumor|_normal|)\\.csv$", "dragen/coverage_overall", NULL,
  "wgs_coverage_metrics(_tumor|_normal|)\\.csv$", "dragen/coverage_metrics", NULL,
  "fragment_length_hist\\.csv$", "dragen/fragment_length_hist", NULL,
  "mapping_metrics\\.csv$", "dragen/mapping_metrics", NULL,
  "ploidy_estimation_metrics\\.csv$", "dragen/ploidy_estimation_metrics", NULL,
  "replay\\.json$", "dragen/replay", NULL,
  "time_metrics\\.csv$", "dragen/time_metrics", NULL,
  "vc_metrics\\.csv$", "dragen/vc_metrics", NULL,
  "multiqc_data\\.json", "multiqc", "MultiqcFile",
  "somatic\\.pcgr\\.json\\.gz$", "pcgr__json", "PcgrJsonFile",
  "somatic\\.pcgr\\.snvs_indels\\.tiers\\.tsv$", "pcgr__tiers", "PcgrTiersFile",
  "chord\\.tsv\\.gz$", "um__chord", "UmChordTsvFile",
  "hrdetect\\.tsv\\.gz$", "um__hrdetect", "UmHrdetectTsvFile",
  "snv_2015\\.tsv\\.gz$", "um__sigs2015_snv", "UmSigsSnvFile",
  "snv_2020\\.tsv\\.gz$", "um__sigs2020_snv", "UmSigsSnvFile",
  "-qc_summary\\.tsv\\.gz$", "um__qcsum", "UmQcSumFile"
)

func_selector <- function(type, tbl = FILE_REGEX) {
  assertthat::assert_that(
    inherits(tbl, "data.frame"),
    all(c("fun", "name") %in% colnames(tbl))
  )
  l <- tbl[["fun"]] |>
    purrr::set_names(tbl[["name"]])
  if (!type %in% names(l)) {
    return(NULL)
  }
  # need this workaround to evaluate string as function
  eval(parse(text = l[[type]]))
}
