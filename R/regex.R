#' Match File to Regex
#'
#' Matches a given file with the regexes found in FILE_REGEX and if there is
#' a match, it returns the 'name' of that match.
#'
#' @param x File to match.
#'
#' @return The 'name' of the matching regex from FILE_REGEX, or NA if there is
#' no match made.
#'
#' @examples
#' match_regex("foo.msi.json.gz")
#' match_regex("foo.fake.tsv")
#'
#' @export
match_regex <- function(x) {
  # d[grepl(FILE_REGEX[["regex"]][i], d[["path"]]), "type"] <- FILE_REGEX[["name"]][i]
  for (i in seq_len(nrow(FILE_REGEX))) {
    if (grepl(FILE_REGEX[["regex"]][i], x)) {
      return(FILE_REGEX[["name"]][i])
    }
  }
  return(NA_character_)
}

FILE_REGEX <- tibble::tribble(
  ~regex, ~name,
  "MetricsOutput\\.tsv$", "tso/metrics_output",
  "AlignCollapseFusionCaller_metrics\\.json\\.gz$", "tso/align_collapse_fusion_caller_metrics",
  "TargetRegionCoverage\\.json\\.gz$", "tso/target_region_coverage",
  "fragment_length_hist\\.json\\.gz$", "tso/fragment_length_hist",
  "msi\\.json\\.gz$", "tso/msi",
  "tmb\\.json\\.gz$", "tso/tmb",
  "TMB_Trace\\.tsv$", "tso/tmb_trace_tsv",
  "_CombinedVariantOutput\\.tsv$", "tso/combined_variant_output",
  "_Fusions\\.json\\.gz$", "tso/fusions_json",
  "_Fusions\\.csv$", "tso/fusions_csv",
  "SampleAnalysisResults\\.json\\.gz$", "tso/sample_analysis_results",
  "wgs_contig_mean_cov_(tumor|normal)\\.csv$", "dragen/contig_mean_cov",
  "wgs_fine_hist_(tumor|normal)\\.csv$", "dragen/fine_hist",
  "wgs_coverage_metrics_(tumor|normal)\\.csv$", "dragen/coverage_metrics",
  "fragment_length_hist\\.csv$", "dragen/fragment_length_hist",
  "mapping_metrics\\.csv$", "dragen/mapping_metrics",
  "ploidy_estimation_metrics\\.csv$", "dragen/ploidy_estimation_metrics",
  "replay\\.json$", "dragen/replay",
  "time_metrics\\.csv$", "dragen/time_metrics",
  "vc_metrics\\.csv$", "dragen/vc_metrics",
  "multiqc_data\\.json", "multiqc",
)
