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
    all(c("regex", "fun") %in% colnames(regexes))
  )
  # d[grepl(FILE_REGEX[["regex"]][i], d[["path"]]), "type"] <- FILE_REGEX[["name"]][i]
  for (i in seq_len(nrow(regexes))) {
    if (grepl(regexes[["regex"]][i], x)) {
      return(regexes[["fun"]][i])
    }
  }
  return(NA_character_)
}

FILE_REGEX <- tibble::tribble(
  ~regex, ~fun,
  "AlignCollapseFusionCaller_metrics\\.json\\.gz$", "TsoAlignCollapseFusionCallerMetricsFile",
  "TargetRegionCoverage\\.json\\.gz$", "TsoTargetRegionCoverageFile",
  "fragment_length_hist\\.json\\.gz$", "TsoFragmentLengthHistFile",
  "msi\\.json\\.gz$", "TsoMsiFile",
  "tmb\\.json\\.gz$", "TsoTmbFile",
  "TMB_Trace\\.tsv$", "TsoTmbTraceTsvFile",
  "_Fusions\\.csv$", "TsoFusionsCsvFile",
  "SampleAnalysisResults\\.json\\.gz$", "TsoSampleAnalysisResultsFile",
  "MergedSmallVariants\\.vcf\\.gz$", "TsoMergedSmallVariantsVcfFile",
  "fastqc_metrics\\.csv$", "NULL",
  "insert-stats\\.tab$", "NULL",
  "sv_metrics\\.csv$", "NULL",
  "trimmer_metrics\\.csv$", "NULL",
  "wgs_contig_mean_cov(_tumor|_normal|)\\.csv$", "ContigMeanCovFile",
  "wgs_fine_hist(_tumor|_normal|)\\.csv$", "FineHistFile",
  "wgs_hist(_tumor|_normal|)\\.csv$", "NULL",
  "wgs_overall_mean_cov(_tumor|_normal|)\\.csv$", "NULL",
  "wgs_coverage_metrics(_tumor|_normal|)\\.csv$", "CoverageMetricsFile",
  "fragment_length_hist\\.csv$", "FragmentLengthHistFile",
  "mapping_metrics\\.csv$", "MappingMetricsFile",
  "ploidy_estimation_metrics\\.csv$", "PloidyEstimationMetricsFile",
  "replay\\.json$", "ReplayFile",
  "time_metrics\\.csv$", "TimeMetricsFile",
  "vc_metrics\\.csv$", "VCMetricsFile",
  "multiqc_data\\.json", "MultiqcFile",
  "somatic\\.pcgr\\.json\\.gz$", "PcgrJsonFile",
  "somatic\\.pcgr\\.snvs_indels\\.tiers\\.tsv$", "PcgrTiersFile",
  "chord\\.tsv\\.gz$", "UmChordTsvFile",
  "hrdetect\\.tsv\\.gz$", "UmHrdetectTsvFile",
  "snv_2015\\.tsv\\.gz$", "UmSigsSnv2015File",
  "snv_2020\\.tsv\\.gz$", "UmSigsSnv2020File",
  "-qc_summary\\.tsv\\.gz$", "UmQcSumFile"
)

func_selector <- function(type, tbl = FILE_REGEX) {
  assertthat::assert_that(
    inherits(tbl, "data.frame"),
    "fun" %in% colnames(tbl)
  )
  l <- tbl[["fun"]] |>
    purrr::set_names(tbl[["fun"]])
  if (!type %in% names(l)) {
    return(NULL)
  }
  # need this workaround to evaluate string as function
  eval(parse(text = l[[type]]))
}
