#' Match File to Regex
#'
#' Matches a given file with the regexes found in `DR_FILE_REGEX` and if there is
#' a match, it returns the 'name' of that match.
#'
#' @param x File to match.
#' @param regexes Tibble with `regex` and `fun`ction name.
#'
#' @return The function corresponding to the matching regex from `DR_FILE_REGEX`, or
#' NA if there is no match made.
#'
#' @examples
#' match_regex("foo.msi.json.gz")
#' match_regex("foo.fake.tsv")
#'
#' @export
match_regex <- function(x, regexes = DR_FILE_REGEX) {
  assertthat::assert_that(
    inherits(regexes, "data.frame"),
    all(c("regex", "fun") %in% colnames(regexes))
  )
  # d[grepl(DR_FILE_REGEX[["regex"]][i], d[["path"]]), "type"] <- DR_FILE_REGEX[["name"]][i]
  for (i in seq_len(nrow(regexes))) {
    if (grepl(regexes[["regex"]][i], x)) {
      return(regexes[["fun"]][i])
    }
  }
  return(NA_character_)
}

DR_FILE_REGEX <- tibble::tribble(
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
  "MergedSmallVariants\\.vcf\\.gz\\.tbi$", "TsoMergedSmallVariantsVcfIndexFile",
  # "MergedSmallVariants\\.genome\\.vcf\\.gz$", "TsoMergedSmallVariantsGenomeVcfFile",
  # "MergedSmallVariants\\.genome\\.vcf\\.gz\\.tbi$", "TsoMergedSmallVariantsGenomeVcfIndexFile",
  "CopyNumberVariants\\.vcf\\.gz$", "TsoCopyNumberVariantsVcfFile",
  "CopyNumberVariants\\.vcf\\.gz\\.tbi$", "TsoCopyNumberVariantsVcfIndexFile",
  "CombinedVariantOutput\\.tsv$", "TsoCombinedVariantOutputFile",
  "fastqc_metrics\\.csv$", "FastqcMetricsFile",
  "sv_metrics\\.csv$", "SvMetricsFile",
  "trimmer_metrics\\.csv$", "TrimmerMetricsFile",
  "wgs_contig_mean_cov(_tumor|_normal|)\\.csv$", "WgsContigMeanCovFile",
  "wgs_fine_hist(_tumor|_normal|)\\.csv$", "WgsFineHistFile",
  "wgs_hist(_tumor|_normal|)\\.csv$", "WgsHistFile",
  "wgs_coverage_metrics(_tumor|_normal|)\\.csv$", "WgsCoverageMetricsFile",
  "fragment_length_hist\\.csv$", "FragmentLengthHistFile",
  "mapping_metrics\\.csv$", "MappingMetricsFile",
  "ploidy_estimation_metrics\\.csv$", "PloidyEstimationMetricsFile",
  "replay\\.json$", "ReplayFile",
  "time_metrics\\.csv$", "TimeMetricsFile",
  "vc_metrics\\.csv$", "VCMetricsFile",
  "multiqc_data\\.json", "MultiqcFile",
  "somatic\\.pcgr\\.json\\.gz$", "PcgrJsonFile",
  "somatic\\.pcgr\\.snvs_indels\\.tiers\\.tsv$", "PcgrTiersFile",
  # "chord\\.tsv\\.gz$", "UmChordTsvFile",
  # "hrdetect\\.tsv\\.gz$", "UmHrdetectTsvFile",
  # "snv_2015\\.tsv\\.gz$", "UmSigsSnvFile",
  # "snv_2020\\.tsv\\.gz$", "UmSigsSnvFile",
  # "-qc_summary\\.tsv\\.gz$", "UmQcSumFile",
  "bcftools_stats\\.txt$", "BcftoolsStatsFile"
)

FILES_DOWNLOAD_BUT_IGNORE <- c(
  "TsoMergedSmallVariantsVcfIndexFile",
  "TsoCopyNumberVariantsVcfIndexFile",
  "TsoMergedSmallVariantsGenomeVcfIndexFile"
)

#' Evaluate dracarys Function
#'
#' This is somewhat a hack for getting a function to evaluate based on a lookup
#' vector. If the function is not found, it returns NULL.
#'
#' @param f Name of function to evaluate.
#' @param v Character vector of strings evaluating to functions. By default,
#' this points to the functions in the DR_FILE_REGEX dracarys tibble.
#'
#' @return Evaluated function.
#' @examples
#' mean_1_to_10 <- dr_func_eval("mean", v = c("mean", "sd"))(1:10)
#' x <- system.file("extdata/tso/sample705.fragment_length_hist.json.gz", package = "dracarys")
#' obj <- dr_func_eval("TsoFragmentLengthHistFile")$new(x)
#' @testexamples
#' expect_equal("sample705.fragment_length_hist.json.gz", obj$bname())
#' expect_equal(mean_1_to_10, base::mean(1:10))
#' expect_null(dr_func_eval("foo"))
#' @export
dr_func_eval <- function(f, v = NULL) {
  v <- v %||% DR_FILE_REGEX[["fun"]]
  if (!f %in% v) {
    return(NULL)
  }
  # evaluate string
  eval(parse(text = f))
}

#' Get dracarys `DR_FILE_REGEX`
#'
#' @return `DR_FILE_REGEX` R tibble object.
#' @export
file_regex_getter <- function() {
  DR_FILE_REGEX
}
