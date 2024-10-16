#' Read DRAGEN CNV Metrics
#'
#' Reads the `cnv_metrics.csv` file output from DRAGEN.
#'
#' @param x Path to file.
#'
#' @return Tibble with metrics.
#' @export
dragen_cnv_metrics_read <- function(x) {
  d0 <- readr::read_lines(x)
  assertthat::assert_that(grepl("SEX GENOTYPER", d0[1]))
  abbrev_nm <- c(
    "Bases in reference genome" = "bases_in_ref_genome_dragen",
    "Average alignment coverage over genome" = "cov_alignment_avg_over_genome_dragen",
    "Number of alignment records" = "n_alignment_records",
    "Number of filtered records (total)" = "n_filtered_records_tot",
    "Number of filtered records (duplicates)" = "n_filtered_records_dup",
    "Number of filtered records (MAPQ)" = "n_filtered_records_mapq",
    "Number of filtered records (unmapped)" = "n_filtered_records_unmap",
    "Coverage MAD" = "coverage_mad",
    "Gene Scaled MAD" = "gene_scaled_mad",
    "Median Bin Count" = "median_bin_count",
    "Number of target intervals" = "n_target_intervals",
    "Number of normal samples" = "n_samp_norm",
    "Number of segments" = "n_seg",
    "Number of amplifications" = "n_amp",
    "Number of deletions" = "n_del",
    "Number of passing amplifications" = "n_amp_pass",
    "Number of passing deletions" = "n_del_pass"
  )

  d1 <- d0 |>
    tibble::as_tibble_col(column_name = "value") |>
    tidyr::separate_wider_delim(
      "value",
      names = c("category", "extra", "var", "count", "pct"),
      delim = ",", too_few = "align_start"
    )
  sexgt <- d1 |>
    dplyr::filter(.data$category == "SEX GENOTYPER") |>
    dplyr::select(sexgt = "count", sexgt_pct = "pct")

  d2 <- d1 |>
    dplyr::filter(!.data$category == "SEX GENOTYPER") |>
    dplyr::mutate(
      count = as.numeric(.data$count),
      pct = round(as.numeric(.data$pct), 2),
      var = dplyr::recode(.data$var, !!!abbrev_nm)
    ) |>
    dplyr::select("var", "count", "pct") |>
    tidyr::pivot_longer(c("count", "pct")) |>
    dplyr::filter(!is.na(.data$value)) |>
    dplyr::mutate(
      name = dplyr::if_else(.data$name == "count", "", "_pct"),
      var = glue("{.data$var}{.data$name}")
    ) |>
    dplyr::select("var", "value") |>
    tidyr::pivot_wider(names_from = "var", values_from = "value")
  res <- dplyr::bind_cols(sexgt, d2)
  return(res)
}

#' Read DRAGEN SV Metrics
#'
#' Reads the `sv_metrics.csv` file output from DRAGEN.
#'
#' @param x Path to file.
#'
#' @return Tibble with metrics.
#' @export
dragen_sv_metrics_read <- function(x) {
  d <- readr::read_lines(x)
  assertthat::assert_that(grepl("SV SUMMARY", d[1]))
  abbrev_nm <- c(
    "Number of deletions (PASS)" = "del",
    "Number of insertions (PASS)" = "ins",
    "Number of duplications (PASS)" = "dup",
    "Number of breakend pairs (PASS)" = "bnd"
  )
  d |>
    tibble::as_tibble_col(column_name = "value") |>
    dplyr::filter(!grepl("Total number of structural variants", .data$value)) |>
    tidyr::separate_wider_delim(
      "value",
      names = c("svsum", "sample", "var", "count", "pct"), delim = ",",
      too_few = "align_start"
    ) |>
    dplyr::mutate(
      count = as.numeric(.data$count),
      var = dplyr::recode(.data$var, !!!abbrev_nm)
    ) |>
    dplyr::select("var", "count") |>
    tidyr::pivot_wider(names_from = "var", values_from = "count")
}

#' Read DRAGEN Trimmer Metrics
#'
#' Reads the `trimmer_metrics.csv` file output from DRAGEN.
#'
#' @param x Path to file.
#'
#' @return Tibble with metrics.
#' @export
dragen_trimmer_metrics_read <- function(x) {
  d <- readr::read_lines(x)
  assertthat::assert_that(grepl("TRIMMER STATISTICS", d[1]))
  abbrev_nm <- c(
    "Total input reads"                              = "reads_tot_input_dragen",
    "Total input bases"                              = "bases_tot_dragen",
    "Total input bases R1"                           = "bases_r1_dragen",
    "Total input bases R2"                           = "bases_r2_dragen",
    "Average input read length"                      = "read_len_avg_dragen",
    "Total trimmed reads"                            = "reads_trimmed_tot_dragen",
    "Total trimmed bases"                            = "bases_trimmed_tot_dragen",
    "Average bases trimmed per read"                 = "bases_trimmed_avg_per_read_dragen",
    "Average bases trimmed per trimmed read"         = "bases_trimmed_avg_per_trimmedread_dragen",
    "Remaining poly-G K-mers R1 3prime"              = "polygkmers3r1_remaining_dragen",
    "Remaining poly-G K-mers R2 3prime"              = "polygkmers3r2_remaining_dragen",
    "Poly-G soft trimmed reads unfiltered R1 3prime" = "polyg_soft_trimmed_reads_unfilt_3r1_dragen",
    "Poly-G soft trimmed reads unfiltered R2 3prime" = "polyg_soft_trimmed_reads_unfilt_3r2_dragen",
    "Poly-G soft trimmed reads filtered R1 3prime"   = "polyg_soft_trimmed_reads_filt_3r1_dragen",
    "Poly-G soft trimmed reads filtered R2 3prime"   = "polyg_soft_trimmed_reads_filt_3r2_dragen",
    "Poly-G soft trimmed bases unfiltered R1 3prime" = "polyg_soft_trimmed_bases_unfilt_3r1_dragen",
    "Poly-G soft trimmed bases unfiltered R2 3prime" = "polyg_soft_trimmed_bases_unfilt_3r2_dragen",
    "Poly-G soft trimmed bases filtered R1 3prime"   = "polyg_soft_trimmed_bases_filt_3r1_dragen",
    "Poly-G soft trimmed bases filtered R2 3prime"   = "polyg_soft_trimmed_bases_filt_3r2_dragen",
    "Total filtered reads"                           = "reads_tot_filt_dragen",
    "Reads filtered for minimum read length R1"      = "reads_filt_minreadlenr1_dragen",
    "Reads filtered for minimum read length R2"      = "reads_filt_minreadlenr2_dragen"
  )

  d |>
    tibble::as_tibble_col(column_name = "value") |>
    tidyr::separate_wider_delim("value", names = c("category", "extra", "var", "count", "pct"), delim = ",", too_few = "align_start") |>
    dplyr::mutate(
      count = as.numeric(.data$count),
      pct = round(as.numeric(.data$pct), 2),
      var = dplyr::recode(.data$var, !!!abbrev_nm)
    ) |>
    dplyr::select("var", "count", "pct") |>
    tidyr::pivot_longer(c("count", "pct")) |>
    dplyr::filter(!is.na(.data$value)) |>
    dplyr::mutate(
      name = dplyr::if_else(.data$name == "count", "", "_pct"),
      var = glue("{.data$var}{.data$name}")
    ) |>
    dplyr::select("var", "value") |>
    tidyr::pivot_wider(names_from = "var", values_from = "value")
}

#' Read DRAGEN VariantCall Metrics
#'
#' Reads the `vc_metrics.csv` file output from DRAGEN.
#'
#' @param x Path to file.
#'
#' @return Tibble with metrics.
#' @export
dragen_vc_metrics_read <- function(x) {
  abbrev_nm1 <- tibble::tribble(
    ~raw, ~clean, ~region,
    "Total", "var_tot_dragen", FALSE,
    "Biallelic", "var_biallelic_dragen", FALSE,
    "Multiallelic", "var_multiallelic_dragen", FALSE,
    "SNPs", "var_snp_dragen", FALSE,
    "Insertions (Hom)", "var_ins_hom_dragen", FALSE,
    "Insertions (Het)", "var_ins_het_dragen", FALSE,
    "Deletions (Hom)", "var_del_hom_dragen", FALSE,
    "Deletions (Het)", "var_del_het_dragen", FALSE,
    "Indels (Het)", "var_indel_het_dragen", FALSE,
    "Chr X number of SNPs over ", "var_snp_x_over_", TRUE,
    "Chr Y number of SNPs over ", "var_snp_y_over_", TRUE,
    "(Chr X SNPs)/(chr Y SNPs) ratio over ", "var_x_over_y_snp_ratio_over_", TRUE,
    "SNP Transitions", "var_snp_transitions_dragen", FALSE,
    "SNP Transversions", "var_snp_transversions_dragen", FALSE,
    "Ti/Tv ratio", "var_ti_tv_ratio_dragen", FALSE,
    "Heterozygous", "var_heterozygous_dragen", FALSE,
    "Homozygous", "var_homozygous_dragen", FALSE,
    "Het/Hom ratio", "var_het_hom_ratio_dragen", FALSE,
    "In dbSNP", "var_in_dbsnp_dragen", FALSE,
    "Not in dbSNP", "var_nin_dbsnp_dragen", FALSE,
    "Percent Callability", "callability_pct_dragen", FALSE,
    "Percent Autosome Callability", "callability_auto_pct_dragen", FALSE,
    "Number of samples", "sample_num_dragen", FALSE,
    "Reads Processed", "reads_processed_dragen", FALSE,
    "Child Sample", "sample_child_dragen", FALSE
  )
  raw <- readr::read_lines(x)
  assertthat::assert_that(grepl("VARIANT CALLER", raw[1]))
  # first detect if this is genome or target region
  res <- raw |>
    tibble::as_tibble_col(column_name = "value") |>
    tidyr::separate_wider_delim(
      "value",
      names = c("category", "sample", "var", "count", "pct"),
      delim = ",", too_few = "align_start"
    )
  reg1 <- NULL
  str1 <- NULL
  tmp <- res |>
    dplyr::filter(grepl("Chr X number of SNPs over ", .data$var)) |>
    dplyr::slice_head(n = 1) |>
    dplyr::pull("var")
  assertthat::assert_that(length(tmp) == 1)
  if (grepl("genome", tmp)) {
    str1 <- "genome"
    reg1 <- "genome"
  } else if (grepl("QC coverage region", tmp)) {
    str1 <- "QC coverage region"
    reg1 <- "qccovreg"
  } else if (grepl("target region", tmp)) {
    str1 <- "target region"
    reg1 <- "targetreg"
  } else {
    cli::cli_abort("Cannot determine the varcall region from: {x}")
  }
  abbrev_nm <- abbrev_nm1 |>
    dplyr::mutate(
      raw = ifelse(.data$region, glue("{.data$raw}{str1}"), .data$raw),
      clean = ifelse(.data$region, glue("{.data$clean}{reg1}_dragen"), .data$clean)
    ) |>
    dplyr::select("raw", "clean") |>
    tibble::deframe()

  d <- res |>
    dplyr::mutate(
      var = dplyr::recode(.data$var, !!!abbrev_nm),
      count = dplyr::na_if(.data$count, "NA"),
      count = as.numeric(.data$count),
      pct = round(as.numeric(.data$pct), 2),
      category = dplyr::case_when(
        grepl("SUMMARY", .data$category) ~ "summary",
        grepl("PREFILTER", .data$category) ~ "prefilter",
        grepl("POSTFILTER", .data$category) ~ "postfilter",
        TRUE ~ "unknown"
      )
    ) |>
    dplyr::filter(.data$category != "summary") |>
    dplyr::select("category", "sample", "var", "count", "pct")
  # pivot
  d |>
    tidyr::pivot_longer(c("count", "pct")) |>
    dplyr::mutate(
      name = dplyr::if_else(.data$name == "count", "", "_pct"),
      var = glue("{.data$var}{.data$name}")
    ) |>
    dplyr::select("category", "sample", "var", "value") |>
    dplyr::filter(!is.na(.data$value)) |>
    tidyr::pivot_wider(names_from = "var", values_from = "value")
}

#' Read DRAGEN Mapping Metrics
#'
#' Reads the `mapping_metrics.csv` file output from DRAGEN.
#'
#' @param x Path to file.
#'
#' @return Tibble with metrics.
#' @export
dragen_mapping_metrics_read <- function(x) {
  abbrev_nm <- c(
    "Total input reads" = "reads_tot_input_dragen",
    "Number of duplicate marked reads" = "reads_num_dupmarked_dragen",
    "Number of duplicate marked and mate reads removed" = "reads_num_dupmarked_mate_reads_removed_dragen",
    "Number of unique reads (excl. duplicate marked reads)" = "reads_num_uniq_dragen",
    "Reads with mate sequenced" = "reads_w_mate_seq_dragen",
    "Reads without mate sequenced" = "reads_wo_mate_seq_dragen",
    "QC-failed reads" = "reads_qcfail_dragen",
    "Mapped reads" = "reads_mapped_dragen",
    "Mapped reads adjusted for filtered mapping" = "reads_mapped_adjfilt_dragen",
    "Mapped reads R1" = "reads_mapped_r1_dragen",
    "Mapped reads R2" = "reads_mapped_r2_dragen",
    "Number of unique & mapped reads (excl. duplicate marked reads)" = "reads_num_uniq_mapped_dragen",
    "Unmapped reads" = "reads_unmapped_dragen",
    "Unmapped reads adjusted for filtered mapping" = "reads_unmapped_adjfilt_dragen",
    "Adjustment of reads matching non-reference decoys" = "reads_match_nonref_decoys_adj_dragen",
    "Singleton reads (itself mapped; mate unmapped)" = "reads_singleton_dragen",
    "Paired reads (itself & mate mapped)" = "reads_paired_dragen",
    "Properly paired reads" = "reads_paired_proper_dragen",
    "Not properly paired reads (discordant)" = "reads_discordant_dragen",
    "Paired reads mapped to different chromosomes" = "reads_paired_mapped_diff_chrom_dragen",
    "Paired reads mapped to different chromosomes (MAPQ>=10)" = "reads_paired_mapped_diff_chrom_mapq10_dragen",
    "Reads with MAPQ [40:inf)" = "reads_mapq_40_inf_dragen",
    "Reads with MAPQ [30:40)" = "reads_mapq_30_40_dragen",
    "Reads with MAPQ [20:30)" = "reads_mapq_20_30_dragen",
    "Reads with MAPQ [10:20)" = "reads_mapq_10_20_dragen",
    "Reads with MAPQ [ 0:10)" = "reads_mapq_0_10_dragen",
    "Reads with MAPQ NA (Unmapped reads)" = "reads_mapq_na_unmapped_dragen",
    "Reads with indel R1" = "reads_indel_r1_dragen",
    "Reads with indel R2" = "reads_indel_r2_dragen",
    "Total bases" = "bases_tot_dragen",
    "Total bases R1" = "bases_tot_r1_dragen",
    "Total bases R2" = "bases_tot_r2_dragen",
    "Mapped bases" = "bases_mapped_dragen",
    "Mapped bases R1" = "bases_mapped_r1_dragen",
    "Mapped bases R2" = "bases_mapped_r2_dragen",
    "Soft-clipped bases" = "bases_softclip_dragen",
    "Soft-clipped bases R1" = "bases_softclip_r1_dragen",
    "Soft-clipped bases R2" = "bases_softclip_r2_dragen",
    "Hard-clipped bases" = "bases_hardclip_dragen",
    "Hard-clipped bases R1" = "bases_hardclip_r1_dragen",
    "Hard-clipped bases R2" = "bases_hardclip_r2_dragen",
    "Mismatched bases R1" = "bases_mismatched_r1_dragen",
    "Mismatched bases R2" = "bases_mismatched_r2_dragen",
    "Mismatched bases R1 (excl. indels)" = "bases_mismatched_r1_noindels_dragen",
    "Mismatched bases R2 (excl. indels)" = "bases_mismatched_r2_noindels_dragen",
    "Q30 bases" = "bases_q30_dragen",
    "Q30 bases R1" = "bases_q30_r1_dragen",
    "Q30 bases R2" = "bases_q30_r2_dragen",
    "Q30 bases (excl. dups & clipped bases)" = "bases_q30_nodups_noclipped_dragen",
    "Total alignments" = "alignments_tot_dragen",
    "Secondary alignments" = "alignments_secondary_dragen",
    "Supplementary (chimeric) alignments" = "alignments_chimeric_dragen",
    "Estimated read length" = "read_len_dragen",
    "Bases in reference genome" = "bases_in_ref_genome_dragen",
    "Bases in target bed [% of genome]" = "bases_in_target_bed_genome_pct_dragen",
    "Insert length: mean" = "insert_len_mean_dragen",
    "Insert length: median" = "insert_len_median_dragen",
    "Insert length: standard deviation" = "insert_len_std_dev_dragen",
    "Provided sex chromosome ploidy" = "ploidy_sex_chrom_provided_dragen",
    "Estimated sample contamination" = "contamination_est_dragen",
    "Estimated sample contamination standard error" = "contamination_stderr_est_dragen",
    "DRAGEN mapping rate [mil. reads/second]" = "mapping_rate_dragen_milreads_per_sec_dragen",
    "Total reads in RG" = "reads_tot_rg_dragen",
    "Mapped reads adjusted for excluded mapping" = "reads_mapped_adjexcl_dragen",
    "Mapped reads adjusted for filtered and excluded mapping" = "reads_mapped_adjfiltexcl_dragen",
    "Unmapped reads adjusted for excluded mapping" = "reads_unmapped_adjexcl_dragen",
    "Unmapped reads adjusted for filtered and excluded mapping" = "reads_unmapped_adjfiltexcl_dragen",
    "Reads mapping to multiple locations" = "reads_map_multiloc_dragen",
    "Adjustment of reads matching filter contigs" = "reads_match_filt_contig_adj_dragen",
    "Reads with splice junction" = "reads_splicejunc_dragen",
    "Average sequenced coverage over genome" = "cov_avg_seq_over_genome_dragen",
    "Filtered rRNA reads" = "reads_rrna_filtered_dragen"
  )
  raw <- readr::read_lines(x)
  assertthat::assert_that(grepl("MAPPING/ALIGNING", raw[1]))
  # split by RG and non-RG
  # tidy
  d <- raw |>
    tibble::as_tibble_col(column_name = "value") |>
    tidyr::separate_wider_delim(
      "value",
      names = c("category", "RG", "var", "count", "pct"),
      delim = ",", too_few = "align_start"
    ) |>
    dplyr::mutate(
      count = dplyr::na_if(.data$count, "NA"),
      count = as.numeric(.data$count),
      pct = as.numeric(.data$pct),
      var = dplyr::recode(.data$var, !!!abbrev_nm),
      RG = dplyr::if_else(.data$RG == "", "Total", .data$RG)
    ) |>
    dplyr::select("RG", "var", "count", "pct")
  # pivot
  d |>
    tidyr::pivot_longer(c("count", "pct")) |>
    dplyr::mutate(
      name = dplyr::if_else(.data$name == "count", "", "_pct"),
      var = glue("{.data$var}{.data$name}")
    ) |>
    dplyr::select("RG", "var", "value") |>
    dplyr::filter(!is.na(.data$value)) |>
    tidyr::pivot_wider(names_from = "var", values_from = "value")
}

#' Read DRAGEN Coverage Metrics
#'
#' Reads the `coverage_metrics.csv` file generated by DRAGEN.
#' @param x Path to file.
#'
#' @return Tibble with metrics.
#' @export
dragen_coverage_metrics_read <- function(x) {
  # all rows except 'Aligned bases' and 'Aligned reads' refer to the region
  abbrev_nm <- tibble::tribble(
    ~raw, ~clean, ~region,
    "Aligned bases", "bases_aligned_tot_dragen", FALSE,
    "Aligned reads", "reads_aligned_tot_dragen", FALSE,
    "Aligned bases in ", "bases_aligned_", TRUE,
    "Average alignment coverage over ", "cov_alignment_avg_over_", TRUE,
    "Uniformity of coverage (PCT > 0.2*mean) over ", "cov_uniformity_pct_gt02mean_", TRUE,
    "Uniformity of coverage (PCT > 0.4*mean) over ", "cov_uniformity_pct_gt04mean_", TRUE,
    "Average chr X coverage over ", "cov_avg_x_over_", TRUE,
    "Average chr Y coverage over ", "cov_avg_y_over_", TRUE,
    "Average mitochondrial coverage over ", "cov_avg_mt_over_", TRUE,
    "Average autosomal coverage over ", "cov_avg_auto_over_", TRUE,
    "Median autosomal coverage over ", "cov_median_auto_over_", TRUE,
    "Mean/Median autosomal coverage ratio over ", "cov_mean_median_auto_ratio_over_", TRUE,
    "Aligned reads in ", "reads_aligned_in_", TRUE
  )
  raw <- readr::read_lines(x)
  assertthat::assert_that(grepl("COVERAGE SUMMARY", raw[1]))
  # first detect if this is genome, QC coverage region, or target region
  res <- raw |>
    tibble::as_tibble_col(column_name = "value") |>
    tidyr::separate_wider_delim(
      "value",
      delim = ",", too_few = "align_start",
      names = c("category", "dummy1", "var", "value", "pct")
    )
  reg1 <- NULL
  str1 <- NULL
  tmp <- res |>
    dplyr::filter(grepl("PCT of .* with coverage ", .data$var)) |>
    dplyr::slice_head(n = 1) |>
    dplyr::pull("var")
  assertthat::assert_that(length(tmp) == 1)
  if (grepl("genome", tmp)) {
    str1 <- "genome"
    reg1 <- "genome"
  } else if (grepl("QC coverage region", tmp)) {
    str1 <- "QC coverage region"
    reg1 <- "qccovreg"
  } else if (grepl("target region", tmp)) {
    str1 <- "target region"
    reg1 <- "targetreg"
  } else {
    cli::cli_abort("Cannot determine the coverage region from: {x}")
  }
  abbrev_nm <- abbrev_nm |>
    dplyr::mutate(
      raw = ifelse(.data$region, glue("{.data$raw}{str1}"), .data$raw),
      clean = ifelse(.data$region, glue("{.data$clean}{reg1}_dragen"), .data$clean)
    ) |>
    dplyr::select("raw", "clean") |>
    tibble::deframe()
  # split to rename the
  # "PCT of genome with coverage [100x: inf)" values
  pat <- glue("PCT of {str1} with coverage ")
  res1 <- res |>
    # pct just shows % for a couple rows which can be
    # calculated from their above values
    dplyr::filter(!grepl(pat, .data$var)) |>
    dplyr::select("var", "value")
  res2 <- res |>
    dplyr::filter(grepl(pat, .data$var)) |>
    dplyr::mutate(
      var = sub(pat, "", .data$var),
      var = gsub("\\[|\\]|\\(|\\)| ", "", .data$var),
      var = gsub("x", "", .data$var),
      var = gsub("inf", "Inf", .data$var)
    ) |>
    tidyr::separate_wider_delim("var", names = c("start", "end"), delim = ":") |>
    dplyr::mutate(var = as.character(glue("cov_pct_{start}_{end}_{reg1}_dragen"))) |>
    dplyr::select("var", "value")
  res <- dplyr::bind_rows(res1, res2) |>
    dplyr::mutate(
      value = dplyr::na_if(.data$value, "NA"),
      value = as.numeric(.data$value),
      var = dplyr::recode(.data$var, !!!abbrev_nm)
    ) |>
    tidyr::pivot_wider(names_from = "var", values_from = "value")
  return(res)
}

#' WgsContigMeanCovFile R6 Class
#'
#' @description
#' Contains methods for reading and displaying contents of
#' the `wgs_contig_mean_cov_<phenotype>.csv` file output from DRAGEN.
#' This file contains the estimated coverage for all contigs, and an autosomal
#' estimated coverage.
#'
#' @examples
#' x1 <- system.file("extdata/wgs/SEQC-II.wgs_contig_mean_cov_normal.csv.gz", package = "dracarys")
#' x2 <- system.file("extdata/wgs/SEQC-II.wgs_contig_mean_cov_tumor.csv.gz", package = "dracarys")
#' cc1 <- WgsContigMeanCovFile$new(x1)
#' cc2 <- WgsContigMeanCovFile$new(x2)
#' d1 <- cc1$read()
#' d2 <- cc2$read()
#' cc1$write(d1, out_dir = tempdir(), prefix = "seqc_n", out_format = "tsv")
#' cc2$write(d2, out_dir = tempdir(), prefix = "seqc_t", out_format = "tsv")
#'
#' cc1$plot(d1)
#' cc2$plot(d2)
#'
#' @export
WgsContigMeanCovFile <- R6::R6Class(
  "WgsContigMeanCovFile",
  inherit = File,
  public = list(
    #' @description
    #' Reads the `wgs_contig_mean_cov_<phenotype>.csv` file output from DRAGEN.
    #'
    #' @param keep_alt Keep the ALT + Mito chromosomes?
    #' @return tibble with the following columns:
    #'   - label: file label.
    #'   - chrom: contig name.
    #'   - n_bases: number of bases aligned to contig (excludes bases from
    #'   duplicate marked reads, reads with MAPQ=0, and clipped bases).
    #'   - coverage: col2 / contig length
    read = function(keep_alt = TRUE) {
      x <- self$path
      readr::read_csv(x, col_names = c("chrom", "n_bases", "coverage"), col_types = "cdd") |>
        dplyr::filter(
          if (!keep_alt) {
            !grepl("chrM|MT|_|Autosomal|HLA-|EBV", .data$chrom)
          } else {
            TRUE
          }
        )
    },

    #' @description
    #' Writes a tidy version of the `wgs_contig_mean_cov_<phenotype>.csv` file output
    #' from DRAGEN.
    #'
    #' @param d Parsed object from `self$read()`.
    #' @param prefix Prefix of output file(s).
    #' @param out_dir Output directory.
    #' @param out_format Format of output file(s).
    #' @param drid dracarys ID to use for the dataset (e.g. `wfrid.123`, `prid.456`).
    write = function(d, out_dir = NULL, prefix, out_format = "tsv", drid = NULL) {
      if (!is.null(out_dir)) {
        prefix <- file.path(out_dir, prefix)
      }
      write_dracarys(obj = d, prefix = prefix, out_format = out_format, drid = drid)
    },


    #' @description Plots the `wgs_contig_mean_cov_<phenotype>.csv` files.
    #' @param d Parsed object from `self$read()`.
    #' @param top_alt_n Number of top covered alt contigs to plot per phenotype.
    #' @return A ggplot2 object with chromosomes on X axis, and coverage on Y axis.
    plot = function(d, top_alt_n = 15) {
      assertthat::assert_that(top_alt_n >= 0)

      # Display chr1-22, X, Y at top (M goes to bottom).
      # Display top 20 of the rest, plus rest as 'other', at bottom
      main_chrom1 <- c(1:22, "X", "Y")
      main_chrom2 <- c(paste0("chr", main_chrom1))
      main_chrom <- c(main_chrom1, main_chrom2, "Autosomal regions")
      min_cvg <- 0.000001

      d <- d |>
        dplyr::mutate(
          panel = dplyr::if_else(.data$chrom %in% main_chrom, "main", "alt"),
          panel = factor(.data$panel, levels = c("main", "alt"))
        ) |>
        dplyr::select("chrom", "coverage", "panel")

      main_panel <- d |> dplyr::filter(.data$panel == "main")
      alt_panel <- d |> dplyr::filter(.data$panel == "alt")
      top_alt <- alt_panel |>
        dplyr::top_n(top_alt_n, wt = .data$coverage) |>
        dplyr::arrange(dplyr::desc(.data$coverage)) |>
        dplyr::pull(.data$chrom) |>
        unique()

      alt_panel2 <- alt_panel |>
        dplyr::mutate(alt_group = dplyr::if_else(.data$chrom %in% top_alt, "top", "bottom"))

      alt_panel_final <- alt_panel2 |>
        dplyr::group_by(.data$alt_group) |>
        dplyr::summarise(mean_cov = mean(.data$coverage)) |>
        dplyr::inner_join(alt_panel2, by = c("alt_group")) |>
        dplyr::mutate(
          chrom = dplyr::if_else(.data$alt_group == "bottom", "OTHER", .data$chrom),
          coverage = dplyr::if_else(.data$alt_group == "bottom", .data$mean_cov, .data$coverage)
        ) |>
        dplyr::distinct() |>
        dplyr::filter(coverage > min_cvg) |>
        dplyr::ungroup() |>
        dplyr::select("chrom", "coverage", "panel")

      chrom_fac_levels <- c(main_chrom, "chrM", "MT", top_alt[!top_alt %in% c("chrM", "MT")], "OTHER")
      d <- dplyr::bind_rows(main_panel, alt_panel_final) |>
        dplyr::mutate(chrom = factor(.data$chrom, levels = chrom_fac_levels))

      d |>
        dplyr::mutate(label = "sampleA") |>
        ggplot2::ggplot(
          ggplot2::aes(
            x = .data$chrom, y = .data$coverage, group = .data$label,
          )
        ) +
        ggplot2::geom_point() +
        ggplot2::geom_line() +
        ggplot2::scale_y_continuous(
          limits = c(0, NA), expand = c(0, 0), labels = scales::comma,
          breaks = scales::pretty_breaks(n = 8)
        ) +
        ggplot2::theme_minimal() +
        ggplot2::labs(title = "Mean Coverage Per Chromosome", colour = "Label") +
        ggplot2::xlab("Chromosome") +
        ggplot2::ylab("Coverage") +
        ggplot2::theme(
          legend.position = "top",
          panel.grid.minor = ggplot2::element_blank(),
          panel.grid.major.y = ggplot2::element_blank(),
          strip.background = ggplot2::element_blank(),
          strip.text.x = ggplot2::element_blank(),
          axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1, size = 6),
          plot.title = ggplot2::element_text(colour = "#2c3e50", size = 14, face = "bold"),
          panel.spacing = ggplot2::unit(2, "lines")
        ) +
        ggplot2::facet_wrap(ggplot2::vars(.data$panel), nrow = 2, scales = "free")
    }
  )
)

#' WgsCoverageMetricsFile R6 Class
#'
#' @description
#' Contains methods for reading and displaying contents of
#' the `wgs_coverage_metrics_<phenotype>.csv` file output from DRAGEN.
#' This file contains read depth of coverage metrics.
#'
#' @examples
#' x1 <- system.file("extdata/wgs/SEQC-II.wgs_coverage_metrics_normal.csv.gz", package = "dracarys")
#' x2 <- system.file("extdata/wgs/SEQC-II.wgs_coverage_metrics_tumor.csv.gz", package = "dracarys")
#' cm1 <- WgsCoverageMetricsFile$new(x1)
#' cm2 <- WgsCoverageMetricsFile$new(x2)
#' d1 <- read(cm1)
#' d2 <- read(cm2)
#' cm1$write(d1, out_dir = tempdir(), prefix = "seqc_n", out_format = "tsv")
#' cm2$write(d2, out_dir = tempdir(), prefix = "seqc_t", out_format = "tsv")
#'
#' @export
WgsCoverageMetricsFile <- R6::R6Class(
  "WgsCoverageMetricsFile",
  inherit = File,
  public = list(
    #' @description
    #' Reads the `wgs_coverage_metrics_<phenotype>.csv` file output from DRAGEN.
    #'
    #' @return tibble with one row and metrics spread across individual columns.
    read = function() {
      abbrev_nm <- c(
        "Aligned bases"                                       = "bases_aligned_dragen",
        "Aligned bases in genome"                             = "bases_aligned_in_genome_dragen",
        "Average alignment coverage over genome"              = "cov_alignment_avg_over_genome_dragen",
        "Uniformity of coverage (PCT > 0.2*mean) over genome" = "cov_uniformity_over_genome_pct_gt02mean_dragen",
        "Uniformity of coverage (PCT > 0.4*mean) over genome" = "cov_uniformity_over_genome_pct_gt04mean_dragen",
        "Average chr X coverage over genome"                  = "cov_avg_x_over_genome_dragen",
        "Average chr Y coverage over genome"                  = "cov_avg_y_over_genome_dragen",
        "Average mitochondrial coverage over genome"          = "cov_avg_mt_over_genome_dragen",
        "Average autosomal coverage over genome"              = "cov_avg_auto_over_genome_dragen",
        "Median autosomal coverage over genome"               = "cov_median_auto_over_genome_dragen",
        "Mean/Median autosomal coverage ratio over genome"    = "cov_mean_median_auto_ratio_over_genome_dragen",
        "Aligned reads"                                       = "reads_aligned_dragen",
        "Aligned reads in genome"                             = "reads_aligned_in_genome_dragen"
      )

      x <- self$path
      raw <- readr::read_lines(x)
      assertthat::assert_that(grepl("COVERAGE SUMMARY", raw[1]))

      res <- raw |>
        tibble::as_tibble_col(column_name = "value") |>
        tidyr::separate_wider_delim(
          "value",
          delim = ",", too_few = "align_start",
          names = c("category", "dummy1", "var", "value", "pct")
        )
      # split to rename the
      # "PCT of genome with coverage [100x: inf)" values
      res1 <- res |>
        # pct just shows 100% for a couple rows
        dplyr::filter(!grepl("PCT of genome with coverage", .data$var)) |>
        dplyr::select("var", "value")
      res2 <- res |>
        dplyr::filter(grepl("PCT of genome with coverage", .data$var)) |>
        dplyr::mutate(
          var = sub("PCT of genome with coverage ", "", .data$var),
          var = gsub("\\[|\\]|\\(|\\)| ", "", .data$var),
          var = gsub("x", "", .data$var),
          var = gsub("inf", "Inf", .data$var)
        ) |>
        tidyr::separate_wider_delim("var", names = c("start", "end"), delim = ":") |>
        dplyr::mutate(var = as.character(glue("cov_genome_pct_{start}_{end}_dragen"))) |>
        dplyr::select("var", "value")
      res <- dplyr::bind_rows(res1, res2) |>
        dplyr::mutate(
          value = dplyr::na_if(.data$value, "NA"),
          value = as.numeric(.data$value),
          var = dplyr::recode(.data$var, !!!abbrev_nm)
        ) |>
        tidyr::pivot_wider(names_from = "var", values_from = "value")
    },
    #' @description
    #' Writes a tidy version of the `wgs_coverage_metrics_<phenotype>.csv` file output
    #' from DRAGEN
    #'
    #' @param d Parsed object from `self$read()`.
    #' @param prefix Prefix of output file(s).
    #' @param out_dir Output directory.
    #' @param out_format Format of output file(s).
    #' @param drid dracarys ID to use for the dataset (e.g. `wfrid.123`, `prid.456`).
    write = function(d, out_dir = NULL, prefix, out_format = "tsv", drid = NULL) {
      if (!is.null(out_dir)) {
        prefix <- file.path(out_dir, prefix)
      }
      write_dracarys(obj = d, prefix = prefix, out_format = out_format, drid = drid)
    }
  )
)

#' WgsFineHistFile R6 Class
#'
#' @description
#' Contains methods for reading and displaying contents of
#' the `wgs_fine_hist_<phenotype>.csv` file output from DRAGEN.
#' This file contains two columns: Depth and Overall.
#' The value in the Depth column ranges from 0 to 1000+ and the Overall
#' column indicates the number of loci covered at the corresponding depth.
#'
#' @examples
#' x1 <- system.file("extdata/wgs/SEQC-II.wgs_fine_hist_normal.csv.gz", package = "dracarys")
#' x2 <- system.file("extdata/wgs/SEQC-II.wgs_fine_hist_tumor.csv.gz", package = "dracarys")
#' ch1 <- WgsFineHistFile$new(x1)
#' ch2 <- WgsFineHistFile$new(x2)
#' d1 <- read(ch1)
#' d2 <- read(ch2)
#' ch1$plot(d1)
#' ch2$plot(d2)
#' ch1$write(d1, out_dir = tempdir(), prefix = "seqc_n", out_format = "tsv")
#' ch2$write(d2, out_dir = tempdir(), prefix = "seqc_t", out_format = "tsv")
#' @export
WgsFineHistFile <- R6::R6Class(
  "WgsFineHistFile",
  inherit = File,
  public = list(
    #' @description
    #' Reads the `wgs_fine_hist_<phenotype>.csv` file output from DRAGEN.
    #' @return tibble with following columns:
    #'   - depth
    #'   - number of loci with given depth
    read = function() {
      x <- self$path
      d <- readr::read_csv(x, col_types = "cd")
      assertthat::assert_that(all(colnames(d) == c("Depth", "Overall")))

      # there's a max Depth of 2000+, so convert to numeric for easier plotting
      d |>
        dplyr::mutate(
          Depth = ifelse(grepl("+", .data$Depth), sub("(\\d*)\\+", "\\1", .data$Depth), .data$Depth),
          Depth = as.integer(.data$Depth)
        ) |>
        dplyr::select(depth = "Depth", n_loci = "Overall")
    },
    #' @description
    #' Writes a tidy version of the `wgs_fine_hist_<phenotype>.csv` file output
    #' from DRAGEN
    #'
    #' @param d Parsed object from `self$read()`.
    #' @param prefix Prefix of output file(s).
    #' @param out_dir Output directory.
    #' @param out_format Format of output file(s).
    #' @param drid dracarys ID to use for the dataset (e.g. `wfrid.123`, `prid.456`).
    write = function(d, out_dir = NULL, prefix, out_format = "tsv", drid = NULL) {
      if (!is.null(out_dir)) {
        prefix <- file.path(out_dir, prefix)
      }
      write_dracarys(obj = d, prefix = prefix, out_format = out_format, drid = drid)
    },

    #' @description Plots the `wgs_fine_hist_<phenotype>.csv` files.
    #' @param d Parsed object from `self$read()`.
    #' @param x_lim X axis range to plot.
    #' @return A ggplot2 object with depth of coverage on X axis,
    #' and number of loci with that depth on Y axis.
    plot = function(d, x_lim = c(0, 300)) {
      assertthat::assert_that(length(x_lim) == 2)
      d |>
        ggplot2::ggplot(ggplot2::aes(x = .data$depth, y = .data$n_loci)) +
        ggplot2::geom_line() +
        ggplot2::coord_cartesian(xlim = x_lim) +
        ggplot2::scale_y_continuous(labels = scales::label_comma()) +
        ggplot2::scale_x_continuous(n.breaks = 8) +
        ggplot2::labs(title = "Coverage Distribution", colour = "Label") +
        ggplot2::xlab("Depth of Coverage") +
        ggplot2::ylab("Number of Loci with Given Coverage") +
        ggplot2::theme_minimal() +
        ggplot2::theme(
          legend.position = c(0.9, 0.9),
          legend.justification = c(1, 1),
          panel.grid.minor = ggplot2::element_blank(),
          panel.grid.major.y = ggplot2::element_blank(),
          plot.title = ggplot2::element_text(colour = "#2c3e50", size = 14, face = "bold")
        )
    }
  )
)

#' FragmentLengthHistFile R6 Class
#'
#' @description
#' Contains methods for reading and plotting contents of
#' the `fragment_length_hist.csv` file output from DRAGEN.
#'
#' @examples
#' x <- system.file("extdata/wgs/SEQC-II.fragment_length_hist.csv.gz", package = "dracarys")
#' fl <- FragmentLengthHistFile$new(x)
#' d <- fl$read() # or read(fl)
#' fl$plot(d) # or plot(fl)
#' fl$write(d |> dplyr::filter(count > 10), out_dir = tempdir(), prefix = "seqc_fl")
#' @export
FragmentLengthHistFile <- R6::R6Class(
  "FragmentLengthHistFile",
  inherit = File,
  public = list(
    #' @description Reads the `fragment_length_hist.csv` file, which contains the
    #' fragment length distribution for each sample.
    #' @return A tibble with the following columns:
    #' - sample: name of sample
    #' - fragmentLength: estimated fragment length
    #' - count: number of reads with estimated fragment length
    read = function() {
      x <- self$path
      d <- readr::read_lines(x)
      assertthat::assert_that(grepl("#Sample", d[1]))

      d |>
        tibble::enframe() |>
        dplyr::mutate(
          sample = dplyr::if_else(
            grepl("#Sample", .data$value),
            sub("#Sample: (.*)", "\\1", .data$value),
            NA_character_
          )
        ) |>
        tidyr::fill("sample", .direction = "down") |>
        dplyr::filter(!grepl("#Sample: |FragmentLength,Count", .data$value)) |>
        tidyr::separate_wider_delim(cols = "value", names = c("fragmentLength", "count"), delim = ",") |>
        dplyr::mutate(
          count = as.numeric(.data$count),
          fragmentLength = as.numeric(.data$fragmentLength)
        ) |>
        dplyr::select("sample", "fragmentLength", "count")
    },
    #' @description
    #' Writes a tidy version of the `fragment_length_hist.csv` file output
    #' from DRAGEN.
    #'
    #' @param d Parsed object from `self$read()`.
    #' @param prefix Prefix of output file(s).
    #' @param out_dir Output directory.
    #' @param out_format Format of output file(s).
    #' @param drid dracarys ID to use for the dataset (e.g. `wfrid.123`, `prid.456`).
    write = function(d, out_dir = NULL, prefix, out_format = "tsv", drid = NULL) {
      if (!is.null(out_dir)) {
        prefix <- file.path(out_dir, prefix)
      }
      write_dracarys(obj = d, prefix = prefix, out_format = out_format, drid = drid)
    },


    #' @description Plots the fragment length distributions as given in the
    #' `fragment_length_hist.csv` file.
    #'
    #' @param d Parsed object from `self$read()`.
    #' @param min_count Minimum read count to be plotted (Default: 10).
    #' @return A ggplot2 plot containing fragment lengths on X axis and read counts
    #'   on Y axis for each sample.
    plot = function(d, min_count = 10) {
      assertthat::assert_that(min_count >= 0)
      d <- d |>
        dplyr::filter(.data$count >= min_count)

      d |>
        ggplot2::ggplot(ggplot2::aes(x = .data$fragmentLength, y = .data$count, colour = sample)) +
        ggplot2::geom_line() +
        ggplot2::labs(title = "Fragment Length Distribution") +
        ggplot2::xlab("Fragment Length (bp)") +
        ggplot2::ylab(glue("Read Count (min: {min_count})")) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
          legend.position = c(0.9, 0.9),
          legend.justification = c(1, 1),
          panel.grid.minor = ggplot2::element_blank(),
          panel.grid.major = ggplot2::element_blank(),
          plot.title = ggplot2::element_text(colour = "#2c3e50", size = 14, face = "bold")
        )
    }
  )
)

#' MappingMetricsFile R6 Class
#'
#' @description
#' Contains methods for reading and displaying contents of
#' the `mapping_metrics.csv` file output from DRAGEN.
#' This file contains mapping and aligning metrics, like the metrics computed by
#' the Samtools Flagstat command. These metrics are available on an aggregate
#' level (over all input data), and on a per read group level. NOTE: we are
#' keeping only the read group level metrics (i.e. removing the aggregate data).
#' Unless explicitly stated, the metrics units are in reads (i.e., not in
#' terms of pairs or alignments).
#'
#' @examples
#' x <- system.file("extdata/wgs/SEQC-II.mapping_metrics.csv.gz", package = "dracarys")
#' mm <- MappingMetricsFile$new(x)
#' d <- mm$read() # or read(mm)
#' mm$write(d, out_dir = tempdir(), prefix = "seqc_mm", out_format = "tsv")
#'
#' @export
MappingMetricsFile <- R6::R6Class(
  "MappingMetricsFile",
  inherit = File,
  public = list(
    #' @description
    #' Reads the `mapping_metrics.csv` file output from DRAGEN.
    #'
    #' @return tibble with one row of X metrics per read group.
    read = function() {
      abbrev_nm <- c(
        "Total input reads" = "reads_tot_input_dragen",
        "Number of duplicate marked reads" = "reads_num_dupmarked_dragen",
        "Number of unique reads (excl. duplicate marked reads)" = "reads_num_uniq_dragen",
        "Reads with mate sequenced" = "reads_w_mate_seq_dragen",
        "Reads without mate sequenced" = "reads_wo_mate_seq_dragen",
        "QC-failed reads" = "reads_qcfail_dragen",
        "Mapped reads adjusted for excluded mapping" = "reads_mapped_adjexcl_dragen",
        "Mapped reads adjusted for filtered and excluded mapping" = "reads_mapped_adjfiltexcl_dragen",
        "Unmapped reads adjusted for excluded mapping" = "reads_unmapped_adjexcl_dragen",
        "Unmapped reads adjusted for filtered and excluded mapping" = "reads_unmapped_adjfiltexcl_dragen",
        "Reads mapping to multiple locations" = "reads_map_multiloc_dragen",
        "Hard-clipped bases R1" = "bases_hardclip_r1_dragen",
        "Hard-clipped bases R2" = "bases_hardclip_r2_dragen",
        "Soft-clipped bases" = "bases_softclip_dragen",
        "Hard-clipped bases" = "bases_hardclip_dragen",
        "Mapped reads" = "reads_mapped_dragen",
        "Mapped reads adjusted for filtered mapping" = "reads_mapped_adjfilt_dragen",
        "Mapped reads R1" = "reads_mapped_r1_dragen",
        "Mapped reads R2" = "reads_mapped_r2_dragen",
        "Number of unique & mapped reads (excl. duplicate marked reads)" = "reads_num_uniq_mapped_dragen",
        "Unmapped reads" = "reads_unmapped_dragen",
        "Unmapped reads adjusted for filtered mapping" = "reads_unmapped_adjfilt_dragen",
        "Adjustment of reads matching non-reference decoys" = "reads_match_nonref_decoys_adj_dragen",
        "Adjustment of reads matching filter contigs" = "reads_match_filt_contig_adj_dragen",
        "Singleton reads (itself mapped; mate unmapped)" = "reads_singleton_dragen",
        "Paired reads (itself & mate mapped)" = "reads_paired_dragen",
        "Properly paired reads" = "reads_paired_proper_dragen",
        "Not properly paired reads (discordant)" = "reads_discordant_dragen",
        "Paired reads mapped to different chromosomes" = "reads_paired_mapped_diff_chrom_dragen",
        "Paired reads mapped to different chromosomes (MAPQ>=10)" = "reads_paired_mapped_diff_chrom_mapq10_dragen",
        "Reads with MAPQ [40:inf)" = "reads_mapq_40_inf_dragen",
        "Reads with MAPQ [30:40)" = "reads_mapq_30_40_dragen",
        "Reads with MAPQ [20:30)" = "reads_mapq_20_30_dragen",
        "Reads with MAPQ [10:20)" = "reads_mapq_10_20_dragen",
        "Reads with MAPQ [ 0:10)" = "reads_mapq_0_10_dragen",
        "Reads with MAPQ NA (Unmapped reads)" = "reads_mapq_na_unmapped_dragen",
        "Reads with indel R1" = "reads_indel_r1_dragen",
        "Reads with indel R2" = "reads_indel_r2_dragen",
        "Reads with splice junction" = "reads_splicejunc_dragen",
        "Total bases" = "bases_tot_dragen",
        "Total bases R1" = "bases_tot_r1_dragen",
        "Total bases R2" = "bases_tot_r2_dragen",
        "Mapped bases" = "bases_mapped_dragen",
        "Mapped bases R1" = "bases_mapped_r1_dragen",
        "Mapped bases R2" = "bases_mapped_r2_dragen",
        "Soft-clipped bases R1" = "bases_softclip_r1_dragen",
        "Soft-clipped bases R2" = "bases_softclip_r2_dragen",
        "Mismatched bases R1" = "bases_mismatched_r1_dragen",
        "Mismatched bases R2" = "bases_mismatched_r2_dragen",
        "Mismatched bases R1 (excl. indels)" = "bases_mismatched_r1_noindels_dragen",
        "Mismatched bases R2 (excl. indels)" = "bases_mismatched_r2_noindels_dragen",
        "Q30 bases" = "bases_q30_dragen",
        "Q30 bases R1" = "bases_q30_r1_dragen",
        "Q30 bases R2" = "bases_q30_r2_dragen",
        "Q30 bases (excl. dups & clipped bases)" = "bases_q30_nodups_noclipped_dragen",
        "Total alignments" = "alignments_tot_dragen",
        "Secondary alignments" = "alignments_secondary_dragen",
        "Supplementary (chimeric) alignments" = "alignments_chimeric_dragen",
        "Estimated read length" = "read_len_dragen",
        "Bases in reference genome" = "bases_in_ref_genome_dragen",
        "Bases in target bed [% of genome]" = "bases_in_target_bed_genome_pct_dragen",
        "Average sequenced coverage over genome" = "cov_avg_seq_over_genome_dragen",
        "Insert length: mean" = "insert_len_mean_dragen",
        "Insert length: median" = "insert_len_median_dragen",
        "Insert length: standard deviation" = "insert_len_std_dev_dragen",
        "Provided sex chromosome ploidy" = "ploidy_sex_chrom_provided_dragen",
        "Estimated sample contamination" = "contamination_est_dragen",
        "DRAGEN mapping rate [mil. reads/second]" = "mapping_rate_dragen_milreads_per_sec_dragen",
        "Number of duplicate marked and mate reads removed" = "reads_num_dupmarked_mate_reads_removed_dragen",
        "Total reads in RG" = "reads_tot_rg_dragen",
        "Filtered rRNA reads" = "reads_rrna_filtered_dragen"
      )
      x <- self$path
      raw <- readr::read_lines(x)
      assertthat::assert_that(grepl("MAPPING/ALIGNING", raw[1]))
      # tidy
      d <- raw |>
        tibble::as_tibble_col(column_name = "value") |>
        tidyr::separate_wider_delim(
          "value",
          names = c("category", "RG", "var", "count", "pct"),
          delim = ",", too_few = "align_start"
        ) |>
        dplyr::filter(.data$RG != "") |>
        dplyr::mutate(
          count = dplyr::na_if(.data$count, "NA"),
          count = as.numeric(.data$count),
          pct = as.numeric(.data$pct),
          var = dplyr::recode(.data$var, !!!abbrev_nm)
        ) |>
        dplyr::select("RG", "var", "count", "pct")
      # pivot
      d |>
        tidyr::pivot_longer(c("count", "pct")) |>
        dplyr::mutate(
          name = dplyr::if_else(.data$name == "count", "", "_pct"),
          var = glue("{.data$var}{.data$name}")
        ) |>
        dplyr::select("RG", "var", "value") |>
        dplyr::filter(!is.na(.data$value)) |>
        tidyr::pivot_wider(names_from = "var", values_from = "value")
    },
    #' @description
    #' Writes a tidy version of the `mapping_metrics.csv` file output
    #' from DRAGEN.
    #'
    #' @param d Parsed object from `self$read()`.
    #' @param prefix Prefix of output file(s).
    #' @param out_dir Output directory.
    #' @param out_format Format of output file(s).
    #' @param drid dracarys ID to use for the dataset (e.g. `wfrid.123`, `prid.456`).
    write = function(d, out_dir = NULL, prefix, out_format = "tsv", drid = NULL) {
      if (!is.null(out_dir)) {
        prefix <- file.path(out_dir, prefix)
      }
      write_dracarys(obj = d, prefix = prefix, out_format = out_format, drid = drid)
    }
  )
)

#' PloidyEstimationMetricsFile R6 Class
#'
#' @description
#' Contains methods for reading contents of
#' the `ploidy_estimation_metrics.csv` file output from DRAGEN.
#'
#' @examples
#' x <- system.file("extdata/wgs/SEQC-II.ploidy_estimation_metrics.csv.gz", package = "dracarys")
#' pem <- PloidyEstimationMetricsFile$new(x)
#' d <- pem$read() # or read(pem)
#' pem$write(d, out_dir = tempdir(), prefix = "seqc_ploidy", out_format = "tsv")
#'
#' @export
PloidyEstimationMetricsFile <- R6::R6Class(
  "PloidyEstimationMetricsFile",
  inherit = File,
  public = list(
    #' @description
    #' Reads the `ploidy_estimation_metrics.csv` file output from DRAGEN.
    #'
    #' @return tibble with one row and metrics spread across individual columns.
    read = function() {
      x <- self$path
      raw <- readr::read_lines(x)
      assertthat::assert_that(grepl("PLOIDY ESTIMATION", raw[1]))
      abbrev_nm <- c(
        "Autosomal median coverage" = "cov_auto_median_dragen",
        "X median coverage" = "cov_x_median_dragen",
        "Y median coverage" = "cov_y_median_dragen",
        "1 median / Autosomal median" = "cov_1_div_auto_median_dragen",
        "2 median / Autosomal median" = "cov_2_div_auto_median_dragen",
        "3 median / Autosomal median" = "cov_3_div_auto_median_dragen",
        "4 median / Autosomal median" = "cov_4_div_auto_median_dragen",
        "5 median / Autosomal median" = "cov_5_div_auto_median_dragen",
        "6 median / Autosomal median" = "cov_6_div_auto_median_dragen",
        "7 median / Autosomal median" = "cov_7_div_auto_median_dragen",
        "8 median / Autosomal median" = "cov_8_div_auto_median_dragen",
        "9 median / Autosomal median" = "cov_9_div_auto_median_dragen",
        "10 median / Autosomal median" = "cov_10_div_auto_median_dragen",
        "11 median / Autosomal median" = "cov_11_div_auto_median_dragen",
        "12 median / Autosomal median" = "cov_12_div_auto_median_dragen",
        "13 median / Autosomal median" = "cov_13_div_auto_median_dragen",
        "14 median / Autosomal median" = "cov_14_div_auto_median_dragen",
        "15 median / Autosomal median" = "cov_15_div_auto_median_dragen",
        "16 median / Autosomal median" = "cov_16_div_auto_median_dragen",
        "17 median / Autosomal median" = "cov_17_div_auto_median_dragen",
        "18 median / Autosomal median" = "cov_18_div_auto_median_dragen",
        "19 median / Autosomal median" = "cov_19_div_auto_median_dragen",
        "20 median / Autosomal median" = "cov_20_div_auto_median_dragen",
        "21 median / Autosomal median" = "cov_21_div_auto_median_dragen",
        "22 median / Autosomal median" = "cov_22_div_auto_median_dragen",
        "X median / Autosomal median" = "cov_x_div_auto_median_dragen",
        "Y median / Autosomal median" = "cov_y_div_auto_median_dragen",
        "Ploidy estimation" = "ploidy_est_dragen"
      )

      d <- raw |>
        tibble::as_tibble_col(column_name = "value") |>
        tidyr::separate_wider_delim("value", names = c("dummy1", "dummy2", "var", "value"), delim = ",") |>
        dplyr::select("var", "value") |>
        dplyr::mutate(
          var = dplyr::recode(.data$var, !!!abbrev_nm)
        ) |>
        tidyr::pivot_wider(names_from = "var", values_from = "value")
      # now convert all except 'Ploidy estimation' to numeric
      cols1 <- colnames(d)[colnames(d) != "ploidy_est_dragen"]
      d |>
        dplyr::mutate(dplyr::across(dplyr::all_of(cols1), as.numeric))
    },
    #' @description
    #' Writes a tidy version of the `ploidy_estimation_metrics.csv` file output
    #' from DRAGEN.
    #'
    #' @param d Parsed object from `self$read()`.
    #' @param prefix Prefix of output file(s).
    #' @param out_dir Output directory.
    #' @param out_format Format of output file(s).
    #' @param drid dracarys ID to use for the dataset (e.g. `wfrid.123`, `prid.456`).
    write = function(d, out_dir = NULL, prefix, out_format = "tsv", drid = NULL) {
      if (!is.null(out_dir)) {
        prefix <- file.path(out_dir, prefix)
      }
      write_dracarys(obj = d, prefix = prefix, out_format = out_format, drid = drid)
    }
  )
)

#' ReplayFile R6 Class
#'
#' @description
#' Contains methods for reading contents of
#' the `replay.json` file output from DRAGEN, which contains the DRAGEN command
#' line, parameters and version for the specific run.
#'
#' @examples
#' x <- system.file("extdata/wgs/SEQC-II-replay.json.gz", package = "dracarys")
#' r <- ReplayFile$new(x)
#' d <- r$read() # or read(r)
#' r$write(d, out_dir = tempdir(), prefix = "seqc_replay", out_format = "tsv")
#' @export
ReplayFile <- R6::R6Class(
  "ReplayFile",
  inherit = File,
  public = list(
    #' @description Reads the `replay.json` file.
    #' @return tibble with one row and metrics spread across individual columns.
    read = function() {
      x <- self$path
      res <- x |>
        jsonlite::read_json(simplifyVector = TRUE) |>
        purrr::map_if(is.data.frame, tibble::as_tibble)

      req_elements <- c("command_line", "hash_table_build", "dragen_config", "system")
      assertthat::assert_that(all(names(res) %in% req_elements))
      res[["system"]] <- res[["system"]] |>
        tibble::as_tibble_row()
      res[["hash_table_build"]] <- res[["hash_table_build"]] |>
        tibble::as_tibble_row()
      # we don't care if the columns are characters, no analysis likely to be done on dragen options
      # (though never say never!)
      res[["dragen_config"]] <- res[["dragen_config"]] |>
        tidyr::pivot_wider(names_from = "name", values_from = "value")

      dplyr::bind_cols(res)
    },
    #' @description
    #' Writes a tidy version of the `replay.json` file output
    #' from DRAGEN.
    #'
    #' @param d Parsed object from `self$read()`.
    #' @param prefix Prefix of output file(s).
    #' @param out_dir Output directory.
    #' @param out_format Format of output file(s).
    #' @param drid dracarys ID to use for the dataset (e.g. `wfrid.123`, `prid.456`).
    write = function(d, out_dir = NULL, prefix, out_format = "tsv", drid = NULL) {
      if (!is.null(out_dir)) {
        prefix <- file.path(out_dir, prefix)
      }
      write_dracarys(obj = d, prefix = prefix, out_format = out_format, drid = drid)
    }
  )
)

#' TimeMetricsFile R6 Class
#'
#' @description
#' Contains methods for reading contents of
#' the `time_metrics.csv` file output from DRAGEN, which contains
#' a breakdown of the run duration for each DRAGEN process.
#'
#' @examples
#' x <- system.file("extdata/wgs/SEQC-II.time_metrics.csv.gz", package = "dracarys")
#' tm <- TimeMetricsFile$new(x)
#' d <- tm$read() # or read(tm)
#' tm$write(d, out_dir = tempdir(), prefix = "seqc_time", out_format = "tsv")
#' @export
TimeMetricsFile <- R6::R6Class(
  "TimeMetricsFile",
  inherit = File,
  public = list(
    #' @description Reads the `time_metrics.csv` file.
    #' @return tibble with one row and metrics spread across individual columns.
    read = function() {
      x <- self$path
      cn <- c("dummy1", "dummy2", "Step", "time_hrs", "time_sec")
      ct <- readr::cols(.default = "c", time_hrs = readr::col_time(format = "%T"), time_sec = "d")
      d <- readr::read_csv(x, col_names = cn, col_types = ct)
      assertthat::assert_that(d$dummy1[1] == "RUN TIME", is.na(d$dummy2[1]))
      assertthat::assert_that(inherits(d$time_hrs, "hms"))
      d |>
        dplyr::mutate(
          Step = tools::toTitleCase(sub("Time ", "", .data$Step)),
          Time = substr(.data$time_hrs, 1, 5)
        ) |>
        dplyr::select("Step", "Time") |>
        tidyr::pivot_wider(names_from = "Step", values_from = "Time") |>
        dplyr::relocate("Total Runtime")
    },
    #' @description
    #' Writes a tidy version of the `time_metrics.csv` file output
    #' from DRAGEN.
    #'
    #' @param d Parsed object from `self$read()`.
    #' @param prefix Prefix of output file(s).
    #' @param out_dir Output directory.
    #' @param out_format Format of output file(s).
    #' @param drid dracarys ID to use for the dataset (e.g. `wfrid.123`, `prid.456`).
    write = function(d, out_dir = NULL, prefix, out_format = "tsv", drid = NULL) {
      if (!is.null(out_dir)) {
        prefix <- file.path(out_dir, prefix)
      }
      write_dracarys(obj = d, prefix = prefix, out_format = out_format, drid = drid)
    }
  )
)

#' Process Multiple TimeMetricsFile Objects
#'
#' Processes multiple TimeMetricsFile objects.
#'
#' @param x Atomic vector with one or more TimeMetricsFile objects.
#' @param id ID for each input, which is used to disambiguate files
#' generated from same samples. Default: index from 1 to length of `x`.
#' @return tibble with the following columns:
#'   - Step: DRAGEN step
#'   - Time: time in HH:MM
#'
#' @examples
#' p <- system.file("extdata/wgs/SEQC-II.time_metrics.csv.gz", package = "dracarys")
#' x <- TimeMetricsFile$new(p)
#' (tm <- time_metrics_process(c(x, x), id = c("run1", "run2")))
#'
#' @testexamples
#' expect_equal(nrow(tm), 2)
#'
#' @export
time_metrics_process <- function(x, id = seq_len(length(x))) {
  assertthat::assert_that(all(purrr::map_lgl(x, ~ inherits(.x, "TimeMetricsFile"))))
  x |>
    purrr::map(read) |>
    purrr::set_names(id) |>
    dplyr::bind_rows(.id = "ID")
}

#' VCMetricsFile R6 Class
#'
#' @description
#' Contains methods for reading and displaying contents of
#' the `vc_metrics.csv` file output from DRAGEN, which contains variant calling metrics.
#' @examples
#' x <- system.file("extdata/wgs/SEQC-II.vc_metrics.csv.gz", package = "dracarys")
#' vm <- VCMetricsFile$new(x)
#' d <- vm$read() # or read(vm)
#' vm$write(d, out_dir = tempdir(), prefix = "seqc_vc", out_format = "tsv")
#'
#' @export
VCMetricsFile <- R6::R6Class(
  "VCMetricsFile",
  inherit = File,
  public = list(
    #' @description
    #' Reads the `vc_metrics.csv` file output from DRAGEN.
    #'
    #' @return tibble with one row and metrics spread across individual columns.
    read = function() {
      abbrev_nm <- c(
        "Total" = "var_tot_dragen",
        "Biallelic" = "var_biallelic_dragen",
        "Multiallelic" = "var_multiallelic_dragen",
        "SNPs" = "var_snp_dragen",
        "Insertions (Hom)" = "var_ins_hom_dragen",
        "Insertions (Het)" = "var_ins_het_dragen",
        "Deletions (Hom)" = "var_del_hom_dragen",
        "Deletions (Het)" = "var_del_het_dragen",
        "Indels (Het)" = "var_indel_het_dragen",
        "Chr X number of SNPs over genome" = "var_snp_x_over_genome_dragen",
        "Chr Y number of SNPs over genome" = "var_snp_y_over_genome_dragen",
        "(Chr X SNPs)/(chr Y SNPs) ratio over genome" = "var_x_over_y_snp_ratio_over_genome_dragen",
        "SNP Transitions" = "var_snp_transitions_dragen",
        "SNP Transversions" = "var_snp_transversions_dragen",
        "Ti/Tv ratio" = "var_ti_tv_ratio_dragen",
        "Heterozygous" = "var_heterozygous_dragen",
        "Homozygous" = "var_homozygous_dragen",
        "Het/Hom ratio" = "var_het_hom_ratio_dragen",
        "In dbSNP" = "var_in_dbsnp_dragen",
        "Not in dbSNP" = "var_nin_dbsnp_dragen",
        "Percent Callability" = "callability_pct_dragen",
        "Percent Autosome Callability" = "callability_auto_pct_dragen",
        "Number of samples" = "sample_num_dragen",
        "Reads Processed" = "reads_processed_dragen",
        "Child Sample" = "sample_child_dragen"
      )
      x <- self$path
      raw <- readr::read_lines(x)
      assertthat::assert_that(grepl("VARIANT CALLER", raw[1]))
      # tidy
      d <- raw |>
        tibble::as_tibble_col(column_name = "value") |>
        tidyr::separate_wider_delim(
          "value",
          names = c("category", "sample", "var", "count", "pct"),
          delim = ",", too_few = "align_start"
        ) |>
        dplyr::mutate(
          var = dplyr::recode(.data$var, !!!abbrev_nm),
          count = dplyr::na_if(.data$count, "NA"),
          count = as.numeric(.data$count),
          pct = round(as.numeric(.data$pct), 2),
          category = dplyr::case_when(
            grepl("SUMMARY", .data$category) ~ "summary",
            grepl("PREFILTER", .data$category) ~ "prefilter",
            grepl("POSTFILTER", .data$category) ~ "postfilter",
            TRUE ~ "unknown"
          )
        ) |>
        dplyr::filter(.data$category != "summary") |>
        dplyr::select("category", "sample", "var", "count", "pct")
      # pivot
      d |>
        tidyr::pivot_longer(c("count", "pct")) |>
        dplyr::mutate(
          name = dplyr::if_else(.data$name == "count", "", "_pct"),
          var = glue("{.data$var}{.data$name}")
        ) |>
        dplyr::select("category", "sample", "var", "value") |>
        dplyr::filter(!is.na(.data$value)) |>
        tidyr::pivot_wider(names_from = "var", values_from = "value")
    },
    #' @description
    #' Writes a tidy version of the `vc_metrics.csv` file output
    #' from DRAGEN.
    #'
    #' @param d Parsed object from `self$read()`.
    #' @param prefix Prefix of output file(s).
    #' @param out_dir Output directory.
    #' @param out_format Format of output file(s).
    #' @param drid dracarys ID to use for the dataset (e.g. `wfrid.123`, `prid.456`).
    write = function(d, out_dir = NULL, prefix, out_format = "tsv", drid = NULL) {
      if (!is.null(out_dir)) {
        prefix <- file.path(out_dir, prefix)
      }
      write_dracarys(obj = d, prefix = prefix, out_format = out_format, drid = drid)
    }
  )
)

#' TrimmerMetricsFile R6 Class
#'
#' @description
#' Contains methods for reading and displaying contents of
#' the `trimmer_metrics.csv` file output from DRAGEN
#'
#' @examples
#' x <- system.file("extdata/wgs/SEQC-II.trimmer_metrics.csv.gz", package = "dracarys")
#' tm <- TrimmerMetricsFile$new(x)
#' d <- tm$read()
#' tm$write(d, out_dir = tempdir(), prefix = "seqc_tm", out_format = "tsv")
#'
#' @export
TrimmerMetricsFile <- R6::R6Class(
  "TrimmerMetricsFile",
  inherit = File,
  public = list(
    #' @description
    #' Reads the `trimmer_metrics.csv` file output from DRAGEN.
    #'
    #' @return tibble with one row and metrics spread across individual columns.
    read = function() {
      x <- self$path
      d <- readr::read_lines(x)
      assertthat::assert_that(grepl("TRIMMER STATISTICS", d[1]))
      abbrev_nm <- c(
        "Total input reads"                              = "reads_tot_input_dragen",
        "Total input bases"                              = "bases_tot_dragen",
        "Total input bases R1"                           = "bases_r1_dragen",
        "Total input bases R2"                           = "bases_r2_dragen",
        "Average input read length"                      = "read_len_avg_dragen",
        "Total trimmed reads"                            = "reads_trimmed_tot_dragen",
        "Total trimmed bases"                            = "bases_trimmed_tot_dragen",
        "Average bases trimmed per read"                 = "bases_trimmed_avg_per_read_dragen",
        "Average bases trimmed per trimmed read"         = "bases_trimmed_avg_per_trimmedread_dragen",
        "Remaining poly-G K-mers R1 3prime"              = "polygkmers3r1_remaining_dragen",
        "Remaining poly-G K-mers R2 3prime"              = "polygkmers3r2_remaining_dragen",
        "Poly-G soft trimmed reads unfiltered R1 3prime" = "polyg_soft_trimmed_reads_unfilt_3r1_dragen",
        "Poly-G soft trimmed reads unfiltered R2 3prime" = "polyg_soft_trimmed_reads_unfilt_3r2_dragen",
        "Poly-G soft trimmed reads filtered R1 3prime"   = "polyg_soft_trimmed_reads_filt_3r1_dragen",
        "Poly-G soft trimmed reads filtered R2 3prime"   = "polyg_soft_trimmed_reads_filt_3r2_dragen",
        "Poly-G soft trimmed bases unfiltered R1 3prime" = "polyg_soft_trimmed_bases_unfilt_3r1_dragen",
        "Poly-G soft trimmed bases unfiltered R2 3prime" = "polyg_soft_trimmed_bases_unfilt_3r2_dragen",
        "Poly-G soft trimmed bases filtered R1 3prime"   = "polyg_soft_trimmed_bases_filt_3r1_dragen",
        "Poly-G soft trimmed bases filtered R2 3prime"   = "polyg_soft_trimmed_bases_filt_3r2_dragen",
        "Total filtered reads"                           = "reads_tot_filt_dragen",
        "Reads filtered for minimum read length R1"      = "reads_filt_minreadlenr1_dragen",
        "Reads filtered for minimum read length R2"      = "reads_filt_minreadlenr2_dragen"
      )

      d |>
        tibble::as_tibble_col(column_name = "value") |>
        tidyr::separate_wider_delim("value", names = c("category", "extra", "var", "count", "pct"), delim = ",", too_few = "align_start") |>
        dplyr::mutate(
          count = as.numeric(.data$count),
          pct = round(as.numeric(.data$pct), 2),
          var = dplyr::recode(.data$var, !!!abbrev_nm)
        ) |>
        dplyr::select("var", "count", "pct") |>
        tidyr::pivot_longer(c("count", "pct")) |>
        dplyr::filter(!is.na(.data$value)) |>
        dplyr::mutate(
          name = dplyr::if_else(.data$name == "count", "", "_pct"),
          var = glue("{.data$var}{.data$name}")
        ) |>
        dplyr::select("var", "value") |>
        tidyr::pivot_wider(names_from = "var", values_from = "value")
    },
    #' @description
    #' Writes a tidy version of the `trimmer_metrics.csv` file output
    #' from DRAGEN.
    #'
    #' @param d Parsed object from `self$read()`.
    #' @param prefix Prefix of output file(s).
    #' @param out_dir Output directory.
    #' @param out_format Format of output file(s).
    #' @param drid dracarys ID to use for the dataset (e.g. `wfrid.123`, `prid.456`).
    write = function(d, out_dir = NULL, prefix, out_format = "tsv", drid = NULL) {
      if (!is.null(out_dir)) {
        prefix <- file.path(out_dir, prefix)
      }
      write_dracarys(obj = d, prefix = prefix, out_format = out_format, drid = drid)
    }
  )
)

#' SvMetricsFile R6 Class
#'
#' @description
#' Contains methods for reading and displaying contents of
#' the `sv_metrics.csv` file output from DRAGEN
#'
#' @examples
#' x <- system.file("extdata/wgs/SEQC-II.sv_metrics.csv.gz", package = "dracarys")
#' sv <- SvMetricsFile$new(x)
#' d <- sv$read()
#' sv$write(d, out_dir = tempdir(), prefix = "seqc_sv", out_format = "tsv")
#'
#' @export
SvMetricsFile <- R6::R6Class(
  "SvMetricsFile",
  inherit = File,
  public = list(
    #' @description
    #' Reads the `sv_metrics.csv` file output from DRAGEN.
    #'
    #' @return tibble with one row and metrics spread across individual columns.
    read = function() {
      x <- self$path
      d <- readr::read_lines(x)
      assertthat::assert_that(grepl("SV SUMMARY", d[1]))
      abbrev_nm <- c(
        "Number of deletions (PASS)" = "n_del",
        "Number of insertions (PASS)" = "n_ins",
        "Number of duplications (PASS)" = "n_dup",
        "Number of breakend pairs (PASS)" = "n_bnd"
      )
      d |>
        tibble::as_tibble_col(column_name = "value") |>
        dplyr::filter(!grepl("Total number of structural variants", .data$value)) |>
        tidyr::separate_wider_delim(
          "value",
          names = c("svsum", "sample", "var", "count", "pct"), delim = ",",
          too_few = "align_start"
        ) |>
        dplyr::mutate(
          count = as.numeric(.data$count),
          pct = round(as.numeric(.data$pct), 2),
          var = dplyr::recode(.data$var, !!!abbrev_nm)
        ) |>
        dplyr::select("var", "count", "pct") |>
        tidyr::pivot_longer(c("count", "pct")) |>
        dplyr::mutate(
          name = dplyr::if_else(.data$name == "count", "", "_pct"),
          var = glue("{.data$var}{.data$name}")
        ) |>
        dplyr::arrange(.data$name, .data$var) |>
        dplyr::select("var", "value") |>
        tidyr::pivot_wider(names_from = "var", values_from = "value")
    },
    #' @description
    #' Writes a tidy version of the `sv_metrics.csv` file output
    #' from DRAGEN.
    #'
    #' @param d Parsed object from `self$read()`.
    #' @param prefix Prefix of output file(s).
    #' @param out_dir Output directory.
    #' @param out_format Format of output file(s).
    #' @param drid dracarys ID to use for the dataset (e.g. `wfrid.123`, `prid.456`).
    write = function(d, out_dir = NULL, prefix, out_format = "tsv", drid = NULL) {
      if (!is.null(out_dir)) {
        prefix <- file.path(out_dir, prefix)
      }
      write_dracarys(obj = d, prefix = prefix, out_format = out_format, drid = drid)
    }
  )
)

#' WgsHistFile R6 Class
#'
#' @description
#' Contains methods for reading and displaying contents of
#' the `wgs_hist.csv` file output from DRAGEN
#'
#' @examples
#' x <- system.file("extdata/wgs/SEQC-II.wgs_hist.csv.gz", package = "dracarys")
#' h <- WgsHistFile$new(x)
#' d <- h$read()
#' h$write(d, out_dir = tempdir(), prefix = "seqc_sv", out_format = "tsv")
#'
#' @export
WgsHistFile <- R6::R6Class(
  "WgsHistFile",
  inherit = File,
  public = list(
    #' @description
    #' Reads the `wgs_hist.csv` file output from DRAGEN.
    #'
    #' @return tibble with one row and metrics spread across individual columns.
    read = function() {
      x <- self$path
      d <- readr::read_csv(x, col_names = c("var", "pct"), col_types = "cd")
      d |>
        dplyr::mutate(
          var = sub("PCT of bases in wgs with coverage ", "", .data$var),
          var = gsub("\\[|\\]|\\(|\\)", "", .data$var),
          var = gsub("x", "", .data$var),
          var = gsub("inf", "Inf", .data$var)
        ) |>
        tidyr::separate_wider_delim("var", names = c("start", "end"), delim = ":") |>
        dplyr::mutate(
          start = as.numeric(.data$start),
          end = as.numeric(.data$end),
          pct = round(.data$pct, 2),
          cumsum = cumsum(.data$pct)
        )
    },
    #' @description
    #' Writes a tidy version of the `wgs_hist.csv` file output
    #' from DRAGEN.
    #'
    #' @param d Parsed object from `self$read()`.
    #' @param prefix Prefix of output file(s).
    #' @param out_dir Output directory.
    #' @param out_format Format of output file(s).
    #' @param drid dracarys ID to use for the dataset (e.g. `wfrid.123`, `prid.456`).
    write = function(d, out_dir = NULL, prefix, out_format = "tsv", drid = NULL) {
      if (!is.null(out_dir)) {
        prefix <- file.path(out_dir, prefix)
      }
      write_dracarys(obj = d, prefix = prefix, out_format = out_format, drid = drid)
    }
  )
)

dragen_subprefix <- function(x, suffix) {
  # L2401290.exon_contig_mean_cov.csv -> exon
  # L2401290.target_bed_contig_mean_cov.csv -> target_bed
  # L2401290.tmb_contig_mean_cov.csv -> tmb
  # L2401290.wgs_contig_mean_cov.csv -> wgs
  bname <- basename(x)
  s1 <- tools::file_path_sans_ext(bname)
  s2 <- sub(".*\\.(.*)", "\\1", s1)
  sub(suffix, "", s2)
}
