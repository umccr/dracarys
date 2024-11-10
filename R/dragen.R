#' DRAGEN Fragment Length Hist Plot
#'
#' Plots the fragment length distributions as given in the
#' `fragment_length_hist` file.
#'
#' @param d Parsed tibble.
#' @param min_count Minimum read count to be plotted (def: 10).
#'
#' @return A ggplot2 plot containing fragment lengths on X axis and read counts
#' on Y axis for each sample.
dragen_fraglenhist_plot <- function(d, min_count = 10) {
  assertthat::assert_that(is.numeric(min_count), min_count >= 0)
  d |>
    dplyr::filter(.data$Count >= min_count) |>
    ggplot2::ggplot(ggplot2::aes(x = .data$FragmentLength, y = .data$Count)) +
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

#' Read DRAGEN UMI Metrics
#'
#' Reads the `umi_metrics.csv` file output from DRAGEN.
#'
#' @param x Path to file.
#'
#' @return Tibble with metrics.
#' @examples
#' \dontrun{
#' dragen_umi_metrics_read(x)
#' }
#' @export
dragen_umi_metrics_read <- function(x) {
  d0 <- readr::read_lines(x)
  assertthat::assert_that(grepl("UMI STATISTICS", d0[1]))
  abbrev_nm <- tibble::tribble(
    ~raw, ~clean, ~target,
    "number of reads", "reads_tot", TRUE,
    "number of reads with valid or correctable umis", "reads_umi_valid_correctable", TRUE,
    "number of reads in discarded families", "reads_discarded_families", TRUE,
    "reads filtered out", "reads_filtered_out", FALSE,
    "reads with all-g umis filtered out", "reads_filtered_out_all_g_umis", FALSE,
    "reads with uncorrectable umis", "reads_uncorrectable_umis", FALSE,
    "total number of families", "families_tot", FALSE,
    "families contextually corrected", "families_contextually_corrected", FALSE,
    "families shifted", "families_shifted", FALSE,
    "families discarded", "families_discarded_tot", TRUE,
    "families discarded by min-support-reads", "families_discarded_minsupportreads", TRUE,
    "families discarded by duplex/simplex", "families_discarded_duplexsimplex", TRUE,
    "families with ambiguous correction", "families_ambiguous_correction", TRUE,
    "duplex families", "duplex_families", TRUE,
    "consensus pairs emitted", "consensus_pairs_emitted", FALSE,
    "mean family depth", "avg_family_depth", TRUE,
    "number of collapsible regions", "collapsible_regions_tot", FALSE,
    "min collapsible region size (num reads)", "collapsible_region_size_min", FALSE,
    "max collapsible region size (num reads)", "collapsible_region_size_max", FALSE,
    "mean collapsible region size (num reads)", "collapsible_region_size_mean", FALSE,
    "collapsible region size standard deviation", "collapsible_region_size_sd", FALSE,
    "histogram of num supporting fragments", "histo_num_supporting_fragments", TRUE,
    "histogram of unique umis per fragment position", "histo_unique_umis", FALSE
  )
  abbrev_nm_target <- abbrev_nm |>
    dplyr::filter(.data$target) |>
    dplyr::mutate(
      raw = as.character(glue("on target {.data$raw}")),
      clean = as.character(glue("{.data$clean}_ontarget"))
    )
  abbrev_nm <- abbrev_nm |>
    dplyr::bind_rows(abbrev_nm_target) |>
    dplyr::select("raw", "clean") |>
    tibble::deframe()
  d1 <- d0 |>
    tibble::as_tibble_col(column_name = "value") |>
    tidyr::separate_wider_delim(
      "value",
      names = c("category", "dummy1", "var", "count", "pct"),
      delim = ",", too_few = "align_start"
    ) |>
    dplyr::mutate(
      var = tolower(.data$var),
      var = dplyr::recode(.data$var, !!!abbrev_nm)
    )
  dirty_names_cleaned(d1$var, abbrev_nm, x)
  hist <- d1 |>
    dplyr::filter(grepl("histo", .data$var)) |>
    dplyr::select(name = "var", "count") |>
    dplyr::mutate(
      count = gsub("\\{|\\}", "", .data$count),
      count = strsplit(.data$count, "\\|")
    ) |>
    tidyr::unnest("count") |>
    dplyr::mutate(count = as.numeric(.data$count)) |>
    tidyr::nest(.by = "name")
  d2 <- d1 |>
    dplyr::filter(!grepl("histo", .data$var)) |>
    dplyr::select("var", "count", "pct") |>
    tidyr::pivot_longer(c("count", "pct"), names_to = "name", values_to = "value") |>
    dplyr::mutate(
      value = as.numeric(.data$value),
      name = dplyr::if_else(.data$name == "count", "", "_pct"),
      var = glue("{.data$var}{.data$name}")
    ) |>
    dplyr::select("var", "value") |>
    dplyr::filter(!is.na(.data$value)) |>
    tidyr::pivot_wider(names_from = "var", values_from = "value")
  res <- list(metrics = d2) |>
    tibble::enframe(name = "name", value = "data") |>
    dplyr::bind_rows(hist) |>
    dplyr::mutate(name = glue("umi_{.data$name}"))
  res
}

#' Read DRAGEN GC Metrics
#'
#' Reads the `gc_metrics.csv` file output from DRAGEN.
#'
#' @param x Path to file.
#'
#' @return Tibble with metrics.
#' @examples
#' \dontrun{
#' dragen_gc_metrics_read(x)
#' }
#'
#' @export
dragen_gc_metrics_read <- function(x) {
  d0 <- readr::read_lines(x)
  assertthat::assert_that(grepl("GC BIAS DETAILS", d0[1]))
  d1 <- d0 |>
    tibble::as_tibble_col(column_name = "value") |>
    tidyr::separate_wider_delim(
      "value",
      names = c("category", "dummy1", "name", "value", "fraction"),
      delim = ",", too_few = "align_start"
    )
  abbrev_nm <- c(
    "Window size" = "window_size",
    "Number of valid windows" = "n_windows_valid",
    "Number of discarded windows" = "n_windows_discard",
    "Average reference GC" = "gc_avg_reference",
    "Mean global coverage" = "cov_avg_global",
    "Normalized coverage at GCs 0-19" = "cov_norm_0to19",
    "Normalized coverage at GCs 20-39" = "cov_norm_20to39",
    "Normalized coverage at GCs 40-59" = "cov_norm_40to59",
    "Normalized coverage at GCs 60-79" = "cov_norm_60to79",
    "Normalized coverage at GCs 80-100" = "cov_norm_80to100",
    "AT Dropout" = "dropout_at",
    "GC Dropout" = "dropout_gc"
  )
  # GC METRICS SUMMARY
  summary <- d1 |>
    dplyr::filter(.data$category == "GC METRICS SUMMARY") |>
    dplyr::select("name", "value") |>
    dplyr::mutate(
      value = as.numeric(.data$value),
      name = dplyr::recode(.data$name, !!!abbrev_nm)
    ) |>
    tidyr::pivot_wider(names_from = "name", values_from = "value")
  dirty_names_cleaned(colnames(summary), abbrev_nm, x)

  # GC BIAS DETAILS
  details <- d1 |>
    dplyr::filter(.data$category == "GC BIAS DETAILS") |>
    dplyr::mutate(
      gc = sub(".* at GC (.*)", "\\1", .data$name),
      name = sub("(.*) at GC .*", "\\1", .data$name),
      name = tolower(.data$name),
      name = sub(" ", "", .data$name),
      value = as.numeric(.data$value),
      fraction = as.numeric(.data$fraction)
    )
  details_wind <- details |>
    dplyr::filter(.data$name == "windows") |>
    dplyr::select("gc", "value")

  details_cov <- details |>
    dplyr::filter(.data$name == "normalizedcoverage") |>
    dplyr::select("gc", "value")

  list(
    gcmetrics_summary = summary,
    gcmetrics_windows = details_wind,
    gcmetrics_coverage = details_cov
  ) |>
    tibble::enframe(name = "name", value = "data")
}

#' Read DRAGEN CNV Metrics
#'
#' Reads the `cnv_metrics.csv` file output from DRAGEN.
#'
#' @param x Path to file.
#'
#' @return Tibble with metrics.
#' @examples
#' \dontrun{
#' dragen_cnv_metrics_read(x)
#' }
#'
#' @export
dragen_cnv_metrics_read <- function(x) {
  d0 <- readr::read_lines(x)
  # first row is sometimes SEX GENOTYPER, others CNV SUMMARY
  assertthat::assert_that(grepl("CNV SUMMARY", d0[2]))
  abbrev_nm <- c(
    "Bases in reference genome" = "bases_in_ref_genome",
    "Average alignment coverage over genome" = "cov_alignment_avg_over_genome",
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
    "Number of passing deletions" = "n_del_pass",
    "Estimated tumor purity" = "purity_tumor",
    "Diploid coverage" = "cov_diploid",
    "Overall ploidy" = "ploidy_overall"
  )
  d1 <- d0 |>
    tibble::as_tibble_col(column_name = "value") |>
    tidyr::separate_wider_delim(
      "value",
      names = c("category", "extra", "var", "count", "pct"),
      delim = ",", too_few = "align_start"
    )
  # in cttso
  sexgt <- d1 |>
    dplyr::filter(.data$category == "SEX GENOTYPER") |>
    dplyr::select(sexgt = "count", sexgt_pct = "pct")
  d2 <- d1 |>
    dplyr::filter(!.data$category == "SEX GENOTYPER") |>
    dplyr::mutate(
      count = as.numeric(.data$count),
      pct = round(as.numeric(.data$pct), 2),
      var = dplyr::recode(.data$var, !!!abbrev_nm)
    )
  dirty_names_cleaned(d2$var, abbrev_nm, x)
  d2 <- d2 |>
    dplyr::select("var", "count", "pct") |>
    tidyr::pivot_longer(c("count", "pct")) |>
    dplyr::filter(!is.na(.data$value)) |>
    dplyr::mutate(
      name = dplyr::if_else(.data$name == "count", "", "_pct"),
      var = glue("{.data$var}{.data$name}")
    ) |>
    dplyr::select("var", "value") |>
    tidyr::pivot_wider(names_from = "var", values_from = "value")
  if (nrow(sexgt) == 0) {
    return(d2)
  }
  dplyr::bind_cols(sexgt, d2)
}

#' Read DRAGEN SV Metrics
#'
#' Reads the `sv_metrics.csv` file output from DRAGEN.
#'
#' @param x Path to file.
#'
#' @return Tibble with metrics.
#' @examples
#' \dontrun{
#' dragen_sv_metrics_read(x)
#' }
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
  res <- d |>
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
  dirty_names_cleaned(colnames(res), abbrev_nm, x)
  res
}

#' Read DRAGEN Trimmer Metrics
#'
#' Reads the `trimmer_metrics.csv` file output from DRAGEN.
#'
#' @param x Path to file.
#'
#' @return Tibble with metrics.
#' @examples
#' \dontrun{
#' dragen_trimmer_metrics_read(x)
#' }
#' @export
dragen_trimmer_metrics_read <- function(x) {
  d0 <- readr::read_lines(x)
  assertthat::assert_that(grepl("TRIMMER STATISTICS", d0[1]))
  abbrev_nm <- c(
    "Total input reads"                              = "reads_tot_input",
    "Total input bases"                              = "bases_tot",
    "Total input bases R1"                           = "bases_r1",
    "Total input bases R2"                           = "bases_r2",
    "Average input read length"                      = "read_len_avg",
    "Total trimmed reads"                            = "reads_trimmed_tot",
    "Total trimmed bases"                            = "bases_trimmed_tot",
    "Average bases trimmed per read"                 = "bases_trimmed_avg_per_read",
    "Average bases trimmed per trimmed read"         = "bases_trimmed_avg_per_trimmedread",
    "Remaining poly-G K-mers R1 3prime"              = "polygkmers3r1_remaining",
    "Remaining poly-G K-mers R2 3prime"              = "polygkmers3r2_remaining",
    "Poly-G soft trimmed reads unfiltered R1 3prime" = "polyg_soft_trimmed_reads_unfilt_3r1",
    "Poly-G soft trimmed reads unfiltered R2 3prime" = "polyg_soft_trimmed_reads_unfilt_3r2",
    "Poly-G soft trimmed reads filtered R1 3prime"   = "polyg_soft_trimmed_reads_filt_3r1",
    "Poly-G soft trimmed reads filtered R2 3prime"   = "polyg_soft_trimmed_reads_filt_3r2",
    "Poly-G soft trimmed bases unfiltered R1 3prime" = "polyg_soft_trimmed_bases_unfilt_3r1",
    "Poly-G soft trimmed bases unfiltered R2 3prime" = "polyg_soft_trimmed_bases_unfilt_3r2",
    "Poly-G soft trimmed bases filtered R1 3prime"   = "polyg_soft_trimmed_bases_filt_3r1",
    "Poly-G soft trimmed bases filtered R2 3prime"   = "polyg_soft_trimmed_bases_filt_3r2",
    "Total filtered reads"                           = "reads_tot_filt",
    "Reads filtered for minimum read length R1"      = "reads_filt_minreadlenr1",
    "Reads filtered for minimum read length R2"      = "reads_filt_minreadlenr2"
  )

  d1 <- d0 |>
    tibble::as_tibble_col(column_name = "value") |>
    tidyr::separate_wider_delim("value", names = c("category", "extra", "var", "count", "pct"), delim = ",", too_few = "align_start") |>
    dplyr::mutate(
      count = as.numeric(.data$count),
      pct = round(as.numeric(.data$pct), 2),
      var = dplyr::recode(.data$var, !!!abbrev_nm)
    ) |>
    dplyr::select("var", "count", "pct")
  dirty_names_cleaned(d1$var, abbrev_nm, x)
  d1 |>
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
#' @examples
#' \dontrun{
#' dragen_vc_metrics_read(x)
#' }
#' @export
dragen_vc_metrics_read <- function(x) {
  abbrev_nm1 <- tibble::tribble(
    ~raw, ~clean, ~region,
    "Total", "var_tot", FALSE,
    "Biallelic", "var_biallelic", FALSE,
    "Multiallelic", "var_multiallelic", FALSE,
    "SNPs", "var_snp", FALSE,
    "Insertions (Hom)", "var_ins_hom", FALSE,
    "Insertions (Het)", "var_ins_het", FALSE,
    "Deletions (Hom)", "var_del_hom", FALSE,
    "Deletions (Het)", "var_del_het", FALSE,
    "Indels (Het)", "var_indel_het", FALSE,
    "Chr X number of SNPs over ", "var_snp_x_over_", TRUE,
    "Chr Y number of SNPs over ", "var_snp_y_over_", TRUE,
    "(Chr X SNPs)/(chr Y SNPs) ratio over ", "var_x_over_y_snp_ratio_over_", TRUE,
    "SNP Transitions", "var_snp_transitions", FALSE,
    "SNP Transversions", "var_snp_transversions", FALSE,
    "Ti/Tv ratio", "var_ti_tv_ratio", FALSE,
    "Heterozygous", "var_heterozygous", FALSE,
    "Homozygous", "var_homozygous", FALSE,
    "Het/Hom ratio", "var_het_hom_ratio", FALSE,
    "In dbSNP", "var_in_dbsnp", FALSE,
    "Not in dbSNP", "var_nin_dbsnp", FALSE,
    "Percent Callability", "callability_pct", FALSE,
    "Percent Autosome Callability", "callability_auto_pct", FALSE,
    "Number of samples", "sample_num", FALSE,
    "Reads Processed", "reads_processed", FALSE,
    "Child Sample", "sample_child", FALSE,
    "Percent QC Region Callability in Region 1", "qc_region_callability_pct_region1", FALSE,
    "Percent QC Region Callability in Region 2", "qc_region_callability_pct_region2", FALSE
  )
  d0 <- readr::read_lines(x)
  assertthat::assert_that(grepl("VARIANT CALLER", d0[1]))
  # first detect if this is genome or target region
  d1 <- d0 |>
    tibble::as_tibble_col(column_name = "value") |>
    tidyr::separate_wider_delim(
      "value",
      names = c("category", "sample", "var", "count", "pct"),
      delim = ",", too_few = "align_start"
    )
  reg1 <- NULL
  str1 <- NULL
  tmp <- d1 |>
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
      clean = ifelse(.data$region, glue("{.data$clean}{reg1}"), .data$clean)
    ) |>
    dplyr::select("raw", "clean") |>
    tibble::deframe()

  d <- d1 |>
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
  dirty_names_cleaned(d$var, abbrev_nm, x)
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
#' @examples
#' \dontrun{
#' dragen_mapping_metrics_read(x)
#' }
#' @export
dragen_mapping_metrics_read <- function(x) {
  abbrev_nm <- c(
    "Total input reads" = "reads_tot_input",
    "Total reads removed by downsampling" = "reads_removed_downsamp",
    "Number of duplicate marked reads" = "reads_num_dupmarked",
    "Number of duplicate marked and mate reads removed" = "reads_num_dupmarked_mate_reads_removed",
    "Number of unique reads (excl. duplicate marked reads)" = "reads_num_uniq",
    "Reads with mate sequenced" = "reads_w_mate_seq",
    "Reads without mate sequenced" = "reads_wo_mate_seq",
    "QC-failed reads" = "reads_qcfail",
    "Mapped reads" = "reads_mapped",
    "Mapped reads adjusted for filtered mapping" = "reads_mapped_adjfilt",
    "Mapped reads R1" = "reads_mapped_r1",
    "Mapped reads R2" = "reads_mapped_r2",
    "Number of unique & mapped reads (excl. duplicate marked reads)" = "reads_num_uniq_mapped",
    "Unmapped reads" = "reads_unmapped",
    "Unmapped reads adjusted for filtered mapping" = "reads_unmapped_adjfilt",
    "Adjustment of reads matching non-reference decoys" = "reads_match_nonref_decoys_adj",
    "Adjustment of reads matching exclude contigs" = "reads_match_excl_contigs",
    "Singleton reads (itself mapped; mate unmapped)" = "reads_singleton",
    "Paired reads (itself & mate mapped)" = "reads_paired",
    "Properly paired reads" = "reads_paired_proper",
    "Not properly paired reads (discordant)" = "reads_discordant",
    "Paired reads mapped to different chromosomes" = "reads_paired_mapped_diff_chrom",
    "Paired reads mapped to different chromosomes (MAPQ>=10)" = "reads_paired_mapped_diff_chrom_mapq10",
    "Reads with MAPQ [40:inf)" = "reads_mapq_40_inf",
    "Reads with MAPQ [30:40)" = "reads_mapq_30_40",
    "Reads with MAPQ [20:30)" = "reads_mapq_20_30",
    "Reads with MAPQ [10:20)" = "reads_mapq_10_20",
    "Reads with MAPQ [ 0:10)" = "reads_mapq_0_10",
    "Reads with MAPQ NA (Unmapped reads)" = "reads_mapq_na_unmapped",
    "Reads with indel R1" = "reads_indel_r1",
    "Reads with indel R2" = "reads_indel_r2",
    "Total bases" = "bases_tot",
    "Total bases removed by downsampling" = "bases_removed_downsamp",
    "Total bases R1" = "bases_tot_r1",
    "Total bases R2" = "bases_tot_r2",
    "Mapped bases" = "bases_mapped",
    "Mapped bases R1" = "bases_mapped_r1",
    "Mapped bases R2" = "bases_mapped_r2",
    "Soft-clipped bases" = "bases_softclip",
    "Soft-clipped bases R1" = "bases_softclip_r1",
    "Soft-clipped bases R2" = "bases_softclip_r2",
    "Hard-clipped bases" = "bases_hardclip",
    "Hard-clipped bases R1" = "bases_hardclip_r1",
    "Hard-clipped bases R2" = "bases_hardclip_r2",
    "Mismatched bases R1" = "bases_mismatched_r1",
    "Mismatched bases R2" = "bases_mismatched_r2",
    "Mismatched bases R1 (excl. indels)" = "bases_mismatched_r1_noindels",
    "Mismatched bases R2 (excl. indels)" = "bases_mismatched_r2_noindels",
    "Q30 bases" = "bases_q30",
    "Q30 bases R1" = "bases_q30_r1",
    "Q30 bases R2" = "bases_q30_r2",
    "Q30 bases (excl. dups & clipped bases)" = "bases_q30_nodups_noclipped",
    "Total alignments" = "alignments_tot",
    "Secondary alignments" = "alignments_secondary",
    "Supplementary (chimeric) alignments" = "alignments_chimeric",
    "Estimated read length" = "read_len",
    "Bases in reference genome" = "bases_in_ref_genome",
    "Bases in target bed [% of genome]" = "bases_in_target_bed_genome_pct",
    "Insert length: mean" = "insert_len_mean",
    "Insert length: median" = "insert_len_median",
    "Insert length: standard deviation" = "insert_len_std_dev",
    "Provided sex chromosome ploidy" = "ploidy_sex_chrom_provided",
    "Estimated sample contamination" = "contamination_est",
    "Estimated sample contamination standard error" = "contamination_stderr_est",
    "DRAGEN mapping rate [mil. reads/second]" = "mapping_rate_dragen_milreads_per_sec",
    "Total reads in RG" = "reads_tot_rg",
    "Mapped reads adjusted for excluded mapping" = "reads_mapped_adjexcl",
    "Mapped reads adjusted for filtered and excluded mapping" = "reads_mapped_adjfiltexcl",
    "Unmapped reads adjusted for excluded mapping" = "reads_unmapped_adjexcl",
    "Unmapped reads adjusted for filtered and excluded mapping" = "reads_unmapped_adjfiltexcl",
    "Reads mapping to multiple locations" = "reads_map_multiloc",
    "Adjustment of reads matching filter contigs" = "reads_match_filt_contig_adj",
    "Reads with splice junction" = "reads_splicejunc",
    "Average sequenced coverage over genome" = "cov_avg_seq_over_genome",
    "Filtered rRNA reads" = "reads_rrna_filtered",
    "Mitochondrial reads excluded" = "reads_mito_excl"
  )
  d0 <- readr::read_lines(x)
  assertthat::assert_that(grepl("MAPPING/ALIGNING", d0[1]))
  # File is separated into two sections, the SUMMARY and the PER RG.
  # Based on what I've seen so far, we can have single samples (where
  # the first column just has MAPPING/ALIGNING) or TUMOR/NORMAL samples (where
  # the first column will have a TUMOR or NORMAL prefix).
  reg1 <- paste0("MAPPING/ALIGNING ", c("SUMMARY", "PER RG"), collapse = "|")
  d <- d0 |>
    tibble::as_tibble_col(column_name = "value") |>
    tidyr::separate_wider_delim(
      "value",
      names = c("dragen_sample", "RG", "var", "count", "pct"),
      delim = ",", too_few = "align_start"
    ) |>
    dplyr::mutate(
      count = dplyr::na_if(.data$count, "NA"),
      count = as.numeric(.data$count),
      pct = as.numeric(.data$pct),
      var = dplyr::recode(.data$var, !!!abbrev_nm),
      RG = dplyr::if_else(.data$RG == "", "Total", .data$RG),
      dragen_sample = sub(reg1, "", .data$dragen_sample) |> trimws(),
      dragen_sample = dplyr::if_else(.data$dragen_sample == "", "SINGLE", .data$dragen_sample)
    ) |>
    dplyr::select("dragen_sample", "RG", "var", "count", "pct")
  dirty_names_cleaned(unique(d$var), abbrev_nm, x)
  # pivot
  d |>
    tidyr::pivot_longer(c("count", "pct")) |>
    dplyr::mutate(
      name = dplyr::if_else(.data$name == "count", "", "_pct"),
      var = glue("{.data$var}{.data$name}")
    ) |>
    dplyr::select("dragen_sample", "RG", "var", "value") |>
    dplyr::filter(!is.na(.data$value)) |>
    tidyr::pivot_wider(names_from = "var", values_from = "value")
}

#' Read DRAGEN Coverage Metrics
#'
#' Reads the `coverage_metrics.csv` file generated by DRAGEN.
#' @param x Path to file.
#'
#' @return Tibble with metrics.
#' @examples
#' \dontrun{
#' dragen_coverage_metrics_read(x)
#' }
#' @export
dragen_coverage_metrics_read <- function(x) {
  # all rows except 'Aligned bases' and 'Aligned reads' refer to the region
  abbrev_nm <- tibble::tribble(
    ~raw, ~clean, ~region,
    "Aligned bases", "bases_aligned_tot", FALSE,
    "Aligned reads", "reads_aligned_tot", FALSE,
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
  d0 <- readr::read_lines(x)
  assertthat::assert_that(grepl("COVERAGE SUMMARY", d0[1]))
  # first detect if this is genome, QC coverage region, or target region
  d1 <- d0 |>
    tibble::as_tibble_col(column_name = "value") |>
    tidyr::separate_wider_delim(
      "value",
      delim = ",", too_few = "align_start",
      names = c("category", "dummy1", "var", "value", "pct")
    )
  reg1 <- NULL
  str1 <- NULL
  tmp <- d1 |>
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
      clean = ifelse(.data$region, glue("{.data$clean}{reg1}"), .data$clean)
    ) |>
    dplyr::select("raw", "clean") |>
    tibble::deframe()
  # split to rename the
  # "PCT of genome with coverage [100x: inf)" values
  pat <- glue("PCT of {str1} with coverage ")
  res1 <- d1 |>
    # pct just shows % for a couple rows which can be
    # calculated from their above values
    dplyr::filter(!grepl(pat, .data$var)) |>
    dplyr::select("var", "value")
  dirty_names_cleaned(res1$var, names(abbrev_nm), x)
  res2 <- d1 |>
    dplyr::filter(grepl(pat, .data$var)) |>
    dplyr::mutate(
      var = sub(pat, "", .data$var),
      var = gsub("\\[|\\]|\\(|\\)| ", "", .data$var),
      var = gsub("x", "", .data$var),
      var = gsub("inf", "Inf", .data$var)
    ) |>
    tidyr::separate_wider_delim("var", names = c("start", "end"), delim = ":") |>
    dplyr::mutate(var = as.character(glue("cov_pct_{start}_{end}_{reg1}"))) |>
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

#' Plot DRAGEN Contig Mean Coverage
#' TODO
#' Plots the `wgs_contig_mean_cov_<phenotype>.csv` files.
#'
#' @param d Parsed tibble.
#' @param top_alt_n Number of top covered alt contigs to plot per phenotype.
#' @return A ggplot2 object with chromosomes on X axis, and coverage on Y axis.
dragen_contig_mean_coverage_plot <- function(d, top_alt_n = 15) {
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
    dplyr::filter(.data$coverage > min_cvg) |>
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

#' Read DRAGEN Ploidy Estimation Metrics
#'
#' Reads the `ploidy_estimation_metrics.csv` file generated by DRAGEN.
#' @param x Path to file.
#'
#' @return Tibble with metrics.
dragen_ploidy_estimation_metrics_read <- function(x) {
  raw <- readr::read_lines(x)
  assertthat::assert_that(grepl("PLOIDY ESTIMATION", raw[1]))
  fun1 <- function(x) {
    setNames(
      as.character(glue("cov_{x}_div_auto_median")),
      as.character(glue("{x} median / Autosomal median"))
    )
  }
  abbrev_nm <- c(
    "Autosomal median coverage" = "cov_auto_median",
    "X median coverage" = "cov_x_median",
    "Y median coverage" = "cov_y_median",
    "Ploidy estimation" = "ploidy_est",
    fun1(c(1:22, "X", "Y"))
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
  cols1 <- colnames(d)[colnames(d) != "ploidy_est"]
  d |>
    dplyr::mutate(dplyr::across(dplyr::all_of(cols1), as.numeric))
}
