#' Dracarys Tidy MultiQC
#'
#' Generate tidier representations of MultiQC JSON output
#' @param json Path to `multiqc_data.json`.
#' @param prefix Prefix for output files.
#' @param outdir Path to output results.
#' @return Generates TSV and Parquet representations of the input
#' MultiQC JSON file.
#' @export
dracarys_tidy_multiqc <- function(json, prefix, outdir) {
  e <- emojifont::emoji
  cli::cli_div(theme = list(
    span.file = list(color = "lightblue"),
    span.emph = list(color = "orange")
  ))
  cli::cli_alert_info("{date_log()} {e('dragon')} Start tidying {.file {json}} {e('fire')}")
  # main dracarys function
  d1 <- multiqc_tidy_json(json)
  ## d2 <- select_column_subset_alignmentqc(d1)
  tsv_out <- file.path(outdir, glue::glue("{prefix}.tsv"))
  parquet_out <- file.path(outdir, glue::glue("{prefix}.parquet"))
  readr::write_tsv(d1, tsv_out)
  arrow::write_parquet(d1, parquet_out)

  cli::cli_alert_success("{date_log()} {e('rocket')} End tidying {.file {json}} {e('comet')}!")
  cli::cli_alert_info("{date_log()} {e('tada')} Path to output directory with results for {.emph {prefix}}: {.file {outdir}}")
}

#' Tidy MultiQC JSON
#'
#' Tidies 'multiqc_data.json' output from MultiQC.
#' Modified from the awesome <https://github.com/multimeric/TidyMultiqc>.
#' @param j Path to `multiqc_data.json`.
#' @return A tidy tibble where each column corresponds to a single metric,
#' and each row corresponds to a single sample.
#' @export
multiqc_tidy_json <- function(j) {
  p <- RJSONIO::fromJSON(j)
  dragen_workflows <- c("alignment", "transcriptome", "somatic", "umccrise")
  config_creation_date <- sub(", ", "_", p$config_creation_date)
  config_title <- p$config_title
  get_workflow <- function(config_title) {
    if (grepl("^SBJ[[:digit:]]+$", config_title)) {
      return("dragen_umccrise")
    } else if (grepl("^UMCCR MultiQC Dragen", config_title)) {
      w <- tolower(sub("UMCCR MultiQC Dragen (.*) Report for .*", "\\1", config_title))
      assertthat::assert_that(w %in% dragen_workflows)
      return(paste0("dragen_", w))
    } else if (grepl("[[:digit:]]+.*-merged$", config_title)) {
      return("bcbio_umccrise")
    } else {
      warning(glue::glue(
        "config_title: '{config_title}'.\n",
        "Unknown which umccr workflow this MultiQC JSON was generated from",
      ))
      return("unknown")
    }
  }
  workflow <- get_workflow(config_title)
  d <- dracarys::multiqc_parse_gen(p)
  if (workflow == "dragen_umccrise") {
    # replace the "NA" strings with NA, else we get a column class error
    # due to trying to bind string ('NA') with numeric.
    # https://stackoverflow.com/questions/35292757/replace-values-in-list
    d <- rapply(d, function(x) ifelse(x == "NA", NA, x), how = "replace")
  }
  d <- d |>
    dplyr::bind_rows(.id = "umccr_id") |>
    dplyr::mutate(
      config_creation_date = config_creation_date,
      umccr_workflow = workflow
    ) |>
    dplyr::select(.data$umccr_id, .data$umccr_workflow, .data$config_creation_date, dplyr::everything())

  if (workflow == "dragen_transcriptome") {
    # discard Ref_X control samples
    wts_ref_samples <- paste0("Ref_", 1:6)
    d <- d |>
      dplyr::filter(!.data$umccr_id %in% wts_ref_samples)
  } else if (workflow %in% c("dragen_umccrise", "bcbio_umccrise")) {
    # discard Alice, Bob etc. control samples
    um_ref_samples <- c("Alice", "Bob", "Chen", "Elon", "Dakota")
    um_ref_samples <- paste0(um_ref_samples, rep(c("_T", "_B", ""), each = length(um_ref_samples)))
    d <- d |>
      dplyr::filter(!.data$umccr_id %in% um_ref_samples)
  }
  d
}

#' Parse MultiQC 'report_general_stats_data' JSON Element
#'
#' Parses MultiQC 'report_general_stats_data' JSON Element. Modified from the
#' awesome <https://github.com/multimeric/TidyMultiqc>.
#' @param p Parsed MultiQC JSON.
#' @return A list.
#' @export
multiqc_parse_gen <- function(p) {
  el <- "report_general_stats_data"
  assertthat::assert_that(inherits(p, "list"), el %in% names(p))
  p[[el]] |> purrr::reduce(~ purrr::list_merge(.x, !!!.y))
}

#' Parse MultiQC 'report_saved_raw_data' JSON Element
#'
#' Parses MultiQC 'report_saved_raw_data' JSON Element. Modified from the
#' awesome <https://github.com/multimeric/TidyMultiqc>.
#' @param p Parsed MultiQC JSON.
#' @param tools2delete Character vector of tools to delete from the parsed JSON.
#' @return A list.
#' @export
multiqc_parse_raw <- function(p, tools2delete = NULL) {
  el <- "report_saved_raw_data"
  assertthat::assert_that(inherits(p, "list"), el %in% names(p))
  tool_nms <- names(p[[el]])
  if (!is.null(tools2delete)) {
    assertthat::assert_that(
      is.vector(tools2delete), !is.list(tools2delete),
      all(tools2delete %in% tool_nms)
    )
    # remove given elements from list
    p[[el]] <- base::within(p[[el]], rm(list = tools2delete))
  }
  # For each tool
  p[[el]] |>
    purrr::imap(function(samples, tool) {
      # Remove the superflous "multiqc_" from the start of the tool name
      # tool <- sub("multiqc_", "", tool)

      # For each sample
      samples |> multiqc_kv_map(function(metrics, sample) {
        # For each metric in the above tool
        list(
          key = sample,
          value = metrics |> multiqc_kv_map(function(mvalue, mname) {
            # Sanitise metric names
            combined_metric <- list(
              # key = paste0(tool, ".", mname),
              key = mname,
              value = mvalue
            )
          })
        )
      })
    }) |>
    purrr::reduce(utils::modifyList)
}

multiqc_kv_map <- function(l, func) {
  mapped <- purrr::imap(l, func) |>
    purrr::set_names(nm = NULL)
  keys <- mapped |> purrr::map_chr("key")
  vals <- mapped |> purrr::map("value")
  vals |> purrr::set_names(keys)
}

MULTIQC_COLUMNS <- tibble::tribble(
  ~workflow, ~raw, ~clean,
  "dragen_alignment_qc", "umccr_id", "umccr_id",
  "dragen_alignment_qc", "umccr_workflow", "umccr_workflow",
  "dragen_alignment_qc", "config_creation_date", "date",
  "dragen_alignment_qc", "Total input reads", "reads_tot_input",
  "dragen_alignment_qc", "Total input reads pct", "reads_tot_input_pct",
  "dragen_alignment_qc", "Number of duplicate marked reads", "reads_num_dupmarked",
  "dragen_alignment_qc", "Number of duplicate marked reads pct", "reads_num_dupmarked_pct",
  "dragen_alignment_qc", "Number of unique reads (excl. duplicate marked reads) pct", "reads_num_unique_pct",
  "dragen_alignment_qc", "Reads with mate sequenced", "reads_w_mate_seq",
  "dragen_alignment_qc", "Reads with mate sequenced pct", "reads_w_mate_seq_pct",
  "dragen_alignment_qc", "Reads without mate sequenced", "reads_wo_mate_seq",
  "dragen_alignment_qc", "Reads without mate sequenced pct", "reads_wo_mate_seq_pct",
  "dragen_alignment_qc", "QC-failed reads", "reads_qcfail",
  "dragen_alignment_qc", "QC-failed reads pct", "reads_qcfail_pct",
  "dragen_alignment_qc", "Mapped reads", "reads_mapped",
  "dragen_alignment_qc", "Mapped reads pct", "reads_mapped_pct",
  "dragen_alignment_qc", "Mapped reads adjusted for filtered mapping", "reads_mapped_adjfilt",
  "dragen_alignment_qc", "Mapped reads adjusted for filtered mapping pct", "reads_mapped_adjfilt_pct",
  "dragen_alignment_qc", "Mapped reads R1", "reads_mapped_r1",
  "dragen_alignment_qc", "Mapped reads R1 pct", "reads_mapped_r1_pct",
  "dragen_alignment_qc", "Mapped reads R2", "reads_mapped_r2",
  "dragen_alignment_qc", "Mapped reads R2 pct", "reads_mapped_r2_pct",
  "dragen_alignment_qc", "Number of unique & mapped reads (excl. duplicate marked reads)", "reads_num_uniq_mapped",
  "dragen_alignment_qc", "Number of unique & mapped reads (excl. duplicate marked reads) pct", "reads_num_uniq_mapped_pct",
  "dragen_alignment_qc", "Unmapped reads", "reads_unmapped",
  "dragen_alignment_qc", "Unmapped reads pct", "reads_unmapped_pct",
  "dragen_alignment_qc", "Unmapped reads adjusted for filtered mapping", "reads_unmapped_adjfilt",
  "dragen_alignment_qc", "Unmapped reads adjusted for filtered mapping pct", "reads_unmapped_adjfilt_pct",
  "dragen_alignment_qc", "Adjustment of reads matching non-reference decoys", "reads_match_nonref_decoys_adj",
  "dragen_alignment_qc", "Adjustment of reads matching non-reference decoys pct", "reads_match_nonref_decoys_adj_pct",
  "dragen_alignment_qc", "Singleton reads (itself mapped; mate unmapped)", "reads_singleton",
  "dragen_alignment_qc", "Singleton reads (itself mapped; mate unmapped) pct", "reads_singleton_pct",
  "dragen_alignment_qc", "Paired reads (itself & mate mapped)", "reads_paired",
  "dragen_alignment_qc", "Paired reads (itself & mate mapped) pct", "reads_paired_pct",
  "dragen_alignment_qc", "Properly paired reads", "reads_paired_proper",
  "dragen_alignment_qc", "Properly paired reads pct", "reads_paired_proper_pct",
  "dragen_alignment_qc", "Not properly paired reads (discordant)", "reads_discordant",
  "dragen_alignment_qc", "Not properly paired reads (discordant) pct", "reads_discordant_pct",
  "dragen_alignment_qc", "Paired reads mapped to different chromosomes", "reads_paired_mapped_diff_chrom",
  "dragen_alignment_qc", "Paired reads mapped to different chromosomes pct", "reads_paired_mapped_diff_chrom_pct",
  "dragen_alignment_qc", "Paired reads mapped to different chromosomes (MAPQ>=10)", "reads_paired_mapped_diff_chrom_mapq10",
  "dragen_alignment_qc", "Paired reads mapped to different chromosomes (MAPQ>=10) pct", "reads_paired_mapped_diff_chrom_mapq10_pct",
  "dragen_alignment_qc", "Reads with MAPQ [40:inf)", "reads_mapq_40_inf",
  "dragen_alignment_qc", "Reads with MAPQ [40:inf) pct", "reads_mapq_40_inf_pct",
  "dragen_alignment_qc", "Reads with MAPQ [30:40)", "reads_mapq_30_40",
  "dragen_alignment_qc", "Reads with MAPQ [30:40) pct", "reads_mapq_30_40_pct",
  "dragen_alignment_qc", "Reads with MAPQ [20:30)", "reads_mapq_20_30",
  "dragen_alignment_qc", "Reads with MAPQ [20:30) pct", "reads_mapq_20_30_pct",
  "dragen_alignment_qc", "Reads with MAPQ [10:20)", "reads_mapq_10_20",
  "dragen_alignment_qc", "Reads with MAPQ [10:20) pct", "reads_mapq_10_20_pct",
  "dragen_alignment_qc", "Reads with MAPQ [ 0:10)", "reads_mapq_0_10",
  "dragen_alignment_qc", "Reads with MAPQ [ 0:10) pct", "reads_mapq_0_10_pct",
  "dragen_alignment_qc", "Reads with MAPQ NA (Unmapped reads)", "reads_mapq_NA_unmapped",
  "dragen_alignment_qc", "Reads with MAPQ NA (Unmapped reads) pct", "reads_mapq_NA_unmapped_pct",
  "dragen_alignment_qc", "Reads with indel R1", "reads_indel_r1",
  "dragen_alignment_qc", "Reads with indel R1 pct", "reads_indel_r1_pct",
  "dragen_alignment_qc", "Reads with indel R2", "reads_indel_r2",
  "dragen_alignment_qc", "Reads with indel R2 pct", "reads_indel_r2_pct",
  "dragen_alignment_qc", "Total bases", "bases_tot",
  "dragen_alignment_qc", "Total bases R1", "bases_tot_r1",
  "dragen_alignment_qc", "Total bases R2", "bases_tot_r2",
  "dragen_alignment_qc", "Mapped bases R1", "bases_mapped_r1",
  "dragen_alignment_qc", "Mapped bases R1 pct", "bases_mapped_r1_pct",
  "dragen_alignment_qc", "Mapped bases R2", "bases_mapped_r2",
  "dragen_alignment_qc", "Mapped bases R2 pct", "bases_mapped_r2_pct",
  "dragen_alignment_qc", "Soft-clipped bases R1", "bases_softclip_r1",
  "dragen_alignment_qc", "Soft-clipped bases R1 pct", "bases_softclip_r1_pct",
  "dragen_alignment_qc", "Soft-clipped bases R2 pct", "bases_softclip_r2_pct",
  "dragen_alignment_qc", "Mismatched bases R1", "bases_mismatched_r1",
  "dragen_alignment_qc", "Mismatched bases R1 pct", "bases_mismatched_r1_pct",
  "dragen_alignment_qc", "Mismatched bases R2 pct", "bases_mismatched_r2_pct",
  "dragen_alignment_qc", "Mismatched bases R1 (excl. indels)", "bases_mismatched_r1_noindels",
  "dragen_alignment_qc", "Mismatched bases R1 (excl. indels) pct", "bases_mismatched_r1_noindels_pct",
  "dragen_alignment_qc", "Mismatched bases R2 (excl. indels) pct", "bases_mismatched_r2_noindels_pct",
  "dragen_alignment_qc", "Q30 bases", "bases_q30",
  "dragen_alignment_qc", "Q30 bases pct", "bases_q30_pct",
  "dragen_alignment_qc", "Q30 bases R1", "bases_q30_r1",
  "dragen_alignment_qc", "Q30 bases R1 pct", "bases_q30_r1_pct",
  "dragen_alignment_qc", "Q30 bases R2", "bases_q30_r2",
  "dragen_alignment_qc", "Q30 bases R2 pct", "bases_q30_r2_pct",
  "dragen_alignment_qc", "Q30 bases (excl. dups & clipped bases)", "bases_q30_nodups_noclipped",
  "dragen_alignment_qc", "Q30 bases (excl. dups & clipped bases) pct", "bases_q30_nodups_noclipped_pct",
  "dragen_alignment_qc", "Total alignments", "alignments_tot",
  "dragen_alignment_qc", "Secondary alignments", "alignments_secondary",
  "dragen_alignment_qc", "Secondary alignments pct", "alignments_secondary_pct",
  "dragen_alignment_qc", "Supplementary (chimeric) alignments", "alignments_chimeric",
  "dragen_alignment_qc", "Estimated read length", "read_len",
  "dragen_alignment_qc", "Bases in reference genome", "bases_in_ref_genome",
  "dragen_alignment_qc", "Bases in target bed [% of genome]", "bases_in_target_bed_genome_pct",
  "dragen_alignment_qc", "Average sequenced coverage over genome", "cov_avg_seq_over_genome",
  "dragen_alignment_qc", "Insert length: mean", "insert_len_mean",
  "dragen_alignment_qc", "Insert length: median", "insert_len_median",
  "dragen_alignment_qc", "Insert length: standard deviation", "insert_len_std_dev",
  "dragen_alignment_qc", "Provided sex chromosome ploidy", "ploidy_sex_chrom_provided",
  "dragen_alignment_qc", "Estimated sample contamination", "contamination_est",
  "dragen_alignment_qc", "DRAGEN mapping rate [mil. reads/second]", "mapping_rate_dragen_milreads_per_sec",
  "dragen_alignment_qc", "Autosomal median coverage", "cov_auto_median",
  "dragen_alignment_qc", "X median coverage", "cov_x_median",
  "dragen_alignment_qc", "Y median coverage", "cov_y_median",
  "dragen_alignment_qc", "1 median / Autosomal median", "cov_1_div_auto_medians",
  "dragen_alignment_qc", "2 median / Autosomal median", "cov_2_div_auto_medians",
  "dragen_alignment_qc", "3 median / Autosomal median", "cov_3_div_auto_medians",
  "dragen_alignment_qc", "4 median / Autosomal median", "cov_4_div_auto_medians",
  "dragen_alignment_qc", "5 median / Autosomal median", "cov_5_div_auto_medians",
  "dragen_alignment_qc", "6 median / Autosomal median", "cov_6_div_auto_medians",
  "dragen_alignment_qc", "7 median / Autosomal median", "cov_7_div_auto_medians",
  "dragen_alignment_qc", "8 median / Autosomal median", "cov_8_div_auto_medians",
  "dragen_alignment_qc", "9 median / Autosomal median", "cov_9_div_auto_medians",
  "dragen_alignment_qc", "10 median / Autosomal median", "cov_10_div_auto_median",
  "dragen_alignment_qc", "11 median / Autosomal median", "cov_11_div_auto_median",
  "dragen_alignment_qc", "12 median / Autosomal median", "cov_12_div_auto_median",
  "dragen_alignment_qc", "13 median / Autosomal median", "cov_13_div_auto_median",
  "dragen_alignment_qc", "14 median / Autosomal median", "cov_14_div_auto_median",
  "dragen_alignment_qc", "15 median / Autosomal median", "cov_15_div_auto_median",
  "dragen_alignment_qc", "16 median / Autosomal median", "cov_16_div_auto_median",
  "dragen_alignment_qc", "17 median / Autosomal median", "cov_17_div_auto_median",
  "dragen_alignment_qc", "18 median / Autosomal median", "cov_18_div_auto_median",
  "dragen_alignment_qc", "19 median / Autosomal median", "cov_19_div_auto_median",
  "dragen_alignment_qc", "20 median / Autosomal median", "cov_20_div_auto_median",
  "dragen_alignment_qc", "21 median / Autosomal median", "cov_21_div_auto_median",
  "dragen_alignment_qc", "22 median / Autosomal median", "cov_22_div_auto_median",
  "dragen_alignment_qc", "X median / Autosomal median", "cov_x_div_auto_median",
  "dragen_alignment_qc", "Y median / Autosomal median", "cov_y_div_auto_median",
  "dragen_alignment_qc", "Ploidy estimation", "ploidy_est",
  "dragen_alignment_qc", "Aligned bases", "bases_aligned",
  "dragen_alignment_qc", "Aligned bases in genome", "bases_aligned_in_genome",
  "dragen_alignment_qc", "Aligned bases in genome pct", "bases_aligned_in_genome_pct",
  "dragen_alignment_qc", "Average alignment coverage over genome", "cov_alignment_avg_over_genome",
  "dragen_alignment_qc", "Uniformity of coverage (PCT > 0.2*mean) over genome", "cov_uniformity_over_genome_pct_gt02mean",
  "dragen_alignment_qc", "PCT of genome with coverage [100x: inf)", "cov_genome_pct_100x_inf",
  "dragen_alignment_qc", "PCT of genome with coverage [ 50x: inf)", "cov_genome_pct_50x_inf",
  "dragen_alignment_qc", "PCT of genome with coverage [ 20x: inf)", "cov_genome_pct_20x_inf",
  "dragen_alignment_qc", "PCT of genome with coverage [ 15x: inf)", "cov_genome_pct_15x_inf",
  "dragen_alignment_qc", "PCT of genome with coverage [ 10x: inf)", "cov_genome_pct_10x_inf",
  "dragen_alignment_qc", "PCT of genome with coverage [  3x: inf)", "cov_genome_pct_3x_inf",
  "dragen_alignment_qc", "PCT of genome with coverage [  1x: inf)", "cov_genome_pct_1x_inf",
  "dragen_alignment_qc", "PCT of genome with coverage [  0x: inf)", "cov_genome_pct_0x_inf",
  "dragen_alignment_qc", "PCT of genome with coverage [ 50x:100x)", "cov_genome_pct_50x_100x",
  "dragen_alignment_qc", "PCT of genome with coverage [ 20x: 50x)", "cov_genome_pct_20x_50x",
  "dragen_alignment_qc", "PCT of genome with coverage [ 15x: 20x)", "cov_genome_pct_15x_20x",
  "dragen_alignment_qc", "PCT of genome with coverage [ 10x: 15x)", "cov_genome_pct_10x_15x",
  "dragen_alignment_qc", "PCT of genome with coverage [  3x: 10x)", "cov_genome_pct_3x_10x",
  "dragen_alignment_qc", "PCT of genome with coverage [  1x:  3x)", "cov_genome_pct_1x_3x",
  "dragen_alignment_qc", "PCT of genome with coverage [  0x:  1x)", "cov_genome_pct_0x_1x",
  "dragen_alignment_qc", "Average chr X coverage over genome", "cov_avg_x_over_genome",
  "dragen_alignment_qc", "Average chr Y coverage over genome", "cov_avg_y_over_genome",
  "dragen_alignment_qc", "Average mitochondrial coverage over genome", "cov_avg_mt_over_genome",
  "dragen_alignment_qc", "Average autosomal coverage over genome", "cov_avg_auto_over_genome",
  "dragen_alignment_qc", "Median autosomal coverage over genome", "cov_median_auto_over_genome",
  "dragen_alignment_qc", "Mean/Median autosomal coverage ratio over genome", "cov_mean_median_auto_ratio_over_genome",
  "dragen_alignment_qc", "Aligned reads", "reads_aligned",
  "dragen_alignment_qc", "Aligned reads in genome", "reads_aligned_in_genome",
  "dragen_alignment_qc", "Aligned reads in genome pct", "reads_aligned_in_genome_pct",
  "dragen_tumor_normal", "Total", "var_tot",
  "dragen_tumor_normal", "Total pct", "var_tot_pct",
  "dragen_tumor_normal", "Biallelic", "var_biallelic",
  "dragen_tumor_normal", "Biallelic pct", "var_biallelic_pct",
  "dragen_tumor_normal", "Multiallelic", "var_multiallelic",
  "dragen_tumor_normal", "Multiallelic pct", "var_multiallelic_pct",
  "dragen_tumor_normal", "SNPs", "var_snp",
  "dragen_tumor_normal", "SNPs pct", "var_snp_pct",
  "dragen_tumor_normal", "Insertions (Hom)", "var_ins_hom",
  "dragen_tumor_normal", "Insertions (Hom) pct", "var_ins_hom_pct",
  "dragen_tumor_normal", "Insertions (Het)", "var_ins_het",
  "dragen_tumor_normal", "Insertions (Het) pct", "var_ins_het_pct",
  "dragen_tumor_normal", "Deletions (Hom)", "var_del_hom",
  "dragen_tumor_normal", "Deletions (Hom) pct", "var_del_hom_pct",
  "dragen_tumor_normal", "Deletions (Het)", "var_del_het",
  "dragen_tumor_normal", "Deletions (Het) pct", "var_del_het_pct",
  "dragen_tumor_normal", "Indels (Het)", "var_indel_het",
  "dragen_tumor_normal", "Indels (Het) pct", "var_indel_het_pct",
  "dragen_tumor_normal", "Chr X number of SNPs over genome", "var_snp_x_over_genome",
  "dragen_tumor_normal", "Chr Y number of SNPs over genome", "var_snp_y_over_genome",
  "dragen_tumor_normal", "(Chr X SNPs)/(chr Y SNPs) ratio over genome", "var_x_over_y_snp_ratio_over_genome",
  "dragen_tumor_normal", "SNP Transitions", "var_snp_transitions",
  "dragen_tumor_normal", "SNP Transversions", "var_snp_transversions",
  "dragen_tumor_normal", "Ti/Tv ratio", "var_ti_tv_ratio",
  "dragen_tumor_normal", "Heterozygous", "var_heterozygous",
  "dragen_tumor_normal", "Homozygous", "var_homozygous",
  "dragen_tumor_normal", "Het/Hom ratio", "var_het_hom_ratio",
  "dragen_tumor_normal", "In dbSNP", "var_in_dbsnp",
  "dragen_tumor_normal", "In dbSNP pct", "var_in_dbsnp_pct",
  "dragen_tumor_normal", "Not in dbSNP", "var_nin_dbsnp",
  "dragen_tumor_normal", "Not in dbSNP pct", "var_nin_dbsnp_pct",
  "dragen_tumor_normal", "Percent Callability", "callability_pct",
  "dragen_tumor_normal", "Percent Autosome Callability", "callability_auto_pct",
  "dragen_tumor_normal", "Insertions", "var_ins_tot",
  "dragen_tumor_normal", "Insertions pct", "var_ins_tot_pct",
  "dragen_tumor_normal", "Deletions", "var_del_tot",
  "dragen_tumor_normal", "Deletions pct", "var_del_tot_pct",
  "dragen_tumor_normal", "Indels", "var_indel_tot",
  "dragen_tumor_normal", "Indels pct", "var_indel_tot_pct",
  "dragen_tumor_normal", "Number of samples", "sample_num",
  "dragen_tumor_normal", "Reads Processed", "reads_processed",
  "dragen_tumor_normal", "Child Sample", "sample_child",
  "dragen_tumor_normal", "Filtered vars", "vars_tot_filt",
  "dragen_tumor_normal", "Filtered vars pct", "vars_tot_filt_pct",
  "dragen_tumor_normal", "Filtered SNPs", "vars_snp_filt",
  "dragen_tumor_normal", "Filtered SNPs pct", "vars_snp_filt_pct",
  "dragen_tumor_normal", "Filtered indels", "vars_indel_filt",
  "dragen_transcriptome", "Adjustment of reads matching filter contigs", "reads_match_filter_contigs_adj",
  "dragen_transcriptome", "Adjustment of reads matching filter contigs pct", "reads_match_filter_contigs_adj_pct",
  "dragen_transcriptome", "Reads with splice junction", "reads_splice_junction",
  "dragen_transcriptome", "Reads with splice junction pct", "reads_splice_junction_pct",
  "dragen_umccrise", "concordance_concordance", "conpair_concordance",
  "dragen_umccrise", "concordance_used_markers", "conpair_concordance_used_markers",
  "dragen_umccrise", "concordance_total_markers", "conpair_concordance_total_markers",
  "dragen_umccrise", "concordance_marker_threshold", "conpair_concordance_marker_threshold",
  "dragen_umccrise", "concordance_min_mapping_quality", "conpair_concordance_min_mapping_quality",
  "dragen_umccrise", "concordance_min_base_quality", "conpair_concordance_min_base_quality",
  "dragen_umccrise", "contamination", "conpair_contamination",
  "dragen_umccrise", "QCStatus", "purple_qc_status",
  "dragen_umccrise", "Method", "purple_method",
  "dragen_umccrise", "CopyNumberSegments", "purple_cn_segs",
  "dragen_umccrise", "UnsupportedCopyNumberSegments", "purple_unsupported_cn_segs",
  "dragen_umccrise", "Purity", "purple_purity",
  "dragen_umccrise", "AmberGender", "amber_gender",
  "dragen_umccrise", "CobaltGender", "cobalt_gender",
  "dragen_umccrise", "DeletedGenes", "purple_deleted_genes",
  "dragen_umccrise", "Contamination", "purple_contamination",
  "dragen_umccrise", "GermlineAberrations", "purple_germline_aberrations",
  "dragen_umccrise", "purity", "purple_purity",
  "dragen_umccrise", "normFactor", "purple_normfactor",
  "dragen_umccrise", "score", "purple_score",
  "dragen_umccrise", "diploidProportion", "purple_diploid_prop",
  "dragen_umccrise", "ploidy", "purple_ploidy",
  "dragen_umccrise", "gender", "purple_gender",
  "dragen_umccrise", "status", "purple_status",
  "dragen_umccrise", "polyclonalProportion", "purple_polyclonal_prop",
  "dragen_umccrise", "minPurity", "purple_min_purity",
  "dragen_umccrise", "maxPurity", "purple_max_purity",
  "dragen_umccrise", "minPloidy", "purple_min_ploidy",
  "dragen_umccrise", "maxPloidy", "purple_max_ploidy",
  "dragen_umccrise", "minDiploidProportion", "purple_min_diploid_prop",
  "dragen_umccrise", "maxDiploidProportion", "purple_max_diploid_prop",
  "dragen_umccrise", "version", "purple_version",
  "dragen_umccrise", "somaticPenalty", "purple_somatic_penalty",
  "dragen_umccrise", "wholeGenomeDuplication", "purple_whole_genome_dup",
  "dragen_umccrise", "msIndelsPerMb", "purple_ms_indels_permb",
  "dragen_umccrise", "msStatus", "purple_ms_status",
  "dragen_umccrise", "tml", "purple_tml",
  "dragen_umccrise", "tmlStatus", "purple_tml_status",
  "dragen_umccrise", "tmbPerMb", "purple_tmb_permb",
  "dragen_umccrise", "tmbStatus", "purple_tmb_status",
  "dragen_umccrise", "svTumorMutationalBurden", "purple_tmb_sv",
  "dragen_umccrise", "1_x_pc", "mosdepth_1x_pc",
  "dragen_umccrise", "5_x_pc", "mosdepth_5x_pc",
  "dragen_umccrise", "10_x_pc", "mosdepth_10x_pc",
  "dragen_umccrise", "30_x_pc", "mosdepth_30x_pc",
  "dragen_umccrise", "50_x_pc", "mosdepth_50x_pc",
  "dragen_umccrise", "median_coverage", "mosdepth_median_cov",
  "dragen_umccrise", "raw_total_sequences", "samtools_raw_total_sequences",
  "dragen_umccrise", "filtered_sequences", "samtools_filtered_sequences",
  "dragen_umccrise", "sequences", "samtools_sequences",
  "dragen_umccrise", "is_sorted", "samtools_is_sorted",
  "dragen_umccrise", "1st_fragments", "samtools_1st_fragments",
  "dragen_umccrise", "last_fragments", "samtools_last_fragments",
  "dragen_umccrise", "reads_mapped", "samtools_reads_mapped",
  "dragen_umccrise", "reads_mapped_and_paired", "samtools_reads_mapped_and_paired",
  "dragen_umccrise", "reads_unmapped", "samtools_reads_unmapped",
  "dragen_umccrise", "reads_properly_paired", "samtools_reads_properly_paired",
  "dragen_umccrise", "reads_paired", "samtools_reads_paired",
  "dragen_umccrise", "reads_duplicated", "reads_num_dupmarked",
  "dragen_umccrise", "reads_duplicated_percent", "reads_num_dupmarked_pct",
  "dragen_umccrise", "reads_MQ0", "samtools_reads_MQ0",
  "dragen_umccrise", "reads_QC_failed", "samtools_reads_QC_failed",
  "dragen_umccrise", "non-primary_alignments", "samtools_non_primary_alignments",
  "dragen_umccrise", "total_length", "samtools_total_length",
  "dragen_umccrise", "total_first_fragment_length", "samtools_total_first_fragment_length",
  "dragen_umccrise", "total_last_fragment_length", "samtools_total_last_fragment_length",
  "dragen_umccrise", "bases_mapped", "samtools_bases_mapped",
  "dragen_umccrise", "bases_mapped_(cigar)", "samtools_bases_mapped_(cigar)",
  "dragen_umccrise", "bases_trimmed", "samtools_bases_trimmed",
  "dragen_umccrise", "bases_duplicated", "samtools_bases_duplicated",
  "dragen_umccrise", "mismatches", "samtools_mismatches",
  "dragen_umccrise", "error_rate", "samtools_error_rate",
  "dragen_umccrise", "average_length", "samtools_average_length",
  "dragen_umccrise", "average_first_fragment_length", "samtools_average_first_fragment_length",
  "dragen_umccrise", "average_last_fragment_length", "samtools_average_last_fragment_length",
  "dragen_umccrise", "maximum_length", "samtools_maximum_length",
  "dragen_umccrise", "maximum_first_fragment_length", "samtools_maximum_first_fragment_length",
  "dragen_umccrise", "maximum_last_fragment_length", "samtools_maximum_last_fragment_length",
  "dragen_umccrise", "average_quality", "samtools_average_quality",
  "dragen_umccrise", "insert_size_average", "samtools_insert_size_average",
  "dragen_umccrise", "insert_size_standard_deviation", "samtools_insert_size_standard_deviation",
  "dragen_umccrise", "inward_oriented_pairs", "samtools_inward_oriented_pairs",
  "dragen_umccrise", "outward_oriented_pairs", "samtools_outward_oriented_pairs",
  "dragen_umccrise", "pairs_with_other_orientation", "samtools_pairs_with_other_orientation",
  "dragen_umccrise", "pairs_on_different_chromosomes", "samtools_pairs_on_different_chromosomes",
  "dragen_umccrise", "percentage_of_properly_paired_reads_(%)", "samtools_reads_properly_paired_pct_round",
  "dragen_umccrise", "reads_properly_paired_percent", "samtools_reads_properly_paired_pct",
  "dragen_umccrise", "reads_mapped_percent", "samtools_reads_mapped_pct",
  "dragen_umccrise", "reads_mapped_and_paired_percent", "samtools_reads_mapped_and_paired_pct",
  "dragen_umccrise", "reads_unmapped_percent", "samtools_reads_unmapped_pct",
  "dragen_umccrise", "reads_paired_percent", "samtools_reads_paired_pct",
  "dragen_umccrise", "reads_MQ0_percent", "samtools_reads_MQ0_pct",
  "dragen_umccrise", "reads_QC_failed_percent", "samtools_reads_QC_failed_pct",
  "dragen_umccrise", "filt_indels", "umccrise_filt_indels",
  "dragen_umccrise", "filt_snps", "umccrise_filt_snps",
  "dragen_umccrise", "filt_vars", "umccrise_filt_vars",
  "dragen_umccrise", "indels", "umccrise_indels",
  "dragen_umccrise", "snps", "umccrise_snps",
  "dragen_umccrise", "viral_content", "oncoviruses_viral_content",
  "dragen_umccrise", "germline", "umccrise_germline",
  "dragen_umccrise", "germline_predispose", "umccrise_germline_predispose",
  "bcbio_umccrise", "concordance", "conpair_concordance",
  "bcbio_umccrise", "concordance_markers", "conpair_concordance_used_over_tot_markers",
  "bcbio_umccrise", "SegmentPass", "purple_segment_pass",
  "bcbio_umccrise", "GenderPass", "purple_gender_pass",
  "bcbio_umccrise", "DeletedGenesPass", "purple_deleted_genes_pass",
  "bcbio_umccrise", "SegmentScore", "purple_segment_score",
  "bcbio_umccrise", "UnsupportedSegments", "purple_unsupported_cn_segs",
  "bcbio_umccrise", "Ploidy", "purple_ploidy",
  "bcbio_umccrise", "%GC", "bcbio_gc_pct",
  "bcbio_umccrise", "Average_insert_size", "insert_len_mean",
  "bcbio_umccrise", "Average_read_length", "read_len",
  "bcbio_umccrise", "Avg_coverage", "bcbio_cov_avg",
  "bcbio_umccrise", "Duplicates", "reads_num_dupmarked",
  "bcbio_umccrise", "Duplicates_pct", "reads_num_dupmarked_pct",
  "bcbio_umccrise", "Mapped_paired_reads", "reads_paired",
  "bcbio_umccrise", "Mapped_reads", "reads_mapped",
  "bcbio_umccrise", "Mapped_reads_pct", "reads_mapped_pct",
  "bcbio_umccrise", "Mapped_unique_reads", "reads_num_uniq_mapped",
  "bcbio_umccrise", "Offtarget_pct", "reads_offtarget_pct",
  "bcbio_umccrise", "Ontarget_pct", "reads_ontarget_pct",
  "bcbio_umccrise", "Ontarget_unique_reads", "reads_ontarget_unique",
  "bcbio_umccrise", "Sequences_flagged_as_poor_quality", "sequences_poor_quality",
  "bcbio_umccrise", "Total_reads", "reads_tot",
  "bcbio_umccrise", "Usable_pct", "reads_usable_pct",
  "bcbio_umccrise", "RiP_pct", "reads_rip_pct",
  "bcbio_umccrise", "non-primary_alignments_percent", "alignments_secondary_pct",
  "bcbio_umccrise", "family_id", "peddy_family_id",
  "bcbio_umccrise", "paternal_id", "peddy_paternal_id",
  "bcbio_umccrise", "maternal_id", "peddy_maternal_id",
  "bcbio_umccrise", "sex", "peddy_sex",
  "bcbio_umccrise", "phenotype", "peddy_phenotype",
  "bcbio_umccrise", "Ethnicity", "peddy_ethnicity",
  "bcbio_umccrise", "duplicates", "peddy_duplicates",
  "bcbio_umccrise", "het_call_rate", "peddy_het_call_rate",
  "bcbio_umccrise", "het_ratio", "peddy_het_ratio",
  "bcbio_umccrise", "het_mean_depth", "peddy_het_mean_depth",
  "bcbio_umccrise", "het_idr_baf", "peddy_het_idr_baf",
  "bcbio_umccrise", "ancestry-prediction", "peddy_ancestry_prediction",
  "bcbio_umccrise", "PC1", "peddy_pc1",
  "bcbio_umccrise", "PC2", "peddy_pc2",
  "bcbio_umccrise", "PC3", "peddy_pc3",
  "bcbio_umccrise", "sex_het_ratio", "peddy_sex_het_ratio",
  "bcbio_umccrise", "sex_fixed", "peddy_sex_fixed",
  "bcbio_umccrise", "depth_outlier_het_check", "peddy_depth_outlier_het_check",
  "bcbio_umccrise", "het_count_het_check", "peddy_het_count_het_check",
  "bcbio_umccrise", "het_ratio_het_check", "peddy_het_ratio_het_check",
  "bcbio_umccrise", "idr_baf_het_check", "peddy_idr_baf_het_check",
  "bcbio_umccrise", "mean_depth_het_check", "peddy_mean_depth_het_check",
  "bcbio_umccrise", "median_depth_het_check", "peddy_median_depth_het_check",
  "bcbio_umccrise", "p10_het_check", "peddy_p10_het_check",
  "bcbio_umccrise", "p90_het_check", "peddy_p90_het_check",
  "bcbio_umccrise", "sampled_sites_het_check", "peddy_sampled_sites_het_check",
  "bcbio_umccrise", "call_rate_het_check", "peddy_call_rate_het_check",
  "bcbio_umccrise", "PC1_het_check", "peddy_pc1_het_check",
  "bcbio_umccrise", "PC2_het_check", "peddy_pc2_het_check",
  "bcbio_umccrise", "PC3_het_check", "peddy_pc3_het_check",
  "bcbio_umccrise", "PC4_het_check", "peddy_pc4_het_check",
  "bcbio_umccrise", "ancestry-prediction_het_check", "peddy_ancestry_prediction_het_check",
  "bcbio_umccrise", "ancestry-prob_het_check", "peddy_ancestry_prob_het_check",
  "bcbio_umccrise", "error_sex_check", "peddy_error_sex_check",
  "bcbio_umccrise", "het_count_sex_check", "peddy_het_count_sex_check",
  "bcbio_umccrise", "het_ratio_sex_check", "peddy_het_ratio_sex_check",
  "bcbio_umccrise", "hom_alt_count_sex_check", "peddy_hom_alt_count_sex_check",
  "bcbio_umccrise", "hom_ref_count_sex_check", "peddy_hom_ref_count_sex_check",
  "bcbio_umccrise", "ped_sex_sex_check", "peddy_ped_sex_sex_check",
  "bcbio_umccrise", "predicted_sex_sex_check", "peddy_predicted_sex_sex_check",
  "bcbio_umccrise", "percent_gc", "gc_pct",
  "bcbio_umccrise", "avg_sequence_length", "avg_sequence_length",
  "bcbio_umccrise", "total_sequences", "sequences_tot",
  "bcbio_umccrise", "percent_duplicates", "reads_num_dupmarked_pct",
  "bcbio_umccrise", "percent_fails", "fastqc_fails_pct",
  "bcbio_umccrise", "Genome", "genome",
  "bcbio_umccrise", "Number_of_variants_before_filter", "snpeff_var_num_before_filter",
  "bcbio_umccrise", "Number_of_known_variants (i.e. non-empty ID)", "snpeff_var_num_known",
  "bcbio_umccrise", "Number_of_known_variants (i.e. non-empty ID)_percent", "snpeff_var_num_known_pct",
  "bcbio_umccrise", "Number_of_effects", "snpeff_effects_num",
  "bcbio_umccrise", "Genome_total_length", "snpeff_genome_tot_len",
  "bcbio_umccrise", "Genome_effective_length", "snpeff_genome_eff_len",
  "bcbio_umccrise", "Change_rate", "snpeff_change_rate",
  "bcbio_umccrise", "HIGH", "snpeff_high",
  "bcbio_umccrise", "HIGH_percent", "snpeff_high_pct",
  "bcbio_umccrise", "LOW", "snpeff_low",
  "bcbio_umccrise", "LOW_percent", "snpeff_low_pct",
  "bcbio_umccrise", "MODERATE", "snpeff_moderate",
  "bcbio_umccrise", "MODERATE_percent", "snpeff_moderate_pct",
  "bcbio_umccrise", "MODIFIER", "snpeff_modifier",
  "bcbio_umccrise", "MODIFIER_percent", "snpeff_modifier_pct",
  "bcbio_umccrise", "MISSENSE", "snpeff_missense",
  "bcbio_umccrise", "MISSENSE_percent", "snpeff_missense_pct",
  "bcbio_umccrise", "NONSENSE", "snpeff_nonsense",
  "bcbio_umccrise", "NONSENSE_percent", "snpeff_nonsense_pct",
  "bcbio_umccrise", "SILENT", "snpeff_silent",
  "bcbio_umccrise", "SILENT_percent", "snpeff_silent_pct",
  "bcbio_umccrise", "Missense_Silent_ratio", "snpeff_missense_silent_ratio",
  "bcbio_umccrise", "Type", "snpeff_type",
  "bcbio_umccrise", "3_prime_UTR_variant", "snpeff_3_prime_utr_variant",
  "bcbio_umccrise", "3_prime_UTR_variant_percent", "snpeff_3_prime_utr_variant_pct",
  "bcbio_umccrise", "5_prime_UTR_premature_start_codon_gain_variant", "snpeff_5_prime_utr_premature_start_codon_gain_variant",
  "bcbio_umccrise", "5_prime_UTR_premature_start_codon_gain_variant_percent", "snpeff_5_prime_utr_premature_start_codon_gain_variant_pct",
  "bcbio_umccrise", "5_prime_UTR_variant", "snpeff_5_prime_utr_variant",
  "bcbio_umccrise", "5_prime_UTR_variant_percent", "snpeff_5_prime_utr_variant_pct",
  "bcbio_umccrise", "TF_binding_site_variant", "snpeff_tf_binding_site_variant",
  "bcbio_umccrise", "TF_binding_site_variant_percent", "snpeff_tf_binding_site_variant_pct",
  "bcbio_umccrise", "conservative_inframe_deletion", "snpeff_conservative_inframe_deletion",
  "bcbio_umccrise", "conservative_inframe_deletion_percent", "snpeff_conservative_inframe_deletion_pct",
  "bcbio_umccrise", "conservative_inframe_insertion", "snpeff_conservative_inframe_insertion",
  "bcbio_umccrise", "conservative_inframe_insertion_percent", "snpeff_conservative_inframe_insertion_pct",
  "bcbio_umccrise", "disruptive_inframe_deletion", "snpeff_disruptive_inframe_deletion",
  "bcbio_umccrise", "disruptive_inframe_deletion_percent", "snpeff_disruptive_inframe_deletion_pct",
  "bcbio_umccrise", "disruptive_inframe_insertion", "snpeff_disruptive_inframe_insertion",
  "bcbio_umccrise", "disruptive_inframe_insertion_percent", "snpeff_disruptive_inframe_insertion_pct",
  "bcbio_umccrise", "downstream_gene_variant", "snpeff_downstream_gene_variant",
  "bcbio_umccrise", "downstream_gene_variant_percent", "snpeff_downstream_gene_variant_pct",
  "bcbio_umccrise", "frameshift_variant", "snpeff_frameshift_variant",
  "bcbio_umccrise", "frameshift_variant_percent", "snpeff_frameshift_variant_pct",
  "bcbio_umccrise", "initiator_codon_variant", "snpeff_initiator_codon_variant",
  "bcbio_umccrise", "initiator_codon_variant_percent", "snpeff_initiator_codon_variant_pct",
  "bcbio_umccrise", "intergenic_region", "snpeff_intergenic_region",
  "bcbio_umccrise", "intergenic_region_percent", "snpeff_intergenic_region_pct",
  "bcbio_umccrise", "intragenic_variant", "snpeff_intragenic_variant",
  "bcbio_umccrise", "intragenic_variant_percent", "snpeff_intragenic_variant_pct",
  "bcbio_umccrise", "intron_variant", "snpeff_intron_variant",
  "bcbio_umccrise", "intron_variant_percent", "snpeff_intron_variant_pct",
  "bcbio_umccrise", "missense_variant", "snpeff_missense_variant",
  "bcbio_umccrise", "missense_variant_percent", "snpeff_missense_variant_pct",
  "bcbio_umccrise", "non_coding_transcript_exon_variant", "snpeff_non_coding_transcript_exon_variant",
  "bcbio_umccrise", "non_coding_transcript_exon_variant_percent", "snpeff_non_coding_transcript_exon_variant_pct",
  "bcbio_umccrise", "non_coding_transcript_variant", "snpeff_non_coding_transcript_variant",
  "bcbio_umccrise", "non_coding_transcript_variant_percent", "snpeff_non_coding_transcript_variant_pct",
  "bcbio_umccrise", "protein_protein_contact", "snpeff_protein_protein_contact",
  "bcbio_umccrise", "protein_protein_contact_percent", "snpeff_protein_protein_contact_pct",
  "bcbio_umccrise", "sequence_feature", "snpeff_sequence_feature",
  "bcbio_umccrise", "sequence_feature_percent", "snpeff_sequence_feature_pct",
  "bcbio_umccrise", "splice_acceptor_variant", "snpeff_splice_acceptor_variant",
  "bcbio_umccrise", "splice_acceptor_variant_percent", "snpeff_splice_acceptor_variant_pct",
  "bcbio_umccrise", "splice_donor_variant", "snpeff_splice_donor_variant",
  "bcbio_umccrise", "splice_donor_variant_percent", "snpeff_splice_donor_variant_pct",
  "bcbio_umccrise", "splice_region_variant", "snpeff_splice_region_variant",
  "bcbio_umccrise", "splice_region_variant_percent", "snpeff_splice_region_variant_pct",
  "bcbio_umccrise", "start_lost", "snpeff_start_lost",
  "bcbio_umccrise", "start_lost_percent", "snpeff_start_lost_pct",
  "bcbio_umccrise", "stop_gained", "snpeff_stop_gained",
  "bcbio_umccrise", "stop_gained_percent", "snpeff_stop_gained_pct",
  "bcbio_umccrise", "stop_lost", "snpeff_stop_lost",
  "bcbio_umccrise", "stop_lost_percent", "snpeff_stop_lost_pct",
  "bcbio_umccrise", "structural_interaction_variant", "snpeff_structural_interaction_variant",
  "bcbio_umccrise", "structural_interaction_variant_percent", "snpeff_structural_interaction_variant_pct",
  "bcbio_umccrise", "synonymous_variant", "snpeff_synonymous_variant",
  "bcbio_umccrise", "synonymous_variant_percent", "snpeff_synonymous_variant_pct",
  "bcbio_umccrise", "upstream_gene_variant", "snpeff_upstream_gene_variant",
  "bcbio_umccrise", "upstream_gene_variant_percent", "snpeff_upstream_gene_variant_pct",
  "bcbio_umccrise", "DOWNSTREAM", "snpeff_downstream",
  "bcbio_umccrise", "DOWNSTREAM_percent", "snpeff_downstream_pct",
  "bcbio_umccrise", "EXON", "snpeff_exon",
  "bcbio_umccrise", "EXON_percent", "snpeff_exon_pct",
  "bcbio_umccrise", "INTERGENIC", "snpeff_intergenic",
  "bcbio_umccrise", "INTERGENIC_percent", "snpeff_intergenic_pct",
  "bcbio_umccrise", "INTRON", "snpeff_intron",
  "bcbio_umccrise", "INTRON_percent", "snpeff_intron_pct",
  "bcbio_umccrise", "MOTIF", "snpeff_motif",
  "bcbio_umccrise", "MOTIF_percent", "snpeff_motif_pct",
  "bcbio_umccrise", "SPLICE_SITE_ACCEPTOR", "snpeff_splice_site_acceptor",
  "bcbio_umccrise", "SPLICE_SITE_ACCEPTOR_percent", "snpeff_splice_site_acceptor_pct",
  "bcbio_umccrise", "SPLICE_SITE_DONOR", "snpeff_splice_site_donor",
  "bcbio_umccrise", "SPLICE_SITE_DONOR_percent", "snpeff_splice_site_donor_pct",
  "bcbio_umccrise", "SPLICE_SITE_REGION", "snpeff_splice_site_region",
  "bcbio_umccrise", "SPLICE_SITE_REGION_percent", "snpeff_splice_site_region_pct",
  "bcbio_umccrise", "TRANSCRIPT", "snpeff_transcript",
  "bcbio_umccrise", "TRANSCRIPT_percent", "snpeff_transcript_pct",
  "bcbio_umccrise", "UPSTREAM", "snpeff_upstream",
  "bcbio_umccrise", "UPSTREAM_percent", "snpeff_upstream_pct",
  "bcbio_umccrise", "UTR_3_PRIME", "snpeff_utr_3_prime",
  "bcbio_umccrise", "UTR_3_PRIME_percent", "snpeff_utr_3_prime_pct",
  "bcbio_umccrise", "UTR_5_PRIME", "snpeff_utr_5_prime",
  "bcbio_umccrise", "UTR_5_PRIME_percent", "snpeff_utr_5_prime_pct",
  "bcbio_umccrise", "Transitions", "snpeff_transitions",
  "bcbio_umccrise", "Transversions", "snpeff_transversions",
  "bcbio_umccrise", "Ts_Tv_ratio", "snpeff_ts_tv_ratio",
  "bcbio_umccrise", "Het", "snpeff_het",
  "bcbio_umccrise", "Hom", "snpeff_hom",
  "bcbio_umccrise", "Missing", "snpeff_missing",
  "bcbio_umccrise", "TFBS_ablation", "snpeff_tfbs_ablation",
  "bcbio_umccrise", "TFBS_ablation_percent", "snpeff_tfbs_ablation_pct",
  "bcbio_umccrise", "bidirectional_gene_fusion", "snpeff_bidirectional_gene_fusion",
  "bcbio_umccrise", "bidirectional_gene_fusion_percent", "snpeff_bidirectional_gene_fusion_pct",
  "bcbio_umccrise", "gene_fusion", "snpeff_gene_fusion",
  "bcbio_umccrise", "gene_fusion_percent", "snpeff_gene_fusion_pct",
  "bcbio_umccrise", "stop_retained_variant", "snpeff_stop_retained_variant",
  "bcbio_umccrise", "stop_retained_variant_percent", "snpeff_stop_retained_variant_pct",
  "bcbio_umccrise", "GENE", "snpeff_gene",
  "bcbio_umccrise", "GENE_percent", "snpeff_gene_pct"
)
