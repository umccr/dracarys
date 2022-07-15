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
  "dragen_alignment_qc", "Estimated read length", "read_length",
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
  "dragen_tumor_normal", "Filtered indels", "vars_indel_filt"
)
