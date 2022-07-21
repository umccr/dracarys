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
  config_creation_date <- sub(", ", "_", p$config_creation_date)
  workflow <- .multiqc_guess_workflow(p)
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

.multiqc_check_cols <- function(d) {
}

.multiqc_guess_workflow <- function(p) {
  assertthat::assert_that(all(c("config_title", "report_data_sources") %in% names(p)))
  config_title <- dplyr::if_else(is.null(p[["config_title"]]), "Unknown", p[["config_title"]])
  ds <- names(p[["report_data_sources"]])
  # bcbio
  if ("bcbio" %in% ds) {
    # wgs, wts, or umccrise?
    if (all(c("Salmon", "STAR", "QualiMap") %in% ds)) {
      return("bcbio_wts")
    } else if (all(c("PURPLE", "Conpair", "mosdepth", "SnpEff") %in% ds)) {
      return("bcbio_umccrise")
    } else if (all(c("Samtools", "Bcftools (somatic)", "Bcftools (germline)") %in% ds)) {
      return("bcbio_wgs")
    } else {
      warning(glue::glue(
        "Unknown which bcbio workflow this MultiQC JSON was generated from",
      ))
      return("bcbio_unknown")
    }
  }

  # dragen
  if ("DRAGEN" %in% ds) {
    dragen_workflows <- c("alignment", "transcriptome", "somatic")
    if (all(c("PURPLE", "Conpair", "mosdepth") %in% ds)) {
      return("dragen_umccrise")
    } else if (grepl("^UMCCR MultiQC Dragen", config_title)) {
      w <- tolower(sub("UMCCR MultiQC Dragen (.*) Report for .*", "\\1", config_title))
      assertthat::assert_that(w %in% dragen_workflows)
      return(paste0("dragen_", w))
    } else {
      warning(glue::glue(
        "config_title: '{config_title}'.\n",
        "Unknown which DRAGEN workflow this MultiQC JSON was generated from",
      ))
      return("dragen_unknown")
    }
  }
  return("UNKNOWN")
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
  "dragen_alignment_qc", "Total input reads", "reads_tot_input_dragen",
  "dragen_alignment_qc", "Total input reads pct", "reads_tot_input_pct_dragen",
  "dragen_alignment_qc", "Number of duplicate marked reads", "reads_num_dupmarked_dragen",
  "dragen_alignment_qc", "Number of duplicate marked reads pct", "reads_num_dupmarked_pct_dragen",
  "dragen_alignment_qc", "Number of unique reads (excl. duplicate marked reads) pct", "reads_num_unique_pct_dragen",
  "dragen_alignment_qc", "Reads with mate sequenced", "reads_w_mate_seq_dragen",
  "dragen_alignment_qc", "Reads with mate sequenced pct", "reads_w_mate_seq_pct_dragen",
  "dragen_alignment_qc", "Reads without mate sequenced", "reads_wo_mate_seq_dragen",
  "dragen_alignment_qc", "Reads without mate sequenced pct", "reads_wo_mate_seq_pct_dragen",
  "dragen_alignment_qc", "QC-failed reads", "reads_qcfail_dragen",
  "dragen_alignment_qc", "QC-failed reads pct", "reads_qcfail_pct_dragen",
  "dragen_alignment_qc", "Mapped reads", "reads_mapped_dragen",
  "dragen_alignment_qc", "Mapped reads pct", "reads_mapped_pct_dragen",
  "dragen_alignment_qc", "Mapped reads adjusted for filtered mapping", "reads_mapped_adjfilt_dragen",
  "dragen_alignment_qc", "Mapped reads adjusted for filtered mapping pct", "reads_mapped_adjfilt_pct_dragen",
  "dragen_alignment_qc", "Mapped reads R1", "reads_mapped_r1_dragen",
  "dragen_alignment_qc", "Mapped reads R1 pct", "reads_mapped_r1_pct_dragen",
  "dragen_alignment_qc", "Mapped reads R2", "reads_mapped_r2_dragen",
  "dragen_alignment_qc", "Mapped reads R2 pct", "reads_mapped_r2_pct_dragen",
  "dragen_alignment_qc", "Number of unique & mapped reads (excl. duplicate marked reads)", "reads_num_uniq_mapped_dragen",
  "dragen_alignment_qc", "Number of unique & mapped reads (excl. duplicate marked reads) pct", "reads_num_uniq_mapped_pct_dragen",
  "dragen_alignment_qc", "Unmapped reads", "reads_unmapped_dragen",
  "dragen_alignment_qc", "Unmapped reads pct", "reads_unmapped_pct_dragen",
  "dragen_alignment_qc", "Unmapped reads adjusted for filtered mapping", "reads_unmapped_adjfilt_dragen",
  "dragen_alignment_qc", "Unmapped reads adjusted for filtered mapping pct", "reads_unmapped_adjfilt_pct_dragen",
  "dragen_alignment_qc", "Adjustment of reads matching non-reference decoys", "reads_match_nonref_decoys_adj_dragen",
  "dragen_alignment_qc", "Adjustment of reads matching non-reference decoys pct", "reads_match_nonref_decoys_adj_pct_dragen",
  "dragen_alignment_qc", "Singleton reads (itself mapped; mate unmapped)", "reads_singleton_dragen",
  "dragen_alignment_qc", "Singleton reads (itself mapped; mate unmapped) pct", "reads_singleton_pct_dragen",
  "dragen_alignment_qc", "Paired reads (itself & mate mapped)", "reads_paired_dragen",
  "dragen_alignment_qc", "Paired reads (itself & mate mapped) pct", "reads_paired_pct_dragen",
  "dragen_alignment_qc", "Properly paired reads", "reads_paired_proper_dragen",
  "dragen_alignment_qc", "Properly paired reads pct", "reads_paired_proper_pct_dragen",
  "dragen_alignment_qc", "Not properly paired reads (discordant)", "reads_discordant_dragen",
  "dragen_alignment_qc", "Not properly paired reads (discordant) pct", "reads_discordant_pct_dragen",
  "dragen_alignment_qc", "Paired reads mapped to different chromosomes", "reads_paired_mapped_diff_chrom_dragen",
  "dragen_alignment_qc", "Paired reads mapped to different chromosomes pct", "reads_paired_mapped_diff_chrom_pct_dragen",
  "dragen_alignment_qc", "Paired reads mapped to different chromosomes (MAPQ>=10)", "reads_paired_mapped_diff_chrom_mapq10_dragen",
  "dragen_alignment_qc", "Paired reads mapped to different chromosomes (MAPQ>=10) pct", "reads_paired_mapped_diff_chrom_mapq10_pct_dragen",
  "dragen_alignment_qc", "Reads with MAPQ [40:inf)", "reads_mapq_40_inf_dragen",
  "dragen_alignment_qc", "Reads with MAPQ [40:inf) pct", "reads_mapq_40_inf_pct_dragen",
  "dragen_alignment_qc", "Reads with MAPQ [30:40)", "reads_mapq_30_40_dragen",
  "dragen_alignment_qc", "Reads with MAPQ [30:40) pct", "reads_mapq_30_40_pct_dragen",
  "dragen_alignment_qc", "Reads with MAPQ [20:30)", "reads_mapq_20_30_dragen",
  "dragen_alignment_qc", "Reads with MAPQ [20:30) pct", "reads_mapq_20_30_pct_dragen",
  "dragen_alignment_qc", "Reads with MAPQ [10:20)", "reads_mapq_10_20_dragen",
  "dragen_alignment_qc", "Reads with MAPQ [10:20) pct", "reads_mapq_10_20_pct_dragen",
  "dragen_alignment_qc", "Reads with MAPQ [ 0:10)", "reads_mapq_0_10_dragen",
  "dragen_alignment_qc", "Reads with MAPQ [ 0:10) pct", "reads_mapq_0_10_pct_dragen",
  "dragen_alignment_qc", "Reads with MAPQ NA (Unmapped reads)", "reads_mapq_NA_unmapped_dragen",
  "dragen_alignment_qc", "Reads with MAPQ NA (Unmapped reads) pct", "reads_mapq_NA_unmapped_pct_dragen",
  "dragen_alignment_qc", "Reads with indel R1", "reads_indel_r1_dragen",
  "dragen_alignment_qc", "Reads with indel R1 pct", "reads_indel_r1_pct_dragen",
  "dragen_alignment_qc", "Reads with indel R2", "reads_indel_r2_dragen",
  "dragen_alignment_qc", "Reads with indel R2 pct", "reads_indel_r2_pct_dragen",
  "dragen_alignment_qc", "Total bases", "bases_tot_dragen",
  "dragen_alignment_qc", "Total bases R1", "bases_tot_r1_dragen",
  "dragen_alignment_qc", "Total bases R2", "bases_tot_r2_dragen",
  "dragen_alignment_qc", "Mapped bases R1", "bases_mapped_r1_dragen",
  "dragen_alignment_qc", "Mapped bases R1 pct", "bases_mapped_r1_pct_dragen",
  "dragen_alignment_qc", "Mapped bases R2", "bases_mapped_r2_dragen",
  "dragen_alignment_qc", "Mapped bases R2 pct", "bases_mapped_r2_pct_dragen",
  "dragen_alignment_qc", "Soft-clipped bases R1", "bases_softclip_r1_dragen",
  "dragen_alignment_qc", "Soft-clipped bases R1 pct", "bases_softclip_r1_pct_dragen",
  "dragen_alignment_qc", "Soft-clipped bases R2 pct", "bases_softclip_r2_pct_dragen",
  "dragen_alignment_qc", "Mismatched bases R1", "bases_mismatched_r1_dragen",
  "dragen_alignment_qc", "Mismatched bases R1 pct", "bases_mismatched_r1_pct_dragen",
  "dragen_alignment_qc", "Mismatched bases R2 pct", "bases_mismatched_r2_pct_dragen",
  "dragen_alignment_qc", "Mismatched bases R1 (excl. indels)", "bases_mismatched_r1_noindels_dragen",
  "dragen_alignment_qc", "Mismatched bases R1 (excl. indels) pct", "bases_mismatched_r1_noindels_pct_dragen",
  "dragen_alignment_qc", "Mismatched bases R2 (excl. indels) pct", "bases_mismatched_r2_noindels_pct_dragen",
  "dragen_alignment_qc", "Q30 bases", "bases_q30_dragen",
  "dragen_alignment_qc", "Q30 bases pct", "bases_q30_pct_dragen",
  "dragen_alignment_qc", "Q30 bases R1", "bases_q30_r1_dragen",
  "dragen_alignment_qc", "Q30 bases R1 pct", "bases_q30_r1_pct_dragen",
  "dragen_alignment_qc", "Q30 bases R2", "bases_q30_r2_dragen",
  "dragen_alignment_qc", "Q30 bases R2 pct", "bases_q30_r2_pct_dragen",
  "dragen_alignment_qc", "Q30 bases (excl. dups & clipped bases)", "bases_q30_nodups_noclipped_dragen",
  "dragen_alignment_qc", "Q30 bases (excl. dups & clipped bases) pct", "bases_q30_nodups_noclipped_pct_dragen",
  "dragen_alignment_qc", "Total alignments", "alignments_tot_dragen",
  "dragen_alignment_qc", "Secondary alignments", "alignments_secondary_dragen",
  "dragen_alignment_qc", "Secondary alignments pct", "alignments_secondary_pct_dragen",
  "dragen_alignment_qc", "Supplementary (chimeric) alignments", "alignments_chimeric_dragen",
  "dragen_alignment_qc", "Estimated read length", "read_len_dragen",
  "dragen_alignment_qc", "Bases in reference genome", "bases_in_ref_genome_dragen",
  "dragen_alignment_qc", "Bases in target bed [% of genome]", "bases_in_target_bed_genome_pct_dragen",
  "dragen_alignment_qc", "Average sequenced coverage over genome", "cov_avg_seq_over_genome_dragen",
  "dragen_alignment_qc", "Insert length: mean", "insert_len_mean_dragen",
  "dragen_alignment_qc", "Insert length: median", "insert_len_median_dragen",
  "dragen_alignment_qc", "Insert length: standard deviation", "insert_len_std_dev_dragen",
  "dragen_alignment_qc", "Provided sex chromosome ploidy", "ploidy_sex_chrom_provided_dragen",
  "dragen_alignment_qc", "Estimated sample contamination", "contamination_est_dragen",
  "dragen_alignment_qc", "DRAGEN mapping rate [mil. reads/second]", "mapping_rate_dragen_milreads_per_sec_dragen",
  "dragen_alignment_qc", "Autosomal median coverage", "cov_auto_median_dragen",
  "dragen_alignment_qc", "X median coverage", "cov_x_median_dragen",
  "dragen_alignment_qc", "Y median coverage", "cov_y_median_dragen",
  "dragen_alignment_qc", "1 median / Autosomal median", "cov_1_div_auto_medians_dragen",
  "dragen_alignment_qc", "2 median / Autosomal median", "cov_2_div_auto_medians_dragen",
  "dragen_alignment_qc", "3 median / Autosomal median", "cov_3_div_auto_medians_dragen",
  "dragen_alignment_qc", "4 median / Autosomal median", "cov_4_div_auto_medians_dragen",
  "dragen_alignment_qc", "5 median / Autosomal median", "cov_5_div_auto_medians_dragen",
  "dragen_alignment_qc", "6 median / Autosomal median", "cov_6_div_auto_medians_dragen",
  "dragen_alignment_qc", "7 median / Autosomal median", "cov_7_div_auto_medians_dragen",
  "dragen_alignment_qc", "8 median / Autosomal median", "cov_8_div_auto_medians_dragen",
  "dragen_alignment_qc", "9 median / Autosomal median", "cov_9_div_auto_medians_dragen",
  "dragen_alignment_qc", "10 median / Autosomal median", "cov_10_div_auto_median_dragen",
  "dragen_alignment_qc", "11 median / Autosomal median", "cov_11_div_auto_median_dragen",
  "dragen_alignment_qc", "12 median / Autosomal median", "cov_12_div_auto_median_dragen",
  "dragen_alignment_qc", "13 median / Autosomal median", "cov_13_div_auto_median_dragen",
  "dragen_alignment_qc", "14 median / Autosomal median", "cov_14_div_auto_median_dragen",
  "dragen_alignment_qc", "15 median / Autosomal median", "cov_15_div_auto_median_dragen",
  "dragen_alignment_qc", "16 median / Autosomal median", "cov_16_div_auto_median_dragen",
  "dragen_alignment_qc", "17 median / Autosomal median", "cov_17_div_auto_median_dragen",
  "dragen_alignment_qc", "18 median / Autosomal median", "cov_18_div_auto_median_dragen",
  "dragen_alignment_qc", "19 median / Autosomal median", "cov_19_div_auto_median_dragen",
  "dragen_alignment_qc", "20 median / Autosomal median", "cov_20_div_auto_median_dragen",
  "dragen_alignment_qc", "21 median / Autosomal median", "cov_21_div_auto_median_dragen",
  "dragen_alignment_qc", "22 median / Autosomal median", "cov_22_div_auto_median_dragen",
  "dragen_alignment_qc", "X median / Autosomal median", "cov_x_div_auto_median_dragen",
  "dragen_alignment_qc", "Y median / Autosomal median", "cov_y_div_auto_median_dragen",
  "dragen_alignment_qc", "Ploidy estimation", "ploidy_est_dragen",
  "dragen_alignment_qc", "Aligned bases", "bases_aligned_dragen",
  "dragen_alignment_qc", "Aligned bases in genome", "bases_aligned_in_genome_dragen",
  "dragen_alignment_qc", "Aligned bases in genome pct", "bases_aligned_in_genome_pct_dragen",
  "dragen_alignment_qc", "Average alignment coverage over genome", "cov_alignment_avg_over_genome_dragen",
  "dragen_alignment_qc", "Uniformity of coverage (PCT > 0.2*mean) over genome", "cov_uniformity_over_genome_pct_gt02mean_dragen",
  "dragen_alignment_qc", "PCT of genome with coverage [100x: inf)", "cov_genome_pct_100x_inf_dragen",
  "dragen_alignment_qc", "PCT of genome with coverage [ 50x: inf)", "cov_genome_pct_50x_inf_dragen",
  "dragen_alignment_qc", "PCT of genome with coverage [ 20x: inf)", "cov_genome_pct_20x_inf_dragen",
  "dragen_alignment_qc", "PCT of genome with coverage [ 15x: inf)", "cov_genome_pct_15x_inf_dragen",
  "dragen_alignment_qc", "PCT of genome with coverage [ 10x: inf)", "cov_genome_pct_10x_inf_dragen",
  "dragen_alignment_qc", "PCT of genome with coverage [  3x: inf)", "cov_genome_pct_3x_inf_dragen",
  "dragen_alignment_qc", "PCT of genome with coverage [  1x: inf)", "cov_genome_pct_1x_inf_dragen",
  "dragen_alignment_qc", "PCT of genome with coverage [  0x: inf)", "cov_genome_pct_0x_inf_dragen",
  "dragen_alignment_qc", "PCT of genome with coverage [ 50x:100x)", "cov_genome_pct_50x_100x_dragen",
  "dragen_alignment_qc", "PCT of genome with coverage [ 20x: 50x)", "cov_genome_pct_20x_50x_dragen",
  "dragen_alignment_qc", "PCT of genome with coverage [ 15x: 20x)", "cov_genome_pct_15x_20x_dragen",
  "dragen_alignment_qc", "PCT of genome with coverage [ 10x: 15x)", "cov_genome_pct_10x_15x_dragen",
  "dragen_alignment_qc", "PCT of genome with coverage [  3x: 10x)", "cov_genome_pct_3x_10x_dragen",
  "dragen_alignment_qc", "PCT of genome with coverage [  1x:  3x)", "cov_genome_pct_1x_3x_dragen",
  "dragen_alignment_qc", "PCT of genome with coverage [  0x:  1x)", "cov_genome_pct_0x_1x_dragen",
  "dragen_alignment_qc", "Average chr X coverage over genome", "cov_avg_x_over_genome_dragen",
  "dragen_alignment_qc", "Average chr Y coverage over genome", "cov_avg_y_over_genome_dragen",
  "dragen_alignment_qc", "Average mitochondrial coverage over genome", "cov_avg_mt_over_genome_dragen",
  "dragen_alignment_qc", "Average autosomal coverage over genome", "cov_avg_auto_over_genome_dragen",
  "dragen_alignment_qc", "Median autosomal coverage over genome", "cov_median_auto_over_genome_dragen",
  "dragen_alignment_qc", "Mean/Median autosomal coverage ratio over genome", "cov_mean_median_auto_ratio_over_genome_dragen",
  "dragen_alignment_qc", "Aligned reads", "reads_aligned_dragen",
  "dragen_alignment_qc", "Aligned reads in genome", "reads_aligned_in_genome_dragen",
  "dragen_alignment_qc", "Aligned reads in genome pct", "reads_aligned_in_genome_pct_dragen",
  "dragen_tumor_normal", "Total", "var_tot_dragen",
  "dragen_tumor_normal", "Total pct", "var_tot_pct_dragen",
  "dragen_tumor_normal", "Biallelic", "var_biallelic_dragen",
  "dragen_tumor_normal", "Biallelic pct", "var_biallelic_pct_dragen",
  "dragen_tumor_normal", "Multiallelic", "var_multiallelic_dragen",
  "dragen_tumor_normal", "Multiallelic pct", "var_multiallelic_pct_dragen",
  "dragen_tumor_normal", "SNPs", "var_snp_dragen",
  "dragen_tumor_normal", "SNPs pct", "var_snp_pct_dragen",
  "dragen_tumor_normal", "Insertions (Hom)", "var_ins_hom_dragen",
  "dragen_tumor_normal", "Insertions (Hom) pct", "var_ins_hom_pct_dragen",
  "dragen_tumor_normal", "Insertions (Het)", "var_ins_het_dragen",
  "dragen_tumor_normal", "Insertions (Het) pct", "var_ins_het_pct_dragen",
  "dragen_tumor_normal", "Deletions (Hom)", "var_del_hom_dragen",
  "dragen_tumor_normal", "Deletions (Hom) pct", "var_del_hom_pct_dragen",
  "dragen_tumor_normal", "Deletions (Het)", "var_del_het_dragen",
  "dragen_tumor_normal", "Deletions (Het) pct", "var_del_het_pct_dragen",
  "dragen_tumor_normal", "Indels (Het)", "var_indel_het_dragen",
  "dragen_tumor_normal", "Indels (Het) pct", "var_indel_het_pct_dragen",
  "dragen_tumor_normal", "Chr X number of SNPs over genome", "var_snp_x_over_genome_dragen",
  "dragen_tumor_normal", "Chr Y number of SNPs over genome", "var_snp_y_over_genome_dragen",
  "dragen_tumor_normal", "(Chr X SNPs)/(chr Y SNPs) ratio over genome", "var_x_over_y_snp_ratio_over_genome_dragen",
  "dragen_tumor_normal", "SNP Transitions", "var_snp_transitions_dragen",
  "dragen_tumor_normal", "SNP Transversions", "var_snp_transversions_dragen",
  "dragen_tumor_normal", "Ti/Tv ratio", "var_ti_tv_ratio_dragen",
  "dragen_tumor_normal", "Heterozygous", "var_heterozygous_dragen",
  "dragen_tumor_normal", "Homozygous", "var_homozygous_dragen",
  "dragen_tumor_normal", "Het/Hom ratio", "var_het_hom_ratio_dragen",
  "dragen_tumor_normal", "In dbSNP", "var_in_dbsnp_dragen",
  "dragen_tumor_normal", "In dbSNP pct", "var_in_dbsnp_pct_dragen",
  "dragen_tumor_normal", "Not in dbSNP", "var_nin_dbsnp_dragen",
  "dragen_tumor_normal", "Not in dbSNP pct", "var_nin_dbsnp_pct_dragen",
  "dragen_tumor_normal", "Percent Callability", "callability_pct_dragen",
  "dragen_tumor_normal", "Percent Autosome Callability", "callability_auto_pct_dragen",
  "dragen_tumor_normal", "Insertions", "var_ins_tot_dragen",
  "dragen_tumor_normal", "Insertions pct", "var_ins_tot_pct_dragen",
  "dragen_tumor_normal", "Deletions", "var_del_tot_dragen",
  "dragen_tumor_normal", "Deletions pct", "var_del_tot_pct_dragen",
  "dragen_tumor_normal", "Indels", "var_indel_tot_dragen",
  "dragen_tumor_normal", "Indels pct", "var_indel_tot_pct_dragen",
  "dragen_tumor_normal", "Number of samples", "sample_num_dragen",
  "dragen_tumor_normal", "Reads Processed", "reads_processed_dragen",
  "dragen_tumor_normal", "Child Sample", "sample_child_dragen",
  "dragen_tumor_normal", "Filtered vars", "vars_tot_filt_dragen",
  "dragen_tumor_normal", "Filtered vars pct", "vars_tot_filt_pct_dragen",
  "dragen_tumor_normal", "Filtered SNPs", "vars_snp_filt_dragen",
  "dragen_tumor_normal", "Filtered SNPs pct", "vars_snp_filt_pct_dragen",
  "dragen_tumor_normal", "Filtered indels", "vars_indel_filt_dragen",
  "dragen_transcriptome", "Adjustment of reads matching filter contigs", "reads_match_filter_contigs_adj_dragen",
  "dragen_transcriptome", "Adjustment of reads matching filter contigs pct", "reads_match_filter_contigs_adj_pct_dragen",
  "dragen_transcriptome", "Reads with splice junction", "reads_splice_junction_dragen",
  "dragen_transcriptome", "Reads with splice junction pct", "reads_splice_junction_pct_dragen",
  "dragen_umccrise", "concordance_concordance", "concordance_conpair",
  "dragen_umccrise", "concordance_used_markers", "concordance_used_markers_conpair",
  "dragen_umccrise", "concordance_total_markers", "concordance_total_markers_conpair",
  "dragen_umccrise", "concordance_marker_threshold", "concordance_marker_threshold_conpair",
  "dragen_umccrise", "concordance_min_mapping_quality", "concordance_min_mapping_quality_conpair",
  "dragen_umccrise", "concordance_min_base_quality", "concordance_min_base_quality_conpair",
  "dragen_umccrise", "contamination", "contamination_conpair",
  "dragen_umccrise", "QCStatus", "qc_status_purple",
  "dragen_umccrise", "Method", "method_purple",
  "dragen_umccrise", "CopyNumberSegments", "cn_segs_purple",
  "dragen_umccrise", "UnsupportedCopyNumberSegments", "unsupported_cn_segs_purple",
  "dragen_umccrise", "Purity", "purity_purple",
  "dragen_umccrise", "AmberGender", "gender_amber",
  "dragen_umccrise", "CobaltGender", "gender_cobalt",
  "dragen_umccrise", "DeletedGenes", "deleted_genes_purple",
  "dragen_umccrise", "Contamination", "contamination_purple",
  "dragen_umccrise", "GermlineAberrations", "germline_aberrations_purple",
  "dragen_umccrise", "purity", "purity_purple",
  "dragen_umccrise", "normFactor", "normfactor_purple",
  "dragen_umccrise", "score", "score_purple",
  "dragen_umccrise", "diploidProportion", "diploid_prop_purple",
  "dragen_umccrise", "ploidy", "ploidy_purple",
  "dragen_umccrise", "gender", "gender_purple",
  "dragen_umccrise", "status", "status_purple",
  "dragen_umccrise", "polyclonalProportion", "polyclonal_prop_purple",
  "dragen_umccrise", "minPurity", "min_purity_purple",
  "dragen_umccrise", "maxPurity", "max_purity_purple",
  "dragen_umccrise", "minPloidy", "min_ploidy_purple",
  "dragen_umccrise", "maxPloidy", "max_ploidy_purple",
  "dragen_umccrise", "minDiploidProportion", "min_diploid_prop_purple",
  "dragen_umccrise", "maxDiploidProportion", "max_diploid_prop_purple",
  "dragen_umccrise", "version", "version_purple",
  "dragen_umccrise", "somaticPenalty", "somatic_penalty_purple",
  "dragen_umccrise", "wholeGenomeDuplication", "whole_genome_dup_purple",
  "dragen_umccrise", "msIndelsPerMb", "ms_indels_permb_purple",
  "dragen_umccrise", "msStatus", "ms_status_purple",
  "dragen_umccrise", "tml", "tml_purple",
  "dragen_umccrise", "tmlStatus", "tml_status_purple",
  "dragen_umccrise", "tmbPerMb", "tmb_permb_purple",
  "dragen_umccrise", "tmbStatus", "tmb_status_purple",
  "dragen_umccrise", "svTumorMutationalBurden", "tmb_sv_purple",
  "dragen_umccrise", "1_x_pc", "cov_1x_pc_mosdepth",
  "dragen_umccrise", "5_x_pc", "cov_5x_pc_mosdepth",
  "dragen_umccrise", "10_x_pc", "cov_10x_pc_mosdepth",
  "dragen_umccrise", "30_x_pc", "cov_30x_pc_mosdepth",
  "dragen_umccrise", "50_x_pc", "cov_50x_pc_mosdepth",
  "dragen_umccrise", "median_coverage", "cov_median_mosdepth",
  "dragen_umccrise", "raw_total_sequences", "raw_total_sequences_samtools",
  "dragen_umccrise", "filtered_sequences", "filtered_sequences_samtools",
  "dragen_umccrise", "sequences", "sequences_samtools",
  "dragen_umccrise", "is_sorted", "is_sorted_samtools",
  "dragen_umccrise", "1st_fragments", "1st_fragments_samtools",
  "dragen_umccrise", "last_fragments", "last_fragments_samtools",
  "dragen_umccrise", "reads_mapped", "reads_mapped_samtools",
  "dragen_umccrise", "reads_mapped_and_paired", "reads_mapped_and_paired_samtools",
  "dragen_umccrise", "reads_unmapped", "reads_unmapped_samtools",
  "dragen_umccrise", "reads_properly_paired", "reads_properly_paired_samtools",
  "dragen_umccrise", "reads_paired", "reads_paired_samtools",
  "dragen_umccrise", "reads_duplicated", "reads_num_dupmarked_samtools",
  "dragen_umccrise", "reads_duplicated_percent", "reads_num_dupmarked_pct_samtools",
  "dragen_umccrise", "reads_MQ0", "reads_mq0_samtools",
  "dragen_umccrise", "reads_QC_failed", "reads_qc_failed_samtools",
  "dragen_umccrise", "non-primary_alignments", "non_primary_alignments_samtools",
  "dragen_umccrise", "total_length", "total_length_samtools",
  "dragen_umccrise", "total_first_fragment_length", "total_first_fragment_length_samtools",
  "dragen_umccrise", "total_last_fragment_length", "total_last_fragment_length_samtools",
  "dragen_umccrise", "bases_mapped", "bases_mapped_samtools",
  "dragen_umccrise", "bases_mapped_(cigar)", "bases_mapped_cigar_samtools",
  "dragen_umccrise", "bases_trimmed", "bases_trimmed_samtools",
  "dragen_umccrise", "bases_duplicated", "bases_duplicated_samtools",
  "dragen_umccrise", "mismatches", "mismatches_samtools",
  "dragen_umccrise", "error_rate", "error_rate_samtools",
  "dragen_umccrise", "average_length", "avg_length_samtools",
  "dragen_umccrise", "average_first_fragment_length", "avg_first_fragment_length_samtools",
  "dragen_umccrise", "average_last_fragment_length", "avg_last_fragment_length_samtools",
  "dragen_umccrise", "maximum_length", "max_length_samtools",
  "dragen_umccrise", "maximum_first_fragment_length", "max_first_fragment_length_samtools",
  "dragen_umccrise", "maximum_last_fragment_length", "max_last_fragment_length_samtools",
  "dragen_umccrise", "average_quality", "avg_quality_samtools",
  "dragen_umccrise", "insert_size_average", "insert_size_avg_samtools",
  "dragen_umccrise", "insert_size_standard_deviation", "insert_size_stddev_samtools",
  "dragen_umccrise", "inward_oriented_pairs", "inward_oriented_pairs_samtools",
  "dragen_umccrise", "outward_oriented_pairs", "outward_oriented_pairs_samtools",
  "dragen_umccrise", "pairs_with_other_orientation", "pairs_with_other_orientation_samtools",
  "dragen_umccrise", "pairs_on_different_chromosomes", "pairs_on_different_chromosomes_samtools",
  "dragen_umccrise", "percentage_of_properly_paired_reads_(%)", "reads_properly_paired_pct_round_samtools",
  "dragen_umccrise", "reads_properly_paired_percent", "reads_properly_paired_pct_samtools",
  "dragen_umccrise", "reads_mapped_percent", "reads_mapped_pct_samtools",
  "dragen_umccrise", "reads_mapped_and_paired_percent", "reads_mapped_and_paired_pct_samtools",
  "dragen_umccrise", "reads_unmapped_percent", "reads_unmapped_pct_samtools",
  "dragen_umccrise", "reads_paired_percent", "reads_paired_pct_samtools",
  "dragen_umccrise", "reads_MQ0_percent", "reads_mq0_pct_samtools",
  "dragen_umccrise", "reads_QC_failed_percent", "reads_qc_failed_pct_samtools",
  "dragen_umccrise", "filt_indels", "filt_indels_umccrise",
  "dragen_umccrise", "filt_snps", "filt_snps_umccrise",
  "dragen_umccrise", "filt_vars", "filt_vars_umccrise",
  "dragen_umccrise", "indels", "indel_umccrise",
  "dragen_umccrise", "snps", "snp_umccrise",
  "dragen_umccrise", "viral_content", "viral_content_umccrise",
  "dragen_umccrise", "germline", "germline_umccrise",
  "dragen_umccrise", "germline_predispose", "germline_predispose_umccrise",
  "bcbio_umccrise", "concordance", "concordance_conpair",
  "bcbio_umccrise", "concordance_markers", "concordance_used_over_tot_markers_conpair",
  "bcbio_umccrise", "SegmentPass", "segment_pass_purple",
  "bcbio_umccrise", "GenderPass", "gender_pass_purple",
  "bcbio_umccrise", "DeletedGenesPass", "deleted_genes_pass_purple",
  "bcbio_umccrise", "SegmentScore", "segment_score_purple",
  "bcbio_umccrise", "UnsupportedSegments", "unsupported_cn_segs_purple",
  "bcbio_umccrise", "Ploidy", "ploidy_purple",
  "bcbio_umccrise", "%GC", "gc_pct_bcbio",
  "bcbio_umccrise", "Average_insert_size", "insert_size_avg_bcbio",
  "bcbio_umccrise", "Average_read_length", "read_len_avg_bcbio",
  "bcbio_umccrise", "Avg_coverage", "cov_avg_bcbio",
  "bcbio_umccrise", "Duplicates", "reads_num_dupmarked_bcbio",
  "bcbio_umccrise", "Duplicates_pct", "reads_num_dupmarked_pct_bcbio",
  "bcbio_umccrise", "Mapped_paired_reads", "reads_paired_bcbio",
  "bcbio_umccrise", "Mapped_reads", "reads_mapped_bcbio",
  "bcbio_umccrise", "Mapped_reads_pct", "reads_mapped_pct_bcbio",
  "bcbio_umccrise", "Mapped_unique_reads", "reads_num_uniq_mapped_bcbio",
  "bcbio_umccrise", "Offtarget_pct", "reads_offtarget_pct_bcbio",
  "bcbio_umccrise", "Ontarget_pct", "reads_ontarget_pct_bcbio",
  "bcbio_umccrise", "Ontarget_unique_reads", "reads_ontarget_unique_bcbio",
  "bcbio_umccrise", "Sequences_flagged_as_poor_quality", "sequences_poor_quality_bcbio",
  "bcbio_umccrise", "Total_reads", "reads_tot_bcbio",
  "bcbio_umccrise", "Usable_pct", "reads_usable_pct_bcbio",
  "bcbio_umccrise", "RiP_pct", "reads_rip_pct_bcbio",
  "bcbio_umccrise", "non-primary_alignments_percent", "alignments_secondary_pct_samtools",
  "bcbio_umccrise", "family_id", "family_id_peddy",
  "bcbio_umccrise", "paternal_id", "paternal_id_peddy",
  "bcbio_umccrise", "maternal_id", "maternal_id_peddy",
  "bcbio_umccrise", "sex", "sex_peddy",
  "bcbio_umccrise", "phenotype", "phenotype_peddy",
  "bcbio_umccrise", "Ethnicity", "ethnicity_peddy",
  "bcbio_umccrise", "duplicates", "duplicates_peddy",
  "bcbio_umccrise", "het_call_rate", "het_call_rate_peddy",
  "bcbio_umccrise", "het_ratio", "het_ratio_peddy",
  "bcbio_umccrise", "het_mean_depth", "het_mean_depth_peddy",
  "bcbio_umccrise", "het_idr_baf", "het_idr_baf_peddy",
  "bcbio_umccrise", "ancestry-prediction", "ancestry_prediction_peddy",
  "bcbio_umccrise", "PC1", "pc1_peddy",
  "bcbio_umccrise", "PC2", "pc2_peddy",
  "bcbio_umccrise", "PC3", "pc3_peddy",
  "bcbio_umccrise", "sex_het_ratio", "sex_het_ratio_peddy",
  "bcbio_umccrise", "sex_fixed", "sex_fixed_peddy",
  "bcbio_umccrise", "depth_outlier_het_check", "depth_outlier_het_check_peddy",
  "bcbio_umccrise", "het_count_het_check", "het_count_het_check_peddy",
  "bcbio_umccrise", "het_ratio_het_check", "het_ratio_het_check_peddy",
  "bcbio_umccrise", "idr_baf_het_check", "idr_baf_het_check_peddy",
  "bcbio_umccrise", "mean_depth_het_check", "mean_depth_het_check_peddy",
  "bcbio_umccrise", "median_depth_het_check", "median_depth_het_check_peddy",
  "bcbio_umccrise", "p10_het_check", "p10_het_check_peddy",
  "bcbio_umccrise", "p90_het_check", "p90_het_check_peddy",
  "bcbio_umccrise", "sampled_sites_het_check", "sampled_sites_het_check_peddy",
  "bcbio_umccrise", "call_rate_het_check", "call_rate_het_check_peddy",
  "bcbio_umccrise", "PC1_het_check", "pc1_het_check_peddy",
  "bcbio_umccrise", "PC2_het_check", "pc2_het_check_peddy",
  "bcbio_umccrise", "PC3_het_check", "pc3_het_check_peddy",
  "bcbio_umccrise", "PC4_het_check", "pc4_het_check_peddy",
  "bcbio_umccrise", "ancestry-prediction_het_check", "ancestry_prediction_het_check_peddy",
  "bcbio_umccrise", "ancestry-prob_het_check", "ancestry_prob_het_check_peddy",
  "bcbio_umccrise", "error_sex_check", "error_sex_check_peddy",
  "bcbio_umccrise", "het_count_sex_check", "het_count_sex_check_peddy",
  "bcbio_umccrise", "het_ratio_sex_check", "het_ratio_sex_check_peddy",
  "bcbio_umccrise", "hom_alt_count_sex_check", "hom_alt_count_sex_check_peddy",
  "bcbio_umccrise", "hom_ref_count_sex_check", "hom_ref_count_sex_check_peddy",
  "bcbio_umccrise", "ped_sex_sex_check", "ped_sex_sex_check_peddy",
  "bcbio_umccrise", "predicted_sex_sex_check", "predicted_sex_sex_check_peddy",
  "bcbio_umccrise", "percent_gc", "gc_pct_fastqc",
  "bcbio_umccrise", "avg_sequence_length", "avg_sequence_length_fastqc",
  "bcbio_umccrise", "total_sequences", "sequences_tot_fastqc",
  "bcbio_umccrise", "percent_duplicates", "reads_num_dupmarked_pct_fastqc",
  "bcbio_umccrise", "percent_fails", "fastqc_fails_pct_fastqc",
  "bcbio_umccrise", "Genome", "genome_snpeff",
  "bcbio_umccrise", "Number_of_variants_before_filter", "var_num_before_filter_snpeff",
  "bcbio_umccrise", "Number_of_known_variants (i.e. non-empty ID)", "var_num_known_snpeff",
  "bcbio_umccrise", "Number_of_known_variants (i.e. non-empty ID)_percent", "var_num_known_pct_snpeff",
  "bcbio_umccrise", "Number_of_effects", "effects_num_snpeff",
  "bcbio_umccrise", "Genome_total_length", "genome_tot_len_snpeff",
  "bcbio_umccrise", "Genome_effective_length", "genome_eff_len_snpeff",
  "bcbio_umccrise", "Change_rate", "change_rate_snpeff",
  "bcbio_umccrise", "HIGH", "high_snpeff",
  "bcbio_umccrise", "HIGH_percent", "high_pct_snpeff",
  "bcbio_umccrise", "LOW", "low_snpeff",
  "bcbio_umccrise", "LOW_percent", "low_pct_snpeff",
  "bcbio_umccrise", "MODERATE", "moderate_snpeff",
  "bcbio_umccrise", "MODERATE_percent", "moderate_pct_snpeff",
  "bcbio_umccrise", "MODIFIER", "modifier_snpeff",
  "bcbio_umccrise", "MODIFIER_percent", "modifier_pct_snpeff",
  "bcbio_umccrise", "MISSENSE", "missense_snpeff",
  "bcbio_umccrise", "MISSENSE_percent", "missense_pct_snpeff",
  "bcbio_umccrise", "NONSENSE", "nonsense_snpeff",
  "bcbio_umccrise", "NONSENSE_percent", "nonsense_pct_snpeff",
  "bcbio_umccrise", "SILENT", "silent_snpeff",
  "bcbio_umccrise", "SILENT_percent", "silent_pct_snpeff",
  "bcbio_umccrise", "Missense_Silent_ratio", "missense_silent_ratio_snpeff",
  "bcbio_umccrise", "Type", "type_snpeff",
  "bcbio_umccrise", "3_prime_UTR_variant", "3_prime_utr_variant_snpeff",
  "bcbio_umccrise", "3_prime_UTR_variant_percent", "3_prime_utr_variant_pct_snpeff",
  "bcbio_umccrise", "5_prime_UTR_premature_start_codon_gain_variant", "5_prime_utr_premature_start_codon_gain_variant_snpeff",
  "bcbio_umccrise", "5_prime_UTR_premature_start_codon_gain_variant_percent", "5_prime_utr_premature_start_codon_gain_variant_pct_snpeff",
  "bcbio_umccrise", "5_prime_UTR_variant", "5_prime_utr_variant_snpeff",
  "bcbio_umccrise", "5_prime_UTR_variant_percent", "5_prime_utr_variant_pct_snpeff",
  "bcbio_umccrise", "TF_binding_site_variant", "tf_binding_site_variant_snpeff",
  "bcbio_umccrise", "TF_binding_site_variant_percent", "tf_binding_site_variant_pct_snpeff",
  "bcbio_umccrise", "conservative_inframe_deletion", "conservative_inframe_deletion_snpeff",
  "bcbio_umccrise", "conservative_inframe_deletion_percent", "conservative_inframe_deletion_pct_snpeff",
  "bcbio_umccrise", "conservative_inframe_insertion", "conservative_inframe_insertion_snpeff",
  "bcbio_umccrise", "conservative_inframe_insertion_percent", "conservative_inframe_insertion_pct_snpeff",
  "bcbio_umccrise", "disruptive_inframe_deletion", "disruptive_inframe_deletion_snpeff",
  "bcbio_umccrise", "disruptive_inframe_deletion_percent", "disruptive_inframe_deletion_pct_snpeff",
  "bcbio_umccrise", "disruptive_inframe_insertion", "disruptive_inframe_insertion_snpeff",
  "bcbio_umccrise", "disruptive_inframe_insertion_percent", "disruptive_inframe_insertion_pct_snpeff",
  "bcbio_umccrise", "downstream_gene_variant", "downstream_gene_variant_snpeff",
  "bcbio_umccrise", "downstream_gene_variant_percent", "downstream_gene_variant_pct_snpeff",
  "bcbio_umccrise", "frameshift_variant", "frameshift_variant_snpeff",
  "bcbio_umccrise", "frameshift_variant_percent", "frameshift_variant_pct_snpeff",
  "bcbio_umccrise", "initiator_codon_variant", "initiator_codon_variant_snpeff",
  "bcbio_umccrise", "initiator_codon_variant_percent", "initiator_codon_variant_pct_snpeff",
  "bcbio_umccrise", "intergenic_region", "intergenic_region_snpeff",
  "bcbio_umccrise", "intergenic_region_percent", "intergenic_region_pct_snpeff",
  "bcbio_umccrise", "intragenic_variant", "intragenic_variant_snpeff",
  "bcbio_umccrise", "intragenic_variant_percent", "intragenic_variant_pct_snpeff",
  "bcbio_umccrise", "intron_variant", "intron_variant_snpeff",
  "bcbio_umccrise", "intron_variant_percent", "intron_variant_pct_snpeff",
  "bcbio_umccrise", "missense_variant", "missense_variant_snpeff",
  "bcbio_umccrise", "missense_variant_percent", "missense_variant_pct_snpeff",
  "bcbio_umccrise", "non_coding_transcript_exon_variant", "non_coding_transcript_exon_variant_snpeff",
  "bcbio_umccrise", "non_coding_transcript_exon_variant_percent", "non_coding_transcript_exon_variant_pct_snpeff",
  "bcbio_umccrise", "non_coding_transcript_variant", "non_coding_transcript_variant_snpeff",
  "bcbio_umccrise", "non_coding_transcript_variant_percent", "non_coding_transcript_variant_pct_snpeff",
  "bcbio_umccrise", "protein_protein_contact", "protein_protein_contact_snpeff",
  "bcbio_umccrise", "protein_protein_contact_percent", "protein_protein_contact_pct_snpeff",
  "bcbio_umccrise", "sequence_feature", "sequence_feature_snpeff",
  "bcbio_umccrise", "sequence_feature_percent", "sequence_feature_pct_snpeff",
  "bcbio_umccrise", "splice_acceptor_variant", "splice_acceptor_variant_snpeff",
  "bcbio_umccrise", "splice_acceptor_variant_percent", "splice_acceptor_variant_pct_snpeff",
  "bcbio_umccrise", "splice_donor_variant", "splice_donor_variant_snpeff",
  "bcbio_umccrise", "splice_donor_variant_percent", "splice_donor_variant_pct_snpeff",
  "bcbio_umccrise", "splice_region_variant", "splice_region_variant_snpeff",
  "bcbio_umccrise", "splice_region_variant_percent", "splice_region_variant_pct_snpeff",
  "bcbio_umccrise", "start_lost", "start_lost_snpeff",
  "bcbio_umccrise", "start_lost_percent", "start_lost_pct_snpeff",
  "bcbio_umccrise", "stop_gained", "stop_gained_snpeff",
  "bcbio_umccrise", "stop_gained_percent", "stop_gained_pct_snpeff",
  "bcbio_umccrise", "stop_lost", "stop_lost_snpeff",
  "bcbio_umccrise", "stop_lost_percent", "stop_lost_pct_snpeff",
  "bcbio_umccrise", "structural_interaction_variant", "structural_interaction_variant_snpeff",
  "bcbio_umccrise", "structural_interaction_variant_percent", "structural_interaction_variant_pct_snpeff",
  "bcbio_umccrise", "synonymous_variant", "synonymous_variant_snpeff",
  "bcbio_umccrise", "synonymous_variant_percent", "synonymous_variant_pct_snpeff",
  "bcbio_umccrise", "upstream_gene_variant", "upstream_gene_variant_snpeff",
  "bcbio_umccrise", "upstream_gene_variant_percent", "upstream_gene_variant_pct_snpeff",
  "bcbio_umccrise", "DOWNSTREAM", "downstream_snpeff",
  "bcbio_umccrise", "DOWNSTREAM_percent", "downstream_pct_snpeff",
  "bcbio_umccrise", "EXON", "exon_snpeff",
  "bcbio_umccrise", "EXON_percent", "exon_pct_snpeff",
  "bcbio_umccrise", "INTERGENIC", "intergenic_snpeff",
  "bcbio_umccrise", "INTERGENIC_percent", "intergenic_pct_snpeff",
  "bcbio_umccrise", "INTRON", "intron_snpeff",
  "bcbio_umccrise", "INTRON_percent", "intron_pct_snpeff",
  "bcbio_umccrise", "MOTIF", "motif_snpeff",
  "bcbio_umccrise", "MOTIF_percent", "motif_pct_snpeff",
  "bcbio_umccrise", "SPLICE_SITE_ACCEPTOR", "splice_site_acceptor_snpeff",
  "bcbio_umccrise", "SPLICE_SITE_ACCEPTOR_percent", "splice_site_acceptor_pct_snpeff",
  "bcbio_umccrise", "SPLICE_SITE_DONOR", "splice_site_donor_snpeff",
  "bcbio_umccrise", "SPLICE_SITE_DONOR_percent", "splice_site_donor_pct_snpeff",
  "bcbio_umccrise", "SPLICE_SITE_REGION", "splice_site_region_snpeff",
  "bcbio_umccrise", "SPLICE_SITE_REGION_percent", "splice_site_region_pct_snpeff",
  "bcbio_umccrise", "TRANSCRIPT", "transcript_snpeff",
  "bcbio_umccrise", "TRANSCRIPT_percent", "transcript_pct_snpeff",
  "bcbio_umccrise", "UPSTREAM", "upstream_snpeff",
  "bcbio_umccrise", "UPSTREAM_percent", "upstream_pct_snpeff",
  "bcbio_umccrise", "UTR_3_PRIME", "utr_3_prime_snpeff",
  "bcbio_umccrise", "UTR_3_PRIME_percent", "utr_3_prime_pct_snpeff",
  "bcbio_umccrise", "UTR_5_PRIME", "utr_5_prime_snpeff",
  "bcbio_umccrise", "UTR_5_PRIME_percent", "utr_5_prime_pct_snpeff",
  "bcbio_umccrise", "Transitions", "transitions_snpeff",
  "bcbio_umccrise", "Transversions", "transversions_snpeff",
  "bcbio_umccrise", "Ts_Tv_ratio", "ts_tv_ratio_snpeff",
  "bcbio_umccrise", "Het", "het_snpeff",
  "bcbio_umccrise", "Hom", "hom_snpeff",
  "bcbio_umccrise", "Missing", "missing_snpeff",
  "bcbio_umccrise", "TFBS_ablation", "tfbs_ablation_snpeff",
  "bcbio_umccrise", "TFBS_ablation_percent", "tfbs_ablation_pct_snpeff",
  "bcbio_umccrise", "bidirectional_gene_fusion", "bidirectional_gene_fusion_snpeff",
  "bcbio_umccrise", "bidirectional_gene_fusion_percent", "bidirectional_gene_fusion_pct_snpeff",
  "bcbio_umccrise", "gene_fusion", "gene_fusion_snpeff",
  "bcbio_umccrise", "gene_fusion_percent", "gene_fusion_pct_snpeff",
  "bcbio_umccrise", "stop_retained_variant", "stop_retained_variant_snpeff",
  "bcbio_umccrise", "stop_retained_variant_percent", "stop_retained_variant_pct_snpeff",
  "bcbio_umccrise", "GENE", "gene_snpeff",
  "bcbio_umccrise", "GENE_percent", "gene_pct_snpeff",
  "bcbio_wts", "rRNA", "rrna_qualimap",
  "bcbio_wts", "rRNA_rate", "rrna_rate_qualimap",
  "bcbio_wts", "5'-3'_bias", "5_3_bias_qualimap1",
  "bcbio_wts", "5_3_bias", "5_3_bias_qualimap2",
  "bcbio_wts", "Intergenic_Rate", "intergenic_rate_qualimap",
  "bcbio_wts", "Intronic_Rate", "intronic_rate_qualimap",
  "bcbio_wts", "Exonic_Rate", "exonic_rate_qualimap",
  "bcbio_wts", "Duplication_Rate_of_Mapped", "mapped_dup_rate_qualimap",
  "bcbio_wts", "total_reads", "total_reads_star",
  "bcbio_wts", "avg_input_read_length", "avg_input_read_length_star",
  "bcbio_wts", "uniquely_mapped", "uniquely_mapped_star",
  "bcbio_wts", "uniquely_mapped_percent", "uniquely_mapped_percent_star",
  "bcbio_wts", "avg_mapped_read_length", "avg_mapped_read_length_star",
  "bcbio_wts", "num_splices", "num_splices_star",
  "bcbio_wts", "num_annotated_splices", "num_annotated_splices_star",
  "bcbio_wts", "num_GTAG_splices", "num_gtag_splices_star",
  "bcbio_wts", "num_GCAG_splices", "num_gcag_splices_star",
  "bcbio_wts", "num_ATAC_splices", "num_atac_splices_star",
  "bcbio_wts", "num_noncanonical_splices", "num_noncanonical_splices_star",
  "bcbio_wts", "mismatch_rate", "mismatch_rate_star",
  "bcbio_wts", "deletion_rate", "deletion_rate_star",
  "bcbio_wts", "deletion_length", "deletion_length_star",
  "bcbio_wts", "insertion_rate", "insertion_rate_star",
  "bcbio_wts", "insertion_length", "insertion_length_star",
  "bcbio_wts", "multimapped", "multimapped_star",
  "bcbio_wts", "multimapped_percent", "multimapped_percent_star",
  "bcbio_wts", "multimapped_toomany", "multimapped_toomany_star",
  "bcbio_wts", "multimapped_toomany_percent", "multimapped_toomany_percent_star",
  "bcbio_wts", "unmapped_mismatches_percent", "unmapped_mismatches_percent_star",
  "bcbio_wts", "unmapped_tooshort_percent", "unmapped_tooshort_percent_star",
  "bcbio_wts", "unmapped_other_percent", "unmapped_other_percent_star",
  "bcbio_wts", "unmapped_mismatches", "unmapped_mismatches_star",
  "bcbio_wts", "unmapped_tooshort", "unmapped_tooshort_star",
  "bcbio_wts", "unmapped_other", "unmapped_other_star",
  "bcbio_wts", "reads_aligned", "reads_aligned_qualimap",
  "bcbio_wgs", "number_of_samples", "num_samples_bcftools",
  "bcbio_wgs", "number_of_records", "num_records_bcftools",
  "bcbio_wgs", "number_of_no-ALTs", "num_noalts_bcftools",
  "bcbio_wgs", "number_of_SNPs", "num_snps_bcftools",
  "bcbio_wgs", "number_of_MNPs", "num_mnps_bcftools",
  "bcbio_wgs", "number_of_indels", "num_indels_bcftools",
  "bcbio_wgs", "number_of_others", "num_others_bcftools",
  "bcbio_wgs", "number_of_multiallelic_sites", "num_multiallelic_sites_bcftools",
  "bcbio_wgs", "number_of_multiallelic_SNP_sites", "num_multiallelic_snp_sites_bcftools",
  "bcbio_wgs", "ts", "ts_bcftools",
  "bcbio_wgs", "tv", "tv_bcftools",
  "bcbio_wgs", "tstv", "tstv_bcftools",
  "bcbio_wgs", "ts_1st_ALT", "ts_1st_alt_bcftools",
  "bcbio_wgs", "tv_1st_ALT", "tv_1st_alt_bcftools",
  "bcbio_wgs", "tstv_1st_ALT", "tstv_1st_alt_bcftools",
  "bcbio_wgs", "substitution_type_A>C", "subst_type_a_c_bcftools",
  "bcbio_wgs", "substitution_type_A>G", "subst_type_a_g_bcftools",
  "bcbio_wgs", "substitution_type_A>T", "subst_type_a_t_bcftools",
  "bcbio_wgs", "substitution_type_C>A", "subst_type_c_a_bcftools",
  "bcbio_wgs", "substitution_type_C>G", "subst_type_c_g_bcftools",
  "bcbio_wgs", "substitution_type_C>T", "subst_type_c_t_bcftools",
  "bcbio_wgs", "substitution_type_G>A", "subst_type_g_a_bcftools",
  "bcbio_wgs", "substitution_type_G>C", "subst_type_g_c_bcftools",
  "bcbio_wgs", "substitution_type_G>T", "subst_type_g_t_bcftools",
  "bcbio_wgs", "substitution_type_T>A", "subst_type_t_a_bcftools",
  "bcbio_wgs", "substitution_type_T>C", "subst_type_t_c_bcftools",
  "bcbio_wgs", "substitution_type_T>G", "subst_type_t_g_bcftools",
  "bcbio_wgs", "variations_hom", "variations_hom_bcftools",
  "bcbio_wgs", "variations_het", "variations_het_bcftools",
  "bcbio_wgs", "shared_hets_ped_check", "shared_hets_ped_check_peddy",
  "bcbio_wgs", "hets_a_ped_check", "hets_a_ped_check_peddy",
  "bcbio_wgs", "ibs2_ped_check", "ibs2_ped_check_peddy",
  "bcbio_wgs", "rel_ped_check", "rel_ped_check_peddy",
  "bcbio_wgs", "hets_b_ped_check", "hets_b_ped_check_peddy",
  "bcbio_wgs", "ibs0_ped_check", "ibs0_ped_check_peddy",
  "bcbio_wgs", "n_ped_check", "n_ped_check_peddy",
  "bcbio_wgs", "pedigree_parents_ped_check", "pedigree_parents_ped_check_peddy",
  "bcbio_wgs", "pedigree_relatedness_ped_check", "pedigree_relatedness_ped_check_peddy",
  "bcbio_wgs", "predicted_parents_ped_check", "predicted_parents_ped_check_peddy",
  "bcbio_wgs", "parent_error_ped_check", "parent_error_ped_check_peddy",
  "bcbio_wgs", "sample_duplication_error_ped_check", "sample_duplication_error_ped_check_peddy",
  "bcbio_wgs", "rel_difference_ped_check", "rel_difference_ped_check_peddy"
)
