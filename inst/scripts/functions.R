library(dplyr, include.only = c("bind_rows", "left_join"))
library(purrr, include.only = c("map", "map_chr", "reduce", "list_merge", "imap", "set_names"))

select_column_subset_alignmentqc <- function(d) {
  cols_to_keep <- c(
    "umccr_subj_id"           = "umccr_subj_id",
    "tot_input_reads"         = "dragen_map_metrics.Total input reads",
    "tot_input_reads_pct"     = "dragen_map_metrics.Total input reads pct",
    "num_dup_reads"           = "dragen_map_metrics.Number of duplicate marked reads",
    "num_dup_reads_pct"       = "dragen_map_metrics.Number of duplicate marked reads pct",
    "num_uniq_reads_pct"      = "dragen_map_metrics.Number of unique reads (excl. duplicate marked reads) pct",
    "mapped_reads"            = "dragen_map_metrics.Mapped reads",
    "mapped_reads_pct"        = "dragen_map_metrics.Mapped reads pct",
    "unmapped_reads"          = "dragen_map_metrics.Unmapped reads",
    "unmapped_reads_pct"      = "dragen_map_metrics.Unmapped reads pct",
    "singleton_reads"         = "dragen_map_metrics.Singleton reads (itself mapped; mate unmapped)",
    "singleton_reads_pct"     = "dragen_map_metrics.Singleton reads (itself mapped; mate unmapped) pct",
    "paired_reads"            = "dragen_map_metrics.Paired reads (itself & mate mapped)",
    "paired_reads_pct"        = "dragen_map_metrics.Paired reads (itself & mate mapped) pct",
    "paired_reads_proper"     = "dragen_map_metrics.Properly paired reads",
    "paired_reads_proper_pct" = "dragen_map_metrics.Properly paired reads pct",
    "read_length"             = "dragen_map_metrics.Estimated read length",
    "insert_length_mean"      = "dragen_map_metrics.Insert length: mean",
    "insert_length_median"    = "dragen_map_metrics.Insert length: median",
    "contamination"           = "dragen_map_metrics.Estimated sample contamination",
    "cov_seqed_avg_genome"    = "dragen_map_metrics.Average sequenced coverage over genome",
    "cov_aligned_avg_genome"  = "dragen_cov_metrics.Average alignment coverage over genome",
    "cov_autosomal_median"    = "dragen_ploidy.Autosomal median coverage",
    "cov_chrx_median"         = "dragen_ploidy.X median coverage",
    "cov_chry_median"         = "dragen_ploidy.Y median coverage",
    "cov50_genome_pct"        = "dragen_cov_metrics.PCT of genome with coverage [ 50x: inf)",
    "ploidy_est"              = "dragen_ploidy.Ploidy estimation"
  )

  assertthat::assert_that(base::all(cols_to_keep %in% base::names(d)))

  d |>
    dplyr::select(dplyr::all_of(cols_to_keep)) |>
    purrr::set_names(base::names(cols_to_keep))
}
