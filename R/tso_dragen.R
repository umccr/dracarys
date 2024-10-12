#' Wf_dragen R6 Class
#'
#' @description
#' Reads and writes tidy versions of files from the `dragen` workflow.
#'
#' @examples
#' \dontrun{
#'
#' #---- Local ----#
#' p <- file.path(
#'   "~/s3/pipeline-prod-cache-503977275616-ap-southeast-2/byob-icav2/production",
#'   "analysis/cttsov2/20240915ff0295ed"
#' )
#' prefix <- "L2401290"
#' t1 <- Wf_tso_ctdna_tumor_only_v2$new(path = p, prefix = prefix)
#' t1$list_files(max_files = 100)
#' t1$list_files_filter_relevant(max_files = 300)
#' d <- t1$download_files(max_files = 100, dryrun = F)
#' d_tidy <- t1$tidy_files(d)
#' d_write <- t1$write(
#'   d_tidy,
#'   outdir = file.path(p, "dracarys_tidy"),
#'   prefix = prefix,
#'   format = "tsv"
#' )
#' }
#' @export
Wf_dragen <- R6::R6Class(
  "Wf_dragen",
  inherit = Wf,
  public = list(
    #' @field prefix The LibraryID prefix of the sample (needed for path lookup).
    prefix = NULL,
    #' @description Create a new Wf_dragen object.
    #' @param path Path to directory with raw workflow results (from S3 or
    #' local filesystem).
    #' @param prefix The LibraryID prefix of the sample (needed for path lookup).
    initialize = function(path = NULL, prefix = NULL) {
      wname <- "dragen"
      pref <- prefix
      reg1 <- tibble::tribble(
        ~regex, ~fun,
        glue("{dc}/{pref}\\-replay\\.json$"), "replay",
        glue("{dc}/{pref}\\.cnv_metrics.csv$"), "DOWNLOAD_ONLY",
        glue("{dc}/{pref}\\.exon_contig_mean_cov\\.csv$"), "contigMeanCov",
        glue("{dc}/{pref}\\.exon_coverage_metrics\\.csv$"), "DOWNLOAD_ONLY",
        glue("{dc}/{pref}\\.exon_fine_hist\\.csv$"), "DOWNLOAD_ONLY",
        glue("{dc}/{pref}\\.exon_hist\\.csv$"), "DOWNLOAD_ONLY",
        glue("{dc}/{pref}\\.exon_overall_mean_cov\\.csv$"), "DOWNLOAD_ONLY",
        glue("{dc}/{pref}\\.fastqc_metrics\\.csv$"), "DOWNLOAD_ONLY",
        glue("{dc}/{pref}\\.fragment_length_hist\\.csv$"), "DOWNLOAD_ONLY",
        glue("{dc}/{pref}\\.gc_metrics\\.csv$"), "DOWNLOAD_ONLY",
        glue("{dc}/{pref}\\.gvcf_metrics\\.csv$"), "DOWNLOAD_ONLY",
        glue("{dc}/{pref}\\.mapping_metrics\\.csv$"), "DOWNLOAD_ONLY",
        glue("{dc}/{pref}\\.microsat_diffs\\.txt$"), "DOWNLOAD_ONLY",
        glue("{dc}/{pref}\\.microsat_output\\.json$"), "DOWNLOAD_ONLY",
        glue("{dc}/{pref}\\.sv_metrics\\.csv$"), "DOWNLOAD_ONLY",
        glue("{dc}/{pref}\\.target_bed_contig_mean_cov\\.csv$"), "contigMeanCov",
        glue("{dc}/{pref}\\.target_bed_coverage_metrics\\.csv$"), "DOWNLOAD_ONLY",
        glue("{dc}/{pref}\\.target_bed_fine_hist\\.csv$"), "DOWNLOAD_ONLY",
        glue("{dc}/{pref}\\.target_bed_hist\\.csv$"), "DOWNLOAD_ONLY",
        glue("{dc}/{pref}\\.target_bed_overall_mean_cov\\.csv$"), "DOWNLOAD_ONLY",
        glue("{dc}/{pref}\\.time_metrics\\.csv$"), "DOWNLOAD_ONLY",
        glue("{dc}/{pref}\\.tmb_contig_mean_cov\\.csv$"), "contigMeanCov",
        glue("{dc}/{pref}\\.tmb_coverage_metrics\\.csv$"), "DOWNLOAD_ONLY",
        glue("{dc}/{pref}\\.tmb_fine_hist\\.csv$"), "DOWNLOAD_ONLY",
        glue("{dc}/{pref}\\.tmb_hist\\.csv$"), "DOWNLOAD_ONLY",
        glue("{dc}/{pref}\\.tmb_overall_mean_cov\\.csv$"), "DOWNLOAD_ONLY",
        glue("{dc}/{pref}\\.trimmer_metrics\\.csv$"), "DOWNLOAD_ONLY",
        glue("{dc}/{pref}\\.umi_metrics\\.csv$"), "DOWNLOAD_ONLY",
        glue("{dc}/{pref}\\.vc_metrics\\.csv$"), "DOWNLOAD_ONLY",
        glue("{dc}/{pref}\\.wgs_contig_mean_cov\\.csv$"), "contigMeanCov",
        glue("{dc}/{pref}\\.wgs_coverage_metrics\\.csv$"), "DOWNLOAD_ONLY",
        glue("{dc}/{pref}\\.wgs_fine_hist\\.csv$"), "DOWNLOAD_ONLY",
        glue("{dc}/{pref}\\.wgs_hist\\.csv$"), "DOWNLOAD_ONLY",
        glue("{dc}/{pref}\\.wgs_overall_mean_cov\\.csv$"), "DOWNLOAD_ONLY"
      )
      regexes <- reg1 |>
        dplyr::mutate(
          fun = paste0("read_", .data$fun),
          fun = ifelse(.data$fun == "read_DOWNLOAD_ONLY", "DOWNLOAD_ONLY", .data$fun)
        )

      super$initialize(path = path, wname = wname, regexes = regexes)
      self$prefix <- prefix
    },
    #' @description Print details about the Workflow.
    #' @param ... (ignored).
    print = function(...) {
      res <- tibble::tribble(
        ~var, ~value,
        "path", self$path,
        "wname", self$wname,
        "filesystem", self$filesystem,
        "prefix", self$prefix
      )
      print(res)
      invisible(self)
    },
    #' @description Read `replay.json` file.
    #' @param x Path to file.
    read_replay = function(x) {
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
      return(dplyr::bind_cols(res))
    },
    #' @description Read `contig_mean_cov.csv` file.
    #' @param x Path to file.
    read_contigMeanCov = function(x) {
      readr::read_csv(x, col_names = c("chrom", "n_bases", "coverage"), col_types = "cdd") |>
        dplyr::filter(
          if (!keep_alt) {
            !grepl("chrM|MT|_|Autosomal|HLA-|EBV", .data$chrom)
          } else {
            TRUE
          }
        )
    },
    #' @description Read `dragen.tsv.gz` cancer report hrd file.
    #' @param x Path to file.
    read_coverageMetrics = function(x) {
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
      return(res)
    },
  ) # end public
) # end Wf_dragen
