#' Wf_dragen R6 Class
#'
#' @description
#' Reads and writes tidy versions of files from the `dragen` workflow.
#'
#' @examples
#' \dontrun{
#'
#' #---- Local ----#
#' prefix <- "L2401290"
#' p <- file.path(
#'   "~/s3/pipeline-prod-cache-503977275616-ap-southeast-2/byob-icav2/production",
#'   "analysis/cttsov2/20240915ff0295ed/Logs_Intermediates/DragenCaller",
#'   prefix
#' )
#' t1 <- Wf_dragen$new(path = p, prefix = prefix)
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
        glue("{pref}\\-replay\\.json$"), "replay",
        glue("{pref}\\.cnv_metrics.csv$"), "DOWNLOAD_ONLY",
        glue("{pref}\\.exon_contig_mean_cov\\.csv$"), "contigMeanCov",
        glue("{pref}\\.exon_coverage_metrics\\.csv$"), "DOWNLOAD_ONLY",
        glue("{pref}\\.exon_fine_hist\\.csv$"), "DOWNLOAD_ONLY",
        glue("{pref}\\.exon_hist\\.csv$"), "DOWNLOAD_ONLY",
        glue("{pref}\\.exon_overall_mean_cov\\.csv$"), "DOWNLOAD_ONLY",
        glue("{pref}\\.fastqc_metrics\\.csv$"), "DOWNLOAD_ONLY",
        glue("{pref}\\.fragment_length_hist\\.csv$"), "DOWNLOAD_ONLY",
        glue("{pref}\\.gc_metrics\\.csv$"), "DOWNLOAD_ONLY",
        glue("{pref}\\.gvcf_metrics\\.csv$"), "DOWNLOAD_ONLY",
        glue("{pref}\\.mapping_metrics\\.csv$"), "DOWNLOAD_ONLY",
        glue("{pref}\\.microsat_diffs\\.txt$"), "DOWNLOAD_ONLY",
        glue("{pref}\\.microsat_output\\.json$"), "DOWNLOAD_ONLY",
        glue("{pref}\\.sv_metrics\\.csv$"), "DOWNLOAD_ONLY",
        glue("{pref}\\.target_bed_contig_mean_cov\\.csv$"), "contigMeanCov",
        glue("{pref}\\.target_bed_coverage_metrics\\.csv$"), "DOWNLOAD_ONLY",
        glue("{pref}\\.target_bed_fine_hist\\.csv$"), "DOWNLOAD_ONLY",
        glue("{pref}\\.target_bed_hist\\.csv$"), "DOWNLOAD_ONLY",
        glue("{pref}\\.target_bed_overall_mean_cov\\.csv$"), "DOWNLOAD_ONLY",
        glue("{pref}\\.time_metrics\\.csv$"), "DOWNLOAD_ONLY",
        glue("{pref}\\.tmb_contig_mean_cov\\.csv$"), "contigMeanCov",
        glue("{pref}\\.tmb_coverage_metrics\\.csv$"), "DOWNLOAD_ONLY",
        glue("{pref}\\.tmb_fine_hist\\.csv$"), "DOWNLOAD_ONLY",
        glue("{pref}\\.tmb_hist\\.csv$"), "DOWNLOAD_ONLY",
        glue("{pref}\\.tmb_overall_mean_cov\\.csv$"), "DOWNLOAD_ONLY",
        glue("{pref}\\.trimmer_metrics\\.csv$"), "DOWNLOAD_ONLY",
        glue("{pref}\\.umi_metrics\\.csv$"), "DOWNLOAD_ONLY",
        glue("{pref}\\.vc_metrics\\.csv$"), "DOWNLOAD_ONLY",
        glue("{pref}\\.wgs_contig_mean_cov\\.csv$"), "contigMeanCov",
        glue("{pref}\\.wgs_coverage_metrics\\.csv$"), "DOWNLOAD_ONLY",
        glue("{pref}\\.wgs_fine_hist\\.csv$"), "DOWNLOAD_ONLY",
        glue("{pref}\\.wgs_hist\\.csv$"), "DOWNLOAD_ONLY",
        glue("{pref}\\.wgs_overall_mean_cov\\.csv$"), "DOWNLOAD_ONLY"
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
    #' @param keep_alt Keep ALT contigs.
    read_contigMeanCov = function(x, keep_alt = FALSE) {
      readr::read_csv(x, col_names = c("chrom", "n_bases", "coverage"), col_types = "cdd") |>
        dplyr::filter(
          if (!keep_alt) {
            !grepl("chrM|MT|_|Autosomal|HLA-|EBV|GL|hs37d5", .data$chrom)
          } else {
            TRUE
          }
        )
    },
    #' @description Read `dragen.tsv.gz` cancer report hrd file.
    #' @param x Path to file.
    read_coverageMetrics = function(x) {
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
        dplyr::mutate(var = as.character(glue("cov_pct_{start}_{end}_{region}_dragen"))) |>
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
  ) # end public
) # end Wf_dragen
