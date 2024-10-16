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
#' d1 <- Wf_dragen$new(path = p, prefix = prefix)
#' d1$list_files(max_files = 100)
#' d1$list_files_filter_relevant(max_files = 300)
#' d <- d1$download_files(max_files = 100, dryrun = F)
#' d_tidy <- d1$tidy_files(d)
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
        glue("{pref}\\.target_bed_contig_mean_cov\\.csv$"), "contigMeanCov",
        glue("{pref}\\.tmb_contig_mean_cov\\.csv$"), "contigMeanCov",
        glue("{pref}\\.wgs_contig_mean_cov\\.csv$"), "contigMeanCov",
        glue("{pref}\\.exon_coverage_metrics\\.csv$"), "coverageMetrics",
        glue("{pref}\\.target_bed_coverage_metrics\\.csv$"), "coverageMetrics",
        glue("{pref}\\.tmb_coverage_metrics\\.csv$"), "coverageMetrics",
        glue("{pref}\\.wgs_coverage_metrics\\.csv$"), "coverageMetrics",
        glue("{pref}\\.exon_fine_hist\\.csv$"), "fineHist",
        glue("{pref}\\.target_bed_fine_hist\\.csv$"), "fineHist",
        glue("{pref}\\.tmb_fine_hist\\.csv$"), "fineHist",
        glue("{pref}\\.wgs_fine_hist\\.csv$"), "fineHist",
        glue("{pref}\\.exon_hist\\.csv$"), "hist",
        glue("{pref}\\.target_bed_hist\\.csv$"), "hist",
        glue("{pref}\\.tmb_hist\\.csv$"), "hist",
        glue("{pref}\\.wgs_hist\\.csv$"), "hist",
        glue("{pref}\\.fastqc_metrics\\.csv$"), "DOWNLOAD_ONLY",
        glue("{pref}\\.fragment_length_hist\\.csv$"), "fragmentLengthHist",
        glue("{pref}\\.gc_metrics\\.csv$"), "DOWNLOAD_ONLY",
        glue("{pref}\\.gvcf_metrics\\.csv$"), "DOWNLOAD_ONLY",
        glue("{pref}\\.mapping_metrics\\.csv$"), "DOWNLOAD_ONLY",
        glue("{pref}\\.microsat_diffs\\.txt$"), "DOWNLOAD_ONLY",
        glue("{pref}\\.microsat_output\\.json$"), "DOWNLOAD_ONLY",
        glue("{pref}\\.sv_metrics\\.csv$"), "DOWNLOAD_ONLY",
        glue("{pref}\\.time_metrics\\.csv$"), "DOWNLOAD_ONLY",
        glue("{pref}\\.trimmer_metrics\\.csv$"), "DOWNLOAD_ONLY",
        glue("{pref}\\.umi_metrics\\.csv$"), "DOWNLOAD_ONLY",
        glue("{pref}\\.vc_metrics\\.csv$"), "DOWNLOAD_ONLY"
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
      dat <- dplyr::bind_cols(res)
      tibble::tibble(name = "replay", data = list(dat))
    },
    #' @description Read `contig_mean_cov.csv` file.
    #' @param x Path to file.
    #' @param keep_alt Keep ALT contigs.
    read_contigMeanCov = function(x, keep_alt = FALSE) {
      subprefix <- dragen_subprefix(x, "_contig_mean_cov")
      dat <- readr::read_csv(x, col_names = c("chrom", "n_bases", "coverage"), col_types = "cdd") |>
        dplyr::filter(
          if (!keep_alt) {
            !grepl("chrM|MT|_|Autosomal|HLA-|EBV|GL|hs37d5", .data$chrom)
          } else {
            TRUE
          }
        )
      tibble::tibble(name = glue("contigmeancov_{subprefix}"), data = list(dat[]))
    },
    #' @description Read `coverage_metrics.csv` file.
    #' @param x Path to file.
    read_coverageMetrics = function(x) {
      subprefix <- dragen_subprefix(x, "_coverage_metrics")
      dat <- dragen_coverage_metrics_read(x)
      tibble::tibble(name = glue("covmetrics_{subprefix}"), data = list(dat))
    },
    #' @description Read `fine_hist.csv` file.
    #' @param x Path to file.
    read_fineHist = function(x) {
      subprefix <- dragen_subprefix(x, "_fine_hist")
      d <- readr::read_csv(x, col_types = "cd")
      assertthat::assert_that(all(colnames(d) == c("Depth", "Overall")))
      # there's a max Depth of 2000+, so convert to numeric for easier plotting
      dat <- d |>
        dplyr::mutate(
          Depth = ifelse(grepl("+", .data$Depth), sub("(\\d*)\\+", "\\1", .data$Depth), .data$Depth),
          Depth = as.integer(.data$Depth)
        ) |>
        dplyr::select(depth = "Depth", n_loci = "Overall")
      tibble::tibble(name = glue("finehist_{subprefix}"), data = list(dat))
    },
    #' @description Read `fragment_length_hist.csv` file.
    #' @param x Path to file.
    read_fragmentLengthHist = function(x) {
      d <- readr::read_lines(x)
      assertthat::assert_that(grepl("#Sample", d[1]))
      dat <- d |>
        tibble::enframe(name = "name", value = "value") |>
        dplyr::filter(!grepl("#Sample: |FragmentLength,Count", .data$value)) |>
        tidyr::separate_wider_delim(cols = "value", names = c("fragmentLength", "count"), delim = ",") |>
        dplyr::mutate(
          count = as.numeric(.data$count),
          fragmentLength = as.numeric(.data$fragmentLength)
        ) |>
        dplyr::select("fragmentLength", "count")
      tibble::tibble(name = "fraglen", data = list(dat))
    },
    #' @description Read `mapping_metrics.csv` file.
    #' @param x Path to file.
    read_mappingMetrics = function(x) {
      dat <- dragen_mapping_metrics_read(x)
      tibble::tibble(name = "mapmetrics", data = list(dat))
    },
    #' @description Read `hist.csv` file.
    read_hist = function(x) {
      subprefix <- dragen_subprefix(x, "_hist")
      d <- readr::read_csv(x, col_names = c("var", "pct"), col_types = "cd")
      dat <- d |>
        dplyr::mutate(
          var = sub("PCT of bases in .* with coverage ", "", .data$var),
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
      tibble::tibble(name = glue("hist_{subprefix}"), data = list(dat))
    }
  ) # end public
) # end Wf_dragen
