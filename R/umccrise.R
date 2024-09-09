#' Wf_umccrise R6 Class
#'
#' @description
#' Reads and writes tidy versions of files from the `umccrise` workflow
#'
#' @examples
#' \dontrun{
#'
#' #---- LOCAL ----#
#' SubjectID <- "SBJ03043"
#' SampleID_tumor <- "PRJ230004"
#' p1_local <- "~/icav1/g/production/analysis_data"
#' p <- file.path(p1_local, "SBJ03043/umccrise/20240830ec648f40/L2300064__L2300063")
#' um1 <- Wf_umccrise$new(path = p, SubjectID = SubjectID, SampleID_tumor = SampleID_tumor)
#' um1$list_files(max_files = 10)
#' um1$list_files_filter_relevant()
#' d <- um1$download_files(max_files = 1000, dryrun = F)
#' d_tidy <- um1$tidy_files(d)
#' d_write <- um1$write(
#'   d_tidy,
#'   outdir = file.path(p, "dracarys_tidy"),
#'   prefix = glue("{SubjectID}__{SampleID_tumor}"),
#'   format = "tsv"
#' )
#'
#' #---- GDS ----#
#' SubjectID <- "SBJ03043"
#' SampleID_tumor <- "PRJ230004"
#' p1_gds <- "gds://production/analysis_data"
#' p <- file.path(p1_gds, "SBJ03043/umccrise/20240830ec648f40/L2300064__L2300063")
#' outdir <- file.path(sub("gds:/", "~/icav1/g", p))
#' token <- Sys.getenv("ICA_ACCESS_TOKEN")
#' um2 <- Wf_umccrise$new(path = p, SubjectID = SubjectID, SampleID_tumor = SampleID_tumor)
#' um2$list_files(max_files = 8)
#' um2$list_files_filter_relevant(ica_token = token, max_files = 500)
#' d <- um2$download_files(
#'   outdir = outdir, ica_token = token,
#'   max_files = 1000, dryrun = F
#' )
#' d_tidy <- um2$tidy_files(d)
#' }
#'
#' @export
Wf_umccrise <- R6::R6Class(
  "Wf_umccrise",
  inherit = Wf,
  public = list(
    #' @field SubjectID The SubjectID of the sample (needed for path lookup).
    #' @field SampleID_tumor The SampleID of the tumor sample (needed for path lookup).
    SubjectID = NULL,
    SampleID_tumor = NULL,
    #' @description Create a new Wf_umccrise object.
    #' @param path Output directory path with results.
    #' @param SubjectID The SubjectID of the sample (needed for path lookup).
    #' @param SampleID_tumor The SampleID of the tumor sample (needed for path lookup).
    initialize = function(path = NULL, SubjectID = NULL, SampleID_tumor = NULL) {
      wname <- "umccrise"
      regexes <- tibble::tribble(
        ~regex, ~fun,
        "-chord\\.tsv\\.gz$", "chordtsv",
        "-hrdetect\\.tsv\\.gz$", "hrdetecttsv",
        "-snv_2015\\.tsv\\.gz$", "sigssnv2015tsv",
        "-snv_2020\\.tsv\\.gz$", "sigssnv2020tsv",
        "-dbs\\.tsv\\.gz$", "sigsdbstsv",
        "-indel\\.tsv\\.gz$", "sigsindeltsv",
        "-qc_summary\\.tsv\\.gz$", "qcsummarytsv",
        "multiqc_conpair.txt", "conpairmultiqc",
        "-somatic\\.pcgr\\.json\\.gz$", "pcgrjson"
      ) |>
        dplyr::mutate(fun = paste0("read_", .data$fun))

      super$initialize(path = path, wname = wname, regexes = regexes)
      self$SubjectID <- SubjectID
      self$SampleID_tumor <- SampleID_tumor
    },
    #' @description Print details about the Workflow.
    #' @param ... (ignored).
    print = function(...) {
      res <- tibble::tribble(
        ~var, ~value,
        "path", self$path,
        "wname", self$wname,
        "filesystem", self$filesystem,
        "SubjectID", self$SubjectID,
        "SampleID_tumor", self$SampleID_tumor
      )
      print(res)
      invisible(self)
    },
    #' @description List dracarys files under given path
    #' @param max_files Max number of files to list (for gds/s3 only).
    #' @param ica_token ICA access token (def: $ICA_ACCESS_TOKEN env var).
    #' @param ... Passed on to the `gds_list_files_filter_relevant` or
    #' the `s3_list_files_filter_relevant` function.
    list_files_filter_relevant = function(max_files = 1000,
                                          ica_token = Sys.getenv("ICA_ACCESS_TOKEN"), ...) {
      path <- self$path
      dir_final <- file.path(path, glue("{self$SubjectID}__{self$SampleID_tumor}"))
      dir_work <- file.path(path, "work", glue("{self$SubjectID}__{self$SampleID_tumor}"))
      dir_work_pcgr <- file.path(dir_work, "pcgr") # for pcgr json
      f1 <- super$list_files_filter_relevant(path = dir_final, max_files = 300, ica_token = ica_token)
      f2 <- super$list_files_filter_relevant(path = dir_work_pcgr, max_files = 50, ica_token = ica_token)
      f_all <- dplyr::bind_rows(f1, f2)
      return(f_all)
    },
    #' @description Download files from GDS/S3 to local filesystem.
    #' @param outdir Path to output directory.
    #' @param ica_token ICA access token (def: $ICA_ACCESS_TOKEN env var).
    #' @param max_files Max number of files to list.
    #' @param dryrun If TRUE, just list the files that will be downloaded (don't
    #' download them).
    #' @param recursive Should files be returned recursively _in and under_ the specified
    #' GDS directory, or _only directly in_ the specified GDS directory (def: TRUE via ICA API).
    #' @param list_filter_fun Function to filter relevant files.
    download_files = function(outdir, ica_token = Sys.getenv("ICA_ACCESS_TOKEN"),
                              max_files = 1000, dryrun = FALSE, recursive = NULL) {
      super$download_files(
        outdir = outdir, ica_token = ica_token, max_files = max_files,
        dryrun = dryrun, recursive = recursive,
        list_filter_fun = self$list_files_filter_relevant
      )
    },
    #' @description Read `pcgr.json.gz` file.
    #' @param x Path to file.
    read_pcgrjson = function(x) {
      j <- read_jsongz_jsonlite(x)
      tmb <-
        j[["content"]][["tmb"]][["variant_statistic"]] %||%
        j[["content"]][["tmb"]][["v_stat"]] %||%
        list(tmb_estimate = NA, n_tmb = NA)
      tmb <- purrr::flatten(tmb) |>
        tibble::as_tibble_row() |>
        dplyr::select("tmb_estimate", "n_tmb")
      msi <- j[["content"]][["msi"]][["prediction"]][["msi_stats"]]
      # handle nulls
      msi <- msi %||% list(fracIndels = NA, predicted_class = NA)
      msi <- purrr::flatten(msi) |>
        tibble::as_tibble_row() |>
        dplyr::select("fracIndels", "predicted_class")
      metrics <- dplyr::bind_cols(msi, tmb)
      return(metrics)
    },
    #' @description Read `chord.tsv.gz` cancer report file.
    #' @param x Path to file.
    read_chordtsv = function(x) {
      ct <- readr::cols_only(
        p_hrd = "d",
        hr_status = "c",
        hrd_type = "c",
        p_BRCA1 = "d",
        p_BRCA2 = "d"
      )
      read_tsvgz(x, col_types = ct)
    },
    #' @description Read `hrdetect.tsv.gz` cancer report file.
    #' @param x Path to file.
    read_hrdetecttsv = function(x) {
      ct <- readr::cols(
        .default = "d",
        sample = "c"
      )
      read_tsvgz(x, col_types = ct) |>
        dplyr::select(-c("sample"))
    },
    #' @description Read signature cancer report file.
    #' @param x Path to file.
    read_sigstsv = function(x) {
      ct <- readr::cols(
        .default = "d",
        Signature = "c"
      )
      read_tsvgz(x, col_types = ct)
    },
    #' @description Read `snv_2015.tsv.gz` sigs cancer report file.
    #' @param x Path to file.
    read_sigssnv2015tsv = function(x) {
      self$read_sigstsv(x)
    },
    #' @description Read `snv_2020.tsv.gz` sigs cancer report file.
    #' @param x Path to file.
    read_sigssnv2020tsv = function(x) {
      self$read_sigstsv(x)
    },
    #' @description Read `dbs.tsv.gz` sigs cancer report file.
    #' @param x Path to file.
    read_sigsdbstsv = function(x) {
      self$read_sigstsv(x)
    },
    #' @description Read `indel.tsv.gz` sigs cancer report file.
    #' @param x Path to file.
    read_sigsindeltsv = function(x) {
      self$read_sigstsv(x)
    },
    #' @description Read `qc_summary.tsv.gz` cancer report file.
    #' @param x Path to file.
    read_qcsummarytsv = function(x) {
      d <- read_tsvgz(x, col_types = readr::cols(.default = "c"))
      d |>
        dplyr::select("variable", "value") |>
        tidyr::pivot_wider(names_from = "variable", values_from = "value") |>
        dplyr::rename(MSI_mb_tmp = "MSI (indels/Mb)") |>
        dplyr::mutate(
          purity_hmf = sub("(.*) \\(.*\\)", "\\1", .data$Purity) |> as.numeric(),
          ploidy_hmf = sub("(.*) \\(.*\\)", "\\1", .data$Ploidy) |> as.numeric(),
          hrd_chord = sub("CHORD: (.*); HRDetect: (.*)", "\\1", .data$HRD) |> as.numeric(),
          hrd_hrdetect = sub("CHORD: (.*); HRDetect: (.*)", "\\2", .data$HRD),
          # handle HRDetect NA
          hrd_hrdetect = ifelse(.data$hrd_hrdetect == "NA", NA_real_, as.numeric(.data$hrd_hrdetect)),
          msi_mb_hmf = sub(".* \\((.*)\\)", "\\1", .data$MSI_mb_tmp) |> as.numeric(),
          contamination_hmf = as.numeric(.data$Contamination),
          deleted_genes_hmf = as.numeric(.data$DeletedGenes),
          msi_hmf = sub("(.*) \\(.*\\)", "\\1", .data$MSI_mb_tmp),
          tmb_hmf = sub("(.*) \\(.*\\)", "\\1", .data$TMB) |> as.numeric(),
          tml_hmf = sub("(.*) \\(.*\\)", "\\1", .data$TML) |> as.numeric(),
          hypermutated = ifelse("Hypermutated" %in% d$variable, .data[["Hypermutated"]], NA) |> as.character(),
          bpi_enabled = ifelse("BPI Enabled" %in% d$variable, .data[["BPI Enabled"]], NA) |> as.character(),
        ) |>
        dplyr::select(
          qc_status_hmf = "QC_Status",
          sex_hmf = "Gender",
          "purity_hmf", "ploidy_hmf", "msi_hmf", "msi_mb_hmf",
          "hrd_chord", "hrd_hrdetect", "contamination_hmf",
          "deleted_genes_hmf", "tmb_hmf", "tml_hmf",
          wgd_hmf = "WGD",
          "hypermutated", "bpi_enabled"
        )
    },
    #' @description Read multiqc_conpair.txt file.
    #' @param x Path to file.
    read_conpairmultiqc = function(x) {
      um_ref_samples <- c("Alice", "Bob", "Chen", "Elon", "Dakota")
      um_ref_samples <- paste0(um_ref_samples, rep(c("_T", "_B", ""), each = length(um_ref_samples)))
      cnames <- list(
        old = c(
          "Sample", "concordance_concordance", "concordance_used_markers",
          "concordance_total_markers", "concordance_marker_threshold",
          "concordance_min_mapping_quality", "concordance_min_base_quality",
          "contamination"
        ),
        new = c(
          "sampleid", "contamination", "concordance", "markers_used",
          "markers_total", "marker_threshold",
          "mapq_min", "baseq_min"
        )
      )
      ctypes <- list(
        old = c("cddddddd"),
        new = c("cddddddd")
      )
      if (!file.exists(x)) {
        return(empty_tbl(cnames$new, ctypes$new))
      }
      d1 <- readr::read_tsv(x, col_types = readr::cols(.default = "d", Sample = "c"))
      assertthat::assert_that(all(colnames(d1) == cnames$old))
      d1 |>
        dplyr::filter(!.data$Sample %in% um_ref_samples) |>
        dplyr::relocate("contamination", .after = "Sample") |>
        rlang::set_names(cnames$new)
    }
  ) # end public
)

#    read = function() {
#      x <- self$path
#      # now return all as list elements
#      list(
#        chord = grep_file(x, "-chord\\.tsv\\.gz$") |> self$read_chordtsv(),
#        hrdetect = grep_file(x, "-hrdetect\\.tsv\\.gz$") |> self$read_hrdetecttsv(),
#        sigs2015 = grep_file(x, "-snv_2015\\.tsv\\.gz$") |> self$read_sigs(),
#        sigs2020 = grep_file(x, "-snv_2020\\.tsv\\.gz$") |> self$read_sigs(),
#        sigsdbs = grep_file(x, "-dbs\\.tsv\\.gz$") |> self$read_sigs(),
#        sigsindel = grep_file(x, "-indel\\.tsv\\.gz$") |> self$read_sigs(),
#        qcsum = grep_file(x, "-qc_summary\\.tsv\\.gz$") |> self$read_qcsummarytsv()
#      )
#    }
