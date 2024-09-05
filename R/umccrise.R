#' Wf_umccrise R6 Class
#'
#' @description
#' Reads and writes tidy versions of files from the `umccrise` workflow
#'
#' @examples
#' \dontrun{
#' token <- Sys.getenv("ICA_ACCESS_TOKEN") |> ica_token_validate()
#' SubjectID <- "SBJ01155"
#' SampleID_tumor <- "PRJ211091"
#' gdsdir1 <- "gds://production/analysis_data/SBJ01155/umccrise/202408300c218043"
#' gdsdir <- file.path(gdsdir1, "L2101566__L2101565")
#' obj <- Wf_umccrise$new(gdsdir)
#' gds_files <- obj$gds_list(
#'   gdsdir = gdsdir, token = token, SubjectID = SubjectID, SampleID_tumor
#' )
#' outdir <- file.path(sub("gds://", "", gdsdir))
#' outdir <- file.path(normalizePath("~/icav1/g"), outdir)
#' out_files <- obj$gds_download(gds_files = gds_files, outdir = outdir, token = token)
#' tidy1 <- obj$tidy(indir = outdir, out_files = out_files)
#' }
#'
#' @export
Wf_umccrise <- R6::R6Class(
  "Wf_umccrise",
  public = list(
    #' @field path Path to the `umccrise` directory.
    #' @field contents Tibble with file path, basename, and size.
    path = NULL,
    contents = NULL,
    #' @description Create a new Wf_umccrise object.
    #' @param path Path to the `umccrise` directory.
    initialize = function(path = NULL) {
      stopifnot(is.character(path), length(path) == 1)
      self$path <- path
    },
    #' @description List Relevant Files In umccrise GDS Directory
    #' @param gdsdir Path to the `umccrise` directory.
    #' @param SubjectID The SubjectID of the sample (used to construct path).
    #' @param SampleID_tumor The SampleID of the tumor sample (used to construct path).
    #' @param token ICA access token.
    gds_list = function(gdsdir, SubjectID, SampleID_tumor, token = Sys.getenv("ICA_ACCESS_TOKEN")) {
      reg1 <- tibble::tribble(
        ~regex, ~fun,
        # "-somatic\\.pcgr\\.snvs_indels\\.tiers\\.tsv$", "PcgrTiersFile",
        "-chord\\.tsv\\.gz$", "UmChordTsvFile",
        "-hrdetect\\.tsv\\.gz$", "UmHrdetectTsvFile",
        "-snv_2015\\.tsv\\.gz$", "UmSigsSnvFile",
        "-snv_2020\\.tsv\\.gz$", "UmSigsSnvFile",
        "-dbs\\.tsv\\.gz$", "UmSigsDbsFile",
        "-indel\\.tsv\\.gz$", "UmSigsIndelFile",
        "-qc_summary\\.tsv\\.gz$", "UmQcSumFile",
        "multiqc_conpair.txt", "UmConpairMultiqc"
      )
      reg2 <- tibble::tribble(
        ~regex, ~fun,
        "-somatic\\.pcgr\\.json\\.gz$", "PcgrJsonFile"
      )
      dir_fin <- file.path(gdsdir, glue("{SubjectID}__{SampleID_tumor}"))
      dir_wrk <- file.path(gdsdir, "work", glue("{SubjectID}__{SampleID_tumor}"))
      dir_wrk_pcgr <- file.path(dir_wrk, "pcgr") # for pcgr json
      f1 <- gds_files_list_filter_relevant(gdsdir = dir_fin, token, page_size = 300, regexes = reg1)
      f2 <- gds_files_list_filter_relevant(gdsdir = dir_wrk_pcgr, token, page_size = 50, regexes = reg2)
      gds_files <- dplyr::bind_rows(f1, f2)
      return(gds_files)
    },

    #' @description GDS File Download via API
    #'
    #' @param gds_files Tibble with bname and file_id for umccrise files.
    #' @param outdir Directory to output files (loosely, not in a structured manner).
    #' @param token ICA access token.
    gds_download = function(gds_files, outdir, token = Sys.getenv("ICA_ACCESS_TOKEN")) {
      assertthat::assert_that(all(c("bname", "file_id") %in% colnames(gds_files)))
      gds_files |>
        dplyr::rowwise() |>
        dplyr::mutate(
          out = file.path(outdir, .data$bname),
          out_dl = gds_file_download_api(.data$file_id, .data$out, token)
        )
    },

    #' @description Tidy up the output files from umccrise
    #'
    #' @param indir Path to the `umccrise` directory.
    #' @param out_files Tibble with file path, basename, and size.
    tidy = function(indir, out_files) {
      obj_canrep <- UmccriseCanRepTables$new(indir)
      canrep_parse <- obj_canrep$read()
      pcgr_json <- out_files |>
        dplyr::filter(.data$type == "PcgrJsonFile") |>
        dplyr::pull("out") |>
        PcgrJsonFile$new() |>
        read()
      conpair_tsv <- out_files |>
        dplyr::filter(.data$type == "UmConpairMultiqc") |>
        dplyr::pull("out") |>
        self$read_conpairmultiqc()
      d <- canrep_parse
      d[["pcgr_json"]] <- pcgr_json[["metrics"]]
      d[["conpair"]] <- conpair_tsv
      d
    },

    #' @description Read multiqc_conpair.txt file.
    #'
    #' @param x (`character(1)`)\cr
    #'   Path to multiqc_conpair.txt file.
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
  )
)

#' UmccriseCanRepTables R6 Class
#'
#' @description
#' Reads and writes tidy versions of files within the `cancer_report_tables` directory
#' output from the `umccrise` workflow.
#'
#' @examples
#' \dontrun{
#' p1 <- "~/icav1/g/production/analysis_data/SBJ01155/umccrise/202408300c218043"
#' p2 <- "L2101566__L2101565"
#' p <- file.path(p1, p2)
#' obj <- UmccriseCanRepTables$new(p)
#' obj$path
#' obj$contents
#' d <- obj$read()
#' obj$write(d, out_dir = tempdir(), prefix = "sampleA", out_format = "tsv")
#' }
#'
#' @export
UmccriseCanRepTables <- R6::R6Class(
  "UmccriseCanRepTables",
  public = list(
    #' @field path Path to the `cancer_report_tables` directory.
    #' @field contents Tibble with file path, basename, and size.
    path = NULL,
    contents = NULL,
    #' @description Create a new UmccriseCanRepTables object.
    #' @param path Path to the `cancer_report_tables` directory.
    initialize = function(path = NULL) {
      stopifnot(is.character(path), length(path) == 1)
      self$path <- normalizePath(path)
      self$contents <- fs::dir_info(path, type = "file", recurse = TRUE) |>
        dplyr::mutate(
          bname = basename(.data$path),
          size = as.character(trimws(.data$size))
        ) |>
        dplyr::select("path", "bname", "size")
    },
    #' @description Print details about the cancer_report_tables directory.
    #' @param ... (ignored).
    print = function(...) {
      bnames <- self$contents |>
        dplyr::mutate(
          low = tolower(.data$bname),
        ) |>
        dplyr::arrange(.data$low) |>
        dplyr::mutate(
          n = dplyr::row_number(),
          bn = glue("{.data$n}. {.data$bname} ({.data$size})")
        ) |>
        dplyr::pull("bn")
      cat("#--- UmccriseCanRepTables ---#\n")
      cat(glue("Path: {self$path}"), "\n")
      cat("Contents:\n")
      cat(bnames, sep = "\n")
      invisible(self)
    },

    #' @description Read `chord.tsv.gz` file output from umccrise.
    #'
    #' @param x (`character(1)`)\cr
    #'   Path to `chord.tsv.gz` file.
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
    #' @description Read `hrdetect.tsv.gz` file output from umccrise.
    #'
    #' @param x (`character(1)`)\cr
    #'   Path to `hrdetect.tsv.gz` file.
    read_hrdetecttsv = function(x) {
      ct <- readr::cols(
        .default = "d",
        sample = "c"
      )
      read_tsvgz(x, col_types = ct) |>
        dplyr::select(-c("sample"))
    },


    #' @description Read `snv_20XX.tsv.gz` file output from umccrise.
    #'
    #' @param x (`character(1)`)\cr
    #'   Path to `snv_20XX.tsv.gz` file.
    read_sigs = function(x) {
      ct <- readr::cols(
        .default = "d",
        Signature = "c"
      )
      read_tsvgz(x, col_types = ct)
    },

    #' @description Read `qc_summary.tsv.gz` file output from umccrise.
    #'
    #' @param x (`character(1)`)\cr
    #'   Path to `qc_summary.tsv.gz` file.
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
    #' @description
    #' Reads contents of `cancer_report_tables` directory output by umccrise.
    #'
    #' @return A list of tibbles.
    #' @export
    read = function() {
      x <- self$path
      # now return all as list elements
      list(
        chord = grep_file(x, "-chord\\.tsv\\.gz$") |> self$read_chordtsv(),
        hrdetect = grep_file(x, "-hrdetect\\.tsv\\.gz$") |> self$read_hrdetecttsv(),
        sigs2015 = grep_file(x, "-snv_2015\\.tsv\\.gz$") |> self$read_sigs(),
        sigs2020 = grep_file(x, "-snv_2020\\.tsv\\.gz$") |> self$read_sigs(),
        sigsdbs = grep_file(x, "-dbs\\.tsv\\.gz$") |> self$read_sigs(),
        sigsindel = grep_file(x, "-indel\\.tsv\\.gz$") |> self$read_sigs(),
        qcsum = grep_file(x, "-qc_summary\\.tsv\\.gz$") |> self$read_qcsummarytsv()
      )
    }
  )
)
