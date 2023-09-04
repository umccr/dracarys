#' UmChordTsvFile R6 Class
#'
#' @description
#' Contains methods for reading and displaying contents of the
#' `chord.tsv.gz` file output from umccrise.
#'
#' @examples
#' \dontrun{
#' x <- "/path/to/chord.tsv.gz"
#' d <- UmChordTsvFile$new(x)
#' d_parsed <- d$read() # or read(d)
#' d$write(d_parsed, out_dir = tempdir(), prefix = "sample705", out_format = "tsv")
#' }
#' @export
UmChordTsvFile <- R6::R6Class(
  "UmChordTsvFile",
  inherit = File,
  public = list(
    #' @description
    #' Reads the `chord.tsv.gz` file output from umccrise.
    #'
    #' @return A tibble.
    read = function() {
      x <- self$path
      ct <- readr::cols_only(
        p_hrd = "d",
        hr_status = "c",
        hrd_type = "c",
        p_BRCA1 = "d",
        p_BRCA2 = "d"
      )
      read_tsvgz(x, col_types = ct)
    },

    #' @description
    #' Writes a tidy version of the `chord.tsv.gz` file output from umccrise.
    #'
    #' @param d Parsed object from `self$read()`.
    #' @param prefix Prefix of output file(s).
    #' @param out_dir Output directory.
    #' @param out_format Format of output file(s).
    #' @param drid dracarys ID to use for the dataset (e.g. `wfrid.123`, `prid.456`).
    write = function(d, out_dir = NULL, prefix, out_format = "tsv", drid = NULL) {
      if (!is.null(out_dir)) {
        prefix <- file.path(out_dir, prefix)
      }
      # prefix2 <- glue("{prefix}_chord")
      write_dracarys(obj = d, prefix = prefix, out_format = out_format, drid = drid)
    }
  )
)

#' UmHrdetectTsvFile R6 Class
#'
#' @description
#' Contains methods for reading and displaying contents of the
#' `hrdetect.tsv.gz` file output from umccrise.
#'
#' @examples
#' \dontrun{
#' x <- "/path/to/hrdetect.tsv.gz"
#' d <- UmHrdetectTsvFile$new(x)
#' d_parsed <- d$read() # or read(d)
#' d$write(d_parsed, out_dir = tempdir(), prefix = "sample705", out_format = "tsv")
#' }
#' @export
UmHrdetectTsvFile <- R6::R6Class(
  "UmHrdetectTsvFile",
  inherit = File,
  public = list(
    #' @description
    #' Reads the `hrdetect.tsv.gz` file output from umccrise.
    #'
    #' @return A tibble.
    read = function() {
      x <- self$path
      ct <- readr::cols(
        .default = "d",
        sample = "c"
      )
      read_tsvgz(x, col_types = ct) |>
        dplyr::select(-c("sample"))
    },

    #' @description
    #' Writes a tidy version of the `hrdetect.tsv.gz` file output from umccrise.
    #'
    #' @param d Parsed object from `self$read()`.
    #' @param prefix Prefix of output file(s).
    #' @param out_dir Output directory.
    #' @param out_format Format of output file(s).
    #' @param drid dracarys ID to use for the dataset (e.g. `wfrid.123`, `prid.456`).
    write = function(d, out_dir, prefix, out_format = "tsv", drid = NULL) {
      if (!is.null(out_dir)) {
        prefix <- file.path(out_dir, prefix)
      }
      # prefix2 <- glue("{prefix}_hrdetect")
      write_dracarys(obj = d, prefix = prefix, out_format = out_format, drid = drid)
    }
  )
)

#' UmSigsSnvFile R6 Class
#'
#' @description
#' Contains methods for reading and displaying contents of the
#' `snv_20XX.tsv.gz` file with SNV signatures output from umccrise.
#'
#' @examples
#' \dontrun{
#' x <- "/path/to/snv_2015.tsv.gz"
#' d <- UmSigsSnvFile$new(x)
#' d_parsed <- d$read() # or read(d)
#' d$write(d_parsed, out_dir = tempdir(), prefix = "sample705", out_format = "tsv")
#' }
#' @export
UmSigsSnvFile <- R6::R6Class(
  "UmSigsSnvFile",
  inherit = File,
  public = list(
    #' @description
    #' Reads the `snv.tsv.gz` file output from umccrise.
    #'
    #' @return A tibble.
    read = function() {
      x <- self$path
      version <- dplyr::if_else(grepl("2015\\.tsv.\\gz", x), "2015", "2020")
      ct <- readr::cols(
        .default = "d",
        Signature = "c"
      )
      read_tsvgz(x, col_types = ct)
    },

    #' @description
    #' Writes a tidy version of the `snv_20XX.tsv.gz` signature file output from umccrise.
    #'
    #' @param d Parsed object from `self$read()`.
    #' @param prefix Prefix of output file(s).
    #' @param out_dir Output directory.
    #' @param out_format Format of output file(s).
    #' @param drid dracarys ID to use for the dataset (e.g. `wfrid.123`, `prid.456`).
    write = function(d, out_dir, prefix, out_format = "tsv") {
      if (!is.null(out_dir)) {
        prefix <- file.path(out_dir, prefix)
      }
      # prefix2 <- glue("{prefix}_sigs_snv")
      write_dracarys(obj = d, prefix = prefix, out_format = out_format, drid = drid)
    }
  )
)

#' UmQcSumFile R6 Class
#'
#' @description
#' Contains methods for reading and displaying contents of the
#' `qc_summary.tsv.gz` file with QC summary metrics output from umccrise.
#'
#' @examples
#' \dontrun{
#' x <- "/path/to/snv_2015.tsv.gz"
#' d <- UmQcSumFile$new(x)
#' d_parsed <- d$read() # or read(d)
#' d$write(d_parsed, out_dir = tempdir(), prefix = "sample705", out_format = "tsv")
#' }
#' @export
UmQcSumFile <- R6::R6Class(
  "UmQcSumFile",
  inherit = File,
  public = list(
    #' @description
    #' Reads the `qc_summary.tsv.gz` file output from umccrise.
    #'
    #' @return A tibble.
    read = function() {
      x <- self$path
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
          hypermutated, bpi_enabled
        )
    },

    #' @description
    #' Writes a tidy version of the `qc_summary.tsv.gz` QC summary file output
    #' from umccrise.
    #'
    #' @param d Parsed object from `self$read()`.
    #' @param prefix Prefix of output file(s).
    #' @param out_dir Output directory.
    #' @param out_format Format of output file(s).
    #' @param drid dracarys ID to use for the dataset (e.g. `wfrid.123`, `prid.456`).
    write = function(d, out_dir, prefix, out_format = "tsv", drid = NULL) {
      if (!is.null(out_dir)) {
        prefix <- file.path(out_dir, prefix)
      }
      # prefix2 <- glue("{prefix}_qc_summary")
      write_dracarys(obj = d, prefix = prefix, out_format = out_format, drid = drid)
    }
  )
)
