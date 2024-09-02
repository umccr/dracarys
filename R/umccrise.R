#' UmccriseCanRepTables R6 Class
#'
#' @description
#' Reads and writes tidy versions of files within the `cancer_report_tables` directory
#' output from the `umccrise` workflow.
#'
#' @examples
#' \dontrun{
#' p1 <- "~/icav1/g/production/analysis_data/SBJ01155/umccrise/202408300c218043"
#' p2 <- "L2101566__L2101565/SBJ01155__PRJ211091/cancer_report_tables"
#' p <- file.path(p1, p2)
#' obj <- UmccriseCanRepTables$new(p)
#' obj$path
#' obj$contents
#' d <- obj$read()
#' b$write(d, out_dir = tempdir(), prefix = "sampleA", out_format = "tsv")
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
      self$contents <- fs::dir_info(path) |>
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
    }
  )
)
