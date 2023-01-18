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
#' d$write(d_parsed, out_dir = tempdir(), prefix = "sample705", out_format = "both")
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
      readr::read_tsv(x, col_types = ct)
    },

    #' @description
    #' Writes a tidy version of the `chord.tsv.gz` file output from umccrise.
    #'
    #' @param d Parsed object from `self$read()`.
    #' @param prefix Prefix of output file(s).
    #' @param out_dir Output directory.
    #' @param out_format Format of output file(s) (one of 'tsv' (def.),
    #' 'parquet', 'both').
    write = function(d, out_dir, prefix, out_format = "tsv") {
      prefix <- file.path(out_dir, prefix)
      prefix2 <- glue("{prefix}_chord")
      write_dracarys(obj = d, prefix = prefix2, out_format = out_format)
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
#' d$write(d_parsed, out_dir = tempdir(), prefix = "sample705", out_format = "both")
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
      readr::read_tsv(x, col_types = ct) |>
        dplyr::select(-c("sample"))
    },

    #' @description
    #' Writes a tidy version of the `hrdetect.tsv.gz` file output from umccrise.
    #'
    #' @param d Parsed object from `self$read()`.
    #' @param prefix Prefix of output file(s).
    #' @param out_dir Output directory.
    #' @param out_format Format of output file(s) (one of 'tsv' (def.),
    #' 'parquet', 'both').
    write = function(d, out_dir, prefix, out_format = "tsv") {
      prefix <- file.path(out_dir, prefix)
      prefix2 <- glue("{prefix}_hrdetect")
      write_dracarys(obj = d, prefix = prefix2, out_format = out_format)
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
#' d$write(d_parsed, out_dir = tempdir(), prefix = "sample705", out_format = "both")
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
      list(
        data = readr::read_tsv(x, col_types = ct),
        version = version
      )
    },

    #' @description
    #' Writes a tidy version of the `snv_20XX.tsv.gz` signature file output from umccrise.
    #'
    #' @param d Parsed object from `self$read()`.
    #' @param prefix Prefix of output file(s).
    #' @param out_dir Output directory.
    #' @param out_format Format of output file(s) (one of 'tsv' (def.),
    #' 'parquet', 'both').
    write = function(d, out_dir, prefix, out_format = "tsv") {
      prefix <- file.path(out_dir, prefix)
      version <- d[["version"]]
      prefix2 <- glue("{prefix}_sigs_snv{version}")
      write_dracarys(obj = d[["data"]], prefix = prefix2, out_format = out_format)
    }
  )
)
