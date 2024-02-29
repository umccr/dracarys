#' BcftoolsStatsFile R6 Class
#'
#' @description
#' Contains methods for reading and displaying contents of the
#' `bcftools_stats.txt` file (QUAL section) output from running `bcftools stats`.
#'
#' @examples
#' \dontrun{
#' x <- "/path/to/bcftools_stats.txt"
#' d <- BcftoolsStatsFile$new(x)
#' d_parsed <- d$read()
#' d$write(d_parsed, out_dir = tempdir(), prefix = "sample123", out_format = "tsv")
#' }
#' @export
BcftoolsStatsFile <- R6::R6Class(
  "BcftoolsStatsFile",
  inherit = File,
  public = list(
    #' @description
    #' Reads the QUAL section from the `bcftools_stats.txt` file.
    #'
    #' @return A tibble.
    read = function() {
      x <- self$path
      ln <- readr::read_lines(x)
      d <- ln[grepl("QUAL", ln)]
      # line1 ignore, line2 is colnames, just clean those up
      cnames <- c("QUAL_dummy", "id", "qual", "snps", "transi", "transv", "indels")
      d[3:length(d)] |>
        I() |>
        readr::read_tsv(col_names = cnames, col_types = readr::cols(.default = "d", "QUAL_dummy" = "c")) |>
        dplyr::select(-"QUAL_dummy")
    },
    #' @description
    #' Writes a tidy version of the `bcftools_stats.txt` (only QUAL section)
    #' file output from TSO.
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
      write_dracarys(obj = d, prefix = prefix, out_format = out_format, drid = drid)
    }
  )
)
