#' PcgrJson R6 Class
#'
#' @description
#' Contains methods for reading and displaying contents of the
#' `pcgr.json.gz` file output from PCGR.
#'
#' @examples
#' \dontrun{
#' x <- "/path/to/pcgr.json.gz"
#' d <- PcgrJsonFile$new(x)
#' d_parsed <- d$read() # or read(d)
#' d$write(d_parsed, out_dir = tempdir(), prefix = "sample705", out_format = "both")
#' }
#' @export
PcgrJsonFile <- R6::R6Class(
  "PcgrJsonFile",
  inherit = File,
  public = list(
    #' @description
    #' Reads the `pcgr.json.gz` file output from PCGR.
    #'
    #' @return List of tibbles.
    read = function() {
      x <- self$path
      j <- jsonlite::read_json(x)
      l2tib <- function(el) {
        purrr::flatten(el) |>
          dplyr::bind_rows() |>
          dplyr::mutate(dplyr::across(dplyr::everything(), ~ as.character(.)))
      }
      dbrel <- j[["metadata"]][["pcgr_db_release"]] |>
        purrr::map(l2tib) |>
        dplyr::bind_rows(.id = "name_tidy") |>
        dplyr::select("name", "name_tidy", "version", "url", "resource_type")
      tmb <- j[["content"]][["tmb"]][["variant_statistic"]] |>
        purrr::flatten() |>
        tibble::as_tibble() |>
        dplyr::select("tmb_estimate", "n_tmb", "tmb_tertile", "target_size_mb")
      msi <- j[["content"]][["msi"]][["prediction"]][["msi_stats"]] |>
        purrr::flatten() |>
        tibble::as_tibble() |>
        dplyr::select("sample_name", "fracIndels", "predicted_class")
      # optimisation required
      snv <- j[["content"]][["snv_indel"]][["variant_set"]][["tsv"]] |>
        purrr::map(l2tib) |>
        dplyr::bind_rows() |>
        readr::type_convert()
      metrics <- dplyr::bind_cols(msi, tmb)

      res <- list(
        metrics = metrics,
        snv = snv,
        dbrel = dbrel
      )
    },

    #' @description
    #' Writes a tidy version of the `pcgr.json.gz` file output from PCGR.
    #'
    #' @param d Parsed object from `self$read()`.
    #' @param prefix Prefix of output file(s).
    #' @param out_dir Output directory.
    #' @param out_format Format of output file(s) (one of 'tsv' (def.),
    #' 'parquet', 'both').
    write = function(d, out_dir, prefix, out_format = "tsv") {
      prefix <- file.path(out_dir, prefix)
      p <- glue("{prefix}_pcgr")
      l <- list(
        meta = list(
          obj = d[["metrics"]],
          pref = glue("{p}_metrics")
        ),
        snv = list(
          obj = d[["snv"]],
          pref = glue("{p}_snv")
        ),
        dbrel = list(
          obj = d[["dbrel"]],
          pref = glue("{p}_dbrel")
        )
      )
      purrr::map(l, function(k) {
        write_dracarys(obj = k[["obj"]], prefix = k[["pref"]], out_format = out_format)
      })
    }
  )
)
