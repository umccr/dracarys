#' PCGR JSON Read
#'
#' @param x Path to file.
#'
#' @return A tibble with: `fracIndels`, `predicted_class`, `tmb_estimate`, `n_tmb`.
#'
#' @examples
#' \dontrun{
#' pcgr_json_read(x)
#' }
#' @export
pcgr_json_read <- function(x) {
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
    dplyr::select(indel_fraction = "fracIndels", "predicted_class")
  metrics <- dplyr::bind_cols(msi, tmb)
  return(metrics)
}
