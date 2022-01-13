#' Replay R6 Class
#'
#' @description
#' Contains methods for reading contents of
#' the `replay.json` file output from DRAGEN, which contains the DRAGEN command
#' line, parameters and version for the specific run.
#'
#' @examples
#' x <- system.file("extdata/SEQC-II-replay.json", package = "dracarys")
#' r <- Replay$new(x)
#' r$read() # or read(r)
#' @export
Replay <- R6::R6Class("Replay", inherit = File, public = list(
  #' @description Reads the `replay.json` file.
  #' @return A list with the following elements:
  #'   - `command_line`: character of DRAGEN command line used.
  #'   - `dragen_config`: tibble of parameters used for the DRAGEN run.
  #'   - `system`: tibble with dragen_version, nodename, and kernel_release.
  #'   - `label`: character of sample label (inferred from file name)
  #'   - `hash_table_build`: tibble with details about the DRAGEN hash table build.
  read = function() {
    x <- self$path
    res <- x |>
      jsonlite::read_json(simplifyVector = TRUE) |>
      purrr::map_if(is.data.frame, tibble::as_tibble)

    req_elements <- c("command_line", "hash_table_build", "dragen_config", "system")
    assertthat::assert_that(all(names(res) %in% req_elements))

    res[["system"]] <- res[["system"]] |>
      tibble::as_tibble() |>
      tidyr::pivot_longer(dplyr::everything())
    res[["hash_table_build"]] <- res[["hash_table_build"]] |>
      tibble::as_tibble() |>
      tidyr::pivot_longer(dplyr::everything())
    res[["label"]] <- sub("-replay.json.*", "", basename(x))

    res
  }
))
