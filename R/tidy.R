#' Tidy UMCCR Results
#'
#' Tidies UMCCR workflow results into a list of tibbles and writes individual
#' tibbles to TSV and/or Parquet format.
#'
#' @param in_dir Directory path to UMCCR workflow results (can be GDS or local).
#' @param prefix Prefix of output file(s).
#' @param out_dir Output directory.
#' @param gds_local_dir If `indir` is a GDS directory, 'recognisable' files
#' will be first downloaded to this directory.
#' @param dryrun Just list the files that will be downloaded (def: FALSE).
#' @param token ICA access token (by default uses $ICA_ACCESS_TOKEN env var).
#' @param out_format Format of output (tsv, parquet, both) (def: tsv).
#' @param pattern Pattern to further filter the returned file type tibble (see
#' `name` column in the `FILE_REGEX` tibble).
#'
#' @return Tibble with path to input file and the resultant tidy object.
#' @examples
#' \dontrun{
#' in_dir <- paste0(
#'   "gds://production/analysis_data/SBJ02858/tso_ctdna_tumor_only/",
#'   "20221104b7ad0b38/L2201560/Results/PRJ222206_L2201560/"
#' )
#' in_dir <- here::here(glue("nogit/tso/2022-12-13/SBJ02858/dracarys_gds_sync"))
#' out_dir <- file.path(in_dir, "../out")
#' gds_local_dir <- NULL
#' prefix <- "SBJ02858"
#' dryrun <- F
#' umccr_tidy(in_dir = in_dir, out_dir = out_dir, prefix = prefix)
#' }
#' @export
umccr_tidy <- function(in_dir, out_dir, prefix, gds_local_dir = NULL, out_format = "tsv",
                       dryrun = FALSE, token = Sys.getenv("ICA_ACCESS_TOKEN"),
                       pattern = NULL) {
  output_format_valid(out_format)
  e <- emojifont::emoji

  if (grepl("^gds://", in_dir)) {
    # in_dir is gds
    gds_local_dir <- gds_local_dir %||% file.path(out_dir, "dracarys_gds_sync")
    pat <- pattern %||% ".*" # keep all recognisable files
    dr_gds_download(
      gdsdir = in_dir, outdir = gds_local_dir, token = token,
      pattern = pat, dryrun = dryrun
    )
    # Use the downloaded results
    in_dir <- gds_local_dir
  } else {
    # in_dir is not gds
    if (!is.null(gds_local_dir)) {
      cli::cli_warn(glue(
        "You have specified the 'gds_local_dir' option to download GDS results, ",
        "but your input directory is local. Ignoring that option."
      ))
    }
  }
  if (dryrun) {
    cli::cli_inform("{e('camel')} You have specified 'dryrun' - just listing files!")
    return(NULL)
  } else {
    # TODO: list recursively to match default ICA API, might change later
    d <- fs::dir_ls(in_dir, recurse = TRUE) |>
      tibble::as_tibble_col(column_name = "path") |>
      dplyr::mutate(
        bname = basename(.data$path),
        type = purrr::map_chr(.data$bname, match_regex)
      ) |>
      dplyr::filter(!is.na(.data$type))

    if (nrow(d) == 0) {
      regex <- FILE_REGEX |>
        dplyr::pull("regex")
      msg <- paste(
        "No UMCCR files for dracarys were found in {.file {in_dir}}.",
        "See current supported regexes: {regex}."
      )
      cli::cli_abort(msg)
    }

    cli::cli_alert_info("{date_log()} {e('dragon')} {.emph {prefix}}: Start tidying UMCCR dir: {.file {in_dir}}")
    res <- d |>
      dplyr::select("type", "path", "bname") |>
      dplyr::rowwise() |>
      dplyr::mutate(
        env = list(func_selector(.data$type)),
        obj = list(.data$env$new(.data$path)),
        has_plot = "plot" %in% names(.data$env[["public_methods"]]),
        obj_parsed = list(.data$obj$read()),
        obj_parsed2 = list(.data$obj$write(.data$obj_parsed, out_dir = out_dir, prefix = prefix, out_format = out_format)),
        plot = ifelse(.data$has_plot, list(.data$obj$plot(.data$obj_parsed)), list(NULL))
      ) |>
      dplyr::select("type", "path", "obj", dat = "obj_parsed", "plot")
    # TODO:
    assertthat::assert_that(
      !any(duplicated(d[["type"]])),
      msg = glue("Aborting - duplicated file types detected in {in_dir}")
    )
    cli::cli_alert_success("{date_log()} {e('tada')} {.emph {prefix}}: UMCCR tidy results at: {.file {out_dir}}")
    return(invisible(res))
  }
}
