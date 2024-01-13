#' Tidy UMCCR Results
#'
#' Tidies UMCCR workflow results into a list of tibbles and writes individual
#' tibbles to TSV, Parquet, SparkDF, or RDS format.
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
#' `name` column in the `DR_FILE_REGEX` tibble).
#'
#' @return Tibble with path to input file and the resultant tidy object.
#' @examples
#' \dontrun{
#' in_dir <- paste0(
#'   "gds://production/analysis_data/SBJ01639/tso_ctdna_tumor_only/",
#'   "202204045ad5743c/L2200214/Results/PRJ220425_L2200214"
#' )
#' o1 <- sub("^gds://", "", in_dir)
#' out_dir <- file.path(fs::path_home(), "icav1/g", o1)
#' # in_dir <- file.path(out_dir, "dracarys_gds_sync")
#' prefix <- "SBJ01639"
#' out_format <- "rds"
#' umccr_tidy(in_dir = in_dir, out_dir = out_dir, prefix = prefix, out_format = out_format, dryrun = T)
#'
#' in_dir <- here::here(glue("nogit/tso/2022-12-13/SBJ02858/dracarys_gds_sync"))
#' out_dir <- file.path(in_dir, "../out")
#' gds_local_dir <- NULL
#' prefix <- "SBJ02858"
#' dryrun <- F
#' umccr_tidy(in_dir = in_dir, out_dir = out_dir, prefix = prefix, dryrun = F)
#' }
#' @export
umccr_tidy <- function(in_dir = NULL, out_dir = NULL, prefix = NULL,
                       gds_local_dir = NULL, out_format = "tsv",
                       dryrun = FALSE, token = Sys.getenv("ICA_ACCESS_TOKEN"),
                       pattern = NULL) {
  assertthat::assert_that(!is.null(in_dir), !is.null(out_dir), !is.null(prefix))
  dr_output_format_valid(out_format)
  e <- emojifont::emoji

  if (grepl("^gds://", in_dir)) {
    # in_dir is gds
    gds_local_dir <- gds_local_dir %||% file.path(out_dir, "dracarys_gds_sync")
    pat <- pattern %||% ".*" # keep all recognisable files
    dr_gds_download(
      gdsdir = in_dir, outdir = gds_local_dir, token = token,
      pattern = pat, dryrun = dryrun
    )
    # Now use the downloaded results
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
      regex <- DR_FILE_REGEX |>
        dplyr::pull("regex") |>
        sort()
      msg <- paste(
        "No UMCCR files for dracarys were found in {.file {in_dir}}.",
        "See current supported regexes: {regex}."
      )
      cli::cli_abort(msg)
    }
    dups <- dup_ftypes(d)
    if (dups$has_dup_ftypes) {
      cli::cli_alert_danger("Aborting - the input dir {.file {in_dir}} contains duplicated file types. See JSON below:")
      print(dups$msg_json)
      cli::cli_abort("Aborting - the input dir {.file {in_dir}} contains duplicated file types. See JSON above!")
    }
  }
  cli::cli_alert_info("{date_log()} {e('dragon')} {.emph {prefix}}: Start tidying UMCCR dir: {.file {in_dir}}")
  res <- d |>
    dplyr::select("type", "path", "bname") |>
    dplyr::filter(!.data$type %in% FILES_DOWNLOAD_BUT_IGNORE) |>
    dplyr::rowwise() |>
    dplyr::mutate(
      env = list(dr_func_eval(.data$type)),
      obj = list(.data$env$new(.data$path)),
      has_plot = "plot" %in% names(.data$env[["public_methods"]]),
      obj_parsed = ifelse(
        .data$type == "MultiqcFile",
        list(.data$obj$read(plot = TRUE, plot_names = "everything")),
        ifelse(
          .data$type %in% c(
            "TsoMergedSmallVariantsVcfFile",
            "TsoMergedSmallVariantsGenomeVcfFile",
            "TsoCopyNumberVariantsVcfFile"
          ),
          list(.data$obj$read(only_pass = FALSE)),
          list(.data$obj$read())
        )
      ),
      obj_parsed2 = list(.data$obj$write(.data$obj_parsed, out_dir = out_dir, prefix = glue("{prefix}_{.data$type}"), out_format = out_format)),
      plot = ifelse(.data$has_plot, list(.data$obj$plot(.data$obj_parsed)), list(NULL))
    ) |>
    dplyr::select("type", "path", "obj", dat = "obj_parsed", "plot")
  cli::cli_alert_success("{date_log()} {e('tada')} {.emph {prefix}}: UMCCR tidy results at: {.file {out_dir}}")
  return(invisible(res))
}

dup_ftypes <- function(d) {
  assertthat::assert_that(all(c("path", "type") %in% colnames(d)))
  dups <- d |>
    dplyr::group_by(.data$type) |>
    dplyr::mutate(type_count = dplyr::n()) |>
    dplyr::ungroup() |>
    dplyr::filter(.data$type_count > 1) |>
    dplyr::arrange(dplyr::desc(.data$type_count), .data$type, .data$path)
  if (nrow(dups) == 0) {
    return(list(has_dup_ftypes = FALSE, msg_json = NULL))
  }
  # show 2 duplicated files per type
  msg_json <- dups |>
    dplyr::group_by(.data$type) |>
    dplyr::slice_head(n = 2) |>
    dplyr::ungroup() |>
    dplyr::select("path", "type", "type_count") |>
    tidyr::nest(two_example_paths = .data$path) |>
    jsonlite::toJSON(pretty = F)
  return(list(has_dup_ftypes = TRUE, msg_json = msg_json))
}
