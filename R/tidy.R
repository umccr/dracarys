#' Tidy Files
#'
#' Tidies files into a tibble with parsed data.
#'
#' @param x Tibble with `localpath` to file and the function `type` to parse it.
#' The function must return a tibble with a `name` column and the tidied `data`
#' as a list-column (see example).
#' @param envir the environment in which to evaluate the function e.g. use `self`
#' when using inside R6 classes.
#'
#' @return Tibble with parsed data in a `data` list-column.
#' @examples
#' \dontrun{
#' p1 <- "~/icav1/g/production/analysis_data/SBJ01155/umccrise/202408300c218043"
#' p2 <- "L2101566__L2101565/SBJ01155__PRJ211091/cancer_report_tables"
#' p <- file.path(p1, p2, "SBJ01155__PRJ211091-qc_summary.tsv.gz")
#' p_dl <- file.path(
#'   p1, "L2101566__L2101565/SBJ01155__PRJ211091/small_variants",
#'   "SBJ01155__PRJ211091-somatic-PASS.vcf.gz"
#' )
#' fun <- function(x) {
#'   d <- readr::read_tsv(x)
#'   tibble::tibble(name = "table1", data = list(d[]))
#' }
#' x <- tibble::tibble(
#'   type = c("fun", "DOWNLOAD_ONLY_foobar"), localpath = c(p, p_dl)
#' )
#' tidy_files(x)
#' }
#'
#' @export
tidy_files <- function(x, envir = parent.frame()) {
  assertthat::assert_that(is.data.frame(x))
  assertthat::assert_that(all(c("type", "localpath") %in% colnames(x)))
  # if there's a DOWNLOAD_ONLY_suffix, extract that suffix and call
  # the DOWNLOAD_ONLY function
  extract_download_suffix <- function(s) {
    sub("DOWNLOAD_ONLY(.*)", "\\1", s)
  }
  x |>
    dplyr::rowwise() |>
    dplyr::mutate(
      data = ifelse(
        !grepl("DOWNLOAD_ONLY", .data$type),
        list(
          dr_func_eval(
            f = .data$type, v = .data$type, envir = envir
          )(.data$localpath)
        ),
        list(
          dr_func_eval(
            f = "DOWNLOAD_ONLY", v = "DOWNLOAD_ONLY", envir = envir
          )(.data$localpath, extract_download_suffix(.data$type))
        )
      )
    ) |>
    dplyr::ungroup() |>
    dplyr::select("data") |>
    tidyr::unnest("data")
}

#' Tidy UMCCR Results
#'
#' Tidies UMCCR workflow results into a list of tibbles and writes individual
#' tibbles to TSV, Parquet, SparkDF, or RDS format.
#'
#' @param in_dir Directory path to UMCCR workflow results (can be S3 or local).
#' @param prefix Prefix of output file(s).
#' @param out_dir Output directory.
#' @param local_dir If `indir` is a S3 directory, 'recognisable' files
#' will be first downloaded to this directory (def: <out_dir>/dracarys_s3_sync).
#' @param dryrun Just list the files that will be downloaded (def: FALSE).
#' @param out_format Format of output (tsv, parquet, both) (def: tsv).
#' @param pattern Pattern to further filter the returned file type tibble (see
#' `name` column in the `DR_FILE_REGEX` tibble).
#' @param regexes Tibble with `regex` and `fun`ction name.
#'
#' @return Tibble with path to input file and the resultant tidy object.
#' @examples
#' \dontrun{
#' in_dir <- file.path(
#'   "s3://umccr-primary-data-prod/UMCCR-Validation/SBJ00596",
#'   "ctTSO/2021-03-17/PTC_SSqCMM05pc_L2100067"
#' )
#' o1 <- sub("s3:/", "~/s3", in_dir)
#' out_dir <- o1
#' out_dir <- file.path(fs::path_home(), "icav1/g", o1)
#' prefix <- "SBJ01639"
#' prefix <- "PTC_SSqCMM05pc_L2100067"
#' out_format <- "rds"
#' umccr_tidy(in_dir = in_dir, out_dir = out_dir, prefix = prefix, out_format = out_format, dryrun = F)
#'
#' in_dir <- here::here(glue("nogit/tso/2022-12-13/SBJ02858/dracarys_gds_sync"))
#' out_dir <- file.path(in_dir, "../out")
#' prefix <- "SBJ02858"
#' dryrun <- F
#' umccr_tidy(in_dir = in_dir, out_dir = out_dir, prefix = prefix, dryrun = F)
#' }
#' @export
umccr_tidy <- function(in_dir = NULL, out_dir = NULL, prefix = NULL,
                       local_dir = NULL, out_format = "tsv",
                       dryrun = FALSE, pattern = NULL, regexes = DR_FILE_REGEX) {
  assertthat::assert_that(!is.null(in_dir), !is.null(out_dir), !is.null(prefix))
  dr_output_format_valid(out_format)
  e <- emojifont::emoji

  if (grepl("^s3://", in_dir)) {
    # in_dir is s3
    cloud_type <- "s3"
    local_dir <- local_dir %||% file.path(out_dir, glue("dracarys_{cloud_type}_sync"))
    pat <- pattern %||% ".*" # keep all recognisable files
    dr_s3_download(
      s3dir = in_dir, outdir = local_dir,
      pattern = pat, dryrun = dryrun
    )
    # Now use the downloaded results
    in_dir <- local_dir
  } else {
    # in_dir is not s3
    if (!is.null(local_dir)) {
      cli::cli_warn(glue(
        "You have specified the 'local_dir' option to download S3 results, ",
        "but your input directory is local. Ignoring that option."
      ))
    }
  }
  if (dryrun) {
    cli::cli_inform("{e('camel')} You have specified 'dryrun' - just listing files!")
    return(NULL)
  } else {
    d <- fs::dir_ls(in_dir, recurse = TRUE) |>
      tibble::as_tibble_col(column_name = "path") |>
      dplyr::mutate(
        bname = basename(.data$path),
        type = purrr::map_chr(.data$bname, \(x) match_regex(x, regexes = regexes))
      ) |>
      dplyr::filter(!is.na(.data$type))

    if (nrow(d) == 0) {
      regex <- regexes |>
        dplyr::pull("regex") |>
        sort()
      msg <- paste(
        "No UMCCR files for dracarys were found in {.file {in_dir}}.",
        "See supported regexes: {regex}."
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
