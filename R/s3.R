#' List Objects in AWS S3 Directory
#'
#' Returns some or all (up to 1,000) of the objects in an S3 directory.
#'
#' @param s3dir S3 directory.
#' @param max_objects Maximum objects returned.
#'
#'
#' @return A tibble with object basename, size, last modified timestamp, and
#' full S3 path.
#' @examples
#' \dontrun{
#' p1 <- "s3://org.umccr.data.oncoanalyser/analysis_data/SBJ05373/sash"
#' p2 <- "20240707becde493/L2401018_L2401017/SBJ05373_MDX240220"
#' s3dir <- file.path(p1, p2, "cancer_report/cancer_report_tables")
#' s3_list_objects_dir(s3dir, max_objects = 15)
#' }
#' @export
s3_list_objects_dir <- function(s3dir, max_objects = 1000) {
  assertthat::assert_that(grepl("^s3://", s3dir))
  bucket <- sub("s3://(.*?)/.*", "\\1", s3dir)
  prefix <- sub("s3://(.*?)/(.*)", "\\2", s3dir)
  s3 <- paws.storage::s3()
  l <- s3$list_objects_v2(Bucket = bucket, Prefix = prefix, MaxKeys = max_objects)
  assertthat::assert_that(all(c("Contents", "KeyCount") %in% names(l)))
  cols_sel <- c("bname", "size", "lastmodified", "path")
  # handle no results
  if (l[["KeyCount"]] == 0) {
    return(empty_tbl(cnames = cols_sel, ctypes = "cccc"))
  }
  d <- l[["Contents"]] |>
    purrr::map(\(x) tibble::tibble(
      Key = x[["Key"]],
      Size = x[["Size"]],
      lastmodified = x[["LastModified"]]
    )) |>
    dplyr::bind_rows() |>
    dplyr::mutate(
      path = glue("s3://{bucket}/{.data$Key}"),
      bname = basename(.data$path),
      size = fs::as_fs_bytes(.data$Size)
    ) |>
    dplyr::select(dplyr::all_of(cols_sel))
  return(d)
}

#' List Relevant Files In AWS S3 Directory
#'
#' Lists relevant files in an AWS S3 directory.
#'
#' @param s3dir S3 directory.
#' @param pattern Pattern to further filter the returned file type tibble.
#' @param max_objects The total number of objects to return.
#' @param presign Include presigned URLs (def: FALSE).
#' @param expiry_sec Number of seconds the presigned URL will be valid for (if generated).
#' @param regexes Tibble with `regex` and `fun`ction name.
#'
#' @return A tibble with file type, basename, file size, date, full path, and presigned URL if requested.
#' @examples
#' \dontrun{
#' s3dir <- "s3://umccr-primary-data-prod/cancer_report_tables"
#' s3_files_list_filter_relevant(s3dir = s3dir, presign = FALSE)
#' }
#' @export
s3_files_list_filter_relevant <- function(s3dir, pattern = NULL, max_objects = 100,
                                          presign = FALSE, expiry_sec = 3600,
                                          regexes = DR_FILE_REGEX) {
  assertthat::assert_that(rlang::is_logical(presign), max_objects <= 1000)
  d_all <- s3_list_objects_dir(s3dir = s3dir, max_objects = max_objects)
  if (nrow(d_all) == 0) {
    return(d_all)
  }
  pattern <- pattern %||% ".*" # keep all recognisable files by default
  d <- d_all |>
    dplyr::rowwise() |>
    dplyr::mutate(
      type = purrr::map_chr(.data$bname, \(x) match_regex(x, regexes))
    ) |>
    dplyr::ungroup() |>
    dplyr::filter(!is.na(.data$type), grepl(pattern, .data$type)) |>
    dplyr::select("type", "bname", "size", "lastmodified", "path")

  if (presign) {
    if (nrow(d) == 0) {
      return(d)
    }
    s3_client <- paws.storage::s3(paws.storage::config(signature_version = "s3v4"))
    d <- d |>
      dplyr::rowwise() |>
      dplyr::mutate(presigned_url = s3_file_presignedurl(
        client = s3_client, s3path = .data$path, expiry_seconds = expiry_sec
      )) |>
      dplyr::ungroup()
  }
  d
}

#' S3 Generate Presigned URL
#'
#' @param client S3 client. Make sure you use `signature_version = "s3v4"` (see example).
#' @param s3path Full path to S3 object.
#' @param expiry_seconds Number of seconds the presigned URL is valid for (3600 = 1 hour).
#'
#' @return An S3 presigned URL.
#' @examples
#' \dontrun{
#' client <- paws.storage::s3(paws.storage::config(signature_version = "s3v4"))
#' s3path <- "s3://bucket1/path/to/file.tsv"
#' s3_file_presignedurl(client, s3path)
#' }
#'
#' @export
s3_file_presignedurl <- function(client, s3path, expiry_seconds = 3600) {
  bucket <- sub("s3://(.*?)/.*", "\\1", s3path)
  prefix <- sub("s3://(.*?)/(.*)", "\\2", s3path)
  client$generate_presigned_url(
    client_method = "get_object",
    params = list(Bucket = bucket, Key = prefix),
    expires_in = expiry_seconds
  )
}

#' Search AWS S3 Objects
#'
#' Searches for the given pattern in the UMCCR `umccr-primary-data-prod` AWS S3
#' bucket.
#'
#' @param pat Pattern to search for (e.g. 'multiqc_data.json').
#' @param rows Max number of rows to return.
#'
#' @return Tibble with S3 path, object size, date modified, id, unique hash.
#'
#' @examples
#' \dontrun{
#' pat <- "qc_summary.tsv.gz"
#' s3_search(pat, 10)
#' }
#' @export
s3_search <- function(pat, rows) {
  au_tz <- "Australia/Melbourne"
  utc_tz <- "UTC"
  base_url <- "https://api.portal.prod.umccr.org/iam/s3"
  url1 <- utils::URLencode(glue("{base_url}?rowsPerPage={rows}&search={pat}"))
  awscurl_cmd <- glue(
    "awscurl '{url1}' ",
    "--header 'Accept: application/json'"
  )
  message(glue("Running {awscurl_cmd}"))
  j <- system(awscurl_cmd, intern = TRUE)
  date_fmt <- "%Y-%m-%dT%H:%M:%S"
  d <- j |>
    jsonlite::fromJSON() |>
    purrr::pluck("results") |>
    tibble::as_tibble()
  d |>
    dplyr::mutate(
      date1 = as.POSIXct(.data$last_modified_date, tz = utc_tz, format = date_fmt),
      date_aest = lubridate::with_tz(.data$date1, tz = au_tz),
      path = glue("s3://{bucket}/{key}"),
      size = fs::as_fs_bytes(.data$size)
    ) |>
    dplyr::select("path", "size", "date_aest", "id", "unique_hash")
}

#' dracarys S3 Download
#'
#' Download only S3 files that can be processed by dracarys.
#'
#' @param s3dir Full path to S3 directory.
#' @param outdir Path to output directory.
#' @param max_objects Maximum objects returned in file listing.
#' @param pattern Pattern to further filter the returned file type tibble.
#' @param regexes Tibble with regex and function name.
#' @param dryrun If TRUE, just list the files that will be downloaded (don't
#' download them).
#' @examples
#' \dontrun{
#' s3dir <- file.path(
#'   "s3://umccr-primary-data-prod/UMCCR-Validation/SBJ00596",
#'   "ctTSO/2021-03-17/PTC_SSqCMM05pc_L2100067"
#' )
#' outdir <- sub("s3:/", "~/s3", s3dir)
#' dr_s3_download(s3dir = s3dir, outdir = outdir, max_objects = 1000, dryrun = F)
#' }
#' @export
dr_s3_download <- function(s3dir, outdir, max_objects = 100, pattern = NULL,
                           regexes = DR_FILE_REGEX, dryrun = FALSE) {
  s3 <- paws.storage::s3()
  e <- emojifont::emoji
  fs::dir_create(outdir)
  d <- s3_files_list_filter_relevant(
    s3dir = s3dir, pattern = NULL, max_objects = max_objects, presign = FALSE, regexes = regexes
  )
  d <- d |>
    dplyr::select("type", "size", "path", "bname") |>
    dplyr::mutate(out = file.path(outdir, .data$bname))

  # download recognisable dracarys files to outdir/{bname}
  if (!dryrun) {
    cli::cli_alert_info("{date_log()} {e('arrow_heading_down')} Downloading files from {.file {s3dir}}")
    d |>
      dplyr::rowwise() |>
      dplyr::mutate(
        s3bucket = sub("s3://(.*?)/.*", "\\1", .data$path),
        s3key = sub("s3://(.*?)/(.*)", "\\2", .data$path),
        dl = list(
          s3$download_file(
            Bucket = .data$s3bucket, Key = .data$s3key, Filename = .data$out
          )
        )
      )
  } else {
    cli::cli_alert_info("{date_log()} {e('camera')} Just list relevant files from {.file {s3dir}}")
    d |>
      dplyr::select("path", "type", "size") |>
      as.data.frame() |>
      print()
  }
}
