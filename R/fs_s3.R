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
#' s3_list_files_dir(s3dir, max_objects = 15)
#' }
#' @export
s3_list_files_dir <- function(s3dir, max_objects = 1000) {
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
      lastmodified = as.character(x[["LastModified"]])
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
#' @inheritParams s3_list_files_dir
#' @param pattern Pattern to further filter the returned file type tibble.
#' @param regexes Tibble with `regex` and `fun`ction name.
#' @param presign Include presigned URLs (def: FALSE).
#' @param expiry_sec Number of seconds the presigned URL will be valid for (if generated).
#'
#' @return A tibble with file type, basename, size, last modified timestamp,
#' full path, and presigned URL if requested.
#' @examples
#' \dontrun{
#' p1 <- "s3://org.umccr.data.oncoanalyser/analysis_data/SBJ05373/sash"
#' p2 <- "20240707becde493/L2401018_L2401017/SBJ05373_MDX240220"
#' s3dir <- file.path(p1, p2)
#' regexes <- tibble::tibble(regex = "multiqc_data\\.json$", fun = "MultiqcJsonFile")
#' s3_list_files_filter_relevant(s3dir = s3dir, regexes = regexes, max_objects = 300)
#' }
#' @export
s3_list_files_filter_relevant <- function(s3dir, pattern = NULL,
                                          regexes = DR_FILE_REGEX, max_objects = 100,
                                          presign = FALSE, expiry_sec = 3600) {
  assertthat::assert_that(rlang::is_logical(presign), max_objects <= 1000)
  d_all <- s3_list_files_dir(s3dir = s3dir, max_objects = max_objects)
  if (nrow(d_all) == 0) {
    return(d_all)
  }
  pattern <- pattern %||% ".*" # keep all recognisable files by default
  cols_sel <- c("type", "bname", "size", "lastmodified", "path")
  d <- d_all |>
    dplyr::rowwise() |>
    dplyr::mutate(
      type = purrr::map_chr(.data$path, \(x) match_regex(x, regexes))
    ) |>
    dplyr::ungroup() |>
    dplyr::filter(!is.na(.data$type), grepl(pattern, .data$type)) |>
    dplyr::select(dplyr::all_of(cols_sel))

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
      dplyr::ungroup() |>
      dplyr::select(dplyr::all_of(c(cols_sel, "presigned_url")))
  }
  d
}

#' dracarys S3 Download
#'
#' Download only S3 files that can be processed by dracarys.
#'
#' @inheritParams s3_list_files_dir
#' @inheritParams s3_list_files_filter_relevant
#' @param outdir Path to output directory.
#' @param dryrun If TRUE, just list the files that will be downloaded (don't
#' download them).
#' @examples
#' \dontrun{
#' p1 <- "s3://org.umccr.data.oncoanalyser/analysis_data/SBJ05373/sash"
#' p2 <- "20240707becde493/L2401018_L2401017/SBJ05373_MDX240220"
#' s3dir <- file.path(p1, p2)
#' regexes <- tibble::tribble(
#'   ~regex, ~fun,
#'   "multiqc_data\\.json$", "MultiqcJsonFile",
#'   "pcgr.*\\.json\\.gz$", "pcgrjson"
#' )
#' outdir <- sub("s3:/", "~/s3", s3dir)
#' dr_s3_download(s3dir = s3dir, outdir = outdir, max_objects = 500, regexes = regexes, dryrun = F)
#' }
#' @export
dr_s3_download <- function(s3dir, outdir, max_objects = 100, pattern = NULL,
                           regexes = DR_FILE_REGEX, dryrun = FALSE) {
  s3 <- paws.storage::s3()
  e <- emojifont::emoji
  fs::dir_create(outdir)
  d <- s3_list_files_filter_relevant(
    s3dir = s3dir, pattern = NULL, regexes = regexes,
    max_objects = max_objects, presign = FALSE
  )
  msg <- glue(
    "S3 input path is: {s3dir}",
    "\nNo relevant files found under there.",
    "\nPlease check that path with `aws s3 ls`, and try to adjust page size."
  )
  assertthat::assert_that(nrow(d) > 0, msg = msg)
  d <- d |>
    dplyr::mutate(
      s3path_minus_s3dir = sub(glue("{s3dir}/"), "", .data$path),
      s3path_minus_s3dir_outdir = file.path(outdir, dirname(.data$s3path_minus_s3dir)) |>
        fs::dir_create() |>
        normalizePath(),
      localpath = file.path(.data$s3path_minus_s3dir_outdir, .data$bname),
      s3path = .data$path
    ) |>
    dplyr::select("type", "bname", "size", "lastmodified", "localpath", "s3path")
  tot_size <- d |>
    dplyr::summarise(tot_size = sum(.data$size)) |>
    dplyr::pull(tot_size)
  # download recognisable dracarys files to outdir/<mirrored-cloud-path>/{bname}
  if (!dryrun) {
    txt <- paste0(
      "{e('arrow_heading_down')} {nrow(d)} files ({tot_size}): {.file {s3dir}}\n"
    )
    cli::cli_alert_info(txt)
    res <- d |>
      dplyr::rowwise() |>
      dplyr::mutate(
        s3bucket = sub("s3://(.*?)/.*", "\\1", .data$s3path),
        s3key = sub("s3://(.*?)/(.*)", "\\2", .data$s3path),
        dl = list(
          s3$download_file(
            Bucket = .data$s3bucket, Key = .data$s3key, Filename = .data$localpath
          )
        ),
        localpath = normalizePath(.data$localpath)
      ) |>
      dplyr::ungroup() |>
      dplyr::select("type", "bname", "size", "lastmodified", "localpath", "s3path")
    return(res)
  } else {
    cli::cli_alert_info("{date_log()} {e('camera')} Just list relevant files from {.file {s3dir}}")
    d |>
      dplyr::select("type", "bname", "size", "lastmodified", "s3path", localpath2be = "localpath") |>
      as.data.frame() |>
      print()
  }
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
s3_file_presignedurl <- function(client, s3path, expiry_seconds = 604800) {
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
