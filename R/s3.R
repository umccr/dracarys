#' List Relevant Files In AWS S3 Directory
#'
#' Lists relevant files in an AWS S3 directory.
#'
#' @param s3dir GDS directory.
#' @param pattern Pattern to further filter the returned file type tibble.
#' @param page_size The size of each page to get in the AWS service call (def: 1000).
#' @param max_items The total number of items to return in the command’s output (def: 1000).
#' @param presign Include presigned URLs (def: FALSE).
#' @param expiry_sec Number of seconds the presigned URL will be valid for (if generated) (def: 43200 (12hrs)).
#'
#' @return A tibble with path, date, file size, file type, and presigned URL if requested.
#' @examples
#' \dontrun{
#' s3dir <- "s3://umccr-primary-data-prod/cancer_report_tables"
#' s3_files_list_filter_relevant(s3dir = s3dir, presign = TRUE)
#' }
#' @export
s3_files_list_filter_relevant <- function(s3dir, pattern = NULL, page_size = 1000, max_items = 1000, presign = FALSE, expiry_sec = 43200) {
  assertthat::assert_that(grepl("^s3://", s3dir), rlang::is_logical(presign))
  pattern <- pattern %||% ".*" # keep all recognisable files by default
  b <- sub("s3://(.*?)/.*", "\\1", s3dir)
  p <- sub("s3://(.*?)/(.*)", "\\2", s3dir)
  cmd <- glue(
    "aws --output json s3api list-objects-v2 --bucket {b} --prefix {p} ",
    "--max-items {max_items} --page-size {page_size}"
  )
  l <- system(cmd, intern = TRUE)
  j <- jsonlite::fromJSON(l)
  assertthat::assert_that("Contents" %in% names(j))
  d <- j[["Contents"]] |>
    tibble::as_tibble() |>
    dplyr::mutate(
      path = glue("s3://{b}/{.data$Key}"),
      date_utc = .data$LastModified,
      size = fs::as_fs_bytes(.data$Size)
    ) |>
    dplyr::rowwise() |>
    dplyr::mutate(
      bname = basename(.data$path),
      type = purrr::map_chr(.data$bname, match_regex)
    ) |>
    dplyr::ungroup() |>
    dplyr::filter(!is.na(.data$type), grepl(pattern, .data$type)) |>
    dplyr::select("path", "date_utc", "size", "type")

  if (presign) {
    d <- d |>
      dplyr::rowwise() |>
      dplyr::mutate(presigned_url = s3_file_presignedurl(.data$path, expiry_seconds = expiry_sec)) |>
      dplyr::ungroup()
  }
  d
}

s3_file_presignedurl <- function(s3path, expiry_seconds = 3600) {
  p <- system(glue("aws s3 presign {s3path} --expires-in {expiry_seconds}"), intern = TRUE)
  p
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
