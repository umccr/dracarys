#' s3_files_list_filter_relevant("s3://umccr-primary-data-prod/Accreditation/ALLOCATE-134131/WGS/2021-07-26/umccrised/ALLOCATE-134131__ALLOCATE-134131_MDx150892_Missing/cancer_report_tables/", presign = TRUE)
s3_files_list_filter_relevant <- function(s3dir, pattern = NULL, page_size = 1000, max_items = 1000, presign = FALSE) {
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
      date1 = .data$LastModified,
      size = fs::as_fs_bytes(.data$Size)
    ) |>
    dplyr::rowwise() |>
    dplyr::mutate(
      bname = basename(.data$path),
      type = purrr::map_chr(.data$bname, match_regex)
    ) |>
    dplyr::ungroup() |>
    dplyr::filter(!is.na(.data$type), grepl(pattern, .data$type)) |>
    dplyr::select(path, date1, size, type)

  if (presign) {
    d <- d |>
      dplyr::rowwise() |>
      dplyr::mutate(presigned_url = s3_file_presignedurl(.data$path)) |>
      dplyr::ungroup()
  }
  d
}

s3_file_presignedurl <- function(s3path, expiry_seconds = 3600) {
  p <- system(glue("aws s3 presign {s3path} --expires-in {expiry_seconds}"), intern = TRUE)
  p
}

# search for files on S3
s3_search <- function(search, rows) {
  au_tz <- "Australia/Melbourne"
  utc_tz <- "UTC"
  base_url <- "https://api.portal.prod.umccr.org/iam/s3"
  url1 <- utils::URLencode(glue::glue("{base_url}?rowsPerPage={rows}&search={search}"))
  awscurl_cmd <- glue::glue(
    "awscurl '{url1}' ",
    "--header 'Accept: application/json'"
  )
  message(glue::glue("Running {awscurl_cmd}"))
  j <- system(awscurl_cmd, intern = TRUE)
  date_fmt <- "%Y-%m-%dT%H:%M:%S"
  d <- j |>
    jsonlite::fromJSON() |>
    purrr::pluck("results") |>
    tibble::as_tibble()
  d |>
    dplyr::mutate(
      date1 = as.POSIXct(.data$last_modified_date, tz = utc_tz, format = date_fmt),
      date_aest = lubridate::with_tz(.data$date1, tz = au_tz)
    ) |>
    dplyr::select(path = key, bucket, size, date_aest, id, unique_hash)
}
