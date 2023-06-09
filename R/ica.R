#' List Relevant Files In GDS Directory
#'
#' Lists relevant files in a GDS directory.
#'
#' @param gdsdir GDS directory.
#' @param token ICA access token.
#' @param pattern Pattern to further filter the returned file type tibble.
#' @param include Use PresignedUrl to include presigned URLs to all files within
#' the GDS directory.
#'
#' @return A tibble with type, bname, size, file_id, path, and presigned URL.
#' @export
gds_files_list_filter_relevant <- function(gdsdir, token, pattern = NULL, include = NULL) {
  pattern <- pattern %||% ".*" # keep all recognisable files by default
  cols_sel <- c("type", "bname", "size", "file_id", "path", "presigned_url")
  d <- dracarys::gds_files_list(gdsdir, token, include = include) |>
    dplyr::rowwise() |>
    dplyr::mutate(type = purrr::map_chr(.data$bname, match_regex)) |>
    dplyr::ungroup() |>
    dplyr::filter(!is.na(.data$type), grepl(pattern, .data$type)) |>
    dplyr::select(dplyr::any_of(cols_sel))
}

#' GDS File Presigned URL
#'
#' Returns presigned URL of given GDS file.
#'
#' @param gds_fileid GDS file ID.
#' @param token ICA access token (def: $ICA_ACCESS_TOKEN env var).
#' @return Presigned URL if valid.
#' @export
gds_file_presignedurl <- function(gds_fileid, token) {
  token <- ica_token_validate(token)
  base_url <- "https://aps2.platform.illumina.com/v1"
  url <- glue("{base_url}/files/{gds_fileid}")
  res <- httr::GET(
    url,
    httr::add_headers(Authorization = glue("Bearer {token}")),
    httr::accept_json()
  )
  presigned_url <- jsonlite::fromJSON(httr::content(x = res, as = "text", encoding = "UTF-8"), simplifyVector = FALSE)[["presignedUrl"]]
  assertthat::assert_that(grepl("^https://stratus-gds-aps2.s3.ap-southeast-2.amazonaws.com", presigned_url))
  presigned_url
}

#' GDS File Download via API
#'
#' @param gds_fileid GDS file ID.
#' @param out_file Path to output file.
#' @param token ICA access token (def: $ICA_ACCESS_TOKEN env var).
#'
#' @examples
#' \dontrun{
#' gds_fileid <- "fil.f9aa2ba7af0c4330095d08dadd2e16b0"
#' out <- tempfile()
#' token <- Sys.getenv("ICA_ACCESS_TOKEN")
#' }
#' @export
gds_file_download_api <- function(gds_fileid, out_file, token) {
  presigned_url <- gds_file_presignedurl(gds_fileid, token)
  # keep quiet instead of logging presigned urls
  status_code <- utils::download.file(url = presigned_url, destfile = out_file, quiet = TRUE)
  assertthat::assert_that(status_code == 0)
  out_file
}

#' GDS File Download via CLI
#'
#' @param gds Full path to GDS file.
#' @param out Path to output file.
#' @param token ICA access token (def: $ICA_ACCESS_TOKEN env var).
#' @export
gds_file_download <- function(gds, out, token = Sys.getenv("ICA_ACCESS_TOKEN")) {
  token <- ica_token_validate(token)
  system(glue("ica files download {gds} {out} --access-token {token}"))
}

#' GDS Files List
#'
#' List files on ICA GDS filesystem.
#'
#' @param gdsdir Full path to GDS directory.
#' @param token ICA access token (def: $ICA_ACCESS_TOKEN env var).
#'
#' @return Tibble with file basename, file size, file full data path, file dir name.
#' @export
gds_files_list <- function(gdsdir, token, include = NULL) {
  token <- ica_token_validate(token)
  assertthat::assert_that(grepl("^gds://", gdsdir))
  gdsdir_original <- gdsdir
  if (!grepl("/$", gdsdir)) {
    gdsdir <- glue("{gdsdir}/")
  }
  base_url <- "https://aps2.platform.illumina.com/v1"
  volname <- sub("gds://(.*?)/.*", "\\1", gdsdir)
  path2 <- sub("gds://(.*?)/(.*)", "\\2", gdsdir)
  query_url <- glue("{base_url}/files?volume.name={volname}&path=/{path2}*&pageSize=100")
  if (!is.null(include)) {
    assertthat::assert_that(include == "PresignedUrl")
    query_url <- glue("{query_url}&include=PresignedUrl")
  }

  res <- httr::GET(
    query_url,
    httr::add_headers(Authorization = glue("Bearer {token}")),
    httr::accept_json()
  )
  j <- jsonlite::fromJSON(httr::content(x = res, type = "text", encoding = "UTF-8"), simplifyVector = FALSE)
  if (j[["itemCount"]] == 0) {
    if (likely_file(gdsdir_original)) {
      cli::cli_abort("{date_log()} ERROR: Is the input directory a file perhaps?\n{.file {gdsdir_original}}")
    }
    msg <- paste0(
      "{date_log()} ERROR: ",
      "No GDS files listed in the input directory. Please confirm you can ",
      "access the following GDS input directory with your token: ",
      "{.file {gdsdir_original}}"
    )
    cli::cli_abort(msg)
  } # endif
  d <- purrr::map_df(j[["items"]], function(x) c(file_id = x[["id"]], path = x[["path"]], size = x[["sizeInBytes"]], presigned_url = x[["presignedUrl"]]))
  d |>
    dplyr::mutate(
      size = fs::as_fs_bytes(.data$size),
      bname = basename(.data$path),
      path = glue("gds://{volname}{.data$path}"),
      dname = basename(dirname(.data$path))
    ) |>
    dplyr::select(dplyr::any_of(c("file_id", "bname", "size", "path", "dname", "presigned_url")))
}

#' List GDS Volumes
#'
#' Lists GDS volumes accessible by the provided ICA token.
#'
#' @param token ICA access token (def: $ICA_ACCESS_TOKEN env var).
#' @param page_size Page size (def: 10).
#'
#' @return A tibble with vol name and vol id.
#' @export
gds_volumes_list <- function(token, page_size = 10) {
  token <- ica_token_validate(token)
  base_url <- "https://aps2.platform.illumina.com/v1"
  query_url <- glue("{base_url}/volumes?pageSize={page_size}")

  res <- httr::GET(
    query_url,
    httr::add_headers(Authorization = glue("Bearer {token}")),
    httr::accept_json()
  )
  j <- jsonlite::fromJSON(httr::content(x = res, type = "text", encoding = "UTF-8"), simplifyVector = FALSE)
  purrr::map_df(j[["items"]], function(x) c(name = x[["name"]], id = x[["id"]]))
}


#' dracarys GDS Download
#'
#' Download only GDS files that can be processed by dracarys.
#'
#' @param gdsdir Full path to GDS directory.
#' @param outdir Path to output directory.
#' @param token ICA access token (def: $ICA_ACCESS_TOKEN env var).
#' @param pattern Pattern to further filter the returned file type tibble.
#' @param dryrun If TRUE, just list the files that will be downloaded (don't
#' download them).
#' @export
dr_gds_download <- function(gdsdir, outdir, token,
                            pattern = NULL, dryrun = FALSE) {
  e <- emojifont::emoji
  fs::dir_create(outdir)
  d <- gds_files_list(gdsdir = gdsdir, token = token) |>
    dplyr::mutate(type = purrr::map_chr(.data$bname, match_regex)) |>
    dplyr::select("file_id", "dname", "type", "size", "path", "bname")

  # download recognisable dracarys files to outdir/{bname}
  pattern <- pattern %||% ".*" # keep all recognisable files
  d_filt <- d |>
    dplyr::filter(!is.na(.data$type), grepl(pattern, .data$type)) |>
    dplyr::mutate(out = file.path(outdir, .data$bname))
  if (!dryrun) {
    cli::cli_alert_info("{date_log()} {e('arrow_heading_down')} Downloading files from {.file {gdsdir}}")
    d_filt |>
      dplyr::rowwise() |>
      dplyr::mutate(out_dl = gds_file_download_api(.data$file_id, .data$out, token))
  } else {
    cli::cli_alert_info("{date_log()} {e('camera')} Just list relevant files from {.file {gdsdir}}")
    d_filt |>
      dplyr::select("path", "type", "size") |>
      as.data.frame() |>
      print()
  }
}

#' Validate ICA access token
#'
#' Validates ICA access token by parsing it and checking its expiration date.
#' @param token ICA access token (def: $ICA_ACCESS_TOKEN env var).
#' @return Returns the token if valid, or else errors out.
#' @export
ica_token_validate <- function(token = Sys.getenv("ICA_ACCESS_TOKEN")) {
  # https://github.com/r-lib/jose/blob/429a46/R/jwt.R#L171
  .ica_token_check_expiration_time <- function(payload) {
    if (length(payload$exp)) {
      stopifnot("exp claim is a number" = is.numeric(payload$exp))
      expdate <- structure(payload$exp, class = c("POSIXct", "POSIXt"))
      if (expdate < (Sys.time() - 60)) {
        stop(paste("Token has expired on", expdate), call. = FALSE)
      }
    }
    if (length(payload$nbf)) {
      stopifnot("nbf claim is a number" = is.numeric(payload$nbf))
      nbfdate <- structure(payload$nbf, class = c("POSIXct", "POSIXt"))
      if (nbfdate > (Sys.time() + 60)) {
        stop(paste("Token is not valid before", nbfdate), call. = FALSE)
      }
    }
  }
  # giving a friendlier error msg in case this isn't even valid jwt
  tmp <- strsplit(token, ".", fixed = TRUE)[[1]]
  msg <- "The input token is not a valid JWT"
  assertthat::assert_that(length(tmp) %in% c(2, 3), msg = msg)
  l <- jose::jwt_split(token)
  .ica_token_check_expiration_time(l[["payload"]])
  token
}

likely_file <- function(x) {
  e <- c(
    "txt", "tsv", "csv", "html", "json", "stdout", "stderr", "stdouterr",
    "log", "vcf", "gz", "bam", "bai"
  )
  tolower(tools::file_ext(x)) %in% e
}
