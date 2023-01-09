#' GDS File Download via API
#'
#' @param gds_fileid GDS file ID.
#' @param out_file Path to output file.
#' @param token ICA access token (by default uses $ICA_ACCESS_TOKEN env var).
#'
#' @examples
#' \dontrun{
#' gds_fileid <- "fil.f9aa2ba7af0c4330095d08dadd2e16b0"
#' out <- tempfile()
#' token <- Sys.getenv("ICA_ACCESS_TOKEN")
#' }
#' @export
gds_file_download_api <- function(gds_fileid, out_file, token) {
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
  # keep quiet instead of logging presigned urls
  utils::download.file(url = presigned_url, destfile = out_file, quiet = TRUE)
  out_file
}

#' GDS File Download via CLI
#'
#' @param gds Full path to GDS file.
#' @param out Path to output file.
#' @param token ICA access token (by default uses $ICA_ACCESS_TOKEN env var).
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
#' @param token ICA access token (by default uses $ICA_ACCESS_TOKEN env var).
#'
#' @return Tibble with file basename, file size, file full data path, file dir name.
#' @export
gds_files_list <- function(gdsdir, token) {
  token <- ica_token_validate(token)
  assertthat::assert_that(grepl("^gds://", gdsdir))
  if (!grepl("/$", gdsdir)) {
    gdsdir <- glue("{gdsdir}/")
  }
  base_url <- "https://aps2.platform.illumina.com/v1"
  volname <- sub("gds://(.*?)/.*", "\\1", gdsdir)
  path2 <- sub("gds://(.*?)/(.*)", "\\2", gdsdir)
  query_url <- glue("{base_url}/files?volume.name={volname}&path=/{path2}*&pageSize=100")

  res <- httr::GET(
    query_url,
    httr::add_headers(Authorization = glue("Bearer {token}")),
    httr::accept_json()
  )
  j <- jsonlite::fromJSON(httr::content(x = res, type = "text", encoding = "UTF-8"), simplifyVector = FALSE)[["items"]]
  d <- purrr::map_df(j, function(x) c(file_id = x[["id"]], path = x[["path"]], size = x[["sizeInBytes"]]))
  d |>
    dplyr::mutate(
      size = fs::as_fs_bytes(.data$size),
      bname = basename(.data$path),
      path = glue("gds://{volname}{.data$path}"),
      dname = basename(dirname(.data$path))
    ) |>
    dplyr::select("file_id", "bname", "size", "path", "dname")
}


#' dracarys GDS Download
#'
#' Download only GDS files that can be processed by dracarys.
#'
#' @param gdsdir Full path to GDS directory.
#' @param outdir Path to output directory.
#' @param token ICA access token (by default uses $ICA_ACCESS_TOKEN env var).
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
  if (is.null(pattern)) {
    pattern <- ".*" # keep all recognisable files
  }
  d_filt <- d |>
    dplyr::filter(!is.na(.data$type), grepl(pattern, .data$type)) |>
    dplyr::mutate(out = file.path(outdir, .data$bname))
  if (!dryrun) {
    cli::cli_alert_info("{date_log()} {e('arrow_heading_down')}  Downloading files from {.file {gdsdir}}")
    d_filt |>
      dplyr::rowwise() |>
      dplyr::mutate(out_dl = gds_file_download_api(.data$file_id, .data$out, token))
  } else {
    cli::cli_alert_info("{date_log()} {e('camera')} Just list relevant files from {.file {gdsdir}}")
    print(d_filt)
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
  l <- jose::jwt_split(token)
  .ica_token_check_expiration_time(l[["payload"]])
  token
}
