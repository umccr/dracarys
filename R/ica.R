#' List Relevant Files In GDS Directory
#'
#' Lists relevant files in a GDS directory.
#'
#' @param gdsdir GDS directory.
#' @param token ICA access token.
#' @param pattern Pattern to further filter the returned file type tibble.
#' @param include_url Include presigned URLs to all files within the GDS directory (def: FALSE).
#' @param page_size Page size (def: 100).
#' @param regexes Tibble with regex and function name.
#'
#' @return A tibble with type, bname, size, file_id, path, and presigned URL.
#' @export
gds_files_list_filter_relevant <- function(gdsdir, token, pattern = NULL, include_url = FALSE, page_size = 100, regexes = DR_FILE_REGEX) {
  pattern <- pattern %||% ".*" # keep all recognisable files by default
  cols_sel <- c("type", "bname", "size", "file_id", "path", "presigned_url")
  d <- dracarys::gds_files_list(gdsdir, token, include_url = include_url, page_size = page_size) |>
    dplyr::rowwise() |>
    dplyr::mutate(type = purrr::map_chr(.data$bname, \(x) match_regex(x, regexes))) |>
    dplyr::ungroup() |>
    dplyr::filter(!is.na(.data$type), grepl(pattern, .data$type)) |>
    dplyr::select(dplyr::any_of(cols_sel))
  d
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
#' @param page_size Page size (def: 10).
#' @param include_url Include presigned URLs to all files within the GDS directory (def: FALSE).
#' @param page_token Page token (def: NULL). Used internally for recursion.
#' @param no_recurse Do not recurse through the file list i.e. just give the first <page_size> items
#' without recursing further down the list using <page_token>.
#' @param recursive Should files be returned recursively _in and under_ the specified
#' GDS directory, or _only directly in_ the specified GDS directory (def: TRUE).
#'
#' @return Tibble with file basename, file size, file full data path, file dir name.
#' @examples
#' \dontrun{
#' gdsdir <- file.path(
#'   "gds://production/primary_data",
#'   "240322_A00130_0290_BH5HLLDSXC/20240323f56ec5a5/WGS_TsqNano"
#' )
#' gdsdir <- file.path(
#'   "gds://bssh.acddbfda498038ed99fa94fe79523959/Runs",
#'   "240322_A00130_0290_BH5HLLDSXC_r.3TbcOsEKZUyetygkqIOXcg/InterOp"
#' )
#' gdsdir <- file.path(
#'   "gds://production/analysis_data/SBJ00699/umccrise",
#'   "202203277dcf8562/L2200352__L2100146/SBJ00699__MDX220105/coverage"
#' )
#' token <- ica_token_validate()
#' page_size <- 11
#' include_url <- TRUE
#' page_token <- NULL
#' no_recurse <- TRUE
#' recursive <- NULL
#' gds_files_list(gdsdir, token, page_size, include_url, no_recurse, page_token, recursive)
#' }
#' @export
gds_files_list <- function(gdsdir, token, page_size = NULL, include_url = FALSE,
                           no_recurse = TRUE, page_token = NULL, recursive = NULL) {
  assertthat::assert_that(is.logical(no_recurse), is.logical(include_url))
  assertthat::assert_that(is.null(recursive) || is.logical(recursive))
  token <- ica_token_validate(token)
  assertthat::assert_that(grepl("^gds://", gdsdir))
  gdsdir_original <- gdsdir
  if (!grepl("/$", gdsdir)) {
    gdsdir <- glue("{gdsdir}/")
  }
  base_url <- "https://aps2.platform.illumina.com/v1"
  volname <- sub("gds://(.*?)/.*", "\\1", gdsdir)
  path2 <- sub("gds://(.*?)/(.*)", "\\2", gdsdir)
  page_size <- ifelse(is.null(page_size), "", glue("&pageSize={page_size}"))
  query_url <- glue("{base_url}/files?volume.name={volname}&path=/{path2}*{page_size}")
  if (include_url) {
    query_url <- glue("{query_url}&include=PresignedUrl")
  }
  if (!is.null(page_token)) {
    query_url <- glue("{query_url}&pageToken={page_token}")
  }
  if (!is.null(recursive)) {
    # without specifying recursive, it's true by default
    recursive <- ifelse(recursive, "true", "false")
    query_url <- glue("{query_url}&recursive={recursive}")
  }
  query_res <- httr::GET(
    query_url,
    httr::add_headers(Authorization = glue("Bearer {token}")),
    httr::accept_json()
  )
  j <- jsonlite::fromJSON(httr::content(x = query_res, type = "text", encoding = "UTF-8"), simplifyVector = FALSE)
  if (j[["itemCount"]] == 0) {
    if (likely_file(gdsdir_original)) {
      cli::cli_abort("{date_log()} ERROR: Is the input directory a file perhaps?\n{.file {gdsdir_original}}")
    }
    # if <somehow> there is a nextPageToken then abort, else continue
    if (!is.null(j[["nextPageToken"]])) {
      msg <- paste0(
        "{date_log()} ERROR: ",
        "No GDS files listed in the input directory. Please confirm you can ",
        "access the following GDS input directory with your token: ",
        "{.file {gdsdir_original}}"
      )
      cli::cli_abort(msg)
    }
  } # endif
  d <- j[["items"]] |>
    purrr::map(\(x) c(file_id = x[["id"]], path = x[["path"]], size = x[["sizeInBytes"]], presigned_url = x[["presignedUrl"]])) |>
    dplyr::bind_rows()
  if (nrow(d) == 0) {
    # We've iterated through all available items, and the next page has 0 items.
    # So dplyr::bind_rows(d, NULL) will return d.
    return(NULL)
  }
  res <- d |>
    dplyr::mutate(
      size = fs::as_fs_bytes(.data$size),
      bname = basename(.data$path),
      path = glue("gds://{volname}{.data$path}"),
      dname = basename(dirname(.data$path))
    ) |>
    dplyr::select(dplyr::any_of(c("file_id", "bname", "size", "path", "dname", "presigned_url")))
  if (!is.null(j[["nextPageToken"]]) && !no_recurse) {
    res2 <- gds_files_list(
      gdsdir = gdsdir, token = token, page_size = NULL,
      include_url = include_url, no_recurse = FALSE, page_token = j[["nextPageToken"]],
      recursive = NULL
    )
    res <- dplyr::bind_rows(res, res2)
  }
  res
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
#' @param page_size Page size (def: 100).
#' @param pattern Pattern to further filter the returned file type tibble.
#' @param dryrun If TRUE, just list the files that will be downloaded (don't
#' download them).
#' @param regexes Tibble with regex and function name.
#' @param recursive Should files be returned recursively _in and under_ the specified
#' GDS directory (TRUE), or _only directly in_ the specified GDS directory (FALSE) (def: TRUE).
#'
#' @export
dr_gds_download <- function(gdsdir, outdir, token, page_size = 100, pattern = NULL,
                            dryrun = FALSE, regexes = DR_FILE_REGEX, recursive = NULL) {
  e <- emojifont::emoji
  fs::dir_create(outdir)
  d <- gds_files_list(
    gdsdir = gdsdir, token = token, page_size = page_size,
    no_recurse = FALSE, recursive = recursive
  ) |>
    dplyr::mutate(type = purrr::map_chr(.data$bname, \(x) match_regex(x, regexes))) |>
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

ica_token_exp <- function(token = Sys.getenv("ICA_ACCESS_TOKEN")) {
  l <- jose::jwt_split(token)
  structure(l$payload$exp, class = c("POSIXct", "POSIXt"))
}

likely_file <- function(x) {
  e <- c(
    "txt", "tsv", "csv", "html", "json", "stdout", "stderr", "stdouterr",
    "log", "vcf", "gz", "bam", "bai"
  )
  tolower(tools::file_ext(x)) %in% e
}
