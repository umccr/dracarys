#' List Files in ICAv1 GDS Directory
#'
#' Lists files in a GDS directory.
#'
#' @param gdsdir Full path to GDS directory.
#' @param token ICA access token (def: $ICA_ACCESS_TOKEN env var).
#' @param page_size Page size (def: 10 via ICA API).
#' @param include_url Include presigned URLs to all files within the GDS directory (def: FALSE via ICA API).
#' @param page_token Page token (def: NULL). Used internally for recursion.
#' @param no_recurse Do not recurse through the file list i.e. just give the first <page_size> items
#' without recursing further down the list using <page_token>.
#' @param recursive Should files be returned recursively _in and under_ the specified
#' GDS directory, or _only directly in_ the specified GDS directory (def: TRUE via ICA API).
#'
#' @return A tibble with file ID, basename, size, last modified timestamp,
#' full GDS path, and presigned URL if requested.
#' @examples
#' \dontrun{
#' gdsdir <- file.path(
#'   "gds://production/analysis_data/SBJ00699/umccrise",
#'   "202203277dcf8562/L2200352__L2100146/SBJ00699__MDX220105/coverage"
#' )
#' token <- ica_token_validate()
#' page_size <- 11
#' include_url <- F
#' page_token <- NULL
#' no_recurse <- TRUE
#' recursive <- NULL
#' gds_list_files_dir(gdsdir, token, page_size, include_url, no_recurse, page_token, recursive)
#' }
#' @export
gds_list_files_dir <- function(gdsdir, token = Sys.getenv("ICA_ACCESS_TOKEN"), page_size = NULL,
                               include_url = FALSE, no_recurse = TRUE, page_token = NULL,
                               recursive = NULL) {
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
    if (gds_likely_file(gdsdir_original)) {
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
    purrr::map(\(x) c(
      file_id = x[["id"]],
      path = x[["path"]],
      size = x[["sizeInBytes"]],
      lastmodified = x[["timeModified"]],
      presigned_url = x[["presignedUrl"]]
    )) |>
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
      path = glue("gds://{volname}{.data$path}")
    ) |>
    dplyr::select(dplyr::any_of(c("bname", "size", "lastmodified", "file_id", "path", "presigned_url")))

  if (!is.null(j[["nextPageToken"]]) && !no_recurse) {
    res2 <- gds_list_files_dir(
      gdsdir = gdsdir, token = token, page_size = NULL,
      include_url = include_url, no_recurse = FALSE, page_token = j[["nextPageToken"]],
      recursive = NULL
    )
    res <- dplyr::bind_rows(res, res2)
  }
  res
}

#' List Relevant Files In ICAv1 GDS Directory
#'
#' Lists relevant files in a GDS directory.
#'
#' @inheritParams gds_list_files_dir
#' @param pattern Pattern to further filter the returned file type tibble.
#' @param regexes Tibble with `regex` and `fun`ction name (see example).
#' @return A tibble with file type, basename, size, last modified timestamp, file_id, full path,
#' and presigned URL if requested.
#' @examples
#' \dontrun{
#' regexes <- tibble::tibble(regex = "multiqc_data\\.json$", fun = "MultiqcJsonFile")
#' gdsdir <- "gds://production/analysis_data/SBJ01155/umccrise/202408300c218043/L2101566__L2101565"
#' gds_list_files_filter_relevant(gdsdir)
#' }
#' @export
gds_list_files_filter_relevant <- function(gdsdir, pattern = NULL, regexes = DR_FILE_REGEX,
                                           token = Sys.getenv("ICA_ACCESS_TOKEN"),
                                           page_size = 100, include_url = FALSE,
                                           no_recurse = TRUE, page_token = NULL,
                                           recursive = NULL) {
  pattern <- pattern %||% ".*" # keep all recognisable files by default
  assertthat::assert_that(all(c("regex", "fun") %in% colnames(regexes)))
  cols_sel <- c("type", "bname", "size", "lastmodified", "file_id", "path", "presigned_url")
  d <- dracarys::gds_list_files_dir(
    gdsdir = gdsdir, token = token, page_size = page_size, include_url = include_url,
    no_recurse = no_recurse, page_token = page_token, recursive = recursive
  ) |>
    dplyr::rowwise() |>
    dplyr::mutate(type = purrr::map_chr(.data$path, \(x) match_regex(x, regexes))) |>
    dplyr::ungroup() |>
    dplyr::filter(!is.na(.data$type), grepl(pattern, .data$type)) |>
    dplyr::select(dplyr::any_of(cols_sel))
  d
}

#' dracarys GDS Download
#'
#' Download only GDS files that can be processed by dracarys.
#'
#' @inheritParams gds_list_files_dir
#' @inheritParams gds_list_files_filter_relevant
#' @param outdir Local output directory.
#' @param dryrun If TRUE, just list the files that will be downloaded (don't
#' download them).
#' @examples
#' \dontrun{
#' gdsdir <- "gds://production/analysis_data/SBJ01155/umccrise/202408300c218043/L2101566__L2101565"
#' outdir <- sub("gds:/", "~/icav1/g", gdsdir)
#' regexes <- tibble::tribble(
#'   ~regex, ~fun,
#'   "multiqc_data\\.json$", "MultiqcJsonFile",
#'   "-somatic\\.pcgr\\.json\\.gz$", "pcgrjson"
#' )
#' dr_gds_download(gdsdir = gdsdir, outdir = outdir, regexes = regexes, dryrun = T)
#' }
#'
#' @export
dr_gds_download <- function(gdsdir, outdir, token = Sys.getenv("ICA_ACCESS_TOKEN"),
                            pattern = NULL, page_size = 100, dryrun = FALSE,
                            regexes = DR_FILE_REGEX, recursive = NULL) {
  e <- emojifont::emoji
  fs::dir_create(outdir)
  d <- gds_list_files_filter_relevant(
    gdsdir = gdsdir, pattern = pattern, regexes = regexes,
    token = token, page_size = page_size, include_url = FALSE,
    no_recurse = FALSE, page_token = NULL,
    recursive = recursive
  )
  d <- d |>
    dplyr::mutate(
      gdspath_minus_gdsdir = sub(glue("{gdsdir}/"), "", .data$path),
      gdspath_minus_gdsdir_outdir = fs::dir_create(
        file.path(outdir, dirname(.data$gdspath_minus_gdsdir))
      ),
      localpath = file.path(.data$gdspath_minus_gdsdir_outdir, .data$bname),
      gdspath = .data$path
    ) |>
    dplyr::select("type", "bname", "size", "lastmodified", "file_id", "localpath", "gdspath")
  # download recognisable dracarys files to outdir/<mirrored-cloud-path>/{bname}
  tot_size <- d |>
    dplyr::summarise(tot_size = sum(.data$size)) |>
    dplyr::pull(tot_size)
  if (!dryrun) {
    txt <- paste0(
      "{e('arrow_heading_down')} {nrow(d)} files ({tot_size}): {.file {gdsdir}}\n"
    )
    cli::cli_alert_info(txt)
    res <- d |>
      dplyr::rowwise() |>
      dplyr::mutate(
        dl = gds_file_download_api(
          gds_fileid = .data$file_id, out_file = .data$localpath, token = token
        ),
        localpath = normalizePath(.data$localpath)
      ) |>
      dplyr::select("type", "bname", "size", "lastmodified", "localpath", "gdspath", "file_id")
    return(res)
  } else {
    cli::cli_alert_info("{date_log()} {e('camera')} Just list relevant files from {.file {gdsdir}}")
    d |>
      dplyr::select("type", "bname", "size", "lastmodified", "gdspath", "file_id", localpath2be = "localpath") |>
      as.data.frame() |>
      print()
  }
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
  normalizePath(out_file)
}

#' GDS File Download via CLI
#'
#' @param gds Full path to GDS file.
#' @param out Path to output file.
#' @param token ICA access token (def: $ICA_ACCESS_TOKEN env var).
#' @export
gds_file_download_cli <- function(gds, out, token = Sys.getenv("ICA_ACCESS_TOKEN")) {
  token <- ica_token_validate(token)
  system(glue("ica files download {gds} {out} --access-token {token}"))
}
