#' List FASTQs In GDS Directory
#'
#' @param gdsdir GDS directory.
#' @param token ICA access token.
#' @param include_url Include presigned URLs to all files within the GDS directory.
#' @param page_size Page size.
#'
#' @return A tibble with type, bname, size, file_id, path, and presigned URL.
#'
#' @examples
#' \dontrun{
#' prim <- "gds://production/primary_data"
#' run <- "240719_A00130_0323_BHMCYHDSXC/202407205bad380d/BiModal_BM-5L"
#' gdsdir <- file.path(prim, run)
#' token <- Sys.getenv("ICA_ACCESS_TOKEN")
#' include_url <- F
#' page_size <- 100
#' gds_files_list_fastq(gdsdir, token, include_url, page_size)
#' }
#' @export
gds_files_list_fastq <- function(gdsdir, token, include_url = FALSE, page_size = 100) {
  fq_regex <- tibble::tribble(
    ~regex, ~fun,
    "fastq\\.gz$", "FASTQ"
  )
  g <- gds_list_files_filter_relevant(
    gdsdir = gdsdir, pattern = NULL, regexes = fq_regex,
    token = token, page_size = page_size, include_url = include_url
  )
  g |>
    dplyr::mutate(
      size_chr = as.character(.data$size),
      size_num = as.numeric(.data$size)
    ) |>
    dplyr::select(
      "type", "bname", "size", "lastmodified", "size_chr", "size_num", "file_id", "path"
    )
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

gds_likely_file <- function(x) {
  e <- c(
    "txt", "tsv", "csv", "html", "json", "stdout", "stderr", "stdouterr",
    "log", "vcf", "gz", "bam", "bai"
  )
  tolower(tools::file_ext(x)) %in% e
}
