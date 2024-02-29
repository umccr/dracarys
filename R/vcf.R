bcftools_installed <- function() {
  system("bcftools -v", ignore.stdout = TRUE) == 0
}

#' Parse VCF with bcftools
#'
#' Parse VCF with bcftools.
#' Uses bcftools under the hood to do the heavy lifting with field splitting,
#' then converts the parsed character vector to a tibble.
#'
#' For VCFs with 0 variants, returns a tibble with 0 rows and proper number of
#' columns.
#'
#' @param vcf VCF with one or more samples.
#' @param only_pass Keep PASS variants only (def: TRUE).
#' @param alias Substitute sample names with S1/S2/... alias (def: TRUE).
#'
#' @return A tibble with all the main, FORMAT, and INFO fields detected in
#' the VCF header as columns.
#' @export
bcftools_parse_vcf <- function(vcf, only_pass = TRUE, alias = TRUE) {
  assertthat::assert_that(is.logical(only_pass), length(only_pass) == 1)
  assertthat::assert_that(bcftools_installed(), msg = "bcftools needs to be on the PATH.")
  if (is_url(vcf)) {
    vcf <- glue("'{vcf}'")
  }
  cmd_header <- glue("bcftools view -h {vcf}")
  h <- system(cmd_header, intern = TRUE)
  main <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER")
  get_samples <- function() {
    final_header_row <- length(h)
    x <- strsplit(h[final_header_row], "\t")[[1]]
    main_len <- length(main) + length(c("INFO", "FORMAT"))
    if (length(x) <= main_len) {
      msg <- paste(c(main, "INFO", "FORMAT"), collapse = ", ")
      stop(
        "What on earth. Check that your VCF has the following VCF columns, followed by sample columns:\n",
        msg, "\nCall me if everything seems fine on 1800-OMG-LOL."
      )
    }
    samples <- x[-(1:main_len)]
    nsamples <- length(samples)
    aliases <- paste0("S", seq_len(nsamples))
    # in case you want the full sample name instead
    if (!alias) {
      aliases <- samples
    }
    list(
      samples = samples,
      n = nsamples,
      aliases = aliases
    )
  }
  samp <- get_samples()
  # splits header sections into tbls, mostly to grab available FORMAT/INFO fields
  split_hdr <- function(pat) {
    h[grepl(pat, h)] |>
      tibble::as_tibble_col(column_name = "x") |>
      dplyr::mutate(x = sub(pat, "", .data$x)) |>
      # Description is likely to have a comma
      tidyr::separate_wider_delim("x", delim = ",", names = c("ID", "Number", "Type", "Description"), too_many = "merge") |>
      dplyr::mutate(
        ID = sub("ID=", "", .data$ID),
        Number = sub("Number=", "", .data$Number),
        Type = sub("Type=", "", .data$Type),
        Description = sub("Description=\\\"(.*)\\\">", "\\1", .data$Description)
      )
  }
  fmt <- split_hdr("##FORMAT=<") |> dplyr::pull(.data$ID)
  info <- split_hdr("##INFO=<") |> dplyr::pull(.data$ID)
  main_cols <- paste0("%", main) |>
    paste(collapse = "\\t")
  info_cols <- paste0("%INFO/", info, collapse = "\\t")
  fmt_cols <- paste0("[\\t", paste0(paste0("%", fmt), collapse = "\\t"), "]\\n")
  q <- paste0(main_cols, "\\t", info_cols, fmt_cols)
  include_pass <- ""
  if (only_pass) {
    include_pass <- "-i 'FILTER=\"PASS\" || FILTER=\".\"'"
  }
  cmd_body <- glue("bcftools query -f \"{q}\" {vcf} {include_pass}")
  # create column names using the main columns, an INFO prefix for the INFO
  # columns, and a S1/2/.._X prefix for the sample columns.
  cnames <- c(
    main,
    paste0("INFO_", info),
    paste0(rep(samp$aliases, each = length(fmt)), "_", fmt)
  )
  # handle empty VCF - fread warns about size 0
  suppressWarnings({
    first_row <- data.table::fread(cmd = cmd_body, sep = "\t", na.strings = ".", nrows = 1)
  })
  if (nrow(first_row) == 0) {
    return(empty_tbl(cnames = cnames))
  }
  data.table::fread(
    cmd = cmd_body,
    header = FALSE,
    sep = "\t",
    col.names = cnames,
    data.table = FALSE,
    na.strings = "."
  ) |>
    tibble::as_tibble()
}

#' Parse VCF regions with bcftools
#'
#' Parses VCF regions with bcftools. The VCF subset is written to a temporary
#' file in the local filesystem, then parsed into a tibble object.
#'
#' @param vcf Path to VCF. Can be S3, http or local. If presigned URL, need to
#' also concatenate the VCF index as in 'vcf_url##idx##vcfi_url'.
#' @param r Character vector of regions to subset (e.g. c('chr1:123-456', 'chr2:789-1000'))
#' @param only_pass Keep PASS variants only (def: TRUE).
#'
#' @return A tibble with all the main, FORMAT, and INFO fields detected in
#' the VCF header as columns, for the regions specified in `r` (if any).
#'
#' @examples
#' \dontrun{
#' vcf_local <- here::here("MergedSmallVariants.vcf.gz")
#' r <- c("chr1:115256529-115256529", "chr2:29443613-29443613")
#' bcftools_parse_vcf_regions(vcf_local, r)
#' }
#' @export
bcftools_parse_vcf_regions <- function(vcf, r, only_pass = TRUE) {
  assertthat::assert_that(is.character(r))
  assertthat::assert_that(bcftools_installed(), msg = "bcftools needs to be on the PATH.")
  if (is_url(vcf)) {
    vcf <- glue("'{vcf}'")
  }
  r <- paste(r, collapse = ",")
  # write to temp then parse into R
  out <- tempfile(fileext = ".vcf")
  cmd <- glue("bcftools view {vcf} -r {r} > {out}")
  system(cmd, intern = TRUE)
  bcftools_parse_vcf(vcf = out, only_pass = only_pass)
}
