tidy_add_args <- function(subp) {
  tidy <- subp$add_parser("tidy", help = "DRAGEN Output Tidying")
  tidy$add_argument("-j", "--json", help = "Path to 'multiqc_data.json'.", required = TRUE)
  tidy$add_argument("-o", "--outdir", help = "Output directory for results.", required = TRUE)
  tidy$add_argument("-p", "--prefix", help = "Prefix name for output files.", required = TRUE)
  tidy$add_argument("-q", "--quiet", help = "Suppress log printing during rendering.", action = "store_true")
}

tidy_parse_args <- function(args) {
  assertthat::assert_that(
    file.exists(args$json), is.character(args$prefix), nchar(args$prefix) > 0,
    is.character(args$outdir), nchar(args$outdir) > 0
  )

  prefix <- args$prefix
  json <- normalizePath(args$json)
  mkdir(args$outdir)
  outdir <- normalizePath(args$outdir)
  quiet <- args$quiet


  if (quiet) {
    suppressMessages(dracarys_tidy_multiqc(json, prefix, outdir))
  } else {
    dracarys_tidy_multiqc(json, prefix, outdir)
  }
}
