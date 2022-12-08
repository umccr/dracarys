multiqc_add_args <- function(subp) {
  multiqc <- subp$add_parser("multiqc", help = "Tidy MultiQC JSON Output")
  multiqc$add_argument("-j", "--json", help = glue("{emoji('poop')} Path to 'multiqc_data.json'."), required = TRUE)
  multiqc$add_argument("-o", "--outdir", help = glue("{emoji('gift')} Output results to this directory."), required = TRUE)
  multiqc$add_argument("-p", "--prefix", help = glue("{emoji('sparkles')} Prefix output files with this string."), required = TRUE)
  multiqc$add_argument("-q", "--quiet", help = glue("{emoji('sleeping')} Shush all the logs."), action = "store_true")
  multiqc$add_argument("-f", "--format", help = glue("{emoji('icecream')} Format of output (default: %(default)s)."), default = "tsv", choices = c("tsv", "parquet", "both"))
}

multiqc_parse_args <- function(args) {
  assert_that(
    file.exists(args$json), is.character(args$prefix), nchar(args$prefix) > 0,
    is.character(args$outdir), nchar(args$outdir) > 0
  )

  prefix <- args$prefix
  json <- normalizePath(args$json)
  mkdir(args$outdir)
  outdir <- normalizePath(args$outdir)
  quiet <- args$quiet
  out_format <- args$format


  if (quiet) {
    suppressMessages(dracarys_multiqc(json, prefix, outdir, out_format))
  } else {
    dracarys_multiqc(json, prefix, outdir, out_format)
  }
}
