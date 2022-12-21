tso_add_args <- function(subp) {
  tso <- subp$add_parser("tso", help = "Tidy TSO500 ctDNA Output")
  tso$add_argument("-i", "--indir", help = glue("{emoji('poop')} Path to directory with TSO ctDNA workflow results."), required = TRUE)
  tso$add_argument("-o", "--outdir", help = glue("{emoji('gift')} Output tidy results to this directory."), required = TRUE)
  tso$add_argument("-q", "--quiet", help = glue("{emoji('sleeping')} Shush all the logs."), action = "store_true")
  tso$add_argument("-n", "--dryrun", help = glue("{emoji('camel')} Dry run."), action = "store_true")
  tso$add_argument("-f", "--format", help = glue("{emoji('icecream')} Format of output (default: %(default)s)."), default = "tsv", choices = c("tsv", "parquet", "rds"))
}

tso_parse_args <- function(args) {
  dir_create(args$outdir)
  outdir <- normalizePath(args$outdir)
  quiet <- args$quiet
  out_format <- args$format


  if (quiet) {
    suppressMessages(dracarys_tso(indir = indir, outdir = outdir, dryrun = quiet, token = token))
  } else {
    dracarys_tso(indir = indir, outdir = outdir, dryrun = dryrun, token = token)
  }
}
