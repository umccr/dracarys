tidy_add_args <- function(subp) {
  tidy <- subp$add_parser("tidy", help = "Tidy UMCCR Workflow Outputs")
  tidy$add_argument("-i", "--in_dir", help = glue("{emoji('snowman')} Directory with untidy UMCCR workflow results. Can be S3 or local."), required = TRUE)
  tidy$add_argument("-o", "--out_dir", help = glue("{emoji('fire')} Directory to output tidy results."), required = TRUE)
  tidy$add_argument("-p", "--prefix", help = glue("{emoji('violin')} Prefix string used for all results."), required = TRUE)
  tidy$add_argument("-l", "--local_dir", help = glue("{emoji('inbox_tray')} If input is an S3 directory, download the recognisable files to this directory. Default: '<out_dir>/dracarys_s3_sync'."))
  tidy$add_argument("-f", "--format", help = glue("{emoji('art')} Format of output. Default: %(default)s."), default = "tsv")
  tidy$add_argument("-n", "--dryrun", help = glue("{emoji('camel')} Dry run - just show files to be tidied."), action = "store_true")
  tidy$add_argument("-q", "--quiet", help = glue("{emoji('sleeping')} Shush all the logs."), action = "store_true")
}

tidy_parse_args <- function(args) {
  # handle dir creation
  fs::dir_create(c(args$out_dir, args$local_dir))
  out_dir <- normalizePath(args$out_dir)
  local_dir <- args$local_dir
  if (!is.null(local_dir)) {
    local_dir <- normalizePath(local_dir)
  }
  if (!is.null(args$rds_dir)) {
    rds_dir <- normalizePath(args$rds_dir)
    rds_path <- file.path(rds_dir, glue("{args$prefix}_dracarys_data.rds"))
  }

  tidy_args <- list(
    in_dir = args$in_dir,
    out_dir = out_dir,
    prefix = args$prefix,
    local_dir = local_dir,
    out_format = args$format,
    dryrun = args$dryrun,
    pattern = args$pattern
  )

  # tidy run
  if (args$quiet) {
    rds_obj <- suppressMessages(do.call(umccr_tidy, tidy_args))
  } else {
    rds_obj <- do.call(umccr_tidy, tidy_args)
  }
}
