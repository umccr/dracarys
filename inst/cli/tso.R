tso_add_args <- function(subp) {
  tso <- subp$add_parser("tso", help = "Tidy TSO500 ctDNA Output")
  tso$add_argument("-i", "--in_dir", help = glue("{emoji('poop')} Directory with TSO500 ctDNA workflow results (GDS or local)."), required = TRUE)
  tso$add_argument("-o", "--out_dir", help = glue("{emoji('gift')} Output tidy results to this directory."), required = TRUE)
  tso$add_argument("-r", "--report_dir", help = glue("{emoji('sparkles')} Output HTML report with tidy RDS object to this directory."))
  tso$add_argument("-p", "--prefix", help = glue("{emoji('dancer')} Prefix string (used for all results)."), required = TRUE)
  tso$add_argument("-g", "--gds_local_dir", help = glue("{emoji('clipboard')} If input is a GDS directory, download the 'recognisable' files to this directory. If not specified, files will be downloaded to 'out_dir/dracarys_gds_sync'."))
  tso$add_argument("--rds_dir", help = glue("{emoji('droplet')} Directory to save RDS object with results from tidy function."))
  tso$add_argument("-f", "--format", help = glue("{emoji('icecream')} Format of output (default: %(default)s)."), default = "tsv", choices = c("tsv", "parquet", "both"))
  tso$add_argument("-n", "--dryrun", help = glue("{emoji('camel')} Dry run (just print tibble with files to be tidied)."), action = "store_true")
  tso$add_argument("-q", "--quiet", help = glue("{emoji('sleeping')} Shush all the logs."), action = "store_true")
  tso$add_argument("--quiet_rmd", help = glue("{emoji('sleeping')} Shush just the Rmd rendering logs."), action = "store_true")
}

tso_parse_args <- function(args) {
  # handle dir creation
  fs::dir_create(c(args$out_dir, args$report_dir, args$gds_local_dir, args$rds_dir))
  out_dir <- normalizePath(args$out_dir)
  gds_local_dir <- args$gds_local_dir
  rds_dir <- tempdir()
  rds_path <- tempfile(fileext = ".rds")
  rds_obj <- NULL
  report_dir <- args$report_dir
  report_run <- FALSE
  if (!is.null(report_dir)) {
    report_dir <- normalizePath(report_dir)
    report_run <- TRUE
  }
  if (!is.null(gds_local_dir)) {
    gds_local_dir <- normalizePath(gds_local_dir)
  }

  if (!is.null(args$rds_dir)) {
    rds_dir <- normalizePath(args$rds_dir)
    rds_path <- file.path(rds_dir, glue("{args$prefix}_dracarys_data.rds"))
  }

  tidy_args <- list(
    in_dir = args$in_dir,
    out_dir = out_dir,
    prefix = args$prefix,
    gds_local_dir = gds_local_dir,
    out_format = args$format,
    dryrun = args$dryrun,
    token = args$token
  )

  # tso_tidy run
  if (args$quiet) {
    rds_obj <- suppressMessages(do.call(tso_tidy, tidy_args))
  } else {
    rds_obj <- do.call(tso_tidy, tidy_args)
  }
  rmd_args <- list(
    tidy_rds = rds_path,
    out_dir = report_dir,
    prefix = args$prefix,
    quiet = args$quiet_rmd || args$quiet
  )

  # tso_rmd run
  if (report_run) {
    readr::write_rds(x = rds_obj, file = rds_path)
    if (args$quiet) {
      suppressMessages(do.call(tso_rmd, rmd_args))
    } else {
      do.call(tso_rmd, rmd_args)
    }
  }
}
