tidy_add_args <- function(subp) {
  tidy <- subp$add_parser("tidy", help="DRAGEN Output Tidying")
  tidy$add_argument("-j", "--json", help = "Path to 'multiqc_data.json'.", required = TRUE)
  tidy$add_argument("-o", "--outdir", help = "Output directory for results.", required = TRUE)
  tidy$add_argument("-p", "--prefix", help = "Prefix name for output files.", required = TRUE)
  tidy$add_argument("-q", "--quiet", help = "Suppress log printing during rendering.", action = "store_true")
}

tidy_parse_args <- function(args) {
  print("FOO")
  cli::cli_h1("{date_log()} Start tidying {args$json}!")
  stopifnot(
    file.exists(args$json), is.character(args$prefix), nchar(args$prefix) > 0,
    is.character(args$outdir), nchar(args$outdir) > 0
  )

  prefix <- args$prefix
  json <- args$json
  mkdir(args$outdir)
  outdir <- normalizePath(args$outdir)

  # main function
  d1 <- dracarys::multiqc_tidy_json(json)
  ## d2 <- select_column_subset_alignmentqc(d1)
  tsv_out <- file.path(outdir, paste0(prefix, ".tsv"))
  parquet_out <- file.path(outdir, paste0(prefix, ".parquet"))
  write_tsv(d1, tsv_out)
  write_parquet(d1, parquet_out)

  cli::cli_h1("{date_log()} End tidying {args$json}")
  cli::cli_alert_info("{date_log()} Path to output directory:\n{outdir}")
}
