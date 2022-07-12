#!/usr/bin/env Rscript

suppressPackageStartupMessages(require(argparse, include.only = "ArgumentParser"))
suppressPackageStartupMessages(require(cli))
suppressPackageStartupMessages(require(dracarys))
suppressPackageStartupMessages(require(glue, include.only = "glue"))
suppressPackageStartupMessages(require(arrow, include.only = "write_parquet"))
suppressPackageStartupMessages(require(readr, include.only = "write_tsv"))

pkg <- "dracarys"
prog_nm <- paste0(pkg, ".R")
version <- as.character(packageVersion(pkg))
p <- ArgumentParser(description = "DRAGEN Output Tidying", prog = prog_nm)
p$add_argument("-v", "--version", action = "version", version = glue("{prog_nm} {version}"))
subparser_name <- "subparser_name"
subp <- p$add_subparsers(help = "sub-command help", dest = subparser_name)


source(system.file("cli/tidy.R", package = pkg))

tidy_add_args(subp)

args <- p$parse_args()
if (length(args$subparser_name) == 0) {
  p$print_help()
} else if (args$subparser_name == "tidy") {
  tidy_parse_args(args)
} else {
  cli_alert_danger("Need to specify 'tidy' in the cli...")
}
