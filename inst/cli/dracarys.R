#!/usr/bin/env Rscript

suppressPackageStartupMessages(require(argparse, include.only = "ArgumentParser"))
suppressPackageStartupMessages(require(arrow, include.only = "write_parquet"))
suppressPackageStartupMessages(require(assertthat, include.only = "assert_that"))
suppressPackageStartupMessages(require(cli))
suppressPackageStartupMessages(require(dracarys))
suppressPackageStartupMessages(require(emojifont, include.only = "emoji"))
suppressPackageStartupMessages(require(fs, include.only = c("dir_create")))
suppressPackageStartupMessages(require(glue, include.only = "glue"))
suppressPackageStartupMessages(require(readr, include.only = "write_tsv"))

pkg <- "dracarys"
prog_nm <- paste0(pkg, ".R")
version <- as.character(packageVersion(pkg))
p <- ArgumentParser(description = glue("{emoji('dragon')} DRAGEN Output Post-Processing {emoji('fire')}"), prog = prog_nm)
p$add_argument("-v", "--version", action = "version", version = glue("{prog_nm} {version}"))
subparser_name <- "subparser_name"
subp <- p$add_subparsers(help = "sub-command help", dest = subparser_name)


source(system.file("cli/multiqc.R", package = pkg))
source(system.file("cli/tso.R", package = pkg))

multiqc_add_args(subp)
tso_add_args(subp)

args <- p$parse_args()
if (length(args$subparser_name) == 0) {
  p$print_help()
} else if (args$subparser_name == "multiqc") {
  multiqc_parse_args(args)
} else if (args$subparser_name == "tso") {
  tso_parse_args(args)
} else {
  cli_alert_danger("Need to specify 'multiqc' in the cli...")
}
