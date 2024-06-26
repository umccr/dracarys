% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/s3.R
\name{dr_s3_download}
\alias{dr_s3_download}
\title{dracarys S3 Download}
\usage{
dr_s3_download(
  s3dir,
  outdir,
  page_size = 100,
  pattern = NULL,
  regexes = DR_FILE_REGEX,
  dryrun = FALSE
)
}
\arguments{
\item{s3dir}{Full path to S3 directory.}

\item{outdir}{Path to output directory.}

\item{page_size}{Page size (def: 100).}

\item{pattern}{Pattern to further filter the returned file type tibble.}

\item{regexes}{Tibble with regex and function name.}

\item{dryrun}{If TRUE, just list the files that will be downloaded (don't
download them).}
}
\description{
Download only S3 files that can be processed by dracarys.
}
\examples{
s3dir <- file.path(
  "s3://umccr-primary-data-prod/UMCCR-Validation/SBJ00596",
  "ctTSO/2021-03-17/PTC_SSqCMM05pc_L2100067"
)
outdir <- sub("s3:/", "~/s3", s3dir)

}
