% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fs_local.R
\name{local_list_files_filter_relevant}
\alias{local_list_files_filter_relevant}
\title{List Relevant Files In Local Directory}
\usage{
local_list_files_filter_relevant(
  localdir,
  regexes = DR_FILE_REGEX,
  max_files = NULL
)
}
\arguments{
\item{localdir}{Path to local directory.}

\item{regexes}{Tibble with \code{regex} and \code{fun}ction name (see example).}

\item{max_files}{Max files returned.}
}
\value{
A tibble with file type, basename, size, last modified timestamp, and
path.
}
\description{
Lists relevant files in a local directory.
}
\examples{
localdir <- system.file("extdata/tso", package = "dracarys")
regexes <- tibble::tibble(regex = "multiqc_data\\\\.json$", fun = "MultiqcFile")
x <- local_list_files_filter_relevant(localdir, regexes)
}
