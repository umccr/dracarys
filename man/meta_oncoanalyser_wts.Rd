% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/portal_meta.R
\name{meta_oncoanalyser_wts}
\alias{meta_oncoanalyser_wts}
\title{Metadata for oncoanalyser_wts workflow}
\usage{
meta_oncoanalyser_wts(pmeta, status = "Succeeded")
}
\arguments{
\item{pmeta}{Path to portal workflows metadata table, or tibble with already parsed data.}

\item{status}{Workflow status to keep (default: Succeeded).}
}
\value{
A tibble with metadata per workflow run.
}
\description{
Metadata for oncoanalyser_wts workflow
}
\examples{
pmeta <- system.file("extdata/portal_meta_top4.csv", package = "dracarys")
(m <- meta_oncoanalyser_wts(pmeta))
}
