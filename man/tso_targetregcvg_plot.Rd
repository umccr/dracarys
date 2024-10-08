% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tso.R
\name{tso_targetregcvg_plot}
\alias{tso_targetregcvg_plot}
\title{Plot TargetRegionCoverage}
\usage{
tso_targetregcvg_plot(d, min_pct = 2)
}
\arguments{
\item{d}{Parsed tibble.}

\item{min_pct}{Minimum percentage to be plotted (def: 2).}
}
\value{
A ggplot2 plot containing read depth on X axis and percentage
covered on Y axis.
}
\description{
Plots stuff from the \code{TargetRegionCoverage.json.gz} file output from the
TSO500 workflow.
}
