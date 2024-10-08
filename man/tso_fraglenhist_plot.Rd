% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tso.R
\name{tso_fraglenhist_plot}
\alias{tso_fraglenhist_plot}
\title{Plot Fragment Length Hist}
\usage{
tso_fraglenhist_plot(d, min_count = 10)
}
\arguments{
\item{d}{Parsed tibble.}

\item{min_count}{Minimum read count to be plotted (def: 10).}
}
\value{
A ggplot2 plot containing fragment lengths on X axis and read counts
on Y axis for each sample.
}
\description{
Plots the fragment length distributions as given in the
\code{fragment_length_hist} file.
}
