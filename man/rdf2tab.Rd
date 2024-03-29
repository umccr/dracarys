% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{rdf2tab}
\alias{rdf2tab}
\title{Write Local R Dataframe to Spark-backed Table}
\usage{
rdf2tab(rdf, outpath, drid = NULL, ...)
}
\arguments{
\item{rdf}{Local R dataframe}

\item{outpath}{Path to output table.}

\item{drid}{dracarys ID to use for the dataset (e.g. \code{wfrid.123}, \code{prid.456})}

\item{...}{Additional arguments for \code{SparkR::write.df}.}
}
\description{
First converts local R dataframe to Spark DataFrame using SparkR,
then \emph{appends} it to the specified table.
}
\examples{
\dontrun{
rdf <- mtcars
rdf2tab(rdf, "dev.wf1.mtcars", drid = "wfr.123")
}
}
