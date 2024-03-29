% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dragen.R
\name{TimeMetricsFile}
\alias{TimeMetricsFile}
\title{TimeMetricsFile R6 Class}
\description{
Contains methods for reading contents of
the \code{time_metrics.csv} file output from DRAGEN, which contains
a breakdown of the run duration for each DRAGEN process.
}
\examples{
x <- system.file("extdata/wgs/SEQC-II.time_metrics.csv.gz", package = "dracarys")
tm <- TimeMetricsFile$new(x)
d <- tm$read() # or read(tm)
tm$write(d, out_dir = tempdir(), prefix = "seqc_time", out_format = "tsv")
}
\section{Super class}{
\code{\link[dracarys:File]{dracarys::File}} -> \code{TimeMetricsFile}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-TimeMetricsFile-read}{\code{TimeMetricsFile$read()}}
\item \href{#method-TimeMetricsFile-write}{\code{TimeMetricsFile$write()}}
\item \href{#method-TimeMetricsFile-clone}{\code{TimeMetricsFile$clone()}}
}
}
\if{html}{\out{
<details open><summary>Inherited methods</summary>
<ul>
<li><span class="pkg-link" data-pkg="dracarys" data-topic="File" data-id="bname"><a href='../../dracarys/html/File.html#method-File-bname'><code>dracarys::File$bname()</code></a></span></li>
<li><span class="pkg-link" data-pkg="dracarys" data-topic="File" data-id="initialize"><a href='../../dracarys/html/File.html#method-File-initialize'><code>dracarys::File$initialize()</code></a></span></li>
<li><span class="pkg-link" data-pkg="dracarys" data-topic="File" data-id="print"><a href='../../dracarys/html/File.html#method-File-print'><code>dracarys::File$print()</code></a></span></li>
<li><span class="pkg-link" data-pkg="dracarys" data-topic="File" data-id="type"><a href='../../dracarys/html/File.html#method-File-type'><code>dracarys::File$type()</code></a></span></li>
</ul>
</details>
}}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-TimeMetricsFile-read"></a>}}
\if{latex}{\out{\hypertarget{method-TimeMetricsFile-read}{}}}
\subsection{Method \code{read()}}{
Reads the \code{time_metrics.csv} file.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{TimeMetricsFile$read()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
tibble with one row and metrics spread across individual columns.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-TimeMetricsFile-write"></a>}}
\if{latex}{\out{\hypertarget{method-TimeMetricsFile-write}{}}}
\subsection{Method \code{write()}}{
Writes a tidy version of the \code{time_metrics.csv} file output
from DRAGEN.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{TimeMetricsFile$write(
  d,
  out_dir = NULL,
  prefix,
  out_format = "tsv",
  drid = NULL
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{d}}{Parsed object from \code{self$read()}.}

\item{\code{out_dir}}{Output directory.}

\item{\code{prefix}}{Prefix of output file(s).}

\item{\code{out_format}}{Format of output file(s).}

\item{\code{drid}}{dracarys ID to use for the dataset (e.g. \code{wfrid.123}, \code{prid.456}).}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-TimeMetricsFile-clone"></a>}}
\if{latex}{\out{\hypertarget{method-TimeMetricsFile-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{TimeMetricsFile$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
