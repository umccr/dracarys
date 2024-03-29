% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dragen.R
\name{MappingMetricsFile}
\alias{MappingMetricsFile}
\title{MappingMetricsFile R6 Class}
\description{
Contains methods for reading and displaying contents of
the \code{mapping_metrics.csv} file output from DRAGEN.
This file contains mapping and aligning metrics, like the metrics computed by
the Samtools Flagstat command. These metrics are available on an aggregate
level (over all input data), and on a per read group level. NOTE: we are
keeping only the read group level metrics (i.e. removing the aggregate data).
Unless explicitly stated, the metrics units are in reads (i.e., not in
terms of pairs or alignments).
}
\examples{
x <- system.file("extdata/wgs/SEQC-II.mapping_metrics.csv.gz", package = "dracarys")
mm <- MappingMetricsFile$new(x)
d <- mm$read() # or read(mm)
mm$write(d, out_dir = tempdir(), prefix = "seqc_mm", out_format = "tsv")

}
\section{Super class}{
\code{\link[dracarys:File]{dracarys::File}} -> \code{MappingMetricsFile}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-MappingMetricsFile-read}{\code{MappingMetricsFile$read()}}
\item \href{#method-MappingMetricsFile-write}{\code{MappingMetricsFile$write()}}
\item \href{#method-MappingMetricsFile-clone}{\code{MappingMetricsFile$clone()}}
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
\if{html}{\out{<a id="method-MappingMetricsFile-read"></a>}}
\if{latex}{\out{\hypertarget{method-MappingMetricsFile-read}{}}}
\subsection{Method \code{read()}}{
Reads the \code{mapping_metrics.csv} file output from DRAGEN.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{MappingMetricsFile$read()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
tibble with one row of X metrics per read group.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-MappingMetricsFile-write"></a>}}
\if{latex}{\out{\hypertarget{method-MappingMetricsFile-write}{}}}
\subsection{Method \code{write()}}{
Writes a tidy version of the \code{mapping_metrics.csv} file output
from DRAGEN.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{MappingMetricsFile$write(
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
\if{html}{\out{<a id="method-MappingMetricsFile-clone"></a>}}
\if{latex}{\out{\hypertarget{method-MappingMetricsFile-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{MappingMetricsFile$clone(deep = FALSE)}\if{html}{\out{</div>}}
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
