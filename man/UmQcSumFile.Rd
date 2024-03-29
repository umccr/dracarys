% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/umccrise.R
\name{UmQcSumFile}
\alias{UmQcSumFile}
\title{UmQcSumFile R6 Class}
\description{
Contains methods for reading and displaying contents of the
\code{qc_summary.tsv.gz} file with QC summary metrics output from umccrise.
}
\examples{
\dontrun{
x <- "/path/to/snv_2015.tsv.gz"
d <- UmQcSumFile$new(x)
d_parsed <- d$read() # or read(d)
d$write(d_parsed, out_dir = tempdir(), prefix = "sample705", out_format = "tsv")
}
}
\section{Super class}{
\code{\link[dracarys:File]{dracarys::File}} -> \code{UmQcSumFile}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-UmQcSumFile-read}{\code{UmQcSumFile$read()}}
\item \href{#method-UmQcSumFile-write}{\code{UmQcSumFile$write()}}
\item \href{#method-UmQcSumFile-clone}{\code{UmQcSumFile$clone()}}
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
\if{html}{\out{<a id="method-UmQcSumFile-read"></a>}}
\if{latex}{\out{\hypertarget{method-UmQcSumFile-read}{}}}
\subsection{Method \code{read()}}{
Reads the \code{qc_summary.tsv.gz} file output from umccrise.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{UmQcSumFile$read()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
A tibble.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-UmQcSumFile-write"></a>}}
\if{latex}{\out{\hypertarget{method-UmQcSumFile-write}{}}}
\subsection{Method \code{write()}}{
Writes a tidy version of the \code{qc_summary.tsv.gz} QC summary file output
from umccrise.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{UmQcSumFile$write(d, out_dir, prefix, out_format = "tsv", drid = NULL)}\if{html}{\out{</div>}}
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
\if{html}{\out{<a id="method-UmQcSumFile-clone"></a>}}
\if{latex}{\out{\hypertarget{method-UmQcSumFile-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{UmQcSumFile$clone(deep = FALSE)}\if{html}{\out{</div>}}
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
