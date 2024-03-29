% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dragen.R
\name{FragmentLengthHistFile}
\alias{FragmentLengthHistFile}
\title{FragmentLengthHistFile R6 Class}
\description{
Contains methods for reading and plotting contents of
the \code{fragment_length_hist.csv} file output from DRAGEN.
}
\examples{
x <- system.file("extdata/wgs/SEQC-II.fragment_length_hist.csv.gz", package = "dracarys")
fl <- FragmentLengthHistFile$new(x)
d <- fl$read() # or read(fl)
fl$plot(d) # or plot(fl)
fl$write(d |> dplyr::filter(count > 10), out_dir = tempdir(), prefix = "seqc_fl")
}
\section{Super class}{
\code{\link[dracarys:File]{dracarys::File}} -> \code{FragmentLengthHistFile}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-FragmentLengthHistFile-read}{\code{FragmentLengthHistFile$read()}}
\item \href{#method-FragmentLengthHistFile-write}{\code{FragmentLengthHistFile$write()}}
\item \href{#method-FragmentLengthHistFile-plot}{\code{FragmentLengthHistFile$plot()}}
\item \href{#method-FragmentLengthHistFile-clone}{\code{FragmentLengthHistFile$clone()}}
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
\if{html}{\out{<a id="method-FragmentLengthHistFile-read"></a>}}
\if{latex}{\out{\hypertarget{method-FragmentLengthHistFile-read}{}}}
\subsection{Method \code{read()}}{
Reads the \code{fragment_length_hist.csv} file, which contains the
fragment length distribution for each sample.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{FragmentLengthHistFile$read()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
A tibble with the following columns:
\itemize{
\item sample: name of sample
\item fragmentLength: estimated fragment length
\item count: number of reads with estimated fragment length
}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-FragmentLengthHistFile-write"></a>}}
\if{latex}{\out{\hypertarget{method-FragmentLengthHistFile-write}{}}}
\subsection{Method \code{write()}}{
Writes a tidy version of the \code{fragment_length_hist.csv} file output
from DRAGEN.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{FragmentLengthHistFile$write(
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
\if{html}{\out{<a id="method-FragmentLengthHistFile-plot"></a>}}
\if{latex}{\out{\hypertarget{method-FragmentLengthHistFile-plot}{}}}
\subsection{Method \code{plot()}}{
Plots the fragment length distributions as given in the
\code{fragment_length_hist.csv} file.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{FragmentLengthHistFile$plot(d, min_count = 10)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{d}}{Parsed object from \code{self$read()}.}

\item{\code{min_count}}{Minimum read count to be plotted (Default: 10).}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A ggplot2 plot containing fragment lengths on X axis and read counts
on Y axis for each sample.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-FragmentLengthHistFile-clone"></a>}}
\if{latex}{\out{\hypertarget{method-FragmentLengthHistFile-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{FragmentLengthHistFile$clone(deep = FALSE)}\if{html}{\out{</div>}}
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
