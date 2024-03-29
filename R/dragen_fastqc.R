#' FastqcMetrics R6 Class
#'
#' @description
#' Contains methods for reading and displaying contents of
#' the `fastqc_metrics.csv` file output from DRAGEN.
#'
#' @examples
#' x <- system.file("extdata/wgs/SEQC-II.fastqc_metrics.csv.gz", package = "dracarys")
#' fq <- FastqcMetricsFile$new(x)
#' d <- fq$read()
#' fq$write(d, out_dir = tempdir(), prefix = "seqc_fq", out_format = "tsv")
#' # fq$plot(d)
#'
#' @export
FastqcMetricsFile <- R6::R6Class(
  "FastqcMetricsFile",
  inherit = File,
  public = list(
    #' @description
    #' Reads the `fastqc_metrics.csv` file output from DRAGEN.
    #'
    #' @return tibble. TODO.
    read = function() {
      x <- self$path
      res <- read_fastqc_metrics(x)
    },
    #' @description
    #' Writes a tidy version of the `fastqc_metrics.csv` file output
    #' from DRAGEN.
    #'
    #' @param d Parsed object from `self$read()`.
    #' @param prefix Prefix of output file(s).
    #' @param out_dir Output directory.
    #' @param out_format Format of output file(s).
    #' @param drid dracarys ID to use for the dataset (e.g. `wfrid.123`, `prid.456`).
    write = function(d, out_dir = NULL, prefix, out_format = "tsv", drid = NULL) {
      if (!is.null(out_dir)) {
        prefix <- file.path(out_dir, prefix)
      }
      d_write <- d |>
        tibble::enframe(name = "section") |>
        dplyr::rowwise() |>
        dplyr::mutate(
          section_low = tolower(.data$section),
          p = glue("{prefix}_{.data$section_low}"),
          out = list(write_dracarys(obj = .data$value, prefix = .data$p, out_format = out_format, drid = drid))
        ) |>
        dplyr::ungroup() |>
        dplyr::select("section", "value") |>
        tibble::deframe()
      invisible(d_write)
    }
  )
)

read_fastqc_metrics <- function(x) {
  # 'SEQUENCE POSITIONS' has an extra field for 'Total Sequence Starts' which
  # can be filtered out since that can be computed by the rest of that section.
  raw <- readr::read_lines(x) |>
    tibble::as_tibble_col(column_name = "raw") |>
    dplyr::filter(!grepl("Total Sequence Starts", .data$raw))
  d <- raw |>
    tidyr::separate_wider_delim(
      "raw",
      delim = ",",
      names = c("section", "mate", "metric", "value"),
    ) |>
    dplyr::mutate(
      value = dplyr::na_if(.data$value, "NA"),
      value = as.numeric(.data$value)
    )

  # there are binned pos e.g. "149-150"
  # the value is an accumulation, so divide by number of grains
  pos_base_cont <- d |>
    dplyr::filter(.data$section == "POSITIONAL BASE CONTENT") |>
    tidyr::separate_wider_delim("metric", delim = " ", names = c("readpos", "pos", "base", "bases")) |>
    dplyr::mutate(
      is_binned = grepl("-", .data$pos), # for debugging
      bin_group = dplyr::row_number()
    ) |>
    tidyr::separate_longer_delim("pos", delim = "-") |>
    dplyr::mutate(
      value = .data$value / dplyr::n(),
      .by = "bin_group"
    ) |>
    dplyr::group_by(.data$mate, .data$pos) |>
    dplyr::mutate(
      pos = as.integer(.data$pos),
      tot = sum(.data$value),
      prop = ifelse(.data$tot == 0, 0, round(.data$value / .data$tot, 3))
    ) |>
    dplyr::ungroup() |>
    dplyr::select("mate", "pos", "base", "prop")

  # there are binned pos e.g. "149-150"
  # the value is an avg, so keep same between grains
  pos_base_mean_qual <- d |>
    dplyr::filter(.data$section == "POSITIONAL BASE MEAN QUALITY") |>
    tidyr::separate_wider_delim("metric", delim = " ", names = c("readpos", "pos", "base", "avg", "qual")) |>
    tidyr::separate_longer_delim("pos", delim = "-") |>
    dplyr::mutate(pos = as.integer(.data$pos), value = round(.data$value, 2)) |>
    dplyr::select("mate", "pos", "base", "value")

  # there are binned pos e.g. "149-150"
  # keep same value between grains
  pos_qual <- d |>
    dplyr::filter(.data$section == "POSITIONAL QUALITY") |>
    tidyr::separate_wider_delim("metric", delim = " ", names = c("readpos", "pos", "pct", "quant", "qv")) |>
    tidyr::separate_longer_delim("pos", delim = "-") |>
    dplyr::mutate(pct = as.integer(sub("%", "", .data$pct)), pos = as.integer(.data$pos)) |>
    dplyr::select("mate", "pos", "pct", "value")

  gc_cont <- d |>
    dplyr::filter(.data$section == "READ GC CONTENT") |>
    tidyr::separate_wider_delim("metric", delim = " ", names = c("pct", "GC", "reads")) |>
    dplyr::mutate(pct = as.integer(sub("%", "", .data$pct))) |>
    dplyr::select("mate", "pct", "value")

  gc_cont_qual <- d |>
    dplyr::filter(.data$section == "READ GC CONTENT QUALITY") |>
    tidyr::separate_wider_delim("metric", delim = " ", names = c("pct", "GC", "reads", "avg", "qual")) |>
    dplyr::mutate(pct = as.integer(sub("%", "", .data$pct)), value = round(.data$value, 2)) |>
    dplyr::select("mate", "pct", "value")

  # there are binned bp e.g. "149-150"
  # the value is an accumulation, so divide by number of grains
  read_len <- d |>
    dplyr::filter(.data$section == "READ LENGTHS") |>
    tidyr::separate_wider_delim("metric", delim = " ", names = c("bp", "len", "reads")) |>
    dplyr::mutate(
      bp = sub("bp", "", .data$bp),
      is_binned = grepl("-", .data$bp), # for debugging
      bin_group = dplyr::row_number()
    ) |>
    tidyr::separate_longer_delim("bp", delim = "-") |>
    dplyr::mutate(
      value = .data$value / dplyr::n(),
      .by = "bin_group"
    ) |>
    dplyr::mutate(bp = as.integer(.data$bp)) |>
    dplyr::select("mate", "bp", "value")

  read_mean_qual <- d |>
    dplyr::filter(.data$section == "READ MEAN QUALITY") |>
    tidyr::separate_wider_delim("metric", delim = " ", names = c("q", "reads")) |>
    dplyr::mutate(q = as.integer(sub("Q", "", .data$q))) |>
    dplyr::select("mate", "q", "value")

  # there are binned bp e.g. "129-130"
  # keep same value between grains
  seq_pos <- d |>
    dplyr::filter(.data$section == "SEQUENCE POSITIONS") |>
    tidyr::separate_wider_delim("metric", delim = " ", names = c("seq", "bp", "starts")) |>
    dplyr::mutate(bp = sub("bp", "", .data$bp)) |>
    tidyr::separate_longer_delim("bp", delim = "-") |>
    dplyr::mutate(bp = as.integer(.data$bp)) |>
    dplyr::select("mate", "seq", "bp", "value")

  list(
    positional_base_content = pos_base_cont,
    positional_base_mean_quality = pos_base_mean_qual,
    positional_quality = pos_qual,
    read_gc_content = gc_cont,
    read_gc_content_quality = gc_cont_qual,
    read_lengths = read_len,
    read_mean_quality = read_mean_qual,
    sequence_positions = seq_pos
  )
}
