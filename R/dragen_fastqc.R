#' DRAGEN FASTQC Metrics
#'
#' Read DRAGEN `fastqc_metrics.csv` file.
#'
#' @param x Path to file.
#'
#' @return Tibble with metrics.
#' @export
dragen_fastqc_metrics_read <- function(x) {
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

  # 1 POSITIONAL BASE CONTENT
  # 2 POSITIONAL BASE MEAN QUALITY
  # 3 POSITIONAL QUALITY
  # 4 READ GC CONTENT
  # 5 READ GC CONTENT QUALITY
  # 6 READ LENGTHS
  # 7 READ MEAN QUALITY
  # 8 SEQUENCE POSITIONS

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
  ) |>
    tibble::enframe(name = "fastqc_name", value = "fastqc_value")
}
