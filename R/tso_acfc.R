#' TsoAlignCollapseFusionCallerMetricsFile R6 Class
#'
#' Contains methods for reading and displaying contents of the
#' `AlignCollapseFusionCaller_metrics.json.gz` file output from TSO.
#'
#' @examples
#' x <- system.file("extdata/tso/sample705.AlignCollapseFusionCaller_metrics.json.gz",
#'   package = "dracarys"
#' )
#' d <- tso_acfc_read(x)
#' p <- tso_acfc_plot(d)
#' @export
tso_acfc_read <- function(x) {
  l2tib <- function(el) {
    # list to tibble and turn cols to char
    fun1 <- function(l) {
      # handle silly NULLs..
      l[["value"]] <- l[["value"]] %||% NA
      l[["value"]] <- ifelse(l[["value"]] == "NA", NA, l[["value"]])
      tibble::as_tibble(l) |>
        dplyr::mutate(dplyr::across(dplyr::everything(), ~ as.character(.)))
    }
    el |>
      purrr::map(fun1) |>
      dplyr::bind_rows()
  }
  j <- read_jsongz_rjsonio(x, simplify = FALSE)
  secs <- list(
    s1 = c("MappingAligningPerRg", "MappingAligningSummary", "TrimmerStatistics", "CoverageSummary"),
    s2 = c("UmiStatistics", "SvSummary", "RunTime")
  )
  s1_new <- c(
    "MappingAligningPerRg" = "acfc_maprg",
    "MappingAligningSummary" = "acfc_map",
    "TrimmerStatistics" = "acfc_trim",
    "CoverageSummary" = "acfc_cvg"
  )
  secs2 <- unlist(secs, use.names = FALSE)
  secs_in_list <- secs2[secs2 %in% names(j)]
  # just extract the following elements if they exist
  j <- j[secs_in_list]
  d <- j |>
    purrr::map(l2tib)

  # Pivot all metrics for easier ingestion,
  # and utilise the multiqc parser to rename dirty columns.
  # Keeping each list section separate for flexibility.
  for (sec in secs[["s1"]]) {
    if (sec %in% names(d)) {
      new_nm <- s1_new[sec]
      d[[new_nm]] <- d[[sec]] |>
        tidyr::pivot_longer(cols = c("value", "percent"), names_to = "name1", values_to = "value1") |>
        dplyr::filter(!is.na(.data$value1)) |>
        dplyr::mutate(
          value = as.numeric(.data$value1),
          name = dplyr::if_else(.data$name1 == "percent", paste0(.data$name, " pct"), .data$name)
        ) |>
        dplyr::select("name", "value") |>
        tidyr::pivot_wider(names_from = "name", values_from = "value") |>
        dplyr::mutate(umccr_workflow = "dragen_ctdna") |>
        multiqc_rename_cols() |>
        dplyr::select(-"umccr_workflow")
      d[[sec]] <- NULL
    }
  }
  if ("UmiStatistics" %in% names(d)) {
    # handle non-hist data
    d[["acfc_umistats"]] <- d[["UmiStatistics"]] |>
      dplyr::filter(!grepl("Hist", .data$name)) |>
      tidyr::pivot_longer(cols = c("value", "percent"), names_to = "name1", values_to = "value1") |>
      dplyr::filter(!is.na(.data$value1)) |>
      dplyr::mutate(
        value = as.numeric(.data$value1),
        name = dplyr::if_else(.data$name1 == "percent", paste0(.data$name, " pct"), .data$name)
      ) |>
      dplyr::select("name", "value") |>
      tidyr::pivot_wider(names_from = "name", values_from = "value") |>
      dplyr::mutate(umccr_workflow = "dragen_ctdna") |>
      multiqc_rename_cols() |>
      dplyr::select(-"umccr_workflow")
    # handle hist data
    d[["acfc_umistatshist"]] <- d[["UmiStatistics"]] |>
      dplyr::filter(grepl("Hist", .data$name)) |>
      dplyr::mutate(
        name = sub("Histogram of ", "", .data$name),
        name = gsub(" ", "_", .data$name),
        value = as.numeric(.data$value)
      ) |>
      dplyr::group_by(name) |>
      dplyr::mutate(num = dplyr::row_number()) |>
      dplyr::ungroup() |>
      dplyr::select(c("name", "num", "value"))
    # now delete
    d[["UmiStatistics"]] <- NULL
  }
  if ("SvSummary" %in% names(d)) {
    d[["acfc_svsum"]] <- d[["SvSummary"]] |>
      dplyr::mutate(
        name = sub("Number of (.*) \\(PASS\\)", "\\1", .data$name),
        name = sub("breakend pairs", "bnd_pairs", .data$name),
        value = as.numeric(.data$value)
      ) |>
      tidyr::pivot_wider(names_from = "name", values_from = "value")
    d[["SvSummary"]] <- NULL
  }
  if ("RunTime" %in% names(d)) {
    # just keep the 'percent' column (number of seconds)
    d[["acfc_runtime"]] <- d[["RunTime"]] |>
      dplyr::mutate(
        seconds = as.numeric(.data$percent),
        name = tools::toTitleCase(sub("Time ", "", .data$name))
      ) |>
      dplyr::select("name", "seconds") |>
      tidyr::pivot_wider(names_from = "name", values_from = "seconds")
    d[["RunTime"]] <- NULL
  }
  tibble::enframe(d, name = "name", value = "data")
}

#' Plot ACFC
#'
#' Generates the UmiStatistics Histogram plots from the
#' `AlignCollapseFusionCaller_metrics.json.gz` file output from TSO.
#'
#' - Histo is the majority from UmiStatistics section, deal with it separately.
#' - Histo of num supporting fragments: Num of families with 0/1/2/3... raw reads.
#' - Histo of unique UMIs per fragment pos: Num of pos with 0/1/2/3... UMI seqs.
#'
#' @param d Parsed object from `tso_acfc_read`.
#' @param max_num Maximum number to display in both plots.
#'
#' @return Both histogram plot objects.
tso_acfc_plot <- function(d, max_num = 15) {
  h <- d |>
    dplyr::filter(.data$name == "acfc_umistatshist") |>
    dplyr::select("data") |>
    tidyr::unnest("data")
  if (nrow(h) == 0) {
    return(
      list(
        p_num_supporting_fragments = NULL,
        p_unique_umis_per_frag_pos = NULL
      )
    )
  }
  # 15 seems like a good cutoff for both plots
  p1 <- h |>
    dplyr::filter(
      .data$name == "num_supporting_fragments",
      .data$num <= max_num
    ) |>
    ggplot2::ggplot(ggplot2::aes(x = num, y = value)) +
    ggplot2::geom_line() +
    ggplot2::scale_y_continuous(
      breaks = scales::pretty_breaks(n = 10),
      labels = scales::label_number(scale_cut = scales::cut_short_scale())
    ) +
    ggplot2::theme_bw() +
    ggplot2::labs(title = "Number of families with 0/1/2/3... raw reads.") +
    ggplot2::xlab("Families") +
    ggplot2::ylab("Reads")
  p2 <- h |>
    dplyr::filter(
      name == "unique_UMIs_per_fragment_position",
      .data$num <= max_num
    ) |>
    ggplot2::ggplot(ggplot2::aes(x = num, y = value)) +
    ggplot2::geom_line() +
    ggplot2::scale_y_continuous(
      breaks = scales::pretty_breaks(n = 10),
      labels = scales::label_number(scale_cut = scales::cut_short_scale())
    ) +
    ggplot2::theme_bw() +
    ggplot2::labs(title = "Number of positions with 0/1/2/3... UMI sequences.") +
    ggplot2::xlab("Positions") +
    ggplot2::ylab("UMI Sequences")
  list(
    p_num_supporting_fragments = p1,
    p_unique_umis_per_frag_pos = p2
  )
}
