#' TsoAlignCollapseFusionCallerMetricsFile R6 Class
#'
#' @description
#' Contains methods for reading and displaying contents of the
#' `AlignCollapseFusionCaller_metrics.json.gz` file output from TSO.
#'
#' @examples
#' x <- system.file("extdata/tso/sample705.AlignCollapseFusionCaller_metrics.json.gz",
#'   package = "dracarys"
#' )
#' m <- TsoAlignCollapseFusionCallerMetricsFile$new(x)
#' d_parsed <- m$read() # or read(m)
#' m$write(d_parsed, out_dir = tempdir(), prefix = "sample705", out_format = c("tsv", "rds"))
#' @export
TsoAlignCollapseFusionCallerMetricsFile <- R6::R6Class(
  "TsoAlignCollapseFusionCallerMetricsFile",
  inherit = File,
  public = list(
    #' @description
    #' Reads the `AlignCollapseFusionCaller_metrics.json.gz` file output from TSO.
    #'
    #' @return tibble with the following columns:
    #' * section: name of original JSON element
    #' * name: name of metric
    #' * value: value of metric
    #' * percent: percentage
    read = function() {
      x <- self$path
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
      d <- j |>
        purrr::map(l2tib)

      # Pivot all metrics for easier ingestion,
      # and utilise the multiqc parser to rename dirty columns.
      # Keeping each list section separate for flexibility.
      secs <- c("MappingAligningPerRg", "MappingAligningSummary", "TrimmerStatistics", "CoverageSummary")
      for (sec in secs) {
        if (sec %in% names(d)) {
          d[[sec]] <- d[[sec]] |>
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
        }
      }
      if ("UmiStatistics" %in% names(d)) {
        # handle non-hist data
        d[["UmiStatisticsMain"]] <- d[["UmiStatistics"]] |>
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
        d[["UmiStatisticsHist"]] <- d[["UmiStatistics"]] |>
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
        d[["SvSummary"]] <- d[["SvSummary"]] |>
          dplyr::mutate(
            name = sub("Number of (.*) \\(PASS\\)", "\\1", .data$name),
            name = sub("breakend pairs", "bnd_pairs", .data$name),
            value = as.numeric(.data$value)
          ) |>
          tidyr::pivot_wider(names_from = "name", values_from = "value")
      }
      if ("RunTime" %in% names(d)) {
        # just keep the 'percent' column (number of seconds)
        d[["RunTime"]] <- d[["RunTime"]] |>
          dplyr::mutate(
            seconds = as.numeric(.data$percent),
            name = tools::toTitleCase(sub("Time ", "", .data$name))
          ) |>
          dplyr::select("name", "seconds") |>
          tidyr::pivot_wider(names_from = "name", values_from = "seconds")
      }
      # keep as list
      d
    },

    #' @description
    #' Writes a tidy version of the `AlignCollapseFusionCaller_metrics.json.gz` file
    #' output from TSO.
    #'
    #' Histo is the majority from UmiStatistics section, write out separately.
    #' Histo of num supporting fragments: Num of families with 0/1/2/3... raw reads.
    #' Histo of unique UMIs per fragment pos: Num of pos with 0/1/2/3... UMI seqs.
    #'
    #' @param d Parsed object from `self$read()`.
    #' @param prefix Prefix of output file(s).
    #' @param out_dir Output directory.
    #' @param out_format Format of output file(s) (one of 'tsv' (def.),
    #' 'parquet', 'both').
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
      # return list to be consistent with the read function
      # since we're writing multiple outputs, return all of them here
      invisible(d_write)
    },
    #' @description
    #' Generates the UmiStatistics Histogram plots from the
    #' `AlignCollapseFusionCaller_metrics.json.gz` file output from TSO.
    #'
    #' - Histo is the majority from UmiStatistics section, deal with it separately.
    #' - Histo of num supporting fragments: Num of families with 0/1/2/3... raw reads.
    #' - Histo of unique UMIs per fragment pos: Num of pos with 0/1/2/3... UMI seqs.
    #' @param d Parsed object from `self$read()`.
    #' @param max_num Maximum number to display in both plots.
    #' @return Both histogram plot objects.
    plot = function(d, max_num = 15) {
      h <- d[["UmiStatisticsHist"]]
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
  )
)
