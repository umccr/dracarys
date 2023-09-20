#' MultiqcFile R6 Class
#'
#' @description
#' Contains methods for reading and displaying contents of the
#' `multiqc_data.json` file output from MultiQC.
#'
#' @examples
#' \dontrun{
#' x <- "/path/to/multiqc_data.json"
#' mqc <- MultiqcFile$new(x)
#' mqc_parsed <- mqc$read() # or read(mqc)
#' mqc$write(mqc_parsed, out_dir = tempdir(), prefix = "sample705", out_format = "both")
#' }
#' @export
MultiqcFile <- R6::R6Class(
  "MultiqcFile",
  inherit = File,
  public = list(
    #' @description
    #' Reads the `multiqc_data.json` file output from MultiQC.
    #'
    #' @return A tidy tibble.
    #'   - label:
    read = function() {
      x <- self$path
      multiqc_tidy_json(x)
    },

    #' @description
    #' Writes a tidy version of the `multiqc_data.json` file output from MultiQC.
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
      # prefix2 <- glue("{prefix}multiqc")
      write_dracarys(obj = d, prefix = prefix, out_format = out_format, drid = drid)
    }
  )
)

#' Tidy MultiQC JSON
#'
#' Tidies 'multiqc_data.json' output from MultiQC.
#' Modified from the awesome <https://github.com/multimeric/TidyMultiqc>.
#' @param j Path to `multiqc_data.json` file.
#' @return A tidy tibble where each column corresponds to a single metric,
#' and each row corresponds to a single sample.
#' @export
multiqc_tidy_json <- function(j) {
  p <- RJSONIO::fromJSON(j)
  cdate <- p[["config_creation_date"]] %||% "UNKNOWN"
  cdate <- multiqc_date_fmt(cdate)
  workflow <- .multiqc_guess_workflow(p)
  d <- dracarys::multiqc_parse_gen(p)
  if (workflow %in% c("dragen_umccrise", "dragen_somatic")) {
    # replace the "NA" strings with NA, else we get a column class error
    # due to trying to bind string ('NA') with numeric.
    # https://stackoverflow.com/questions/35292757/replace-values-in-list
    d <- rapply(d, function(x) ifelse(x == "NA", NA, x), how = "replace")
  }
  if (workflow == "bcl_convert") {
    # general_stats_data is empty
    # use raw instead
    return(dracarys::multiqc_parse_raw(p))
  }
  d <- d |>
    # get rid of duplicated elements - see umccr/dracarys#96
    purrr::map(\(x) {
      x[which(duplicated(names(x)))] <- NULL
      x
    }) |>
    dplyr::bind_rows(.id = "umccr_id") |>
    dplyr::mutate(
      config_creation_date = cdate,
      umccr_workflow = workflow
    ) |>
    dplyr::select("umccr_id", "umccr_workflow", "config_creation_date", dplyr::everything())

  if (workflow == "dragen_transcriptome") {
    # discard Ref_X control samples
    wts_ref_samples <- paste0("Ref_", 1:6)
    d <- d |>
      dplyr::filter(!.data$umccr_id %in% wts_ref_samples)
  } else if (workflow %in% c("dragen_umccrise", "bcbio_umccrise")) {
    # discard Alice, Bob etc. control samples
    um_ref_samples <- c("Alice", "Bob", "Chen", "Elon", "Dakota")
    um_ref_samples <- paste0(um_ref_samples, rep(c("_T", "_B", ""), each = length(um_ref_samples)))
    d <- d |>
      dplyr::filter(!.data$umccr_id %in% um_ref_samples)
  }
  return(multiqc_rename_cols(d))
}

multiqc_rename_cols <- function(d) {
  umccr_workflows <- c(
    "dragen_alignment", "dragen_somatic",
    "dragen_transcriptome", "dragen_umccrise",
    "dragen_ctdna",
    "bcbio_umccrise", "bcbio_wgs", "bcbio_wts"
  )

  w <- unique(d[["umccr_workflow"]])
  assertthat::assert_that(length(w) == 1)
  # if known workflow rename columns
  if (w %in% umccr_workflows) {
    m <- MULTIQC_COLUMNS |>
      dplyr::filter(.data$workflow == w)
    # warn if any metrics aren't accounted for
    missing_col <- names(d)[!names(d) %in% m[["raw"]]]
    if (length(missing_col) > 0) {
      warning(
        "The following metrics in the input MultiQC JSON are unknown: ",
        paste(missing_col, collapse = ", ")
      )
    }
    # grab only those m rows of interest since not all json files contain
    # the same metrics (e.g. tool versions change).
    m <- m |>
      dplyr::filter(.data$raw %in% names(d))
    # awesomeness
    rename_vec <- purrr::set_names(m[["raw"]], m[["clean"]])
    d <- dplyr::rename(d, !!!rename_vec)
  }
  return(d)
}

.multiqc_guess_workflow <- function(p) {
  assertthat::assert_that(all(c("config_title", "report_data_sources") %in% names(p)))
  config_title <- p[["config_title"]] %||% "Unknown"
  ds <- names(p[["report_data_sources"]])
  # bcbio
  if ("bcbio" %in% ds) {
    # wgs, wts, or umccrise?
    if (all(c("Salmon", "STAR", "QualiMap") %in% ds)) {
      return("bcbio_wts")
    } else if (all(c("PURPLE", "Conpair", "mosdepth") %in% ds)) {
      return("bcbio_umccrise")
    } else if (all(c("Samtools", "Bcftools (somatic)", "Bcftools (germline)") %in% ds)) {
      return("bcbio_wgs")
    } else {
      warning(glue(
        "Unknown which bcbio workflow this MultiQC JSON was generated from",
      ))
      return("bcbio_unknown")
    }
  }

  # dragen
  if ("DRAGEN" %in% ds) {
    dragen_workflows <- c("alignment", "transcriptome", "somatic", "ctdna")
    if (all(c("PURPLE", "Conpair", "mosdepth") %in% ds)) {
      return("dragen_umccrise")
    } else if (grepl("^UMCCR MultiQC Dragen", config_title)) {
      w <- tolower(sub("UMCCR MultiQC Dragen (.*) Report for .*", "\\1", config_title))
      if (w == "somatic and germline") {
        w <- "somatic"
      }
      assertthat::assert_that(w %in% dragen_workflows)
      return(paste0("dragen_", w))
    } else if (grepl("^UMCCR MultiQC ctDNA", config_title)) {
      return("dragen_ctdna")
    } else {
      warning(glue(
        "config_title: '{config_title}'.\n",
        "Unknown which DRAGEN workflow this MultiQC JSON was generated from",
      ))
      return("dragen_unknown")
    }
  }
  if (length(ds) == 1 && ds == "bclconvert") {
    return("bcl_convert")
  }
  return("UNKNOWN")
}

#' Parse MultiQC 'report_general_stats_data' JSON Element
#'
#' Parses MultiQC 'report_general_stats_data' JSON Element. Modified from the
#' awesome <https://github.com/multimeric/TidyMultiqc>.
#' @param p Parsed MultiQC JSON.
#' @return A list.
#' @export
multiqc_parse_gen <- function(p) {
  el <- "report_general_stats_data"
  assertthat::assert_that(inherits(p, "list"), el %in% names(p))
  if (length(p[[el]]) == 0) {
    return(list(NA))
  }
  p[[el]] |> purrr::reduce(~ purrr::list_merge(.x, !!!.y))
}

#' Parse MultiQC 'report_saved_raw_data' JSON Element
#'
#' Parses MultiQC 'report_saved_raw_data' JSON Element.
#' @param p Parsed MultiQC JSON.
#' @return A list.
#' @export
multiqc_parse_raw <- function(p) {
  x <- p[["report_saved_raw_data"]]
  tool_nms <- names(x)
  res <- list()
  # some elements can be empty dicts
  replace_empty_with_na <- function(x) {
    ifelse(length(x) == 0, NA_real_, x)
  }
  for (tool in tool_nms) {
    tool_d <- x[[tool]]
    sample_nms <- names(tool_d)
    for (sample in sample_nms) {
      d <- tool_d[[sample]]
      if (is.list(d)) {
        d <- purrr::map(d, replace_empty_with_na)
      }
      res[[tool]][[sample]] <- tibble::as_tibble_row(d)
    }
    res[[tool]] <- res[[tool]] |>
      dplyr::bind_rows(.id = "multiqc_sample")
  }
  res |>
    purrr::map(\(x) tidyr::nest(x, .by = "multiqc_sample")) |>
    dplyr::bind_rows(.id = "multiqc_tool")
}

# From https://github.com/multimeric/TidyMultiqc
multiqc_kv_map <- function(l, func, map_keys = FALSE) {
  mapper <- ifelse(map_keys, purrr::imap, purrr::map)
  mapped <- mapper(l, func) |> purrr::set_names(nm = NULL)
  keys <- mapped |> purrr::map_chr("key")
  vals <- mapped |> purrr::map("value")
  vals |> purrr::set_names(keys)
}

# Tibble with three columns: workflow name, raw and cleaned multiqc metric name.
# For consistency, make sure that names are in snake_case and that the tool that
# generates them is added as the suffix e.g. awesome_metric_tool1.
MULTIQC_COLUMNS <- system.file("extdata/multiqc_column_map.tsv", package = "dracarys") |>
  readr::read_tsv(col_types = "ccc")


#' Append New MultiQC Workflow Columns
#'
#' @param d Tidy MultiQC tibble with raw names (i.e. pre-rename).
#' @param w New workflow name.
#' @param m Path to the 'inst/extdata/multiqc_column_map.tsv' dracarys file.
#'
#' @export
multiqc_column_map_append <- function(d, w, m) {
  assertthat::assert_that(
    inherits(d, "data.frame"),
    inherits(w, "character"),
    file.exists(m)
  )
  new_workflow <- w
  tibble::tibble(new_workflow = new_workflow, raw = colnames(d)) |>
    dplyr::left_join(MULTIQC_COLUMNS, by = "raw", multiple = "first") |>
    dplyr::distinct(.data$new_workflow, .data$raw, .data$clean) |>
    dplyr::rename(workflow = "new_workflow") |>
    readr::write_tsv(append = TRUE, file = m)
}

#' Format MultiQC Config Date
#'
#' @param cdate String with config date in "YYYY-MM-DD, HH:MM UTC" format.
#' @return Properly formatted datetime string in "YYYY-MM-DDTHH:MM:SS" format.
#' If the input string was formatted differently, it will return NA.
#'
#' @examples
#' cdate <- "2023-04-07, 09:09 UTC"
#' (res1 <- multiqc_date_fmt(cdate))
#' (res2 <- multiqc_date_fmt("2023-04-07"))
#' (res3 <- multiqc_date_fmt("UNKNOWN"))
#' @testexamples
#' expect_equal(res1, "2023-04-07T09:09:00")
#' expect_equal(res2, NA_character_)
#' expect_equal(res3, NA_character_)
#' @export
multiqc_date_fmt <- function(cdate) {
  # Input: "2023-04-07, 09:09 UTC"
  # "YYYY-MM-DD, HH:MM UTC"
  # Output: 2023-04-07T09:09:00
  res <- sub(" UTC", "", cdate) |>
    strptime(format = "%Y-%m-%d, %H:%M", tz = "UTC") |>
    format("%Y-%m-%dT%H:%M:%S")
  return(res)
}

#' Parse Plot Data from MultiQC JSON
#'
#' @param j Path to `multiqc_data.json` file.
#' @param plot_names Names of plots to parse. Use "everything" if you wantz all
#' the plotz.
#'
#' @return Nested tibble with plot name and result as list column (use `tidyr::unnest` to access).
#'
#' @examples
#' \dontrun{
#' j <- here::here("nogit/warehouse/cttso/2023/05/wfr.ec3748db414f422381e419d367f73eec/dracarys_gds_sync/multiqc_data.json")
#' multiqc_parse_plots(j, plot_names = c("fastqc_gc_content_mean_sequence_quality_plot", "dragen_coverage_per_contig"))
#' multiqc_parse_plots(j, plot_names = NULL)
#' }
#' @export
multiqc_parse_plots <- function(j, plot_names = NULL) {
  parsed <- RJSONIO::fromJSON(j)
  funcs <- c(
    bar_graph = "multiqc_parse_bargraph_plot",
    xy_line = "multiqc_parse_xyline_plot"
  )
  assertthat::assert_that(
    inherits(parsed, "list"),
    "report_plot_data" %in% names(parsed)
  )
  plot_data <- parsed[["report_plot_data"]]
  # I want all the plots
  if (!is.null(plot_names) & ("everything" %in% plot_names)) {
    plot_names <- names(plot_data)
  }
  # loop over each plot item and run
  # corresponding function based on plot type.
  final_cols <- c("plot_nm", "plot_res")
  res <- multiqc_list_plots(plot_data) |>
    # keep only specified plots; NULL plot_names will discard everything.
    dplyr::filter(
      .data$plot_nm %in% plot_names,
      .data$plot_type %in% names(funcs)
    )
  if (nrow(res) == 0) {
    return(empty_tbl(cnames = final_cols))
  }
  res |>
    dplyr::rowwise() |>
    dplyr::mutate(
      input = list(plot_data[[.data$plot_nm]]),
      plot_res = list(dr_func_eval(f = funcs[.data$plot_type], v = funcs)(.data$input))
    ) |>
    dplyr::ungroup() |>
    dplyr::select(dplyr::all_of(final_cols))
}

#' MultiQC Extract XY Line Plot Data
#'
#' Extracts `xy_line` data for the given plot name from a MultiQC JSON.
#'
#' Note that the `dragen_coverage_per_contig/non_main_contig` plots do not conform
#' to the normal structure of `xy_line` plots, so we handle those separately.
#'
#' @param dat Parsed JSON list element with specific plot data to extract.
#'
#' @return Tibble with name and data list-column.
#' @export
multiqc_parse_xyline_plot <- function(dat) {
  assertthat::assert_that(dat[["plot_type"]] == "xy_line")
  if (dat[["config"]][["id"]] %in% c("dragen_coverage_per_contig", "dragen_coverage_per_non_main_contig")) {
    return(multiqc_parse_xyline_plot_contig_cvg(dat))
  }
  dat[["datasets"]][[1]] |>
    purrr::map(
      function(nm_data) {
        name <- nm_data[["name"]]
        d <- nm_data[["data"]] |>
          # handle NULLs in raw JSON
          purrr::map(~ tibble::tibble(x = unlist(.x[1]) %||% NA_real_, y = unlist(.x[2]) %||% NA_real_)) |>
          dplyr::bind_rows()
        tibble::tibble(name = name, x = d[["x"]], y = d[["y"]])
      }
    ) |>
    dplyr::bind_rows() |>
    tidyr::nest(dat = c(.data$x, .data$y))
}

#' MultiQC Extract XY Line Plot Data for DRAGEN Contig Coverage
#'
#' @param dat Parsed JSON list element with specific plot data to extract.
#'
#' @return Tibble with name and data list-column.
#' @export
multiqc_parse_xyline_plot_contig_cvg <- function(dat) {
  assertthat::assert_that(dat[["plot_type"]] == "xy_line")
  assertthat::assert_that(
    dat[["config"]][["id"]] %in% c("dragen_coverage_per_contig", "dragen_coverage_per_non_main_contig")
  )
  contig <- dat[["config"]][["categories"]]
  dat[["datasets"]][[1]] |>
    purrr::map(
      function(nm_data) {
        name <- nm_data[["name"]]
        val <- nm_data[["data"]]
        assertthat::assert_that(length(val) == length(contig))
        tibble::tibble(name = name, contig = contig, cvg = val)
      }
    ) |>
    dplyr::bind_rows() |>
    tidyr::nest(dat = c(.data$contig, .data$cvg))
}

#' MultiQC Extract Bar Graph Data
#'
#' Extracts `bar_graph` data for the given plot name from a MultiQC JSON.
#' Each `bar_graph` element in the JSON object contains:
#'
#' - a `samples` array of N subarrays (of strings) of length A.
#' - a `datasets` array of N subarrays (of objects) of length B.
#' Each of the subarrays contains name (single string) and data (array of length
#' equal to the length of the index-corresponding `samples` subarray).
#'
#' @param dat Parsed JSON list element with specific plot data to extract.
#'
#' @return Tibble with name and data list-column.
#' @examples
#' \dontrun{
#' j1 <- here::here("nogit/bcl_convert/multiqc_data.json")
#' j2 <- here::here("nogit/warehouse/wgs_tumor_normal/SBJ03197/2023-04-30_0813ce/dracarys_gds_sync/multiqc_data.json")
#' j <- j1
#' j <- j2
#' multiqc_list_plots(j)
#' parsed <- RJSONIO::fromJSON(j)
#' dat <- parsed$report_plot_data$mapping_dup_percentage_plot
#' dat <- parsed$report_plot_data$time_metrics_plot
#' dat <- parsed$report_plot_data$bclconvert_lane_counts
#' dat <- parsed$report_plot_data$bclconvert_sample_counts
#' multiqc_parse_bargraph_plot(dat)
#' }
#' @export
multiqc_parse_bargraph_plot <- function(dat) {
  assertthat::assert_that(
    dat[["plot_type"]] == "bar_graph",
    all(c("samples", "datasets") %in% names(dat))
  )
  samp <- dat[["samples"]]
  ds_labs <- rep(NA, length(samp))
  if (!is.null(dat[["config"]][["data_labels"]])) {
    ds_labs <- dat[["config"]][["data_labels"]] |>
      purrr::map_chr("name", .default = NA)
  }
  d <- NULL # binding x to NULL gives us x
  # for each dataset subarray
  for (i in seq_len(length(dat[["datasets"]]))) {
    ds <- dat[["datasets"]][[i]]
    # for each object in the subarray
    for (y in seq_len(length(ds))) {
      assertthat::assert_that(length(ds[[y]][["data"]]) == length(samp[[i]]))
      d <- tibble::tibble(
        label = ds_labs[[i]],
        sample = samp[[i]],
        name = ds[[y]][["name"]],
        data = ds[[y]][["data"]]
      ) |>
        dplyr::bind_rows(d)
    }
  }
  d |>
    tidyr::pivot_wider(names_from = "name", values_from = "data")
}

multiqc_list_plots <- function(x) {
  # handle both parsed and raw json input
  if (rlang::is_bare_atomic(x) && file.exists(x)) {
    x <- RJSONIO::fromJSON(x)[["report_plot_data"]]
  }
  assertthat::assert_that(inherits(x, "list"))
  # keep unsorted as in the original json
  purrr::map_chr(x, "plot_type") |>
    tibble::enframe(name = "plot_nm", value = "plot_type")
}
