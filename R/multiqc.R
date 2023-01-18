#' Dracarys MultiQC Tidy
#'
#' Generate tidier representations of MultiQC JSON output
#' @param json Path to `multiqc_data.json`.
#' @param prefix Prefix for output files.
#' @param outdir Path to output results.
#' @param out_format Format of output (tsv, parquet, both) (def: tsv).
#' @return Generates TSV and/or Parquet representations of the input
#' MultiQC JSON file.
#' @export
dracarys_multiqc <- function(json, prefix, outdir, out_format = "tsv") {
  output_format_valid(out_format)
  e <- emojifont::emoji
  cli::cli_div(theme = list(
    span.file = list(color = "lightblue"),
    span.emph = list(color = "orange")
  ))
  cli::cli_alert_info("{date_log()} {e('dragon')} Start tidying {.file {json}} {e('fire')}")
  # main dracarys function
  d1 <- multiqc_tidy_json(json)
  if (out_format %in% c("tsv", "both")) {
    tsv_out <- file.path(outdir, glue("{prefix}.tsv"))
    readr::write_tsv(d1, tsv_out)
  }
  if (out_format %in% c("parquet", "both")) {
    parquet_out <- file.path(outdir, glue("{prefix}.parquet"))
    arrow::write_parquet(d1, parquet_out)
  }

  cli::cli_alert_success("{date_log()} {e('rocket')} End tidying {.file {json}} {e('comet')}!")
  cli::cli_alert_info("{date_log()} {e('tada')} Path to output directory with results for {.emph {prefix}}: {.file {outdir}}")
}

#' Tidy MultiQC JSON
#'
#' Tidies 'multiqc_data.json' output from MultiQC.
#' Modified from the awesome <https://github.com/multimeric/TidyMultiqc>.
#' @param j Path to `multiqc_data.json`.
#' @return A tidy tibble where each column corresponds to a single metric,
#' and each row corresponds to a single sample.
#' @export
multiqc_tidy_json <- function(j) {
  p <- RJSONIO::fromJSON(j)
  cdate <- p[["config_creation_date"]] %||% "UNKNOWN"
  cdate <- sub(", ", "_", cdate)
  workflow <- .multiqc_guess_workflow(p)
  d <- dracarys::multiqc_parse_gen(p)
  if (workflow == "dragen_umccrise") {
    # replace the "NA" strings with NA, else we get a column class error
    # due to trying to bind string ('NA') with numeric.
    # https://stackoverflow.com/questions/35292757/replace-values-in-list
    d <- rapply(d, function(x) ifelse(x == "NA", NA, x), how = "replace")
  }
  d <- d |>
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
  return(.multiqc_rename_cols(d))
}

.multiqc_rename_cols <- function(d) {
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
  p[[el]] |> purrr::reduce(~ purrr::list_merge(.x, !!!.y))
}

#' Parse MultiQC 'report_saved_raw_data' JSON Element
#'
#' Parses MultiQC 'report_saved_raw_data' JSON Element. Modified from the
#' awesome <https://github.com/multimeric/TidyMultiqc>.
#' @param p Parsed MultiQC JSON.
#' @param tools2delete Character vector of tools to delete from the parsed JSON.
#' @return A list.
#' @export
multiqc_parse_raw <- function(p, tools2delete = NULL) {
  el <- "report_saved_raw_data"
  assertthat::assert_that(inherits(p, "list"), el %in% names(p))
  tool_nms <- names(p[[el]])
  if (!is.null(tools2delete)) {
    assertthat::assert_that(
      is.vector(tools2delete), !is.list(tools2delete),
      all(tools2delete %in% tool_nms)
    )
    # remove given elements from list
    p[[el]] <- base::within(p[[el]], rm(list = tools2delete))
  }
  # For each tool
  p[[el]] |>
    purrr::imap(function(samples, tool) {
      # For each sample
      samples |> multiqc_kv_map(function(metrics, sample) {
        # For each metric in the above tool
        list(
          key = sample,
          value = metrics |> multiqc_kv_map(function(mvalue, mname) {
            # Sanitise metric names
            combined_metric <- list(
              # key = paste0(tool, ".", mname),
              key = mname,
              value = mvalue
            )
          })
        )
      })
    }) |>
    purrr::reduce(utils::modifyList)
}

multiqc_kv_map <- function(l, func) {
  mapped <- purrr::imap(l, func) |>
    purrr::set_names(nm = NULL)
  keys <- mapped |> purrr::map_chr("key")
  vals <- mapped |> purrr::map("value")
  vals |> purrr::set_names(keys)
}

# Tibble with three columns: workflow name, raw and cleaned multiqc metric name.
# For consistency, make sure that names are in snake_case and that the tool that
# generates them is added as the suffix e.g. awesome_metric_tool1.
MULTIQC_COLUMNS <- system.file("extdata/multiqc_column_map.tsv", package = "dracarys") |>
  readr::read_tsv(col_types = "ccc")
