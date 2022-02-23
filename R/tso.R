#' Read TSO SampleAnalysisResults JSON
#'
#' @param x Path to '_SampleAnalysisResults.json.gz' file.
#'
#' @return
#'
#' @examples
#' x <- here::here("nogit/tso/PRJ210016_L2100355_SampleAnalysisResults.json.gz")
tso_read_sample_analysis_results <- function(x) {
  d <- jsonlite::read_json(x)[["data"]]
  # samp_info <- d[["sampleInformation"]] # leave as-is
  # soft_conf <- d[["softwareConfiguration"]] # simplify nirvanaVersionList$dataSources
  # soft_conf_nvl <- soft_conf$nirvanaVersionList[[1]]
  # soft_conf_ds <- dplyr::bind_rows(soft_conf_nvl$dataSources)
  # samp_met <- d[["sampleMetrics"]]
  vsmall <- d[["variants"]][["smallVariants"]]

  # grab 3 elements from there
  v <- tibble::tibble(v = vsmall[1:3])
  v |>
    unnest_wider(v) |>
    unnest_longer(nirvana) |>
    unnest_wider(nirvana) |>
    unnest(transcripts) |>
    unnest_wider(transcripts) |>
    unnest_longer(consequence)

  n <- purrr::map(v, function(el) {
    l <- el[["nirvana"]][[1]]
    if (length(l$transcripts) == 0) {
      l$transcripts <- list(NA)
    }
    l
  })
}

#' Read TSO tmb JSON file
#'
#' @param x Path to 'tmb.json.gz' file.
#'
#' @return A tibble.
#'
#' @examples
#' x <- here::here("nogit/tso/tmb/MDX200144_L2100176_rerun.tmb.json.gz")
tso_read_tmb <- function(x) {
  j <- jsonlite::read_json(x)
  j$Settings <- NULL
  tibble::as_tibble(j)
}

d <-
  tibble(
    x = list.files(here::here("nogit/tso/tmb"),
                   pattern = "tmb\\.json\\.gz$",
                   recursive = TRUE, full.names = TRUE)) |>
  rowwise() |>
  mutate(sample = sub(".tmb.json.gz", "", basename(x)),
         y = list(tso_read_tmb(x))) |>
  unnest(y) |>
  select(-x) |>
  pivot_longer(TmbPerMb:CodingRegionSizeMb)

theme_set(theme_bw())
d |>
  ggplot(aes(x = "", y = value, label = sample, colour = sample)) +
  geom_violin(fill = "transparent", colour = "grey80", alpha = 0.4) +
  # geom_point() +
  ggforce::geom_sina(aes(
    group = name,
    colour = sample,
  ), seed = 42) +
  facet_wrap(~name, scales = "free")

d |>
  ggplot(aes(x = value)) +
  geom_histogram(fill = "purple", colour = "black") +
  facet_wrap(~name, scales = "free")

# arrow::read_parquet(here::here("nogit/tso/tmb/56.TMB.parquet"))
