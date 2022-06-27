require(tidyverse)
require(arrow)

# Summarise and plot TMB results

d <-
  tibble(
    x = list.files(here::here("nogit/tso"),
                   pattern = "tmb\\.json\\.gz$",
                   recursive = TRUE, full.names = TRUE)) |>
  rowwise() |>
  mutate(
    obj = list(TsoTmbFile$new(x)),
    sample = sub(".tmb.json.gz", "", obj$bname()),
    y = list(read(obj))) |>
  unnest(y) |>
  select(-c(x, obj)) |>
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

arrow::read_parquet(here::here("nogit/tso/tmb/56.TMB.parquet"))


