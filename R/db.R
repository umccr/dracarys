#' Test DB
#'
#' @examples
#' \dontrun{
#' path <- file.path(
#'   "~/s3/pipeline-prod-cache-503977275616-ap-southeast-2/byob-icav2/production",
#'   "analysis/wgts-qc/20241123ffa837c4/L2401621_dragen_alignment"
#' )
#' path <- file.path(
#'   "~/s3/pipeline-prod-cache-503977275616-ap-southeast-2/byob-icav2/production",
#'   "analysis/cttsov2/20250117c5b9baa8"
#' )
#' prefix <- "L2500039" # tsov2
#' prefix <- "L2401621" # alignqc
#' outdir <- path
#' max_files <- 1000
#' prid <- "abcd1234"
#' prid <- "efgh5678"
#' dbname <- "nemo"
#' dbuser <- "orcabus"
#' d <- db_test(
#'   path = path, prefix = prefix, outdir = outdir,
#'   prid = prid, max_files = max_files,
#'   dbname = dbname, dbuser = dbuser
#' )
#' }
#'
#' @export
db_test <- function(
  path,
  prefix,
  outdir,
  prid,
  max_files,
  dbname = "nemo",
  dbuser = "orcabus"
) {
  # TODO: add workflow type dispatcher
  obj <- Wf_dragen$new(path = path, prefix = prefix)
  d_dl <- obj$download_files(
    outdir = outdir,
    max_files = max_files
  )
  d_tidy <- obj$tidy_files(d_dl) |>
    # add ID column at the beginning
    dplyr::mutate(
      data = purrr::map(
        .data$data,
        \(x)
          tibble::add_column(
            x,
            dracarys_id = as.character(prid),
            .before = 1
          )
      )
    )
  con <- DBI::dbConnect(
    drv = RPostgres::Postgres(),
    dbname = dbname,
    user = dbuser
  )
  # now write each table to db
  fin <- d_tidy |>
    dplyr::rowwise() |>
    dplyr::mutate(
      write_tbl = list(
        DBI::dbWriteTable(
          conn = con,
          name = .data$name,
          value = .data$data,
          append = T,
          overwrite = F
        )
      )
    ) |>
    dplyr::ungroup()
  DBI::dbDisconnect(con)
  return(fin)
}

#' Generate Django Models from a Schema Tibble
#'
#' @param schema_tbl A tibble containing schema information (from `s1$schema()`)
#' @param output_file Name of the Python file to write the Django models (default: "models.py")
#' @return Writes Django models to a file
#' @examples
#' schema <- tibble::tibble(
#'   table_name = "foo",
#'   column_name = c("dracarysId", "foo_num", "foo_dbl", "foo_chr", "foo_int", "foo_date"),
#'   data_type = c("TextField", "FloatField", "FloatField", "TextField", "IntegerField", "DateField"),
#'   primary_key = rep(c(TRUE, FALSE), c(1L, 5L)),
#' )
#' out <- tempfile()
#' schema2django(s1$schema(), out)
schema2django <- function(schema, out) {
  assertthat::assert_that(tibble::is_tibble(schema))
  assertthat::assert_that(
    all(
      c("table_name", "column_name", "data_type", "primary_key") %in%
        colnames(schema)
    )
  )

  # Generate model definitions
  models <- schema |>
    # dplyr::group_by(.data$table_name) |>
    dplyr::mutate(
      fields = glue(
        "    {.data$column_name} = models.{.data$data_type}({dplyr::if_else(.data$primary_key, 'primary_key=True', '')})"
      )
    ) |>
    dplyr::mutate(code = paste(.data$fields, collapse = "\n")) |>
    dplyr::select("table_name", "code") |>
    dplyr::distinct() |>
    dplyr::mutate(
      code = paste0(
        glue("class {.data$table_name}(models.Model):"),
        "\n",
        .data$code,
        collapse = "\n"
      )
    ) |>
    dplyr::pull(code)

  # Write to file
  model_code <- paste0(
    "from django.db import models\n\n",
    paste0(models, collapse = "\n\n")
  )
  readr::write_lines(model_code, out)
}
