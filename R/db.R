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
#'   "analysis/cttsov2/20250308f448d4e0"
#' )
#' prefix <- "L2500183" # tsov2
#' prefix <- "L2401621" # alignqc
#' outdir <- path
#' max_files <- 1000
#' prid <- "20250308f448d4e0"
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
  obj <- Wf_cttsov2$new(path = path, prefix = prefix)
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
  res <- d_tidy |>
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
  return(res)
}
