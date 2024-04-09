#' Read Google LIMS
#'
#' Reads UMCCR's Google LIMS spreadsheet.
#'
#' @return Tibble with all columns and rows from the Google LIMS spreadsheet.
#' @export
glims_excel_read <- function() {
  lims_key <- googledrive::drive_find("^Google LIMS$", shared_drive = "LIMS")$id
  lims <- lims_key |>
    googlesheets4::read_sheet("Sheet1", na = c(".", "", "-"), col_types = "c")
  lims |> readr::type_convert(col_types = readr::cols(.default = "c", Timestamp = "T"))
}
