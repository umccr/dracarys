#' File R6 Class
#'
#' @description File is a base R6 class representing a TSV/CSV/JSON output from
#' a DRAGEN workflow.
#'
#' A File has a path, a basename, a type, and a default read method for its type.
#'
#' @examples
#' F1 <- File$new(readr::readr_example("mtcars.csv"))
#' F1$read(col_types = readr::cols("double"))
#' F1$bname()
#'
#' @export
File <- R6::R6Class("File", public = list(
  #' @field path Name or full path of the file.
  path = NULL,

  #' @description Create a new File object.
  #' @param path Name or full path of the file.
  initialize = function(path = NULL) {
    stopifnot(is.character(path), length(path) == 1)
    self$path <- normalizePath(path)
  },

  #' @description Basename of the file.
  #' @return Basename of the file as a character vector.
  bname = function() {
    basename(self$path)
  },

  #' @description Get the type of file.
  #' @return String describing the type of file (CSV, TSV, JSON or OTHER).
  type = function() {
    nm <- self$path
    dplyr::case_when(
      grepl("\\.json", nm) ~ "JSON",
      grepl("\\.csv", nm) ~ "CSV",
      grepl("\\.tsv", nm) ~ "TSV",
      TRUE ~ "OTHER"
    )
  },

  #' @description Print details about the File.
  #' @param ... (ignored).
  print = function(...) {
    cat("#--- File ---#\n")
    cat(glue("Path: {self$path}"), "\n")
    cat(glue("Basename: {self$bname()}"), "\n")
    cat(glue("Type: {self$type()}"), "\n")
    invisible(self)
  },

  #' @description Read the file based on its type.
  #' @param ... Arguments passed on to appropriate read_* function.
  read = function(...) {
    x <- self$path
    t <- self$type()
    possible_types <- c("CSV", "TSV", "JSON", "OTHER")
    assertthat::assert_that(t %in% possible_types)
    if (t == "CSV") {
      readr::read_csv(x, ...)
    } else if (t == "TSV") {
      readr::read_tsv(x, ...)
    } else if (t == "JSON") {
      jsonlite::read_json(x, ...)
    } else {
      stop(glue("Don't know how to read file of type {t}."))
    }
  }
))

#' Read File object
#'
#' @description Read the file based on its type.
#'
#' @param x Object of class `File`.
#' @param ... Arguments passed on to appropriate read_* function.
#' @export
read.File <- function(x, ...) {
  assertthat::assert_that(inherits(x, "File"))
  x$read(...)
}
