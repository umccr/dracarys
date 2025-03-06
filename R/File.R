#' File R6 Class
#'
#' @description File is a base R6 class representing a TSV/CSV/JSON output from
#' a UMCCR workflow.
#'
#' A File has a path, a basename, a type, and can be a presigned URL.
#'
#' @examples
#' F1 <- File$new(readr::readr_example("mtcars.csv"))
#' (bname_f1 <- F1$bname())
#' (F2 <- File$new("https://stratus-gds-aps2/foo/bar/baz.csv?bla"))
#'
#' @testexamples
#' expect_true(inherits(F1, c("File", "R6")))
#' expect_equal(bname_f1, "mtcars.csv")
#' expect_equal(F2$bname(), "baz.csv")
#' expect_equal(F2$type(), NA_character_)
#' expect_equal(F2$is_url, TRUE)
#'
#' @export
File <- R6::R6Class(
  "File",
  public = list(
    #' @field path Name or full path of the file.
    #' @field is_url Is the file a presigned URL?
    path = NULL,
    is_url = NULL,

    #' @description Create a new File object.
    #' @param is_url Is the file a presigned URL?
    #' @param path Name or full path of the file.
    initialize = function(path = NULL, is_url = NULL) {
      stopifnot(is.character(path), length(path) == 1)
      self$is_url <- is_url(path)
      self$path <- base::ifelse(is_url(path), path, normalizePath(path))
    },

    #' @description Basename of the file.
    #' @return Basename of the file as a character vector.
    bname = function() {
      x <- self$path
      if (is_url(x)) {
        x <- strsplit(self$path, "\\?")[[1]][1]
      }
      basename(x)
    },

    #' @description Get the type of file.
    #' @return String describing the specific type of dracarys file (NA if not a dracarys-recognised file).
    type = function() {
      nm <- self$bname()
      match_regex(nm)
    },

    #' @description Print details about the File.
    #' @param ... (ignored).
    print = function(...) {
      cat("#--- File ---#\n")
      cat(glue("Path: {self$path}"), "\n")
      cat(glue("Basename: {self$bname()}"), "\n")
      cat(glue("Type: {self$type()}"), "\n")
      cat(glue("isURL: {self$is_url}"), "\n")
      invisible(self)
    }
  )
)

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

is_url <- function(x) {
  grepl("(http|https)://[a-zA-Z0-9./?=_%:-]*", x)
}
