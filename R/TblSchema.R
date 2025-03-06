#' TblSchema R6 Class
#'
#' @description Encapsulates a tibble with metadata such as column types
#' and primary keys. It allows schema inference, and export to YAML for
#' integration with Django.
#'
#' @examples
#' tbl <- tibble::tibble(
#'   dracarysId = "abcd1234",
#'   foo_num = 123,
#'   foo_dbl = 3.14,
#'   foo_chr = "foobar",
#'   foo_int = 35L,
#'   foo_date = Sys.Date()
#' )
#' s1 <- TblSchema$new(tbl = tbl, name = "Foo", pk = "dracarysId")
#'
#' @testexamples
#' expect_true(inherits(s1, c("TblSchema", "R6")))
#'
#' @export
TblSchema <- R6::R6Class(
  "TblSchema",
  private = list(
    .tbl = NULL,
    .name = NULL,
    .pk = NULL,

    # returns vector of django field data types for the given tbl
    .djangofields = function() {
      sql_types <- DBI::dbDataType(DBI::ANSI(), private$.tbl)
      sql2django <- c(
        TEXT = "TextField",
        DOUBLE = "FloatField",
        INT = "IntegerField",
        DATE = "DateField",
        TIMESTAMP = "DateTimeField"
      )
      assertthat::assert_that(all(sql_types %in% names(sql2django)))
      sql2django[sql_types] |>
        purrr::set_names(names(sql_types))
    }
  ),
  public = list(
    #' @description Create a new TblSchema object.
    #' @param tbl Tibble object.
    #' @param name Name of table.
    #' @param pk Character vector of primary key(s).
    initialize = function(tbl = NULL, name = NULL, pk = NULL) {
      assertthat::assert_that(tibble::is_tibble(tbl))
      assertthat::assert_that(rlang::is_character(name, n = 1))
      assertthat::assert_that(is.null(pk) | all(pk %in% colnames(tbl)))
      private$.tbl <- tbl
      private$.name <- name
      private$.pk <- pk
    },
    #' @description
    #' Schema of tbl.
    schema = function() {
      private$.djangofields() |>
        tibble::enframe(name = "column_name", value = "data_type") |>
        dplyr::mutate(
          primary_key = .data$column_name %in% private$.pk
        ) |>
        tibble::add_column(table_name = private$.name, .before = 1)
    },
    #' @description Print details about the File.
    #' @param ... (ignored).
    print = function(...) {
      cat("#--- TblSchema ---#\n")
      cat(glue("Table Name: {private$.name}"), "\n")
      cat(glue("Primary Key(s): {private$.pk}"), "\n")
      cat(glue("Number of Columns: {ncol(private$.tbl)}"), "\n")
      invisible(self)
    }
  )
)
