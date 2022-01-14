#' Generic Reading
#'
#' @description Generic function for reading dracarys objects.
#' @param x Object of the respective class.
#' @param ... Additional arguments.
#'
#' @return Possibly a tibble or a list, depending on the object.
#' @export
read <- function(x, ...) {
  UseMethod("read")
}
