#' Print current timestamp for logging
#'
#' @return Current timestamp as character.
#' @export
date_log <- function() {
  as.character(glue::glue('[{format(Sys.time(), "%Y-%m-%dT%H:%M:%S%Z")}]'))
}
