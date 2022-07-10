#' Print current timestamp for logging
#'
#' @return Current timestamp as character.
#' @export
date_log <- function() {
  as.character(paste0("[", as.POSIXct(Sys.time()), "]"))
}

#' Create directory
#'
#' @param d Directory to create.
#'
#' @export
mkdir <- function(d) {
  if (!dir.exists(d)) {
    dir.create(d, recursive = TRUE)
  }
}

# copy recursively
cpdir <- function(from, to) {
  mkdir(to)
  file.copy(from = from, to = to, recursive = TRUE)
}
