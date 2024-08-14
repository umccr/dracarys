#' Workflow R6 Class
#'
#' @description Workflow is a base R6 class representing a bioinformatic
#' workflow run from a UMCCR workflow manager.
#'
#' @examples
#' p1 <- system.file("extdata/portaldb_workflow_top4.rds", package = "rportal") |>
#'   readRDS() |>
#'   dplyr::filter(type_name == "umccrise") |>
#'   dplyr::slice(1)
#' w <- Wf$new(
#'   prid = p1$portal_run_id, type = p1$type_name, start = p1$start, end = p1$end,
#'   status = p1$end_status, input = p1$input, output = p1$output
#' )
#' w
#' @export
Wf <- R6::R6Class(
  "Wf",
  public = list(
    #' @field prid Portal run ID.
    #' @field type Workflow type.
    #' @field start Workflow start datetime.
    #' @field end Workflow end datetime.
    #' @field status Workflow end status.
    #' @field input Workflow input JSON string.
    #' @field output Workflow output JSON string.
    prid = NULL,
    type = NULL,
    start = NULL,
    end = NULL,
    status = NULL,
    input = NULL,
    output = NULL,
    #' @description Create a new Workflow object.
    #' @param prid Portal run ID.
    #' @param type Workflow type.
    #' @param start Workflow start datetime.
    #' @param end Workflow end datetime.
    #' @param status Workflow end status.
    #' @param input Workflow input JSON string.
    #' @param output Workflow output JSON string.
    initialize = function(prid = NULL, type = NULL, start = NULL, end = NULL,
                          status = NULL, input = NULL, output = NULL) {
      types <- c(
        "bcl_convert",
        "tso_ctdna_tumor_only",
        "wgs_alignment_qc",
        "wts_alignment_qc",
        "wts_tumor_only",
        "wgs_tumor_normal",
        "umccrise",
        "rnasum",
        "star_alignment",
        "oncoanalyser_wts",
        "oncoanalyser_wgs",
        "oncoanalyser_wgts_existing_both",
        "sash"
      )
      assertthat::assert_that(
        type %in% types
      )
      self$prid <- prid
      self$type <- type
      self$start <- start
      self$end <- end
      self$status <- status
      self$input <- input
      self$output <- output
    },
    #' @description Print details about the Workflow.
    #' @param ... (ignored).
    print = function(...) {
      res <- tibble::tribble(
        ~var, ~value,
        "prid", self$prid,
        "type", self$type,
        "start", as.character(self$start),
        "end", as.character(self$end),
        "status", self$status,
      )
      print(res)
      invisible(self)
    }
  ) # end public
)
