library(shiny)
library(DT)

prid_choices <- c("efgh5678", "20250308f448d4e0")
ui <- fluidPage(
  selectInput(
    "prid",
    "PortalRunId",
    choices = prid_choices
  ),
  "Available Tables",
  DT::DTOutput("table_list")
  # DT::DTOutput("table")
)

server <- function(input, output, session) {
  con <- DBI::dbConnect(
    drv = RPostgres::Postgres(),
    dbname = "nemo",
    user = "orcabus"
  )
  prid <- reactive(input$prid)
  output$table_list <- DBI::dbListTables(con) |>
    tibble::as_tibble_col() |>
    DT::renderDT()
  # qc <- dplyr::tbl(con, "tso_sar_qc") |>
  #   dplyr::as_tibble() |>
  #   dplyr::filter(dracarys_id == prid())
  # output$table <- DT::renderDT(qc)
}

shinyApp(ui, server)
