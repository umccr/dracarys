awsvault_profile <- function(p) {
  assertthat::assert_that(system("command -v aws-vault") == 0, msg = "aws-vault needs to be installed - see https://github.com/99designs/aws-vault")
  creds <- file.path(tempdir(), "creds.tsv")
  system(glue("unset AWS_VAULT && aws-vault exec {p} -- env | grep AWS > {creds}"))
  x <- readr::read_tsv(creds, col_names = "var_eq_value", col_types = "c") |>
    tidyr::separate_wider_delim("var_eq_value", delim = "=", names = c("var", "value"), too_many = "merge") |>
    tidyr::pivot_wider(names_from = "var", values_from = "value") |>
    as.list()

  Sys.setenv(AWS_VAULT = p)
  Sys.setenv(AWS_ACCESS_KEY_ID = x[["AWS_ACCESS_KEY_ID"]])
  Sys.setenv(AWS_SECRET_ACCESS_KEY = x[["AWS_SECRET_ACCESS_KEY"]])
  Sys.setenv(AWS_SESSION_TOKEN = x[["AWS_SESSION_TOKEN"]])
  Sys.setenv(AWS_DEFAULT_REGION = x[["AWS_DEFAULT_REGION"]])
  Sys.setenv(AWS_REGION = x[["AWS_REGION"]])
  cat(glue::glue("Switched to AWS Profile: {p}"))
}
