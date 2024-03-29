# Generated by roxytest: do not edit by hand!

# File R/File.R: @testexamples

test_that("Function File() @ L21", {
  
  F1 <- File$new(readr::readr_example("mtcars.csv"))
  (bname_f1 <- F1$bname())
  (F2 <- File$new("https://stratus-gds-aps2/foo/bar/baz.csv?bla"))
  
  expect_true(inherits(F1, c("File", "R6")))
  expect_equal(bname_f1, "mtcars.csv")
  expect_equal(F2$bname(), "baz.csv")
  expect_equal(F2$type(), NA_character_)
  expect_equal(F2$is_url, TRUE)
})

