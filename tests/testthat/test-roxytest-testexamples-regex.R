# Generated by roxytest: do not edit by hand!

# File R/regex.R: @testexamples

test_that("Function dr_func_eval() @ L96", {
  
  mean_1_to_10 <- dr_func_eval("mean", v = c("mean", "sd"))(1:10)
  expect_equal(mean_1_to_10, base::mean(1:10))
  expect_null(dr_func_eval("foo"))
})

