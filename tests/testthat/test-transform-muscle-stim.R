context("transforming muscle_stim objects")

test_data <- read.csv("manual_import.csv")
names(test_data) <- c("Stim", "Force", "Position")
test_workloop <- as_muscle_stim(test_data, "workloop", 1000)

test_that("transformation functions work", {
  expect_error(fix_GR(1), "should be of class `muscle_stim`")
  expect_error(fix_GR(test_workloop, "a"), "must be numeric")
  expect_equal(attr(fix_GR(test_workloop, 5), "gear_ratio"), 5)
  expect_error(invert_position(1), "should be of class `muscle_stim`")
  expect_true(attr(invert_position(test_workloop), "position_inverted"))
  expect_equal(invert_position(test_workloop)$Position, c(-5, -6))
})

test_that("trapezoidal integration works", {
  expect_error(trapezoidal_integration("a"), "first argument")
  expect_error(trapezoidal_integration(1, "a"), "second argument")
  expect_error(trapezoidal_integration(2, 1:2), "lengths of the variable")
  expect_equal(trapezoidal_integration(test_workloop$Time,
                                       test_workloop$Force),
               0.0035)
})
