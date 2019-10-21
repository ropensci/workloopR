context("timing isometric data")

tetanus_example <- read_ddf(system.file("extdata", "tetanus.ddf", package = "workloopR"))
twitch_example <- read_ddf(system.file("extdata", "twitch.ddf", package = "workloopR"))

test_that("isometric timing stops when it should", {
  expect_error(isometric_timing(1), "from an isometric experiment!")
  expect_error(isometric_timing(tetanus_example, "a"), "rising set points")
  expect_error(isometric_timing(tetanus_example, 101), "rising set points")
  expect_error(isometric_timing(twitch_example, 10, T), "relaxing set points")
  expect_error(isometric_timing(twitch_example, 10, -1), "relaxing set points")
})

isometric_timing(twitch_example, 43, 34)
#      file_id time_stim force_stim time_peak force_peak
# 1 twitch.ddf    0.1002    224.067    0.1141   412.4495
#   time_rising_43 force_rising_43 time_relaxing_34
# 1         0.1066        307.1295           0.1375
#   force_relaxing_34
# 1          288.4205


test_that("isometric timing calculates the correct values", {
  expect_equal(isometric_timing(tetanus_example),
    data.frame(
      file_id = "tetanus.ddf",
      time_stim = 0.1002,
      force_stim = 141.0045,
      time_peak = 0.1364,
      force_peak = 2626.289,
      time_rising_10 = 0.1097,
      force_rising_10 = 395.9985,
      time_rising_90 = 0.1282,
      force_rising_90 = 2383.387,
      stringsAsFactors = F
    ),
    tolerance = 1e-4
  )
  expect_equal(isometric_timing(twitch_example, 43, 34),
    data.frame(
      file_id = "twitch.ddf",
      time_stim = 0.1002,
      force_stim = 224.067,
      time_peak = 0.1141,
      force_peak = 412.4495,
      time_rising_43 = 0.1066,
      force_rising_43 = 307.1295,
      time_relaxing_34 = 0.1375,
      force_relaxing_34 = 288.4205,
      stringsAsFactors = F
    ),
    tolerance = 1e-4
  )
})
