context("printing muscle_stim objects")

workloop_example<-read_ddf(system.file("extdata","workloop.ddf",package = 'workloopR'))
tetanus_example<-read_ddf(system.file("extdata","tetanus.ddf",package = 'workloopR'))
twitch_example<-read_ddf(system.file("extdata","twitch.ddf",package = 'workloopR'))

test_that("print works for basic experiment types", {
  expect_output(print(workloop_example),"# Workloop Data: .* with 3238 more rows")
  expect_output(print(tetanus_example,10),"# Tetanus Data: .* with 7491 more rows")
  expect_output(print(twitch_example,20),"# Twitch Data: .* with 3981 more rows")
})

selected_cycles<-select_cycles(workloop_example,"lo",2:4)
analyzed_workloop<-analyze_workloop(selected_cycles)

test_that("print works for transformed and analyzed objects",{
  expect_output(print(selected_cycles),"# Workloop Data: .* with 1065 more rows")
  expect_output(print(analyzed_workloop),"File ID: workloop.ddf\nCycles: 3 cycles kept out of 6\nMean Work: 0.00266 J\nMean Power: 0.0807 W")
})

context("summarizing muscle_stim objects")

test_that("summarizing works for basic experiment types", {
  expect_output(summary(workloop_example),"# Workloop Data:.*File ID: workloop.ddf.*Cycle Frequency: 28Hz.*L0-to-L0")
  expect_output(summary(tetanus_example),"# Tetanus Data:.*data.frame Columns:.*Stimulus Length: 0.081s")
  expect_output(summary(twitch_example),"# Twitch Data:.*Stimulus Offset: 0.1s")
  expect_output(summary(invert_position(workloop_example)),"Position is inverted")
})

test_that("summarizing works for transformed and analyzed objects",{
  expect_output(summary(selected_cycles),"# Workloop Data:.*  Cycle \\(letters\\)")
  expect_output(summary(analyzed_workloop),"# Workloop Data:.*  Cycle        Work  Net_Power")
})
summary(analyzed_workloop)
