context("reading and analyzing single ddf files")

test_that("read_analyze_wl stops when it should", {
  expect_error(read_analyze_wl(system.file("extdata","twitch.ddf",package = 'workloopR')),"does not appear to contain data from a workloop experiment")
  expect_warning(read_analyze_wl(system.file("extdata","workloop.ddf",package = 'workloopR'),invalid_arg=4),"One or more provided attributes")
})

test_that("read_analyze_wl passes arguments appropriately", {
  example_analyzed_wl<-read_analyze_wl(system.file("extdata","workloop.ddf",package = 'workloopR'),file_id="test_id",cycle_def="t2t",GR=2)
  expect_equal(attr(example_analyzed_wl,"file_id"),"test_id")
  expect_equal(attr(example_analyzed_wl,"cycle_def"),"t2t")
  expect_equal(attr(example_analyzed_wl,"gear_ratio"),2)
})

context("reading and analyzing ddf files by directory")

test_that("read_ddf_dir stops when it should",{
  expect_error(read_analyze_wl_dir("non_existant_folder"),"No files matching the pattern")
  expect_warning(read_analyze_wl_dir(system.file("extdata/wl_duration_trials",package="workloopR"),
                                     sort_by='non-existant-attribute'),
                                     "The provided sort_by argument is not a valid attribute")
})

analyzed_wl_dir_example<-read_analyze_wl_dir(system.file("extdata/wl_duration_trials",package="workloopR"),cycle_def="p2p",keep_cycles=2:4,sort_by="file_id")

test_that("read_analyze_wl_dir reads in data correctly", {
  expect_length(analyzed_wl_dir_example,4)
  expect_true(all(sapply(analyzed_wl_dir_example,function(x) class(x)==c("analyzed_workloop","list"))))
})
