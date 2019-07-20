context("importing ddf files")

test_that("read_ddf stops when it should", {
  expect_error(read_ddf(),"A filename is required")
  expect_error(read_ddf(system.file("extdata","twitch001.ddf",package = 'workloopR')),"File * not found!")
  expect_error(read_ddf(system.file("CITATION",package = 'workloopR')),"DMC Datafile")
  expect_error(read_ddf("invalid_exp.ddf"),"Could not parse experiment type")
  expect_warning(read_ddf("bad_units.ddf"),"Non-standard units detected")
  expect_warning(read_ddf("non-numeric.ddf"),"includes non-numeric data")
})

workloop_example<-read_ddf(system.file("extdata","workloop.ddf",package = 'workloopR'))
tetanus_example<-read_ddf(system.file("extdata","tetanus.ddf",package = 'workloopR'))
twitch_example<-read_ddf(system.file("extdata","twitch.ddf",package = 'workloopR'))

test_that("read_ddf parses experiments correctly", {
  expect_equal(class(workloop_example),c("workloop","muscle_stim","data.frame"))
  expect_equal(class(tetanus_example),c("tetanus","isometric","muscle_stim","data.frame"))
  expect_equal(class(twitch_example),c("twitch","isometric","muscle_stim","data.frame"))
})

test_that("read_ddf reads data in correctly",{
  expect_equal(names(workloop_example),c("Time","Position","Force","Stim"))
  expect_equal(sapply(workloop_example,sum),c(Time=526.339,Position=1622.828,Force=3267965.800,Stim=48.000),tolerance=1e-3)
  expect_equal(sapply(tetanus_example,sum),c(Time=2813.625,Position=4448.732,Force=2210595.621,Stim=44.000),tolerance=1e-3)
  expect_equal(sapply(twitch_example,sum),c(Time=800.6001,Position=-24008.8626,Force=950311.5720,Stim=2.0000),tolerance=1e-3)
})

test_that("read_ddf reads in attributes correctly",{
  expect_false(any(is.na(attributes(workloop_example))))
  expect_length(attributes(workloop_example),22)

  expect_false(any(is.na(attributes(tetanus_example))))
  expect_length(attributes(tetanus_example),17)

  expect_false(any(is.na(attributes(twitch_example))))
  expect_length(attributes(twitch_example),15)
})

context("importing ddf files by directory")

test_that("read_ddf_dir stops when it should",{
  expect_error(read_ddf_dir("non_existant_folder"),"No files matching the pattern")
  expect_warning(read_ddf_dir(system.file("extdata/wl_duration_trials",package="workloopR"),
                              sort_by='non-existant-attribute'),
                 "The provided sort_by argument is not a valid attribute")
})

workloop_dir_example<-read_ddf_dir(system.file("extdata/wl_duration_trials",package = 'workloopR'),sort_by="file_id")

test_that("read_ddf_dir parses experiments correctly",{
  expect_true(all(sapply(workloop_dir_example,function(x) class(x)==c("workloop","muscle_stim","data.frame"))))
})

test_that("read_ddf_dir reads data in correctly",{
  expect_equal(names(workloop_dir_example[[1]]),c("Time","Position","Force","Stim"))
  expect_equal(sapply(workloop_dir_example,sum),c(1404208.3,890138.2,2253317.5,1250923.5),tolerance=1e-3)
})

test_that("read_ddf_dir reads in attributes correctly",{
  expect_false(any(sapply(workloop_dir_example,function(x)any(is.na(attributes(x))))))
  expect_equal(mean(sapply(workloop_dir_example,function(x)length(attributes(x)))),22)
})

context("pulling metadata")

test_that("get wl metadata pulls return is correct dimensions", {
  expect_equal(names(get_wl_metadata(".")), c("size","isdir","mode","mtime","ctime","atime","uid","gid","uname","grname","exp_names"))
  expect_equal(nrow(get_wl_metadata(".")),4)
})
