context("importing manually")

test_data<-read.csv("manual_import.csv")

test_that("as_muscle_stim stops when it should",{
  expect_error(as_muscle_stim(test_data),"Please specify the experiment type")
  expect_error(as_muscle_stim(test_data,"invalid type"),"Invalid experiment type")
  expect_error(as_muscle_stim(test_data,"workloop"),"Position, Force, Stim")
})

# fix names
names(test_data)<-c("Stim","Force","Position")

test_that("sample frequency and/or time are inferred correctly",{
  expect_error(as_muscle_stim(test_data,"workloop"),"sampling frequency")
  expect_equal(attr(as_muscle_stim(cbind(test_data,Time=c(0,1e-3)),"workloop"),"sample_frequency"),1000)
  expect_equal(as_muscle_stim(test_data,"workloop",1000)$Time,c(0,1e-3))
})

test_workloop<-as_muscle_stim(test_data,"workloop",1000)
test_tetanus<-as_muscle_stim(test_data,"tetanus",1000)
test_twitch<-as_muscle_stim(test_data,"twitch",1000)

test_that("data is read in correctly",{
  expect_equal(sapply(test_workloop,sum),c(Stim=3,Force=7,Position=11,Time=1e-3))
})

test_that("attributes are handled correctly",{
  # Correct number of attributes made
  expect_length(attributes(test_workloop),22)
  expect_length(attributes(test_tetanus),17)
  expect_length(attributes(test_twitch),15)

  # Correct number of attributes are NA by default
  expect_equal(sum(is.na(attributes(test_workloop))),16)
  expect_equal(sum(is.na(attributes(test_tetanus))),12)
  expect_equal(sum(is.na(attributes(test_twitch))),10)

  # Passing inappropriate/invalid attributes
  expect_warning(as_muscle_stim(test_data,"workloop",1000,invalid_arg=T),"do not match known attr")
  expect_warning(as_muscle_stim(test_data,"twitch",1000,cycle_frequency=T),"do not match known attr")

  # Pass valid attributes
  expect_equal(attr(as_muscle_stim(test_data,"workloop",1000,cycle_frequency=30),"cycle_frequency"),30)
})
