context("selecting cycles")

workloop_example<-read_ddf(system.file("extdata","workloop.ddf",package = 'workloopR'))

test_that("select_cycles stops when it should", {
  expect_error(select_cycles(1),"Input data should be of class")
  expect_warning(select_cycles(workloop_example),"Cycle definition not supplied")
  expect_error(select_cycles(workloop_example,"a"),"Invalid cycle definition")
  expect_error(select_cycles(workloop_example,"p2p","a"),"keep_cycles should be numeric")
  expect_warning(select_cycles(workloop_example,"p2p",-1),"keep_cycles argument includes cycles that don't exist")
  attr(workloop_example,"cycle_frequency")<-NA
  expect_error(select_cycles(workloop_example,"p2p"),"Length-out cycle frequency is needed")
})

attr(workloop_example,"cycle_frequency")<-28

test_that("cycles are selected properly",{
  expect_length(unique(select_cycles(workloop_example,"p2p",1:4)$Cycle),4)
  expect_equal(attr(select_cycles(workloop_example,"lo",1),"cycle_def"),"lo")
  expect_equal(attr(select_cycles(workloop_example,"p2p",1),"cycle_def"),"p2p")
  expect_equal(attr(select_cycles(workloop_example,"t2t",1),"cycle_def"),"t2t")
})

context("analyzing workloops")

selected_cycles<-select_cycles(workloop_example,"lo",2:4)

test_that("analyze_workloop stops when it should", {
  expect_error(analyze_workloop(1),"Input data should be of class")
  expect_error(analyze_workloop(workloop_example),"The Cycle column")
  expect_error(analyze_workloop(selected_cycles,F,"a"),"Gear ratio")
})

test_that("simplified analyze_workloop output is accurate",{
  expect_equal(analyze_workloop(selected_cycles,T),
    data.frame(
      Cycle=c("A","B","C"),
      Work=c(0.002263076,0.002708145,0.003007098),
      Net_Power=c(0.06903262,0.08220278,0.09085517),
      row.names=c("a","b","c")
    ),
    tolerance=1e-4
    )
})

analyzed_workloop<-analyze_workloop(selected_cycles)

test_that("full analyze_workloop output is accurate",{
  expect_equal(class(analyzed_workloop),c("analyzed_workloop","list"))
  expect_true(all(sapply(analyzed_workloop,function(x)class(x)==c("workloop","muscle_stim","data.frame"))))
  expect_true(all(sapply(analyzed_workloop,function(x)names(x)==c("Time","Position","Force","Stim","Cycle","Inst_Velocity","Filt_Velocity","Inst_Power","Percent_of_Cycle"))))
  expect_equal(suppressWarnings(sapply(analyzed_workloop[[1]],mean,na.rm=T)),
               c(Time=6.570000e-02,
                 Position=4.985869e-01,
                 Force=1.348556e+03,
                 Stim=2.240896e-02,
                 Cycle=NA,
                 Inst_Velocity=1.322893e-03,
                 Filt_Velocity=1.123992e-02,
                 Inst_Power=6.903262e-02,
                 Percent_of_Cycle=5.000000e+01),
               tolerance=1e-3)
})

context("summarizing multiple analyzed workloop trials")

analyzed_wl_list<-read_analyze_wl_dir(system.file("extdata/wl_duration_trials",package="workloopR"),phase_from_peak=TRUE,cycle_def="p2p",keep_cycles=2:4,sort_by="file_id")
analyzed_wl_list[[5]]<-workloop_example

test_that("summarize wl trials stops when it should", {
  expect_error(summarize_wl_trials(workloop_example),"Please provide a list")
  expect_error(summarize_wl_trials(analyzed_wl_list),"includes elements that are not analyzed")
})

analyzed_wl_list[[5]]<-NULL
swl_output<-summarize_wl_trials(analyzed_wl_list)

test_that("summarize wl trials pulls the correct data",{
  # drop mtime column, which will depend on the installation datetime of package
  expect_equal(swl_output[,-7],
               data.frame(
                 File_ID=c("01_4pulse.ddf","02_2pulse.ddf","03_6pulse.ddf","04_4pulse.ddf"),
                 Cycle_Frequency=rep(28,4),
                 Amplitude=rep(3.15,4),
                 Phase=c(-24.36,-24.64,-24.92,-24.64),
                 Stimulus_Pulses=c(4,2,6,4),
                 Stimulus_Frequency=rep(300,4),
                 Mean_Work=c(0.0028362363,0.0009686570,-0.0001310863,0.0024082708),
                 Mean_Power=c(0.078967198,0.026247519,-0.004017894,0.066959552)),
               tolerance=1e-4)
})

context("time correcting for muscle degradation")

test_that("time_correct stops when it should",{
  expect_error(time_correct(1),"provide a data.frame")
  expect_error(time_correct(analyze_workloop(selected_cycles,T)),"provide summarized workloop trial")
})

test_that("time_correct calculates the correct values",{
  # standardize mtime for consistent time correction
  swl_output$mtime<-1:4
  # only check the newly calculated columns
  expect_equal(time_correct(swl_output)[,10:11],
               data.frame(
                 Time_Corrected_Work=c(0.002836236,0.001111312,0.000154224,0.002836236),
                 Time_Corrected_Power=c(0.078967198,0.030250068,0.003987203,0.078967198)),
               tolerance=1e-4)
})
