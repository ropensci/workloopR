# custom functions
# all written by Vikram B. Baliga (vbaliga@zoology.ubc.ca) and Shreeram
# Senthivasan
# last updated: 2019-07-20


############################ trapezoidal integration ###########################

#' Approximate the definite integral via the trapezoidal rule
#'
#' Mostly meant for internal use in our analysis functions, but made available
#' for other use cases. Accordingly, it does not strictly rely on objects of
#' class \code{muscle_stim}.
#'
#' @param x a variable, e.g. vector of positions
#' @param f integrand, e.g. vector of forces
#'
#' @details In the functions \code{analyze_workloop()}, \code{read_analyze_wl()}
#' , and \code{read_analyze_wl_dir()}, work is calculated as the difference
#' between the integral of the upper curve and the integral of the lower curve
#' of a work loop.
#'
#' @references Atkinson, Kendall E. (1989), An Introduction to Numerical
#' Analysis (2nd ed.), New York: John Wiley & Sons
#'
#' @author Vikram B. Baliga
#'
#' @seealso
#' \code{\link{analyze_workloop}},
#' \code{\link{read_analyze_wl}},
#' \code{\link{read_analyze_wl_dir}}
#'
#' @examples
#'
#' # create a circle centered at (x = 10, y = 20) with radius 2
#' t <- seq(0, 2*pi, length = 1000)
#' coords <- t(rbind(10 + sin(t)*2, 20 + cos(t)*2))
#'
#' # sanity check: does it look like a circle?
#' #plot(coords, asp = 1)
#'
#' # use the function to get the area
#' trapezoidal_integration(coords[,1],coords[,2])
#'
#' # does it match (pi * r^2)?
#' 3.14159265358 * (2^2) # very close
#'
#'
#' @export
trapezoidal_integration <- function(x,
                                    f)
{
  if(!is.numeric(x))
    stop('The variable (first argument) is not numeric.')
  if(!is.numeric(f))
    stop('The integrand (second argument) is not numeric.')
  if (length(x) != length(f))
    stop('The lengths of the variable and the integrand are not equal.')

  # obtain length of variable of integration and integrand
  n=length(x)
  # integrate using the trapezoidal rule
  integral=0.5*sum((x[2:n]-x[1:(n-1)])*(f[2:n]+f[1:(n-1)]))
  # return the definite integral
  return(integral)
}


########################### work loop data analysis ##########################

#' Analyze work loop object to compute work and power output
#'
#' Compute work and power output from a work loop experiment on a per-cycle
#' basis.
#'
#' @param x A \code{workloop} object of class \code{muscle_stim} that has been
#' passed through \code{select_cycles}. See Details.
#' @param simplify Logical. If \code{FALSE}, the full analyzed workloop
#' object is returned. If \code{TRUE} a simpler table of net work and power
#' (by cycle) is returned.
#' @param GR Gear ratio, set to 1 by default
#' @param M Velocity multiplier, set adjust the sign of velocity. This parameter
#'  should generally be either -1 (the default) or 1.
#' @param vel_bf Critical frequency (scalar) for low-pass filtering of velocity
#' via \code{signal::butter()}
#' @param ... Additional arguments potentially passed down from
#' \code{read_analyze_wl()} or \code{read_analyze_wl_dir()}
#'
#' @details Please note that \code{select_cycles()} must be run on data prior to
#' using this function. This function relies on the input \code{muscle_stim}
#' object being organized by cycle number.
#'
#' The \code{muscle_stim} object (\code{x}) must be a \code{workloop},
#' preferably read in by one of our data import functions. Please see
#' documentation for \code{as_muscle_stim()} if you need to manually construct
#' a \code{muscle_stim} object from a non .ddf source.
#'
#' The gear ratio (GR) and velocity multiplier (M) parameters can help correct
#' for issues related to the magnitude and sign of data collection. By default,
#' they are set to apply no gear ratio adjustment and to positivize velocity.
#' Instanteous velocity is often noisy and the \code{vel_bf} parameter allows
#' for low-pass filtering of velocity data. See \code{signal::butter()} and
#' \code{signal::filtfilt()} for details of how filtering is achieved.
#'
#' Please also be careful with units! Se Warning section below.
#'
#' @section Warning:
#' Most systems we have enountered record Position data in millimeters
#' and Force in millinewtons, and therefore this function assumes data are
#' recorded in those units. Through a series of internal conversions, this
#' function computes velocity in meters/sec, work in Joules, and power in
#' Watts. If your raw data do not originate in millimeters and millinewtons,
#' please transform your data accordingly and ignore what you see in the
#' attribute \code{units}.
#'
#' @return
#' The function returns a \code{list} of class \code{analyzed_workloop}
#' that provides instantaneous velocity, a smoothed velocity, and computes work,
#'  instantaneous power, and net power from a work loop experiment. All data are
#'  organized by the cycle number and important metadata are stored as
#'  Attributes.
#'
#' Within the \code{list}, each entry is labeled by cycle and includes:
#' \item{Time}{Time, in sec}
#' \item{Position}{Length change of the muscle, corrected for gear ratio, in mm}
#' \item{Force}{Force, corrected for gear ratio, in mN}
#' \item{Stim}{When stimulation occurs, on a binary scale}
#' \item{Cycle}{Cycle ID, as a letter}
#' \item{Inst_velocity}{Instanteous velocity, computed from \code{Position}
#' change, reported in meters/sec}
#' \item{Filt_velocity}{Instaneous velocity, after low-pass filtering, again in
#' meter/sec}
#' \item{Inst_Power}{Instantaneous power, a product of \code{Force} and
#' \code{Filt_velocity}, reported in J}
#' \item{Percent_of_Cycle}{The percent of that particular cycle which has
#' elapsed}
#'
#' In addition, the following information is stored in the
#' \code{analyzed_workloop} object's attributes:
#' \item{stimulus_frequency}{Frequency at which stimulus pulses occurred}
#' \item{cycle_frequency}{Frequency of oscillations (assuming sine wave
#' trajectory)}
#' \item{total_cycles}{Total number of oscillatory cycles (assuming sine wave
#' trajectory) that the muscle experienced.}
#' \item{cycle_def}{Specifies what part of the cycle is understood as the
#' beginning and end. There are currently three options:
#' 'lo' for L0-to-L0;
#' 'p2p' for peak-to-peak; and
#' 't2t' for trough-to-trough}
#' \item{amplitude}{Amplitude of length change (assuming sine wave
#' trajectory)}
#' \item{phase}{Phase of the oscillatory cycle (in percent) at which stimulation
#' occurred. Somewhat experimental, please use with caution}
#' \item{position_inverted}{Logical; whether position inversion has been
#' applied)}
#' \item{units}{The units of measurement for each column in the object after
#' running this function. See Warning}
#' \item{sample_frequency}{Frequency at which samples were collected}
#' \item{header}{Additional information from the header}
#' \item{units_table}{Units from each Channel of the original ddf file}
#' \item{protocol_table}{Protocol in tabular format; taken from the original
#' ddf file}
#' \item{stim_table}{Specific info on stimulus protocol; taken from the original
#' ddf file}
#' \item{stimulus_pulses}{Number of sequential pulses within a stimulation
#' train}
#' \item{stimulus_offset}{Timing offset at which stimulus began}
#' \item{gear_ratio}{Gear ratio applied by this function}
#' \item{file_id}{Filename}
#' \item{mtime}{Time at which file was last modified}
#' \item{retained_cycles}{Which cycles were retained, as numerics}
#' \item{summary}{Simple table showing work (in J) and net power (in W) for each
#'  cycle}
#'
#' @author Vikram B. Baliga and Shreeram Senthivasan
#'
#' @family data analyses
#' @family workloop functions
#'
#' @references Josephson RK. 1985. Mechanical Power output from Striated Muscle
#'  during Cyclic Contraction. Journal of Experimental Biology 114: 493-512.
#'
#' @examples
#'
#' library(workloopR)
#'
#' # import the workloop.ddf file included in workloopR
#' wl_dat <-read_ddf(system.file("extdata", "workloop.ddf",
#'                               package = 'workloopR'))
#'
#' # select cycles 3 through 5 via the peak-to-peak definition
#' wl_selected <- select_cycles(wl_dat, cycle_def = "p2p", keep_cycles = 3:5)
#'
#' # run the analysis function and get the full object
#' wl_analyzed <- analyze_workloop(wl_selected, GR = 2)
#'
#' # print methods give a short summary
#' print(wl_analyzed)
#'
#' # summary provides a bit more detail
#' summary(wl_analyzed)
#'
#' # run the analysis but get the simplified version
#' wl_analyzed_simple <- analyze_workloop(wl_selected, simplify = TRUE, GR = 2)
#'
#' @seealso
#' \code{\link{read_ddf}},
#' \code{\link{read_analyze_wl}},
#' \code{\link{select_cycles}}
#'
#' @export
#'
analyze_workloop <- function(x,
                             simplify = FALSE,
                             GR = 1,
                             M = -1,
                             vel_bf = 0.05,
                             ...){
  if(!any(class(x)=="workloop"))
    stop("Input data should be of class `workloop`")
  if(!any(names(x)=='Cycle'))
    stop('The Cycle column is missing with no default. Please use select_cycles() to generate this column or check that the column is named correctly.')
  if(!is.numeric(GR))
    stop('Gear ratio (GR) must be numeric')
  if(!is.numeric(M))
    stop('Velocity multiplier (M) must be numeric and is recommended to be either -1 or 1.')

  # transform variables
  x<-fix_GR(x,GR)

  # first chop up the data by cycle:
  cycle_names<-unique(x$Cycle)
  x_by_cycle<-lapply(cycle_names,function(cycle)x[x$Cycle==cycle,])

  # create a percent cycle index column
  percent_of_cycle<-lapply(x_by_cycle,function(x)seq(0,100,100/(nrow(x)-1)))

  # work is calculated as the path integral of Force with respect to Position
  # (displacement)
  # Position and Force are each divided by 1000 to convert mm to meters and mN
  # to N prior to taking the integral. This ensures that the integral reports
  # work in J. The negative is used to match conventions for work
  work<-lapply(x_by_cycle,function(x)-trapezoidal_integration(x$Position/1000,
                                                              x$Force/1000))
  names(work)<-cycle_names

  # velocity is the instantanous change in length (i.e. position) multiplied by sampling frequency
  # the result is divided by 1000 to convert to m/s and multiplied by the velocity multiplier, M
  velocity<-lapply(x_by_cycle,function(x)(x$Position-c(NA,utils::head(x$Position,-1)))*attributes(x)$sample_frequency/1000*M)

  # apply a butterworth filter to velocity to smooth it out a bit
  buttah<-signal::butter(2,vel_bf)
  filt_velocity<-lapply(velocity,function(v)c(NA,signal::filtfilt(buttah,v[-1])))

  # instantaneous power is calculated as the product of instantaneous velocity
  # and force. However since velocity is calculated between two time points,
  # corresponding pairs of force measurements are averaged first
  # the result is divided by 1000 to convert mW to W
  instant_power<-mapply(function(x,v)x$Force*v/1000,x_by_cycle,filt_velocity,
                        SIMPLIFY=FALSE)

  # net power is simply the mean of all instantaneous power
  net_power<-lapply(instant_power,mean,na.rm=TRUE)

  # Early escape for simplified output
  summary_table<-data.frame(
      Cycle=paste0(toupper(cycle_names)),
      Work=unlist(work),
      Net_Power=unlist(net_power)
    )
  if(simplify) return(summary_table)

  # combine everything into one useful object
  result<-mapply(
    function(x,v,filt_v,w,inst_p,net_p,perc){
      x$Inst_Velocity<-v
      x$Filt_Velocity<-filt_v
      x$Inst_Power<-inst_p
      x$Percent_of_Cycle<-perc
      attr(x,"work")<-w
      attr(x,"net_power")<-net_p
      if(!all(is.na(attr(x,"units"))))
        attr(x,"units")<-c(attr(x,"units"),"m/s","m/s","W")
      x
    },
    x_by_cycle,
    velocity,
    filt_velocity,
    work,
    instant_power,
    net_power,
    percent_of_cycle,
    SIMPLIFY=FALSE)
  attr(x,"row.names")<-attr(x,"names")<-NULL
  attributes(result)<-attributes(x)
  attr(result,"summary")<-summary_table
  class(result)<-c("analyzed_workloop","list")

  stats::setNames(result,paste0("cycle_",cycle_names))
}


################################# time correct #################################

#' Time correction for work loop experiments
#'
#' Correct for potential degradation of muscle over time.
#'
#' @param x A \code{data.frame} with summary data, e.g. an object created by
#' \code{summarize_wl_trials()}.
#'
#' @details This function assumes that across a batch of successive trials, the
#' stimulation parameters for the first and final trials are identical. If not,
#' DO NOT USE. Decline in power output is therefore assumed to be a linear
#' function of time. Accordingly, the difference between the final and first
#' trial's (absolute) power output is used to 'correct' trials that occur in
#' between, with explicit consideration of run order and time elapsed (via
#' mtime). A similar correction procedure is applied to work.
#'
#' @return A \code{data.frame} that additionally contains:
#' \item{Time_Corrected_Work }{Time corrected work output, transformed from
#'  \code{$Mean_Work}}
#' \item{Time_Corrected_Power }{Time corrected net power output, transformed
#' from \code{$Mean_Power}}
#'
#' And new attributes:
#' \item{power_difference }{Difference in mass-specific net power output
#' between the final and first trials.}
#' \item{time_difference }{Difference in mtime between the final and first
#' trials.}
#' \item{time_correction_rate }{Overall rate; \code{power_difference} divided
#'  by \code{time_difference}.}
#'
#' @author Vikram B. Baliga and Shreeram Senthivasan
#'
#' @family workloop functions
#' @family batch analyses
#'
#' @examples
#'
#' library(workloopR)
#'
#' # batch read and analyze files included with workloopR
#' analyzed_wls <- read_analyze_wl_dir(system.file("extdata/wl_duration_trials",
#'                                                 package = 'workloopR'),
#'                                     cycle_def = "p2p", keep_cycles = 2:4)
#'
#' # now summarize
#' summarized_wls <- summarize_wl_trials(analyzed_wls)
#'
#'
#' # mtimes within the package are not accurate, so we'll supply
#' # our own vector of mtimes
#' summarized_wls$mtime <- read.csv(
#'                           system.file(
#'                             "extdata/wl_duration_trials/ddfmtimes.csv",
#'                             package="workloopR"))$mtime
#'
#' # now time correct
#' timecor_wls <- time_correct(summarized_wls)
#' timecor_wls
#'
#' # or on your own directory
#' #my_meta <- read_analyze_wl_dir("./my/file/path/")
#' #my_summaries <- summarize_wl_trials(my_meta)
#' #my_timecors <- time_correct(my_summaries)
#'
#' @seealso
#' \code{\link{summarize_wl_trials}}
#'
#' @export
time_correct <- function(x){
  if(class(x)[[1]]!="data.frame")
    stop("Please provide a data.frame of summarized workloop trial data generated by summarize_wl_trials")
  if(!all(c("Mean_Work","Mean_Power","mtime") %in% names(x)))
    stop("Please provide summarized workloop trial data generated by summarize_wl_trials")

  x$Time_Corrected_Work <-
    x$Mean_Work-(utils::tail(x$Mean_Work,1)-utils::head(x$Mean_Work,1)) /
    (utils::tail(x$mtime,1)-utils::head(x$mtime,1))*(x$mtime-utils::head(x$mtime,1))
  x$Time_Corrected_Power <-
    x$Mean_Power-(utils::tail(x$Mean_Power,1)-utils::head(x$Mean_Power,1)) /
    (utils::tail(x$mtime,1)-utils::head(x$mtime,1))*(x$mtime-utils::head(x$mtime,1))
  attr(x,"power_difference") <-
    utils::tail(x$Mean_Power,1)-utils::head(x$Mean_Power,1)
  attr(x,"time_difference") <-
    utils::tail(x$mtime,1)-utils::head(x$mtime,1)
  attr(x,"time_correction_rate") <-
    attr(x,"power_difference") / attr(x,"time_difference")
  x
}


############################## isometric timing ################################

#' Compute timing and magnitude of force in isometric trials
#'
#' Calculate timing and magnitude of force at stimluation, peak force, and
#' various parts of the rising (force development) and relaxation (falling)
#' phases of the twitch.
#'
#' @param x A \code{muscle_stim} object that contains data from an isometric
#'  twitch trial, ideally created via \code{read_ddf}.
#' @param rising Set points of the rising phase to be described.
#'  By default: 10\% and 90\%.
#' @param relaxing Set points of the relaxation phase to be described.
#'  By default: 90\% and 50\%.
#'
#' @details The \code{data.frame} (x) must have time series data organized in
#' columns. Generally, it is preferred that you use a \code{muscle_stim} object
#' imported by \code{read_ddf()}.
#'
#' The \code{rising} and \code{relaxing} arguments allow for the user to supply
#' numeric vectors of any length. By default, these arguments are
#' \code{rising = c(10, 90)} and \code{relaxing  = c(90, 50)}. Numbers in each
#' of these correspond to percent values and capture time and force at that
#' percent of the corresponding curve. These values can be replaced by those
#' that the user specifies and do not necessarily need to have length = 2. But
#' please note that 0 and 100 should not be used, e.g.
#' \code{rising = seq(10, 90, 5)} works, but \code{rising = seq(0, 100, 5)}
#' does not.
#'
#' @return A \code{data.frame} with the following metrics as columns:
#' \item{file_ID }{File ID}
#' \item{time_stim}{Time between beginning of data collection and when
#' stimulation occurs}
#' \item{force_stim}{Magnitude of force at the onset of stimulation}
#' \item{time_peak}{Absolute time of peak force, i.e. time between beginning of
#' data collection and when peak force occurs}
#' \item{force_peak}{Magnitude of peak force}
#' \item{time_rising_X}{Time between beginning of data collection and X\% of
#'  force development}
#' \item{force_rising_X}{Magnitude of force at X\% of force development}
#' \item{time_relaxing_X}{Time between beginning of data collection and X\% of
#'  force relaxation}
#' \item{force_relaxing_X}{Magnitude of force at X\% of relaxation}
#'
#'
#' @references Ahn AN, and Full RJ. 2002. A motor and a brake: two leg extensor
#' muscles acting at the same joint manage energy differently in a running
#' insect. Journal of Experimental Biology 205, 379-389.
#'
#' @author Vikram B. Baliga
#'
#' @examples
#'
#' library(workloopR)
#'
#' # import the twitch.ddf file included in workloopR
#' twitch_dat <-read_ddf(system.file("extdata", "twitch.ddf",
#'                                   package = 'workloopR'))
#'
#' # run isometric_timing() to get info on twitch kinetics
#' # we'll use different set points than the defaults
#' analyze_twitch <- isometric_timing(twitch_dat,
#'                                 rising = c(25, 50, 75),
#'                                 relaxing = c(75, 50, 25))
#'
#' # see the results
#' analyze_twitch
#'
#' @family data analyses
#' @family twitch functions
#'
#' @export
isometric_timing <- function(x,
                             rising = c(10, 90),
                             relaxing = c(90, 50)){
  # check input data
  if(!("isometric" %in% class(x)))
    stop("Please ensure that your data is from an isometric experiment!")
  if("tetanus" %in% class(x))
    relaxing=c()

  # check that set points are numeric between 0 and 100
  if(any(!is.numeric(rising) | rising < 0 | rising > 100))
    stop("Please ensure that all rising set points are numeric values between 0 and 100.")
  if(any(!is.numeric(relaxing) | relaxing < 0 | relaxing > 100))
    stop("Please ensure that all relaxing set points are numeric values between 0 and 100.")

  # convert precents to proportions for easier math
  rising<-rising/100
  relaxing<-relaxing/100

  # find position of peak force and stimulus in dataset
  stim_row<-which.max(x$Stim)
  pf_row<-which.max(x$Force)

  # get force and timing for peak force and stim
  main_results<-data.frame(
    'file_id'=attr(x,"file_id"),
    'time_stim'=x$Time[stim_row],
    'force_stim'=x$Force[stim_row],
    'time_peak'=x$Time[pf_row],
    'force_peak'=x$Force[pf_row],
    stringsAsFactors=FALSE)

  # calculate absolute force at optional set points
  rising_forces<-rising*(x$Force[pf_row]-x$Force[stim_row])+x$Force[stim_row]
  relaxing_forces<-relaxing*(x$Force[pf_row]-x$Force[stim_row])+x$Force[stim_row]

  # calculate corresponding position in dataset
  rising_row<-sapply(rising_forces,function(i)utils::head(which(x$Force>i),1))
  relaxing_row<-sapply(relaxing_forces,function(i)utils::tail(which(x$Force>i),1))

  # extract time and force at these positions, bind together into a vector
  set_point_results<-
    c(unlist(lapply(rising_row, function(i) c(x$Time[i],x$Force[i]))),
      unlist(lapply(relaxing_row, function(i) c(x$Time[i],x$Force[i]))))

  # add names and convert to data.frame
  names(set_point_results)<-
    c(sapply(rising*100, function(i) c(paste0("time_rising_",i),
                                       paste0("force_rising_",i))),
      sapply(relaxing*100, function(i) c(paste0("time_relaxing_",i),
                                         paste0("force_relaxing_",i))))
  set_point_results<-data.frame(as.list(set_point_results))

  # return both result
  cbind(main_results,set_point_results)
}

