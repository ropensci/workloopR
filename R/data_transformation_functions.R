# custom functions
# all written by Vikram B. Baliga (vbaliga@zoology.ubc.ca) and Shreeram
# Senthivasan
# last updated: 2019-07-20


############################### select cycles ###############################

#' Select cycles from a work loop object
#'
#' Retain data from a work loop experiment based on position cycle
#'
#' @param x A \code{workloop} object (see Details for how it should be
#' organized)
#' @param cycle_def A string specifying how cycles should be defined; one of:
#' "lo", "p2p", or "t2t". See Details more info
#' @param keep_cycles The indices of the cycles to keep. Include 0 to keep data
#' identified as being outside complete cycles
#' @param bworth_order Filter order for low-pass filtering of \code{Position}
#'  via \code{signal::butter()} prior to finding L0
#' @param bworth_freq Critical frequency (scalar) for low-pass filtering of
#' \code{Position} via \code{signal::butter()} prior to finding L0
#' @param ... Additional arguments passed to/from other functions that make use
#' of \code{select_cycles()}
#'
#' @details \code{select_cycles()} subsets data from a workloop trial by
#' position cycle. The \code{cycle_def} argument is used to specify which part
#' of the cycle is understood as the beginning and end. There are currently
#' three options: \cr
#' 'lo' for L0-to-L0; \cr
#' 'p2p' for peak-to-peak; and \cr
#' 't2t' for trough-to-trough \cr
#'
#' Peaks are identified using \code{pracma::findpeaks()}. L0 points on the
#' rising edge are found by finding the midpoints between troughs and the
#' following peak. However the first and last extrema and L0 points may be
#' misidentified by this method. Please plot your \code{Position} cycles to
#' ensure the edge cases are identified correctly.
#'
#' The \code{keep_cycles} argument is used to determine which cycles (as
#' defined by \code{cycle_def} should be retained in the final dataset. Zero
#' is the index assigned to all data points that are determined to be outside
#' a complete cycle.
#'
#' The \code{muscle_stim} object (\code{x}) must be a \code{workloop},
#' preferably read in by one of our data import functions. Please see
#' documentation for \code{as_muscle_stim()} if you need to manually construct
#' a \code{muscle_stim} object from another source.
#'
#' @return A \code{workloop} object with rows subsetted by the chosen position
#' cycles. A \code{Cycle} column is appended to denote which cycle each time
#' point is associated with. Finally, all attributes from the input
#' \code{workloop} object are retained and one new attribute is added to
#' record which cycles from the original data were retained.
#'
#' @author Vikram B. Baliga and Shreeram Senthivasan
#'
#' @family data transformations
#' @family workloop functions
#'
#' @examples
#'
#' library(workloopR)
#'
#' # import the workloop.ddf file included in workloopR
#' wl_dat <-read_ddf(system.file("extdata", "workloop.ddf",
#'                               package = 'workloopR'),
#'                   phase_from_peak = TRUE)
#'
#' # select cycles 3 through 5 via the peak-to-peak definition
#' wl_selected <- select_cycles(wl_dat, cycle_def = "p2p", keep_cycles = 3:5)
#'
#' # are only three cycles present, running peak-to-peak?
#' #plot(wl_selected$Position)
#'
#' # are the cycles of (approximately) the same length?
#' summary(as.factor(wl_selected$Cycle))
#'
#' @seealso
#' \code{\link{analyze_workloop}},
#' \code{\link{read_analyze_wl}},
#' \code{\link{read_analyze_wl_dir}}
#'
#' @export
select_cycles <- function(x,
                          cycle_def,
                          keep_cycles = 4:6,
                          bworth_order = 2,
                          bworth_freq = 0.05,
                          ...)
{
  if(!any(class(x) == "workloop"))
    stop("Input data should be of class `workloop`")
  if(!is.numeric(keep_cycles))
    stop("keep_cycles should be numeric")
  if(missing(cycle_def)){
    cycle_def<-"lo"
    warning("Cycle definition not supplied! Defaulting to L0-to-L0")
  }
  if(is.na(attr(x,"cycle_frequency")))
    stop("Length-out cycle frequency is needed to identify cycles! Please set the `cycle_frequency` attribute accordingly.")

  # get cycle frequency and sample frequency
  cyc_freq<-attr(x,"cycle_frequency")
  samp_freq<-attr(x,"sample_frequency")

  # Use butterworth-filtered position data to identify peaks
  bworth<-signal::butter(bworth_order,bworth_freq)
  smPos<-signal::filtfilt(bworth,x$Position)

  # Calculate minimum number of ups before a peak from proportion ups
  qf<-floor(0.25*(1/cyc_freq)*samp_freq)-1
  peaks<-stats::setNames(data.frame(pracma::findpeaks(smPos,nups=qf)[,2:4]),
                         c("peak","start","end"))

  switch(cycle_def,
    # L0-to-Lo assumes position cycle starts and ends on an L0
    # Most L0 are found by averaging indices of a peak and previous trough
    "lo"={splits<-round((peaks$start+peaks$peak)/2)
    # The first L0 is the lowest point before first peak
          splits[1]<-peaks$start[1]
    # The last L0 is the last peak
          splits[length(splits)]<-peaks$peak[nrow(peaks)]
          splits<-c(0,splits,nrow(x))},
    "p2p"=splits<-c(0,peaks$peak,nrow(x)),
    "t2t"=splits<-c(0,peaks$start,utils::tail(peaks$end,1),nrow(x)),
    stop("Invalid cycle definition! Please select one of:\n
      'lo':  L0-to-L0
      'p2p': peak-to-peak
      't2t': trough-to-trough")
  )
  splits<-(splits-c(NA,utils::head(splits,-1)))[-1]

  cycle<-unlist(sapply(seq_along(splits), function(i)rep(i-1,splits[i])))
  x$Cycle<-replace(cycle,cycle==max(cycle),0)

  # Update cycle definition and total cycles
  attr(x,"cycle_def")<-cycle_def
  attr(x,"total_cycles")<-max(x$Cycle)

  # Subset by keep_cycles and rename cycles by letters
  if(any(keep_cycles<0 | keep_cycles>max(x$Cycle)))
    warning("The keep_cycles argument includes cycles that don't exist (negative or greater than total_cycles). These are being ignored.")
  x<-x[x$Cycle %in% keep_cycles,]
  x$Cycle<-letters[as.factor(x$Cycle)]
  if(!all(is.na(attr(x,"units")))) attr(x,"units")<-c(attr(x,"units"),"letters")
  attr(x,"retained_cycles")<-keep_cycles
  x
}

############################## position inversion ##############################

#' Invert the position data
#'
#' Multiply instantaneous position by -1.
#'
#' @param x A \code{muscle_stim} object
#'
#' @details The \code{muscle_stim} object can be of any type, including
#' \code{workloop}, \code{twitch}, or \code{tetanus}.
#'
#' If you have manually constructed the object via \code{as_muscle_stim()},
#' the \code{muscle_stim} object should have a column entitled \code{Position}.
#' Other columns and attributes are welcome and will be passed along unchanged.
#'
#' @return A \code{workloop} object with inverted position. The
#' \code{position_inverted} attribute is set to \code{TRUE} and all others are
#' retained.
#'
#' @author Vikram B. Baliga
#'
#' @family data transformations
#' @family workloop functions
#' @family twitch functions
#' @family tetanus functions
#'
#' @examples
#'
#' library(workloopR)
#'
#' # import the workloop.ddf file included in workloopR
#' wl_dat <-read_ddf(system.file("extdata", "workloop.ddf",
#'                               package = 'workloopR'),
#'                   phase_from_peak = TRUE)
#'
#' # invert the sign of Position
#' wl_fixed <- invert_position(wl_dat)
#'
#' # quick check:
#' max(wl_fixed$Position)/min(wl_dat$Position) # -1
#'
#' @export
invert_position <- function(x)
{
  if(!any(class(x) == "muscle_stim"))
    stop("Input data should be of class `muscle_stim`")
  x$Position<-x$Position*-1
  attr(x,"position_inverted")<-TRUE
  return(x)
  }


############################## gear ratio correction ###########################

#' Adjust for the gear ratio of a motor arm
#'
#' Fix a discrepancy between the gear ratio of the motor arm used and the gear
#' ratio recorded by software.
#'
#' @param x A \code{muscle_stim} object
#' @param GR Gear ratio, set to 1 by default
#'
#' @details The \code{muscle_stim} object can be of any type, including
#' \code{workloop}, \code{twitch}, or \code{tetanus}.
#'
#' If you have manually constructed the object via \code{as_muscle_stim()},
#' the \code{muscle_stim} object should have columns as follows: \cr
#' \code{Position}: length change of the muscle; \cr
#' \code{Force}: force \cr
#'
#' @return An object of the same class(es) as the input (\code{x}). The function
#'  will multiply \code{Position} by (1/GR) and multiply \code{Force} by GR,
#'  returning an object with new values in \code{$Position} and \code{$Force}.
#'  Other columns and attributes are welcome and will simply be passed on
#'  unchanged into the resulting object.
#'
#' @author Vikram B. Baliga
#'
#' @family data transformations
#' @family workloop functions
#' @family twitch functions
#' @family tetanus functions
#'
#' @examples
#'
#' library(workloopR)
#'
#' # import the workloop.ddf file included in workloopR
#' wl_dat <-read_ddf(system.file("extdata", "workloop.ddf",
#'                               package = 'workloopR'),
#'                   phase_from_peak = TRUE)
#'
#' # apply a gear ratio correction of 2
#' # this will multiply Force by 2 and divide Position by 2
#' wl_fixed <- fix_GR(wl_dat, GR = 2)
#'
#' # quick check:
#' max(wl_fixed$Force)/max(wl_dat$Force)       #5592.578 / 2796.289 = 2
#' max(wl_fixed$Position)/max(wl_dat$Position) #1.832262 / 3.664524 = 0.5
#'
#' @seealso
#' \code{\link{analyze_workloop}},
#' \code{\link{read_analyze_wl}},
#' \code{\link{read_analyze_wl_dir}}
#'
#' @export
fix_GR <- function(x,
                   GR = 1)
{
  # Check that x is correct type of object
  if (!any(class(x) == "muscle_stim"))
    stop("Input data should be of class `muscle_stim`")

  # check that gear ratio is numeric
  if (!is.numeric(GR))
  {
    stop('Gear ratio (GR) must be numeric')
  }

  x$Position<-x$Position*(1/GR)

  x$Force<-x$Force*GR

  attr(x,"gear_ratio")<-attr(x,"gear_ratio")*GR
  if("workloop" %in% class(x))
    if(!is.na(attr(x,"amplitude")))
      attr(x,"amplitude")<-attr(x,"amplitude")*(1/GR)
  return(x)
}
