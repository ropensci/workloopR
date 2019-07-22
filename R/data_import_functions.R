# custom functions
# all written by Vikram B. Baliga (vbaliga@zoology.ubc.ca) and Shreeram
# Senthivasan
# last updated: 2019-07-20


########################## workloop object functions ###########################
# Object constructor
# See read_ddf for details

# Top level class for all objects created by read_ddf
#' @noRd
muscle_stim <-
  function(data,
           units,
           sample_frequency,
           header,
           units_table,
           protocol_table,
           stim_table,
           file_id,
           mtime,
           ...){
    attr(data, "units") <- units
    attr(data, "sample_frequency") <- sample_frequency
    attr(data, "header") <- header
    attr(data, "units_table") <- units_table
    attr(data, "protocol_table") <- protocol_table
    attr(data, "stim_table") <- stim_table
    attr(data, "stimulus_pulses") <- stim_table$pulses[1]
    attr(data, "stimulus_offset") <- stim_table$offset[1]
    attr(data, "stimulus_width") <- stim_table$width[1]
    attr(data, "gear_ratio") <- 1
    attr(data, "file_id") <- file_id
    attr(data, "mtime") <- mtime
    class(data) <- c(class(data), "muscle_stim", "data.frame")
    data
  }

# Classes for specific trial type
#' @noRd
workloop <-
  function(data,
           stim_table,
           cycle_table,
           sample_frequency,
           phase_from_peak,
           ...) {
    attr(data, "stimulus_frequency") <- stim_table$frequency[1]
    attr(data, "cycle_frequency") <- cycle_table$frequency[1]
    attr(data, "total_cycles") <- cycle_table$total_cycles[1]
    attr(data, "cycle_def") <- "lo"
    attr(data, "amplitude") <- cycle_table$amplitude[1]

    # Calculate Phase
    phase<-(which.max(data$Stim)-which.max(data$Position))/sample_frequency*stim_table$cycle_frequency[1]
    if(!phase_from_peak)
      phase<-phase+0.25
    # convert 0-1 scale to -50 to +50
    attr(data,"phase")<-(((phase+0.5)%%1)-0.5)*100

    attr(data,"position_inverted")<-FALSE
    class(data)<-c("workloop")
    muscle_stim(data=data,
                stim_table=stim_table,
                sample_frequency=sample_frequency,
                ...)
  }

#' @noRd
tetanus<-
  function(data,stim_table,...){
    attr(data,"stimulus_frequency")<-stim_table$frequency[1]
    attr(data,"stimulus_length")<-stim_table$length[1]
    class(data)<-c("tetanus","isometric")
    muscle_stim(data=data,
                stim_table=stim_table,
                ...)
  }

#' @noRd
twitch<-
  function(data,...){
    class(data)<-c("twitch","isometric")
    muscle_stim(data=data,...)
  }

# Generate type-specific output header for workloop objects
#' @noRd
print_muscle_stim_header<-function(x,include_time=TRUE){
  type<-paste(toupper(substring(class(x)[1],1,1)),
              substring(class(x)[1],2),sep="")
  if(include_time)
    cat(paste0("# ",type," Data: ",
               ncol(x)-1," channels recorded over ",
               nrow(x)/attr(x,"sample_frequency"),"s\n"))
  else
    cat(paste0("# ",type," Data:\n\n"))
}

# Print method

#' @export
print.muscle_stim <- function(x, n = 6, ...){
  print_muscle_stim_header(x)
  cat(paste0("File ID: ",attr(x,"file_id"),"\n\n"))
  class(x)<-"data.frame"
  print(utils::head(x,n=n))
  if(n < nrow(x))
    cat(paste0("# \u2026 with ",nrow(x)-n," more rows\n"))
}

#' @export
print.analyzed_workloop <- function(x, n = 6, ...){
  cat(paste0("File ID: ",
             attr(x,"file_id")))
  cat(paste0("\nCycles: ",
             length(attr(x,"retained_cycles")),
             " cycles kept out of ",
             attr(x,"total_cycles")))
  cat(paste0("\nMean Work: ",
             round(mean(attr(x,"summary")$Work),5),
             " J"))
  cat(paste0("\nMean Power: ",
             round(mean(attr(x,"summary")$Net_Power),5),
             " W\n\n"))
}

# Summary method

#' @export
summary.muscle_stim <- function(object, ...){
  print_muscle_stim_header(object,...)
  cat(paste0("\nFile ID: ",attr(object,"file_id")))
  cat(paste0("\nMod Time (mtime): ",attr(object,"mtime")))
  cat(paste0("\nSample Frequency: ",attr(object,"sample_frequency"),"Hz\n\n"))
  cat(paste0("data.frame Columns: \n"))
  for(i in 2:ncol(object))
    cat(paste0("  ",colnames(object)[i]," (",attr(object,"units")[i],")\n"))
  cat(paste0("\nStimulus Offset: ",attr(object,"stimulus_offset"),"s\n"))
  cat(paste0("Stimulus Frequency: ",attr(object,"stimulus_frequency"),"Hz\n"))
  cat(paste0("Stimulus Width: ",attr(object,"stimulus_width"),"ms\n"))
  cat(paste0("Stimulus Pulses: ",attr(object,"stimulus_pulses"),"\n"))
  cat(paste0("Gear Ratio: ",attr(object,"gear_ratio"),"\n"))
}

#' @export
summary.workloop <- function(object, ...){
  NextMethod()
  cat(paste0("\nCycle Frequency: ",attr(object,"cycle_frequency"),"Hz\n"))
  cat(paste0("Total Cycles (",
             switch(attr(object,"cycle_def"),
                    "lo"="L0-to-L0",
                    "p2p"="peak-to-peak",
                    "t2t"="trough-to-trough",
                    "undefined"),
             "): ",
             attr(object,"total_cycles"),"\n"))
  if(!is.null(attr(object,"retained_cycles")))
     cat(paste0("Cycles Retained: ",
                length(attr(object,"retained_cycles")),
                "\n"))
  cat(paste0("Amplitude: ",
             attr(object,"amplitude"),
             attr(object,"units")[grep("Position",colnames(object))],"\n\n"))
  if(attr(object,"position_inverted"))
    cat("\nPlease note that Position is inverted!\n\n")
}

#' @export
summary.tetanus <- function(object, ...){
  NextMethod()
  cat(paste0("Stimulus Length: ",attr(object,"stimulus_length"),"s\n\n"))
}

#' @export
summary.analyzed_workloop <- function(object, ...){
  summary(object[[1]],include_time=FALSE)
  cat("\n")
  print(attr(object,"summary"))
}

######################### read non-ddf work loop files #########################

#' Create your own muscle_stim object
#'
#' For use when data are not stored in .ddf format and you would like
#' to create a \code{muscle_stim} object that can be used by other workloopR
#' functions.
#'
#' @param x A \code{data.frame}. See Details for how it should be organized.
#' @param type Experiment type; must be one of: "workloop", "tetanus", or
#' "twitch."
#' @param sample_frequency Numeric value of the frequency at which samples were
#' recorded; must be in Hz. Please format as numeric, e.g. \code{10000} works
#' but \code{10000 Hz} does not
#' @param ... Additional arguments that can be passed in as attributes. See
#' Details.
#'
#' @details \code{muscle_stim} objects, which are required by (nearly) all
#' workloopR functions, are automatically created via \code{read_ddf()}. Should
#' you have data that are stored in a format other than .ddf, use this function
#' to create your own object of class \code{muscle_stim}.
#'
#' The input \code{x} must be a \code{data.frame} that contains time series
#' of numeric data collected from an experiment. Each row must correspond to a
#' sample, and these columns (exact title matches) must be included: \cr
#' "Time" - time, recorded in seconds \cr
#' "Position" - instantaneous position of the muscle, preferably in millimeters \cr
#' "Force" - force, preferably in millinewtons \cr
#' "Stim" - whether stimulation has occurred. All entries must be either 0 (no
#' stimulus) or 1 (stimulus ocurrence).
#'
#' Additional arguments can be provided via \code{...}. For all experiment
#' types, the following attributes are appropriate: \cr
#' "units","header", "units_table",
#' "protocol_table", "stim_table",
#' "stimulus_pulses", "stimulus_offset",
#' "stimulus_width", "gear_ratio",
#' "file_id", or "mtime".
#'
#' Please ensure that further attributes are appropriate to your experiment
#' type.
#'
#' For workloops, these include:
#' "stimulus_frequency", "cycle_frequency",
#' "total_cycles", "cycle_def",
#' "amplitude", "phase",
#' and "position_inverted"
#'
#' For twitches or tetanic trials:
#' "stimulus_frequency", and "stimulus_length"
#'
#'
#' @inherit read_ddf return
#'
#' @examples
#'
#' library(workloopR)
#'
#' # import the workloop.ddf file included in workloopR
#' wl_dat <-read_ddf(system.file("extdata", "workloop.ddf", package = 'workloopR'))
#'
#' # see how this object is organized - this will give you a sense
#' # of how your inputs to `as_muscle_stim()` should be arranged:
#' #head(wl_dat)
#' #str(wl_dat)
#' # formatting of attributes:
#' #names(attributes(wl_dat))
#' #str(attributes(wl_dat))
#'
#' @author Shreeram Senthivasan
#'
#' @family data import functions
#'
#' @seealso
#' \code{\link{read_ddf}}
#'
#' @export
as_muscle_stim <- function(x,
                           type,
                           sample_frequency,
                           ...){
  # Check for missing information
  if(missing(type))stop("Please specify the experiment type! The type argument should be one of: workloop, tetanus, or twitch.")
  if(!(type %in% c("workloop","tetanus","twitch"))|length(type)!=1)
    stop("Invalid experiment type! The type argument should be one of: workloop, tetanus, or twitch.")
  if(!all(c("Position","Force","Stim") %in% names(x)))
    stop("Couldn't find one or more of the following necessary columns: Position, Force, Stim. Please ensure that the columns match the naming conventions.")
  if(missing(sample_frequency)&!("Time" %in% names(x)))
    stop("Insufficient information to infer the sampling frequency. Please provide a value for the sample_frequency argument or include a column named `Time` in the dataframe.")

  # Consolidate time / sample frequency information
  if(!missing(sample_frequency))
    x$Time=(1:nrow(x)-1)/sample_frequency
  else
    sample_frequency<-1/(x$Time[2]-x$Time[1])

  # Generate a list of acceptable attributes given experiment type
  valid_args <-
    c(
      "units",
      "header",
      "units_table",
      "protocol_table",
      "stim_table",
      "stimulus_pulses",
      "stimulus_offset",
      "stimulus_width",
      "gear_ratio",
      "file_id",
      "mtime"
    )
  switch(
    type,
    "workloop" = valid_args <-
      c(
        valid_args,
        "stimulus_frequency",
        "cycle_frequency",
        "total_cycles",
        "cycle_def",
        "amplitude",
        "phase",
        "position_inverted"
      ),
    "tetanus" = valid_args <-
      c(valid_args, "stimulus_frequency", "stimulus_length")
  )

  # Check for invalid attributes and assign valids
  args<-list(...)
  if(!all(names(args) %in% valid_args))
    warning("One or more provided attributes do not match known attributes. These will attributes will not be assigned.")
  for(i in intersect(names(args),valid_args))
    attr(x,i)<-args[[i]]
  for(i in setdiff(valid_args,names(args)))
    attr(x,i)<-NA
  attr(x,"sample_frequency")<-sample_frequency
  if(is.na(attr(x,"gear_ratio"))) attr(x,"gear_ratio")<-1
  if(type=="workloop")
    if(is.na(attr(x,"position_inverted"))) attr(x,"position_inverted")<-FALSE

  # Assign classes and return
  class(x)<-c("muscle_stim","data.frame")
  switch(type,
         "workloop"=class(x)<-c("workloop",class(x)),
         "tetanus"=class(x)<-c("tetanus","isometric",class(x)),
         "twitch"=class(x)<-c("twitch","isometric",class(x)))
  x
}

########################## read_ddf files - work loops #########################

#' Import work loop or isometric data from .ddf files
#'
#' \code{read_ddf} reads in workloop, twitch, or tetanus experiment data from
#' .ddf files.
#'
#' @param filename A .ddf file that contains data from a single workloop,
#' twitch, or tetanus experiment
#' @param file_id A string identifying the experiment. The filename is used by
#' default.
#' @param rename_cols List consisting of a vector of indices of columns to
#' rename and a vector of new column names. See Details.
#' @param skip_cols Numeric vector of column indices to skip. See Details.
#' @param phase_from_peak Logical, indicating whether percent phase of
#' stimulation should be recorded relative to peak length or relative to L0
#' (default)
#' @param ... Additional arguments passed to/from other functions that work
#' with \code{read_ddf()}
#'
#' @details Read in a .ddf file that contains data from an experiment. If
#' position and force do not correspond to columns 2 and 3 (respectively),
#' replace "2" and "3" within \code{rename_cols} accordingly. Similarly,
#' \code{skip_cols = 4:11} should be adjusted if more than 11 columns are
#' present and/or columns 4:11 contain important data.
#'
#' Please note that there is no correction for gear ratio or further
#' manipulation of data. See \code{fix_GR} to adjust gear ratio. Gear ratio can
#' also be adjusted prior to analyses within the \code{analyze_workloop()}
#' function, the data import all-in-one function \code{read_analyze_wl()}, or
#' the batch analysis all-in-one \code{read_analyze_wl_dir()}.
#'
#' Please also note that organization of data within the .ddf file is assumed to
#' conform to that used by Aurora Scientific's Dynamic Muscle Control and
#' Analysis Software. YMMV if using a .ddf file from another source. The
#' \code{as_muscle_stim()} function can be used to generate \code{muscle_stim}
#' objects if data are imported via another function. Please feel free to
#' contact us with any issues or requests.
#'
#'
#' @return An object of class \code{workloop}, \code{twitch}, or \code{tetanus},
#' all of which inherit class \code{muscle_stim}. These objects behave like
#' \code{data.frames} in most situations but also store metadata from the ddf
#' as attributes.
#'
#' The \code{muscle_stim} object's columns contain:
#' \item{Time}{Time}
#' \item{Position}{Length change of the muscle, uncorrected for gear ratio}
#' \item{Force}{Force, uncorrected for gear ratio}
#' \item{Stim}{When stimulation occurs, on a binary scale}
#'
#' In addtion, the following information is stored in the \code{data.frame}'s
#' attributes:
#' \item{sample_frequency}{Frequency at which samples were collected}
#' \item{pulses}{Number of sequential pulses within a stimulation train}
#' \item{total_cycles_lo}{Total number of oscillatory cycles (assuming sine
#' wave trajectory) that the muscle experienced. Cycles are defined with respect
#' to initial muscle length (L0-to-L0 as opposed to peak-to-peak).}
#' \item{amplitude}{amplitude of length change (again, assuming sine wave
#' trajectory)}
#' \item{cycle_frequency}{Frequency of oscillations (again, assuming sine wave
#' trajectory)}
#' \item{units}{The units of measurement for each column in the
#' \code{data.frame}. This might be the most important attribute so please check
#'  that it makes sense!}
#'
#' @author Vikram B. Baliga and Shreeram Senthivasan
#'
#'
#' @family data import functions
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
#' # or import your own file
#' #my_dat <- read_ddf("./my/file/path/myfile.ddf")
#'
#' @export
read_ddf <-
  function(filename,
           file_id = NA,
           rename_cols = list(c(2, 3), c("Position", "Force")),
           skip_cols = 4:11,
           phase_from_peak = FALSE,
           ...)
  {
    # Import and checks
    if(missing(filename)) stop("A filename is required")
    if(!file.exists(filename)) stop(paste0("File ",filename," not found!"))
    f<-file(filename,"r")
    if(!grepl("DMC.*Data File",readLines(f,1))){
      close(f)
      stop("The input file does not appear to be a DMC Datafile (ddf)")
    }
    if(is.na(file_id)) file_id<-basename(filename)

    # get metadata
    mtime<-file.info(filename)$mtime

    # Setup for reading in file
    header<-c()
    units_table<-c()
    protocol_table<-c()

    # Read in Header
    while(!grepl("Calibration Data",(l<-readLines(f,1))))
      header<-c(header,l)
    sample_frequency<-as.numeric(sub(".*: ","",header[1]))

    # Read in Calibration Table
    while(!grepl("Comments",(l<-readLines(f,1))))
      units_table<-c(units_table,l)
    units_table<-t(utils::read.table(text=units_table,
                               row.names=1,
                               sep="\t",
                               stringsAsFactors=FALSE))
    rownames(units_table)<-c()
    colnames(units_table)<-sub(" .*","",colnames(units_table))
    units_table<-data.frame(units_table,stringsAsFactors=FALSE)
    units_table[3:5]<-lapply(units_table[3:5],as.numeric)
    units<-c("s",units_table$Units[-skip_cols+1],"TTL")
    if(!all(units %in% c("s","mm","mN","TTL")))
      warning("Non-standard units detected in ddf file! Please note that calculations currently assume raw data are in seconds, millimeters, and millinewtons.")

    # Read in Protocol Array
    while(!grepl("Protocol",readLines(f,1))) {}
    readLines(f,1) # Discard empty line
    while((l<-readLines(f,1))!="")
      protocol_table<-c(protocol_table,l)
    protocol_table<-utils::read.table(text=protocol_table,
                                 sep="\t",
                                 stringsAsFactors=FALSE,
                                 col.names=c("Wait.s",
                                             "Then.action",
                                             "On.port",
                                             "Units",
                                             "Parameters"))

    # Read in data
    while(!grepl("Test Data",(l<-readLines(f,1)))){}
    readLines(f,1)
    dataz<-utils::read.table(text=readLines(f),
                             header=TRUE,
                             sep="\t",
                             stringsAsFactors=FALSE)
    if(any(!apply(dataz,2,is.numeric))) warning("The ddf file includes non-numeric data. Please ensure that this is intentional before proceeding.")
    close(f)

    # Parse file type
    read_filetype.ddf<-NULL
    switch(
           grep("Stim",protocol_table[[2]],value=TRUE)[1],
           "Stimulus-Train"=read_filetype.ddf<-read_wl.ddf,
           "Stimulus-Twitch"=read_filetype.ddf<-read_twitch.ddf,
           "Stimulus-Tetanus"=read_filetype.ddf<-read_tetanus.ddf,
            stop("Could not parse experiment type (workloop, twitch, or tetanus)! Please ensure that the protocol section of the ddf header includes a label with one of the following: Stimulus-Train, Stimulus-Twitch, or Stimulus-Tetanus.")
    )
    read_filetype.ddf(file_id=file_id,
                     mtime=mtime,
                     header=header,
                     units_table=units_table,
                     units=units,
                     protocol_table=protocol_table,
                     raw_data=dataz,
                     sample_frequency=sample_frequency,
                     rename_cols=rename_cols,
                     skip_cols=skip_cols,
                     phase_from_peak=phase_from_peak)
  }

#' Import a batch of work loop or isometric data files from a directory
#'
#' Uses \code{read_ddf()} to read in workloop, twitch, or tetanus experiment
#' data from multiple .ddf files.
#'
#' @param filepath Path where files are stored. Should be in the same folder.
#' @param pattern Regex pattern for identifying relevant files in the filepath.
#' @param sort_by Metadata by which files should be sorted to be in the correct
#' run order. Defaults to \code{mtime}, which is time of last modification of
#' files.
#' @param ... Additional arguments to be passed to \code{read_ddf()}.
#'
#' @inherit read_ddf details
#'
#' @return A list of objects of class \code{workloop}, \code{twitch}, or
#' \code{tetanus}, all of which inherit class \code{muscle_stim}. These objects
#' behave like \code{data.frames} in most situations but also store metadata
#' from the ddf as attributes.
#'
#' Each \code{muscle_stim} object's columns contain:
#' \item{Time}{Time}
#' \item{Position}{Length change of the muscle, uncorrected for gear ratio}
#' \item{Force}{Force, uncorrected for gear ratio}
#' \item{Stim}{When stimulation occurs, on a binary scale}
#'
#' In addtion, the following information is stored in each \code{data.frame}'s
#' attributes:
#' \item{sample_frequency}{Frequency at which samples were collected}
#' \item{pulses}{Number of sequential pulses within a stimulation train}
#' \item{total_cycles_lo}{Total number of oscillatory cycles (assuming sine
#' wave trajectory) that the muscle experienced. Cycles are defined with respect
#' to initial muscle length (L0-to-L0 as opposed to peak-to-peak).}
#' \item{amplitude}{amplitude of length change (again, assuming sine wave
#' trajectory)}
#' \item{cycle_frequency}{Frequency of oscillations (again, assuming sine wave
#' trajectory)}
#' \item{units}{The units of measurement for each column in the
#' \code{data.frame}. This might be the most important attribute so please check
#'  that it makes sense!}
#'
#' @author Vikram B. Baliga and Shreeram Senthivasan
#'
#'
#' @family data import functions
#'
#' @examples
#'
#' library(workloopR)
#'
#' # import a set of twitch .ddf files included in workloopR
#' workloop_dat <-read_ddf_dir(system.file("extdata/wl_duration_trials",
#'                  package = 'workloopR'))
#'
#' # or import your own file
#' #my_dat <- read_ddf_dir("./my/file/path/")
#'
#' @export
read_ddf_dir <- function(filepath,
                         pattern = "*.ddf",
                         sort_by = "mtime",
                         ...){
  # Generate list of filenames
  filename_list<-list.files(path=filepath,pattern=pattern,full.names=TRUE)
  if(length(filename_list)==0) stop("No files matching the pattern found at the given directory!")

  # Generate list of muscle_stim objects
  ms_list<-lapply(filename_list,function(i) read_ddf(i,...))

  # Sort list, likely by modification time
  if(is.null(attr(ms_list[[1]],sort_by))){
    warning("The provided sort_by argument is not a valid attribute. Defaulting to `mtime`.")
    sort_by<-"mtime"
  }
  ms_list<-ms_list[order(sapply(ms_list,function(i)attr(i,sort_by)))]

  ms_list
}

############################# rescale data matrix ##############################
# Rescales data in ddf files using the scale and offset parameters
#' @noRd
rescale_data<-
  function(dataz,
           unitz,
           sample_frequency,
           rename_cols,
           skip_cols){
    rescaled<-mapply(function(raw,offset,scale){
                       if(!is.numeric(raw))
                         return(raw)
                       else
                         (raw+offset)*scale},
                     dataz[unitz$Channel],
                     unitz$Offset,
                     unitz$Scale)

    rescaled<-data.frame(Time=(1:nrow(dataz))/sample_frequency,
                  rescaled,
                  Stim=dataz$Stim)

    #rename columns, if desired
    if(!is.null(rename_cols))
      names(rescaled)[rename_cols[[1]]]<-rename_cols[[2]]

    rescaled[,-skip_cols]
  }

########################## read_ddf files - workloop ###########################
#' @noRd
read_wl.ddf<-
  function(raw_data,
           units_table,
           protocol_table,
           sample_frequency,
           rename_cols,
           skip_cols,
           ...)
  {
    #get info on experimental parameters
    stim_table<-
      utils::read.table(
        text=protocol_table[grepl("Stim",protocol_table$Then.action),"Units"],
        sep=",",
        col.names=c("offset","frequency","width","pulses","cycle_frequency")
      )
    cycle_table<-
      utils::read.table(
        text=protocol_table[grepl("Sine",protocol_table$Then.action),"Units"],
        sep=",",
        col.names=c("frequency","amplitude","total_cycles")
      )

    #use scale (and maybe offset) to convert Volts into units
    rescaled_data<-rescale_data(raw_data,
                                units_table,
                                sample_frequency,
                                rename_cols,
                                skip_cols)

    #construct and return workloop object
    workloop(data=rescaled_data,
             sample_frequency=sample_frequency,
             units_table=units_table,
             protocol_table=protocol_table,
             stim_table=stim_table,
             cycle_table=cycle_table,
             ...)
  }


############################ read_ddf files - twitch ###########################
#' @noRd
read_twitch.ddf<-
  function(raw_data,
           units_table,
           protocol_table,
           sample_frequency,
           rename_cols=list(c(2,3),c("Position","Force")),
           skip_cols=4:11,
           ...)
  {
    #get info on experimental parameters
    stim_table<-
      utils::read.table(
        text=protocol_table[grepl("Stim",protocol_table$Then.action),"Units"],
        sep=",",
        col.names=c("offset","width")
      )
    stim_table$pulses<-rep(1,nrow(stim_table))

    #use scale (and maybe offset) to convert Volts into units
    rescaled_data<-rescale_data(raw_data,
                                units_table,
                                sample_frequency,
                                rename_cols,
                                skip_cols)

    #construct and return workloop object
    twitch(data=rescaled_data,
           sample_frequency=sample_frequency,
           units_table=units_table,
           protocol_table=protocol_table,
           stim_table=stim_table,
           ...)
  }


########################## read_ddf files - tetanus ##########################
#' @noRd
read_tetanus.ddf<-
  function(raw_data,
           units_table,
           protocol_table,
           sample_frequency,
           rename_cols=list(c(2,3),c("Position","Force")),
           skip_cols=4:11,
           ...)
  {
    #get info on experimental parameters
    stim_table <-
      utils::read.table(
        text = protocol_table[grepl("Stim",protocol_table$Then.action),"Units"],
        sep=",",
        col.names=c("offset","frequency","width","length")
      )
    stim_table$pulses<-
      as.integer(floor(stim_table$frequency * stim_table$length))

    #use scale (and maybe offset) to convert Volts into units
    rescaled_data<-rescale_data(raw_data,
                                units_table,
                                sample_frequency,
                                rename_cols,
                                skip_cols)

    #construct and return workloop object
    tetanus(data=rescaled_data,
            sample_frequency=sample_frequency,
            units_table=units_table,
            protocol_table=protocol_table,
            stim_table=stim_table,
            ...)
  }


###################### work loop reading and data extraction ###################

#' All-in-one import function for work loop files
#'
#' \code{read_analyze_wl()} is an all-in-one function to read in a work loop
#' file, select cycles, and compute work and power output.
#'
#' @param filename A .ddf file that contains data from a
#' single workloop experiment
#' @param ... Additional arguments to be passed to \code{read_ddf()},
#' \code{select_cycles()},
#' or \code{analyze_workloop()}.
#'
#' @details Please be careful with units! See Warnings below. This function
#' combines \code{read_ddf()} with \code{select_cycles()} and then ultimately
#' \code{analyze_workloop()} into one handy function.
#'
#' As detailed in these three functions, possible arguments include: \cr
#' \code{cycle_def} - used to specify which part of the cycle is understood as
#' the beginning and end. There are currently three options: 'lo' for L0-to-L0;
#' 'p2p' for peak-to-peak; and 't2t' for trough-to-trough \cr
#' \code{bworth_order} - Filter order for low-pass filtering of \code{Position}
#'  via \code{signal::butter} prior to finding peak lengths. Default: 2. \cr
#' \code{bworth_freq} - Critical frequency (scalar) for low-pass filtering of
#' \code{Position} via \code{signal::butter} prior to finding peak lengths.
#' Default: 0.05. \cr
#' \code{keep_cycles} - Which cycles should be retained. Default: 4:6. \cr
#' \code{GR} - Gear ratio. Default: 1. \cr
#' \code{M} - Velocity multiplier used to positivize velocity; should be either
#' -1 or 1. Default: -1. \cr
#' \code{vel_bf} - Critical frequency (scalar) for low-pass filtering of
#' velocity via \code{signal::butter}. Default: 0.05. \cr
#'
#' The gear ratio (GR) and velocity multiplier (M) parameters can help correct
#' for issues related to the magnitude and sign of data collection. By
#' default, they are set to apply no gear ratio adjustment and to positivize
#' velocity. Instanteous velocity is often noisy and the \code{vel_bf} parameter
#' allows for low-pass filtering of velocity data. See \code{signal::butter()}
#' and \code{signal::filtfilt()} for details of how filtering is achieved.
#'
#' @inherit analyze_workloop return
#' @inheritSection analyze_workloop Warning
#'
#' @references Josephson RK. 1985. Mechanical Power output from Striated Muscle
#'  during Cyclic Contraction. Journal of Experimental Biology 114: 493-512.
#'
#' @author Vikram B. Baliga
#'
#' @family data analyses
#' @family data import functions
#' @family workloop functions
#'
#' @examples
#'
#' library(workloopR)
#'
#' # import the workloop.ddf file included in workloopR and analyze with
#' # a gear ratio correction of 2 and cycle definition of peak-to-peak
#' wl_dat <- read_analyze_wl(system.file("extdata", "workloop.ddf",
#'                                       package = 'workloopR'),
#'                           phase_from_peak = TRUE,
#'                           GR = 2, cycle_def = "p2p")
#'
#'
#' @seealso
#' \code{\link{read_ddf}},
#' \code{\link{select_cycles}}
#' \code{\link{analyze_workloop}}
#'
#' @export
read_analyze_wl <- function(filename,
                            ...){
  valid_args<-c("file_id","rename_cols","skip_cols","phase_from_peak","cycle_def","keep_cycles","bworth_order","bworth_freq","simplify","GR","M","vel_bf")
  arg_names<-names(list(...))
  if(!all(arg_names %in% valid_args)) warning("One or more provided attributes do not match known attributes. These will attributes will not be assigned.")

  fulldata<-read_ddf(filename,...)
  if(!("workloop" %in% class(fulldata)))
    stop(paste0("The provided file ",filename," does not appear to contain data from a workloop experiment!"))
  analyze_workloop(select_cycles(fulldata,...),...)
}


###################### file info for sequence of work loops ####################

#' Get file info for a sequence of experiment files
#'
#' Grab metadata from files stored in the same folder (e.g. a sequence of trials
#'  in an experiment).
#'
#' @param filepath Path where files are stored. Should be in the same folder.
#' @param pattern Regex pattern for identifying relevant files in the filepath.
#'
#' @details If several files (e.g. successive trials from one experiment) are
#' stored in one folder, use this function to obtain metadata in a list
#' format. Runs \code{file.info} from base R to extract info from files.
#'
#' This function is not truly considered to be part of the batch analysis
#' pipeline;
#' see \code{read_analyze_wl_dir()} for a similar function that not
#' only grabs metadata but also imports & analyzes files. Instead,
#' \code{get_wl_metadata()} is meant to be a handy function to investigate
#' metadata issues that arise if running \code{read_analyze_wl_dir()} goes awry.
#'
#' Unlike \code{read_analyze_wl_dir()}, this function does not necessarily need
#' files to all be work loops. Any filetype is welcome (as long as the Regex
#' \code{pattern} argument makes sense).
#'
#' @family data import functions
#' @family workloop functions
#' @family batch analyses
#'
#' @seealso
#' \code{\link{summarize_wl_trials}}
#'
#' @author Vikram B. Baliga
#'
#' @examples
#'
#' library(workloopR)
#'
#' # get file info for files included with workloopR
#' wl_meta <- get_wl_metadata(system.file("extdata/wl_duration_trials",
#'                                        package = 'workloopR'))
#'
#' # or on your own directory
#' #my_meta <- get_wl_metadata("./my/file/path/")
#'
#' @export
get_wl_metadata <- function(filepath,
                            pattern = "*.ddf"){
  exp_list<-file.info(list.files(path=filepath,pattern=pattern,
                                 full.names=TRUE,recursive=TRUE))
  exp_list$exp_names<-rownames(exp_list)
  # re-order by run order, using time stamps
  exp_list<-exp_list[with(exp_list, order(as.POSIXct(mtime))), ]
  return(exp_list)
}


###################### read and analyze sequence of work loops #################

#' Read and analyze work loop files from a directory
#'
#' All-in-one function to import multiple workloop .ddf files from a directory,
#' sort them by mtime, analyze them, and store the resulting objects in an
#' ordered list.
#'
#' @param filepath Directory in which files are located
#' @param pattern Regular expression used to specify files of interest. Defaults
#' to all .ddf files within filepath.
#' @param sort_by Metadata by which files should be sorted to be in the correct
#' run order. Defaults to \code{mtime}, which is time of last modification of
#' files.
#' @param ... Additional arguments to be passed to \code{read_analyze_wl()},
#' \code{analyze_workloop()}, \code{select_cycles()}, or \code{read_ddf()}.
#'
#' @details Work loop data files will be imported and then arranged in the order
#' in which they were run (assuming run order is reflected in \code{mtime}).
#' Chiefly used in conjunction with \code{summarize_wl_trials()} and
#' \code{time_correct()} if time correction is desired.
#'
#' @return
#' A list containing \code{analyzed_workloop} objects, one for each file that is
#' imported and subsequently analyzed. The list is sorted according to the
#' \code{sort_by} parameter, which by default uses the time of last modification
#' of each file's contents (mtime).
#'
#' @inheritSection analyze_workloop Warning
#'
#' @references Josephson RK. 1985. Mechanical Power output from Striated Muscle
#'  during Cyclic Contraction. Journal of Experimental Biology 114: 493-512.
#'
#' @seealso
#' \code{\link{read_analyze_wl}},
#' \code{\link{get_wl_metadata}},
#' \code{\link{summarize_wl_trials}},
#' \code{\link{time_correct}}
#'
#' @author Shreeram Senthivasan
#'
#' @family data analyses
#' @family data import functions
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
#'                                     phase_from_peak = TRUE,
#'                                     cycle_def = "p2p", keep_cycles = 2:4)
#'
#' # or on your own directory
#' #my_analyzed_wls <- read_analyze_wl_dir("./my/file/path/")
#'
#' @export
read_analyze_wl_dir <- function(filepath,
                                pattern = "*.ddf",
                                sort_by = "mtime",
                                ...){
  # Generate list of filenames
  filename_list<-list.files(path=filepath,pattern=pattern,full.names=TRUE)
  if(length(filename_list)==0) stop("No files matching the pattern found at the given directory!")

  # Generate list of analyzed workloop objects
  wl_list<-lapply(filename_list,function(i) read_analyze_wl(i,...))

  # Sort list, likely by modification time
  if(is.null(attr(wl_list[[1]],sort_by))){
    warning("The provided sort_by argument is not a valid attribute. Defaulting to `mtime`.")
    sort_by<-"mtime"
  }
  wl_list<-wl_list[order(sapply(wl_list,function(i)attr(i,sort_by)))]
}

######################### summarize sequence of work loops #####################

#' Summarize work loop files
#'
#' Summarize important info from work loop files stored in the same folder
#' (e.g. a sequence of trials in an experiment) including experimental
#' parameters, run order, and \code{mtime}.
#'
#' @param wl_list List of \code{analyzed_workloop} objects, preferably one
#' created by \code{read_analyze_wl_dir()}.
#'
#' @details If several files (e.g. successive trials from one experiment) are
#' stored in one folder, use this function to obtain summary stats and
#' metadata and other parameters. This function requires a list of
#' \code{analyze_workloop} objects, which can be readily obtained by first
#' running \code{read_analyze_wl_dir()} on a specified directory.
#'
#' @return
#' A \code{data.frame} of information about the collection of workloop files.
#' Columns include:
#' \item{File_ID }{Name of the file}
#' \item{Cycle_Frequency }{Frequency of Position change}
#' \item{Amplitude }{amplitude of Position change}
#' \item{Phase }{Phase of the oscillatory cycle (in percent) at which
#' stimulation occurred. Somewhat experimental, please use with caution}
#' \item{Stimulus_Pulses }{Number of stimulation pulses}
#' \item{mtime }{Time at which file's contents were last changed (\code{mtime})}
#' \item{Mean_Work }{Mean work output from the selected cycles}
#' \item{Mean_Power }{Net power output from the selected cycles}
#'
#' @references Josephson RK. 1985. Mechanical Power output from Striated Muscle
#'  during Cyclic Contraction. Journal of Experimental Biology 114: 493-512.
#'
#' @seealso
#' \code{\link{read_analyze_wl_dir}},
#' \code{\link{get_wl_metadata}},
#' \code{\link{time_correct}}
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
#'                                                package = 'workloopR'),
#'                                     phase_from_peak = TRUE,
#'                                     cycle_def = "p2p",
#'                                     keep_cycles = 2:4
#'                                     )
#'
#' # now summarize
#' summarized_wls <- summarize_wl_trials(analyzed_wls)
#'
#' # or on your own directory
#' #my_meta <- read_analyze_wl_dir("./my/file/path/")
#' #my_summaries <- summarize_wl_trials(my_meta)
#'
#'
#' @export
summarize_wl_trials <- function(wl_list){
  if(class(wl_list)[[1]]!="list")
     stop("Please provide a list of analyzed workloop objects")
  if(!all(sapply(wl_list,function(x) 'analyzed_workloop' %in% class(x))))
    stop("The provided list includes elements that are not analyzed workloop objects")

  data.frame(
    File_ID = sapply(wl_list,function(i)attr(i,"file_id")),
    Cycle_Frequency = sapply(wl_list,function(i)attr(i,"cycle_frequency")),
    Amplitude = sapply(wl_list,function(i)attr(i,"amplitude")),
    Phase = sapply(wl_list,function(i)attr(i,"phase")),
    Stimulus_Pulses = sapply(wl_list,function(i)attr(i,"stimulus_pulses")),
    Stimulus_Frequency = sapply(wl_list,function(i)attr(i,"stimulus_frequency")),
    mtime = sapply(wl_list,function(i)attr(i,"mtime")),
    Mean_Work = sapply(wl_list,function(i)mean(attr(i,"summary")$Work)),
    Mean_Power = sapply(wl_list,function(i)mean(attr(i,"summary")$Net_Power))
  )
}

