# custom functions
# all written by Vikram B. Baliga (vbaliga@zoology.ubc.ca) and Shreeram
# Senthivasan
# last updated: 2019-10-22



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
#' "Position" - instantaneous position of the muscle,
#'   preferably in millimeters \cr
#' "Force" - force, preferably in millinewtons \cr
#' "Stim" - whether stimulation has occurred. All entries must be either 0 (no
#' stimulus) or 1 (stimulus occurrence).
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
#' wl_dat <-read_ddf(system.file("extdata", "workloop.ddf",
#'                               package = 'workloopR'))
#'
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
                           ...) {
  # Check for missing information
  if (missing(type)) stop("Please specify the experiment type!
                        \nThe type argument should be one of:
                        \nworkloop, tetanus, or twitch.")
  if (!(type %in% c("workloop", "tetanus", "twitch")) | length(type) != 1) {
    stop("Invalid experiment type!
         \nThe type argument should be one of:
         \nworkloop, tetanus, or twitch.")
  }
  if (!all(c("Position", "Force", "Stim") %in% names(x))) {
    stop("Couldn't find one or more of the following necessary columns:
         \nPosition, Force, Stim.
         \nPlease ensure that the columns match the naming conventions.")
  }
  if (missing(sample_frequency) & !("Time" %in% names(x))) {
    stop("Insufficient information to infer the sampling frequency.
         \nPlease provide a value for the sample_frequency argument
         \nor include a column named `Time` in the dataframe.")
  }

  # Consolidate time / sample frequency information
  if (!missing(sample_frequency)) {
    x$Time <- (seq_len(nrow(x)) - 1) / sample_frequency
  } else {
    sample_frequency <- 1 / (x$Time[2] - x$Time[1])
  }

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
  args <- list(...)
  if (!all(names(args) %in% valid_args)) {
    warning("One or more provided attributes do not match known attributes.
            \nThese attributes will not be assigned.")
  }
  for (i in intersect(names(args), valid_args)) {
    attr(x, i) <- args[[i]]
  }
  for (i in setdiff(valid_args, names(args))) {
    attr(x, i) <- NA
  }
  attr(x, "sample_frequency") <- sample_frequency
  if (is.na(attr(x, "gear_ratio"))) attr(x, "gear_ratio") <- 1
  if (type == "workloop") {
    if (is.na(attr(x, "position_inverted")))
      attr(x, "position_inverted") <- FALSE
  }

  # Assign classes and return
  class(x) <- c("muscle_stim", "data.frame")
  switch(type,
    "workloop" = class(x) <- c("workloop", class(x)),
    "tetanus" = class(x) <- c("tetanus", "isometric", class(x)),
    "twitch" = class(x) <- c("twitch", "isometric", class(x))
  )
  return(x)
}


########################## read_ddf files - work loops #########################

#' Import work loop or isometric data from .ddf files
#'
#' \code{read_ddf} reads in workloop, twitch, or tetanus experiment data from
#' .ddf files.
#'
#' @param file_name A .ddf file that contains data from a single workloop,
#' twitch, or tetanus experiment
#' @param file_id A string identifying the experiment. The file name is used by
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
#' In addition, the following information is stored in the \code{data.frame}'s
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
#'
#' @export
read_ddf <-
  function(file_name,
           file_id = NA,
           rename_cols = list(c(2, 3), c("Position", "Force")),
           skip_cols = 4:11,
           phase_from_peak = FALSE,
           ...) {
    # Import and checks
    if (missing(file_name)) stop("A file_name is required")
    if (!file.exists(file_name)) stop(paste0("File ", file_name, " not found!"))
    f <- file(file_name, "r")
    if (!grepl("DMC.*Data File", readLines(f, 1))) {
      close(f)
      stop("The input file does not appear to be a DMC Datafile (ddf)")
    }
    if (is.na(file_id)) file_id <- basename(file_name)

    # get metadata
    mtime <- file.info(file_name)$mtime

    # Setup for reading in file
    header <- c()
    units_table <- c()
    protocol_table <- c()

    # Read in Header
    while (!grepl("Calibration Data", (l <- readLines(f, 1)))) {
      header <- c(header, l)
    }
    sample_frequency <- as.numeric(sub(".*: ", "", header[1]))

    # Read in Calibration Table
    while (!grepl("Comments", (l <- readLines(f, 1)))) {
      units_table <- c(units_table, l)
    }
    units_table <- t(utils::read.table(
      text = units_table,
      row.names = 1,
      sep = "\t",
      stringsAsFactors = FALSE
    ))
    rownames(units_table) <- c()
    colnames(units_table) <- sub(" .*", "", colnames(units_table))
    units_table <- data.frame(units_table, stringsAsFactors = FALSE)
    units_table[3:5] <- lapply(units_table[3:5], as.numeric)
    units <- c("s", units_table$Units[-skip_cols + 1], "TTL")
    if (!all(units %in% c("s", "mm", "mN", "TTL"))) {
      warning("Non-standard units detected in ddf file!
              \nPlease note that calculations currently assume raw data
              are in seconds, millimeters, and millinewtons.")
    }

    # Read in Protocol Array
    while (!grepl("Protocol", readLines(f, 1))) {}
    readLines(f, 1) # Discard empty line
    while ((l <- readLines(f, 1)) != "") {
      protocol_table <- c(protocol_table, l)
    }
    protocol_table <- utils::read.table(
      text = protocol_table,
      sep = "\t",
      stringsAsFactors = FALSE,
      col.names = c(
        "Wait.s",
        "Then.action",
        "On.port",
        "Units",
        "Parameters"
      )
    )

    # Read in data
    while (!grepl("Test Data", (l <- readLines(f, 1)))) {}
    readLines(f, 1)
    dataz <- utils::read.table(
      text = readLines(f),
      header = TRUE,
      sep = "\t",
      stringsAsFactors = FALSE
    )
    if (any(!apply(dataz, 2, is.numeric))) {
      warning("The ddf file includes non-numeric data.
              \nPlease ensure that this is intentional before proceeding.")
    }
    close(f)

    # Parse file type
    read_filetype.ddf <- NULL
    switch(
      grep("Stim", protocol_table[[2]], value = TRUE)[1],
      "Stimulus-Train" = read_filetype.ddf <- read_wl_ddf,
      "Stimulus-Twitch" = read_filetype.ddf <- read_twitch_ddf,
      "Stimulus-Tetanus" = read_filetype.ddf <- read_tetanus_ddf,
      stop("Could not parse experiment type (workloop, twitch, or tetanus)!
               \nPlease ensure that the protocol section of the ddf header
               includes a label with one of the following:
               \nStimulus-Train, Stimulus-Twitch, or Stimulus-Tetanus.")
    )
    return(read_filetype.ddf(
      file_id = file_id,
      mtime = mtime,
      header = header,
      units_table = units_table,
      units = units,
      protocol_table = protocol_table,
      raw_data = dataz,
      sample_frequency = sample_frequency,
      rename_cols = rename_cols,
      skip_cols = skip_cols,
      phase_from_peak = phase_from_peak
    ))
  }


############################# rescale data matrix ##############################
# Rescales data in ddf files using the scale and offset parameters
#' @noRd
rescale_data <-
  function(dataz,
           unitz,
           sample_frequency,
           rename_cols,
           skip_cols) {
    rescaled <- mapply(
      function(raw, offset, scale) {
        if (!is.numeric(raw)) {
          return(raw)
        } else {
          (raw + offset) * scale
        }
      },
      dataz[unitz$Channel],
      unitz$Offset,
      unitz$Scale
    )

    rescaled <- data.frame(
      Time = (seq_len(nrow(dataz))) / sample_frequency,
      rescaled,
      Stim = dataz$Stim
    )

    # rename columns, if desired
    if (!is.null(rename_cols)) {
      names(rescaled)[rename_cols[[1]]] <- rename_cols[[2]]
    }

    return(rescaled[, -skip_cols])
  }


########################## read_ddf files - workloop ###########################
#' @noRd
read_wl_ddf <-
  function(raw_data,
           units_table,
           protocol_table,
           sample_frequency,
           rename_cols,
           skip_cols,
           ...) {
    # get info on experimental parameters
    stim_table <-
      utils::read.table(
        text = protocol_table[grepl("Stim", protocol_table$Then.action),
                              "Units"],
        sep = ",",
        col.names = c("offset",
                      "frequency",
                      "width",
                      "pulses",
                      "cycle_frequency")
      )
    cycle_table <-
      utils::read.table(
        text = protocol_table[grepl("Sine", protocol_table$Then.action),
                              "Units"],
        sep = ",",
        col.names = c("frequency", "amplitude", "total_cycles")
      )

    # use scale (and maybe offset) to convert Volts into units
    rescaled_data <- rescale_data(
      raw_data,
      units_table,
      sample_frequency,
      rename_cols,
      skip_cols
    )

    # construct and return workloop object
    return(workloop(
      data = rescaled_data,
      sample_frequency = sample_frequency,
      units_table = units_table,
      protocol_table = protocol_table,
      stim_table = stim_table,
      cycle_table = cycle_table,
      ...
    ))
  }


############################ read_ddf files - twitch ###########################
#' @noRd
read_twitch_ddf <-
  function(raw_data,
           units_table,
           protocol_table,
           sample_frequency,
           rename_cols = list(c(2, 3), c("Position", "Force")),
           skip_cols = 4:11,
           ...) {
    # get info on experimental parameters
    stim_table <-
      utils::read.table(
        text = protocol_table[grepl("Stim", protocol_table$Then.action),
                              "Units"],
        sep = ",",
        col.names = c("offset", "width")
      )
    stim_table$pulses <- rep(1, nrow(stim_table))

    # use scale (and maybe offset) to convert Volts into units
    rescaled_data <- rescale_data(
      raw_data,
      units_table,
      sample_frequency,
      rename_cols,
      skip_cols
    )

    # construct and return workloop object
    return(twitch(
      data = rescaled_data,
      sample_frequency = sample_frequency,
      units_table = units_table,
      protocol_table = protocol_table,
      stim_table = stim_table,
      ...
    ))
  }


########################## read_ddf files - tetanus ##########################
#' @noRd
read_tetanus_ddf <-
  function(raw_data,
           units_table,
           protocol_table,
           sample_frequency,
           rename_cols = list(c(2, 3), c("Position", "Force")),
           skip_cols = 4:11,
           ...) {
    # get info on experimental parameters
    stim_table <-
      utils::read.table(
        text = protocol_table[grepl("Stim", protocol_table$Then.action),
                              "Units"],
        sep = ",",
        col.names = c("offset", "frequency", "width", "length")
      )
    stim_table$pulses <-
      as.integer(floor(stim_table$frequency * stim_table$length))

    # use scale (and maybe offset) to convert Volts into units
    rescaled_data <- rescale_data(
      raw_data,
      units_table,
      sample_frequency,
      rename_cols,
      skip_cols
    )

    # construct and return workloop object
    return(tetanus(
      data = rescaled_data,
      sample_frequency = sample_frequency,
      units_table = units_table,
      protocol_table = protocol_table,
      stim_table = stim_table,
      ...
    ))
  }
