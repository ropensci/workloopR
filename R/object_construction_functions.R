# custom functions
# all written by Vikram B. Baliga (vbaliga@zoology.ubc.ca) and Shreeram
# Senthivasan
# last updated: 2019-10-22



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
           ...) {
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
    return(data)
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
    phase <-
      (
        which.max(data$Stim) - which.max(data$Position)
      ) / sample_frequency * stim_table$cycle_frequency[1]
    if (!phase_from_peak) {
      phase <- phase + 0.25
    }
    # convert 0-1 scale to -50 to +50
    attr(data, "phase") <- (((phase + 0.5) %% 1) - 0.5) * 100

    attr(data, "position_inverted") <- FALSE
    class(data) <- c("workloop")
    return(muscle_stim(
      data = data,
      stim_table = stim_table,
      sample_frequency = sample_frequency,
      ...
    ))
  }

#' @noRd
tetanus <-
  function(data, stim_table, ...) {
    attr(data, "stimulus_frequency") <- stim_table$frequency[1]
    attr(data, "stimulus_length") <- stim_table$length[1]
    class(data) <- c("tetanus", "isometric")
    return(muscle_stim(
      data = data,
      stim_table = stim_table,
      ...
    ))
  }

#' @noRd
twitch <-
  function(data, ...) {
    class(data) <- c("twitch", "isometric")
    return(muscle_stim(data = data, ...))
  }
