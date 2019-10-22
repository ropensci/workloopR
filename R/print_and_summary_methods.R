# custom functions
# all written by Vikram B. Baliga (vbaliga@zoology.ubc.ca) and Shreeram
# Senthivasan
# last updated: 2019-10-22



# Generate type-specific output header for workloop objects
#' @noRd
print_muscle_stim_header <- function(x, include_time = TRUE) {
  type <- paste(toupper(substring(class(x)[1], 1, 1)),
    substring(class(x)[1], 2),
    sep = ""
  )
  if (include_time) {
    cat(paste0(
      "# ", type, " Data: ",
      ncol(x) - 1, " channels recorded over ",
      nrow(x) / attr(x, "sample_frequency"), "s\n"
    ))
  } else {
    cat(paste0("# ", type, " Data:\n\n"))
  }
}

# Print method

#' @export
print.muscle_stim <- function(x, n = 6, ...) {
  print_muscle_stim_header(x)
  cat(paste0("File ID: ", attr(x, "file_id"), "\n\n"))
  class(x) <- "data.frame"
  print(utils::head(x, n = n))
  if (n < nrow(x)) {
    cat(paste0("# \u2026 with ", nrow(x) - n, " more rows\n"))
  }
}

#' @export
print.analyzed_workloop <- function(x, n = 6, ...) {
  cat(paste0(
    "File ID: ",
    attr(x, "file_id")
  ))
  cat(paste0(
    "\nCycles: ",
    length(attr(x, "retained_cycles")),
    " cycles kept out of ",
    attr(x, "total_cycles")
  ))
  cat(paste0(
    "\nMean Work: ",
    round(mean(attr(x, "summary")$Work), 5),
    " J"
  ))
  cat(paste0(
    "\nMean Power: ",
    round(mean(attr(x, "summary")$Net_Power), 5),
    " W\n\n"
  ))
}

# Summary method

#' @export
summary.muscle_stim <- function(object, ...) {
  print_muscle_stim_header(object, ...)
  cat(paste0("\nFile ID: ", attr(object, "file_id")))
  cat(paste0("\nMod Time (mtime): ", attr(object, "mtime")))
  cat(paste0("\nSample Frequency: ", attr(object, "sample_frequency"),
             "Hz\n\n"))
  cat(paste0("data.frame Columns: \n"))
  for (i in 2:ncol(object)) {
    cat(paste0("  ", colnames(object)[i], " (", attr(object, "units")[i],
               ")\n"))
  }
  cat(paste0("\nStimulus Offset: ", attr(object, "stimulus_offset"), "s\n"))
  cat(paste0("Stimulus Frequency: ", attr(object, "stimulus_frequency"),
             "Hz\n"))
  cat(paste0("Stimulus Width: ", attr(object, "stimulus_width"), "ms\n"))
  cat(paste0("Stimulus Pulses: ", attr(object, "stimulus_pulses"), "\n"))
  cat(paste0("Gear Ratio: ", attr(object, "gear_ratio"), "\n"))
}

#' @export
summary.workloop <- function(object, ...) {
  NextMethod()
  cat(paste0("\nCycle Frequency: ", attr(object, "cycle_frequency"), "Hz\n"))
  cat(paste0(
    "Total Cycles (",
    switch(attr(object, "cycle_def"),
      "lo" = "L0-to-L0",
      "p2p" = "peak-to-peak",
      "t2t" = "trough-to-trough",
      "undefined"
    ),
    "): ",
    attr(object, "total_cycles"), "\n"
  ))
  if (!is.null(attr(object, "retained_cycles"))) {
    cat(paste0(
      "Cycles Retained: ",
      length(attr(object, "retained_cycles")),
      "\n"
    ))
  }
  cat(paste0(
    "Amplitude: ",
    attr(object, "amplitude"),
    attr(object, "units")[grep("Position", colnames(object))], "\n\n"
  ))
  if (attr(object, "position_inverted")) {
    cat("\nPlease note that Position is inverted!\n\n")
  }
}

#' @export
summary.tetanus <- function(object, ...) {
  NextMethod()
  cat(paste0("Stimulus Length: ", attr(object, "stimulus_length"), "s\n\n"))
}

#' @export
summary.analyzed_workloop <- function(object, ...) {
  summary(object[[1]], include_time = FALSE)
  cat("\n")
  print(attr(object, "summary"))
}

