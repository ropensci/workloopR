# custom functions
# all written by Vikram B. Baliga (vbaliga@zoology.ubc.ca) and Shreeram
# Senthivasan
# last updated: 2019-10-22



#' Import a batch of work loop or isometric data files from a directory
#'
#' Uses \code{read_ddf()} to read in workloop, twitch, or tetanus experiment
#' data from multiple .ddf files.
#'
#' @param file_path Path where files are stored. Should be in the same folder.
#' @param pattern Regex pattern for identifying relevant files in the file_path.
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
#' In addition, the following information is stored in each \code{data.frame}'s
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
#' \dontrun{
#' my_dat <- read_ddf_dir("./my/file/path/")
#' }
#'
#' @export
read_ddf_dir <- function(file_path,
                         pattern = "*.ddf",
                         sort_by = "mtime",
                         ...) {
  # Generate list of file_names
  file_name_list <- list.files(path = file_path,
                               pattern = pattern,
                               full.names = TRUE)
  if (length(file_name_list) == 0) {
    stop("No files matching the pattern found at the given directory!")
  }

  # Generate list of muscle_stim objects
  ms_list <- lapply(file_name_list, function(i) read_ddf(i, ...))

  # Sort list, likely by modification time
  if (is.null(attr(ms_list[[1]], sort_by))) {
    warning("The provided sort_by argument is not a valid attribute.
            \nDefaulting to `mtime`.")
    sort_by <- "mtime"
  }
  ms_list <- ms_list[order(unlist(lapply(ms_list, function(i)
    attr(i, sort_by))))]

  return(ms_list)
}


###################### file info for sequence of work loops ####################

#' Get file info for a sequence of experiment files
#'
#' Grab metadata from files stored in the same folder (e.g. a sequence of trials
#'  in an experiment).
#'
#' @param file_path Path where files are stored. Should be in the same folder.
#' @param pattern Regex pattern for identifying relevant files in the file_path.
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
#' files to all be work loops. Any file type is welcome (as long as the Regex
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
#' \dontrun{
#' my_meta <- get_wl_metadata("./my/file/path/")
#' }
#'
#' @export
get_wl_metadata <- function(file_path,
                            pattern = "*.ddf") {
  exp_list <- file.info(list.files(
    path = file_path, pattern = pattern,
    full.names = TRUE, recursive = TRUE
  ))
  exp_list$exp_names <- rownames(exp_list)
  # re-order by run order, using time stamps
  exp_list <- exp_list[with(exp_list, order(as.POSIXct(mtime))), ]
  return(exp_list)
}


###################### read and analyze sequence of work loops #################

#' Read and analyze work loop files from a directory
#'
#' All-in-one function to import multiple workloop .ddf files from a directory,
#' sort them by mtime, analyze them, and store the resulting objects in an
#' ordered list.
#'
#' @param file_path Directory in which files are located
#' @param pattern Regular expression used to specify files of interest. Defaults
#' to all .ddf files within file_path
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
#' \dontrun{
#' my_analyzed_wls <- read_analyze_wl_dir("./my/file/path/")
#' }
#'
#' @export
read_analyze_wl_dir <- function(file_path,
                                pattern = "*.ddf",
                                sort_by = "mtime",
                                ...) {
  # Generate list of file_names
  file_name_list <- list.files(path = file_path,
                               pattern = pattern,
                               full.names = TRUE)
  if (length(file_name_list) == 0) {
    stop("No files matching the pattern found at the given directory!")
  }

  # Generate list of analyzed workloop objects
  wl_list <- lapply(file_name_list, function(i) read_analyze_wl(i, ...))

  # Sort list, likely by modification time
  if (is.null(attr(wl_list[[1]], sort_by))) {
    warning("The provided sort_by argument is not a valid attribute.
            \nDefaulting to `mtime`.")
    sort_by <- "mtime"
  }
  return(wl_list <- wl_list[order(unlist(lapply(wl_list, function(i)
    attr(i, sort_by))))])
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
#' \dontrun{
#' my_meta <- read_analyze_wl_dir("./my/file/path/")
#' my_summaries <- summarize_wl_trials(my_meta)
#' }
#'
#' @export
summarize_wl_trials <- function(wl_list) {
  if (class(wl_list)[[1]] != "list") {
    stop("Please provide a list of analyzed workloop objects")
  }
  if (!all(unlist(lapply(wl_list,
                         function(x) "analyzed_workloop" %in% class(x))))) {
    stop("The provided list includes elements that are
         not analyzed workloop objects")
  }

  summarized <- data.frame(
    File_ID = vapply(wl_list, function(i) attr(i, "file_id"),
                     character(1)),
    Cycle_Frequency = vapply(wl_list, function(i) attr(i, "cycle_frequency"),
                             numeric(1)),
    Amplitude = vapply(wl_list, function(i) attr(i, "amplitude"),
                       numeric(1)),
    Phase = vapply(wl_list, function(i) attr(i, "phase"),
                   numeric(1)),
    Stimulus_Pulses = vapply(wl_list, function(i) attr(i, "stimulus_pulses"),
                             numeric(1)),
    Stimulus_Frequency = vapply(wl_list, function(i) {
      attr(i, "stimulus_frequency")
    }, numeric(1)),
    mtime = vapply(wl_list, function(i) attr(i, "mtime"),
                   numeric(1)),
    Mean_Work = vapply(wl_list, function(i) mean(attr(i, "summary")$Work),
                       numeric(1)),
    Mean_Power = vapply(wl_list, function(i) mean(attr(i, "summary")$Net_Power),
                        numeric(1))
  )

  return(summarized)
}
