% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analysis_functions.R
\name{time_correct}
\alias{time_correct}
\title{Time correction for work loop experiments}
\usage{
time_correct(x)
}
\arguments{
\item{x}{A \code{data.frame} with summary data, e.g. an object created by
\code{summarize_wl_trials()}.}
}
\value{
A \code{data.frame} that additionally contains:
\item{Time_Corrected_Work }{Time corrected work output, transformed from
 \code{$Mean_Work}}
\item{Time_Corrected_Power }{Time corrected net power output, transformed
from \code{$Mean_Power}}

And new attributes:
\item{power_difference }{Difference in mass-specific net power output
between the final and first trials.}
\item{time_difference }{Difference in mtime between the final and first
trials.}
\item{time_correction_rate }{Overall rate; \code{power_difference} divided
 by \code{time_difference}.}
}
\description{
Correct for potential degradation of muscle over time.
}
\details{
This function assumes that across a batch of successive trials, the
stimulation parameters for the first and final trials are identical. If not,
DO NOT USE. Decline in power output is therefore assumed to be a linear
function of time. Accordingly, the difference between the final and first
trial's (absolute) power output is used to 'correct' trials that occur in
between, with explicit consideration of run order and time elapsed (via
mtime). A similar correction procedure is applied to work.
}
\examples{

library(workloopR)

# batch read and analyze files included with workloopR
analyzed_wls <- read_analyze_wl_dir(system.file("extdata/wl_duration_trials",
                                                package = 'workloopR'),
                                    phase_from_peak = TRUE,
                                    cycle_def = "p2p", keep_cycles = 2:4)

# now summarize
summarized_wls <- summarize_wl_trials(analyzed_wls)


# mtimes within the package are not accurate, so we'll supply
# our own vector of mtimes
summarized_wls$mtime <- read.csv(
                          system.file(
                            "extdata/wl_duration_trials/ddfmtimes.csv",
                            package="workloopR"))$mtime

# now time correct
timecor_wls <- time_correct(summarized_wls)
timecor_wls


}
\seealso{
\code{\link{summarize_wl_trials}}

Other workloop functions: 
\code{\link{analyze_workloop}()},
\code{\link{fix_GR}()},
\code{\link{get_wl_metadata}()},
\code{\link{invert_position}()},
\code{\link{read_analyze_wl_dir}()},
\code{\link{read_analyze_wl}()},
\code{\link{select_cycles}()},
\code{\link{summarize_wl_trials}()}

Other batch analyses: 
\code{\link{get_wl_metadata}()},
\code{\link{read_analyze_wl_dir}()},
\code{\link{summarize_wl_trials}()}
}
\author{
Vikram B. Baliga and Shreeram Senthivasan
}
\concept{batch analyses}
\concept{workloop functions}
