% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_transformation_functions.R
\name{fix_GR}
\alias{fix_GR}
\title{Adjust for the gear ratio of a motor arm}
\usage{
fix_GR(x, GR = 1)
}
\arguments{
\item{x}{A \code{muscle_stim} object}

\item{GR}{Gear ratio, set to 1 by default}
}
\value{
An object of the same class(es) as the input (\code{x}). The function
 will multiply \code{Position} by (1/GR) and multiply \code{Force} by GR,
 returning an object with new values in \code{$Position} and \code{$Force}.
 Other columns and attributes are welcome and will simply be passed on
 unchanged into the resulting object.
}
\description{
Fix a discrepancy between the gear ratio of the motor arm used and the gear
ratio recorded by software.
}
\details{
The \code{muscle_stim} object can be of any type, including
\code{workloop}, \code{twitch}, or \code{tetanus}.

If you have manually constructed the object via \code{as_muscle_stim()},
the \code{muscle_stim} object should have columns as follows: \cr
\code{Position}: length change of the muscle; \cr
\code{Force}: force \cr
}
\examples{

library(workloopR)

# import the workloop.ddf file included in workloopR
wl_dat <-read_ddf(system.file("extdata", "workloop.ddf",
                              package = 'workloopR'),
                  phase_from_peak = TRUE)

# apply a gear ratio correction of 2
# this will multiply Force by 2 and divide Position by 2
wl_fixed <- fix_GR(wl_dat, GR = 2)

# quick check:
max(wl_fixed$Force) / max(wl_dat$Force) # 5592.578 / 2796.289 = 2
max(wl_fixed$Position) / max(wl_dat$Position) # 1.832262 / 3.664524 = 0.5

}
\seealso{
\code{\link{analyze_workloop}},
\code{\link{read_analyze_wl}},
\code{\link{read_analyze_wl_dir}}

Other data transformations: 
\code{\link{invert_position}()},
\code{\link{select_cycles}()}

Other workloop functions: 
\code{\link{analyze_workloop}()},
\code{\link{get_wl_metadata}()},
\code{\link{invert_position}()},
\code{\link{read_analyze_wl_dir}()},
\code{\link{read_analyze_wl}()},
\code{\link{select_cycles}()},
\code{\link{summarize_wl_trials}()},
\code{\link{time_correct}()}

Other twitch functions: 
\code{\link{invert_position}()},
\code{\link{isometric_timing}()}

Other tetanus functions: 
\code{\link{invert_position}()}
}
\author{
Vikram B. Baliga
}
\concept{data transformations}
\concept{tetanus functions}
\concept{twitch functions}
\concept{workloop functions}
