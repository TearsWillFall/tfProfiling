% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_motif_coverage.R
\name{plot_motif_coverage}
\alias{plot_motif_coverage}
\title{Plot Mean Depth Coverage for TFBS}
\usage{
plot_motif_coverage(
  data = "",
  trend_line = TRUE,
  tf_name = "",
  sample_name = "",
  output_dir = ""
)
}
\arguments{
\item{data}{Data.frame with mean depth coverage data.}

\item{trend_line}{Draw trend line. Default TRUE}

\item{tf_name}{Transcription Factor name}

\item{sample_name}{Sample name}

\item{output_dir}{Directory to output results.}
}
\description{
This function takes a DATA.FRAME with Mean Depth Coverage data as input, and plots it
}
