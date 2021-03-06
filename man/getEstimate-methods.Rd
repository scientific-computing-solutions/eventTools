% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lagTimeEstimate.R
\docType{methods}
\name{getEstimate}
\alias{getEstimate}
\alias{getEstimate,LagTimeEstimate-method}
\title{Methods retrieving the optimal value for the lag time}
\usage{
getEstimate(object, smoothing = TRUE, ...)

\S4method{getEstimate}{LagTimeEstimate}(object, smoothing = TRUE, ...)
}
\arguments{
\item{object}{An EventData object}

\item{smoothing}{Whether a spline should be used to connect the points or not}

\item{...}{Additional arguments to be passed to the method}
}
\description{
Methods retrieving the optimal value for the lag time

Method returning the optimal value based on AIC or BIC
}

