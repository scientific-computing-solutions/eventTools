% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/KMRightCensorData.R
\name{KMRightCensor}
\alias{KMRightCensor}
\title{KMRightCensor the Event Data}
\usage{
KMRightCensor(object, time.cut)
}
\arguments{
\item{object}{The EventData object}

\item{time.cut}{The time in days when to cut the data}
}
\value{
\code{EventData} object with the data censored at the appropriate point
}
\description{
Creates a new EventData from an existing one where each subject
is censored at subject time T (not study time)
}

