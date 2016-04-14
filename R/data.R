##' A simulated clinical trial with 979 patients. A lag time of 5.5 months has been assumed.
##' HR(after lag)=0.25, Control median=5.8. Survival function Weibull with shape 1.2. No dropouts.
##' Linear recruitment. 
##'
##' @name lag.data
##'
##' @format A data frame with 979 rows and 10 variables:
##' \describe{
##'   \item{subject}{Subject identifier}
##'   \item{site}{Site (center/hospital) identifier}
##'   \item{arm}{Treatment arm}
##'   \item{randDate}{Randomization date}
##'   \item{eventDate}{Event date}
##'   \item{time}{Time to event or withdrawal from study}
##'   \item{hasEvent}{Whether subject has had event or not}
##'   \item{lastDate}{Last date}
##'   \item{withdrawn}{If subject has withdrawn or not.}
##'   \item{eventType}{Type of event}
##' }
##' 
##' @examples
##' data(lag.data)
NULL