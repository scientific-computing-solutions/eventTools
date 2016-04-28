# Some common functionality and constants

##' Sets the smoothing parameter, usally in (0,1]. 
##' 
##' @return The spline resolution
##' @export 
splineResolution <- function() {
  getOption('eventTools.spline.resolution', 0.5 )
}

##' Checks that the change point is ok
##' @param x EventModel object
##' @param y EventModelExtended object
checkValidTimeCut <- function( x, y ) {
  if( y@time.cut < max(x@event.data@subject.data$time, na.rm=T) ) {
    warning( "The change point is less than the largest time in the right censored model!" )
  } 
}
