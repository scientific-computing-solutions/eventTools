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
  if( ( y@time.cut - max(x@event.data@subject.data$time, na.rm=T) ) > 2 ) {
    warning( paste0( "The change point should equal the largest time in the right censored model!", 
                     "That or all patients have had events! Regardless, this situation may result", 
                     "undefined behaiviour!" ) )
  } 
}
