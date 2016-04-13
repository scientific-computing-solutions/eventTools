# Some common functionality and constants

##' Sets the smoothing parameter, usally in (0,1]. 
##' 
##' @return The spline resolution
##' @export 
splineResolution <- function() {
  getOption('eventTools.spline.resolution', 0.5 )
}

