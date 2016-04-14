# This file contains the public functions associated with
# the fitted Weibull/loglogistic model for the predict from data simulations  
##' @import eventPrediction methods
NULL

##' An S4 class containing a fitted survival model
##' of an Eventdata object
##' @slot model An S3 "weibreg" object
##' @slot event.data An EventData object used to fit a Weibull survial model
##' @slot time.cut The change point
##' @slot simParams The \code{FromDataSimParam} object which contains the information
##' needed to generate the survial times
##' @export
setClass("EventModelExtended", 
   slots = list( model="list",  #Would prefer to put "weibreg" here but doesn't work?
                 event.data="EventData",
                 time.cut="numeric",
                 simParams="FromDataSimParam" )
)


##' Displays information about EventModelExtended objects
##' @name show
##' @param object EventModelExtended object
##' @rdname show-methods
##' @aliases show,EventModelExtended-method
##' @export
setMethod("show",
          "EventModelExtended",
  function(object) {
    print(object@model[[1]])
})



##' Plots the EventModelExtended object
##' @name plot
##' @param ... Additional arguments to be passed to the method
##' @rdname plot-methods
##' @aliases plot,EventModelExtended,missing-method
##' @export
setMethod( "plot", signature( x="EventModelExtended", y="missing" ),
  function(x, xlab=paste("Time in study [",units,"]",sep=""),
                ylab="", main="", ylim=NULL, xlim=NULL, ...) { 
      warning( paste0( "This function is not yet defined, please use: ",
                "plot( right.censored.model, left.truncated.model )" ) )
  }
)


##' Plots the EventModel,EventModelExtended object
##' @name plot
##' @param units Scale for the x-axis. "Days", "Months" or "Years"
##' @rdname plot-methods
##' @aliases plot,EventModel,EventModelExtended-method
##' @export
## Basically joins the curves where the right-cens model ends. Are there otherways to do this?
setMethod( "plot", signature( x="EventModel", y="EventModelExtended" ),
 function( x, y, units="Days", xlab=paste("Time in study [",units,"]",sep=""),
           ylab="", main="", ylim=NULL, xlim=NULL, ...) { 
   xscale <- eventPrediction:::GetScaleForKM(units,daysinyear)
   daysinyear <- standarddaysinyear()
   
   KM <- survfit(Surv( time, has.event) ~ 1, data=y@event.data@subject.data,...)
   
   plot(KM, lwd=c(2,1,1), col=c("red", "black", "black" ), 
        xlab=xlab, ylab=ylab,
        xlim=xlim, ylim=ylim,
        main=main,xscale=xscale)
   
   npts = 100
   xx <- seq( y@time.cut, max(y@event.data@subject.data$time, na.rm=T), length = npts )
   y.end.right.cens <- pweibull( y@time.cut,   # Make sure curve starts where left curve ends.
                                 scale = 1/x@simParams@parameters$rate, 
                                 shape = x@simParams@parameters$shape, 
                                 lower.tail = FALSE)
   y.new <- pweibull( xx, 
                      scale = 1/y@simParams@parameters$rate, 
                      shape = y@simParams@parameters$shape, 
                      lower.tail = FALSE) 
   # Connect start point of left-truncated curve with 
   # the end of the right-censored curve
   y.new <- y.new - y.new[ 1 ] + y.end.right.cens
   lines( xx, y.new, lwd=3, col="orange" )
   
   y.right <- seq(.99,.01,by=-.01)
   right.model <- predict(x@model, type="quantile", p=rev(y.right))[1,]/xscale
   idx <- which( seq(.99,.01,by=-.01) > y.end.right.cens )
   lines( right.model[idx], y.right[idx], col="blue", type="l", lwd=3 )
   abline( v = y@time.cut, lty=3 )
   pos <- if(is.null(xlim)) 0.75*max(KM$time/xscale) else xlim[1] + 0.75*(xlim[2]-xlim[1])
   legend( pos, 1, c( "Data", "Model" ), col=c( "red", "orange" ), lty=c(1,1) )
 }
)


