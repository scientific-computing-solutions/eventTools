# This file contains the public functions associated with
# the fitted Weibull/loglogistic model for the predict from data simulations  
##' @import eventPrediction methods
NULL

##' An S4 class containing a fitted survival model
##' of an Eventdata object
##' @slot model An S3 "weibreg" object
##' @slot event.data An EventData object used to fit a Weibull survial model
##' @slot simParams The \code{FromDataSimParam} object which contains the information
##' needed to generate the survial times
##' @export
setClass("EventModelExtended", 
         slots = list( model="list",  #Would prefer to put "weibreg" here but doesn't work?
                       event.data="EventData",
                       time.cut="numeric",
                       simParams="FromDataSimParam" ) #,
      #   validity = function(object) { T },
      #   contains = "eventPrediction"
)


##' @export
setMethod("show",
          "EventModelExtended",
          function(object) {
            print(object@model[[1]])
          })




##' @export
setMethod( "plot",
           signature( x="EventModelExtended", y="missing" ),
           function(x, units="Days", xlab=paste("Time in study [",units,"]",sep=""),
                    ylab="", main="", ylim=NULL, xlim=NULL, ...) { 
             
             xscale <- eventPrediction:::GetScaleForKM(units,daysinyear)
             daysinyear <- standarddaysinyear()
             
             KM <- survfit(Surv( time, has.event) ~ 1, data=x@event.data@subject.data,...)
             
             plot(KM, lwd=c(2,1,1), col=c("red", "black", "black" ), 
                  xlab=xlab, ylab=ylab,
                  xlim=xlim, ylim=ylim,
                  main=main,xscale=xscale)
             
             npts = 100
             xx <- seq( x@time.cut, max(x@event.data@subject.data$time, na.rm=T), length = npts )
             y <- pweibull( xx, 
                            scale = 1/x@simParams@parameters$rate, 
                            shape = x@simParams@parameters$shape, 
                            lower.tail = FALSE)
             lines( xx, y, lwd=3, col="orange" )
             
             pos <- if(is.null(xlim)) 0.75*max(KM$time/xscale) else xlim[1] + 0.75*(xlim[2]-xlim[1])
             
             legend(pos, 1, c( "Data", "Model" ), col=c( "red", "orange" ), lty=c(1,1))
             
           }
)



