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


##' Help-function to get a segment for plotting
getSegment <- function( x.start, x.end, y.start, npts, shape, scale ){
  xx <- seq( x.start, x.end, length = npts )
  yy <- pweibull( xx,   # Make sure curve starts where left curve ends.
                  scale = scale, 
                  shape = shape, 
                  lower.tail = FALSE )
  yy <- yy - yy[ 1 ] + y.start
  list( x = xx, y = yy )
}

##' Help-function that plots a segment
plotSegment <- function( seg, col, xscale ){
  lines( seg$x/xscale, seg$y, lwd=3, col=col )
  abline( v = max( seg$x ), lty=3 )
}


##' Plots the EventModelExtended object
##' @rdname plot-methods
##' @aliases plot,EventModelExtended,missing-method
setMethod( "plot", signature( x="EventModelExtended", y="missing" ),
  function(x, xlab=paste("Time in study [",units,"]",sep=""),
                ylab="", main="", ylim=NULL, xlim=NULL, ...) { 
      warning( paste0( "This function is not yet defined, please use: ",
                "plot( right.censored.model, left.truncated.model )" ) )
  }
)





##' Plots the EventModel,EventModelExtended object
##' @rdname plot-methods
##' @param units Scale for the x-axis. "Days", "Months" or "Years"
##' @param ... Additional arguments to be passed to the method
##' @aliases plot,EventModel,EventModelExtended-method
## Basically joins the curves where the right-cens model ends. Are there otherways to do this?
setMethod( "plot", signature( x="EventModel", y="EventModelExtended" ),
 function( x, y, units="Days", xlab=paste("Time in study [",units,"]",sep=""),
           ylab="", main="", ylim=NULL, xlim=NULL, ...) { 
   daysinyear <- standarddaysinyear()
   xscale <- eventPrediction:::GetScaleForKM(units,daysinyear)

   checkValidTimeCut( x, y )
   
   KM <- survfit(Surv( time, has.event) ~ 1, data=y@event.data@subject.data,...)
   
   plot(KM, lwd=c(2,1,1), col=c("red", "black", "black" ), 
        xlab=xlab, ylab=ylab,
        xlim=xlim, ylim=ylim,
        main=main,xscale=xscale)
   
   subject.data.2 <- y@event.data@subject.data
   chP1 <- y@time.cut
   chP2 <- max( subject.data.2$time )

   ## Piece 1 
   p1 <- getSegment( x.start=0, 
                     x.end=chP1,
                     y.start = 1,
                     npts = 100, 
                     shape = x@simParams@parameters$shape, 
                     scale = 1/x@simParams@parameters$rate )
   plotSegment( p1, "blue", xscale )
   
   
   ## Piece 2
   p2 <- getSegment( x.start=chP1, 
                     x.end=chP2,
                     y.start = tail( p1$y, 1 ),
                     npts = 100, 
                     shape = y@simParams@parameters$shape, 
                     scale = 1/y@simParams@parameters$rate )
   plotSegment( p2, "darkgreen", xscale )
   
   pos <- if(is.null(xlim)) 0.75*max(KM$time/xscale) else xlim[1] + 0.75*(xlim[2]-xlim[1])
   legend( pos, 1, c( "Data", "Model (right)", "Model (left)" ), col=c( "red", "blue", "darkgreen" ), lty=c(1,1) )
 }
)


##' Plots an arbitrary number of change points where the EventModel and a list 
##' of EventModelExtended objects are given
##' @rdname plot-methods
##' @param units Scale for the x-axis. "Days", "Months" or "Years"
##' @param ... Additional arguments to be passed to the method
##' @aliases plot,EventModel,list-method
## Basically joins the curves where the right-cens model ends. Are there otherways to do this?
setMethod( "plot", signature( x="EventModel", y="list" ),
           function( x, y, z, units="Days", xlab=paste("Time in study [",units,"]",sep=""),
                     ylab="", main="", ylim=NULL, xlim=NULL, ...) { 
             
             # Check that all entries in list are eventModelExtended
             if( !all( sapply( y, function(i) class( i )[1]=="EventModelExtended" ) ) ){
               stop( "All elements in list y need to be of class EventModelExtended" )
             }
             
             x.models <- append( list( x ), y ) 
             
             daysinyear <- standarddaysinyear()
             xscale <- eventPrediction:::GetScaleForKM(units,daysinyear)
          
             N.seg <- length(x.models)
             # Change points including boarders
             chPoints <- c( 0,  
                            sapply( 2:length(x.models), function(i) x.models[[i]]@time.cut ),
                            max( x.models[[ N.seg ]]@event.data@subject.data$time ) ) 
             
             
             for( i in 1:( length( x.models ) - 1 ) ) {
               checkValidTimeCut( x.models[[i]], x.models[[i+1]] )
             }
             
             KM <- survfit( Surv( time, has.event) ~ 1, 
                            data=x.models[[ N.seg ]]@event.data@subject.data,
                            ... )
             
             plot( KM, lwd=c(2,1,1), 
                   col=c("red", "black", "black" ), 
                   xlab=xlab, ylab=ylab,
                   xlim=xlim, ylim=ylim,
                   main=main,xscale=xscale)
             
             
             my.colors <- rep( c("blue", "darkgreen", "orange"), 100 )
             
             y.start <- 1
             for( i in 1:N.seg ){
               ## Piece i 
               pI <- getSegment( x.start=chPoints[ i ], 
                                 x.end=chPoints[ i+1 ],
                                 y.start = y.start,
                                 npts = 100, 
                                 shape = x.models[[i]]@simParams@parameters$shape, 
                                 scale = 1/x.models[[i]]@simParams@parameters$rate )
               plotSegment( pI, my.colors[i], xscale )
               y.start <- tail( pI$y, 1 )
             }
             
             pos <- if(is.null(xlim)) 0.75*max(KM$time/xscale) else xlim[1] + 0.75*(xlim[2]-xlim[1])
             my.text <- c( "Data", 
                           paste0( "Model (P", 1:N.seg, ")" ) )
             legend( pos, 1, my.text,  cex=0.5,
                     col=c( "red", my.colors[1:N.seg] ), lty=rep(1,4) )
           }
)



