# This file contains the public functions associated with
# the fitted Weibull/loglogistic model for the predict from data simulations  
##' @import eventPrediction methods
NULL

# Register the runjags and mcmc as a formally defined classes
setOldClass( "runjags" )
setOldClass( "mcmc" )

##' An S4 class containing a fitted survival model obtained from JAGS
##' of an Eventdata object
##' @slot model A "runjags" object
##' @slot event.data An EventData object used to fit a Weibull survial model
##' @slot shape.median Shape median of sampled distribution. 
##' @slot scale.1.median Median of the first of the scale parameters of sampled distribution.
##' @slot scale.2.median Median of rthe second of the scale parameters of sampled distribution.
##' @slot pi.median Median of the mixture coefficient of the sampled distribution. 
##' @slot mcmc.object The MCMC object containing the sampled distributions.  
##' @slot seed Seed used stored for reproducibility.
##' @slot jags.model The jags model used as a character string. 
##' 
##' @export
setClass( "EventModelBayesian",
   slots = list( model="runjags", 
                 event.data = "EventData",
                 shape.median = "numeric",
                 scale.1.median = "numeric",
                 scale.2.median = "numeric",
                 pi.median = "numeric",
                 mcmc.object = "mcmc",
                 seed = "numeric", 
                 jags.model = "character"
                 )
)


##' Displays information about EventModelBayesian objects
##' @name show
##' @rdname show-methods
##' @aliases show,EventModelBayesian-method
##' @export
setMethod("show",
          "EventModelBayesian",
  function(object) {
    cat( paste0( 
      "Parameters:\n",
      "Burnin (adapt+burnin); ", object@model$burnin,"\n",
      "Samples: ", object@model$sample, "\n",
      "Seed: ", object@seed, "\n",
      "Estimates:\n",
      "Shape.hat (median): ", object@shape.median, "\n",
      "Scale.hat[1] (median): ", object@scale.1.median, "\n",
      "Scale.hat[2] (median): ", object@scale.2.median, "\n",
      "pi.mixture.hat (median): ", object@pi.median, "\n" ) )
})



##' Plots the EventModelBayesian object
##' @rdname plot-methods
##' @param x EventModel, EventModelExtended, EventModelBayesian or LagTimeEstimate 
##' @param y For missing or EventModelExtended objects only
##' @param z dummy [EventModel only]
##' @param units Scale for the x-axis. "Days", "Months" or "Years"
##' @param xlab X-label titel
##' @param ylab Y-label titel
##' @param main Plot titel
##' @param ylim Y-range
##' @param xlim X-range 
##' @param smoothing [LagTimeEstimate objects only] If true will add a smoothed 
##' line to the curve and use that for when retreiving the optimal time value.
##' @param ... Additional arguments to be passed to the method
##' @aliases plot,EventModelBayesian,missing-method
setMethod( "plot", signature( x="EventModelBayesian", y="missing" ),
  function(x, units="Days", xlab=paste("Time in study [",units,"]",sep=""),
                ylab="", main="", ylim=NULL, xlim=NULL, ...) { 
    daysinyear <- standarddaysinyear()
    xscale <- eventPrediction:::GetScaleForKM(units,daysinyear)
    
    KM <- survfit(Surv( time, has.event) ~ 1, data=x@event.data@subject.data,...)
    
    plot(KM, lwd=c(2,1,1), col=c("red", "black", "black" ), 
         cex=0.5, conf.int=F,
         xlab=xlab, ylab=ylab,
         xlim=xlim, ylim=ylim,
         main=main,xscale=xscale)
    
    npts = 100
    xx <- seq( 0, max(x@event.data@subject.data$time, na.rm=T), 
               length = npts )
    
    # Calculate mix coefficients so do not get mixup between arms
    N.col <- ncol(x@mcmc.object)
    N.row <- nrow(x@mcmc.object)
    tmp <- sapply( seq_len(N.row), function(i){sum( x@mcmc.object[i,7:N.col] )/N.col} )
    mix.coef <- median( tmp )
    
    y.new <- ((1-mix.coef))*pweibull( xx, 
                        scale = x@scale.1.median, 
                        shape = x@shape.median, 
                        lower.tail = FALSE ) + 
                    mix.coef*pweibull( xx, 
                               scale = x@scale.2.median, 
                               shape = x@shape.median, 
                               lower.tail = FALSE )
    lines( xx/xscale, y.new, lwd=3, col="darkgreen" )
    
    legend( "topright", 1, c( "Data", "Model" ), col=c( "red", "darkgreen" ), 
            lty=c(1,1) )
  }
)


##' Plots the traces for each parameter of the fitted object
##' @name tracePlot
##' @param object EventModelBayesian object
##' @param ... Additional arguments to be passed to the method
##' @rdname tracePlot-methods
##' @aliases tracePlot,EventModelBayesian-method
##' @export
setGeneric( "tracePlot", 
            function( object, ... ) 
            standardGeneric("tracePlot") )


##' @rdname tracePlot-methods
##' @export
setMethod( "tracePlot", signature( object="EventModelBayesian" ),
           function( object ){
             idx <- which( names( object@mcmc.object[1,] ) %in% 
                             c( "shape", "pi", "scale[1]", "scale[2]" ) )
             xyplot( object@mcmc.object[,idx] )             
})

##' Plots the densities for each parameter of the fitted object
##' @name densityPlot
##' @param object EventModelBayesian object
##' @param ... Additional arguments to be passed to the method
##' @rdname densityPlot-methods
##' @aliases densityPlot,EventModelBayesian-method
##' @export
setGeneric( "densityPlot", 
            function( object, ... ) 
              standardGeneric("densityPlot") )

##' @rdname densityPlot-methods
##' @export
setMethod( "densityPlot", signature( object="EventModelBayesian" ),
           function( object ){
             idx <- which( names( object@mcmc.object[1,] ) %in% 
                             c( "shape", "pi", "scale[1]", "scale[2]" ) )
             densityplot( object@mcmc.object[,idx], 
               panel=function(x,...){
                 panel.densityplot(x,...)
                 panel.abline(v=median(x), col.line="red", lwd=3 )
                 # Maximum density value
                 a <- density(x)
                 panel.abline(v=a$x[ which.max( a$y ) ], col.line="gray", lwd=3)
               })
})


