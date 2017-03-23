# Functions to estimate the change-point in time when HR is changes

##' An S4 class containing an LagTimeEstimate Object
##' @slot event.data An EventData object used to fit a Weibull survial model.
##' @slot times time vector with the values at which XIC was calculated.
##' @slot XIC vector of AIC or BIC values of fits at different times.
##' @slot criterion AIC or BIC (AIC)
##' @export
setClass("LagTimeEstimate", 
         slots = list( event.data="EventData",
                       times="numeric",
                       XIC="numeric",
                       criterion="character" )
)

##' Methods to estimate lag times
##' @rdname estimateLagTime-methods
##' @name estimateLagTimeG
##' @title GMethod to estimate lag times
##' @param object EventData object
##' @param t.start Start time for search interval [Integer, months]
##' @param t.stop end time for search interval [Integer, months]
##' @param dt Time increment for search [Integer, months]
##' @param dist Survival distribution, see survival package ("exponential")
##' @param criterion Estimation criterion for the point estimate, AIC or BIC (AIC).
##' @export
setGeneric( "estimateLagTime",
            function( object,  t.start, t.stop, dt, dist="exponential", 
                      criterion="AIC" )
              standardGeneric( "estimateLagTime") 
)


##' Methods to estimate lag times using survival package
##' @rdname estimateLagTime-methods
##' @name estimateLagTime
##' @aliases estimateLagTime,EventData-method
##' @export
setMethod( "estimateLagTime",signature=c( "EventData" ),
 function( object, t.start, t.stop, dt, dist="exponential", criterion="AIC" ) {
   if( t.start < 0 ){ stop( "Start time needs to be greater or equal to zero" ) }
   if( t.stop < t.start ){ stop( "Stop time needs to be greater than start time" ) }
   if( dt < 0 ){ stop( "Time increment (dt) needs to be positive" ) }
   
   lagT <- seq( t.start, t.stop, dt )*standarddaysinyear()/12
   xic <- rep( NA, length( lagT ))
   for( i in seq_along( lagT ) ){
     tmp.data <-  object@subject.data
     tmp.data$ind <- tmp.data$time > lagT[i]
     model <- survreg( Surv( time, has.event )~ind, 
                       data=tmp.data, 
                       dist=dist )
     xic[i] <- ifelse( criterion=="AIC", AIC( model ), BIC( model ) )
   }

   new( "LagTimeEstimate",
        event.data=object,
        times = lagT, 
        XIC = xic, 
        criterion = criterion ) 
})


##' Methods retrieving the optimal value for the lag time
##' @param object An EventData object 
##' @param smoothing Whether a spline should be used to connect the points or not
##' @param ... Additional arguments to be passed to the method
##' @rdname getEstimate-methods
##' @export
setGeneric( "getEstimate",
            function( object,  smoothing=TRUE, ... )
              standardGeneric( "getEstimate") 
)


##' Method returning the optimal value based on AIC or BIC
##' @rdname getEstimate-methods
##' @export
setMethod( "getEstimate",signature=c( "LagTimeEstimate" ),
           function( object, smoothing = TRUE, ... ) {
    if( is.null( object@times ) ){ return( NULL ) }
     if( smoothing ) {
       my.smooth <- smooth.spline( object@XIC~object@times, spar=splineResolution() )
       my.smooth$x[ which.min( my.smooth$y ) ]
      } 
     else {
       object@times[ which.min( object@XIC ) ]
      }
})



##' Plot method for LagTimeEstimate objects
##' @name plot
##' @rdname plot-methods
##' @aliases plot,LagTimeEstimate,missing-method
##' @export
setMethod( "plot",
   signature( x="LagTimeEstimate", y="missing" ),
   function( x, smoothing=TRUE, ... ) { 
     main <- if( !is.null( getEstimate(x, smoothing) ) ){  
       paste0( "Change point estimate = ", round( getEstimate(x), 2 ) )
     }
     { NULL }
     
     my.smooth <- if( smoothing == TRUE ) {
       smooth.spline( x@XIC~x@times, spar=splineResolution() )
     }
     { NULL }
     
     plot( x@times, x@XIC, type="l", col="black", main=main, 
           xlab="Time [Days]", ylab=x@criterion,... )
     abline( v=getEstimate(x), col="black", lty=2 )
     if( !is.null( my.smooth)  ){
       lines( my.smooth$x, my.smooth$y, col="darkgreen", lwd=2 )
       legend( "topright", c( "Smoothed", "Values" ), col=c( "darkgreen", "black" ), lty=c(1,1) )
     }

   }
)


