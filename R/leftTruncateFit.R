##' @importFrom eventPrediction FromDataParam
NULL

##' LeftFit the Event Data
##' 
##' Creates an EventModel based on the left truncated data
##' @param object The EventData object
##' @param time.cut The time in days when to cut the data
##' @param ... Additional arguments to be passed to the method
##' @rdname LeftFit-methods
##' @name LeftFit
##' @return \code{EventData} object with the data censored at the appropriate point
##' @export
setGeneric( "LeftFit", 
            function( object, time.cut, ... ) 
  standardGeneric("LeftFit") )


##' LeftFit the Event Data
##' 
##' @param object The EventData object
##' @param time.cut The time in days when to cut the data
##' @param dist Distribution (weibull)
##' @param ... Additional arguments to be passed to the method
##' @export
setMethod( "LeftFit", signature=c( object="EventData", time.cut="numeric" ),
           function( object, time.cut, dist="weibull", ... ){
  if( !dist %in% c( "weibull" ) ){
     stop( "dist must be weibull")
  }
  
  if( nrow( object@subject.data ) == 0 ) stop( "Empty data frame!" ) 

  subject.data <- object@subject.data
  fit.data <- subject.data[ subject.data$time > time.cut, ]
  
  if( sum( fit.data$has.event )==0 ){
    stop("Cannot fit a model to a dataset with no events, time.cut too high?")
  }

  fit.data$time.left <- time.cut
  model <- tryCatch( weibreg( Surv( time.left, time, has.event ) ~ 1, 
                              data = fit.data ), 
                     error=function(e){
                       warning( "Regression did not converge, exponential distribution will be used!")
                       weibreg( Surv( time.left, time, has.event ) ~ 1, 
                                data = fit.data, shape=1 )
                     })

  if( !is.null( model )) {
    shape <- 1
    if( !is.na(model$coefficients[2]) ) {
      shape <- as.numeric( exp( model$coefficients[2] ) ) 
    }
    scale <- as.numeric( exp( model$coefficients[1] ) ) 
    sigma <- model$var
    if( length( sigma ) == 1 ){
      # Fix for exponential distribution (Shape has no uncertainty)
      sigma <- matrix( c( sigma, 0, 0, 0 ), nrow=2 )
    }
    rate <- 1/scale
    # Just for consistency with survreg
    sigma[2] <- -sigma[2]
    sigma[3] <- -sigma[3]
    new( "EventModelExtended",
         model=list( model ),
         event.data=object,
         time.cut=time.cut,
         simParams=FromDataParam( type=dist,
                                  rate=rate,
                                  shape=shape,
                                  sigma=sigma ) ) 
  } else {
    NULL
  }
})


