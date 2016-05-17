##' Conditional function for piece wise model
##' @param t.conditional time vector to condition on  
##' @param params paramter object with study settings
##' @param HR HR
rcpweibull <- function( t.conditional, params, HR ) {
  params <- as.list( params )

  # Generate times for patients in the first time-interval
  # using the right-censored model
  idx <- which( t.conditional < params$time.cut.1 )
  my.times <- rep( NA, length( t.conditional ) )
  my.times[ idx ] <- t.conditional[ idx ]*params$rate
  my.times[ idx ] <- ((my.times[ idx ]^params$shape + 
                         rexp(length(my.times[ idx ])))^(1/params$shape))/params$rate

  # Loop over all other time intervals
  
  # All times generated from previous model that exceeds time.cut.i needs 
  # to be resimulated with model i conditioned on time.cut.i
  var.time.cut <- paste0( "time.cut.", params$N+1 )
  params[[ var.time.cut ]] = Inf
  for( i in seq_len( params$N ) ) {
    var.time.cut <- paste0( "time.cut.", i )
    var.time.cut.next <- paste0( "time.cut.", i+1 )
    var.rate <- paste0( "rate.lft.", i )
    var.shape <- paste0( "shape.lft.", i )
    
    idx2 <- which( my.times > params[[ var.time.cut ]] )

    # All times simulated that exceeds time.cut.i will be set to time.cut.i and
    # simulated from the left-truncated model
    t.conditional[ idx2 ] = params[[ var.time.cut ]]
    my.times[ idx2 ] = NA

    # Simulate all with undefined times that are less than nexd change-point
    idx3 <- which( is.na( my.times ) & t.conditional < params[[ var.time.cut.next ]] ) 
  
    my.times[ idx3 ] <- t.conditional[ idx3 ]*params[[ var.rate ]]
    my.times[ idx3 ] <- ((my.times[ idx3 ]^params[[ var.shape ]] + 
                          rexp(length(my.times[ idx3 ])))^(1/params[[ var.shape ]] ))/params[[ var.rate ]]
  }
  my.times 
}


##' Simulate method for use when predicting from piecewise Weibull
##' @param object The model fitted on Right censored data
##' @param leftobjects The model fitted on Left truncated/Right censored data. Either a single
##' or a list of many EventModelExtended objects. 
##' @param ... Additional arguments to be passed to the method
##' @rdname simulatePW-methods
##' @export 
setGeneric( "simulatePW",
  function( object, leftobjects, ... )
           standardGeneric( "simulatePW")
)



# ##' @rdname simulatePW-methods
# ##' @export
setMethod( "simulatePW",signature=c( "EventModel", "EventModelExtended" ),
  function( object, leftobjects, ... ){
    # Currently only Weibull allowed
    if( object@simParams@type != "weibull" ||
        leftobjects@simParams@type != "weibull" ){
      stop( "Only weibull is currenlty allowed for the piecewise simulation!")
    }
    checkValidTimeCut(  object, leftobjects )
  simulate.Internal( leftobjects@event.data, 
                     object@simParams, 
                     list( leftobjects@simParams ), 
                     list( leftobjects@time.cut ), ... )
})

##' @rdname simulatePW-methods
##' @export
setMethod( "simulatePW",signature=c( "EventModel", "list" ),
           function( object, leftobjects, ... ){
             leftobjects <- as.list( leftobjects )
             # Currently only Weibull allowed
             if( any( sapply( leftobjects, function(x) class(x) != "EventModelExtended" ) ) ){
               stop( "All leftobjects need to be of type EventModelExtended!")
             }
             if( object@simParams@type != "weibull" ||
                 any( sapply( leftobjects, function(x) x@simParams@type ) != "weibull" ) ){
               stop( "Only weibull is currenlty allowed for the piecewise simulation!")
             }
             
             x.models <- append( list( object ), leftobjects ) 
             simParamsList <- list()
             timeCutList <- list()
             for( i in 1:( length( x.models ) - 1 ) ) {
               checkValidTimeCut( x.models[[i]], x.models[[i+1]] )
               simParamsList <- append( simParamsList, x.models[[i+1]]@simParams )
               timeCutList   <- append( timeCutList, x.models[[i+1]]@time.cut    )
             }
             
             simulate.Internal( x.models[[ length( x.models ) ]]@event.data, 
                                object@simParams, 
                                simParamsList, 
                                timeCutList, ... )
           })



##' Only used internally
##' @param data EventData object 
##' @param SimParamsRgt Parameters for fitted KM right censored model
##' @param SimParamsList List of parameters for Fitted left truncated model
##' @param timeCutList ;ist of change points, should equal length of SimParamsList
##' @param accrualGenerator accrualGenerator object 
##' @param Naccrual Number of subjects to recruit
##' @param Nsim Number of simulations
##' @param seed Random seed
##' @param limit Used to set CI limit
##' @param longlagsettings LonLagSettings object
##' @param dropout Dropout object 
simulate.Internal <- function( data, 
                               SimParamsRgt,
                               SimParamsList, 
                               timeCutList, 
                               accrualGenerator=NULL, Naccrual=0, 
                               Nsim=1e4, 
                               seed=NULL, limit=0.05, 
                               longlagsettings=NULL,
                               dropout=NULL ){
  
  # These are not applicable here
  r <- NULL
  HR <- NULL
  # Validate the arguments
  eventPrediction:::validate.simulate.arguments(accrualGenerator,Naccrual,Nsim,seed,
                              limit,longlagsettings,HR,r,data)  
  if( !all( sapply( SimParamsList, function(i) class( i )[1]=="FromDataSimParam" ) ) ){
    stop( "All elements in SimParamsList need to be of class FromDataSimParam" )
  }
  if( length( SimParamsList ) != length( timeCutList ) ){
    stop( "Length of SimParamsList should equal timeCutList" )
  }
  if( !all( sapply( timeCutList, function(i) length(i)>0 ) ) ){
    stop( "Negative change point(s) in timeCutList!" )
  }
  
  #calculate the dropout rate and shape for drop out
  dropoutctrlSpec <- eventPrediction:::CtrlSpecFromList(dropout,eventtext="",1)[[1]]
  dropout.shape <- if(is.null(dropout) || is.null(dropout$shape)) 1 else dropout$shape
  dropout.rate <- log(2)^(1/dropout.shape)/dropoutctrlSpec@median
  
  #set seed to be used
  if(!is.null(seed)) set.seed(seed)
  
  #pre-process data to deal with subjects censored
  #a long time in the past
  indat <- eventPrediction:::DealWithReportingLag(data@subject.data,longlagsettings)      
  
  #create matrix of subject recruitment times including additional
  #accrual we have a matrix with 1 row per simulation, 1 column per subject
  rec.details <- eventPrediction:::CalculateAccrualTimes(Naccrual,Nsim,indat$rand.date,accrualGenerator) 
  
  #calculate quantiles from the recruitment details matrix for storing in output 
  recQuantiles <- eventPrediction:::SimQOutputFromMatrix(rec.details,limit,Nsim)
  
  #subset the recruitment details to get the new subjects        
  newrecs <- if(Naccrual!= 0) rec.details[,(ncol(rec.details)-Naccrual+1):ncol(rec.details)]
  else NULL
  
  #generate the simulation specific parameters
  #e.g. rate and shape Weibull parameters used for each simulation 
  singleSimParamsRgt <- SimParamsRgt@generateParameterFunction( Nsim )
  singleSimParams <- singleSimParamsRgt
  
  for( i in seq_len( length( timeCutList ) ) ){  
    singleSimParamsLft <- SimParamsList[[ i ]]@generateParameterFunction( Nsim )
    singleSimParams <- cbind( singleSimParams, 
                              singleSimParamsLft[,2], 
                              singleSimParamsLft[,3], 
                              timeCutList[[ i ]] )
    colnames( singleSimParams )[(length(singleSimParams[1,])-2):length(singleSimParams[1,])] <- 
      c(  paste0( "rate.lft.", i ),
          paste0( "shape.lft.", i ),
          paste0( "time.cut.", i ) )
  }
  singleSimParams <- cbind( singleSimParams, N=rep( length( timeCutList ), nrow( singleSimParams ) ) )
  
  #perform the simulations    
  outcomes <-apply( singleSimParams, 1, eventPrediction:::PerformOneSimulation,
                    number.subjects=nrow(indat),
                    Naccrual=Naccrual,
                    indat=indat,newrecs=newrecs,HR=HR,r=r,
                    dropout.rate=dropout.rate, 
                    dropout.shape=dropout.shape, followup=data@followup,
                    conditionalFunction=rcpweibull )
  
  #post process the output
  event.type <- sapply( outcomes,function(x){x$event.type})
  if(class(event.type)=="numeric") event.type <- matrix(event.type,ncol=Nsim)
  times <- sapply(outcomes,function(x){x$event.times})
  if(class(times)=="numeric") times <- matrix(times,ncol=Nsim)
  
  #calculate the quantiles for the events and dropouts
  event.times <- t(times)
  eventQuantiles <- eventPrediction:::SimQOutputFromMatrix(event.times,limit,Nsim,event.type=t(event.type),non.inf.event.type=0)
  dropoutQuantiles <- eventPrediction:::SimQOutputFromMatrix(event.times,limit,Nsim,event.type=t(event.type),non.inf.event.type=1)
  
  #use a dummy AccrualGenerator if one not given
  if(is.null(accrualGenerator)) 
    accrualGenerator <-  new("AccrualGenerator",f=function(N){NULL},model="NONE",text="NONE")
  
  return(new("FromDataResults", 
             eventQuantiles = eventQuantiles,  
             recQuantiles=recQuantiles,
             limit = limit,
             event.data = data,
             accrualGenerator=accrualGenerator,
             Naccrual=Naccrual,
             time.pred.data=eventPrediction:::EmptyPredictionDF(),
             event.pred.data=eventPrediction:::EmptyPredictionDF(),
             singleSimDetails=eventPrediction:::SingleSimDetails(event.type=event.type,event.times=times,rec.times=t(rec.details)),
             dropout.shape=dropout.shape,
             dropout.rate=dropout.rate,
             dropoutQuantiles=dropoutQuantiles,
             simParams=SimParamsRgt ## Do we need this?
  ))  
}