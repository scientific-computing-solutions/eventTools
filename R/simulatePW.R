##' Conditional function for piece wise model 
##' @param params paramter object with study settings
##' @param hr HR
rcpwweibull <- function( t.conditional, params, HR ) {
  params <- as.list( params )
  
  idx <- which( t.conditional < params$time.cut )
  
  # HR is a vector why rate becomes a vector, not 
  # sure why this is necessary, need to check!
  params$rate <- params$rate*(HR)^{1/params$shape}
  params$rate.lft <- params$rate.lft*(HR)^{1/params$shape}
  
  my.times <- rep( NA, length( t.conditional ) )
  my.times[ idx ] <- t.conditional[ idx ]*params$rate[ idx ]
  my.times[ idx ] <- ((my.times[ idx ]^params$shape + 
                         rexp(length(my.times[ idx ])))^(1/params$shape))/params$rate[ idx ]

  # All times generated from model 1 that exceeds the cutoff needs 
  # to be resimulated with model 2 conditioned on time.cut
  idx2 <- which( my.times > params$time.cut )

  # All times simulated that exceeds time.cut will be set to time.cut and
  # simulated from the left-truncated model
  t.conditional[ idx2 ] = params$time.cut 
  my.times[ idx2 ] = NA

  # Simulate all with undefined times
  idx3 <- which( is.na( my.times ) ) 
  
  my.times[ idx3 ] <- t.conditional[ idx3 ]*params$rate.lft[ idx3 ]
  my.times[ idx3 ] <- ((my.times[ idx3 ]^params$shape.lft + 
                          rexp(length(my.times[ idx3 ])))^(1/params$shape.lft ))/params$rate.lft[ idx3 ]
  my.times 
}


##' @export
setGeneric( "simulatePW",
  function( object, object2, ... )
           standardGeneric( "simulatePW")
)


##' @export
setMethod( "simulatePW",signature=c( "EventModel", "EventModelExtended" ),
  function( object, object2, ... ){
    # Currently only Weibull allowed
    if( object@simParams@type != "weibull" ||
        object2@simParams@type != "weibull" ){
      stop( "Only weibull is currenlty allowed for the piecewise simulation!")
    }
  simulate.Internal( object2@event.data, object@simParams, object2@simParams, object2@time.cut, ... )
})

simulate.Internal <- function( data, SimParamsRgt, SimParamsLft, time.cut, 
                               accrualGenerator=NULL,Naccrual=0, Nsim=1e4, 
                               seed=NULL, limit=0.05, 
                               longlagsettings=NULL,HR=NULL,r=NULL,
                               dropout=NULL){
  #validate the arguments
  eventPrediction:::validate.simulate.arguments(accrualGenerator,Naccrual,Nsim,seed,
                              limit,longlagsettings,HR,r,data)  


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
  singleSimParamsLft   <- SimParamsLft@generateParameterFunction( Nsim )
  singleSimParams <- cbind( singleSimParamsRgt, rate.lft=singleSimParamsLft[,2], shape.lft=singleSimParamsLft[,3], 
                            time.cut = time.cut )
  
  #perform the simulations    
  outcomes <-apply( singleSimParams, 1, eventPrediction:::PerformOneSimulation,
                    number.subjects=nrow(indat),
                    Naccrual=Naccrual,
                    indat=indat,newrecs=newrecs,HR=HR,r=r,
                    dropout.rate=dropout.rate, 
                    dropout.shape=dropout.shape, followup=data@followup,
                    conditionalFunction= rcpwweibull )
  
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