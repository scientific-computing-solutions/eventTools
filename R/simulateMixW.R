##' Simulate method for use when predicting from piecewise Weibull
##' @param object The fitted mixture model
##' @param ... Additional arguments to be passed to the method
##' @rdname SimulateMixW-methods
##' @export 
setGeneric( "SimulateMixW",
            function( object, ... )
              standardGeneric( "SimulateMixW")
)



##' @rdname SimulateMixW-methods
##' @export
setMethod( "SimulateMixW",signature=c( "EventModelBayesian" ),
           function( object, ... ){
             df <- as.data.frame( object@mcmc.object )
             params <- data.frame( 
               df[, grep( "shape", names(df)) ],
               df[, grep( "scale.1", names(df)) ],
               df[, grep( "scale.2", names(df)) ],
               df[, grep( "subgroup*", names(df)) ] )
             names(params) <- c( "shape", 
                                 "scale.1", 
                                 "scale.2", 
                                 paste0( "subject.alloc.", 1:(ncol(params)-3)) )
             simulate.Mixture( object@event.data, params, ... )
           })

##' Only used internally
##' @param data EventData object 
##' @param params Parameters from jags-sampling
##' @param accrualGenerator accrualGenerator object
##' @param Naccrual Number of subjects to recruit
##' @param limit Used to set CI limit
##' @param longlagsettings LonLagSettings object
##' @param dropout Dropout object 
simulate.Mixture <- function( data, params, 
                              accrualGenerator=NULL,Naccrual=0, 
                              limit=0.05, 
                              longlagsettings=NULL,
                              dropout=NULL ){
  message( "!! This code has not been fully reviewed - please report any inconsistencies !!" )
  # THese are not applicable here
  r <- NULL
  HR <- NULL
  seed <- 1234 #Not used (only dummy here)
  
  Nsim <- nrow( params )
  
  #validate the arguments
  eventPrediction:::validate.simulate.arguments(accrualGenerator,Naccrual,Nsim,seed,
                                                limit,longlagsettings,HR,r,data)  
  
  
  
  #calculate the dropout rate and shape for drop out
  dropoutctrlSpec <- eventPrediction:::CtrlSpecFromList(dropout,eventtext="",1)[[1]]
  dropout.shape <- if(is.null(dropout) || is.null(dropout$shape)) 1 else dropout$shape
  dropout.rate <- log(2)^(1/dropout.shape)/dropoutctrlSpec@median
  
  
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

  #perform the simulations    
  outcomes <- lapply( 1:Nsim,  function(i) {
                    PerformOneMixtureSimulation( params[i,],
                      number.subjects=nrow(indat),
                      Naccrual=Naccrual,
                      indat=indat,newrecs=newrecs,HR=HR,r=r,
                      dropout.rate=dropout.rate, 
                      dropout.shape=dropout.shape, followup=data@followup,
                      conditionalFunction=eventPrediction:::rcweibull ) } )
  
  #post process the output
  event.type <- sapply( outcomes,function(x){x$event.type})
  if(class(event.type)=="numeric") event.type <- matrix(event.type,ncol=Nsim)
  times <- sapply(outcomes,function(x){x$event.times})
  if(class(times)=="numeric") times <- matrix(times,ncol=Nsim)
  
  #calculate the quantiles for the events and dropouts
  event.times <- t( times )
  eventQuantiles <- eventPrediction:::SimQOutputFromMatrix(event.times,limit,Nsim,event.type=t(event.type),non.inf.event.type=0)
  dropoutQuantiles <- eventPrediction:::SimQOutputFromMatrix(event.times,limit,Nsim,event.type=t(event.type),non.inf.event.type=1)
  
  #use a dummy AccrualGenerator if one not given
  if(is.null(accrualGenerator)) 
    accrualGenerator <-  new("AccrualGenerator",f=function(N){NULL},model="NONE",text="NONE")
  
  # This is necessary to create the FromDataResult Object. 
  dummy <- eventPrediction::FromDataParam(
    type="weibull",
    rate=1e-15,
    shape=1e-15,
    sigma=matrix(rep(0, 4), nrow=2) )
  
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
             simParams=dummy ## This is never used but required...
  ))  
}


# Function to perform a single simulation
# 
# @param row A vector containing 3 elements first the rate for events hazard, 
# second the shape for the event hazard and third the index of the simulation
# @param number.subjects The number of existing subjects in the data frame
# @param Naccrual The number of additional subjects recruited
# @param indat The subject.data slot of the EventData object (after the long lagsettings have been applied)
# @param newrecs A matrix of recruitment times newrecs[row[[3]],] is a vector of recruitment times
# needed for simulated subjects for the single simulation
# @param HR The hazard ratio: an advanced option which allows two arm trials to be simulated. This replicates the 
# Predict from parameters functionality but uses the recruitment times found in \code{data}. See the vignette for
# further details    
# @param r The allocation ratio: see \code{HR} argument.
# @param dropout.rate Weibull rate parameter for subject drop out
# @param dropout.shape Weibull shape parameter for subject drop out
# @param followup The subject follow up time (Inf if followed until event/withdrawn) 
# @param conditionalFunction The conditionalFunction slot of a \code{FromDataSimParam} object
# @return A list with two elements both vectors of length 
# Naccrual+number.subjects. event.type which contains 0 if event occurs
# 1 if subject dropped out or 2 if subject completed follow up period
# event.times the dates (as numeric) of subjects leaving the trial due to follow up
PerformOneMixtureSimulation <- function(row,number.subjects,Naccrual,
                                 indat, newrecs,HR,r, dropout.rate, 
                                 dropout.shape, followup,conditionalFunction){
  #first set all subjects to having an event
  #0 = event, 1 = dropout, 2 = follow up
  event.type <- rep(0,number.subjects+Naccrual)
  
  
  #And the structure for the time of leaving the study
  retVal <- structure(rep(NA_real_, number.subjects+Naccrual), class="Date")
  
  #subjects who have already left the trial
  existing.events.index=which(indat$has.event==1 | indat$withdrawn==1 | indat$censored.at.follow.up==1)
  
  #set their date and event type correct
  retVal[existing.events.index] <- eventPrediction:::LastDate(indat)[existing.events.index]
  event.type[existing.events.index] <- ifelse(indat$has.event[existing.events.index]==1, 0,
                                              ifelse(indat$censored.at.follow.up[existing.events.index]==1,2,1))
  

  # Divide subjects based on subgroup assignments
  idx <- as.logical( row[4:length(row)] )
  
  #get subjects who have not had an event/withdrawn
  w.1 <- which( !idx & indat$has.event==0 & indat$withdrawn==0 & indat$censored.at.follow.up==0)
  w.2 <- which( idx &  indat$has.event==0 & indat$withdrawn==0 & indat$censored.at.follow.up==0)
 
  lower.bounds.1 <- indat$time[w.1]
  lower.bounds.2 <- indat$time[w.2]
  rand.date.1 <- indat$rand.date[w.1]
  rand.date.2 <- indat$rand.date[w.2]
  
  #add to them accrued subjects
  # Add a fraction of patients to each subgroup based on 
  # the assignments
  if(Naccrual > 0){
    f.rand <- sum( idx ) / length( idx )
    N.1 <- ceiling(  f.rand * Naccrual )
    N.2 <- floor(  (1-f.rand) * Naccrual )
    lower.bounds.1 <- c(lower.bounds.1,rep(0,N.1))
    lower.bounds.2 <- c(lower.bounds.2,rep(0,N.2))
    nw.tmp <- newrecs[as.integer( row[names(row) == "Id"] ),]
    idx.1 <- sample( seq_len(Naccrual), N.1 )
    idx.2 <- !(seq_len(Naccrual) %in% idx.1)
    nw.1 <- nw.tmp[idx.1]
    nw.2 <- nw.tmp[idx.2]
    
    rand.date.1 <- if( N.1 > 1 ) c( rand.date.1,as.Date( nw.1, origin="1970-01-01" ) )
    else c(rand.date.1,as.Date( nw.1, origin="1970-01-01" ) )
    rand.date.2 <- if( N.1 > 1 ) c( rand.date.2,as.Date( nw.2, origin="1970-01-01" ) )
    else c(rand.date.2,as.Date( nw.2, origin="1970-01-01" ) )
    
    w.1 <- c(w.1,(number.subjects+1):(number.subjects+N.1) )
    w.2 <- c(w.2,(number.subjects+N.1+1):(number.subjects+N.1+N.2))
  }
  
  row <- as.data.frame(row)

  # Perhaps not optimally coded but decided to keep it this way due to the large
  # number of parameters. 
  # Simulate times from first mixture model
  HR.1 <- if(is.null(HR)) rep(1,length(lower.bounds.1)) else eventPrediction:::GetHRs(HR,r,length(lower.bounds.1))
  row.1 <- data.frame( shape = row$shape, rate = 1/row$scale.1 )
  times.1 <- conditionalFunction(t.conditional=lower.bounds.1,params=row.1,HR=HR.1) 

  # Simulate times from second mixture model
  HR.2 <- if(is.null(HR)) rep(1,length(lower.bounds.2)) else eventPrediction:::GetHRs(HR,r,length(lower.bounds.2))
  row.2 <- data.frame( shape = row$shape, rate = 1/row$scale.2 )
  times.2 <- conditionalFunction(t.conditional=lower.bounds.2,params=row.2,HR=HR.2) 
      
  
  #if competing risks dropout then deal with this
  if(dropout.rate != 0){
    dropout.times.1 <- eventPrediction:::rcweibull( lower.bounds.1,params=list(shape=dropout.shape,rate=dropout.rate),
                                  HR=rep(1,length(lower.bounds.1)))
    event.type[w.1] <- ifelse(times.1 < dropout.times.1,0,1)
    times.1 <- pmin( times.1, dropout.times.1 )
    
    dropout.times.2 <- eventPrediction:::rcweibull( lower.bounds.2,params=list(shape=dropout.shape,rate=dropout.rate),
                                  HR=rep( 1,length(lower.bounds.2)))
    event.type[w.2] <- ifelse( times.2 < dropout.times.2, 0, 1 )
    times.2 <- pmin( times.2, dropout.times.2 )
  }
  
  #if anytimes are > followup then set them as followup
  lost.to.followup.1 <- which(times.1 > followup)
  times.1[lost.to.followup.1] <- followup
  event.type[w.1[lost.to.followup.1]] <- 2
  
  lost.to.followup.2 <- which(times.2 > followup)
  times.2[lost.to.followup.2] <- followup
  event.type[w.2[lost.to.followup.2]] <- 2
  
  # Put the times into the return vector
  retVal[w.1] <- eventPrediction:::internal.Date(pmax(1,round(times.1)),rand.date.1)   
  retVal[w.2] <- eventPrediction:::internal.Date(pmax(1,round(times.2)),rand.date.2) 
    
  return(list(event.type=event.type,
              event.times=retVal))
  
}