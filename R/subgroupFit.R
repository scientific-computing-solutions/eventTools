##' @importFrom eventPrediction FromDataParam
NULL

##' SubgroupFit the Event Data
##' 
##' Creates an EventModel based on the left truncated data
##' @param object The EventData object
##' @param ... Additional arguments to be passed to the method
##' @rdname SubgroupFit-methods
##' @name SubgroupFit
##' @return \code{EventData} object with the data censored at the appropriate point
##' @export
setGeneric( "SubgroupFit", 
            function( object, time.cut, ... ) 
              standardGeneric("SubgroupFit") )


##' SubgroupFit the Event Data
##' 
##' @param object The EventData object
##' @param time.cut The time in days when to cut the data
##' @param dist Distribution (weibull)
##' @param ... Additional arguments to be passed to the method
##' @export
setMethod( "SubgroupFit", signature=c( object="EventData" ),
           function( object, 
                     dist.1="weibull", dist.2="weibull", ... ){
             if(!dist.1 %in% c("weibull","loglogistic")){
               stop("dist.1 must be weibull or loglogistic")
             }
             if(!dist.2 %in% c("weibull","loglogistic")){
               stop("dist.2 must be weibull or loglogistic")
             }
             
             if(nrow(object@subject.data)==0)stop("Empty data frame!") 
             if(sum(object@subject.data$has.event)==0){
               stop("Cannot fit a model to a dataset with no events")
             }
             
             #subjects with time = 0 are set to NA for the model fitting
             #so they are ignored
             indat <- object@subject.data
             indat$subgroup <- as.factor( indat$subgroup )
             indat$time <- ifelse(indat$time==0,NA,indat$time)
             subgr <- levels( indat$subgroup )
             s1 <- subgr[1]
             s2 <- subgr[2]
             indat.1 <- indat[ indat$subgroup==s1, ]
             indat.2 <- indat[ indat$subgroup==s2, ]
             
             if(nrow(indat.1)==0)stop( paste0( "Second subgroup empty (", s1, ")" ) ) 
             if(sum(indat.1$has.event)==0){
               stop(paste0( "Cannot fit a model to a dataset with no events (", s1, ")" ) )
             }
             if(nrow(indat.2)==0)stop( paste0( "Second subgroup empty (", s2, ")" ) ) 
             if(sum(indat.2$has.event)==0){
               stop(paste0( "Cannot fit a model to a dataset with no events (", s2, ")" ) )
             }
             
             model.1 <-survreg(Surv(time, has.event) ~ 1, data=indat.1, dist=dist.1, y = TRUE)
             model.2 <-survreg(Surv(time, has.event) ~ 1, data=indat.2, dist=dist.2, y = TRUE)
             new("EventModelSubgroup",
                 model.s1=model.1,
                 model.s2=model.2,
                 s1=s1,
                 s2=s2,
                 event.data=object,
                 simParams.s1=eventPrediction:::FromDataParam(object=model.1,type=dist.1),
                 simParams.s2=eventPrediction:::FromDataParam(object=model.2,type=dist.2))
           })




