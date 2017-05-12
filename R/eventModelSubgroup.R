# This file contains the public functions associated with
# the fitted Weibull/loglogistic model for the predict from data simulatisons  
##' @import eventPrediction methods
NULL

##' An S4 class containing a fitted survival model
##' of an Eventdata object
##' @slot model.1 An S3 "survreg" object
##' @slot model.2 An S3 "survreg" object
##' @slot s1 Name of subgroup 1
##' @slot s2 Name of subgroup 2
##' @slot event.data An EventData object used to fit a Weibull survial model
##' @slot simParams The \code{FromDataSimParam} object which contains the information
##' needed to generate the survial times
##' @export
setClass("EventModelSubgroup", 
   slots = list( model.s1="survreg",
                 model.s2="survreg",
                 s1="character",
                 s2="character",
                 event.data="EventData",
                 simParams.s1="FromDataSimParam",
                 simParams.s2="FromDataSimParam")
)


##' Displays information about EventModelSubgroup objects
##' @name show
##' @param object EventModelSubgroup object
##' @rdname show-methods
##' @aliases show,EventModelSubgroup-method
##' @export
setMethod("show",
          "EventModelSubgroup",
  function(object) {
    cat( paste0( "*** Model for subgroup \"", s1, "\" ***\n" ) )
    print( object@model.s1 )
    cat( paste0(  "\n*** Model for subgroup \"", s2, "\" ***\n" ) )
    print( object@model.s2 ) 
})



##' Plots the EventModelSubgroup object
##' @rdname plot-methods
##' @aliases plot,EventModelSubgroup,missing-method
setMethod( "plot", signature( x="EventModelSubgroup", y="missing" ),
  function(x, units="Days", xlab=paste("Time in study [",units,"]",sep=""),
                ylab="", main="", ylim=NULL, xlim=NULL, ...) { 
    xscale <- eventPrediction:::GetScaleForKM(units,daysinyear)
    daysinyear <- standarddaysinyear()
    
    indat <- x@event.data@subject.data
    indat$subgroup <- as.factor( indat$subgroup )
    indat$time <- ifelse(indat$time==0,NA,indat$time)
    subgr <- levels( indat$subgroup )
    s1 <- subgr[1]
    s2 <- subgr[2]
    indat.1 <- indat[ indat$subgroup==s1, ]
    indat.2 <- indat[ indat$subgroup==s2, ]
    
    KM.1 <- survfit(Surv(time, has.event) ~ 1, data=indat.1,...)
    KM.2 <- survfit(Surv(time, has.event) ~ 1, data=indat.2,...)
    

    plot(KM.1, lwd=c(2,1,1), col=c("red", "black", "black" ), 
         xlab=xlab, ylab=ylab,
         xlim=xlim, ylim=ylim,
         main=main,xscale=xscale)
    lines(predict(x@model.s1, type="quantile", p=seq(.01,.99,by=.01))[1,]/xscale, 
          seq(.99,.01,by=-.01), col="brown", type="l", lwd=3)
    
    lines(KM.2, lwd=c(2,1,1), col=c("blue", "black", "black" ))
    lines(predict(x@model.s2, type="quantile", p=seq(.01,.99,by=.01))[1,]/xscale, 
          seq(.99,.01,by=-.01), col="orange", type="l", lwd=3)
    legend("topright", 1, c( "Data", paste0( "Model-", s1 ), paste0( "Model-", s2 ) ), 
           col=c( "red", "brown", "orange" ), lty=c(1,1))
    
  }
)

##' Plots the traces for each parameter of the fitted object
##' @name params
##' @param object EventModelSubgroup object
##' @param ... Additional arguments to be passed to the method
##' @rdname params-methods
##' @aliases params,EventModelSubgroup-method
##' @export
setGeneric( "params", 
            function( object, ... ) 
              standardGeneric("params") )


##' @rdname params-methods
##' @export
setMethod( "params", signature( object="EventModelSubgroup" ),
           function( object ){
             rr <- list() 
             rr[[ object@s1 ]] <- object@simParams.s1@parameters 
             rr[[ object@s2 ]] <- object@simParams.s2@parameters
             rr
           })




