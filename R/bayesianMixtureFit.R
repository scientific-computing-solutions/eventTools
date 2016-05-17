##' @importFrom eventPrediction FromDataParam
NULL

##' An S4 class containing a fitted survival model obtained from JAGS
##' of an Eventdata object
##' @slot shape.init Initial value for shape parameter.
##' @slot shape.min The shape parameter prior is uniformly distributed between 
##' [shape.min, shape.max]. 
##' @slot shape.max The shape parameter prior is uniformly distributed between 
##' [shape.min, shape.max]. 
##' @slot ctrl.median Control median [months]
##' @slot hr Hazard ratio
##' @slot mixture.min The minimum mixture coefficient, e.g. randomization 
##' balance (boundary of uniform distribution), prior is uniformly distributed between 
##' [mixture.min, mixture.max]
##' @slot mixture.max The maximum mixture coefficient, e.g. randomization 
##' balance (boundary of uniform distribution), prior is uniformly distributed between 
##' [mixture.min, mixture.max]
##' @slot scale.sd.min The minimum standard deviation of the log(scale) (0.5), 
##' prior is uniformly distributed between [scale.sd.min, scale.sd.max]
##' @slot scale.sd.max The maximum standard deviation of the log(scale) (2), 
##' prior is uniformly distributed between [scale.sd.min, scale.sd.max]
##' @export
setClass( "PriorAndInitsMixture",
          slots = list( shape.init="numeric",
                        shape.min="numeric",
                        shape.max="numeric",
                        ctrl.median="numeric",
                        hr="numeric",                        
                        mixture.min="numeric",
                        mixture.max="numeric",
                        scale.sd.min="numeric",
                        scale.sd.max="numeric"
          )
)


##' BayesianMixtureFit of the Event Data
##' 
##' Creates an EventModel based on the fitting a Bayesian MCMC mixture model using
##' run.jags and jags.
##' @param object The EventData object
##' @param prior.init Initial values to parameters and prior values, see 
##' \code{PriorAndInitValues}
##' @param ... Additional arguments to be passed to the method
##' @rdname BayesianMixtureFit-methods
##' @name BayesianMixtureFit
##' @return \code{EventData} object with the data censored at the appropriate point
##' @export
setGeneric( "BayesianMixtureFit", 
            function( object, prior.init, ... ) 
  standardGeneric("BayesianMixtureFit") )


##' BayesianMixtureFit the Event Data
##' 
##' @param dist Distribution (only Weibull available for mixture). Can "force" into
##' exponential by specifying a tight boundary around 1 for the shape parameter. 
##' @param N.samples The number of samples to take (2000). This will also determine
##' the number of simulations used in the prediction step.
##' @param N.burnin Number of burnin iterations before sampling starts (1000). Before
##' these 1000 adaptive iterations are tun to enhance the sampling efficiency. 
##' @param N.thin Thining interval to use. May reduce auto-correlation.
##' @param N.chains Number of chins to use. IF running on a parallel machine, one may
##' match this with the number of processors (N.proc).
##' @param N.proc Number of processors to use if running on a cluster.
##' @param parallel Use more than one processor (FALSE)
##' @param seed Random seed, integer (NULL)
##' @rdname BayesianMixtureFit-methods
##' @export
setMethod( "BayesianMixtureFit", 
           signature=c( object="EventData", prior.init="PriorAndInitsMixture" ),
           function( object, 
                     prior.init,
                     dist="weibull", 
                     N.samples=2000, 
                     N.burnin=1000,
                     N.chains=1, 
                     N.thin=1,
                     seed=NULL,
                     parallel=FALSE,
                     N.proc=NULL,
                     ... ){
  if( !dist %in% c( "weibull" ) ){
     stop( "dist must be weibull")
  }
  
  if( nrow( object@subject.data ) == 0 ) stop( "Empty data frame!" ) 

  subject.data <- object@subject.data

  # Censored data must be recorded as NA
  censored = subject.data$has.event==0
  is.censored = censored*1
  times  <- subject.data$time
  
  # Unique version needed for dinterval
  times.cens <- subject.data$time
  times.na   <- subject.data$time
  times.na[censored] = NA
  
  jags.data = list( times = times,
                    times.cens = times.cens,
                    times.na = times.na,
                    is.censored = is.censored,
                    N = nrow(subject.data) )
  
  dayspermonth <- standarddaysinyear()/12
  ctrl.rate <- log(2)^(1/prior.init@shape.init)/prior.init@ctrl.median
  exp.rate  <- log(2)^(1/prior.init@shape.init)/
    (prior.init@ctrl.median/prior.init@hr^(1/prior.init@shape.init))
  ctrl.rate <- ctrl.rate/dayspermonth
  exp.rate  <- exp.rate/dayspermonth
  
  my.model <- "
  model
  {
    for( i in 1:N )
    {
    is.censored[i] ~ dinterval( times.na[i], times.cens[i] )
    subgroup[i] ~ dbern( p )
    times.na[i] ~ dweib( shape, lambda[ 1+subgroup[i] ] )
    surv[i] <- 1 - pweib( times[i], shape, lambda[ 1+subgroup[i] ] )
    }"
    
    my.model <- paste0( my.model,
    "shape ~ dunif( ", prior.init@shape.min ,",", prior.init@shape.max , ")\n",
    "p ~ dunif( ", prior.init@mixture.min, ", ", prior.init@mixture.max, " )\n",
    "log.scale[1] ~ dnorm( log(", ctrl.rate , "), sd1 )\n", 
    "log.scale[2] ~ dnorm( log(", exp.rate, "), sd2 )\n",
    "sd1 ~ dunif(", prior.init@scale.sd.min, ",", prior.init@scale.sd.max,")\n",
    "sd2 ~ dunif(", prior.init@scale.sd.min, ",", prior.init@scale.sd.max,")\n", 
    "for( k in 1:2 ) {\n",
    "  lambda[k] <- pow( scale[k], -shape )\n", 
    "  scale[k] <- exp( log.scale[k] )\n",
    "  mu[k] <- scale[k]*exp(loggam(1+1/(shape)))\n",
    "}\n",
    "}" )
  
  jags.init <- rep( NA, length( times.na ) )
  jags.init[ censored ] <- times.cens[ censored ] + 1
  
  shape <- exp( log( prior.init@shape.min ) + log( prior.init@shape.max ) )
  jags.inits <- function() {
    list( 'shape' = shape, 
          'log.scale' = c( log(ctrl.rate) + rnorm(1), 
                          log(exp.rate) + rnorm(1) ), 
          'times.na' = jags.init )
  }
  
   if( !is.null( seed ) ) {
     set.seed( seed )
     jags.inits <- function() {
       list( 'shape' = shape, 
             'log.scale' = c( max( log(ctrl.rate) + rnorm(1), 0 ), 
                              max( log(exp.rate) + rnorm(1), 0 ) ), 
             'times.na' = jags.init,
             '.RNG.name'= "base::Wichmann-Hill",
             '.RNG.seed'=seed )
     }
   } 
  
  model.fit <- NULL
  if( !parallel ){
    model.fit <- run.jags( data = jags.data,
                           model = my.model,
                           monitor = c( "mu", "scale", "shape", "p", "subgroup" ),
                           sample =  N.samples, 
                           burnin = N.burnin,
                           n.chains =  N.chains,
                           thin = N.thin,
                           summarise = FALSE,
                           inits = jags.inits )
  } else {
    if( is.null( N.proc ) || N.proc <= 1 ) {
      stop( "N.proc needs to be defined and greater than 1!")
    }
    library( parallel )
    cl <- makeCluster( N.chains )
    model.fit <- run.jags( data = jags.data,
                           model = my.model,
                           monitor = c( "mu", "scale", "shape", "p", "subgroup" ),
                           sample = N.samples, 
                           burnin = N.burnin,
                           method="rjparallel",
                           n.chains = N.chains,
                           thin = N.thin,
                           summarise = FALSE,
                           inits = jags.inits,
                           cl = cl )
    stopCluster( cl )
  }
  
  # Retrieve the mcmc object
  x <- suppressWarnings( as.mcmc( model.fit ) ) 
  
  shape.hat = median( x[,names( x[1,] )=="shape"] )
  scale.1.hat = median( x[,names( x[1,] )=="scale[1]"] )
  scale.2.hat = median( x[,names( x[1,] )=="scale[2]"] )
  p.hat = median( x[,names( x[1,] )=="p"] )
  seed = ifelse( !is.null(seed), seed, NA  )
  
  new( "EventModelBayesian",
       model=model.fit,
       event.data=object,
       shape.median = shape.hat,
       scale.1.median = scale.1.hat,
       scale.2.median = scale.2.hat,
       p.median = p.hat,
       mcmc.object = x, 
       seed=as.numeric(seed) )
})




##' PriorAndInitsMixture constructor 
##' 
##' Creates a PriorAndInitsMixture object. These are guesses of the original 
##' parameters based on (protocol) assumptions.
##' @param shape.min The shape parameter prior is uniformly distributed between 
##' [shape.min, shape.max]. 
##' @param shape.max The shape parameter prior is uniformly distributed between 
##' [shape.min, shape.max]. 
##' @param ctrl.median Control median [months]
##' @param hr Hazard ratio
##' @param mixture.min The minimum mixture coefficient, e.g. randomization 
##' balance (boundary of uniform distribution), prior is uniformly distributed between 
##' [mixture.min, mixture.max]
##' @param mixture.max The maximum mixture coefficient, e.g. randomization 
##' balance (boundary of uniform distribution), prior is uniformly distributed between 
##' [mixture.min, mixture.max]
##' @param scale.sd.min The minimum standard deviation of the log(scale) (0.5), 
##' prior is uniformly distributed between [scale.sd.min, scale.sd.max]
##' @param scale.sd.max The maximum standard deviation of the log(scale) (2), 
##' prior is uniformly distributed between [scale.sd.min, scale.sd.max]
##' @export
PriorAndInitValues <- function( shape.init,
                                shape.min,
                                shape.max,
                                ctrl.median,
                                hr,                        
                                mixture.min,
                                mixture.max,
                                scale.sd.min = 0.9,
                                scale.sd.max = 1.1 ){
  if( shape.min > shape.max || shape.min < 0 ) stop( "Invalid range of shape-parameter!" )
  if( mixture.min > mixture.max || mixture.min < 0 ) stop( "Invalid range of mixture-coefficient!" )
  if( ctrl.median < 0 ) stop( "Invalid control median!" )
  if( hr < 0 ) stop( "Invalid hazard ratio!" )
  if( scale.sd.min > scale.sd.max || scale.sd.min < 0 ) stop( "Invalid range of shape-parameter!" )
  
  new( "PriorAndInitsMixture",
       shape.init=shape.init,
       shape.min=shape.min,
       shape.max=shape.max,
       ctrl.median=ctrl.median,
       hr=hr,                        
       mixture.min=mixture.min,
       mixture.max=mixture.max,
       scale.sd.min=scale.sd.min,
       scale.sd.max=scale.sd.max)
}



