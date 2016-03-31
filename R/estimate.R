# Functions to estimate the cut-point in time when HR changes


# From: http://stackoverflow.com/questions/15874214/piecewise-function-fitting-with-nls-in-r
findPiecewiseLinearCut <- function( x, y ) {
  f <- function( Cx ) 
  {
    lhs <- function(x) ifelse(x < Cx,Cx-x,0)
    rhs <- function(x) ifelse(x < Cx,0,x-Cx)
    
    fit <- tryCatch( lm(y ~ lhs(x) + rhs(x)), error=function(e) NULL )
    if( !is.null( fit )) {
      c(summary(fit)$r.squared, 
        summary(fit)$coef[1], 
        summary(fit)$coef[2],
        summary(fit)$coef[3])
    }
    else{
      c( NA, NA, NA, NA )
    }
  }
  
  r2 <- function(x) -(f(x)[1])
  
  res <- optimize(r2,interval=c(min(x),max(x)))
  c(res$minimum,f(res$minimum))
}

##' @export
setGeneric( "piecewiseLinearEstimateT",
            function( object, ... )
              standardGeneric( "piecewiseLinearEstimateT")
)

##' @export
setMethod( "piecewiseLinearEstimateT",signature=c( "EventData" ),
           function( object, ... ) {
  KM <- survfit(Surv( time, has.event) ~ 1, data=object@subject.data)
  x <- log( KM$time )
  y <- log( -log( KM$surv ) )
  
  
  res <- findPiecewiseLinearCut( x, y )  
  
  best_Cx <- res[1]
  coef1 <- res[3]
  coef2 <- res[4]
  coef3 <- res[5]
  plot( x, y, main=paste0( "log(T.hat)=", round( best_Cx, 2 ), 
                           ", T.hat=", round( exp(best_Cx), 2 ) ) )
  abline( coef1+best_Cx*coef2,-coef2, col="red", lwd=3 ) #lhs  
  abline( coef1-best_Cx*coef3, coef3, col="blue", lwd=3 )  #rs
  abline( v=best_Cx, lty=2, lwd=3 )
})


setGeneric( "piecewiseLinearBootT",
            function( object, N, ... )
              standardGeneric( "piecewiseLinearBootT")
)

##' @export
setMethod( "piecewiseLinearBootT",signature=c( "EventData" ),
           function( object, N=100, ... ) {
  best_Cx <- rep( NA, N )          
  for( i in 1:N ){
    idx <- sample( seq_len(nrow(object@subject.data)), replace=T ) 
    boot.sample <-object@subject.data[ idx, ] 
    
    KM <- survfit(Surv( time, has.event) ~ 1, data=boot.sample )
    x <- log( KM$time )
    y <- log( -log( KM$surv ) )
    suppressWarnings(
      res <- findPiecewiseLinearCut( x, y )  
    )
    best_Cx[i] <- res[1]
  }
  best_Cx
})
