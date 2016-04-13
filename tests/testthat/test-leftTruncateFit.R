context("simulatePW")

test_that("leftFit",{
  
  # Test that estimate works for simulated data
  N <- 1e5
  shape <- 2 
  scale <- 500
  
  set.seed( 20160413 )
  event.data <- data.frame( 
    subject=paste0( "Subj_", 1:N ),
    rand.date=rep( Sys.Date(), N ),
    withdrawn=rep( FALSE, N ),
    has.event=c( rep( TRUE, N/2 ), rep( FALSE, N/2 ) ),
    time=rweibull( N, shape=shape, scale=scale )
    )
    
  my.data <- EventData( data=event.data, 
                        subject="subject", 
                        rand.date="rand.date", 
                        has.event="has.event", 
                        withdrawn="withdrawn", 
                        time="time" )
  
  my.fit <- fit( my.data )
  
  estLagT <- 200
  left.trunc.fit <- LeftFit( my.data, estLagT )
  expect_true( ( left.trunc.fit@simParams@parameters$rate - my.fit@simParams@parameters$rate )/
                 left.trunc.fit@simParams@parameters$rate < 0.01 )
  expect_true( left.trunc.fit@simParams@parameters$shape - my.fit@simParams@parameters$shape < 0.01 )
})

