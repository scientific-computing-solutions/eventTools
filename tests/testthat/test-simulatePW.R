context("simulatePW")

test_that("validateArgs",{
  data( lag.data )
  my.data <- EventData( data=lag.data,
                        subject="subject",
                        rand.date="randDate",
                        has.event="hasEvent",
                        withdrawn="withdrawn",
                        time=list( event.date="eventDate",
                                   last.date="lastDate" ) )
  
  estLagT <- 160 
  my.fit <- fit( my.data )
  right.cens.data <- KMRightCensor( my.data, estLagT )
  right.cens.fit <- fit( right.cens.data )
  left.trunc.fit <- LeftFit( my.data, estLagT )
  
  expect_error( simulatePW( right.cens.fit, left.trunc.fit, Nsim=-6 ) )
  expect_error( simulatePW( right.cens.fit, left.trunc.fit, Nsim=0 ) )
  expect_error( simulatePW( right.cens.fit, left.trunc.fit, Nsim="sd" ) )
  expect_error( expect_warning( simulatePW( right.cens.fit, left.trunc.fit, Nsim="we" ) ) )
  expect_error( simulatePW( simulatePW( left.trunc.fit, right.cens.fit, Nsim=10 ) ) )
  expect_error( simulatePW( my.fit, Nsim=10 ) ) 
  expect_error( simulatePW( right.cens.fit, left.trunc.fit, Nsim=10, fix.shape=1.0 ) )
  expect_error( simulatePW( right.cens.fit, left.trunc.fit, Nsim=10, Naccrual=3 ) )
  
  ag <- Generate.PoissonAccrual(start.date="2016-04-04",rate=1)
  expect_warning( simulatePW( right.cens.fit, left.trunc.fit, accrualGenerator=ag, Nsim=10 ) )
  expect_warning( simulatePW( right.cens.fit, left.trunc.fit, accrualGenerator=ag, Nsim=10, Naccrual=0 ) )
  expect_error( simulatePW( right.cens.fit, left.trunc.fit, accrualGenerator=g, Nsim=10, Naccrual=10) )
  expect_error( simulatePW( right.cens.fit, left.trunc.fit, accrualGenerator=ag, Nsim=10, Naccrual=-4) )
  
  expect_error( simulatePW( right.cens.fit, left.trunc.fit, Nsim=10,longlagsettings="err") )
  expect_error( simulatePW( right.cens.fit, left.trunc.fit, Nsim=10,longlagsettings=ag ) )
})



test_that("deterministic_bits",{
  set.seed( 11 )
  
  data( lag.data )
  e <- lag.data[lag.data$hasEvent==1,]
  e <- rbind( e,lag.data[700,] )

  NROWS = nrow(e)
  
  my.data <- EventData(data=e,
                       subject="subject",
                       rand.date="randDate",
                       has.event="hasEvent",
                       withdrawn="withdrawn",
                       time=list( event.date="eventDate",
                                  last.date="lastDate" )) 
  
  estLagT <- 160 
  right.cens.data <- KMRightCensor( my.data, estLagT )
  right.cens.fit <- fit( right.cens.data )
  left.trunc.fit <- LeftFit( my.data, estLagT )
  results <- simulatePW( right.cens.fit, left.trunc.fit, Nsim=50, limit=0.25 )
  
  
  expect_equal(NROWS,length(results@eventQuantiles@median))
  expect_equal(NROWS,length(results@eventQuantiles@upper))
  expect_equal(NROWS,length(results@eventQuantiles@lower))
  expect_equal(NROWS,length(results@recQuantiles@median))
  expect_equal(NROWS,length(results@recQuantiles@upper))
  expect_equal(NROWS,length(results@recQuantiles@lower))
  expect_equal(my.data,results@event.data)
  
  expect_equal(0.25,results@limit)
  expect_equal(0,results@Naccrual)
  
  expect_equal(sort(my.data@subject.data$rand.date),results@recQuantiles@median)
  expect_equal(sort(my.data@subject.data$rand.date),results@recQuantiles@upper)
  expect_equal(sort(my.data@subject.data$rand.date),results@recQuantiles@lower)
  expect_equal(sort(as.Date(e$eventDate[1:(NROWS-1)])),results@eventQuantiles@median[1:(NROWS-1)])
  expect_equal(sort(as.Date(e$eventDate[1:(NROWS-1)])),results@eventQuantiles@lower[1:(NROWS-1)])
  expect_equal(sort(as.Date(e$eventDate[1:(NROWS-1)])),results@eventQuantiles@upper[1:(NROWS-1)])

  my.accrual <- Generate.PoissonAccrual(start.date="2015-01-01",rate=1)
  
  expect_warning(results <- simulatePW( right.cens.fit, left.trunc.fit, Nsim=50,
                                        limit=0.25,
                                        accrualGenerator=my.accrual, 
                                        Naccrual=99 )) # (Recruiting in the past)
  
  expect_equal( NROWS+99, length( results@eventQuantiles@median ) )
  expect_equal( NROWS+99, length(results@eventQuantiles@lower))
  expect_equal( NROWS+99, length(results@eventQuantiles@upper))
  expect_equal( NROWS+99,length(results@recQuantiles@median))
  
  expect_equal(results@eventQuantiles@median,sort(results@eventQuantiles@median))
  expect_equal(results@recQuantiles@median,sort(results@recQuantiles@median))
  expect_true(all(results@eventQuantiles@median <= results@eventQuantiles@upper))
  expect_true(all(results@eventQuantiles@median >= results@eventQuantiles@lower))
})


test_that("conditional_function",{
 params <- list( rate = log(2)/3,
                 shape = 2, 
                 rate.lft.1 = log(2)/3,
                 shape.lft.1 = 2, 
                 time.cut=12.01, 
                 N=1 )
 t.conditional <- seq( 1, 12, 1 )
 
 # Test that times generated before cut-point equals rcweibull
 set.seed(10)
 rcpw <- eventTools:::rcpweibull( t.conditional, params, 1 )
 set.seed(10)
 rcw <- eventPrediction:::rcweibull( t.conditional, params, 1 )
 idx <- which( rcpw < 12 )
 expect_equal( rcw[idx], rcpw[idx] )
 
 # Test that times generated after cut-point equals rcweibull (left defined with same parameters)
 t.conditional <- seq( 13, 26, 1 )
 set.seed(10)
 rcpw <- eventTools:::rcpweibull( t.conditional, params, 1 )
 set.seed(10)
 rcw <- eventPrediction:::rcweibull( t.conditional, params, 1 )
 expect_equal( rcw, rcpw )
})

