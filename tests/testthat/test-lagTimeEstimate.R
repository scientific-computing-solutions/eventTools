context("simulatePW")

test_that("lagTime",{
  
  data( lag.data )
  my.data <- EventData( data=lag.data,
                        subject="subject",
                        rand.date="randDate",
                        has.event="hasEvent",
                        withdrawn="withdrawn",
                        time=list( event.date="eventDate",
                                   last.date="lastDate" ) )
  
  N.out <- length( seq( 0.5, 20, 0.5 ) )
  est.obj <- estimateLagTime( my.data, t.start=0.5, t.stop=20, dt=0.5 )
  expect_equal( nrow( est.obj@event.data@subject.data ), nrow( my.data@subject.data ) )
  expect_equal( est.obj@criterion, "AIC" )
  expect_equal( length( est.obj@times ), N.out )
  expect_equal( length( est.obj@XIC ), N.out )
  expect_true( is.numeric( getEstimate( est.obj ) ) )
  est.obj <- estimateLagTime( my.data, t.start=0.5, t.stop=20, dt=0.5, crit="BIC" )
  expect_equal( est.obj@criterion, "BIC" )
  expect_equal( length( est.obj@times ), N.out )
  expect_equal( length( est.obj@XIC ), N.out )
  expect_true( is.numeric( getEstimate( est.obj ) ) )
})

