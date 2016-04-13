context("simulatePW")

test_that("cuttingData",{
  data( lag.data )
  my.data <- EventData( data=lag.data,
                        subject="subject",
                        rand.date="randDate",
                        has.event="hasEvent",
                        withdrawn="withdrawn",
                        time=list( event.date="eventDate",
                                   last.date="lastDate" ) )
  
  estLagT <- 200
  right.cens.data <- KMRightCensor( my.data, estLagT )
  expect_equal( nrow( right.cens.data@subject.data ), nrow( my.data@subject.data ) )
  expect_equal( sum( right.cens.data@subject.data$time > estLagT ), 0 )
  expect_true( all( subset( right.cens.data@subject.data, 
                            has.event==F & 
                            withdrawn==F & 
                            rand.date <= max( right.cens.data@subject.data$rand.date ) - estLagT  )$time == estLagT ) )
})

