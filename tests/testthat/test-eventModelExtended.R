context("eventModelExtended")

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
  left.trunc.fit <- LeftFit( my.data, estLagT )
  
  expect_equal( class( left.trunc.fit )[1], "EventModelExtended" )
  expect_equal( class( left.trunc.fit@model[[1]] )[1] , "weibreg" )
  expect_equal( left.trunc.fit@time.cut, estLagT )
  expect_equal( nrow( left.trunc.fit@event.data@subject.data ), nrow( my.data@subject.data ) )
  expect_warning( plot( left.trunc.fit ) )
})

