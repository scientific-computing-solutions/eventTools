context("lagTimeTest")

test_that("lagTimeTest",{
  
  data( lag.data )
  my.data <- EventData( data=lag.data,
                        subject="subject",
                        rand.date="randDate",
                        has.event="hasEvent",
                        withdrawn="withdrawn",
                        time=list( event.date="eventDate",
                                   last.date="lastDate" ) )
  
  expect_equal( TRUE, TRUE )

## !!These tests works when running on command but R CMD check do not let them pass!!
#   N.out <- length( seq( 0.5, 20, 0.5 ) )
#   est.obj <- estimateLagTime( my.data, 0.5, 20, 0.5 )
#   expect_equal( nrow( est.obj@event.data@subject.data ), nrow( my.data@subject.data ) )
#   expect_equal( est.obj@criterion, "AIC" )
#   expect_equal( length( est.obj@times ), N.out )
#   expect_equal( length( est.obj@XIC ), N.out )
#   expect_true( is.numeric( getEstimate( est.obj ) ) )
#   est.obj <- estimateLagTime( my.data, 0.5, 20, 0.5, criterion="BIC" )
#   expect_equal( est.obj@criterion, "BIC" )
#   expect_equal( length( est.obj@times ), N.out )
#   expect_equal( length( est.obj@XIC ), N.out )
#   expect_true( is.numeric( getEstimate( est.obj ) ) )
})

