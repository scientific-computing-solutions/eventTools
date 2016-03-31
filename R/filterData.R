#' This file contains functions to cut data based on the KM curve, 
#' i.e. based on the subject times. 

##' KMRightCensor the Event Data
##' 
##' Creates a new EventData from an existing one where each subject
##' is censored at subject time T (not study time)
##' @param object The EventData object
##' @param time.cut The time in days when to cut the data
##' @rdname KMRightCensor-methods
##' @name KMRightCensor
##' @return \code{EventData} object with the data censored at the appropriate point
##' @export
setGeneric( "KMRightCensor", function( object, time.cut ) 
  standardGeneric("KMRightCensor") )


##' @rdname KMRightCensor-methods
##' @name KMRightCensor
##' @export
setMethod( "KMRightCensor", "EventData", function( object, time.cut ){
  
  subject.data <- object@subject.data
  
  if( time.cut <= 0 ){
    warning( "The subject time (time.cut) must be greater than 0!" )
  }
  
  # No events after time.cut are considered
  subject.data$has.event[ subject.data$time > time.cut ] <- 0

  # No withdrawns after time.cut are considered
  subject.data$withdrawn[ subject.data$time > time.cut] <- 0
  
  # No event.types after time.cut 
  subject.data$event.type[ subject.data$time > time.cut ] <- NA
  
  # No times after time.cut
  subject.data$time[ subject.data$time >= time.cut ] <- time.cut
  
  EventData(
    data=subject.data,
    subject="subject",
    rand.date="rand.date",
    has.event="has.event",
    withdrawn="withdrawn",
    time="time",
    site="site",
    event.type="event.type",
    followup=object@followup
  )
})

# ##' KMLeftTruncate the Event Data
# ##' 
# ##' Creates a new EventData from an existing one where each subject
# ##' is truncated at subject time T (not study time)
# ##' @param object The EventData object
# ##' @param time.cut The time in days when to left-truncate the data
# ##' @rdname KMLeftTruncate-methods
# ##' @name KMLeftTruncate
# ##' @return \code{EventData} object with the data censored at the appropriate point
# ##' @export
# setGeneric( "KMLeftTruncate", function( object, time.cut ) 
#   standardGeneric("KMLeftTruncate") )
# 
# 
# ##' @rdname KMLeftTruncate-methods
# ##' @name KMLeftTruncate
# ##' @export
# setMethod( "KMLeftTruncate", "EventData", function( object, time.cut ){
#   
#   subject.data <- object@subject.data
#   
#   if( time.cut <= 0 ){
#     warning( "The subject time (time.cut) must be greater than 0!" )
#   }
#   
#   subject.data <- subject.data[ subject.data$time > time.cut, ]
#   if( nrow( subject.data ) == 0 ){
#     warning( "No subjects left, please use a lower cut value (time-cut)!" )
#   }  
#   
#   EventData(
#     data=subject.data,
#     subject="subject",
#     rand.date="rand.date",
#     has.event="has.event",
#     withdrawn="withdrawn",
#     time="time",
#     site="site",
#     event.type="event.type",
#     followup=object@followup
#   )
# })
# 



