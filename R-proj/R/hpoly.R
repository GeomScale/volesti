#' The title for my S4 class that extends \code{"character"} class.
#'
#' Some details about this class and my plans for it in the body.
#'
#' \describe{
#'    \item{myslot1}{A logical keeping track of something.}
#'
#'    \item{myslot2}{An integer specifying something else.}
#' 
#'    \item{myslot3}{A data.frame holding some data.}
#'  }
#' @name mynewclass-class
#' @rdname mynewclass-class
#' @exportClass mynewclass
Hpoly <- setClass (
  # Class name
  "Hpoly",
  
  # Defining slot type
  representation (
    A = "matrix",
    b = "numeric",
    volume = "numeric",
    type = "numeric",
    dimension = "numeric"
    ),
  
  # Initializing slots
  prototype = list(
    A = as.matrix(0),
    b = as.numeric(NULL),
    volume = as.numeric(NULL),
    type = as.numeric(NULL),
    dimension = as.numeric(NULL)
  )
)