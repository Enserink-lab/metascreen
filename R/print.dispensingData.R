#' Printing a dispensing data set
#'
#' @description
#' The function \emph{\code{print}} extracts the dispensing data from an object of class \code{dispensingData}.
#'  
#' @param x an object of class 'dispensingData'.
#' @param ... further arguments passed to or from other methods.
#' 
#' @details This function extends the generic function \code{\link{print}} for objects of class 'dispensingData'. 
#' It can be used, whenever the dispensing data is needed for downstream analysis. 
#' It allows to quickly extract the dispensing data from an S3 object of class \code{dispensingData} and return it as a data.frame, without prior knowledge of the data structure.
#' 
#' @seealso \code{\link{generateDispensingData}}
#'
#' @return Returns dispensing data from an object of class 'dispensingData' as a data.frame.
#'
#' @examples
#' \donttest{\dontrun{
#' print(dispensingData)
#' }}
#'
#' @keywords drug screen drug combination dispensing print
#'
#' 
#' @export

print.dispensingData <- function(x, ...){
  
  # Check, if the dispensing data has been provided as an object of class S3:dispensingData
  if(missing(x)){stop("Dispensing data missing! Please provide a dispensing data set.", call. = TRUE)}
  if(class(x) != "dispensingData"){stop("Dispensing data not of class 'dispensingData'!", call. = TRUE)}
  
  # Assign object to internal variable
  dispensingData = x
  
  
  print(dispensingData[["output"]])
  
}