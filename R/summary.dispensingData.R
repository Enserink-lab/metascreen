#' Summarizing a dispensing data set
#'
#' @description
#' The function \emph{\code{summary}} generates a concise summary for a dispensing data set. The function provides insight into the scale of an experimental set-up.  
#' 
#' @param object an object of class 'dispensingData'.
#' @param ... further arguments passed to or from other methods.
#' @method summary dispensingData
#' @S3method summary dispensingData
#' 
#' @details This function extends the generic function \code{\link{summary}} for objects of class 'dispensingData'.
#' It returns the number of unique drug treatments, such as drug combinations (if, a drug combination screen is being carried out) and single drug treatments, as well as the number of drugs and controls used. 
#' It also provides the number of doses, and in cases where the number of doses differs between single drug treatments and drug combinations, each number is shown respectively (sing/comb). 
#' A more important number in the planning of a drug screen is the amount of plates an experimental set-up is leading up to, which will also be estimated by this function.
#'  
#' @seealso \code{\link{generateDispensingData}}
#'
#' @return Returns summary statistics in form of a matrix for a given dispensing data set.
#'
#' @examples
#' \donttest{\dontrun{
#' summary(dispensingData)
#' }}
#'
#' @keywords drug screen drug combination dispensing summary
#'
#'
#' @export

summary.dispensingData <- function(object, ...){
  
  # Check, if the dispensing data has been provided as an object of class S3:dispensingData
  if(missing(object)){stop("Dispensing data missing! Please provide a dispensing data set.", call. = TRUE)}
  if(class(object) != "dispensingData"){stop("Dispensing data not of class 'dispensingData'!", call. = TRUE)}
  
  # Assign object to internal variable
  dispensingData = object
  
  
  # Printing a summary for a given dispensing set
  message('\n', "Summary for dispensing ID: ", dispensingData[["dispensingID"]])
  cat("Number of drug treatments:",               format(dispensingData[["summary"]][["treatments"]], big.mark=" "), fill = TRUE)
  cat("Number of unique drug combinations:",      format(dispensingData[["summary"]][["combinations"]], big.mark=" "), fill = TRUE)
  cat("Number of unique single drug treatments:", format(dispensingData[["summary"]][["single"]], big.mark=" "), fill = TRUE)
  cat("Number of excluded wells per plate:",      format(dispensingData[["summary"]][["excluded"]], big.mark=" "), fill = TRUE)
  cat("Number of controls per plate:",            format(dispensingData[["summary"]][["controls"]], big.mark=" "), fill = TRUE)
  cat("Number of total plates:",                  format(dispensingData[["summary"]][["plates"]], big.mark=" "), fill = TRUE)
  cat("Number of drugs: ", dispensingData[["summary"]][["drugs"]], ",  Number of doses: ", dispensingData[["summary"]][["doses"]][["single"]], "/", dispensingData[["summary"]][["doses"]][["combination"]], sep = "", fill = TRUE)
  
  .summary <- matrix(c(as.numeric(dispensingData[["dispensingID"]]),
                       as.numeric(dispensingData[["summary"]][["treatments"]]),
                       as.numeric(dispensingData[["summary"]][["combinations"]]),
                       as.numeric(dispensingData[["summary"]][["single"]]),
                       as.numeric(dispensingData[["summary"]][["excluded"]]),
                       as.numeric(dispensingData[["summary"]][["controls"]]),
                       as.numeric(dispensingData[["summary"]][["plates"]]),
                       as.numeric(dispensingData[["summary"]][["drugs"]]),
                       as.numeric(dispensingData[["summary"]][["doses"]][["single"]]),
                       as.numeric(dispensingData[["summary"]][["doses"]][["combination"]])), 
                     nrow = 10, ncol = 1, byrow = TRUE, 
                     dimnames = list(c("dispensing ID",
                                       "Drug treatments",
                                       "Drug combinations",
                                       "Single drug treatments",
                                       "Exluded wells p. plate",
                                       "Controls p. plate",
                                       "Number of plates",
                                       "Number of drugs",
                                       "Number of sing. drug doses",
                                       "Number of comb. drug doses"), ""))
  
  invisible(.summary)
  
}
