#' Generate a list of doses for each drug
#'
#' @description
#' An complementary component of the modular library \pkg{metascreen}, invaluable for the initial set-up of a drug sensitivity screen.
#' \emph{\code{generateListofDoses}} is a function that converts a list of doses from wide to long-format.
#' A list of doses is needed for any consecutive steps that facilitate the setup of a drug sensitivity screen as well as the subsequent downstream analysis of data obtained from such a screen.
#'
#' @param data \code{data frame}; a list of doses in wide-format. the data frame includes a multiple columns of doses for a given dose range at which a drug is being tested.
#' @param .doseIdentifier \code{character}; an identifier for columns that contain any dose of a given drug. The default identifier is 'dose'. 
#' @param .dropCol \code{logical}; if TRUE, then non-essential columns are being dropped. The default is FALSE, in which all columns are retained.
#'
#' @details This function is useful for converting a dose list from a wide data format, which is considered more human readable, to a long data format. The long-format is more machine friendly, 
#' with a primary focus on efficiency, rather than readability. The use of repetitive and redundant columns is being avoided, while an unique link between each data attribute (columns) is being maintained for each data entry (row). 
#' 
#' The provided data should contain a column with drug names, one or more columns with doses, and a column with the unit. 
#' While only columns that contain the provided dose identifier will be converted into a long-format, other columns remain untouched. 
#' The data frame provided can contain an unlimited number of columns of doses, however only those that contain the identifier will be used.
#' 
#' If \emph{.doseIdentifier} is not provided, any columns that contained 'dose' in their name will be selected.
#' 
#' The final list of doses should contain the following essential columns: a column with the drug names, a column with the doses and a column with the unit. 
#' Those columns should be labeled as \emph{drug}, \emph{dose} and \emph{unit}, or at least contain these labels in their column name. 
#' Any other columns can be dropped by setting \emph{.dropCol} to TRUE.
#' 
#' @return Returns an object of class "data.frame", with an essential column of ['drug'] names, ['doses'] and ['units'].
#'
#' @references Robert Hanes et al.
#' @author Robert Hanes
#' @note Version: 2021.03
#'
#' @examples
#' \donttest{\dontrun{
#' transformListofDoses(listofDoses, .dropCol = TRUE)
#' }}
#'
#' @keywords drug screen list drugs doses
#'
#'
#' @importFrom stats reshape
#' 
#' @export

generateListofDoses <- function(data, .doseIdentifier = "dose", .dropCol = FALSE){
  
  # Check, if arguments are provided and in the proper format
  if(missing(data)){stop("Data missing! Please provide a list of doses.", call. = TRUE)}

  listofDoses <- data

  # Convert list from wide to long-format
  finalListofDoses <- stats::reshape(listofDoses, varying = list(grepl(.doseIdentifier, names(listofDoses), ignore.case = TRUE)), v.names = c("Dose"), timevar = NULL, direction = "long", new.row.names = NULL)
  finalListofDoses <- within(finalListofDoses, rm("id"))
  rownames(finalListofDoses) <- NULL
  
  # Drop non-essential columns, if requested
  if(.dropCol == TRUE){ finalListofDoses <- finalListofDoses[grepl("Drug|Dose|Unit", names(finalListofDoses), ignore.case = TRUE)]}
  

  # Rearrange essential columns, while leaving other columns unchanged 
  # get the column order and replace matching columns with new order 
  .col <- seq_along(finalListofDoses)
  .col[.col %in% match(c("Drug", "Dose", "Unit"), names(finalListofDoses))] <- match(c("Drug", "Dose", "Unit"), names(finalListofDoses))
  finalListofDoses <- finalListofDoses[.col]
  
  return(finalListofDoses)
  
}
