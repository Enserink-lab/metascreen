#' Generate a list of rows, columns and wells for a given plate format
#'
#' @description
#' A complementary component of the modular library \pkg{metascreen}, valuable for a wide range of analytical task in a drug sensitivity screen. 
#' \emph{\code{baseplate}} is a function that returns a list of rows, columns and wells associated with a given plate type. 
#' This tool is helpful for cross-referencing plates and plate layouts from different sources, i.e. different plate readers and other machines that use microplates. 
#' It proved to be essential in the process of consolidating data sets from different sources, i.e. dispensing files and raw measurement files with different row, column or well designations. 
#' \emph{baseplate} was designed with the objective in mind to address the incompatibility of different formats.
#'
#' @param plateType \code{integer}; a predefined numeric value designating the plate type based on either 6, 12, 24, 48, 96, 384 or 1536 wells.
#' @param zeroBase \code{logical}; if TRUE, then row and column numbering follow a zero-based numbering.
#'  system. The default is FALSE, in which the row and column numbering corresponds to the actual numbering on the plate.
#' @param leadingZero logical; if TRUE, then numeric column labels will have a leading zero based on a double digit format. The default is FALSE, in which column labeling corresponds to the actual numbering on the plate.
#'
#' @details The value for \emph{plateType} is restricted to numbers that represent the number of wells for a given plate type. Any other number or value will lead to an error message and fail to generate a list.
#' The following plate types are supported: 6, 12, 24, 48, 96, 384 and 1536
#'
#' If \emph{zeroBase} is TRUE, the \strong{numbering} of rows, columns and wells will start with 0, otherwise if \emph{
#' zeroBase} is FALSE, the numbering starts with 1.
#'
#' If \emph{leadingZero} is TRUE, the \strong{labeling} of columns and wells will have a 0 preceding single digit numbers, i.e. A01, instead of A1.
#'
#' @return Returns an object of S3 class "baseplate", with a list of \emph{rows}, \emph{columns} and \emph{wells}.
#' \preformatted{
#'  Attributes:
#'    $names
#'  [1] "no." "id."
#' }
#'
#' Each list contains two objects: \emph{no.}, an object of type "numeric" containing the row, column and well \strong{numbers} and \emph{id.}, an object of type "character" containing the row, column and well \strong{labels}.
#'
#' @references Robert Hanes et al.
#' @author Robert Hanes
#' @note Version: 2021.03
#'
#' @examples
#' baseplate(96)
#' baseplate(384, leadingZero = TRUE)
#' baseplate(1536, zeroBase = TRUE, leadingZero = FALSE)
#'
#' @keywords drug screen well microplate plate
#'
#'
#' @export

baseplate <- function(plateType, zeroBase=FALSE, leadingZero=FALSE){

  # Check, if arguments are in the proper format
  if(!is.numeric(plateType)){stop("Plate type needs to be a numeric object.\nPlease select a plate with 6, 12, 24, 48, 96, 384 or 1536 wells.", call. = TRUE)}
  if(!(plateType %in% c(6, 12, 24, 48, 96, 384, 1536))){stop("Plate type not supported. Please select a plate with 6, 12, 24, 48, 96, 384 or 1536 wells.", call. = TRUE)}


  # Calculate the number of columns and rows based on the plate type
  # The plate matrix is based on columns and rows at a ratio of 3:2
  # This part calculates the closest equal (lowest distance between them) factors of a number in which the modulus
  # between the two factors is zero. Example: 12 = 4 x 3, but not 6 x 2 (distance between 6 and 2 is greater than of 4 and 3)
  .noRow <- as.integer(sqrt(plateType))
  while (plateType %% .noRow != 0) {
    .noRow <- .noRow-1
  }
  .noCol <- plateType / .noRow


  # Check if zero-based numbering should be used
  if(zeroBase == TRUE){ s = 0 }else{ s = 1}

  # Generate row labels for the given plate type and set attributes designating the number attr(,"no.") and name attr(,"id.")
  .listofRows <- c(LETTERS, sapply(LETTERS, function(x) paste0(x, LETTERS)))[1:.noRow]
  attr(.listofRows, "no.") <- seq(s, .noRow-(1-s), 1)
  attr(.listofRows, "id.") <- as.vector(.listofRows)

  # Generate column numbers for the given plate type and set attributes designating the number attr(,"no.") and name attr(,"id.")
  .listofColumns <- seq(1, .noCol, 1)
  attr(.listofColumns, "no.") <- seq(s, .noCol-(1-s), 1)
  attr(.listofColumns, "id.") <- as.vector(.listofColumns)

  # Generate wells from the row and column labels and set attributes designating the number attr(,"no.") and name attr(,"id.")
  .listofWells <- as.character(sapply(.listofRows, FUN = function(x) paste(x, .listofColumns, sep = "")))
  attr(.listofWells, "no.") <- seq(s, plateType-(1-s), 1)
  attr(.listofWells, "id.") <- as.vector(.listofWells)



  # Check if leading zeros should be added to the column numbers and well labels
  if(leadingZero == TRUE){
    # Add leading zero with fixed width to column numbers: 1 to 01, but not 10 to 010
    .listofColumns <- paste0(
      # Extract the well number and add leading zero to one digit numbers
      sprintf("%02d", as.numeric(gsub("[^[:digit:]]", "0\\1", .listofColumns)))
    )
    # Set attributes designating the number and names of columns
    attr(.listofColumns, "no.") <- seq(s, .noCol-(1-s), 1)
    attr(.listofColumns, "id.") <- as.vector(.listofColumns)
    # Add leading zero with fixed width to well numbers: A1 to A01, but not A10 to A010
    .listofWells <- paste0(
      # Extract the well letter and keep it
      gsub("[[:digit:]]", "", .listofWells),
      # Extract the well number and add leading zero to one digit numbers
      sprintf("%02d", as.numeric(gsub("[^[:digit:]]", "0\\1", .listofWells)))
    )
    # Set attributes designating the number and names of wells
    attr(.listofWells, "no.") <- seq(s, plateType-(1-s), 1)
    attr(.listofWells, "id.") <- as.vector(.listofWells)
  }


  # Create an object of S3 class "baseplate", with a list from the attributes of all three vectors
  .baseplate <- list(rows = attributes(.listofRows), columns = attributes(.listofColumns), wells = attributes(.listofWells))


  # Set the class of the object to be returned
  class(.baseplate) <- "baseplate"

  return(.baseplate)

}
