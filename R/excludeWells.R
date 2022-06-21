#' Exclude wells on 6-, 12-, 24-, 48-, 96-, 384-, and 1536-well microplates from an experimental setup
#'
#' @description
#' A complementary component of the modular library \pkg{metascreen}, valuable for a wide range of analytical task in a drug sensitivity screen.
#' \emph{\code{excludeWells}} is a helpful function that returns a list of wells to be excluded from a 6, 12, 24, 48, 96, 384 or 1536-well microplate in an experimental setup of a high-throughput single or drug combination sensitivity screen (DSS).
#'
#' @param plateType \code{integer}; a predefined numeric value designating the plate type based on either 6, 12, 24, 48, 96, 384 or 1536 wells.
#' @param wells \code{vector}; a set of wells, row or column designations. Wells and rows need to be an object of type "character", while columns need to be "numeric".
#' @param outer.wells \code{logical}; if TRUE, the function will pre-select the outer wells for exclusion.
#'
#' @details The value for \emph{plateType} is restricted to numbers that represent the number of wells for a given plate type. Any other number or value will lead to an error message and fail to generate a list.
#' The following plate types are supported: 6, 12, 24, 48, 96, 384 and 1536
#'
#' If \emph{wells} is provided, a combination of designated rows, columns and individual wells will be marked for exclusion. Both, rows and columns can be specified by their names, which are  designated letters for rows (A, B, C, ..) and designated numbers for columns (1, 2, 3, ..). Alternatively, in addition to entire rows and columns, individual wells (A1, B2, C3, ..) can be marked for exclusion.
#'
#' If \emph{outer.wells} is TRUE, all outer wells of given plate type will be marked for exclusion. Note: allowing outer.wells to be TRUE for a 6-well plate will consequently mark all wells for exclusion on that plate!
#'
#' Along with the \emph{plateType} at least one additional argument needs to be provided, either \emph{wells} and/or \emph{outer.wells}.
#' Designations that are not found for a given plate type will be ignored and produce a warning.
#'
#' @return A vector of wells on a given microplate to be excluded from an experimental setup.
#'
#' Additionally, the plate layout for the given plate is plotted on screen, while highlighting all the excluded wells.
#'
#' @references Robert Hanes et al.
#' @author Robert Hanes
#' @note Version: 2021.02
#'
#' @examples
#' excludeWells(96, outer.wells = TRUE)
#' excludeWells(384, wells = c("A", "E", "M", 8, 12, "H3"))
#' excludeWells(1536, wells = c("AE", 12, "H3"), outer.wells = TRUE)
#'
#' @keywords drug screen well microplate plate
#'
#'
#' @importFrom ggplot2 ggplot aes alpha element_line element_rect element_text geom_hline geom_raster geom_text geom_tile geom_vline scale_fill_manual scale_x_discrete scale_y_discrete coord_fixed theme xlab ylab
#' @importFrom utils head tail
#'
#' @export

excludeWells <- function(plateType, wells=NULL, outer.wells=FALSE){

  # Check, if arguments are in the proper format
  if(!is.numeric(plateType)){stop("Plate type needs to be a numeric object.\nPlease select a plate with 6, 12, 24, 48, 96, 384 or 1536 wells.", call. = TRUE)}
  if(!(plateType %in% c(6, 12, 24, 48, 96, 384, 1536))){stop("Plate type not supported. Please select a plate with 6, 12, 24, 48, 96, 384 or 1536 wells.", call. = TRUE)}
  if(!is.null(wells) & !is.vector(wells)){stop("Wells need to be provided as a vector.", call. = TRUE)}
  if(is.null(wells) & outer.wells == FALSE){return(message("No wells have been selected! Please validate your input."))}


  # Calculate the number of columns and rows based on the plate type
  # The plate matrix is based on columns and rows at a ratio of 3:2
  # This part calculates the closest equal (lowest distance between them) factors of a number in which the modulus
  # between the two factors is zero. Example: 12 = 4 x 3, but not 6 x 2 (distance between 6 and 2 is greater than of 4 and 3)
  .noRow <- as.integer(sqrt(plateType))
  while (plateType %% .noRow != 0) {
    .noRow <- .noRow-1
  }
  .noCol <- plateType / .noRow


  # Define row and column names based on the plate type
  .listofRows <- c(LETTERS, sapply(LETTERS, function(x) paste0(x, LETTERS)))[1:.noRow]
  .listofColumns <- seq(1, .noCol, 1)
  .listofWells <- as.character(sapply(.listofRows, FUN = function(x) paste(x, .listofColumns, sep = "")))



  # Start selecting wells based on user input
  # Select outer wells based on plate type
  if(outer.wells == TRUE){

    # Select outer wells from the first and last column and row
    .ouw <- unique(c(as.character(sapply(c(head(.listofRows, 1), tail(.listofRows, 1)), FUN = function(x) paste(x, 1:length(.listofColumns), sep = ""))),
                     as.character(sapply(c(head(.listofColumns, 1), tail(.listofColumns, 1)), FUN = function(x) paste(.listofRows, x, sep = "")))))

  }else{.ouw = NULL}


  # Select rows, columns and individual wells based on plate type
  if(!is.null(wells)){

    # Remove empty strings and potential NAs from provided wells
    # wells <- na.omit(wells)
    wells <- wells[!wells %in% c("", NA)]

    # Select for only rows, columns and wells
    .rows <- wells[!grepl('[^A-Za-z]', wells)]
    .cols <- wells[grepl('^[[:digit:]]', wells)]
    .wells <- wells[grepl('[A-Za-z][0-9]', wells)]
    # .cols <- wells[!grepl('[^0-9]', wells)]

    # Prompt a warning for rows, columns and wells that are not found for a given plate type
    if(TRUE %in% !(.rows %in% .listofRows) | TRUE %in% !(.cols %in% .listofColumns) | TRUE %in% !(.wells %in% .listofWells)){warning(paste(c(
    if(TRUE %in% !(.rows %in% .listofRows)){paste("Invalid rows have been specified for a ", plateType, "-well plate. ", "The following rows have been omited: ", paste(.rows[!(.rows %in% .listofRows)], collapse = ", "), sep = "")},
    if(TRUE %in% !(.cols %in% .listofColumns)){paste("Invalid columns have been specified for a ", plateType, "-well plate. ", "The following columns have been omited: ", paste(.cols[!(.cols %in% .listofColumns)], collapse = ", "), sep = "")},
    if(TRUE %in% !(.wells %in% .listofWells)){paste("Invalid wells have been specified for a ", plateType, "-well plate. ", "The following wells have been omited: ", paste(.wells[!(.wells %in% .listofWells)], collapse = ", "), sep = "")}
    ), collapse = "\n"), call. = FALSE, noBreaks. = TRUE)}


    # Generate a list of wells based on the given rows, columns and individual wells
    .exrw <- as.character(sapply(.rows[.rows %in% .listofRows], FUN = function(x) paste(x, .listofColumns, sep = "")))
    .excw <- as.character(sapply(.cols[.cols %in% .listofColumns], FUN = function(x) paste(.listofRows, x, sep = "")))
    .exw <-  as.character(.wells[.wells %in% .listofWells])

  }else{.exrw <- .excw <- .exw <- NULL}

  # Combine all wells together and remove duplicate wells
  .exwells <- unique(c(.ouw, .exrw, .excw, .exw))

  # Sort wells according to their well number
  .exwells <- .listofWells[sort(match(.exwells, .listofWells))]


  g <- ggplot2::ggplot(data.frame(Wells = .listofWells, Rows = gsub("([0-9])", "", .listofWells), Columns = gsub("([A-Z])", "", .listofWells), Excluded = .exwells[match(.listofWells, .exwells)]), aes(x = factor(gsub("([A-Z])", "", .listofWells), .listofColumns), y = factor(gsub("([0-9])", "", .listofWells), rev(.listofRows)))) +
    geom_raster(aes(fill = Excluded), na.rm=TRUE) +
    geom_tile(aes(fill = Excluded, width=1, height=1), color="lightgrey", linetype = "dotted", size=0.3) +
    scale_fill_manual(values = rep("#FF0066", 1, length(.exwells)), na.value = NA) +
    # scale_fill_gradient(low = "#FF0066", high = "#FF0066", na.value = NA) +
    # scale_fill_hue("Drugs", l=60, h = c(350, 350), na.value = NA) +
    geom_text(aes(label = paste(ifelse(is.na(.exwells[match(.listofWells, .exwells)]),"",.exwells[match(.listofWells, .exwells)]), sep = "\n")), fontface = "plain", color ="white", hjust = 0.5, size = 2.6) +
    geom_hline(yintercept=seq(1.5, length(unique(gsub("([0-9])", "", .listofWells)))-0.5, 1), lwd=0.2, colour="#FFFFFF", linetype = "dotted") +
    geom_vline(xintercept=seq(1.5, length(unique(gsub("([A-Z])", "", .listofWells)))-0.5, 1), lwd=0.2, colour="#FFFFFF", linetype = "dotted") +
    #scale_x_discrete(position = "top") +
    scale_x_discrete(expand = c(0,0), position = "top") +
    scale_y_discrete(expand = c(0,0)) +
    # Set the ratio of the plot reflecting the ratio of the plate: 3:2
    coord_fixed(ratio = 2/3) +
    xlab(NULL) +
    ylab(NULL) +
    # labs(caption = paste("Dispensing Plate: ", destinationPlateID, ".", plateNumber, sep = "")) +
    theme(plot.title = element_text(hjust = 0.5),
          plot.caption = element_text(size = 7),
          text = element_text(size = 12),
          axis.line = element_line(color='black'),
          # panel.grid.major = element_line(color="lightgrey", linetype = "dotted"),
          # panel.grid.minor = element_line(color="lightgrey", linetype = "dotted"),
          panel.background = element_rect(fill=alpha('white', 1)),
          panel.border = element_rect(colour = "black", fill=NA, size=1),
          legend.position = "none",
          # Set the ratio of the plot reflecting the ratio of the plate: 3:2
          # aspect.ratio = 2/3
    )

  print(g)

  # Set the class of the object to be returned
  class(.exwells) <- "exwells"

  return(.exwells)

}
