#' Plotting a set of dispensing plates
#'
#' @description
#' The function \emph{\code{plot}} generates individual plots for each plate of a given dispensing data set. 
#' 
#' @param x an object of class 'dispensingData'.
#' @param .saveto string; path to a folder location where the object is saved to.
#' @param ... further arguments passed to or from other methods.
#' 
#' @details This function extends the generic function \code{\link{plot}} for objects of class 'dispensingData'. 
#' It is used to visualize the dispensing data of each \emph{destination plate} from a given dispensing set.
#' The individual drug treatments are arranged according to the dispensing instructions and portrayed in the corresponding well-plate format.
#' The plots will be save to the provided folder location. If the folder does not exist, it will be created, if the folder location is not provided, 
#' the plots will be saved to the current working directory.
#' 
#' @seealso \code{\link{generateDispensingData}}
#'
#' @examples
#' \donttest{\dontrun{
#' plot(dispensingData, .saveto = "../myDispensing/plots")
#' }}
#'
#' @keywords drug screen drug combination dispensing plot
#'
#'
#' @importFrom ggplot2 ggplot ggsave labs scale_fill_hue
#' @importFrom stats aggregate
#' 
#' @export

plot.dispensingData <- function(x, .saveto, ...){
  
  # Check, if the dispensing data has been provided as an object of class S3:dispensingData
  if(missing(x)){stop("Dispensing data missing! Please provide a dispensing data set.", call. = TRUE)}
  if(class(x) != "dispensingData"){stop("Dispensing data not of class 'dispensingData'!", call. = TRUE)}
  
  # Check, if a folder location has been provided
  if(missing(.saveto)){ .saveto <- getwd() }
  # Create folder, if it does not exist
  if(!file.exists(file.path(.saveto))){ dir.create(file.path(.saveto), showWarnings = FALSE, recursive = TRUE) }
  if(!utils::file_test("-d", .saveto)){stop("in '.saveto'. Argument needs to be a valid folder location.", call. = TRUE)}
  
  # Assign object to internal variable
  dispensingData = x
  
  message("> Reading dispensing data for ", dispensingData$dispensingID, ".", appendLF = FALSE)
  
  
  # Plotting plate layout for a given dispensing set
  # Remove untreated and the additional DMSO from the single drug treatments
  dispPlot <- subset(dispensingData[["output"]], Drug != "Untreated")
  dispPlot <- subset(dispPlot, !Drug %in% dispensingData[["dataList"]][["ctrlList"]]$NAME | Drug %in% dispensingData[["dataList"]][["ctrlList"]]$NAME & Drug.Concentration %in% dispensingData[["dataList"]][["ctrlList"]]$DOSE)

  
  # Select only essential columns: drug name, dose and destination well from the dispensing data
  # Convert drug names to capitalized four letter symbols for controls and three letter symbols for drug names
  dispPlot <- dispPlot[grepl("Combination.ID|Drug|Drug.Concentration|Destination.Well|Plate.Number", names(dispPlot), ignore.case = TRUE)]
  dispPlot$Drug <- gsub("[^[:alnum:]]", "", dispPlot$Drug)
  dispPlot$Drug <- with(dispPlot, ifelse(Drug %in% dispensingData[["dataList"]][["ctrlList"]]$NAME, Drug, toupper(Drug)))
  dispPlot$Drug <- with(dispPlot, ifelse(Drug %in% dispensingData[["dataList"]][["ctrlList"]]$NAME, strtrim(Drug, 4), strtrim(Drug, 3)))
  # Round numbers to three digits
  dispPlot$Drug.Concentration <- with(dispPlot, as.numeric(format(round(Drug.Concentration, 3), nsmall = 3, scientific = FALSE, drop0trailing = TRUE)))
  # Create labels for each well
  dispPlot <- within(dispPlot, Drug <- ifelse(Drug %in% dispensingData[["dataList"]][["ctrlList"]]$NAME, Drug, paste(Drug, Drug.Concentration, sep = " ")))
  # Remove duplicate labels with same drug dispensing, such as single drug treatments
  dispPlot <- dispPlot[!duplicated(dispPlot[c("Drug", "Drug.Concentration", "Destination.Well", "Plate.Number")]), ]
  
  # Combine labels, in which two drugs are dispensed into the same well
  dispPlot <- merge(dispPlot[,c("Combination.ID", "Drug.Concentration", "Destination.Well", "Plate.Number")],
                                   aggregate(Drug~Combination.ID, data = dispPlot, FUN = function(x) paste(x, collapse="\n")),
                                   by = "Combination.ID", all = TRUE, sort = FALSE)
  
  # Remove duplicate labels from combining drug treatments
  dispPlot <- dispPlot[!duplicated(dispPlot[c("Combination.ID", "Drug", "Destination.Well", "Plate.Number")]), ]
  
  # Convert the plate number to a numeric value
  dispPlot$Plate.Number <- as.numeric(regmatches(dispPlot$Plate.Number, gregexpr("\\d+", dispPlot$Plate.Number)))
  
  # Rename columns and select only columns needed for plotting
  names(dispPlot) <- sub("Destination.Well", "Well", names(dispPlot))
  names(dispPlot) <- sub("Plate.Number", "Plate", names(dispPlot))
  
  dispPlot <- dispPlot[,c("Drug", "Well", "Plate")]
  
  
  # Add remaining wells to obtain a complete plate layout with a total set of wells
  # Create a list of wells for a given plate format
  baseplate <- baseplate(dispensingData[["plateFormat"]], zeroBase = FALSE, leadingZero = FALSE)
  dispPlot <- merge(dispPlot, setNames(data.frame(rep(baseplate[["wells"]]$id., times = max(dispPlot$Plate)), rep(1:max(dispPlot$Plate), each = dispensingData$plateFormat), stringsAsFactors = FALSE), c("Well", "Plate")), 
                    by = c("Well", "Plate"), all = TRUE,  sort = FALSE)
  # Sort by well after converting it to a factor
  dispPlot$Well <- factor(dispPlot$Well, baseplate[["wells"]]$id.)
  dispPlot <- dispPlot[with(dispPlot, order(Plate, Well)),]

  
  # Split wells into row letters and column numbers
  dispPlot$Row <- as.character(regmatches(dispPlot$Well, gregexpr("[[:alpha:]]+", dispPlot$Well)))
  dispPlot$Column <- as.numeric(regmatches(dispPlot$Well, gregexpr("[[:digit:]]+", dispPlot$Well)))
  
  
  # Plotting each plate to file
  
  for(plate in unique(dispPlot$Plate)){
    
    dfs.plot <- subset(dispPlot, Plate == plate)
    
    ggplot(dfs.plot, aes(x = factor(Column), y = factor(Row, rev(unique(Row))))) + 
      geom_raster(aes(fill = Drug), na.rm=TRUE) +
      #scale_fill_discrete(na.value = NA) +
      scale_fill_hue("Drugs", l=80, h = c(270, 360), na.value = NA) +
      geom_text(aes(label = paste(ifelse(is.na(Drug),"",Drug), sep = "\n")), color ="black", hjust = 0.5, size = 2.7) +
      geom_hline(yintercept=seq(1.5, length(unique(dfs.plot$Row))-0.5, 1), lwd=0.2, colour="black") +
      geom_vline(xintercept=seq(1.5, length(unique(dfs.plot$Column))-0.5, 1), lwd=0.2, colour="black") +
      scale_x_discrete(expand = c(0,0), position = "top") +
      scale_y_discrete(expand = c(0,0)) +
      coord_fixed(ratio = 1) +
      xlab(NULL) +
      ylab(NULL) +
      labs(caption = paste("Dispensing Plate: ", dispensingData$dispensingID, ":", plate, sep = "")) +
      theme(plot.title = element_text(hjust = 0.5),
            plot.caption = element_text(size = 7),
            text = element_text(size = 12),
            axis.line = element_line(color='black'),
            panel.grid.major = element_line(color='lightgrey', linetype = "dotted"),
            panel.grid.minor = element_line(color='lightgrey', linetype = "dotted"),
            panel.background = element_rect(fill=alpha('white', 1)),
            panel.border = element_rect(colour = "black", fill=NA, size=1),
            legend.position = "none"
      )
    
    
    message('\r', "> Plotting plate no. ", plate, " / ", max(dispPlot$Plate), ".", strrep(" ", 100), appendLF = FALSE)
    ggsave(filename = file.path(.saveto, paste(dispensingData$dispensingID, "-Plate", plate, ".png", sep = "")), device = "png", width = 840, height = 840/(3/2), units = "mm", dpi = 300, limitsize = FALSE)
    
    if(plate == max(dispPlot$Plate)){
      cat('\r', "Finished plotting a total of ", plate, " plates for ", dispensingData$dispensingID, ".", strrep(" ", 100), '\n', sep = "")
      rm(plate, dfs.plot)
    }
    
  }
}