#' Plotting raw measurement files
#'
#' @description
#' The function \emph{\code{plot}} generates individual plots for each plate of a raw measurement file. 
#' 
#' @param x an object of class 'rawMeasurements'.
#' @param .saveto string; path to a folder location where the object is saved to.
#' @param ... further arguments passed to or from other methods.
#' 
#' @details This function extends the generic function \code{\link{plot}} for objects of class 'rawMeasurements'. 
#' It is used to visualize the raw measurement of each \emph{measurement plate}.
#' The color range corresponds to the range of the measured signal across all plates. 
#' The plots will be save to the provided folder location. If the folder does not exist, it will be created, if the folder location is not provided, 
#' the plots will be saved to the current working directory.
#' 
#' @seealso \code{\link{readFiles}}
#'
#' @examples
#' \donttest{\dontrun{
#' plot(rawMeasurements, .saveto = "path/to/folder/")
#' }}
#'
#' @keywords drug screen raw measurement plot
#'
#'
#' @importFrom ggplot2 ggplot ggsave labs scale_fill_gradient
#' @importFrom utils tail
#' 
#' @export

plot.readFiles <- function(x, .saveto, ...){
  
  # Check, if the dispensing data has been provided as an object of class S3:dispensingData
  if(missing(x)){stop("Data missing! Please provide a set of raw measurements.", call. = TRUE)}
  if("rawMeasurements" %in% class(x)){stop("Provided data not of class 'rawMeasurements'!", call. = TRUE)}
  
  # Check, if a folder location has been provided
  if(missing(.saveto)){ .saveto <- getwd() }
  # Create folder, if it does not exist
  if(!file.exists(file.path(.saveto))){ dir.create(file.path(.saveto), showWarnings = FALSE, recursive = TRUE) }
  if(!utils::file_test("-d", .saveto)){stop("in '.saveto'. Argument needs to be a valid folder location.", call. = TRUE)}
  
  # Assign object to internal variable
  rfs = x
  
  
  # Retrieve the minimum and maximum CPS measurement
  maxCPSLimit <- max(sapply(rfs, function(x) max(x$CPS)))
  minCPSLimit <- min(sapply(rfs, function(x) min(x$CPS)))

  
  
  for (plate in names(rfs)){
    
    match(plate, names(rfs))
    
    cat('\r', "[", match(plate, names(rfs)), "/", length(rfs), "]", " > Reading raw measurements for plate ", plate, ".", strrep(" ", 100), sep = "")
    
    data <- rfs[[plate]]
    
    data$Row <- as.character(with(data, gsub("\\d", "", Well)))
    data$Column <- as.numeric(with(data, gsub("\\D", "", Well)))
    data$Row <- factor(data$Row, levels = with(data, unique(Row[order(nchar(Row), Row)])))
    
    
    ggplot(data, aes(x = factor(Column), y = factor(Row, rev(levels(Row))))) + 
      geom_raster(aes(fill = CPS), na.rm=TRUE) +
      scale_fill_gradient(low = "#FFCCCC", high = "#FF0033", limits = c(minCPSLimit, maxCPSLimit)) +
      geom_hline(yintercept=seq(1.5, length(unique(data$Row))-0.5, 1), lwd=0.2, colour="black") +
      geom_vline(xintercept=seq(1.5, length(unique(data$Column))-0.5, 1), lwd=0.2, colour="black") +
      scale_x_discrete(position = "top") +
      coord_fixed(ratio = 1) +
      xlab(NULL) +
      ylab(NULL) +
      labs(caption = paste("Raw Measurement Plate: ", plate, sep = "")) +
      theme(plot.title = element_text(hjust = 0.5),
            plot.caption = element_text(size = 10),
            text = element_text(size = 26),
            axis.line = element_line(color='black'),
            panel.grid.major = element_line(color='lightgrey', linetype = "dotted"),
            panel.grid.minor = element_line(color='lightgrey', linetype = "dotted"),
            panel.background = element_rect(fill=alpha('white', 1)),
            panel.border = element_rect(colour = "black", fill=NA, size=1),
            legend.position = "none")
    

    cat('\r', "[", match(plate, names(rfs)), "/", length(rfs), "]", " > Plotting raw measurements for plate ", plate, ".", strrep(" ", 100), sep = "")
    ggsave(filename = file.path(.saveto, paste("Plate", plate, ".png", sep = "")), device = "png", width = 840, height = 840/(3/2), units = "mm", dpi = 300, limitsize = FALSE)
    
    if(plate == utils::tail(names(rfs),n = 1L)){
      cat('\r', "Finished plotting raw measurements for plate ", plate, ".", strrep(" ", 100), '\n', sep = "")
      rm(data, rfs, x)
    }
    
  }
  
}
