#' Reading raw measurement files
#'
#' @description
#' An essential component of the modular library \pkg{metascreen}, and imperative for the post-experimental processing of data from a drug sensitivity screen.
#' \emph{\code{readFiles}} is a function that reads raw data files from an experimental assay. It supports the most commonly used file formats from a selection of different read-out methods and machines.
#' 
#' @param .readfrom \code{character}; a path to a folder location containing the raw data files to import.
#' @param .fileformat \code{character}; a string or list determining the file format of the raw data files.
#' @param .format \code{character}; a predefined identifier, which will decide the format and layout of the raw data files (see Details, for more information).
#' 
#' @details The function \code{readFiles} is used to read raw measurement files from selected plate readers. It supports the most common text-based file types (.txt, .csv) in various formats.
#' Those formats are machine specific and follow a proprietary layout. This function is capable of detecting the export-format and identifying the used text delimiter, such as comma-, semicolon- or tab-separated. 
#' It will automatically detect the plate format used and find the raw data of interest.
#' 
#' As of now, files from the following machines are supported:
#' 
#' \tabular{rlrl}{
#' \verb{  } \tab Type:             \tab \verb{ } \tab Multimode Plate Reader  \cr
#' \verb{  } \tab Manufacturer:     \tab \verb{ } \tab PerkinElmer, Inc. \cr
#' \verb{  } \tab Instrument Name:  \tab \verb{ } \tab EnVision multimode plate reader \cr
#' \verb{  } \tab Instrument Model: \tab \verb{ } \tab \strong{EnVision} 21xx \cr
#' \verb{  } \tab                   \tab \verb{ } \tab \cr
#' }
#' \tabular{rlrl}{
#' \verb{  } \tab Type:             \tab \verb{ } \tab Multimode Plate Reader  \cr
#' \verb{  } \tab Manufacturer:     \tab \verb{ } \tab PerkinElmer, Inc. \cr
#' \verb{  } \tab Instrument Name:  \tab \verb{ } \tab VICTOR X multimode plate reader \cr
#' \verb{  } \tab Instrument Model: \tab \verb{ } \tab \strong{VICTOR} X \cr
#' }
#'
#' @return Returns a list with a data frame for each input file.
#'
#' @references Robert Hanes et al.
#' @author Robert Hanes
#' @note Version: 2021.03
#'
#' @examples
#' \donttest{\dontrun{
#' # Read a set of output files from an EnVision plate reader
#' readFiles(.readfrom = "path/to/folder/", .fileformat = ".csv", .format="EnVision")
#' }}
#'
#' @keywords drug screen drug read file plate reader
#' 
#' @importFrom utils read.delim
#' 
#' @export

readFiles <- function(.readfrom, .fileformat, .format){
  
  
  # Check, if a folder location has been provided to read file from
  if(missing(.readfrom)){ message("Folder location not specified. Saving to default working directory:\n", getwd()) }
  if(missing(.readfrom)){ .readfrom <- getwd() }
  # Create folder, if it does not exist
  if(!file.exists(file.path(.readfrom))){ dir.create(file.path(.readfrom), showWarnings = FALSE, recursive = TRUE) }
  if(!utils::file_test("-d", .readfrom)){stop("in '.readfrom'. Argument needs to be a valid folder location.", call. = TRUE)}
  
  
  # Check, if the format have been specified
  if(missing(.format)){
    message("Dispensing format not sepecified. Using default: 'EnVision'")
  } else if(all(.format != "EnVision", !grepl("envision", .format, ignore.case = TRUE)) &
            all(.format != "VICTOR", !grepl("victor", .format, ignore.case = TRUE))){
    stop("in '.format'. Data format not supported! Please select a valid data format. Supported formats: 'EnVision', 'VICTOR'. \nFor more information and help, see 'R Documentation'.", call. = TRUE)
  }
  
  # Check, if the file formats have been specified
  if(missing(.fileformat)){
    message("File format not sepecified. Using default: '.csv', '.txt'")
    .fileformat <- c(".csv", ".txt")
  }
  
  
  
  # Create a list of files to read and import. Remove directories from the list.
  fileList <- setdiff(list.files(paste(.readfrom, sep = ""), pattern = paste("\\", .fileformat, "$", sep = "", collapse = " || "), full.names = TRUE, ignore.case = TRUE, include.dirs = FALSE, recursive = FALSE),
                      list.dirs(paste(.readfrom, sep = ""), full.names = TRUE, recursive = FALSE))
  
  
  # SECTION A: READING ALL DATA FILES ##############################
  # Import all measurement files into a list of data frames
  
  
  rfs = list()
  
  for (i in fileList) {
    
    filename <- gsub("(.*\\/)([^.]+)(\\.[[:alnum:]]+$)", "\\2", i)
    fileextension <- gsub("(.*\\/)([^.]+)(\\.[[:alnum:]]+$)", "\\3", i)
    
    print(paste("Reading file: ", basename(i), sep=""))
    
    
    # Identify text deliminator for each file
    filedelim <- ifelse(grepl(",", readLines(i, n=2)[2]), ",", 
                        ifelse(grepl(";", readLines(i, n=2)[2]), ";",
                               ifelse(grepl("\t", readLines(i, n=2)[2]), "\t")))
    
    print(paste(" Detected text delimiter: ", ifelse(filedelim == ",", "comma-separated [,]",
                                                     ifelse(filedelim == ";", "semicolon-separated [;]",
                                                            ifelse(filedelim == "\t", "tab-separated [ ]"))), sep=""))
    
    
    # Read out the exported format for a given file
    exportformat <- unlist(lapply(strsplit(grep("Export format", readLines(i), value = TRUE),
                                           filedelim), function(x) x[!x %in% ""]))[2]
    
    print(paste(" Detected file export format: ", exportformat, sep=""))
    
    
    
    # SECTION B: EXTRACT INFORMATION FROM DATA FILES ##############################
    # Identify plate type and define the number of rows and columns based on plate type
    
    plateType <- as.numeric(unlist(lapply(strsplit(grep("Number of the wells in the plate", readLines(i), value = TRUE),
                                                   filedelim), function(x) x[!x %in% ""]))[2])
    
    print(paste(" Detected plate type: ", plateType, "-well format", sep=""))
    cat("\n")
    
    
    # Check if the measurement section contains headers
    headers <- length(grep("Results for", readLines(i), value = FALSE)) != 0
    
    
    # Determine the first line of raw measurements
    switch(exportformat, 
           Plate = {
             # set the first line of measurements based on export format 'Plate'
             firstline <- ifelse(headers == FALSE,
                                 grep("Background information", readLines(i), value = FALSE)+3,
                                 grep("Background information", readLines(i), value = FALSE)+3)
           },
           Plate4 = {
             # set the first line of measurements based on export format 'Plate4'
             firstline <- ifelse(headers == FALSE,
                                 grep("Background information", readLines(i), value = FALSE)[2]+3,
                                 grep("Results for", readLines(i), value = FALSE))
           },
           {
             # default:
             message("Export format '", exportformat, "' not supported.", appendLF = FALSE)
           }
    )
    
    
    # Estimate number of columns and rows for a given plate type
      nofRow <- as.integer(sqrt(plateType))
      while (plateType %% nofRow != 0) {
        nofRow <- nofRow-1
      }
      nofCol <- plateType / nofRow
      noWells <- nofCol*nofRow

    
    
      # Read files into a list of individual data frames for each cell line and plate (barcode)
      if(headers == FALSE){
        readings <- utils::read.delim(file.path(i), check.names=FALSE, header=FALSE,
                               stringsAsFactors=FALSE, strip.white = TRUE, 
                               sep = filedelim, dec=".", skip=firstline, nrows=nofRow)[1:nofRow, 1:nofCol]
        colnames(readings) <- seq(1, nofCol, 1)
        rownames(readings) <- c(LETTERS, sapply(LETTERS, function(x) paste0(x, LETTERS)))[1:nofRow]
      }else{
        readings <- utils::read.delim(file.path(i), check.names=FALSE, header=TRUE, row.names = 1,
                               stringsAsFactors=FALSE, strip.white = TRUE, 
                               sep = filedelim, dec=".", skip=firstline, nrows=nofRow)[1:nofRow, 1:nofCol]
      }
    
    
    
      # Converting the format of the data from a matrix layout to a list
      rfs[[filename]] <- data.frame(Well = paste(rownames(readings)[row(readings)], colnames(readings)[col(readings)], sep = ""), CPS = unlist(readings))
      
      
      rm(i, filename, fileextension, readings, filedelim, firstline, headers, exportformat, nofCol, nofRow, noWells)
    
  }
  
  
  class(rfs) <- c("readFiles", "rawMeasurements")
  
  return(rfs)

}