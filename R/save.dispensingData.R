#' Saving a set of dispensing files
#'
#' @description
#' The function \emph{\code{save}} creates a set of dispensing files, which can be used to provide instructions to dispensing robots on how to dispense a set of drugs, controls and samples.
#' Those files serve as a blueprint for printing a set of \emph{destination plates} from a set of \emph{source plates} based on the dispensing layout provided.
#' 
#' @param x an object of class 'dispensingData'.
#' @param .saveto \code{character}; a path to a folder location where the object is saved to.
#' @param .sets \code{numeric}; a value determining the number of sets to create, or logical; if TRUE, the number of sets will depend on the labels provided.
#' @param .labels \code{vector}; a set of labels to be used for each individual set.
#' @param .split \code{logical}; if TRUE, the function will split each set into individual files. The default is FALSE, in which a single file is created.
#' @param .format \code{character}; a predefined identifier, which will decide  the format and layout of the dispensing file (see Details, for more information).
#' @method save dispensingData
#' @S3method save dispensingData
#' 
#' @details This function extends the generic function \code{\link{save}} for objects of class 'dispensingData'. It exports a set of dispensing files for a given number of sets. The possibility to generate multiple identical sets from a single dispensing layout, is required, whenever a drug screen is being performed on multiple samples, e.g. cell lines.
#' In certain instances, an experiment requires to be repeated, either to obtain experimental replicates, or due to a failed experiment. In either case, a single dispensing layout can be replicated based to an unlimited number of sets as needed. 
#' 
#' The labels can be assigned to each dispensing set either through a predefined selection of labels by using either a set of alphabetic or numeric values. 
#' This can be accomplished by setting \code{.format} = "alphabetic", or \code{.format} = "numeric", respectively. Otherwise, the labels have to be explicitly stated for each individual set.
#' 
#' Please note, that if the number of sets is not provided, or explicitly set as \code{FALSE}, any provided labels will be ignored. If the argument is set as TRUE, instead of a numeric value, the number of sets will depend on the number of labels provided. 
#' In that case, for each label a set will be created. If the number of sets is large than the number of labels provided, the function will result in an error, while if the number of sets is smaller than the number of labels provided, only the first \emph{n} labels will be used.
#'  
#' Also in cases, where only one set is being requested, the \code{.split} argument will be ignored, since the dispensing file cannot be split with a single data set.
#' 
#' 
#' Furthermore, depending on how the dispensing files are read and used for dispensing, each set can be saved as an individual file, or as a single file containing all the individual replicates.
#' 
#' The format of the dispensing file, is device specific and depends on the machine that is being used to carry out the dispensing. As of now, the following machines are supported:
#' 
#' \tabular{rlrl}{
#' \verb{  } \tab Type:             \tab \verb{ } \tab Echo Acoustic Liquid Handlers \cr
#' \verb{  } \tab Manufacturer:     \tab \verb{ } \tab Labcyte Inc., Beckman Coulter Life Sciences \cr
#' \verb{  } \tab Instrument Name:  \tab \verb{ } \tab E5XX-1366 \cr
#' \verb{  } \tab Instrument Model: \tab \verb{ } \tab Echo 550 \cr
#' }
#' 
#' 
#' To save the dispensing file in the format supported by the \emph{Echo Acoustic Liquid Handler}, set the argument using the instrument name or model \code{.format} = "E5XX-1366", or \code{.format} = "Echo 550".
#' If another model of the \emph{Echo Acoustic Liquid Handler} is used that supports the same dispensing file format, simply use the instrument group \code{.format} = "Echo".
#' 
#' Alternatively, the complete dispensing data can be saved, which allows machines not listed or currently supported in this package to pick individual instructions for dispensing. 
#' This is accomplished by setting the argument \code{.format} = "full" or \code{.format} = "complete".
#' 
#' The dispensing files will be save as a comma-separated values (csv) file, in compliance with the standard outlined in RFC 4180 (2005). For more information, see \url{https://www.ietf.org/rfc/rfc4180.txt}
#' 
#'  
#' @seealso \code{\link{generateDispensingData}}
#'
#' @examples
#' \donttest{\dontrun{
#' # Save a set of twenty-six dispensing files labeled from 'A to Z' in the format as required by an 'Echo Acoustic Liquid Handler'
#' save(dispensingData, .saveto = "../myDispensing/files", .sets = 26, .labels = "alphabetic", .split = TRUE, .format = "E5XX-1366")
#' 
#' Save a single dispensing files with three sets labeled individually containing the full dispensing information
#' save(dispensingData, .saveto = "../myDispensing/files", .sets = 3, .labels = c("REP1", "REP2", "CTRL"), .split = FALSE, .format = "full")
#' }}
#'
#' @keywords drug screen drug combination dispensing save
#' 
#' 
#' @importFrom utils write.csv
#' 
#' 
#' @export

save.dispensingData <- function(x, .saveto, .sets, .labels, .split = FALSE, .format){
  
  # Check, if the dispensing data has been provided as an object of class S3:dispensingData
  if(missing(x)){stop("Dispensing data missing! Please provide a dispensing data set.", call. = TRUE)}
  if(class(x) != "dispensingData"){stop("Dispensing data not of class 'dispensingData'!", call. = TRUE)}
  
  # Check, if a folder location has been provided
  if(missing(.saveto)){ message("Folder location not specified. Saving to default working directory:\n", getwd()) }
  if(missing(.saveto)){ .saveto <- getwd() }
  # Create folder, if it does not exist
  if(!file.exists(file.path(.saveto))){ dir.create(file.path(.saveto), showWarnings = FALSE, recursive = TRUE) }
  if(!utils::file_test("-d", .saveto)){stop("in '.saveto'. Argument needs to be a valid folder location.", call. = TRUE)}
  
  
  # Check, if sets are requested and how many sets to generate
  if(missing(.sets)){ .sets = FALSE } 
  if(all(.sets == TRUE, missing(.labels))){
    stop("in '.labels'. Dispensing sets requested, but no labels provided. Please provide labels from which to generate a number of sets.", call. = TRUE)
  }
  
  if(is.numeric(.sets)){
    
    if(length(.sets) > 1){ 
      message("Too many arguments have been provided for the dispensing sets. Selecting only the first element.")
      .sets <- .sets[1]
    } else if(.sets < 0){ stop("in '.sets'. Number of sets cannot be negative! Please provide a valid number of sets.", call. = TRUE) } 
    else if(.sets == 0){ .sets = FALSE } 
    else if(grepl(".", .sets, fixed = TRUE)){ stop("in '.sets'. Number of sets cannot have a decimal point and must be a whole number! Please provide a valid number of sets.", call. = TRUE) }
    
  } else if(!any(is.numeric(.sets), is.logical(.sets))) {
    stop("in '.sets'. Argument not recognized!\nPlease provide valid number of sets, or provide labels from which to generate a number of sets.", call. = TRUE)
  }
  

  
  # Check, if labels have been provided
  # and if they are in the proper combination with the dispensing set flag
  if(all(missing(.labels), is.numeric(.sets))){
    stop("Labels missing! Dispensing sets requested, but no labels provided. Please provide labels, which to be assigned to each disepnsing set.")
  # Check, if sets are requested, but labels are missing
  } else if(all(missing(.labels), .sets)){
    stop("Labels missing! Please provide labels from which to generate a number of sets." )
  }
  
  # This sections checks, whether custom labels should be used
  # Check, if labels are requested as alphabetic or numeric
  if(all(grepl("alpha", .labels), is.numeric(.sets))){
  .labels <- c(LETTERS, sapply(LETTERS, function(x) paste0(x, LETTERS)))[1:.sets]
  } else if(all(grepl("num|digit", .labels),  is.numeric(.sets))){
  .labels <- 1:.sets
  } else if(all(grepl("alpha|num|digit", .labels), !is.numeric(.sets), .sets != FALSE)) {
    stop("Number of sets missing! The selected argument for labels is only valid with a numeric dispensing set value. \nPlease provide a valid number of sets, or provide labels from which to generate a number of sets." )
  }
  
  # Check, if the number of labels corresponds with the number of sets
  if(all(is.vector(.labels), is.numeric(.sets))){
    if(length(.labels) < .sets){
      stop(.sets-length(.labels), " label(s) missing. A total of ", .sets, " sets have been requested, but only ", length(.labels), " labels have been provided.")
    } else if(length(.labels) > .sets){ 
      message(length(.labels)-.sets, " more label(s) provided than neccessary. Using the first ", length(.labels)-(length(.labels)-.sets), " labels.")
      .labels <- .labels[1:.sets]
    }
  # Check, if number of sets has been provided through a logical argument
  } else if(all(is.vector(.labels), .sets)){
    .sets <- length(.labels)
  } else if(all(is.logical(.labels), any(.sets, is.numeric(.sets)))){
    if(.labels){
      stop(.sets, " sets requested, but no labels provided.\nPlease provide labels for each set.")
    }
  }
  
  # Ignore splitting flag, if sets are set to FALSE
  if(all(any(.sets == FALSE, .sets == 1), .split == TRUE)) { 
    message("Argument '.split' ignored. Cannot split a single data set.")
    .split = FALSE 
  }
  
  # Ignore splitting flag, if sets are set to FALSE
  if(all(.sets == FALSE, any(is.vector(.labels), .labels == TRUE))) { 
    message("Argument '.labels' ignored. Dispensing sets have not been requested.")
    .split = FALSE 
  }

  
  
  # Check, if the format have been specified
  if(missing(.format)){
    message("Dispensing format not sepecified. Using default: 'Echo 550:E5XX-1366'")
    .format <- "Echo"
  } else if(any(.format == "E5XX-1366", grepl("Echo", .format, ignore.case = TRUE))){
    .format <- "Echo"
  } else if(grepl("full|complete|auto", .format, ignore.case = TRUE)){
    .format <- "full"
  } else {
    stop("in '.format'. Dispensing format not supported! Please select a valid dispensing format.\nFor more information and help, see 'R Documentation'.", call. = TRUE)
  }
  
  
  # Assign object to internal variable
  dispensingData = x
  
  
  message('\n', "> Reading dispensing data...", appendLF = FALSE)
  
  # Extract the dispensing data and remove untreated controls
  .dispensingData <- dispensingData[["output"]]
  .dispensingData <- subset(.dispensingData, !Drug %in% dispensingData[["origData"]][[".addUntreated"]]$name)
  .dispensingData <- transform(.dispensingData, Plate.Number = as.numeric(regmatches(Plate.Number, gregexpr("[[:digit:]]+", Plate.Number))))
  rownames(.dispensingData) <- NULL
  
  # Replicate the dispensing data based on the number of sets
  # .dispensingData <- cbind(.dispensingData, Dispensing.Set = rep(1:.sets, each = nrow(.dispensingData)))
  # Note: Instead of replicating the data set, the dispensing data is saved multiple times over based on the number of sets requested

  # Assign the provided labels to each set
  # if(.sets == TRUE){
  #   .dispensingData <- transform(.dispensingData, Dispensing.Set = .labels[.dispensingData$Dispensing.Set])
  # }
  
  
  # Build the destination plate barcode from the destination plate id, the dispensing set and the plate number
  # .dispensingData <- transform(.dispensingData, Destination.Plate.Barcode = paste(Destination.Plate.Barcode, Dispensing.Set, Plate.Number, sep = ""))
  

  
  # Select columns from the complete dispensing data set based on the dispensing format selected
  .select <- function(.format){
    switch(.format, 
           Echo = {
             c("Drug", "Source.Plate.Barcode", "Source.Well", "Transfer.Volume",  "Destination.Well", "Destination.Plate.Barcode")
           },
           full = {
             # selecting full data set, no columns will be removed
             names(.dispensingData)
           },
           {
             # default: based on Echo 550:E5XX-1366
             c("Drug", "Source.Plate.Barcode", "Source.Well", "Transfer.Volume",  "Destination.Well", "Destination.Plate.Barcode")
           }
    )
  }
  
  
  # If splitting the data set into individual file is not requested, the dispensing data is exported as a single data file
  if(.split){
    
    for (s in 1:.sets) {
      # Assign the corresponding labels to each set
      .dispensingFile <- transform(.dispensingData, Dispensing.Set = .labels[s])
      # Build the destination plate barcode from the destination plate id, the dispensing set and the plate number
      .dispensingFile <- transform(.dispensingFile, Destination.Plate.Barcode = paste(Destination.Plate.Barcode, Dispensing.Set, Plate.Number, sep = ""))
      # Select only columns based on the selected dispensing format
      .dispensingFile <- .dispensingFile[.select(.format)]
      
      
      message('\r', "> ", "[", s, "/", .sets, "]", " Saving dispensing set: ", .labels[s], appendLF = FALSE)
      
      # Export the corresponding set of the dispensing data to an individual .csv file
      write.csv(.dispensingFile, file = file.path(.saveto, paste("Dispensing file", "_", dispensingData$dispensingID, .labels[s], "_V2", ".csv", sep = "")), row.names = FALSE, quote = FALSE)
    } 
    
    cat('\r', "Finished saving all dispensing files for ", dispensingData$dispensingID, ".", strrep(" ", 100), '\n', sep = "")
    
  } else {
    
    # Replicate dispensing data, if requested and assign a label to each set
    if(.sets != FALSE){

      # Replicate the dispensing data based on the number of sets
      .dispensingData <- cbind(.dispensingData, Dispensing.Set = rep(1:.sets, each = nrow(.dispensingData)))
      # Note: Instead of replicating the data set, the dispensing data could also be saved multiple times over based on the number of sets requested
      
      # Assign the provided labels to each set
      .dispensingData <- transform(.dispensingData, Dispensing.Set = .labels[.dispensingData$Dispensing.Set])
      
      # Build the destination plate barcode from the destination plate id, the dispensing set and the plate number
      .dispensingData <- transform(.dispensingData, Destination.Plate.Barcode = paste(Destination.Plate.Barcode, Dispensing.Set, Plate.Number, sep = ""))
      
    }
    
    # Select only columns based on the selected dispensing format
    .dispensingFile <- .dispensingData[.select(.format)]
    
    message('\r', "> Saving dispensing file... ", strrep(" ", 100), appendLF = FALSE)
    
    # Export the dispensing data to a single .csv file
    write.csv(.dispensingFile, file = file.path(.saveto, paste("Dispensing file", "_", dispensingData$dispensingID, "_V2", ".csv", sep = "")), row.names = FALSE, quote = FALSE)
    
    cat('\r', "Finished saving all dispensing files for ", dispensingData$dispensingID, ".", strrep(" ", 100), '\n', sep = "")
    
  }
  
  
}