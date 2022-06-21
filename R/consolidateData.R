#' Consolidating raw measurements with dispensing data 
#'
#' @description
#' An essential component of the modular library \pkg{metascreen}, and imperative for the post-experimental processing of data from a drug sensitivity screen.
#' \emph{\code{consolidateData}} is a function that consolidates raw measurements from an experimental assay with the dispensing data. It is essential for building a reference data set for downstream analysis.
#' 
#' @param dispensingData an object of class 'dispensingData'.
#' @param rawMeasurements an object of class 'rawMeasurement'.
#' @param .barcodeReference \code{data frame}; a list of samples for each set of plates.
#' 
#' @details The function \code{consolidateData} is used to build a reference data set by mapping raw measurements (\emph{rawMeasurements}) to the dispensing data (\emph{dispensingData}). Only raw measurements from samples listed in the 
#' \emph{.barcodeReference} will be consolidated. The raw data from all other sets will be ignored by default.
#' 
#' The \emph{.barcodeReference} is a list of samples used for a given set of experiments in a drug sensitivity screen. Each sample is listed in reference to the \emph{plate id} (unique to a drug screen), the label of each \emph{set} and the \emph{number of plates} for a given set.
#' 
#' @return Returns an object of class S3:consolidatedData.
#'
#' @examples
#' \donttest{\dontrun{
#' # Consolidate raw measurements with dispensing data
#' consolidateData(dispensingData, rawMeasurements, .barcodeReference)
#' }}
#'
#' @keywords drug screen consolidate raw measurement dispensing
#' 
#' 
#' @export

consolidateData <- function(dispensingData, rawMeasurements, .barcodeReference){
 
  
  # Check, if the dispensing data has been provided as an object of class S3:dispensingData
  if(missing(dispensingData)){stop("Dispensing data missing! Please provide a dispensing data set.", call. = TRUE)}
  if(class(dispensingData) != "dispensingData"){
    stop("Dispensing data not of class 'dispensingData'!", call. = TRUE)
    } else {
      # Assign object to internal variable
      .dispensingData <- dispensingData[["output"]]
    }
  
  # Check, if the dispensing data has been provided as an object of class S3:rawMeasurements
  if(missing(rawMeasurements)){stop("Raw measurements missing! Please provide a list with raw measueremnts for each plate.", call. = TRUE)}
  if(!"rawMeasurements" %in% class(rawMeasurements)){stop("Raw measurement data not of class 'rawMeasurements'!", call. = TRUE)}
   
  
  message(" > Processing raw measurements...", appendLF = FALSE)

  # Replicate the dispensing data based on the number of sets
  .dispensingData <- cbind(.dispensingData, Dispensing.Set = rep(.barcodeReference[["Set"]], each = nrow(.dispensingData)), stringsAsFactors = FALSE)

  # Build the destination plate barcode from the destination plate id, the dispensing set and the plate number
  .dispensingData <- transform(.dispensingData, Destination.Plate.Barcode = paste(Destination.Plate.Barcode, Dispensing.Set, Plate.Number, sep = ""))
  
  

  # Extend data set containing raw measurements with adding additional columns in order to be able to consolidate data sets
  
  dfs = list()
  
  for (i in names(rfs)) {
    
    plateBarcode <- i
    destinationPlateID <- unlist(regmatches(plateBarcode, gregexpr("\\d+", plateBarcode)))[1]
    plateSet <- gsub("\\d+", "", plateBarcode)
    plateNumber <- unlist(regmatches(plateBarcode, gregexpr("\\d+", plateBarcode)))[2]

    # Check, if plate belongs to a set found in the barcode reference list
    # Note: Skip sets that should not be analyzed
    if(!(plateSet %in% .barcodeReference[["Set"]])){next}
    
    # Fetch sample id
    sampleName <- .barcodeReference$Sample[match(paste(destinationPlateID, plateSet, sep = ""), paste(.barcodeReference$PlateID, .barcodeReference$Set, sep = ""))]
    
    
    # Create a data frame based on the raw data set
    dfs[[plateBarcode]] <- rfs[[plateBarcode]]
    
    # Adding additional columns to the data set 
    dfs[[plateBarcode]]$Sample <- sampleName
    dfs[[plateBarcode]]$Destination.Plate.Barcode <- plateBarcode
    dfs[[plateBarcode]]$Plate.Number <- plateNumber
    
    
    # Check, if well numbers have a leading zero in the raw data readings
    # If so, convert well numbers to standard format: A01 to A1
    levels(dfs[[plateBarcode]]$Well) <- gsub("(?<![0-9])0+", "", levels(dfs[[plateBarcode]]$Well), perl = TRUE)

    # Reorder rows, and columns
    dfs[[plateBarcode]] <- dfs[[plateBarcode]][with(dfs[[plateBarcode]], order(Plate.Number, Well)),]
    dfs[[plateBarcode]] <- dfs[[plateBarcode]][,c("Sample", "Destination.Plate.Barcode", "Plate.Number", "Well", grep("CPS", names(dfs[[plateBarcode]]), ignore.case = TRUE, value = TRUE))]
    

    rm(plateBarcode, destinationPlateID, plateSet, plateNumber, sampleName, i)
    
  }
  
  
  message('\r', " > Consolidating data sets...", strrep(" ", 100), appendLF = FALSE)
  
  # Combine all individual raw measurement data sets into one single data set
  # rename, if necessary, the column with the raw measurements
  finalMeasurementData <- do.call(rbind, setNames(dfs, NULL))
  names(finalMeasurementData)[names(finalMeasurementData) == grep("CPS", names(finalMeasurementData), value = TRUE)] <- "CPS"
  
  # Consolidate the dispensing data with the raw measurements
  clDataSet <- merge(.dispensingData, finalMeasurementData, by.x = c("Destination.Plate.Barcode", "Plate.Number", "Destination.Well"), by.y = c("Destination.Plate.Barcode", "Plate.Number", "Well"))
  
  # Set factors for individual columns to allow proper sorting
  clDataSet$Dispensing.Set <- factor(clDataSet$Dispensing.Set, levels = with(clDataSet, unique(Dispensing.Set[order(nchar(Dispensing.Set), Dispensing.Set)])))
  clDataSet$Destination.Well <- factor(clDataSet$Destination.Well, levels = with(finalMeasurementData, unique(Well[order(nchar(gsub("[[:digit:]]", "", Well)), gsub("[[:digit:]]", "", Well), as.numeric(gsub("[[:alpha:]]", "", Well)))])))
  
  # Reorder rows, and columns
  clDataSet <- clDataSet[with(clDataSet, order(Dispensing.Set, Plate.Number, Destination.Well)),]
  clDataSet <- clDataSet[,c("Sample", "Destination.Plate.Barcode", "Plate.Number", "Dispensing.Set", "Combination.ID", "Drug", "CAS.number", "Drug.Concentration", "Unit", "Transfer.Volume", "Source.Plate.Barcode",
                                                "Source.Well", "Destination.Well", "CPS" )]

  
  
  cat('\r', "Finished processing data. Data sets succesfully consolidated.")
  
  
  # Create an object of S3 class with a list from the attributes of all three vectors
  object <- list(consolidated = clDataSet, measurements = finalMeasurementData, dispensingData = dispensingData)
  class(object) <- "consolidatedData"

  return(object)
  
}
