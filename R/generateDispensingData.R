#' Generate a dispensing data set
#'
#' @description
#' An essential component of the modular library \pkg{metascreen}, imperative for the initial set-up of a drug sensitivity screen.
#' \emph{\code{generateDispensingData}} is a function that generates a dispensing data set from a list of drug treatments. 
#' It allows to probe a hypothetical dispensing prior to generating a dispensing data set or any dispensing files, providing a summary with essential information about the potential dispensing.
#' 
#' @param listofCombinations \code{data frame}; a list of drug treatments in long-format. The data frame needs to include a column with the drug name, dose and unit for each drug. For help generating a list of combinations, see \emph{combineDrugs}.
#' @param listofDrugs \code{data frame}; a list of drugs with a unique drug id, name and CAS number.
#' @param listofDoses \code{data frame}; a list of drugs and doses in long-format with essential columns that include the drug names, all doses and the units. For help generating a list of doses, see \emph{generateListofDoses}.
#' @param listofVolumes \code{data frame}; a list of volumes for each drug and dose in wide- or long-format with essential columns that include the drug names, the dispensing volume for each dose and the units.
#' @param listofCtrls \code{data frame}; a list of controls to include in the screen. The data frame needs to include a column with the name, dose, unit, dispensing volume and the CAS number for each of the controls.
#' @param listofStockConcentrations \code{data frame}; a list of stock concentrations for each of the drugs with essential columns that include a unique drug id, drug name, the stock concentration of the drugs and the unit.
#' @param sourcePlate \code{data frame}; the source plate with the drugs and doses to dispense from. The data frame needs to include a column with a unique plate id, the well, the drug name, drug concentration and unit. For help importing a source plate from a .PlateMap file, see \emph{importPlateMap}.
#' @param listofExWells \code{vector}; a list of wells to excluded from dispensing to. For help generating a list of excluded wells, see \emph{excludeWells}.
#' @param .ctrlReplicates \code{numeric}; a value providing the number of replicates for each of the controls.
#' @param .addUntreated \code{list}; a list containing the name and the number of replicates for the untreated control. The list needs to contain the named elements \emph{name} and \emph{replicates}, respectively. If no untreated controls are to be used, the argument can be left out or explicitly set to FALSE.
#' @param .finalWellVolume \code{numeric}; a value providing the final volume across all wells found on the destination (treatment) plate.
#' @param .plateFormat \code{integer}; a predefined numeric value designating the plate format based on either 6, 12, 24, 48, 96, 384 or 1536 wells.
#' @param .destinationPlateID \code{character}; a unique id for a set of destination (treatment) plates.
#' @param .randomizeDispensing \code{logical}; if TRUE, then drug treatments will be randomized across all plates. If FALSE, drug treatments will be dispensed in sequential order. The default is TRUE, in which the drug treatments are randomized.
#' @param .probeDispensing \code{logical}; if TRUE, then the dispensing will be probed only and no dispensing data will be generated. The default is FALSE, in which probing is being omitted.
#'
#' @details The function \code{generateDispensingData} is used to accomplish at least two fundamental objectives: (a) to generate dispensing files, which are used to dispense a list of drugs from a \strong{source plate} onto a \strong{destination plate}, and/or (b) to map the dispensing to the readout of a drug screen.
#' 
#' The dispensing data basically contains instructions from where to pick individual drugs and where to dispense them to. In most drug sensitivity screens, drugs are dispensed from a \emph{source plate} onto a \emph{destination plate}. The \emph{destination plate} is the plate, in which an assay, and hence the actual experiment, is conducted.
#' At the final stage of the experiment, the assay will be terminated and the drug treatments evaluated. To map the readout of the assay to the actual drug treatments, the dispensing data is needed.
#' In order to successfully generate a dispensing set, the provision of a number of lists is required. Those lists are essential for building a complete dispensing data set.
#' 
#' A list of treatments, such as the \emph{listofCombinations} for a drug combination screen, contains a number of drug combinations, which have been previously generated from a list of drugs and doses. The list of combinations can be generated using another function in this package, \code{\link{combineDrugs}}.
#' 
#' Furthermore, a number of additional list are needed. The \emph{listofDrugs} and the \emph{listofCtrls} serve as a reference list, both containing unique identifiers (id, and cas numbers) for each of the drugs used in the drug screen.
#' The \emph{listofDoses} provides the doses at which each drug is being dispensed. This list is a precursor of the \emph{listofCombinations}. In most cases, data is only available in a human readable format, which often is provided in a long-format.
#' In order to increase efficiency, the data is converted into a machine friendly, long-format. The function \emph{transformListofDoses} found in this package was build with this purpose in mind and is used to convert a list of \emph{listofDoses} from wide- to long-format.
#' 
#' The same is being done with the \emph{listofVolumes}. A list that can be provided in a wide-format, in which the function will convert it to a long-format.
#' 
#' The \emph{listofStockConcentrations} contains another reference list with the stock concentrations of all drugs provided. Please note, that for all reference lists, a set of unique identifiers is required. All drugs must have an unique id and should have a unique CAS number.
#' 
#' The \emph{sourcePlate} is another essential component, in order to be able to generate a dispensing set. It serves as a map to the \emph{source plate} from which individual drugs are picked and dispensed onto the \emph{destination plate}. The source plate contains therefore a list of wells with drugs at different concentrations. 
#' Please note: as of current, the doses of drugs on the destination plate can only be mapped to the concentrations of drugs on the source plate though the volumes provided. The well on the source plate, containing the corresponding drug concentration for dispensing, is identified through: \code{((final volume / transfer volume) * drug concentration)}.
#' In order generate a source plate from a plate map, such as from the proprietary .PlateMap file format, the function \code{\link{importPlateMap}} can be used.
#' 
#' All wells in the \emph{listofExWells} will be excluded from dispensing into. This list can be conveniently generated using the function \code{excludeWells()} from this package.
#' 
#' The argument \emph{.ctrlReplicates} allows to set the number of replicates for all the controls that have been provided in the \emph{listofCtrls}. If an untreated control is needed, this can be done through the argument \emph{.addUntreated}, which has to be a list containing the name and the number of replictes for the untreated control, i.e. \code{.addUntreated = list(name = "Untreated", replicates = 8)}
#' 
#' The argument \emph{.finalWellVolume} provides the final and total volume in all wells on the \emph{destination plate} that has been dispensed into. This argument is important in mapping the stock concentration on the \emph{source plate} to the final drug concetration to the \emph{destination plate}.
#' 
#' The argument \emph{.plateFormat} defines the format of the \emph{destination plate}, in which the assay is performed. This value is representative of the number of wells found on the microplate.
#' 
#' The argument \emph{.destinationPlateID} allows to provide a unique dispensing, or experimental, identifier from which a unique \emph{destination plate} barcode is being generated for each plate and set.
#' 
#' The argument \emph{.randomizeDispensing} offers the possibility to randomize drug treatments across all plates of a dispensing set. Randomization offers the benefit of reducing spatial variance by randomizing the noise across a larger sample size. If drug treatments are dispensed in sequential order, the technical noise as a result of the experimental variation will impact the entire dose range of a drug treatment. 
#' This can cause the result to be skewed and potentially result in a false positive or negative interpretation of the response. Hence, due to the critical importance of randomization for quality assurance, the default is set to TRUE. If the argument is set to FASLE, the drug treatments will be dispensed in sequential order as found in the \emph{listofCombinations}.
#' 
#' 
#' An important feature of this function is the possibility to probe a dispensing based on the data provided prior to generating a dispensing file. This can be crucial, if there are limitation in place in form of equipment, material, time or any other resources, in which only a limited number of drug treatments, plates or sets can be experimentally carried out.
#' This argument will generate a short summary and allow the user to estimate the number of drugs, controls, and plates among other things required for a given experimental set-up and make necessary changes, if the limitations require it. In order to probe a dispensing, only a limited number of essential data is required, which are 
#' \emph{listofCombinations}, \emph{listofCtrls}, \emph{listofExWells}, \emph{.ctrlReplicates}, \emph{.addUntreated}, and \emph{.plateFormat}
#' 
#' 
#' The function \code{generateDispensingData} generate an object of class S3, which can be used to extract the dispensing data set with the function \code{\link{print}}, to generate a dispensing summary with \code{\link{summary}}, plot the dispensing data for each plate with \code{\link{plot}} or save each dispensing set to file with \code{\link{save}}.
#' 
#' @return Returns an object of class S3:dispensingData.
#'
#' @references Robert Hanes et al.
#' @author Robert Hanes
#' @note Version: 2021.03
#'
#' @examples
#' \donttest{\dontrun{
#' library(metascreen)
#' 
#' # Generate a dispensing data set
#' generateDispensingData(listofCombinations, listofDrugs, listofDoses, listofVolumes, listofCtrls, 
#' listofStockConcentrations, sourcePlate, listofExWells, .ctrlReplicates = 8, 
#' .addUntreated = list(name = "Untreated", replicates = 8), .finalWellVolume = 5, 
#' .plateFormat = 1536, .destinationPlateID = "0521", .randomizeDispensing = TRUE, 
#' .probeDispensing = FALSE)
#' 
#' # Probe a potential dispensing for summary statistics
#' generateDispensingData(listofCombinations, listofCtrls, listofExWells, .ctrlReplicates = 8, 
#' .addUntreated = list(name = "Untreated", replicates = 8), .plateFormat = 1536, 
#' .probeDispensing = TRUE)
#' }}
#'
#' @keywords drug screen drug combination dispensing
#'
#'
#' @importFrom stats setNames
#' 
#' @export

generateDispensingData <- function(listofCombinations, listofDrugs, listofDoses, listofVolumes, listofCtrls, listofStockConcentrations, sourcePlate,
                                   listofExWells, .ctrlReplicates, .addUntreated, .finalWellVolume, .plateFormat, .destinationPlateID, 
                                   .randomizeDispensing = TRUE, .probeDispensing = FALSE){
  
  
  
  # Check, if arguments are provided and in the proper format
  if(missing(listofCombinations)){stop("Data missing! Please provide a list of combinations.", call. = TRUE)}
  if(!all(sapply(c("drug.1", "dose.1", "unit.1", "drug.2", "dose.2", "unit.2"), function(x) any(grepl(x, names(listofCombinations), ignore.case = TRUE))))){stop("in 'listofCombinations'. List of combinations has missing data! Please provide a data frame containing a column for the drug name ['Drug.1, Drug.2'], the dose ['Dose.1, Dose.2'] and the unit ['Unit.1, Unit.2'] for each drug pair.", call. = TRUE)
  } else {
    colnames(listofCombinations)[grepl("drug.1", names(listofCombinations), ignore.case = TRUE)] <- "Drug.1"
    colnames(listofCombinations)[grepl("dose.1", names(listofCombinations), ignore.case = TRUE)] <- "Dose.1"
    colnames(listofCombinations)[grepl("unit.1", names(listofCombinations), ignore.case = TRUE)] <- "Unit.1"
    colnames(listofCombinations)[grepl("drug.2", names(listofCombinations), ignore.case = TRUE)] <- "Drug.2"
    colnames(listofCombinations)[grepl("dose.2", names(listofCombinations), ignore.case = TRUE)] <- "Dose.2"
    colnames(listofCombinations)[grepl("unit.2", names(listofCombinations), ignore.case = TRUE)] <- "Unit.2"
  }
  

  if(missing(listofCtrls)){stop("Data missing! Please provide a list of controls.", call. = TRUE)}
  if(!is.list(listofCtrls)){
    stop("in 'listofCtrls'. Argument needs to be an object of class data.frame! Please provide a data frame containing essential columns with the name of the control ['Name'], the CAS number ['CAS number'], dose ['Dose'], unit ['Unit'] and the volume ['Volume'] at atwhich to dispense the control.", call. = TRUE)
  } else if(!all(sapply(c("name", "cas", "dose", "unit", "vol"), function(x) any(grepl(x, names(listofCtrls), ignore.case = TRUE))))){stop("in 'listofCtrls'. List of controls has missing data! Please provide a data frame containing essential columns with the name of the control ['Name'], the CAS number ['CAS number'], dose ['Dose'], unit ['Unit'] and the volume ['Volume'] at atwhich to dispense the control.", call. = TRUE)
    colnames(listofCtrls)[grepl("name", names(listofCtrls), ignore.case = TRUE)] <- "NAME"
    colnames(listofCtrls)[grepl("number", names(listofCtrls), ignore.case = TRUE) & !grepl("cas", names(listofCtrls), ignore.case = TRUE)] <- "ID"
    colnames(listofCtrls)[grepl("cas", names(listofCtrls), ignore.case = TRUE)] <- "CAS_NUMBER"
    colnames(listofCtrls)[grepl("dose", names(listofCtrls), ignore.case = TRUE)] <- "DOSE"
    colnames(listofCtrls)[grepl("unit", names(listofCtrls), ignore.case = TRUE)] <- "UNIT"
    colnames(listofCtrls)[grepl("vol", names(listofCtrls), ignore.case = TRUE)] <- "VOLUME"
  }
  
  
  if(missing(listofExWells)){ listofExWells = NULL }
  if(missing(.ctrlReplicates)){
    message("Number of control replicates not specified. Using default: 1\nControls won't be replicated.")
    .ctrlReplicates = 1
    }
  
  
  # Handle arguments for the untreated control
  if (any(missing(.addUntreated), isFALSE(.addUntreated))){ 
    .addUntreated = NULL 
  } else if (isTRUE(.addUntreated)) { 
    message("Untreated control not specified. Using default name: 'Untreated' and default number of replicates: 1")
    .addUntreated = list(name = "Untreated", replicates = 1) 
  } else if (class(.addUntreated) != "list") { stop("in '.addUntreated'. Argument needs to be a list! Please provide a list with a name and the number of replicates for the untreated control.", call. = TRUE) 
  } else if (length(.addUntreated) == 0){ stop("in '.addUntreated'. List cannot be empty! Please provide a name and the number of replicates for the untreated control.", call. = TRUE)
  } else if (length(.addUntreated) == 1){ stop("in '.addUntreated'. List is missing an element! Please provide both a name and the number of replicates for the untreated control.", call. = TRUE)
  } else if (!is.character(.addUntreated[[1]])){ stop("in '.addUntreated'. First element of the list needs to be a character! Please provide a valid name for the untreated control.", call. = TRUE)
  } else if (!is.numeric(.addUntreated[[2]])){ stop("in '.addUntreated'. Second element of the list needs to be a number! Please provide a valid number of replicates for the untreated control.", call. = TRUE)
  } else if (.addUntreated[[2]] == 0){
    .addUntreated = NULL
    message("Replicates for untreated control set to zero. No untreated control will be used.")
  } else if (.addUntreated[[2]] < 0){ stop("in '.addUntreated'. Second element of the list cannot be negative! Please provide a valid number of replicates for the untreated control.", call. = TRUE)
  } else if (.addUntreated[[2]] < 1){ stop("in '.addUntreated'. Second element of the list cannot be less than one! Please provide a valid number of replicates for the untreated control.", call. = TRUE)
  } else {
    .addUntreated <- setNames(.addUntreated[1:2], c("name", "replicates"))
    .addUntreated[["name"]] <- as.character(.addUntreated[[1]])
    .addUntreated[["replicates"]] <- as.integer(round(.addUntreated[[2]]))
  }
  
  
  # Check, if the plate format has been provided
  if(missing(.plateFormat)){stop("Data missing! Please provide a plate format specifying the number of wells of a microplate to be dispensed into.", call. = TRUE)}
  if(!is.numeric(.plateFormat)){stop("in '.plateFormat'. Plate type needs to be a numeric object.\nPlease select a plate with 6, 12, 24, 48, 96, 384 or 1536 wells.", call. = TRUE)}
  if(!(.plateFormat %in% c(6, 12, 24, 48, 96, 384, 1536))){stop("Plate type not supported. Please select a plate with 6, 12, 24, 48, 96, 384 or 1536 wells.", call. = TRUE)}
  

  
  # SECTION A: PROBE DISPENSING #######################################################################################
  # Probe dispensing with a minimum set of provided data
  if (.probeDispensing){
  
    # Calculate the number of controls
    .noCtrls <- nrow(listofCtrls) * .ctrlReplicates + ifelse(is.null(.addUntreated), 0, .addUntreated$replicates)
    # Calculate the number of plates:
    # by number of wells needed / number of wells available (w/o excluded wells and wells for controls)
    .noPlates <- ceiling(nrow(listofCombinations) / (.plateFormat-length(listofExWells)-.noCtrls))
    
    # Estimate unique number of drug concentrations for single drug treatments and in combination
    .noDoses <- list(single = as.numeric(unique(with(subset(listofCombinations, Drug.1 == Drug.2), by(subset(listofCombinations, Drug.1 == Drug.2), INDICES = Drug.1, FUN = function(x){ length(unique(x$Dose.1)) })))), 
                     combination = unique(with(subset(listofCombinations, Drug.1 != Drug.2), apply(cbind(by(subset(listofCombinations, Drug.1 != Drug.2), INDICES = Drug.1, FUN = function(x){ length(unique(x$Dose.1)) }),
                                                                                                         by(subset(listofCombinations, Drug.1 != Drug.2), INDICES = Drug.2, FUN = function(x){ length(unique(x$Dose.2)) })), 1, FUN = max))) )
    
    # Print the dispensing summary for a potential dispensing
    message("Predicted dispensing summary:")
    cat("Number of drug treatments:",               format(nrow(listofCombinations), big.mark=" "), fill = TRUE)
    cat("Number of unique drug combinations:",      format(nrow(unique(subset(listofCombinations, Drug.1 != Drug.2))), big.mark=" "), fill = TRUE)
    cat("Number of unique single drug treatments:", format(nrow(unique(subset(listofCombinations, Drug.1 == Drug.2))), big.mark=" "), fill = TRUE)
    cat("Number of excluded wells per plate:",      format(length(listofExWells), big.mark=" "), fill = TRUE)
    cat("Number of controls per plate:",            format(.noCtrls, big.mark=" "), fill = TRUE)
    cat("Number of total plates:",                  format(.noPlates, big.mark=" "), fill = TRUE)
    cat("Number of drugs: ", with(listofCombinations, length(unique(c(Drug.1, Drug.2)))), ",  Number of doses: ", ifelse(!length(.noDoses$single), 0, .noDoses$single), "/", ifelse(!length(.noDoses$combination), 0, .noDoses$combination), sep = "", fill = TRUE)
    
    return(cat(""))
    
  }
  
  
  
  # Check further arguments required to generate the dispensing data set
  if(missing(listofDrugs)){stop("Data missing! Please provide a list of drugs.", call. = TRUE)}
  if(!is.list(listofDrugs)){
    stop("in 'listofDrugs'. Argument needs to be an object of class data.frame! Please provide a data frame containing a column with the drug name ['Name'], a unique drug number ['Number'] and the CAS number ['CAS number'].", call. = TRUE)
  } else if (!all(sapply(c("name", "number", "cas"), function(x) any(grepl(x, names(listofDrugs), ignore.case = TRUE))))){stop("in 'listofDrugs'. List of drugs has missing data! Please provide a data frame containing a column with the drug name ['Name'], a unique drug number ['Number'] and the CAS number ['CAS number'].", call. = TRUE)
  } else {
    colnames(listofDrugs)[grepl("name", names(listofDrugs), ignore.case = TRUE)] <- "NAME"
    colnames(listofDrugs)[grepl("number", names(listofDrugs), ignore.case = TRUE) & !grepl("cas", names(listofDrugs), ignore.case = TRUE)] <- "ID"
    colnames(listofDrugs)[grepl("cas", names(listofDrugs), ignore.case = TRUE)] <- "CAS_NUMBER"
  }
  
  
  if(missing(listofDoses)){stop("Data missing! Please provide a list of doses.", call. = TRUE)}
  if(!is.list(listofDoses)){
    stop("in 'listofDoses'. Argument needs to be an object of class data.frame! Please provide a data frame containing essential columns with the drug name ['Drug'], the drug dose ['Dose'] and the unit ['Unit'].", call. = TRUE)
  } else if(!all(sapply(c("drug", "dose", "unit"), function(x) any(grepl(x, names(listofDoses), ignore.case = TRUE))))){stop("in 'listofDoses'. List of doses has missing data! Please provide a data frame containing essential columns with the drug name ['Drug'], drug dose ['Dose] and the unit ['Unit'].", call. = TRUE)
  } else {
    colnames(listofDoses)[grepl("drug", names(listofDoses), ignore.case = TRUE)] <- "Drug"
    colnames(listofDoses)[grepl("dose", names(listofDoses), ignore.case = TRUE)] <- "Dose"
    colnames(listofDoses)[grepl("unit", names(listofDoses), ignore.case = TRUE)] <- "Unit"
  }
  
  
  if(missing(listofVolumes)){stop("Data missing! Please provide a list of volumes.", call. = TRUE)}
  if(!is.list(listofVolumes)){
    stop("in 'listofVolumes'. Argument needs to be an object of class data.frame! Please provide a data frame containing essential columns with the drug name ['Drug'], the volume ['Volume'] for each of the doses at which this drug is to be dispensed and the unit ['Unit'] of the dispensing volume.", call. = TRUE)
  } else if(!all(sapply(c("drug", "vol", "unit"), function(x) any(grepl(x, names(listofVolumes), ignore.case = TRUE))))){stop("in 'listofVolumes'. List of volumes has missing data! Please provide a data frame containing essential columns with the drug name ['Drug'], the volume ['Volume'] for each of the doses at which a drug is to be dispensed and the unit ['Unit'] of the dispensing volume.", call. = TRUE)
  } else if(sum(grepl("Vol", names(listofVolumes), ignore.case = TRUE)) != max(with(listofDoses, ave(seq_along(Drug), Drug, FUN=length)))){stop("in 'listofVolumes'. List is missing volumes for at least one dose! Please provide a column with the volume ['Volume'] for each of the doses at which a drug is to be dispensed.", call. = TRUE)
    colnames(listofVolumes)[grepl("drug", names(listofVolumes), ignore.case = TRUE)] <- "Drug"
    colnames(listofVolumes)[grepl("unit", names(listofVolumes), ignore.case = TRUE)] <- "Unit"
  }

  
  if(missing(listofStockConcentrations)){stop("Data missing! Please provide a list of drugs.", call. = TRUE)}
  if(!is.list(listofStockConcentrations)){
    stop("in 'listofStockConcentrations'. Argument needs to be an object of class data.frame! Please provide a data frame containing a column with the drug name ['Name'], a unique drug number ['Number'], the stock concentration of that drug and the concentration unit ['Unit'].", call. = TRUE)
  } else if (!all(sapply(c("name", "number", "conc", "unit"), function(x) any(grepl(x, names(listofStockConcentrations), ignore.case = TRUE))))){stop("in 'listofStockConcentrations'. List of drugs has missing data! Please provide a data frame containing a column with the drug name ['Name'], a unique drug number ['Number'], the stock concentration of that drug and the concentration unit ['Unit'].", call. = TRUE)
  } else {
    colnames(listofStockConcentrations)[grepl("name", names(listofStockConcentrations), ignore.case = TRUE)] <- "NAME"
    colnames(listofStockConcentrations)[grepl("number", names(listofStockConcentrations), ignore.case = TRUE) & !grepl("CAS", names(listofStockConcentrations), ignore.case = TRUE)] <- "ID"
    colnames(listofStockConcentrations)[grepl("conc", names(listofStockConcentrations), ignore.case = TRUE)] <- "CONCENTRATION"
    colnames(listofStockConcentrations)[grepl("unit", names(listofStockConcentrations), ignore.case = TRUE)] <- "UNIT"
  }
  

  if(missing(sourcePlate)){stop("Data missing! Please provide a source plate.", call. = TRUE)}
  if(!is.list(sourcePlate)){
    stop("in 'sourcePlate'. Argument needs to be an object of class data.frame! Please provide a data frame containing a column with the source plate id ['PlateID'], source plate well ['Well'], the drug ['Drug'] found in that well, the concentration of that drug ['Concentration'], the concentration unit ['Unit'].", call. = TRUE)
  } else if (!all(sapply(c("plateid", "well", "drug", "dose|conc", "unit"), function(x) any(grepl(x, names(sourcePlate), ignore.case = TRUE))))){stop("in 'sourcePlate'. List of drugs has missing data! Please provide a data frame containing a column with the source plate id ['PlateID'], source plate well ['Well'], the drug ['Drug'] found in that well, the concentration of that drug ['Concentration'], the concentration unit ['Unit'].", call. = TRUE)
  } else {
    colnames(sourcePlate)[grepl("plateid", names(sourcePlate), ignore.case = TRUE)] <- "PlateID"
    colnames(sourcePlate)[grepl("well", names(sourcePlate), ignore.case = TRUE)] <- "Well"
    colnames(sourcePlate)[grepl("drug", names(sourcePlate), ignore.case = TRUE)] <- "Drug"
    colnames(sourcePlate)[grepl("dose|conc", names(sourcePlate), ignore.case = TRUE)] <- "Concentration"
    colnames(sourcePlate)[grepl("unit", names(sourcePlate), ignore.case = TRUE)] <- "Unit"
  }
  
  
  # Check, if a final volume has been provided 
  if(missing(.finalWellVolume)){stop("Data missing! Please provide the final volume found in all wells.", call. = TRUE)}
  
  # Check, if a destination plate ID has been provided
  if(missing(.destinationPlateID)){stop("Data missing! Please provide a unique destination plate ID.", call. = TRUE)}
  
  # Check, if dispensing should be randomized
  if(missing(.randomizeDispensing)){message("Randomization argument not provided. Using default: Dispensing NOT randomized.")}

    
  
  # SECTION B: CREATE DIFFERENT SETS OF LISTS #########################################################################
  
  # Add the corresponding dose number to each dose of a drug
  doseList <- listofDoses
  doseList$DoseID <- with(listofDoses, rev(ave(seq_along(Drug), Drug, FUN=seq_along)))
  names(doseList) <- sub("Unit", "Unit (Dose)", names(doseList))
  
  # Convert list of volumes from wide to long-format
  volumeList <- stats::reshape(listofVolumes, varying = list(grepl("Vol", names(listofVolumes), ignore.case = TRUE)), v.names = c("Volume"), timevar = NULL, direction = "long", new.row.names = NULL)
  volumeList <- within(volumeList, rm("id"))
  rownames(volumeList) <- NULL
  volumeList$DoseID <- with(volumeList, rev(ave(seq_along(Drug), Drug, FUN=seq_along)))
  names(volumeList) <- sub("Unit", "Unit(Volume)", names(volumeList))
  
  # Create a reference data set with doses and dispensing volumes for each drug
  refVolDoseList <- merge(doseList, volumeList, by = intersect(names(doseList), names(volumeList)), all = TRUE)
  
  
  # Create a list of drugs, controls and stocks
  # Create a list with drug names
  drugList <- as.list(listofDrugs[,c("ID", "NAME", "CAS_NUMBER")])
  # Name the list of drugs based on their batch ID
  drugList$NAME <- setNames(drugList$NAME, drugList$ID)
  
  # Create a list of controls
  ctrlList <- as.list(listofCtrls[,c("NAME", "CAS_NUMBER", "DOSE", "UNIT", "VOLUME")])
  # Name the list of controls based on their name
  ctrlList$NAME <- setNames(ctrlList$NAME, ctrlList$ID)
  
  # Create a list of stock concentrations
  stocklList <- as.list(listofStockConcentrations[,c("ID", "NAME", "CONCENTRATION", "UNIT")])
  # Name the list of drugs based on their batch ID
  stocklList$NAME <- setNames(stocklList$NAME, stocklList$ID)
  

  
  

  
  #### SECTION C: EMPTY PLATE SET-UP ##################################################################################
  # Creates a plate with only the wells available for dispensing: w/o excluded wells
  
  # Define the row and column names based on plate type
  # Create a list of row and column names that match the raw labeling of the source plate
  baseplate <- baseplate(.plateFormat, zeroBase = FALSE, leadingZero = FALSE)
  

  # Create a reference list with only the wells that are available for dispensing
  dispWells <- data.frame(Number = baseplate[["wells"]]$no., Well = baseplate[["wells"]]$id., stringsAsFactors = FALSE)
  dispWells <- dispWells[!dispWells$Well %in% listofExWells,]
  
  
  
  #### SECTION D: GENERATE DISPENSING DATA ############################################################################
  
  # Create a column with a unique combination ID
  # this ID is supposed to track and identify each dispensing element throughout the generation of the dispensing file,
  # such as association of each drug with the source plate
  
  # Calculate the number of controls
  .noCtrls <- nrow(listofCtrls) * .ctrlReplicates + ifelse(is.null(.addUntreated), 0, .addUntreated$replicates)
  # Estimate the number of plates needed
  # Number of combinations divided by available number of wells per plate w/o outer wells
  .noPlates <- ceiling(nrow(listofCombinations) / (.plateFormat-length(listofExWells)-.noCtrls))

  
  # Check before proceeding, if any of the drug names contain a ":" which might cause an unexpected behavior in the section below
  if(any(grepl(":", drugList$NAME)) == TRUE){ 
    warning("One of the drug names contains an invalid character.", call. = FALSE, immediate. = TRUE)
    message("The following drugs contain incompatible characters:\n")
    message(paste(grep(":", drugList$NAME, value = TRUE), collapse="\n"), appendLF = FALSE)
  }else{
    message("Message: All drug names have passed the quality check.")
  }
  
  
  cat(" > Generating dispensing data...")
  
  # Reformat the doses from numeric to character
  # listofCombinations <- rapply(listofCombinations, f = format, scientific = FALSE, drop0trailing = TRUE, classes = c("numeric"), how = "replace")
  
  # Create a dispensing data set from the list of unique drug combinations 
  # Merge name, dose and unit for each drug pair
  dispensingData <- data.frame(Drug.1 = do.call(paste, c(listofCombinations[grep("1", names(listofCombinations), value=TRUE)], sep = ":")),
                               Drug.2 = do.call(paste, c(listofCombinations[grep("2", names(listofCombinations), value=TRUE)], sep = ":")),
                               stringsAsFactors = FALSE)
  
  
  # Distribute a set of combinations to a number of plates based on the number of available wells per plate
  dispensingData$Plate.Number <- rep(1:.noPlates, each=(.plateFormat-length(listofExWells)-.noCtrls), length.out=nrow(listofCombinations))
  # Randomize the drug treatments also across the plates, if requested
  dispensingData$Plate.Number <- if(.randomizeDispensing == TRUE){
    sample(dispensingData$Plate.Number)
    }else{ dispensingData$Plate.Number }
  # Arrange the data frame by plate number
  # The order of plates is not kept numerically (A1, A2, ... A24), but alphanumerically (A1, A10, A11 ... A9)
  # This can be changed by sorting plate numbers by factor
  dispensingData$Plate.Number <- factor(dispensingData$Plate.Number, paste(unique(gsub("[[:digit:]]", "", dispensingData$Plate.Number)), sort(order(unique(dispensingData$Plate.Number))), sep = ""))
  dispensingData <- dispensingData[order(dispensingData[["Plate.Number"]]),]

  # Add dispensing volume before separating drug combination to avoid double dispensing of the drug
  # Build data set by adding additional columns
  dispensingData[c("Combination ID", "Transfer Volume", "CAS number", "Source Plate Barcode", "Source Well", "Destination Well")] <- NA
  names(dispensingData) <- make.names(names(dispensingData), unique = TRUE)
  dispensingData[c("Destination.Plate.Barcode")] <- .destinationPlateID

  # Add a set of controls to each plate
  # NOTE: Each control is dispensed twice in combination with itself to obtain the same amount per well such as for the drug combination
  dispensingData <- rbind(dispensingData, data.frame(Drug.1 = paste(rep(ctrlList$NAME, .noPlates, each = .ctrlReplicates), rep(ctrlList$DOSE, .noPlates, each = .ctrlReplicates), rep(ctrlList$UNIT, .noPlates, each = .ctrlReplicates), sep = ":"),
                                                     Drug.2 = paste(rep(ctrlList$NAME, .noPlates, each = .ctrlReplicates), rep(ctrlList$DOSE, .noPlates, each = .ctrlReplicates), rep(ctrlList$UNIT, .noPlates, each = .ctrlReplicates), sep = ":"),
                                                     Plate.Number = rep(1:.noPlates, each = .ctrlReplicates*length(ctrlList$NAME)),
                                                     Combination.ID = NA, 
                                                     Transfer.Volume = NA,
                                                     CAS.number = rep(ctrlList$CAS_NUMBER, .noPlates, each = .ctrlReplicates),
                                                     Source.Plate.Barcode = NA,
                                                     Source.Well = NA,
                                                     Destination.Well = NA,
                                                     Destination.Plate.Barcode = .destinationPlateID,
                                                     check.names = TRUE, stringsAsFactors = FALSE))
  # Add untreated controls to each plate
  if(!is.null(.addUntreated)){
    dispensingData <- rbind(dispensingData, data.frame(Drug.1 = paste(rep(.addUntreated$name, .noPlates, each = .addUntreated$replicates), rep(0, .noPlates, each = .addUntreated$replicates), rep(NA, .noPlates, each = .addUntreated$replicates), sep = ":"),
                                                       Drug.2 = paste(rep(.addUntreated$name, .noPlates, each = .addUntreated$replicates), rep(0, .noPlates, each = .addUntreated$replicates), rep(NA, .noPlates, each = .addUntreated$replicates), sep = ":"),
                                                       Plate.Number = rep(1:.noPlates, each = .addUntreated$replicates),
                                                       Combination.ID = NA, 
                                                       Transfer.Volume = NA,
                                                       CAS.number = NA,
                                                       Source.Plate.Barcode = NA,
                                                       Source.Well = NA,
                                                       Destination.Well = NA,
                                                       Destination.Plate.Barcode = .destinationPlateID,
                                                       check.names = TRUE, stringsAsFactors = FALSE))
  }
  
  # Arrange dispensing data by plate number
  dispensingData <- dispensingData[order(dispensingData$Plate.Number),]
  
  # Add unique ID to each drug treatment
  dispensingData$Combination.ID <- seq(1, nrow(listofCombinations)+(.noCtrls*.noPlates), 1)
  # Randomize the list of combinations
  dispensingData$Combination.ID <- if(.randomizeDispensing == TRUE){
    sample(dispensingData$Combination.ID)
    }else{ dispensingData$Combination.ID }
  dispensingData <- dispensingData[order(dispensingData$Combination.ID),]

  
  # Assign destination wells to each drug treatment, either systematically or randomly
  # Split data set by plate and assign a well for each drug treatment to be dispensed into
  if(.randomizeDispensing == TRUE){
    dispensingData <- do.call(rbind, setNames(lapply(split(dispensingData, dispensingData$Plate.Number), function(x){ x$Destination.Well = sample(dispWells[["Well"]][1:nrow(x)]); x}), NULL))
  }else{
    dispensingData <- do.call(rbind, setNames(lapply(split(dispensingData, dispensingData$Plate.Number), function(x){ x$Destination.Well = dispWells[["Well"]][1:nrow(x)]; x}), NULL))
  }
  rownames(dispensingData) <- NULL

  # Separate drug pairs into individual rows retaining the unique drug combination number
  # Convert drug pairs from wide to long-format, in which each drug has its own row
  dispensingData <- stats::reshape(dispensingData, varying = list(grepl("Drug", names(dispensingData), ignore.case = TRUE)), v.names = c("Drug"), timevar = "Drug.Number", times = c("Drug.1", "Drug.2"), direction = "long", new.row.names = NULL)
  dispensingData <- within(dispensingData, rm("id"))
  # dispensingData <- dispensingData[order(dispensingData[["Combination.ID"]]),]
  rownames(dispensingData) <- NULL
  # Replace initial drug column with individual columns for the drug name, dose and unit 
  dispensingData <- cbind(within(dispensingData, rm("Drug")), setNames(data.frame(do.call('rbind', strsplit(as.character(dispensingData$Drug),':', fixed=TRUE)), stringsAsFactors = FALSE), c("Drug", "Drug.Concentration", "Unit")))
  dispensingData <- transform(dispensingData, Drug.Concentration = as.numeric(Drug.Concentration))
  dispensingData <- replace(dispensingData, dispensingData == "NA", NA)
  
  # Remove duplicate dispensing from single drug combinations
  # Note: Avoiding dispensing of the same drug into the same well twice
  dispensingData <- dispensingData[!duplicated(dispensingData[c("Drug", "Drug.Concentration", "Combination.ID", "Plate.Number", "Destination.Well")]), ]
  
  # Fill additional columns that are needed for the dispensing file
  dispensingData[["CAS.number"]] <- ifelse(is.na(dispensingData$CAS.number), drugList$CAS_NUMBER[match(dispensingData$Drug, drugList$NAME)], dispensingData$CAS.number)
  dispensingData[["Source.Plate.Barcode"]] <- sourcePlate$PlateID[match(dispensingData$Drug, sourcePlate$Drug)]
  
  # Calculate the transfer volume based on the stock and final concentration in relation to the total volume per well 
  # (`Transfer Volume` = finalVolume / ((stocklList$CONCENTRATION[match(`Sample Name`, stocklList$NAME)] * 1000) / `Drug Concentration`))
  # OR
  
  # Add transfer volume based on a reference list of dispensing volumes, either by merging the reference volume-dose data set or by
  # dispensingData <- merge(dispensingData, refVolDoseList[c("Drug", "Dose", "Volume")], by.x = c("Drug", "Drug.Concentration"), by.y = c("Drug", "Dose"), all.x = TRUE, all.y = FALSE, sort = FALSE)
  # matching the volume to drug and dose
  dispensingData <-  transform(dispensingData, Transfer.Volume = refVolDoseList$Volume[match(paste(Drug, as.numeric(format(round(Drug.Concentration, 12), nsmall = 12))), paste(refVolDoseList$Drug, as.numeric(format(round(refVolDoseList$Dose, 12), nsmall = 12))))])
  dispensingData <-  transform(dispensingData, Transfer.Volume = ifelse(Drug %in% ctrlList$NAME, ctrlList$VOLUME[match(paste(Drug, Drug.Concentration), paste(ctrlList$NAME, ctrlList$DOSE))], Transfer.Volume))
  
  # Reorder rows, and columns
  dispensingData <- dispensingData[with(dispensingData, order(Drug.Number, Combination.ID)),]
  dispensingData <- dispensingData[,c("Combination.ID", "Drug", "CAS.number", "Drug.Concentration", "Unit", "Transfer.Volume", "Source.Plate.Barcode",
                                      "Source.Well", "Destination.Well", "Destination.Plate.Barcode", "Plate.Number")]

  
  # Adding DMSO to each well that has only a single drug dispensed
  # Note: this is done to have the same DMSO concentration in each well,
  # since a combination treatment will have twice the DMSO concentration than a single drug treatment
  
  # select only single drug treatments
  addDMSO <- dispensingData[with(dispensingData, ave(Combination.ID, Combination.ID, FUN = length)) == 1, ]
  addDMSO <- subset(addDMSO, !(Drug %in% c(ctrlList$NAME, if(!is.null(.addUntreated)){.addUntreated$name})))
  # remove treatments that already have the full concentration of DMSO
  addDMSO <- subset(addDMSO, !(Transfer.Volume == 10))
  # add DMSO to those wells
  addDMSO <- transform(addDMSO, Drug = "DMSO", CAS.number = ctrlList$CAS_NUMBER[match("DMSO", ctrlList$NAME)],
                       Drug.Concentration = ctrlList$DOSE[match("DMSO", ctrlList$NAME)]/2,
                       Unit = ctrlList$UNIT[match("DMSO", ctrlList$NAME)],
                       Source.Plate.Barcode = sourcePlate$PlateID[match("DMSO", sourcePlate$Drug)])

  # Add additional DMSO dispensing to the main dispensing data set
  dispensingData <- rbind(dispensingData, addDMSO)
  rm(addDMSO)
  
  
  
  
  # Adding source well to the dispensing data: this part is done separately, but can be included above for convenience
  # NOTE: The stock concentration on the source plate might be in some cases different from the needed concentration on the destination plate
  # The final concentration used for the destination plate is related to the final volume of 5 μl and a dispensing volume of 2.5 nl
  # This relates in cases of combinations to a total volume per well of 5 μM divided by the volume (2.5 nl) of a single drug dispensed times the final concentration desired:
  # (((finalWellVolume) / transferVolume) * `Drug Concentration`)
  
  # The code below does not differentiate by source plate. This might lead to problems, if there are multiple source wells across different source plates for the same drug
  # Separate the controls and the treatments and assign the source well individually for each data set and merge both sets back together
  # Split the data set containing only the controls by each individual control and distribute multiple available source wells equally among the dispensings
  # Note: Each set of controls can only be on a single source plate and not across multiple plates
  # finalDispensingData <- do.call(rbind, setNames(c(lapply(split(dispensingData[dispensingData$Drug %in% ctrlList$NAME,], dispensingData$Drug[dispensingData$Drug %in% ctrlList$NAME]), 
  #                                                         function(x){ x$Source.Well = rep(sourcePlate$Well[sourcePlate$Drug == unique(x$Drug)], length.out = nrow(x)); x}),
  #                                                  # Extract only non-control treatments and assign the source well matching the name and drug concentration related to the volume 
  #                                                  lapply(list(subset(dispensingData, !Drug %in% ctrlList$NAME)), function(x){ 
  #                                                    x$Source.Well = sapply(1:nrow(x), function(y){ sourcePlate$Well[grep(x[y,"Drug"], sourcePlate$Drug, fixed = TRUE)][match(as.numeric(format(round(((finalVolume / x[y,"Transfer.Volume"]) * x[y,"Drug.Concentration"]), 9), nsmall = 9)), 
  #                                                                                                                                                                             as.numeric(format(round(sourcePlate$Concentration[sourcePlate$Drug == x[y,"Drug"]], 9), nsmall = 9)))] })
  #                                                    ; x})), NULL))

  # Split data set by source plate
  finalDispensingData <- split(dispensingData, dispensingData$Source.Plate.Barcode)
  # Assign the source well separately for the controls and the drug treatments,
  # since there might be multiple source wells for each control.
  for (sP in names(finalDispensingData)){
    
    # Select only the control treatments and assign the source well for each individual control
    ctrlDispensingData <- subset(finalDispensingData[[sP]], Drug %in% ctrlList$NAME)
    if (nrow(ctrlDispensingData) != 0){
      ctrlDispensingData <- by(ctrlDispensingData, INDICES = ctrlDispensingData$Drug, FUN = function(x){ x$Source.Well = rep(sourcePlate$Well[sourcePlate$Drug == unique(x$Drug)], length.out = nrow(x)); x  })
      ctrlDispensingData <- do.call(rbind, setNames(ctrlDispensingData, NULL))
    }
    
    # Select only the drug treatments and assign the source well for each individual drug (row-wise)
    drugDispensingData <- subset(finalDispensingData[[sP]], !Drug %in% ctrlList$NAME)
    if (nrow(drugDispensingData) != 0){
      drugDispensingData$Source.Well <- sapply(1:nrow(drugDispensingData), function(x){ 
        # Map the source well by matching the drug name and dose rowwise
        sourcePlate$Well[grep(drugDispensingData[x,"Drug"], sourcePlate$Drug, fixed = TRUE)][mapply(function(x, y) isTRUE(all.equal(x, y, tolerance = 5e-6)), ((.finalWellVolume / drugDispensingData[x,"Transfer.Volume"]) * drugDispensingData[x,"Drug.Concentration"]), sourcePlate$Concentration[grep(drugDispensingData[x,"Drug"], sourcePlate$Drug, fixed = TRUE)])]
      })
    }
    
    # Merge both data sets together
    finalDispensingData[[sP]] <- rbind(ctrlDispensingData, drugDispensingData)
    rm(sP, ctrlDispensingData, drugDispensingData)
  }
  
  finalDispensingData <- do.call(rbind, setNames(finalDispensingData, NULL))
  finalDispensingData <- rbind(finalDispensingData, subset(dispensingData, Drug == .addUntreated$name))
  
  # Reorder rows by Combination.ID
  finalDispensingData <- finalDispensingData[with(finalDispensingData, order(Combination.ID)),]
  rownames(finalDispensingData) <- NULL

  

  # Check before proceeding if any of the essential data for dispensing is missing before generating and exporting the dispensing file
  if(any(is.na(finalDispensingData[finalDispensingData$Drug != .addUntreated$name,!(names(finalDispensingData) == "CAS.number")]))){
    warning("One of the essential columns contains missing data.", call. = FALSE, immediate. = TRUE)
    message("Columns: ", paste(names(finalDispensingData)[sapply(finalDispensingData, function(x)any(is.na(x)))], collapse = ", "))
  }else{
    cat('\r', "Finished generating dispensing data. The dispensing data was successfully generated.")
  }
  
  
  
  #### SECTION E: CREATE OBJECT OF CLASS S3 ###########################################################################
  
  
  # Generate a dispensing summary for the dispensing file
  # Calculate the number of controls
  .noCtrls <- nrow(listofCtrls) * .ctrlReplicates + ifelse(is.null(.addUntreated), 0, .addUntreated$replicates)
  # Calculate the number of plates:
  # by number of wells needed / number of wells available (w/o excluded wells and wells for controls)
  .noPlates <- ceiling(nrow(listofCombinations) / (.plateFormat-length(listofExWells)-.noCtrls))
  
  # Estimate unique number of drug concentrations for single drug treatments and in combination
  .noDoses <- list(single = as.numeric(unique(with(subset(listofCombinations, Drug.1 == Drug.2), by(subset(listofCombinations, Drug.1 == Drug.2), INDICES = Drug.1, FUN = function(x){ length(unique(x$Dose.1)) })))), 
                   combination = unique(with(subset(listofCombinations, Drug.1 != Drug.2), apply(cbind(by(subset(listofCombinations, Drug.1 != Drug.2), INDICES = Drug.1, FUN = function(x){ length(unique(x$Dose.1)) }),
                                                                                                       by(subset(listofCombinations, Drug.1 != Drug.2), INDICES = Drug.2, FUN = function(x){ length(unique(x$Dose.2)) })), 1, FUN = max))) )

  .noTreatments <-     as.numeric(nrow(listofCombinations))
  .noCombTreatments <- as.numeric(nrow(unique(subset(listofCombinations, Drug.1 != Drug.2))))
  .noSingTreatments <- as.numeric(nrow(unique(subset(listofCombinations, Drug.1 == Drug.2))))
  .noDrugs <-          as.numeric(with(listofCombinations, length(unique(c(Drug.1, Drug.2)))))
  .noExWells <-        as.numeric(length(listofExWells))
  
  .summary <- list(controls     = .noCtrls,
                   plates       = .noPlates,
                   drugs        = .noDrugs,
                   doses        = .noDoses,
                   treatments   = .noTreatments,
                   combinations = .noCombTreatments,
                   single       = .noSingTreatments,
                   excluded     = .noExWells)
  class(.summary) <- "summary.dispensingData"

  
  # Create a list of data used to generate the dispensing data set
  .data <- list(listofCombinations        = listofCombinations,
                listofDrugs               = listofDrugs,
                listofDoses               = listofDoses,
                listofVolumes             = listofVolumes,
                listofCtrls               = listofCtrls,
                listofStockConcentrations = listofStockConcentrations,
                sourcePlate               = sourcePlate,
                listofExWells             = listofExWells,
                .ctrlReplicates           = .ctrlReplicates,
                .addUntreated             = .addUntreated,
                .finalWellVolume          = .finalWellVolume,
                .plateFormat              = .plateFormat,
                .destinationPlateID       = .destinationPlateID,
                .randomizeDispensing      = .randomizeDispensing)
  
  # Create a list of data lists
  .lists <- list(drugList = drugList,
                 doseList = doseList,
                 ctrlList = ctrlList)

  
  # Create an object of S3 class with a list from the attributes of all three vectors
  object <- list(output = finalDispensingData, origData = .data, dataList = .lists, summary = .summary, 
                 plateFormat = .plateFormat, dispensingID = .destinationPlateID, randomized = .randomizeDispensing)
  # Set the class of the object to be returned
  class(object) <- "dispensingData"
  

  return(object)

}
