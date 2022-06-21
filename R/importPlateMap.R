#' Import plate maps from IncuCyte® .PlateMap files, Microsoft Excel or .CSV files
#'
#' @description
#' A component of the package \pkg{PlateMap}, useful for import, conversion and modification of plate maps to be used in a drug sensitivity screen.
#' The function \emph{\code{importPlateMap}} can be used for importing \emph{IncuCyte® .PlateMap} files, either as such or converted into a source plate format. Alternatively, plate maps can be imported from Microsoft Excel, or as comma-separated values (csv) files. 
#' The source plate used in a drug sensitivity screen contains drugs at stock concentrations from which desired concentrations are dispensed onto destination plates, which are used to run the drug sensitivity screen on. 
#' The source plate serves as a reference for the location of drugs at their stock concentration used in dispensing files.
#'
#' @param importFile \code{character}; a character vector with either the path to a folder containing multiple or individual plate map files. 
#' @param .fileFormat \code{character}; a predefined character vector stating the file format of the file(s) to import: .PlateMap, .xlsx or .csv
#' @param .plateMapGrp \code{list}; a list allowing to group individual components of the plate map and to generate a barcode for individual plates through the labels provided.
#' @param .sourcePlateConv \code{logical}; if TRUE, then the plate map will be converted to a source plate format. The default is FALSE, in which the plate map is imported in it's original format.
#'
#'
#' @details Drug sensitivity screens that use a library of drugs stored across multiple microplates, require the import of multiple plate maps, one for each microplate. 
#' Those plate maps can be directly imported from a comma separated values (.csv) file, or by utilizing a third-party software, such as the IncuCyte® Plate Map Editor. 
#' Third-party applications and tools usually work with proprietary file formats, which are not directly usable with R. The function \emph{importPlateMap} from the package \pkg{metascreen} is able to read and import the 
#' proprietary .PlateMap file format. This allows the user to design source plates using the IncuCyte® Plate Map Editor and import the plate map into their pipeline.
#' 
#' The \emph{IncuCyte® .PlateMap} files are in a proprietary file format of the \strong{IncuCyte® Plate Map Editor}.
#' 
#' The function prints the progress of the plate map conversion on screen, along with a short summary of all compounds and cell lines that have been found.
#' 
#' The plate maps can be imported either as single files, or from multiple files. It is recommended to specify the file format \emph{\code{.fileFormat}} of the file to import, especially if the folder from which the files are imported contain different file formats.
#' Alternatively, the function is able to detect any of the supported file formats automatically and import them. Note: for the time being, it is not possible to import or combine plate maps from different file formats at once, but rather one file format at a time.  
#' 
#' Please note that this function allows the import plate maps with multiple components, such as \emph{compounds}, \emph{cells} and/or \emph{other}. However, doing so, requires a specific file nomenclature, if the plate maps are being imported from individual .csv files.
#' 
#' \strong{compounds}, or drugs need to come from files that are named or contain at least the word \emph{drug} in the file name. The corresponding file with the drug concentrations, needs to be named or contain the word \emph{concentrations}.
#' If multiple drugs are found in individual wells, such as drug combinations, the drug pairs need to be imported from individual .csv files and labeled with a numeric suffix. The same applies to the concentration files.
#' 
#' For example: If two drugs are to be combined in the same well of a single plate, the files need to be named drugs-1.csv, drugs-2.csv, concentrations-1.csv and concentrations-2.csv.
#' 
#' This can be extended with any number of drugs. If no concentrations are found, only the drugs will be imported. However, it is not possible to import concentrations without the corresponding drug files.
#' 
#' \strong{cells}, or cell lines can be imported from .csv files named or containing the word \emph{cells}. If multiple cell lines are to be used, the files need to be labeled with a numeric suffix.
#' For example: cells-1.csv, cells-2.csv ...
#' 
#' 
#' The \emph{.plateMapGrp} allows to group individual components where necessary. A case scenario would be, in which multiple cell lines need to be imported that share the same drug treatment, but are found across different plates, or in which multiple cell lines are found across a single plate.
#' The benefit of grouping is that it does not require to import the same compounds multiple times for each cell line, but rather will replicate components that are not grouped for each grouped component.
#' 
#' For example: (a) If the same drugs and drug concentrations are to be used for multiple cell lines, the cells can be grouped, while the drugs left ungrouped. This will replicate the same drug treatment for each of the cell lines.
#' (b) If multiple cell lines are found across a single plate. in case where multiple cell lines are spread across a number of plates, all cell lines found on an individual plate can be assigned to the same group. The drug treatment will then be replicated across that group of cell lines.
#' Please note that only one component can be grouped at a time.
#' 
#' The grouping can be extended for the generation of plate barcodes. Instead of assigning individual grouping labels, a unique barcode can be assigned to each group, denoting the plate map to a single plate.
#' 
#' This does not apply to plate maps imported from .PlateMap files. Due to the internal standard of the file format, multiple components can be imported simultaneously through a single file. However, the name of the .PlateMap file will denote the plate id, or plate barcode.
#' 
#' 
#' The option \emph{\code{.sourcePlateConv}} will automatically convert the plate map to the requirement of a source plate. If this option is not set, the full data set is returned, including empty wells, or missing components.
#' If the plate map is being converted to a source plate format, the empty data will be removed as well as components that are not compounds, which can then be used as a reference for the generation of dispensing files using \code{generateDispensingData()}. 
#' 
#' @return Returns an object of class "data.frame", with essential columns containing the plate id ['plateid'] for each source plate, the location ['well'] of the drugs ['drug'] and their concentrations ['concentration'] along with the units ['unit'].
#'
#' @examples
#' \donttest{\dontrun{
#' # Import a .PlateMap file and convert plate map into a source plate format
#' sourcePlate <- importPlateMap("myPlateMapFile.PlateMap", .sourcePlateConv=TRUE)
#' 
#' # Importing plate map files from folder
#' importPlateMap("path/to/folder/", sourcePlateConv=FALSE)
#' 
#' # Importing plate maps from  .PlateMap, .xlsx or .csv files files
#' importPlateMap("path/to/folder/", .fileFormat = "PlateMap", sourcePlateConv=FALSE)
#' importPlateMap("path/to/folder/", .fileFormat = ".xlsx", sourcePlateConv=FALSE)
#' importPlateMap("path/to/folder/", .fileFormat = "csv", sourcePlateConv=FALSE)
#' 
#' # Grouping of cell lines across single and multiple plates
#' # Note: First two cell lines are grouped on the same plate, while the last two cell lines on another plate
#' importPlateMap("path/to/folder/", .fileFormat = ".csv", 
#'                .plateMapGrp = list(type = "Cells", name = c("A375", "WM1366", "WM45.1", "FEMXV"), 
#'                                    group = c("A", "A", "B", "B")))
#'                                    
#' # Extending labels with plate barcodes, same set two plates
#' importPlateMap("path/to/folder/", .fileFormat = ".csv", 
#'                .plateMapGrp = list(type = "Cells", name = c("A375", "WM1366", "WM45.1", "FEMXV"), 
#'                                    group = c("0921A1", "0921A1", "0921A2", "0921A2")))
#' 
#' # Extending labels with plate barcodes, two different sets of plates
#' importPlateMap("path/to/folder/", .fileFormat = ".csv", 
#'                .plateMapGrp = list(type = "Cells", name = c("A375", "WM1366", "WM45.1", "FEMXV"), 
#'                                    group = c("0921A1", "0921A1", "0921B1", "0921B1")))
#' 
#' }}
#' 
#' @keywords incucyte plate map platemap
#' 
#' 
#' @importFrom XML xmlParse xpathSApply xmlGetAttr xmlSize
#' @importFrom utils file_test type.convert read.csv
#' @importFrom stats setNames
#'
#' @export

importPlateMap <- function(importFile, .fileFormat, .plateMapGrp = list(type, name, group), .sourcePlateConv=FALSE){

  
  # Check, if arguments are provided and in the proper format
  if(missing(importFile)){stop("File missing! Please provide a plate map file to import.", call. = TRUE)}
  
  if(missing(.fileFormat)){
    message("missing: '.fileFormat'. File format not provided! Checking file format of provided files...", appendLF = TRUE)
  } else if (grepl("PlateMap|xlsx|csv", .fileFormat)) {
    
    # Format .fileFormat
    if(grepl("IncuCyte|PlateMap", .fileFormat, ignore.case = TRUE)){
      .fileFormat <- paste(".", gsub("\\.", "", "PlateMap"), sep = "")
    } else if(grepl("Microsoft|Excel|xlsx", .fileFormat, ignore.case = TRUE)){
      .fileFormat <- paste(".", gsub("\\.", "", ".xlsx"), sep = "")
    } else if(grepl("csv", .fileFormat, ignore.case = TRUE)){
      .fileFormat <- paste(".", gsub("\\.", "", tolower(.fileFormat)), sep = "")
    }
    
  } else {
    stop("in '.fileFormat'. File format not supported! Please select one of the following file formats: .PlateMap, .xlsx, or .csv", call. = TRUE)
    }
  
  
  # Check whether a file or file path was provided
  # In case of a file, check whether the file format is supported
  # otherwise in case of a directory, search for all supported files within the directory
  if(utils::file_test("-f", importFile)){
    if(!grepl("PlateMap|xlsx|csv", tail(unlist(strsplit(basename(importFile), ".", fixed = TRUE)), 1), ignore.case = TRUE)){
      stop("in 'importFile'. File format not supported! Please provide a file in the following format: .PlateMap, .xlsx, or .csv", call. = TRUE)
    } else { fileList <- importFile }
  } else if(utils::file_test("-d", importFile)){
    
    # Fetch all files inside the directory
    fileList <- list.files(file.path(importFile), full.names = TRUE, ignore.case = TRUE, recursive = TRUE)
    
    # Check, if any of the files are in the supported file format, or 
    # if multiple file formats are found, while not file format has been specified
    if(all(!grepl("PlateMap|xlsx|csv", sapply(strsplit(fileList, split="\\."), tail, 1L), ignore.case = TRUE))){
      stop("in 'importFile'. File not found! Please provide a file in the following format: .PlateMap, .xlsx, or .csv", call. = TRUE)
    } else if (all(length(unique(sapply(strsplit(fileList, split="\\."), tail, 1L))) != 1, missing(.fileFormat))) {
      stop("in 'importFile'. Multiple file formats found inside directory. Please provide a file format to import.", call. = TRUE)
    } else if (missing(.fileFormat)) { .fileFormat <- paste(".", tolower(unique(sapply(strsplit(fileList, split="\\."), tail, 1L))), sep = "") }
    
    fileList <- list.files(file.path(importFile), pattern = .fileFormat, full.names = TRUE, ignore.case = TRUE, recursive = TRUE)
    
  } else {
    stop("File or directory does not exist! Please provide a file or a file path to one or more supported files.", call. = TRUE)
  }
  
  
  # Check, if a plate map group has been provided
  if (missing(.plateMapGrp)){ 
    .plateMapGrp = FALSE 
  } else if (class(.plateMapGrp) != "list") { stop("in '.plateMapGrp'. Argument needs to be a list! Please provide a list with the type and name of individual components to group by, along with the labels for those groups.", call. = TRUE) 
  } else if (length(.plateMapGrp) == 0){ stop("in '.plateMapGrp'. List cannot be empty! Please provide a list with the type and name of individual components to group by, along with the labels for those groups.", call. = TRUE)
  } else if (!("type" %in% names(.plateMapGrp))){ stop("in '.plateMapGrp'. List is missing an element! Please provide the 'type' of compnent to group by: 'cell', 'compound', or 'other'", call. = TRUE)
  } else if (!("name" %in% names(.plateMapGrp))){ stop("in '.plateMapGrp'. List is missing an element! Please provide the 'name' of each compnent to group by.", call. = TRUE)
  } else if (!("group" %in% names(.plateMapGrp))){ stop("in '.plateMapGrp'. List is missing an element! Please provide the 'group' labels for each of the components.", call. = TRUE)
  } else if (length(.plateMapGrp[["name"]]) > length(.plateMapGrp[["group"]])){ stop("in '.plateMapGrp'. Labels missing for one or more names. Please provide for each name a corresponding label.", call. = TRUE)
  } else if (grepl("cell|comp|other", .plateMapGrp[["type"]], ignore.case = TRUE)) {
    if(grepl("cell", .plateMapGrp[["type"]], ignore.case = TRUE)){
      .plateMapGrp[["type"]] <- "Cell"
    } else if(grepl("comp", .plateMapGrp[["type"]], ignore.case = TRUE)){
      .plateMapGrp[["type"]] <- "Compound"
    } else if(grepl("other", .plateMapGrp[["type"]], ignore.case = TRUE)){
      .plateMapGrp[["type"]] <- "Other"
    }
  }
  
  
  
  cat("\nA total of", length(fileList), "files have been found:", paste(basename(fileList), collapse = ", "), "\n\n")
  
  
  # Function for importing IncuCyte® .PlateMap files
  .PlateMap <- function(fileList){
    
    for(platemapFile in fileList) {
      
      # Import XML file into object an convert it to a list
      xmlData <- XML::xmlParse(file=file.path(platemapFile))
      # xmlData <- xmlToList(xmlData, addAttributes = TRUE, simplify = FALSE)
      
      
      # Display the compounds and cell lines found in the plate map
      # check, if any compounds are found and how many
      if(xpathSApply(xmlData, sprintf("boolean(//*/referenceItems/referenceItem[@type='Compound'])"))){
        cat("A total of ", length(xpathSApply(xmlData, "//*/referenceItems/referenceItem[@type='Compound']")), " compounds found in ", basename(platemapFile), ".\n", sep = "")
        cat("Compounds found: ", 
            paste(sort(xpathSApply(xmlData, sprintf("//*/referenceItems/referenceItem[@type='Compound']"), xmlGetAttr, "displayName", default = NA)), collapse = ", "), "\n", sep = "")
      } else { cat("No compounds found in ", basename(platemapFile), "\n", sep = "") }
      # check, if any cell types are found and how many
      if(xpathSApply(xmlData, sprintf("boolean(//*/referenceItems/referenceItem[@type='CellType'])"))){
        cat("A total of ", length(xpathSApply(xmlData, "//*/referenceItems/referenceItem[@type='CellType']")), " cell lines found in ", basename(platemapFile), ".\n", sep = "")
        cat("Cell lines found: ", 
            paste(sort(xpathSApply(xmlData, sprintf("//*/referenceItems/referenceItem[@type='CellType']"), xmlGetAttr, "displayName", default = NA)), collapse = ", "), "\n", sep = "")
      } else { cat("No cell lines found in ", basename(platemapFile), ".\n", sep = "") }
      cat('\n')
      
      
      platemap = data.frame()
      # Extract individual attributes from the plate map into a single data frame
      for (well in 1:xpathSApply(xmlData, "//wellStore/wells", xmlSize)){
        cat('\r', "> Reading PlateMap file. Processing well ", well, " of a ", xpathSApply(xmlData, "//wellStore/wells", xmlSize), "-well plate.", sep = "")
        platemap <- rbind(platemap, data.frame(No. = well, Row = xpathSApply(xmlData, sprintf("//*/well[%d]/@row", well)), Column = xpathSApply(xmlData, sprintf("//*/well[%d]/@col", well)),
                                               
                                               # Check if node exists and assign its attribute value
                                               Type = if(!xpathSApply(xmlData, sprintf("boolean(//*/well[%d]/items/wellItem)", well))){NA}else{xpathSApply(xmlData, sprintf("//*/well[%d]/items/wellItem/referenceItem", well), xmlGetAttr, "type", default = NA)},
                                               Name = if(!xpathSApply(xmlData, sprintf("boolean(//*/well[%d]/items/wellItem)", well))){NA}else{xpathSApply(xmlData, sprintf("//*/well[%d]/items/wellItem/referenceItem", well), xmlGetAttr, "displayName", default = NA)},
                                               Description = if(!xpathSApply(xmlData, sprintf("boolean(//*/well[%d]/items/wellItem)", well))){NA}else{xpathSApply(xmlData, sprintf("//*/well[%d]/items/wellItem/referenceItem", well), xmlGetAttr, "description", default = NA)},
                                               
                                               # Extract attribute from node otherwise assign NA
                                               Concentration = if(!xpathSApply(xmlData, sprintf("boolean(//*/well[%d]/items/wellItem)", well))){NA}else{xpathSApply(xmlData, sprintf("//*/well[%d]/items/wellItem", well), xmlGetAttr, "concentration", default = NA)},
                                               Unit = if(!xpathSApply(xmlData, sprintf("boolean(//*/well[%d]/items/wellItem)", well))){NA}else{xpathSApply(xmlData, sprintf("//*/well[%d]/items/wellItem", well), xmlGetAttr, "concentrationUnits", default = NA)},
                                               Passage = if(!xpathSApply(xmlData, sprintf("boolean(//*/well[%d]/items/wellItem)", well))){NA}else{xpathSApply(xmlData, sprintf("//*/well[%d]/items/wellItem", well), xmlGetAttr, "passage", default = NA)},
                                               Density = if(!xpathSApply(xmlData, sprintf("boolean(//*/well[%d]/items/wellItem)", well))){NA}else{xpathSApply(xmlData, sprintf("//*/well[%d]/items/wellItem", well), xmlGetAttr, "seedingDensity", default = NA)},
                                               
                                               colorArgb = if(!xpathSApply(xmlData, sprintf("boolean(//*/well[%d]/items/wellItem)", well))){NA}else{xpathSApply(xmlData, sprintf("//*/well[%d]/items/wellItem/referenceItem", well), xmlGetAttr, "colorArgb", default = NA)},
                                               
                                               row.names = NULL, check.names = FALSE, stringsAsFactors = FALSE), stringsAsFactors = FALSE)
        
        if(well == xpathSApply(xmlData, "//wellStore/wells", xmlSize)){
          cat('\r' , strrep(" ", 100), '\n', sep = "")
          rm(well)
        }
        
      }
      
      
      # Modify individual attributes post-extraction
      # correct the data type of each column
      platemap$No. <- as.numeric(platemap$No.)
      # correct for zero base offset
      platemap$Row <- as.numeric(platemap$Row) + 1
      # convert row numbers to names
      platemap$Row <- c(LETTERS, sapply(LETTERS, function(x) paste0(x, LETTERS)))[platemap$Row]
      platemap$Column <- as.numeric(platemap$Column) + 1
      # Map well numbers to row names and column numbers
      platemap$Well <- paste(platemap$Row, platemap$Column, sep = "")
      platemap$Concentration <- as.numeric(platemap$Concentration)
      platemap$Passage <- as.numeric(platemap$Passage)
      platemap$Density <- as.numeric(platemap$Density)
      platemap$colorArgb <- as.numeric(platemap$colorArgb)
      # Add new column with plate map id
      platemap$PlateID <- gsub(".PlateMap", "", basename(platemapFile), ignore.case = TRUE)
      # reorder columns
      platemap <- platemap[,c('No.', 'Well', 'Row', 'Column', 'Type', 'Name', 'Description', 'Concentration', 'Unit', 'Passage', 'Density', 'PlateID')]
      
      if(!exists(".sourcePlate", inherits = FALSE)){.sourcePlate = data.frame()}
      # Combine individual source plates together
      .sourcePlate <- rbind(.sourcePlate, platemap)
      
    }
    
    cat('\r' , "Finished importing PlateMap file(s).", strrep(" ", 100), '\n', sep = "")
    
    return(.sourcePlate)
    
  }
  
  
  
  # Function for importing Microsoft Excel, .xlsx files
  .xlsx <- function(){
    print("Work in progress.")
  }
  
  
  
  # Function for importing comma-separated values, .csv files
  .csv <- function(fileList){
    
    dfs = list()
    
    # Fetch one file at a time
    for(platemapFile in fileList) {
      
      # Check, if file contains compound data: drugs, concentrations, units
      if(grepl("drug", basename(platemapFile), ignore.case = TRUE)){
        
        cat("Reading plate map from file:", basename(platemapFile), sep = " ")
        
        # Import file containing drug labels
        compd.drug <- utils::read.csv(file=file.path(platemapFile), check.names=FALSE, header=TRUE, stringsAsFactors=FALSE, na.strings="", sep=",", dec=".", row.names = 1, skip=0)
        rownames(compd.drug) <- toupper(rownames(compd.drug))
        
        # Converting plate format from a matrix layout to a list or data frame
        compd.drug <- data.frame(Well = paste(rownames(compd.drug)[row(compd.drug)], colnames(compd.drug)[col(compd.drug)], sep = ""), Name = unlist(compd.drug), stringsAsFactors = FALSE)
        
        # Set factors for individual columns to allow proper sorting
        compd.drug$Well <- factor(compd.drug$Well, levels = with(compd.drug, unique(Well[order(nchar(gsub("[[:digit:]]", "", Well)), gsub("[[:digit:]]", "", Well), as.numeric(gsub("[[:alpha:]]", "", Well)))])))
        compd.drug <- compd.drug[with(compd.drug, order(Well)),]
        
        
        # Extend list with additional information
        compd.drug$No. <- 1:nrow(compd.drug)
        compd.drug$Row <- with(compd.drug, gsub("\\d", "", Well))
        compd.drug$Column <- with(compd.drug, gsub("\\D", "", Well))
        compd.drug$Type <- ifelse(is.na(compd.drug$Name), NA, "Compound")
        compd.drug$PlateID <- NA
        
        
        # Extract suffix from file name
        .suffix <- tail(unlist(regmatches(basename(platemapFile), regexec(".*\\D(\\d+).*", basename(platemapFile)))), 1L)
        .suffix <- as.character(ifelse(!length(.suffix), 0, .suffix))
        
        
        
        # Select the matching file with the drug concentrations
        platemapFile <- file.path(dirname(platemapFile), gsub("drug", "concentration", basename(platemapFile)))
        
        if(utils::file_test("-f", platemapFile)){
          cat(" Found matching concetrations file:", basename(platemapFile), '\n', sep = " ")
        } else {
          cat('\n')
          message("Warning: Matching concentrations file not found.", '\n', sep = " ")
          dfs[["compound"]][[.suffix]] <- compd.drug
          next
        }
        
        
        # Import file containing drug concentrations
        compd.conc <- utils::read.csv(file=file.path(platemapFile), check.names=FALSE, header=TRUE, stringsAsFactors=FALSE, na.strings="", sep=",", dec=".", row.names = 1, skip=0)
        rownames(compd.conc) <- toupper(rownames(compd.conc))
        
        # Converting plate format from a matrix layout to a list or data frame
        compd.conc <- data.frame(Well = paste(rownames(compd.conc)[ row(compd.conc)], colnames(compd.conc)[col(compd.conc)], sep = ""), Concentration = unlist(compd.conc), stringsAsFactors = FALSE)
        
        # Separate dose and unit
        compd.conc <- setNames(data.frame(compd.conc$Well, as.numeric(regmatches(compd.conc$Concentration, regexec("[[:digit:]]+\\.*[[:digit:]]*", compd.conc$Concentration))), 
                                          as.character(gsub("([.,-])|[[:digit:]]|\\s","", compd.conc$Concentration)), stringsAsFactors = FALSE), paste(c("Well", "Concentration", "Unit"), sep = ""))
        
        # Merge the list of drugs with the list of doses
        compound <- merge(compd.drug, compd.conc, by = intersect(names(compd.drug), names(compd.conc)), all = TRUE)
        compound <- compound[with(compound, order(Well)),]
        
        dfs[["compound"]][[.suffix]] <- compound
        
        
      }
      
      
      # Check, if file contains cell data: cell line name, cell number, cell passage number 
      if(grepl("cell", basename(platemapFile), ignore.case = TRUE)){
        
        cat("Reading plate map from file:", basename(platemapFile), '\n', sep = " ")
        
        
        # Import file containing drug labels
        cell.type <- utils::read.csv(file=file.path(platemapFile), check.names=FALSE, header=TRUE, stringsAsFactors=FALSE, na.strings="", sep=",", dec=".", row.names = 1, skip=0)
        rownames(cell.type) <- toupper(rownames(cell.type))
        
        # Converting plate format from a matrix layout to a list or data frame
        cell.type <- data.frame(Well = paste(rownames(cell.type)[row(cell.type)], colnames(cell.type)[col(cell.type)], sep = ""), Name = unlist(cell.type), stringsAsFactors = FALSE)
        
        # Set factors for individual columns to allow proper sorting
        cell.type$Well <- factor(cell.type$Well, levels = with(cell.type, unique(Well[order(nchar(gsub("[[:digit:]]", "", Well)), gsub("[[:digit:]]", "", Well), as.numeric(gsub("[[:alpha:]]", "", Well)))])))
        cell.type <- cell.type[with(cell.type, order(Well)),]
        
        
        # Extend list with additional information
        cell.type$No. <- 1:nrow(cell.type)
        cell.type$Row <- with(cell.type, gsub("\\d", "", Well))
        cell.type$Column <- with(cell.type, gsub("\\D", "", Well))
        cell.type$Type <- ifelse(is.na(cell.type$Name), NA, "Cell")
        cell.type$Passage <- cell.type$Density <- NA
        
        
        # Extract suffix from file name
        .suffix <- tail(unlist(regmatches(basename(platemapFile), regexec(".*\\D(\\d+).*", basename(platemapFile)))), 1L)
        .suffix <- ifelse(!length(.suffix), 0, .suffix)
        
        dfs[["cell-type"]][[.suffix]] <- cell.type
        
      }
      
    }
    
    cat('\n' , "Finished importing plate map from csv file(s).", strrep(" ", 100), '\n', sep = "")
    
    
    # Merge all (sub-lists) components of the same category 
    dfs <- lapply(dfs, function(x) do.call(rbind, setNames(x, NULL)))
    
    # Merge the list of each category to build a complete plate map
    # :rbind data frames with different columns
    platemap <- do.call(rbind, c(lapply(dfs, function(x) data.frame(c(x, sapply(setdiff(unique(unlist(lapply(dfs, names))), names(x)), function(y) NA)))), make.row.names=FALSE))
    
    # Sort rows and columns by name
    platemap <- platemap[with(platemap, order(Well)),]
    platemap <- platemap[,intersect(c("No.", "Well", "Row", "Column", "Type", "Name", "Description", "Concentration", "Unit", "Passage", "Density", "PlateID"), names(platemap))]
    rownames(platemap) <- NULL
    
    return(platemap)
    
  }
  
  
  
  # Select the function to import the plate map based on the file format
  .sourcePlate <- switch(.fileFormat, .PlateMap = .PlateMap(fileList), .xlsx = .xlsx(fileList), .csv = .csv(fileList))
  
  

    # Assign individual plat maps to a specific groups, where individual components are found on one plate, but not on another
  if(!is.logical(.plateMapGrp)){
    for(name in .plateMapGrp[["name"]]){
      .sourcePlate$PlateID <- with(.sourcePlate, replace(PlateID, (Type == .plateMapGrp[["type"]] & Name == name ), .plateMapGrp[["group"]][match(name, .plateMapGrp[["name"]])] ))
    }
    
    # Replicate each un-grouped entry  for each group
    .sourcePlate <- rbind(.sourcePlate, .sourcePlate[rep(rownames(.sourcePlate[!is.na(.sourcePlate$Name) & is.na(.sourcePlate$PlateID),]), each = length(unique(.plateMapGrp[["group"]]))-1), ])
    # Assign the group to each individual entry
    .sourcePlate[!is.na(.sourcePlate$Name) & is.na(.sourcePlate$PlateID),]$PlateID <- rep(unique(.plateMapGrp[["group"]]), each = nrow(.sourcePlate[!is.na(.sourcePlate$Name) & is.na(.sourcePlate$PlateID),])/length(unique(.plateMapGrp[["group"]])))
    
    # Sort rows and columns by name
    .sourcePlate <- .sourcePlate[with(.sourcePlate, order(Well, Type, PlateID)),]
    .sourcePlate <- .sourcePlate[,intersect(c("No.", "Well", "Row", "Column", "Type", "Name", "Description", "Concentration", "Unit", "Passage", "Density", "PlateID"), names(.sourcePlate))]
    rownames(.sourcePlate) <- NULL
  }
  
  
  # Select only essential columns for a source plate and drop all other columns
  if(.sourcePlateConv == TRUE){ 
    
    # Filter out empty well, in which the names is NA, rename column containing drug names
    # and select only essential columns
    .sourcePlate <- subset(.sourcePlate, !is.na(Name) & Type == "Compound")
    names(.sourcePlate) <- sub("Name", "Drug", names(.sourcePlate))
    .sourcePlate <- .sourcePlate[grepl("Well|Drug|Concentration|Unit|PlateID", names(.sourcePlate), ignore.case = TRUE)]
    
  }
  
  
  # Remove duplicate wells that are empty
  .sourcePlate <- unique(.sourcePlate)
  
  
  return(.sourcePlate)
  
}
