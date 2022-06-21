#' Quality control (QC) for drug sensitivity screens 
#'
#' @description
#' An essential component of the modular library \pkg{metascreen}, and imperative for the quality assurance of a drug sensitivity screen.
#' \emph{\code{qc}} is a function that offers a set of quality control analysis tools assessing the variance of individual controls, reporting the z'-factor between the positive and negative controls,
#' as well as looking at the signal of empty and untreated wells.
#' 
#' @param consolidatedData an object of class 'consolidatedData'.
#' @param .saveto string; path to a folder location where the object is saved to.
#' @param .ctrls \code{vector}; a list of controls to analyze.
#' @param .qcMethod \code{character}; a set of predefined quality control methods.
#' 
#' @details The function \code{qc} is used to assess the quality of a drug sensitivity screen by looking at the variance and signal distribution between individual controls.
#' 
#' With the parameter \emph{.qcMethod}, it is possible to choose between individual quality assessments. At the moment the following quality control methods are available:
#' 
#' \strong{variance} : assessing the variance between individual controls both, across all plates, as well as by individual plate
#' \strong{emptywells} : assessing the signal of empty and untreated wells, this will also include any excluded wells
#' \strong{firstcolumn} : assessing the signal of wells in the first column of each plate
#' \strong{zprime} : assessing the z'-factor based on the distribution between the positive and negative control. 
#' 
#' Note: If \emph{zprime} is selected as a method, only two controls can be used for the assessment. The z'-factor is being calculated for each set of controls on each plate.
#' In any case where more than two controls are provided with \emph{.ctrls}, only the first two controls will be used as positive and negative controls. Please set the order of controls accordingly.
#'  
#' The z'-factor is a metric providing a measure of quality for high-throughput screens. The calculations are based on the z'-factor described in the
#' original paper by Zhang, 1999, J Biomol Screen (see references).
#'  
#'  
#' It is possible to run multiple QC methods at once, or all by simply specifying \strong{all} in \emph{.qcMethod}.
#' 
#' Any plots generated during the QC analysis will be saved to the location provided with \emph{.saveto}.
#'   
#' If a set of control plates need to be assessed, this can be achieved with the functions \code{\link{readFiles}} and consequently plotting them with \code{\link{plot.readFiles}}.
#'  
#' @return Returns an object of class S3:controlData. THe object will contain the data and plots generated with the analysis.
#' 
#' @references Zhang, J. H., Chung, T. D., Oldenburg, K. R. (1999): A Simple Statistical Parameter for Use in Evaluation and Validation of High Throughput Screening Assays. J Biomol Screening 4 (2), S. 67–73. DOI: 10.1177/108705719900400206
#' 
#' @seealso \code{\link{consolidateData}}
#'
#' @examples
#' \donttest{\dontrun{
#' # Assessing the variance between individual controls
#' qc(consolidatedData, .saveto = "path/to/folder/", .ctrls = c("BzCl", "DMSO", "Untreated"), .qcMethod = "variance")
#' 
#' # Assessing the signal of empty wells across plates
#' qc(consolidatedData, .saveto = "path/to/folder/", .qcMethod = "emptywells")
#' 
#' # Assessing the z'-factor for individual plates by looking at the distribution between the positive and negative controls
#' qc(consolidatedData, .saveto = "path/to/folder/", .ctrls = c("BzCl", "DMSO"), .qcMethod = "zprime")
#' 
#' # Running multiple QC assessments
#' qc(consolidatedData, .saveto = "path/to/folder/", .ctrls = c("BzCl", "DMSO", "Untreated"), .qcMethod = c("variance", "emptywells", "firstcolumn"))
#' 
#' # Running all QC assessments
#' # Note: If zprime is being used, only the first two controls will be used as the positive and negative control. Set the order of controls in that case accordingly.
#' qc(consolidatedData, .saveto = "path/to/folder/", .ctrls = c("BzCl", "DMSO", "Untreated"), .qcMethod = "all")
#' 
#' # Assessing control plates
#' data <- readFiles(.readfrom, .fileformat = c(".csv", ".txt"), .format = "EnVision")
#' plot(data, .saveto = "path/to/folder/")
#' }}
#'
#' @keywords drug screen consolidate raw measurement dispensing
#' 
#' @importFrom ggplot2 ggplot aes facet_grid geom_boxplot geom_point ggtitle scale_fill_continuous scale_y_continuous
#' @importFrom ggrepel geom_text_repel
#' @importFrom gridExtra grid.arrange
#' @importFrom utils write.csv2
#' @importFrom stats IQR quantile
#' 
#' @export

qc <- function(consolidatedData, .saveto, .ctrls,  .qcMethod){
 
  
  # Check, if the data has been provided as an object of class S3:consolidatedData
  if(missing(consolidatedData)){stop("Data missing! Please provide a consolidated data set.", call. = TRUE)}
  if(class(consolidatedData) != "consolidatedData"){
    stop("Provided data not of class 'consolidatedData'!", call. = TRUE)
  } else {
    # Assign object to internal variable
    ctrlData <- consolidatedData[["consolidated"]]
  }
  
  # Check, if a folder location has been provided
  if(missing(.saveto)){ message("Folder location not specified. Saving to default working directory:\n", getwd()) }
  if(missing(.saveto)){ .saveto <- getwd() }
  # Create folder, if it does not exist
  if(!file.exists(file.path(.saveto))){ dir.create(file.path(.saveto, "quality"), showWarnings = FALSE, recursive = TRUE) }
  if(!utils::file_test("-d", .saveto)){stop("in '.saveto'. Argument needs to be a valid folder location.", call. = TRUE)}
  

  # Check, if a set of controls have been selected, otherwise retrieve controls from the consolidated dispensing data
  if(missing(.ctrls)){
    if(is.null(consolidatedData[["dispensingData"]][["dataList"]][["ctrlList"]])){
      stop("in '.ctrls'. Please select a set of controls to analyse.", call. = TRUE)
    }else{
      .controls <- c(consolidatedData[["dispensingData"]][["dataList"]][["ctrlList"]]$NAME,
                    consolidatedData[["dispensingData"]][["origData"]][[".addUntreated"]]$name)
    }
  }
  else if(any(!.ctrls %in% ctrlData$Drug)){ stop("in '.ctrls'. One or more controls missing in th data set. Only the controls ", gsub(",([^,]*)$"," and\\1", paste(.ctrls[.ctrls %in% ctrlData$Drug], collapse = ", ")), " were found.", call. = TRUE) }
  else { .controls <- .ctrls[.ctrls %in% ctrlData$Drug] }
  
  
  # Check, if a qc method has been provided
  if(any(missing(.qcMethod), all(!.qcMethod %in% c("variance", "emptywells", "firstcolumn", "zprime", "all")))){
    stop("Quality control method missing! Please provide a valid QC method: 'variance', 'emptywells', 'firstcolumn', 'zprime' or 'all'.", call. = TRUE)
  }else if("all" %in% .qcMethod){ 
    .qcMethod <- c("variance", "emptywells", "firstcolumn", "zprime") 
  }else{ 
    .qcMethod <- .qcMethod[.qcMethod %in% c("variance", "emptywells", "firstcolumn", "zprime")] 
  }
  
  
  
  # Remove treatments, in which one of the controls has been added to a single drug treatment
  ctrlData <- do.call(rbind, setNames(lapply(split(ctrlData, ctrlData$Sample), function(x) x[!(duplicated(x["Combination.ID"]) | duplicated(x["Combination.ID"], fromLast = TRUE)), ]), NULL))

  # Get the data for only the control treatments
  ctrlData <- subset(ctrlData, Drug %in% .controls)
  
  # Set factors for individual columns to allow proper sorting
  ctrlData$Plate.Number <- factor(as.numeric(ctrlData$Plate.Number), levels = order(unique(ctrlData$Plate.Number)))
  ctrlData <- ctrlData[with(ctrlData, order(Dispensing.Set, Plate.Number, Drug)),]
  
  
  ctrlDFS <- lapply(split(ctrlData, ctrlData$Sample), list)
  ctrlDFS <- lapply(ctrlDFS, setNames, "data")
  

  
  # Function for analyzing variance of controls
  .variance <- function(data){
    
    ctrlDFS = list()
    
    # Run analysis for each sample
    for(samplename in unique(ctrlData$Sample)){
      
      cat(" > Running QC for ", samplename, ".", sep = "")
      
      
      # Subset data set by sample
      ctrlDFS[[samplename]][["variance"]][["data"]] <- data[[samplename]][["data"]]
      
      
      # Creating a function to mark outliers
      is_outlier <- function(x) {
        return(x < stats::quantile(x, 0.45) - 1 * stats::IQR(x) | x > stats::quantile(x, 0.55) + 1 * stats::IQR(x))
      }
      
      # Set the limit for the y-axis
      yScaleLimit <- max(ctrlDFS[[samplename]][["variance"]][["data"]]$CPS)
      yScaleLimit <- ceiling(yScaleLimit/10^(nchar(yScaleLimit)-1))*10^(nchar(yScaleLimit)-1)
      
      
      # and split data set into individual control treatments
      ctrlDFS[[samplename]][["variance"]][["data"]] <- split(ctrlDFS[[samplename]][["variance"]][["data"]], ctrlDFS[[samplename]][["variance"]][["data"]]$Drug)
      
      
      # Create folder structure for data and plots to be saved to
      if(!file.exists(file.path(.saveto, "quality", samplename, "variance"))){ dir.create(file.path(.saveto, "quality", samplename, "variance"), showWarnings = FALSE, recursive = TRUE) }
      if(!file.exists(file.path(.saveto, "quality", samplename, "data"))){ dir.create(file.path(.saveto, "quality", samplename, "data"), showWarnings = FALSE, recursive = TRUE) }
      
      
      # Save raw measurements for each of the controls to a .csv file
      for(ctrl in .controls){
        utils::write.csv2(ctrlDFS[[samplename]][["variance"]][["data"]][[ctrl]], file = file.path(.saveto, "quality", samplename, "data", paste(ctrl, "-rawCPS.csv", sep ="")), row.names = FALSE, quote = FALSE)
      }
      
      
      
      # Identifying outliers and labeling them with the corresponding well number
      ctrlDFS[[samplename]][["variance"]][["data"]] = lapply(ctrlDFS[[samplename]][["variance"]][["data"]], function(x){ x$.outlier = with(x, ifelse(is_outlier(CPS), as.character(Destination.Well), NA)); return(x)} )
      
      
      
      
      # -----------------------------------------------------------------------------------------------------------------------------------
      # Plotting variance for selected controls as boxplot
      for(ctrl in .controls){
        cat('\r', " > Plotting ", ctrl, " variance for ", samplename, ".", strrep(" ", 100), sep = "")
        
        ctrlDFS[[samplename]][["variance"]][["boxplot"]][[ctrl]] <- 
          ggplot2::ggplot(ctrlDFS[[samplename]][["variance"]][["data"]][[ctrl]], aes(x = Sample, y = CPS)) +
          geom_point(aes(colour = Sample)) +
          # geom_text(aes(label = .outlier), na.rm = TRUE, hjust = -0.4) +
          geom_boxplot(alpha = 0.2, notchwidth = 2) +
          scale_y_continuous(name = "CPS", limits = c(0, yScaleLimit)) +
          # facet_grid(Sample~.) +
          xlab(NULL) +
          ylab("CPS") +
          ggtitle(paste(ctrl, sep = "")) +
          theme(plot.title = element_text(hjust = 0.5), legend.position = "none")
        
        # ggsave(filename = file.path(.saveto, "quality", samplename, "controls", paste(ctrl, "variance.png", sep = "-")), device = "png", width = 224, height = 127, units = "mm", dpi = 300, limitsize = FALSE)
        
      }
      
      
      g <- gridExtra::grid.arrange(grobs = ctrlDFS[[samplename]][["variance"]][["boxplot"]], nrow = 1)
      ggsave(filename = file.path(.saveto, "quality", samplename, "variance", paste(gsub(",([^,]*)$"," &\\1", paste(names(ctrlDFS[[samplename]][["variance"]][["boxplot"]]), collapse = ", ")), "variance.png", sep = " ")), g, device = "png", width = 224, height = 127, units = "mm", dpi = 300, limitsize = FALSE)
      

      
      
      # -----------------------------------------------------------------------------------------------------------------------------------
      # Plotting variance for selected controls as boxplot by plate
      for(ctrl in .controls){
        cat('\r', " > Plotting \'", ctrl, "\' variance by plate for ", samplename, ".", sep = "")
        
        ctrlDFS[[samplename]][["variance"]][["boxplot-byplate"]][[ctrl]] <- 
          ggplot2::ggplot(ctrlDFS[[samplename]][["variance"]][["data"]][[ctrl]], aes(x = Plate.Number, y = CPS)) +
          geom_point(aes(colour = Plate.Number)) +
          # geom_text(aes(label = .outlier), na.rm = TRUE, hjust = -0.4) +
          geom_text_repel(aes(label = .outlier), na.rm = TRUE, force = 1, direction = "both", point.padding = unit(0.1, "lines"), box.padding = unit(0.1, "lines"),
                          segment.size = 0.005, segment.color = "black", nudge_x = 0.2, nudge_y = 0, hjust = 0, size = 1.5) +
          geom_boxplot(alpha = 0.2, notchwidth = 2) +
          scale_y_continuous(name = "CPS", limits = c(0, yScaleLimit)) +
          facet_grid(Sample~.) +
          xlab(NULL) +
          ylab("CPS") +
          ggtitle(paste(ctrl, sep = "")) +
          labs(tag = paste("Plate:", sep = " ")) +
          theme(plot.title = element_text(hjust = 0.5), legend.position = "none",
                text = element_text(size = 10), axis.text.x = element_text(angle = 0, vjust = 1, hjust=0.5),
                plot.tag = element_text(hjust = 1, size = 8), plot.tag.position = c(0.06, 0.015))
        
        ggsave(filename = file.path(.saveto, "quality", samplename, "variance", paste(ctrl, "byPlate.png", sep = "-")), device = "png", width = 224, height = 127, units = "mm", dpi = 300, limitsize = FALSE)
        
      }
      

      # Remove outlier annotation from the composite plots
      ctrlDFS[[samplename]][["variance"]][["boxplot-byplate"]] <- lapply(ctrlDFS[[samplename]][["variance"]][["boxplot-byplate"]], function(x){ x$layers[[2]] <- NULL; return(x) })
      
      g <- gridExtra::grid.arrange(grobs = ctrlDFS[[samplename]][["variance"]][["boxplot-byplate"]], nrow = 3)
      ggsave(filename = file.path(.saveto, "quality", samplename, "variance", paste(gsub(",([^,]*)$"," &\\1", paste(names(ctrlDFS[[samplename]][["variance"]][["boxplot"]]), collapse = ", ")), "-variance byPlate.png", sep = " ")), g, device = "png", width = 224, height = 64*3, units = "mm", dpi = 300, limitsize = FALSE)
      
      
      cat('\r', "Finished plotting ", gsub(",([^,]*)$"," &\\1", paste(names(ctrlDFS[[samplename]][["variance"]][["boxplot"]]), collapse = ", ")), " variance for ", samplename, ".", strrep(" ", 100), '\n', sep = "")
      

    }
    
    return(ctrlDFS)
    
  }
  
  
  
  .zprime <- function(data){
    
    ctrlDFS = list()
    
    # Run analysis for each sample
    for(samplename in unique(ctrlData$Sample)){
      
      cat(" > Calculating z'-factor score for ", samplename, ".", sep = "")
      
      
      # Subset and split data set by positive and negative control
      ctrlDFS[[samplename]][["z-factor"]][["data"]] <- subset(data[[samplename]][["data"]], Drug %in% .controls[1:2])
      
      
      # Create a data set for the calculation of the z'-factor
      .zdata <- ctrlDFS[[samplename]][["z-factor"]][["data"]]

      
      # Function for calculating z'-factor based on
      # Zhang, J. H., Chung, T. D., Oldenburg, K. R. (1999): A Simple Statistical Parameter for Use in Evaluation and Validation 
      # of High Throughput Screening Assays. J Biomol Screening 4 (2), S. 67–73. DOI: 10.1177/108705719900400206
      zprime <- function(a, b){
        1-3*(sd(a)+sd(b))/abs(mean(a)-mean(b))
      }

      # Call zprime() function and run calculations
      .zdata <- do.call(rbind, setNames(lapply(split(.zdata, .zdata$Plate.Number), function(x){ x[[paste(".", "score", sep = "")]] <- zprime(x$CPS[which(x$Drug == .controls[1])], x$CPS[which(x$Drug == .controls[2])]); return(x) }), NULL))
      
      # Select specific rows and columns
      .zdata <- .zdata[!grepl("Dispensing.Set|Combination.ID|Drug|CAS.number|Drug.Concentration|Unit|Transfer.Volume|Source.Plate.Barcode|Source.Well|Destination.Well|CPS|.outlier", names(.zdata), ignore.case = TRUE)]
      .zdata <- .zdata[!duplicated(.zdata),]
      
      

      # Create folder structure for data and plots to be saved to
      if(!file.exists(file.path(.saveto, "quality", samplename, "z-factor"))){ dir.create(file.path(.saveto, "quality", samplename, "z-factor"), showWarnings = FALSE, recursive = TRUE) }
      
      
      # Plot the z-factor scores for each plate    
      cat('\r', " > Plotting z'-factor score for ", samplename, ".", strrep(" ", 100), sep = "")
      
      
      ctrlDFS[[samplename]][["z-factor"]][["plot"]] <- 
        ggplot2::ggplot(.zdata, aes(x = Plate.Number, y = .zdata[,grepl("score", names(.zdata))])) +
        geom_point(aes(colour = .zdata[,grepl("score", names(.zdata))]), size = 7, alpha = .zdata[,grepl("score", names(.zdata))]) +
        scale_y_continuous(name = "z-factor", limits = c(0, 1)) +
        geom_hline(aes(yintercept=0.5), lwd=0.2, colour="red") +
        facet_grid(Sample ~ .) +
        xlab(NULL) +
        ylab("z-factor") +
        labs(tag = paste("Plate:", sep = " ")) +
        ggtitle(paste("Z'factor score (", samplename, ")", sep = "")) +
        theme(plot.title = element_text(hjust = 0.5), legend.position = "none",
              text = element_text(size = 26),
              axis.text.x = element_text(angle = 0, vjust = 1, hjust=1),
              plot.tag = element_text(hjust = 1, size = 21), plot.tag.position = c(0.025, 0.01))
      
      ggsave(filename = file.path(.saveto, "quality", samplename, "z-factor", paste(samplename, " z-factor.png", sep = "")), device = "png", width = 840, height = 567, units = "mm", dpi = 300, limitsize = FALSE)
      
      
      cat('\r', "Finished plotting z'factor for ", samplename, ".", strrep(" ", 100), '\n', sep = "")
      
      ctrlDFS[[samplename]][["z-factor"]][["data"]] <- .zdata
      
    }
    
    return(ctrlDFS)
    
  }
  
  
  
  .emptywells <- function(data){
    
    ctrlDFS = list()
    
    # Run analysis for each sample
    for(samplename in unique(ctrlData$Sample)){
      
      cat(" > Analyzing empty & excluded wells for ", samplename, ".", sep = "")
      
      
      # Create folder structure for data and plots to be saved to
      if(!file.exists(file.path(.saveto, "quality", samplename, "empty-wells", "plates"))){ dir.create(file.path(.saveto, "quality", samplename, "empty-wells", "plates"), showWarnings = FALSE, recursive = TRUE) }
      
      
      for(platenumber in unique(ctrlData$Plate.Number)){
        
        # Subset data by sample and plate number
        .exwells <- subset(consolidatedData[["measurements"]], Sample == samplename & Plate.Number == platenumber)
        # Subset to only the excluded wells
        # .exwells <- subset(.exwells, Well %in% consolidatedData[["dispensingData"]][["origData"]]$listofExWells)
        # Subset to excluded and empty wells, such as on the last plate
        .exwells <- subset(.exwells, !Well %in% unique(as.character(subset(consolidatedData[["consolidated"]], Sample == samplename & Plate.Number == platenumber)$Destination.Well)))

        # Separate wells into rows and columns
        .exwells$Row <- as.character(with(.exwells, gsub("\\d", "", Well)))
        .exwells$Column <- as.numeric(with(.exwells, gsub("\\D", "", Well)))
        .exwells$Row <- factor(.exwells$Row, levels = with(.exwells, unique(Row[order(nchar(Row), Row)])))
        
        
        cat('\r', " > Plotting empty & excluded wells on plate ", platenumber, " for ", samplename, ".", strrep(" ", 100), sep = "")
        
        ctrlDFS[[samplename]][["empty-wells"]][[paste("plate", platenumber, sep = "-")]] <- 
          ggplot2::ggplot(.exwells, aes(x = factor(Column), y = factor(Row, rev(levels(Row))))) + 
            geom_raster(aes(fill = CPS), na.rm=TRUE) +
            scale_fill_continuous(low = "#FFCCFF", high = "#FF0033") +
            geom_hline(yintercept=seq(1.5, length(unique(.exwells$Row))-0.5, 1), lwd=0.2, colour="black") +
            geom_vline(xintercept=seq(1.5, length(unique(.exwells$Column))-0.5, 1), lwd=0.2, colour="black") +
            scale_x_discrete(position = "top") +
            facet_grid(Sample~.) +
            coord_fixed(ratio = 1) +
            xlab(NULL) +
            ylab(NULL) +
            labs(caption = paste("Dispensing Plate:", .exwells$Destination.Plate.Barcode, sep = " ")) +
            theme(plot.title = element_text(hjust = 0.5),
                  plot.caption = element_text(size = 12),
                  text = element_text(size = 26),
                  axis.line = element_line(color='black'),
                  panel.grid.major = element_line(color='lightgrey', linetype = "dotted"),
                  panel.grid.minor = element_line(color='lightgrey', linetype = "dotted"),
                  panel.background = element_rect(fill=alpha('white', 1)),
                  panel.border = element_rect(colour = "black", fill=NA, size=1),
                  legend.position = "none")
        
        
        ggsave(filename = file.path(.saveto, "quality", samplename, "empty-wells", "plates", paste(samplename, " Plate-", platenumber, " empty-wells.png", sep = "")), device = "png", width = 840, height = 840/(3/2), units = "mm", dpi = 300, limitsize = FALSE)
        
      
      }
      
      cat('\r', "Finished plotting empty & excluded wells for ", samplename, ".", strrep(" ", 100), '\n', sep = "")
      
      
    }
  
    return(ctrlDFS)
    
  }
  
  
  
  .firstcolumn <- function(data){
    
    ctrlDFS = list()
    
    # Run analysis for each sample
    for(samplename in unique(ctrlData$Sample)){
      
      cat(" > Analyzing wells of the first column for ", samplename, ".", sep = "")
      

      # Create folder structure for data and plots to be saved to
      if(!file.exists(file.path(.saveto, "quality", samplename, "empty-wells", "plates"))){ dir.create(file.path(.saveto, "quality", samplename, "empty-wells", "plates"), showWarnings = FALSE, recursive = TRUE) }
      

      # Subset data by sample and plate number
      .exwells <- subset(consolidatedData[["measurements"]], Sample == samplename)
      # Subset to only the excluded wells
      .exwells <- subset(.exwells, !Well %in% consolidatedData[["dispensingData"]][["origData"]]$listofExWells)

      # Select only wells found in the first column
      .exwells <- .exwells[which(as.numeric(gsub("\\D", "", .exwells$Well)) == 1),]
      
      # Separate wells into rows and columns
      .exwells$Row <- as.character(with(.exwells, gsub("\\d", "", Well)))
      .exwells$Column <- as.numeric(with(.exwells, gsub("\\D", "", Well)))
      .exwells$Row <- factor(.exwells$Row, levels = with(.exwells, unique(Row[order(nchar(Row), Row)])))
      
      # Set factors for individual columns to allow proper sorting
      .exwells$Plate.Number <- factor(as.numeric(.exwells$Plate.Number), levels = order(unique(.exwells$Plate.Number)))
      
      
      cat('\r', " > Plotting wells of the first column for ", samplename, ".", strrep(" ", 100), sep = "")
      
      ctrlDFS[[samplename]][["first-column"]] <- 
        ggplot2::ggplot(.exwells, aes(x = Plate.Number, y = factor(Row, rev(levels(Row))))) + 
          geom_raster(aes(fill = CPS), na.rm=TRUE) +
          scale_fill_continuous(low = "#FFCCFF", high = "#FF0033") +
          # scale_fill_discrete(na.value = NA) +
          # scale_fill_hue("Drugs", l=80, h = c(270, 360), na.value = NA) +
          # geom_text(aes(label = paste(ifelse(is.na(ctrlData$Drug),"",ctrlData$Drug), sep = "\n")), color ="black", hjust = 0.5, size = 2.7) +
          geom_hline(yintercept=seq(1.5, length(unique(.exwells$Row))-0.5, 1), lwd=0.2, colour="black") +
          geom_vline(xintercept=seq(1.5, length(unique(.exwells$Plate.Number))-0.5, 1), lwd=0.2, colour="black") +
          scale_x_discrete(position = "bottom") +
          facet_grid(Sample~.) +
          coord_fixed(ratio = 2/3) +
          xlab(NULL) +
          ylab(NULL) +
          labs(tag = paste("Plate:", sep = " ")) +
          theme(plot.title = element_text(hjust = 0.5),
                plot.caption = element_text(size = 12),
                text = element_text(size = 26),
                axis.line = element_line(color='black'),
                panel.grid.major = element_line(color='lightgrey', linetype = "dotted"),
                panel.grid.minor = element_line(color='lightgrey', linetype = "dotted"),
                panel.background = element_rect(fill=alpha('white', 1)),
                panel.border = element_rect(colour = "black", fill=NA, size=1),
                legend.position = "none",
                plot.tag = element_text(hjust = 1, size = 21), 
                plot.tag.position = c(0.015, 0.01))
      
      ggsave(filename = file.path(.saveto, "quality", samplename, "empty-wells", paste(samplename, " empty-wells (first column).png", sep = "")), device = "png", width = 840, height = 840/(3/2), units = "mm", dpi = 300, limitsize = FALSE)
     
      
      cat('\r', "Finished plotting wells of the first column for ", samplename, ".", strrep(" ", 100), '\n', sep = "")
      
    }
    
    return(ctrlDFS)
            
  }
  
  
  
  # Run all selected analysis methods
  for(.method in .qcMethod){
    
    # Select the function to run the analysis with and attach list of the new analysis to existing list
    # list from previously run qc analysis will be replaced 
    ctrlDFS <- utils::modifyList(ctrlDFS, switch(.method, variance = .variance(ctrlDFS), emptywells = .emptywells(ctrlDFS), firstcolumn = .firstcolumn(ctrlDFS), zprime = .zprime(ctrlDFS)))
    # alternatively, keep and attach each new list to existing list
    # ctrlDFS <- mapply(c, ctrlDFS, switch(.method, variance = .variance(ctrlDFS), emptywells = .emptywells(ctrlDFS), firstcolumn = .firstcolumn(ctrlDFS), zprime = .zprime(ctrlDFS)), SIMPLIFY=FALSE)
    
  }
  
  class(ctrlDFS) <- "controlData"
  
  return(ctrlDFS)
  
}