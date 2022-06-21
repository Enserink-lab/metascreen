#' Run dose-response model
#'
#' @description
#' An essential component of the modular library \pkg{metascreen}, and imperative for the analysis of the experimental data from a drug sensitivity screen.
#' \emph{\code{runDRM}} is a function that runs the dose-response model for the analysis of single drug treatments. It performs curve fitting using a \emph{four-parameter log-logistic function (LL.4)}, 
#' and estimates the EC10, EC50 and EC90. The function also plots the single drug curves by viability and inhibition for each single drug.
#' 
#' The dose-response analysis is performed with the function \code{\link[drc]{drm}} from the R-package \emph{drc} as described in the
#' original paper by Ritz, 2015, PLOS ONE (see references).
#'  
#' @param processedData an object of class 'processedData'.
#' @param .saveto string; path to a folder location where the results are saved to.
#' @param .plot logical; if TRUE, then fitted dose-response will be plotted. The default is FALSE, in which no plots will be generated.
#' 
#' 
#' @details The function \code{runDRM} is used to run the first of a series of analyses, specifically, it looks at the dose response of each single drug. 
#' It performs curve fitting based on a four-parameter log-logistic model and estimates the EC50s, along with the EC10 and EC90. The curve fitting is performed with lower limits set to 0 and 
#' the upper limit to 1. In cases where the curve fitting fails, such as due to lack of response resulting in a flat line, or due to outliers, or any other irregularities, the fitting is performed without restrictions to either one of the limits.
#' 
#' All estimates of the EC50s is exported as a table in a csv file format. Plots are exported as composite of all drugs for each sample as individual png image files.
#' 
#' Files are saved either to the specified location or the default working environment, with the corresponding folder structure: 'results/EC50s' and 'results/graphs/single drug response'
#' 
#' 
#' @return Returns an object of type list with the dose-response model, a summary and the estimated ECs.
#' 
#' @references Ritz, C., Baty, F., Streibig, J. C., Gerhard, D. (2015) Dose-Response Analysis Using R. PLOS ONE 10 (12), e0146021. DOI: 10.1371/journal.pone.0146021
#' 
#' @seealso \code{\link[drc]{drm}}
#'
#' @examples
#' \donttest{\dontrun{
#' # Process consolidated data
#' runDRM(processedData, .saveto = "path/to/folder/")
#' }}
#'
#' @keywords drug screen analysis single curve fitting log logistic drc drm
#' 
#' @importFrom utils capture.output packageVersion
#' @importFrom drc drm LL.4 drmc ED
#' @importFrom grDevices dev.off png
#' @importFrom graphics grid legend mtext par plot points title
#' 
#' @export

runDRM <- function(processedData, .saveto, .plot){
  
  # Check, if the data has been provided as an object of class S3:processedData
  if(missing(processedData)){stop("Data missing! Please provide a normalized data set.", call. = TRUE)}
  if(class(processedData) != "processedData"){
    stop("Provided data not of class 'processedData'!", call. = TRUE)
  } else {
    # Extract all the necessary data from the class S3 object
    analysisData <- processedData[["analysisData"]]
    synDataset <- processedData[["splitDataset"]]
    synDRM <- processedData[["doserespMatrix"]]
  }
  
  
  # Check, if a folder location has been provided
  if(missing(.saveto)){ .saveto <- getwd() }
  # Create folder, if it does not exist
  if(!file.exists(file.path(.saveto))){ dir.create(file.path(.saveto), showWarnings = FALSE, recursive = TRUE) }
  if(!utils::file_test("-d", .saveto)){stop("in '.saveto'. Argument needs to be a valid folder location.", call. = TRUE)}
  
  
  # Check, if dose response models should be plotted
  if(missing(.plot)){ 
    .plot <- FALSE 
  } else if(!is.logical(.plot)){
    stop("in '.plot'. Argument needs to be a logical value, either TRUE or FALSE.", call. = TRUE)
  }
  
  
  
  # SECTION A: DOSE RESPONSE MODEL AND CURVE FITTING ##########################################################
  
  cat('\n', "Running dose-response (drm) analysis.", strrep(" ", 100), '\n\n', sep = "")
  
  
  efs <- list()
  
  for(samplename in names(synDataset)){
    
    for(drugname in sort(unique(synDataset[[samplename]][["singleDrugResponseData"]]$Drug))){
      
      nocalculations <- length(synDataset) * length(unique(synDataset[[samplename]][["singleDrugResponseData"]]$Drug))
      noutput <- 2
      
      if(exists("p")){ p <- p+noutput }else{ p <- 1 }
      
      cat('\r', "[", p, "/", nocalculations*noutput, "] ", "Fitting curve for ", samplename, " and ", drugname, ".", strrep(" ", 100), sep = "")
      

      # Run dose response model and curve fitting with R-package drc based on
      # Ritz, C., Baty, F., Streibig, J. C., Gerhard, D. (2015) Dose-Response Analysis Using R. PLOS ONE 10 (12), e0146021. DOI: 10.1371/journal.pone.0146021
      # package version: drc v3.0-1
      
      # Note: viability used for viability plotting 
      # and for summary statistics (ED50, Std. Error , p-value)
      utils::capture.output(type="message",
      efs[[samplename]][[drugname]][["drm"]] <- tryCatch({drc::drm(Viability ~ Drug.Concentration, curveid = Drug, data = synDataset[[samplename]][["singleDrugResponseData"]],
                                                              subset = Drug == drugname, fct = drc::LL.4(fixed = c(NA,NA,NA,NA), names = c("Slope", "Lower Limit", "Upper Limit", "ED50")),
                                                              lowerl = c(-Inf, 0, 0, 0), upperl = c(Inf, 1, 1, Inf),
                                                              separate = TRUE, control = drc::drmc(useD = FALSE))}, error=function(e){
                                                                # If model fails with lower and upper limits, same model is run again without limits 
                                                                # and with derivatives for estimation, this usually works for most of the problematic dose response curves
                                                                drc::drm(Viability ~ Drug.Concentration, curveid = Drug, data = synDataset[[samplename]][["singleDrugResponseData"]],
                                                                    subset = Drug == drugname, fct = drc::LL.4(fixed = c(NA,NA,NA,NA), names = c("Slope", "Lower Limit", "Upper Limit", "ED50")),
                                                                    separate = TRUE, control = drc::drmc(useD = TRUE))
                                                              }) )
      
      
      # Calculate summary and EC50 based on the curve fitting data
      cat('\r', "[", p+1, "/", nocalculations*noutput, "] ", "Estimating the ED50 for ", samplename, " and ", drugname, ".", strrep(" ", 100), sep = "")
      
      efs[[samplename]][[drugname]][["summary"]] <- suppressWarnings(summary(efs[[samplename]][[drugname]][["drm"]]))
      efs[[samplename]][[drugname]][["ED"]] <- suppressWarnings(drc::ED(efs[[samplename]][[drugname]][["drm"]], c(10, 50, 90), interval = "delta", type = "relative", bound = FALSE, display = FALSE))
      
      if(p == nocalculations*noutput){cat('\r', "Finished fitting dose response model for all samples.", sep = "")}
      
    }
    
  }
  
  
  
  # SECTION B: EC50 CONSOLIDATION #############################################################################
  # Extract EC50s and merge into a single data set for each sample
  
  EC50DSList <- list()
  
  for(samplename in names(efs)){
    
    for(drugname in names(efs[[samplename]])){
      
      cat('\r', " > Extracting ECs for ", samplename, ":", drugname, ".", strrep(" ", 100), sep = "")
      
      
      # Extract the ED50 value from efs for each cell line and drug
      EC50DSList[[samplename]][[drugname]] <- as.data.frame(efs[[samplename]][[drugname]][["ED"]], stringsAsFactors = FALSE)
      # Add columns with cell line and drug name
      EC50DSList[[samplename]][[drugname]]$Sample <- samplename
      EC50DSList[[samplename]][[drugname]]$Drug <- drugname
      
      # Convert row names to columns
      EC50DSList[[samplename]][[drugname]] <- cbind(ED = rownames(EC50DSList[[samplename]][[drugname]]), data.frame(EC50DSList[[samplename]][[drugname]], row.names=NULL, check.names = FALSE), stringsAsFactors = FALSE)
      # Clean-up of ED labels
      EC50DSList[[samplename]][[drugname]]$ED <- with(EC50DSList[[samplename]][[drugname]], sapply(strsplit(ED, ":"), `[`, 3))
      # Select the order of columns
      EC50DSList[[samplename]][[drugname]] <- EC50DSList[[samplename]][[drugname]][c("Sample", "Drug", "ED", "Estimate", "Std. Error")]

    }
    
    
    EC50DSList[[samplename]] <- do.call(rbind, setNames(EC50DSList[[samplename]], NULL))
    
    # Create subfolder if it does not exist
    if(!file.exists(file.path(.saveto, "results", "EC50s"))){ dir.create(file.path(.saveto, "results", "EC50s"), showWarnings = FALSE, recursive = TRUE) }

    # Export the data set prior EC50 calculation
    write.csv2(EC50DSList[[samplename]], file = file.path(.saveto, "results", "EC50s", paste(samplename, "EC50s.csv", sep = "_")), row.names = FALSE, quote = FALSE)
    
    cat('\r', "Finished estimating ECs for ", samplename, ".", strrep(" ", 100), '\n', sep = "")
    
    
    if(samplename == utils::tail(names(efs), n = 1L)){
      message('\n', "Estimated ECs saved to: ", '\n', file.path(.saveto, "results", "EC50s"), strrep(" ", 100), sep = "")
      rm(samplename, drugname)
    }
  
  }
  
  
  if(isTRUE(.plot)){
    # SECTION C: PLOTTING DOSE RESPONSE CURVES ##################################################################
    
    # Plotting Dose Response curves based on viability.
    # Plot dose response curves as a composite.
    
    for(samplename in names(efs)){
      
      cat('\r', "Plotting single drug viability curves for ", samplename, ".", strrep(" ", 100), sep = "")
      
      # Create subfolder if it does not exist
      if(!file.exists(file.path(.saveto, "results", "graphs/single drug response/by viability"))){ dir.create(file.path(.saveto, "results", "graphs/single drug response", "by viability"), showWarnings = FALSE, recursive = TRUE) }
      
      
      grDevices::png(filename = file.path(.saveto, "results", "graphs/single drug response", "by viability", paste(samplename, "dose-response (viability).png", sep = " ")), 
                     width = 7920, height = 4080, units = "px", pointsize = 24)
      op <- graphics::par(mfrow = c(7, 10), oma = c(4,0,7,0))
      
      for(drugname in sort(names(efs[[samplename]]))){
        
        cat('\r', " > Plotting dose-response curve: ", samplename, ": ", drugname, strrep(" ", 100), sep = "")
        
        EC50 <- efs[[samplename]][[drugname]][["ED"]][paste("e", drugname, "50", sep = ":"), "Estimate"]
        
        # plot all points (replicates) by type = all, or average of points by type = average
        graphics::plot(efs[[samplename]][[drugname]][["drm"]], level = drugname, log = "x", type = "all", main = drugname, xlab = "Concentration (\u00b5M)", ylab = "Viability",
                       pch=1, cex=1, lty=1, lwd=2, ylim=c(0,1), yaxs="i", xaxs="r", bty="n")
        
        graphics::points(EC50, efs[[samplename]][[drugname]][["summary"]][["coefficients"]]["Upper Limit:(Intercept)", "Estimate"]-((
          efs[[samplename]][[drugname]][["summary"]][["coefficients"]]["Upper Limit:(Intercept)", "Estimate"]-
            efs[[samplename]][[drugname]][["summary"]][["coefficients"]]["Lower Limit:(Intercept)", "Estimate"])/2),
          type="p", pch=16, cex = 1.0, col="red")
        
        graphics::grid(10,10, lwd = 2)
        graphics::legend("topright", ifelse(EC50 > max(efs[[samplename]][[drugname]][["drm"]][["dataList"]][["dose"]]) | EC50 < min(efs[[samplename]][[drugname]][["drm"]][["dataList"]][["dose"]]) 
                                            # Not printing the EC50, if outside the dose range, or if response is linear
                                            | efs[[samplename]][[drugname]][["drm"]][["coefficients"]]["Lower Limit:(Intercept)"] == efs[[samplename]][[drugname]][["drm"]][["coefficients"]]["Upper Limit:(Intercept)"],
                                            "", paste("EC50:", format(EC50, digits = 3), "\u00b5M", sep = " ")), cex = 1.0, col = "black", bty = "n")
        
      }
      
      cat('\r', "Saving dose-response plots for ", samplename, ".", strrep(" ", 100), sep = "")
      
      graphics::title(paste(samplename, ": Dose Response (Viability)", sep = ""), outer = TRUE, cex.main = 4)
      graphics::mtext(paste("  Package: ", "drc", " (v", utils::packageVersion("drc"), ")   ",  sep = ""), side = 1, line = 2, adj = 1, cex = 1, col = "black", outer = TRUE) 
      grDevices::dev.off()
      
      cat('\r', "Finished plotting dose-response curves for ", samplename, ".", strrep(" ", 100), '\n', sep = "")
      
      
      if(samplename == utils::tail(names(efs), n = 1L)){
        message('\n', "Dose-response plots saved to: ", '\n', file.path(.saveto, "results", "graphs/single drug response/by viability"), strrep(" ", 100), sep = "")
        rm(samplename, drugname, EC50, op)
      }
      
    }
    
    
    
    # Plotting Dose Response curves based on inhibition.
    # Plot dose response curves as a composite.
    
    for (samplename in names(synDataset)){
      
      cat('\r', "Plotting single drug viability curves for ", samplename, ".", strrep(" ", 100), sep = "")
      
      # Create subfolder if it does not exist
      if(!file.exists(file.path(.saveto, "results", "graphs/single drug response/by inhibition"))){ dir.create(file.path(.saveto, "results", "graphs/single drug response", "by inhibition"), showWarnings = FALSE, recursive = TRUE) }
      
      
      grDevices::png(filename = file.path(.saveto, "results", "graphs/single drug response/by inhibition", paste(samplename, "dose-response (inhibition).png", sep = " ")),
                     width = 7920, height = 4080, units = "px", pointsize = 24)
      op <- graphics::par(mfrow = c(7, 10), oma = c(4,0,7,0))
      
      for(drugname in sort(unique(synDataset[[samplename]][["singleDrugResponseData"]]$Drug))){
        
        cat('\r', " > Plotting dose-response curve: ", samplename, ": ", drugname, strrep(" ", 100), sep = "")
        
        .drm.data <- subset(synDataset[[samplename]][["singleDrugResponseData"]], Drug == drugname)
        
        
        # Run dose response model and curve fitting with R-package drc based on
        # Ritz, C., Baty, F., Streibig, J. C., Gerhard, D. (2015) Dose-Response Analysis Using R. PLOS ONE 10 (12), e0146021. DOI: 10.1371/journal.pone.0146021
        # package version: drc v3.0-1
        
        utils::capture.output(type="message",
                              .drm.model <- tryCatch({drc::drm(Inhibition ~ Drug.Concentration, curveid = Drug, data = .drm.data,
                                                               subset = Drug == drugname, fct = drc::LL.4(fixed = c(NA,NA,NA,NA), names = c("Slope", "Lower Limit", "Upper Limit", "ED50")),
                                                               lowerl = c(-Inf, 0, 0, 0), upperl = c(Inf, 1, 1, Inf),
                                                               separate = TRUE, control = drc::drmc(useD = FALSE))}, error=function(e){
                                                                 # If model fails with lower and upper limits, same model is run again without limits 
                                                                 # and with derivatives for estimation, this usually works for most of the problematic dose response curves
                                                                 drc::drm(Inhibition ~ Drug.Concentration, curveid = Drug, data = .drm.data,
                                                                          subset = Drug == drugname, fct = drc::LL.4(fixed = c(NA,NA,NA,NA), names = c("Slope", "Lower Limit", "Upper Limit", "ED50")),
                                                                          separate = TRUE, control = drc::drmc(useD = TRUE))
                                                               }) )
        
        EC50 <- .drm.model[["coefficients"]]["ED50:(Intercept)"]
        
        graphics::par(bg="white", fg="black", col="black", col.axis="black", col.lab="black", col.main="black", col.sub="black") 
        graphics::plot(.drm.model, log = "x", main=drugname, type = "all", pch=1, cex=1, lty=1, lwd=2, 
                       ylim=c(0,1), axes = TRUE, yaxs="i", xaxs="r", bty="n", xlab="Concentration (\u00b5M)", ylab="Inhibition")
        # axis(2, at = seq(0,1,by=.2), labels = as.character(1*seq(1,0, by=-.2)), tick = TRUE)
        graphics::points(EC50, .drm.model[["coefficients"]]["Upper Limit:(Intercept)"]-((
          .drm.model[["coefficients"]]["Upper Limit:(Intercept)"]-.drm.model[["coefficients"]]["Lower Limit:(Intercept)"])/2),
          type="p", pch=16, cex = 1.0, col="red")
        graphics::grid(10,10, lwd = 2)
        graphics::legend("topleft", ifelse(EC50 > max(.drm.model[["dataList"]][["dose"]]) | EC50 < min(.drm.model[["dataList"]][["dose"]]) 
                                           # Not printing the EC50, if outside the dose range, or if response is linear
                                           | .drm.model[["coefficients"]]["Lower Limit:(Intercept)"] == .drm.model[["coefficients"]]["Upper Limit:(Intercept)"], 
                                           "", paste("EC50:", format(EC50, digits = 3), "\u00b5M", sep = " ")), cex = 1.0, col = "black", bty = "n")
        
      }
      
      cat('\r', "Saving dose-response plots for ", samplename, ".", strrep(" ", 100), sep = "")
      
      graphics::title(paste(samplename, ": Dose Response (Inhibition)", sep = ""), outer = TRUE, cex.main = 4)
      graphics::mtext(paste("  Package: ", "drc", " (v", utils::packageVersion("drc"), ")   ",  sep = ""), side = 1, line = 2, adj = 1, cex = 1, col = "black", outer = TRUE) 
      grDevices::dev.off()
      
      cat('\r', "Finished plotting dose-response curves for ", samplename, ".", strrep(" ", 100), '\n', sep = "")
      
      
      if(samplename == utils::tail(names(efs), n = 1L)){
        message('\n', "Dose-response plots saved to: ", '\n', file.path(.saveto, "results", "graphs/single drug response/by inhibition"), strrep(" ", 100), sep = "")
        rm(samplename, drugname, .drm.data, op)
      }
    }
    
  }
  
  
  # Return object of class DRM
  class(efs) <- "drm"
  return(efs)

}