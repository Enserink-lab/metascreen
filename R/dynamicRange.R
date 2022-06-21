#' Assessing the dynamic drug-activity range
#'
#' @description
#' An complementary component of the modular library \pkg{metascreen} assessing the dynamic drug-activity range for individual drug-dose responses.
#' \emph{\code{dynamicRange}} is a function that estimates the dynamic drug-activity range (DDAR) across a number of doses for each drug response and generates a set of plots.
#' 
#' @param data an object of class 'processedData' or 'drm'.
#' @param .saveto string; path to a folder location where the results are saved to.
#' 
#' 
#' @details The function \code{dynamicRange} requires either one of the two data sets, either an object of class 'processedData' or an object of class 'drm'.
#' 
#' \code{dynamicRange} is used to provide an assessment of the drug activity range across the selected range of doses. It is an essential indicator, especially for
#' drug combination screens, in which only a selected range of doses are combined to assess synergies. If doses are combined at which drugs enfold their full inhibitory potential, it won't leave enough room for potential combinatory effects. 
#' 
#' The dynamic range is considered the dose-response range between the ED10 and ED90. This function will indicate the expected and the observed drug activity range for each drug and sample. 
#' It will generate plots based on the fitted dose-response models as well as unfitted curves.
#' 
#' Files are saved either to the specified location or the default working environment, with the corresponding folder structure: 'results/Dynamic Drug Activity Range' and 'results/graphs/dynamic range'
#' 
#' @examples
#' \donttest{\dontrun{
#' # Run dose-response model and estimate the dynamic range
#' dynamicRange(processedData, .saveto = "path/to/folder/")
#' 
#' # Estimate the dynamic range based on the provided dose-response models
#' dynamicRange(doseRespModel, .saveto = "path/to/folder/")
#' }}
#'
#' @keywords drug screen analysis dose response curve matrix
#' 
#' @importFrom utils capture.output tail
#' @importFrom stats ave predict
#' @importFrom drc drm LL.4 drmc ED
#' @importFrom grDevices dev.off png rainbow rgb
#' @importFrom graphics abline axis box grid legend lines mtext par plot plot.new points rect title
#' 
#' @export

dynamicRange <- function(data, .saveto){
  
  # Check, if the data has been provided as an object of class S3:processedData
  if(missing(data)){stop("Data missing! Please provide either a data set as an object of 'processedData' or 'drm'.", call. = TRUE)}
  
  # Extract the required data from the class S3 object
  if(class(data) == "processedData"){
    analysisData <- data[["analysisData"]]
    synDataset <- data[["splitDataset"]]
    synDRM <- data[["doserespMatrix"]]
  } else if(class(data) == "drm") {
    efs <- data
  } else {
    stop("in 'data'. Argument needs to be an object of 'processedData' or 'drm'.", call. = TRUE)
  }
  
  # Check, if a folder location has been provided
  if(missing(.saveto)){ .saveto <- getwd() }
  # Create folder, if it does not exist
  if(!file.exists(file.path(.saveto))){ dir.create(file.path(.saveto), showWarnings = FALSE, recursive = TRUE) }
  if(!utils::file_test("-d", .saveto)){stop("in '.saveto'. Argument needs to be a valid folder location.", call. = TRUE)}
  
  
  
  if(class(data) == "processedData"){
    
    # SECTION A: DOSE RESPONSE MODEL AND CURVE FITTING ##########################################################
    
    cat('\n', "Running dose-response (drm) analysis.", strrep(" ", 100), '\n\n', sep = "")
    
    
    efs <- list()
    
    for(samplename in names(synDataset)){
      
      for(drugname in sort(unique(synDataset[[samplename]][["singleDrugResponseData"]]$Drug))){
        
        nocalculations <- length(synDataset) * length(unique(synDataset[[samplename]][["singleDrugResponseData"]]$Drug))
        noutput <- 2
        
        if(exists("p")){ p <- p+noutput }else{ p <- 1 }
        
        cat('\r', "[", p, "/", nocalculations*noutput, "] ", "Fitting curve for ", samplename, " and ", drugname, ".", strrep(" ", 50), sep = "")
        
        
        # Run dose response model and curve fitting
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
        cat('\r', "[", p+1, "/", nocalculations*noutput, "] ", "Estimating the ED50 for ", samplename, " and ", drugname, ".", strrep(" ", 50), sep = "")
        
        efs[[samplename]][[drugname]][["summary"]] <- suppressWarnings(summary(efs[[samplename]][[drugname]][["drm"]]))
        efs[[samplename]][[drugname]][["ED"]] <- suppressWarnings(drc::ED(efs[[samplename]][[drugname]][["drm"]], c(10, 50, 90), interval = "delta", type = "relative", bound = FALSE, display = FALSE))
        
        if(p == (nocalculations*noutput)-1){cat('\r', "Finished fitting dose response model for all samples.", strrep(" ", 100), '\n\n', sep = "")}
        
      }
      
    }
    
  }
  
  
  
  # Extracting Data from the Analysis. Exporting a list of ECs.
  # Extract the dynamic range for each drug and sample based on the EC10 and EC90.
  
  cat("Extracting Dynmaic Drug Activity Range (DDAR).", '\n', sep = "")
  
  
  ECDRList <- list()
  
  for(samplename in names(efs)){
    
    for(drugname in names(efs[[samplename]])){
      
      cat('\r', "Extracting dynmaic drug-activity range (DDAR) for ", samplename, ":", drugname, ".", sep = "")
      
      
      # Extract the all ED values for each drug and sample
      ECDRList[[samplename]][[drugname]] <- cbind(ED = rownames(efs[[samplename]][[drugname]][["ED"]]), data.frame(efs[[samplename]][[drugname]][["ED"]], row.names = NULL, check.names = FALSE), stringsAsFactors = FALSE)
      ECDRList[[samplename]][[drugname]] <- ECDRList[[samplename]][[drugname]][c("ED", "Estimate")]
      # Clean-up of ED labels
      ECDRList[[samplename]][[drugname]]$ED <- with(ECDRList[[samplename]][[drugname]] , sapply(strsplit(ED, ":"), `[`, 3))
      # Add columns with cell line and drug name
      ECDRList[[samplename]][[drugname]]$Sample <- samplename
      ECDRList[[samplename]][[drugname]]$Drug <- drugname
      
      
      # Convert list from long to wide-format
      ECDRList[[samplename]][[drugname]] <- stats::reshape(ECDRList[[samplename]][[drugname]], idvar = c("ED", "Estimate"), v.names = "Estimate", timevar = "ED", direction = "wide", new.row.names = NULL, sep = "")
      names(ECDRList[[samplename]][[drugname]]) <- gsub("Estimate", "ED", names(ECDRList[[samplename]][[drugname]]))
      # Combine rows dropping empty columns
      ECDRList[[samplename]][[drugname]] <- as.data.frame(lapply(ECDRList[[samplename]][[drugname]], function(x) x[!is.na(x)][1]), stringsAsFactors = FALSE)
 
    }
    
    ECDRList[[samplename]] <- do.call(rbind, setNames(ECDRList[[samplename]], NULL)) 
    

    # Create subfolder if it does not exist
    if(!file.exists(file.path(.saveto, "results", "Dynamic Drug Activity Range"))){ dir.create(file.path(.saveto, "results", "Dynamic Drug Activity Range"), showWarnings = FALSE, recursive = TRUE) }
    
    
    if(samplename == utils::tail(names(efs), n = 1)){
      ECDRList <- do.call(rbind, setNames(ECDRList, NULL))
      
      # Calculate the minimum, mean and maximum Efor each ED
      ECDRList$minED10 <- stats::ave(x = ECDRList$ED10, ECDRList$Drug, FUN = min)
      ECDRList$maxED90 <- stats::ave(x = ECDRList$ED90, ECDRList$Drug, FUN = max)
      ECDRList$avgED10 <- stats::ave(x = ECDRList$ED10, ECDRList$Drug, FUN = mean)
      ECDRList$avgED90 <- stats::ave(x = ECDRList$ED90, ECDRList$Drug, FUN = mean)
      ECDRList$medED10 <- stats::ave(x = ECDRList$ED10, ECDRList$Drug, FUN = median)
      ECDRList$medED90 <- stats::ave(x = ECDRList$ED90, ECDRList$Drug, FUN = median)
      
      ECDRList$Unit <- unique(efs[[samplename]][[drugname]][["drm"]][["origData"]]$Unit)
      
      ECDRList <- ECDRList[with(ECDRList, order(Drug)),]

      
      # Export the data set prior EC50 calculation
      write.csv2(ECDRList, file = file.path(.saveto, "results", "Dynamic Drug Activity Range", "Dynmaic Drug Activity Range (DDAR).csv"), row.names = FALSE, quote = FALSE)

    }
    
    if(samplename == utils::tail(names(efs), n = 1L)){
      cat('\r', "Finished extracting dynamic drug-activity range (DDAR).", strrep(" ", 100), '\n\n', sep = "")
      rm(samplename, drugname)
    }
    
  }
 
  
  
  # Plotting dose response curves with the dynamic range.
  # Plot curves as composite for each cell line with all the drugs and export as png.
  
  cat("Plotting dynamic drug-activity range (DDAR).", '\n', sep = "")
  
  for(samplename in names(efs)){
    
    cat('\r', " > Plotting dynamic range based on single drug viability curves for ", samplename, ".", strrep(" ", 100), sep = "")
    
    # Create subfolder if it does not exist
    if(!file.exists(file.path(.saveto, "results", "graphs/dynamic range", "by sample"))){ dir.create(file.path(.saveto, "results", "graphs/dynamic range", "by sample"), showWarnings = FALSE, recursive = TRUE) }
    
    

    grDevices::png(filename = file.path(.saveto, "results", "graphs/dynamic range", "by sample", paste(samplename, "dynamic range (viability).png", sep = " ")), 
        width = 7920, height = 4080, units = "px", pointsize = 24)
    op <- graphics::par(mfrow = c(7, 10), oma = c(4,0,7,0))
    
    for(drugname in sort(names(efs[[samplename]]))){
      
      drm.data <- subset(efs[[samplename]][[drugname]][["drm"]][["origData"]], Drug == drugname)
      drm.model <- efs[[samplename]][[drugname]][["drm"]]
      
      EC10 <- efs[[samplename]][[drugname]][["ED"]][paste("e", drugname, "10", sep = ":"), "Estimate"]
      EC50 <- efs[[samplename]][[drugname]][["ED"]][paste("e", drugname, "50", sep = ":"), "Estimate"]
      EC90 <- efs[[samplename]][[drugname]][["ED"]][paste("e", drugname, "90", sep = ":"), "Estimate"]
      
      # Plot all points (replicates) by type = all, or average of points by type = average
      graphics::par(bg="white", fg="black", col="black", col.axis="black", col.lab="black", col.main="black", col.sub="black") 
      graphics::plot(efs[[samplename]][[drugname]][["drm"]], level = drugname, log = "x", type = "average", main = drugname, 
                     xlab = paste("Concentration (", unique(efs[[samplename]][[drugname]][["drm"]][["origData"]]$Unit), ")", sep = ""), ylab = "Viability", axes = FALSE
           # set an individual scale for plotting
           , ylim=c(0, 1), pch=1, cex = 1.0
      )
      # graphics::grid(10,10, lwd = 1)
      
      if(samplename == head(names(efs), n = 1)){axis(1, at=unique(drm.data$Drug.Concentration), labels=format(unique(drm.data$Drug.Concentration), trim = TRUE, digits = 2, scientific = FALSE, drop0trailing = TRUE))}
      if(samplename == head(names(efs), n = 1)){axis(2)}
      
      # Indicate the measured EC50
      graphics::points(EC50, efs[[samplename]][[drugname]][["summary"]][["coefficients"]]["Upper Limit:(Intercept)", "Estimate"]-((
        efs[[samplename]][[drugname]][["summary"]][["coefficients"]]["Upper Limit:(Intercept)", "Estimate"]-
          efs[[samplename]][[drugname]][["summary"]][["coefficients"]]["Lower Limit:(Intercept)", "Estimate"])/2),
        type="p", pch=16, cex = 1.0, col="red")
      
      # Indicate the observed dynamic range between the EC10 and EC90
      graphics::abline(v = EC10, col="grey", lwd = 1.2, lty = 3)
      graphics::abline(v = EC90, col="grey", lwd = 1.2, lty = 3)
      rect(EC10,par("usr")[3],EC90,par("usr")[4],col=rgb(0.89,0.89,0.89,alpha=0.3),lty=0)
      
      graphics::legend("bottom", c(paste("EC10: ", format(EC10, digits = 3), "  EC50: ", format(EC50, digits = 3), "  EC90: ", format(EC90, digits = 3), sep = "")),
                       inset=c(0,0.99), xpd=TRUE, horiz=TRUE, cex = 0.9, col = "black", bty = "n")
      

      # Indicating the doses used for a given drug and the theoretical dynamic range
      for(n in 1:length(unique(drm.data$Drug.Concentration))){
        graphics::points(sort(unique(drm.data$Drug.Concentration))[n], 
               suppressWarnings(drc::PR(efs[[samplename]][[drugname]]$drm, sort(unique(drm.data$Drug.Concentration))[n])), 
               type="p", pch=16, cex = 1.2, col=rgb(0.0,0.0,0.0,alpha=1))
        
        graphics::abline(v = sort(unique(drm.data$Drug.Concentration))[n], col=rgb(0.0,0.0,0.0,alpha=1), lwd = 1.2, lty = 3)
        
      }
      
      
      #   # Import a list of maximum possible doses
      #   listofDoses <- read.csv2(file=file.path(libDirectory, "listofDoses.csv"), check.names=FALSE, header=TRUE, stringsAsFactors=FALSE, na.strings="", sep=";", dec=",", skip=0)
      #   # Section checking if the EC50 is above the highest dose
      #   if(EC50 > as.numeric(listofDoses$`6th Dose`[which(listofDoses$Drug == drugname)])){
      #     graphics::box("plot", lwd = 2.4, col="orange")}
      # 
      #   # Section checking if the EC10 is below the lowest dose
      #   if(EC10 < as.numeric(listofDoses$`1st Dose`[which(listofDoses$Drug == drugname)])){
      #     graphics::box("plot", lwd = 2.4, col="orange")}
      
      # Section checking if the Std. Error is larger than the EC50
      if(is.na(efs[[samplename]][[drugname]][["summary"]][["coefficients"]]["ED50:(Intercept)", "Std. Error"]) &
         EC50 > as.numeric(max(sort(unique(drm.data$Drug.Concentration))))){
        graphics::box("plot", lwd = 2.4, col="orange")}
      
      # 
      #   # Section checking if there is less than 25% killing
      #   if(efs[[samplename]][[drugname]][["summary"]][["coefficients"]]["Upper Limit:(Intercept)", "Estimate"] -
      #      efs[[samplename]][[drugname]][["summary"]][["coefficients"]]["Lower Limit:(Intercept)", "Estimate"] < 0.2){
      #     graphics::box("plot", lwd = 2.4, col="red")}
      
      # Section checking if the slope is positive or negative
      if(efs[[samplename]][[drugname]][["summary"]][["coefficients"]]["Slope:(Intercept)", "Estimate"] < 0.0){
        graphics::box("plot", lwd = 2.4, col="red")}
      
      # Section checking if the EC50 is outside the dose range
      if(EC50 > max(efs[[samplename]][[drugname]][["drm"]][["dataList"]][["dose"]]) | EC50 < min(efs[[samplename]][[drugname]][["drm"]][["dataList"]][["dose"]])){
        graphics::box("plot", lwd = 2.4, col="orange")}
      
      # Section checking if response is linear
      if(efs[[samplename]][[drugname]][["drm"]][["coefficients"]]["Lower Limit:(Intercept)"] == efs[[samplename]][[drugname]][["drm"]][["coefficients"]]["Upper Limit:(Intercept)"]){
        graphics::box("plot", lwd = 2.4, col="red")}
      
      
    }
    
    plot.new()
    graphics::legend("topleft", c("observed dynamic range (DR)"), fill=c(rgb(0.89,0.89,0.89,alpha=0.3)), xpd = TRUE, horiz = FALSE, inset = c(0, 0.1), bty = "n", cex = 1.5)
    graphics::legend("topleft", c("EC50"), col=rgb(1,0,0,alpha=1), xpd = TRUE, horiz = FALSE, inset = c(0.0125, 0.25), pch=16, bty = "n", cex = 1.5)
    graphics::title(paste(samplename, ": Dynamic Drug Activity Range", sep = ""), outer = TRUE, cex.main = 4)
    graphics::mtext(paste("  Package: ", "drc", " (v", packageVersion("drc"), ")   ",  sep = ""), side = 1, line = 2, adj = 1, cex = 1, col = "black", outer = TRUE) 
    grDevices::dev.off()
    

    if(samplename == utils::tail(names(efs), n = 1L)){
      cat('\r', "Finished plotting dynamic drug-activity range (DDAR).", strrep(" ", 100), '\n\n', sep = "")
      rm(samplename, drugname, EC10, EC50, EC90, n, op)
    }
    
  }
  
  
  
  
  # Plotting the dynamic range for all samples.
  # Plot curves as composite with all cell lines for each of the drugs and export as png.
  
  cat("Plotting dynamic drug-activity range (DDAR) based on single drug viability curves as composite for all samples.", '\n', sep = "")
  
  
  grDevices::png(filename = file.path(.saveto, "results", "graphs/dynamic range", "dynamic range (DR) (all cell lines) graph.png"), 
      width = 7920, height = 4080, units = "px", pointsize = 24)
  op <- graphics::par(mfrow = c(7, 10), oma = c(4,0,7,0))
  
  for(drugname in sort(names(efs[[1]]))){
    
    cat('\r', " > Plotting dynamic range based on single drug viability curves for ", drugname, ".", strrep(" ", 100), sep = "")
    
    
    for(samplename in names(efs)){
      
      drm.data <- subset(efs[[samplename]][[drugname]][["drm"]][["origData"]], Drug == drugname)
      drm.model <- efs[[samplename]][[drugname]][["drm"]]
      
      # Plot all points (replicates) by type = all, or average of points by type = average
      graphics::plot(drm.model, log = "x", type = "none", main = drugname, xlab = paste("Concentration (", unique(drm.data$Unit), ")", sep = ""), ylab = "Viability", axes = FALSE,
           lty = 1, lwd = 0.5, cex = 1, pch = 16, add = ifelse(samplename == head(names(efs), n = 1), FALSE, TRUE),
           # plot the full concentration range possible
           xlim=c(0, max(drm.data$Drug.Concentration)),
           # set an individual scale for plotting
           ylim=c(0, 1), col = FALSE, legend = FALSE)
      
      if(samplename == head(names(efs), n = 1)){axis(1, at=unique(drm.data$Drug.Concentration), labels=format(unique(drm.data$Drug.Concentration), trim = TRUE, digits = 2, scientific = FALSE, drop0trailing = TRUE))}
      if(samplename == head(names(efs), n = 1)){axis(2)}
      
    }
    
    # Indicate the observed dynamic range between the EC10 and EC90
    # Note, if LL2.4 is used, the scale is logarithmic: use log(ED50)
    graphics::abline(v = unique(ECDRList$avgED10[which(ECDRList$Drug == drugname)]), col="grey", lwd = 1.2, lty = 3)
    graphics::abline(v = unique(ECDRList$avgED90[which(ECDRList$Drug == drugname)]), col="grey", lwd = 1.2, lty = 3)
    rect(unique(ECDRList$avgED10[which(ECDRList$Drug == drugname)]),
         par("usr")[3],
         unique(ECDRList$avgED90[which(ECDRList$Drug == drugname)]),
         par("usr")[4],col=rgb(0.00,0.89,0.00,alpha=0.1),lty=0)
    
    # graphics::legend("bottomleft", adj = c(0.05, 0.5), paste("DR: EC10 (median): ", format(unique(ECDRList$medED10[which(ECDRList$Drug == drugname)]), digits = 3),
    #                      " ", "EC90 (median): ", format(unique(ECDRList$medED90[which(ECDRList$Drug == drugname)]), digits = 3), sep = ""),
    #        cex = 0.8, col = "black", bty = "n")
    
    
    
    # Indicating the doses used for a given drug and the theoretical dynamic range
    for(n in 1:length(unique(drm.data$Drug.Concentration))){
      graphics::abline(v = sort(unique(drm.data$Drug.Concentration))[n], col=rgb(0.0,0.0,0.0,alpha=1), lwd = 1.2, lty = 3)
    }
    
    rect(sort(unique(drm.data$Drug.Concentration))[2],
         par("usr")[3],
         sort(unique(drm.data$Drug.Concentration))[length(sort(unique(drm.data$Drug.Concentration)))-1],
         par("usr")[4],col=rgb(1.00,0.00,0.00,alpha=0.05),lty=0)
    
    

    if(drugname == utils::tail(sort(names(efs[[1]])), n = 1)){
      plot.new()
      graphics::legend("topleft", c("observed dynamic range (DR)", "expected dynamic range (DR)"), fill=c(rgb(0.00,0.89,0.00,alpha=0.1), rgb(1.00,0.00,0.00,alpha=0.05)), xpd = TRUE, horiz = FALSE, inset = c(0, 0.1), bty = "n", cex = 1.5)
      graphics::legend("topleft", c("max. Dose"), col=c(rgb(1,0,0,alpha=1)), xpd = TRUE, horiz = FALSE, inset = c(0, 0.4), lwd = 1.2, lty = 3, bty = "n", cex = 1.5)
      graphics::title(paste("Dynamic Drug Activity Range (all cell lines)", sep = " "), outer = TRUE, cex.main = 4)
      graphics::mtext(paste("  Package: ", "drc", " (v", packageVersion("drc"), ")   ",  sep = ""), side = 1, line = 2, adj = 1.0, cex = 1, col = "black", outer = TRUE) 
      grDevices::dev.off()
      
      cat('\r', "Finished plotting dynamic drug-activity range (DDAR) as composite.", strrep(" ", 100), '\n', sep = "")
      rm(drm.data, samplename, drugname, n, op)
    }
    
  }
  
  
  
  
  # Plotting the dynamic range for all drugs
  # Plotting individual dose response curves for all cell line and by each drug
  
  cat("Plotting dynamic drug-activity range (DDAR) based on single drug viability curves by drug.", '\n', sep = "")
  
  for(drugname in sort(names(efs[[1]]))){
    
    cat('\r', " > Plotting dynamic range based on single drug viability curves for ", drugname, ".", strrep(" ", 100), sep = "")
    
    # Create subfolder if it does not exist
    if(!file.exists(file.path(.saveto, "results", "graphs/dynamic range", "by drug"))){ dir.create(file.path(.saveto, "results", "graphs/dynamic range", "by drug"), showWarnings = FALSE, recursive = TRUE) }


    grDevices::png(filename = file.path(.saveto, "results", "graphs/dynamic range", "by drug", paste(gsub("/", " ", drugname), " dynamic range (viability)", ".png", sep = "")),
        width = 1920*0.8, height = 1080*0.8, units = "px", pointsize = 24)
    
    graphics::par(mar=c(5.1, 4.1, 3.1, 8.1), xpd=TRUE)
    
    for(samplename in names(efs)){
      
      drm.data <- subset(efs[[samplename]][[drugname]][["drm"]][["origData"]], Drug == drugname)
      drm.model <- efs[[samplename]][[drugname]][["drm"]]

      
      # Plot all points (replicates) by type = all, or average of points by type = average
      graphics::plot(drm.model, log = "x", type = "none", main = drugname, xlab = paste("Concentration (", unique(drm.data$Unit), ")", sep = ""), ylab = "Viability", axes=FALSE,
           lty = match(samplename, names(efs)), lwd = 1.5, cex = 1, cex.axis = 0.8, cex.lab = 0.8, pch = 16, bty="l",
           add = ifelse(samplename == head(names(efs), n = 1), FALSE, TRUE),
           # plot the full concentration range possible
           xlim=c(0, max(drm.data$Drug.Concentration)),
           # set an individual scale for plotting
           ylim=c(0, 1), col = grDevices::rainbow(18)[match(samplename, names(efs))], legend = FALSE
      )
      
      if(samplename == head(names(efs), n = 1)){axis(1, at=unique(drm.data$Drug.Concentration), labels=format(unique(drm.data$Drug.Concentration), trim = TRUE, digits = 2, scientific = FALSE, drop0trailing = TRUE))}
      if(samplename == head(names(efs), n = 1)){axis(2)}
      
    }
    
    graphics::legend("right", inset = c(-0.15, 0), legend = unique(names(efs)), xpd = TRUE, 
           horiz = FALSE, col = grDevices::rainbow(18), lty = 1:18, lwd = 1,5, cex = 0.8, bty = "n")
    
    # indicate the dynamic range
    # graphics::abline(v = unique(ECDRList$medED10[which(ECDRList$Drug == drugname)]), col="grey", lwd = 1.2, lty = 3)
    # graphics::abline(v = unique(ECDRList$medED90[which(ECDRList$Drug == drugname)]), col="grey", lwd = 1.2, lty = 3)
    # rect(unique(ECDRList$medED10[which(ECDRList$Drug == drugname)]),
    #      par("usr")[3],
    #      unique(ECDRList$medED90[which(ECDRList$Drug == drugname)]),
    #      par("usr")[4],col=rgb(0.00,0.89,0.00,alpha=0.1),lty=0)
    
    
    
    # Indicating the doses used for a given drug and the theoretical dynamic range
    for(n in 1:length(unique(drm.data$Drug.Concentration))){
      graphics::abline(v = sort(unique(drm.data$Drug.Concentration))[n], col=rgb(0.0,0.0,0.0,alpha=0.8), lwd = 0.5, lty = 3, xpd=FALSE)
    }
    
    rect(sort(unique(drm.data$Drug.Concentration))[2],
         par("usr")[3],
         sort(unique(drm.data$Drug.Concentration))[length(sort(unique(drm.data$Drug.Concentration)))-1],
         par("usr")[4],col=rgb(1.00,0.00,0.00,alpha=0.05),lty=0)

    
    grDevices::dev.off()
    
    if(drugname == utils::tail(sort(names(efs[[1]])), n = 1)){
      cat('\r', "Finished plotting dynamic drug-activity range (DDAR) by drug.", strrep(" ", 100), '\n', sep = "")
      rm(drm.data, drm.model, samplename, drugname, n)
    }
  }
  
  
  
  
  
  
  # Plotting dose response curves manually
  # Note: drc is only used to generate the base of the plots, 
  # the points and lines are then added manually
  
  cat("Plotting dynamic drug-activity range (DDAR) based on unfitted single drug viability curves.", '\n', sep = "")
  
  
  grDevices::png(filename = file.path(.saveto, "results", "graphs/dynamic range", "dynamic range (DR) (all cell lines) graph (unfitted wo points).png"), 
      width = 7920, height = 4080, units = "px", pointsize = 24)
  op <- graphics::par(mfrow = c(7, 10), oma = c(4,0,7,0))
  
  for(drugname in sort(names(efs[[1]]))){
    
    cat('\r', " > Plotting dynamic range based on single drug viability curves for ", drugname, ".", strrep(" ", 100), sep = "")
    
    
    for(samplename in names(efs)){
      
      drm.data <- subset(efs[[samplename]][[drugname]][["drm"]][["origData"]], Drug == drugname)
      drm.data$`Viability (Mean)` <- with(drm.data, stats::ave(Viability, Drug.Concentration, FUN = mean))
      drm.data$`Viability (Mean)` <- with(drm.data, ifelse(`Viability (Mean)` > 1, 1, `Viability (Mean)`))
      drm.model <- efs[[samplename]][[drugname]][["drm"]]

      
      # Plot all points (replicates) by type = all, or average of points by type = average
      graphics::plot(drm.model, log = "x", type = "none", main = drugname, xlab = paste("Concentration (", unique(drm.data$Unit), ")", sep = ""), ylab = "Viability", axes = FALSE,
           lty = 1, lwd = 0, cex = 1, pch = 1, add = ifelse(samplename == head(names(efs), n = 1), FALSE, TRUE),
           # plot the full concentration range possible
           xlim=c(0, max(drm.data$Drug.Concentration)),
           # set an individual scale for plotting
           ylim=c(0, 1),  bty="n", col = FALSE, legend = FALSE)
      
      if(samplename == head(names(efs), n = 1)){axis(1, at=unique(drm.data$Drug.Concentration), labels=format(unique(drm.data$Drug.Concentration), trim = TRUE, digits = 2, scientific = FALSE, drop0trailing = TRUE))}
      if(samplename == head(names(efs), n = 1)){axis(2)}
      # graphics::grid(nx = NULL, ny = nx, col = "lightgrey", lty = "dotted", lwd = par("lwd"), equilogs = TRUE)
      # graphics::abline(v=outer((1:10),(10^(-5:5))), col="grey85", lty = "dotted", lwd = 0.5)
      # graphics::abline(h=seq(0,1,0.1), col="grey85", lty = "dotted", lwd = 0.5)
      # graphics::points(drm.data$Drug.Concentration, drm.data$`Viability (Mean)`, type="p", pch=16, cex = 0.6, col="#FF0066")
      lines(drm.data$Drug.Concentration[order(drm.data$Drug.Concentration)], drm.data$`Viability (Mean)`[order(drm.data$Drug.Concentration)], lwd = 0.6)
      # lines(drm.data$Drug.Concentration[order(drm.data$Drug.Concentration)], stats::predict(lo, drm.data$Drug.Concentration[order(drm.data$Drug.Concentration)]), lwd = 0.5)
      # lines(lo, lwd = 0.5)
    }
    
    if(drugname == utils::tail(sort(names(efs[[1]])), n = 1)){
      plot.new()
      graphics::title(paste("Dynamic Drug Activity Range (all cell lines)", sep = " "), outer = TRUE, cex.main = 4)
      grDevices::dev.off()
      
      cat('\r', "Finished plotting dynamic drug-activity range (DDAR) based on unfitted single drug viability curves.", strrep(" ", 100), '\n', sep = "")
      rm(samplename, drugname, drm.data, drm.model, op)
    }
    
  }
  
  
  cat("Plotting dynamic drug-activity range (DDAR) based on unfitted single drug viability curves with points.", '\n', sep = "")
  
  
  grDevices::png(filename = file.path(.saveto, "results", "graphs/dynamic range", "dynamic range (DR) (all cell lines) graph (unfitted w points).png"), 
      width = 7920, height = 4080, units = "px", pointsize = 24)
  op <- graphics::par(mfrow = c(7, 10), oma = c(4,0,7,0))
  
  for(drugname in sort(names(efs[[1]]))){
    
    cat('\r', " > Plotting dynamic range based on single drug viability curves for ", drugname, ".", strrep(" ", 100), sep = "")
    
    
    for(samplename in names(efs)){
      
      drm.data <- subset(efs[[samplename]][[drugname]][["drm"]][["origData"]], Drug == drugname)
      drm.data$`Viability (Mean)` <- with(drm.data, stats::ave(Viability, Drug.Concentration, FUN = mean))
      drm.data$`Viability (Mean)` <- with(drm.data, ifelse(`Viability (Mean)` > 1, 1, `Viability (Mean)`))
      drm.model <- efs[[samplename]][[drugname]][["drm"]]
      
      
      # Plot all points (replicates) by type = all, or average of points by type = average
      graphics::plot(drm.model, log = "x", type = "none", main = drugname, xlab = paste("Concentration (", unique(drm.data$Unit), ")", sep = ""), ylab = "Viability", axes = FALSE,
           lty = 1, lwd = 0, cex = 1, pch = 1, add = ifelse(samplename == head(names(efs), n = 1), FALSE, TRUE),
           # plot the full concentration range possible
           xlim=c(0, max(drm.data$Drug.Concentration)),
           # set an individual scale for plotting
           ylim=c(0, 1),  bty="n", col = FALSE, legend = FALSE)
      
      if(samplename == head(names(efs), n = 1)){axis(1, at=unique(drm.data$Drug.Concentration), labels=format(unique(drm.data$Drug.Concentration), trim = TRUE, digits = 2, scientific = FALSE, drop0trailing = TRUE))}
      if(samplename == head(names(efs), n = 1)){axis(2)}
      graphics::grid(nx = NULL, ny = "", col = "lightgrey", lty = "dotted", lwd = par("lwd"), equilogs = TRUE)
      graphics::abline(v=outer((1:10),(10^(-5:5))), col="grey80", lty = "dotted", lwd = 0.3)
      graphics::abline(h=seq(0,1,0.1), col="grey80", lty = "dotted", lwd = 0.3)
      lines(drm.data$Drug.Concentration[order(drm.data$Drug.Concentration)], drm.data$`Viability (Mean)`[order(drm.data$Drug.Concentration)], lwd = 0.6)
      graphics::points(drm.data$Drug.Concentration, drm.data$`Viability (Mean)`, type="p", pch=16, cex = 0.6, col="#FF0066")
      
      # lines(drm.data$Drug.Concentration[order(drm.data$Drug.Concentration)], stats::predict(lo, drm.data$Drug.Concentration[order(drm.data$Drug.Concentration)]), lwd = 0.5)
      # lines(lo, lwd = 0.5)
    }
    
    if(drugname == utils::tail(sort(names(efs[[1]])), n = 1)){
      plot.new()
      graphics::title(paste("Dynamic Drug Activity Range (all cell lines)", sep = " "), outer = TRUE, cex.main = 4)
      grDevices::dev.off()
      
      cat('\r', "Finished plotting dynamic drug-activity range (DDAR) based on unfitted single drug viability curves with points.", strrep(" ", 100), '\n', sep = "")
      rm(samplename, drugname, drm.data, drm.model, op)
    }
    
  }
  
  
  message("Drug activity range plots saved to: ", '\n', file.path(.saveto, "results", "graphs/dynamic range"), strrep(" ", 100), sep = "")
  
   
  # Return object of class DRM
  class(ECDRList) <- "ECDR"
  return(ECDRList)
  
}