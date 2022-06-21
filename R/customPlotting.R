#' Custom plotting of drug-dose response
#'
#' @description
#' An complementary component of the modular library \pkg{metascreen} providing a set of custom plots for the visualization of drug-dose responses.
#' \emph{\code{customPlotting}} is a function that generates a set of plots that provide an overview of the dose-response between drugs and samples.
#' Furthermore, it plots the dose-response matrix for drug combinations along with the dose-response curves for each individual drug pair.
#' 
#' @param processedData an object of class 'processedData'.
#' @param doseRespModel an object of class 'drm'.
#' @param .saveto string; path to a folder location where the results are saved to.
#' 
#' 
#' @details The function \code{customPlotting} is used to provide an overview of the responses of all treatments between individual drugs and samples. This is of particularly benefit for large drug screens
#' in which a large number of drugs and samples have been screened.
#' 
#' Files are saved either to the specified location or the default working environment, with the corresponding folder structure: 'results/graphs/custom plots'
#' 
#' @examples
#' \donttest{\dontrun{
#' # Plot dose response curves and matrix
#' customPlotting(processedData, doseRespModel, .saveto = "path/to/folder/")
#' }}
#'
#' @keywords drug screen analysis dose response curve matrix
#' 
#' @importFrom utils capture.output tail
#' @importFrom stats ave predict
#' @importFrom drc drm LL.4 drmc ED
#' @importFrom ggplot2 ggplot ggsave labs
#' @importFrom grid gpar textGrob unit viewport
#' @importFrom gridExtra grid.arrange
#' 
#' @export

customPlotting <- function(processedData, doseRespModel, .saveto){
  
  # Check, if the data has been provided as an object of class S3:processedData
  if(missing(processedData)){stop("Data missing! Please provide a processed data set.", call. = TRUE)}
  if(class(processedData) != "processedData"){
    stop("First argument needs to be data of class 'processedData'!", call. = TRUE)
  } else {
    # Extract all the necessary data from the class S3 object
    analysisData <- processedData[["analysisData"]]
    synDataset <- processedData[["splitDataset"]]
    synDRM <- processedData[["doserespMatrix"]]
  }
  
  # Check, if the data has been provided as an object of class S3:drm
  if(missing(doseRespModel)){stop("Data missing! Please provide a data set with the dose response models.", call. = TRUE)}
  if(class(doseRespModel) != "drm"){
    stop("Second argument needs to be data of class 'drm'!", call. = TRUE)
  } else {
    # Extract all the necessary data from the class S3 object
    efs <- doseRespModel
  }
  
  
  # Check, if a folder location has been provided
  if(missing(.saveto)){ .saveto <- getwd() }
  # Create folder, if it does not exist
  if(!file.exists(file.path(.saveto))){ dir.create(file.path(.saveto), showWarnings = FALSE, recursive = TRUE) }
  if(!utils::file_test("-d", .saveto)){stop("in '.saveto'. Argument needs to be a valid folder location.", call. = TRUE)}
  
  
  
  if(!exists("synPlots")){synPlots = list()}
  
  # SECTION A: PLOTTING D-R CURVES AS COMPOSITE BY SAMPLE ##################################################
  # Plotting all single dose-response curves
  
  cat("Plotting dose-response curves by sample.", '\n', sep = "")
  
  for(samplename in names(efs)){
    
    # Skip samples that have already been plotted
    # if(file_test("-f", file.path(.saveto, "results", "graphs/custom plots/dose-response curves", "by sample", paste(samplename, "dose-response curves.png", sep = " ")))){next}
    
    cat('\r', "Plotting dose-resonse curves for ", samplename, ".", strrep(" ", 100), sep = "")
    
    
    # Create subfolder if it does not exist
    if(!file.exists(file.path(.saveto, "results", "graphs/custom plots/dose-response curves", "by sample"))){ dir.create(file.path(.saveto, "results", "graphs/custom plots/dose-response curves", "by sample"), showWarnings = FALSE, recursive = TRUE) }
    
    
    
    # Select all drugs sorted according to the list of drugs
    for(drugname in sort(names(efs[[samplename]]))){
      
      cat('\r', " > Plotting dose-resonse curve: ", samplename, ": ", drugname, strrep(" ", 100), sep = "")
      
      # Filter data set by drug, select columns drug, dose and inhibition and order data by dose
      .drm.data <- subset(efs[[samplename]][[drugname]][["drm"]][["origData"]], Drug == drugname)
      .drm.data <- .drm.data[c("Drug", "Drug.Concentration", "Unit", "Inhibition")]
      names(.drm.data)[match("Drug.Concentration", names(.drm.data))] <- "Dose"
      .drm.data <- .drm.data[with(.drm.data, order(Dose)),]
      
      # Calculate the median inhibition for each dose
      .drm.data$Median <- ave(x = .drm.data$Inhibition, .drm.data$Dose, FUN = median)
      
      # Extract the dose response model
      .drm.model <- efs[[samplename]][[drugname]][["drm"]]
      
      # Run predictions to estimate the confidence interval
      .drm.prediction <- expand.grid(Dose=exp(seq(log(min(.drm.data$Dose)), log(max(.drm.data$Dose)), length=100)))
      
      .drm.prediction$p <- 1-suppressWarnings(stats::predict(.drm.model, newdata=.drm.prediction, interval="confidence")[,1])
      .drm.prediction$pmin <- 1-suppressWarnings(stats::predict(.drm.model, newdata=.drm.prediction, interval="confidence")[,2])
      .drm.prediction$pmax <- 1-suppressWarnings(stats::predict(.drm.model, newdata=.drm.prediction, interval="confidence")[,3])
      
      # Don't plot predictions lower than the minimum inhibition or higher than 1, 
      # keeping the scale of the plot to the actual data rather than the predictions.
      # .drm.prediction <- transform(.drm.prediction, pmin = ifelse(pmin < min(.drm.data$Inhibition), min(.drm.data$Inhibition)-0.05, pmin))
      .drm.prediction <- transform(.drm.prediction, pmin = ifelse(pmin < 0, 0, pmin), pmax = ifelse(pmax > 1, 1, pmax))
      
      # Replace pmax and pmin predictions that are NAs with prediction numbers
      .drm.prediction <- transform(.drm.prediction, pmin = ifelse(is.na(pmin), p, pmin), pmax = ifelse(is.na(pmax), p, pmax))
      
      
      synPlots[["dose.response.curve"]][[samplename]][[drugname]] <- ggplot2::ggplot(data = .drm.data, aes(x = Dose, y = Inhibition)) +
        geom_point(data = .drm.data, aes(x = Dose, y = Median), shape = 16, size = 0.8) +
        geom_point(shape = 1, size = 0.8) +
        ggplot2::geom_line(data = .drm.prediction, aes(x = Dose, y = p), size = 0.3) +
        # geom_ribbon(data = .drm.prediction, aes(x = Dose, y = p, ymin = pmin, ymax = pmax), alpha = 0.05) +
        ggplot2::scale_x_continuous(name = paste("Concentration (", unique(.drm.data$Unit), ")", sep = ""), breaks=sort(unique(.drm.data$Dose)),
                           labels=format(round(sort(unique(.drm.data$Dose)), 4), nsmall = 4, scientific = FALSE, drop0trailing = TRUE)) +
        ggplot2::scale_y_continuous(name = "Inhibition", breaks=c(-1.00, -0.75, -0.50, -0.25, 0.00, 0.25, 0.50, 0.75, 1.00),
                           labels=c("-1.00", "-0.75", "-0.50", "-0.25", "0.00", "0.25", "0.50", "0.75", "1.00"),
                           limits = c(NA, 1),  expand = c(0,0)) +
        xlab(paste("Concentration (", unique(.drm.data$Unit), ")", sep = "")) +
        ylab("Inhibition") +
        labs(title = drugname) +
        ggplot2::annotate(x=sort(unique(.drm.data$Dose))[2], xend=sort(unique(.drm.data$Dose))[5], y=min(.drm.data$Inhibition)-0.05, yend=min(.drm.data$Inhibition)-0.05, 
                 colour="black", lwd=0.5, geom="segment") +
        ggplot2::coord_trans(x="log",ylim = c(min(.drm.data$Inhibition)-0.05, 1), clip = "on") +
        theme(plot.title = element_text(color="black", size=12, face="bold", hjust = 0.5),
              axis.title = element_text(color="black", size=12, face="plain", hjust = 0.5),
              axis.title.x = element_text(margin = ggplot2::margin(t = 10, r = 0, b = 0, l = 0)),
              axis.title.y = element_text(margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 0)),
              legend.title = element_text(color="black", size=10, face="plain", hjust = 0.0),
              legend.text = element_text(color="black", size=10, face="plain", hjust = 0.0),
              plot.caption = element_text(size = 10),
              text = element_text(size = 12),
              axis.line = element_line(color='black', size=0),
              panel.grid.major = element_line(color='lightgrey', linetype = "dotted"),
              panel.grid.minor = element_line(color='lightgrey', linetype = "dotted"),
              panel.background = element_rect(fill=alpha('white', 2)),
              panel.border = element_rect(colour = "black", fill=NA, size=0),
              legend.position = "right"
        )
      

      rm(.drm.data, .drm.model, .drm.prediction)
      
    }
    
    cat('\r', " > Saving dose-response curves for ", samplename, ".", strrep(" ", 100), sep = "")

    g <- gridExtra::grid.arrange(grobs=synPlots[["dose.response.curve"]][[samplename]], vp=grid::viewport(width=1, height=1), nrow = 7, ncol = 10, top = grid::textGrob(paste(samplename, ": Dose-Response Curves\n", sep = ""), gp=grid::gpar(fontsize = 16, fontface = "bold")))
    
    ggplot2::ggsave(filename = file.path(.saveto, "results", "graphs/custom plots/dose-response curves", "by sample", paste(samplename, "dose-response curves.png", sep = " ")),
           g, device = "png", width = 670, height = 345, units = "mm", dpi = 300, limitsize = FALSE)
    
    cat('\r', "Finished plotting dose-response curves for ", samplename, ".", strrep(" ", 100), '\n', sep = "")
    
    
    if(samplename == utils::tail(names(efs), n = 1L)){
      message('\n', "Dose-response curves saved to: ", '\n', file.path(.saveto, "results", "graphs/custom plots/dose-response curves", "by sample"), strrep(" ", 100), sep = "")
      rm(samplename, drugname, g)
    }
    
  }
  
  
  
  
  
  
  # SECTION B: PLOTTING D-R CURVES AS COMPOSITE BY DRUG ##################################################
  # Plotting all single dose-response curves
  
  cat("Plotting dose-response curves by drug", '\n', sep = "")
  
  for (drugname in sort(names(efs[[1]]))){
    
    # Skip samples that have already been plotted
    # if(file_test("-f", file.path(.saveto, "results", "graphs/custom plots/dose-response curves", "by drug", paste(drugname, "dose-response curves.png", sep = " ")))){next}
    
    cat('\r', " > Plotting dose-resonse curves for ", drugname, ".", strrep(" ", 100), sep = "")

    # Create subfolder if it does not exist
    if(!file.exists(file.path(.saveto, "results", "graphs/custom plots/dose-response curves", "by drug"))){ dir.create(file.path(.saveto, "results", "graphs/custom plots/dose-response curves", "by drug"), showWarnings = FALSE, recursive = TRUE) }
    

    
    # Pick plots from list for each cell line and assemble them by drug
    g <- gridExtra::grid.arrange(grobs=lapply(names(synPlots[["dose.response.curve"]]), function(x){synPlots[["dose.response.curve"]][[x]][[drugname]]+labs(title = x)}), vp=grid::viewport(width=1, height=1), nrow = 7, ncol = 9, top = grid::textGrob(paste(drugname, ": Dose-Response Curves\n", sep = ""), gp=grid::gpar(fontsize = 16, fontface = "bold")))
    
    ggplot2::ggsave(filename = file.path(.saveto, "results", "graphs/custom plots/dose-response curves", "by drug", paste(gsub("/", " ", drugname), "dose-response curves.png", sep = " ")),
           g, device = "png", width = 670*(10/9), height = 345*(10/9), units = "mm", dpi = 300, limitsize = FALSE)
    
    
    cat('\r', "Finished plotting dose-response curves for ", drugname, ".", strrep(" ", 100), '\n', sep = "")
    
    
    if(drugname == utils::tail(sort(names(efs[[1]])), n = 1L)){
      message('\n', "Dose-response curves saved to: ", '\n', file.path(.saveto, "results", "graphs/custom plots/dose-response curves", "by drug"), strrep(" ", 100), sep = "")
      rm(drugname, g)
    }
    
  }
  
  
  
  
  # SECTION C: Plotting Composite D-R CURVES + MATRIX ###############################
  # Plotting dose response and synergies as composites

  # Plot drug response matrix and export as png
  # Generating plots for each sample
  
  cat("Plotting dose-response curves and matrix.", '\n', sep = "")
  
  for(samplename in names(synDRM)){
    
    cat('\r', "Plotting dose-resonse curves and matrix for ", samplename, ".", strrep(" ", 100), sep = "")
    
    
    for (drugname in names(synDRM[[samplename]])){
      
      # Create subfolder if it does not exist
      if(!file.exists(file.path(.saveto, "results", "graphs/custom plots/dose-response", samplename, gsub("/", " ", drugname)))){ dir.create(file.path(.saveto, "results", "graphs/custom plots/dose-response", samplename, gsub("/", " ", drugname)), showWarnings = FALSE, recursive = TRUE) }
      
      
      for (drugpair in names(synDRM[[samplename]][[drugname]])){
        
        # Skip samples that have already been plotted
        # if(file_test("-f", file.path(.saveto, "results", "graphs/custom plots/dose-response", samplename, gsub("/", " ", drugname), paste(samplename, " ", gsub("/", " ", drugname), " + ", gsub("/", " ", drugpair), ".png", sep = "")))){next}
        
        cat('\r', " > Plotting dose-resonse for ", samplename, ": ", drugname, " + ", drugpair, strrep(" ", 100), sep = "")

        
        
        # PLOTTING DRUG RESPONSE CURVES
        for(i in c(drugname, drugpair)){
          
          # Filter data set by drug, select columns drug, dose and inhibition and order data by dose
          .drm.data <- subset(efs[[samplename]][[i]][["drm"]][["origData"]], Drug == i)
          .drm.data <- .drm.data[c("Drug", "Drug.Concentration", "Unit", "Inhibition")]
          names(.drm.data)[match("Drug.Concentration", names(.drm.data))] <- "Dose"
          .drm.data <- .drm.data[with(.drm.data, order(Dose)),]
          
          # Calculate the median inhibition for each dose
          .drm.data$Median <- ave(x = .drm.data$Inhibition, .drm.data$Dose, FUN = median)
          
          # Extract the dose response model
          .drm.model <- efs[[samplename]][[i]][["drm"]]
          
          # Run predictions to estimate the confidence interval
          .drm.prediction <- expand.grid(Dose=exp(seq(log(min(.drm.data$Dose)), log(max(.drm.data$Dose)), length=100)))
          
          .drm.prediction$p <- 1-suppressWarnings(stats::predict(.drm.model, newdata=.drm.prediction, interval="confidence")[,1])
          .drm.prediction$pmin <- 1-suppressWarnings(stats::predict(.drm.model, newdata=.drm.prediction, interval="confidence")[,2])
          .drm.prediction$pmax <- 1-suppressWarnings(stats::predict(.drm.model, newdata=.drm.prediction, interval="confidence")[,3])
          
          # Don't plot predictions lower than the minimum inhibition or higher than 1, 
          # keeping the scale of the plot to the actual data rather than the predictions.
          # .drm.prediction <- transform(.drm.prediction, pmin = ifelse(pmin < min(.drm.data$Inhibition), min(.drm.data$Inhibition)-0.05, pmin))
          # .drm.prediction <- transform(.drm.prediction, pmin = ifelse(pmin < 0, 0, pmin), pmax = ifelse(pmax > 1, 1, pmax))
          
          # Replace pmax and pmin predictions that are NAs with prediction numbers
          # .drm.prediction <- transform(.drm.prediction, pmin = ifelse(is.na(pmin), p, pmin), pmax = ifelse(is.na(pmax), p, pmax))
          

          .textscalefactor = 1.45
          
          synPlots[["dose.response.curve"]][[i]] <- ggplot2::ggplot(data = .drm.data, aes(x = Dose, y = Inhibition)) +
            ggplot2::geom_point(data = .drm.data, aes(x = Dose, y = Median), shape = 16) +
            ggplot2::geom_point(shape = 1) +
            ggplot2::geom_line(data=.drm.prediction, aes(x=Dose, y=p)) +
            # geom_ribbon(data=.drm.prediction, aes(x=Dose, y=p, ymin=pmin, ymax=pmax), alpha=0.05) +
            ggplot2::scale_x_continuous(name = paste("Concentration (", unique(.drm.data$Unit), ")", sep = ""), breaks=sort(unique(.drm.data$Dose)),
                               labels=format(round(sort(unique(.drm.data$Dose)), 4), nsmall = 4, scientific = FALSE, drop0trailing = TRUE)) +
            ggplot2::scale_y_continuous(name = "Inhibition", breaks=c(-1.00, -0.75, -0.50, -0.25, 0.00, 0.25, 0.50, 0.75, 1.00),
                               labels=c("-1.00", "-0.75", "-0.50", "-0.25", "0.00", "0.25", "0.50", "0.75", "1.00"),
                               limits = c(NA, 1),  expand = c(0,0)) +
            xlab(paste("Concentration (", unique(.drm.data$Unit), ")", sep = "")) +
            ylab("Inhibition") +
            labs(title = paste("Dose-Response curve:", i, sep = " ")) +
            ggplot2::annotate(x=sort(unique(.drm.data$Dose))[2], xend=sort(unique(.drm.data$Dose))[5], y=min(.drm.data$Inhibition)-0.05, yend=min(.drm.data$Inhibition)-0.05, 
                     colour="black", lwd=0.5, geom="segment") +
            ggplot2::coord_trans(x="log",ylim = c(min(.drm.data$Inhibition)-0.05, 1), clip = "on") +
            theme(plot.title = element_text(color="black", size=12*.textscalefactor, face="bold", hjust = 0.0, margin = grid::unit(c(0,0,12,0), "mm")),
                  axis.title = element_text(color="black", size=12*.textscalefactor, face="plain", hjust = 0.5),
                  axis.title.x = element_text(margin = ggplot2::margin(t = 10, r = 0, b = 0, l = 0)),
                  axis.title.y = element_text(margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 0)),
                  legend.title = element_text(color="black", size=10*.textscalefactor, face="plain", hjust = 0.0),
                  legend.text = element_text(color="black", size=10*.textscalefactor, face="plain", hjust = 0.0),
                  plot.caption = element_text(size = 10),
                  text = element_text(size = 12*.textscalefactor),
                  axis.line = element_line(color='black', size=0),
                  panel.grid.major = element_line(color='lightgrey', linetype = "dotted"),
                  panel.grid.minor = element_line(color='lightgrey', linetype = "dotted"),
                  panel.background = element_rect(fill=alpha('white', 2)),
                  panel.border = element_rect(colour = "black", fill=NA, size=0),
                  plot.margin = grid::unit(c(2,2,2,2), "mm"),
                  legend.position = "right"
            )
          
          rm(.drm.data, .drm.model, .drm.prediction, i)
          
        }
        
        
        # PLOTTING DRUG RESPONSE MATRIX
        .dose.response.matrix <- synDRM[[samplename]][[drugname]][[drugpair]][["dfs"]]
        
        .textscalefactor = 1.45
        
        synPlots[["dose.response.matrix"]][[drugname]][[drugpair]] <- ggplot2::ggplot(.dose.response.matrix, aes(x = factor(sort(as.numeric(format(round(ConcCol, 4), nsmall = 4)))), y = factor(as.numeric(format(round(ConcRow, 4), nsmall = 4))))) + 
          geom_raster(aes(fill = Response), na.rm=TRUE) +
          #scale_fill_hue("Value", l=80, h = c(270, 360), na.value = NA) +
          ggplot2::scale_fill_gradient2(low = "#00CC66", mid = "#FFFFFF", high = "#FF0066",
                               # low = "#EDC9AF", mid = "#FFFFFF", high = "#4B573E",
                               space = "Lab", na.value = "grey95", guide = "colourbar",
                               aesthetics = "fill", limits = c(-1,1)) +
          ggplot2::geom_text(aes(label = format(round(Response, 2), nsmall = 2), color = ifelse(is.na(Response), TRUE, FALSE)), hjust = 0.5, size = 3*.textscalefactor) +
          ggplot2::geom_hline(yintercept=seq(1.5, length(unique(.dose.response.matrix$Row))-0.5, 1), lwd=0.1, colour="black") +
          ggplot2::geom_vline(xintercept=seq(1.5, length(unique(.dose.response.matrix$Column))-0.5, 1), lwd=0.1, colour="black") +
          ggplot2::scale_x_discrete(position = "bottom") +
          ggplot2::scale_color_manual(values = c("black", "grey75"), guide = FALSE) +
          xlab(.dose.response.matrix$DrugCol) +
          ylab(.dose.response.matrix$DrugRow) +
          labs(title = "Dose-Response Matrix", fill = "Inhibition") +
          theme(plot.title = element_text(color="black", size=12*.textscalefactor, face="bold", hjust = 0.0),
                axis.title = element_text(color="black", size=12*.textscalefactor, face="plain", hjust = 0.5),
                axis.title.x = element_text(margin = ggplot2::margin(t = 10, r = 0, b = 0, l = 0)),
                axis.title.y = element_text(margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 0)),
                legend.title = element_text(color="black", size=10*.textscalefactor, face="plain", hjust = 0.0),
                legend.text = element_text(color="black", size=10*.textscalefactor, face="plain", hjust = 0.0),
                legend.text.align = 1,
                plot.caption = element_text(size = 10*.textscalefactor),
                text = element_text(size = 12*.textscalefactor),
                axis.line = element_line(color='black'),
                panel.grid.major = element_line(color='lightgrey', linetype = "dotted"),
                panel.grid.minor = element_line(color='lightgrey', linetype = "dotted"),
                panel.background = element_rect(fill=alpha('white', 1)),
                panel.border = element_rect(colour = "black", fill=NA, size=1),
                legend.position = "right"
          )
        
        
        g1 <- gridExtra::grid.arrange(synPlots[["dose.response.curve"]][[drugpair]], ggplot()+theme(panel.background = element_rect(fill=alpha('white', 1))), synPlots[["dose.response.curve"]][[drugname]], nrow = 3, heights=c(1,0.1,1))
        g <- gridExtra::grid.arrange(g1, synPlots[["dose.response.matrix"]][[drugname]][[drugpair]], nrow = 1)
        
        ggplot2::ggsave(filename = file.path(.saveto, "results", "graphs/custom plots/dose-response", samplename, gsub("/", " ", drugname), paste(samplename, " ", gsub("/", " ", drugname), " + ", gsub("/", " ", drugpair), ".png", sep = "")),
               g, device = "png", width = 245*2, height = 245, units = "mm", dpi = 300, limitsize = FALSE)
        
        rm(.dose.response.matrix, g, g1)
        
      }
      
    }
    
    
    cat('\r', "Finished plotting dose-response curves for ", drugname, ".", strrep(" ", 100), '\n', sep = "")
    
    if(samplename == utils::tail(names(synDRM), n = 1L)){
      message('\n', "Dose-response plots saved to: ", '\n', file.path(.saveto, "results", "graphs/custom plots/dose-response", samplename, gsub("/", " ", drugname)), strrep(" ", 100), sep = "")
      rm(samplename, drugname, drugpair)
    }
    
  }

  
}
