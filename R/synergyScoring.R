#' Scoring drug interactions
#'
#' @description
#' An essential component of the modular library \pkg{metascreen}, and imperative for the analysis of the experimental data from a drug sensitivity screen.
#' \emph{\code{synergyScoring}} is a function that scores drug interactions and ranks them based on their individual synergy scores.
#' 
#' @param bayesdata an object of class 'bayesdata'.
#' @param .saveoutput logical; if TRUE, then the data will be saved to file. The default is FALSE, in which the output is not saved.
#' @param .plot logical; if TRUE, then drug interactions scores will be plotted. The default is FALSE, in which no plots will be generated.
#' @param .saveto string; path to a folder location where the results are saved to.
#' 
#' 
#' @details The function \code{synergyScoring} is used to score and rank drug interactions between drugs of a combination treatment.
#' if plotting is set to TRUE, the function will generate three types of plots: (1) one in which all drug synergy and antagonism scores are plotted for each sample, (2) another one 
#' in which only the synergy score is plotted by drug fore each sample and (3) one in which only synergy and antagonism scores are plotted for each individual drug.
#' 
#' Files are saved either to the specified location or the default working environment, with the corresponding folder structure: 'results/synergyscores/bayesynergy' and 'results/graphs/synergyscores/bayesynergy'
#' 
#' 
#' @return Returns an class S3 object of type list with ranked drug combinations and additional statistics.
#' 
#' 
#' @examples
#' \donttest{\dontrun{
#' # Score drug interactions and save both the output and plots
#' synergyScoring(bayesdata, .saveoutput = TRUE, .plot = TRUE, .saveto = "path/to/folder/")
#' 
#' # Run synergy scoring without generating plots or saving the output
#' synergyScoring(bayesdata)
#' }}
#'
#' @keywords drug screen analysis dose response curve matrix
#' 
#' @importFrom stats mad sd
#' @importFrom ggplot2 annotate element_blank geom_abline ggplot ggsave guide_legend guides margin scale_alpha scale_colour_gradient2 scale_fill_gradient2 scale_size_continuous scale_x_continuous
#' 
#' @export

synergyScoring <- function(bayesdata, .saveoutput, .plot, .saveto){
  
  # Check, if the data has been provided as an object of class S3:bayesdata
  if(missing(bayesdata)){stop("Data missing! Please provide a data set as an object of class 'bayesdata' ", call. = TRUE)}
  if(class(bayesdata) != "bayesdata"){
    stop("Provided data not of class 'bayesdata'!", call. = TRUE)
  } else {
    bayDFS <- bayesdata
    synergyMethod <- "bayesynergy"
  }
  
  
  # Check, if output should be saved
  if(missing(.saveoutput)){ .saveoutput <- FALSE } 
  else if(!is.logical(.saveoutput)){
    stop("in '.saveoutput'. Argument needs to be a logical value, either TRUE or FALSE.", call. = TRUE)
  }
  
  # Check, if plots should be generated
  if(missing(.plot)){ .plot <- FALSE } 
  else if(!is.logical(.plot)){
    stop("in '.plot'. Argument needs to be a logical value, either TRUE or FALSE.", call. = TRUE)
  }
  
  # Check, if a folder location has been provided
  if(missing(.saveto)){ .saveto <- getwd() }
  # Create folder, if it does not exist
  if(!file.exists(file.path(.saveto))){ dir.create(file.path(.saveto), showWarnings = FALSE, recursive = TRUE) }
  if(!utils::file_test("-d", .saveto)){ stop("in '.saveto'. Argument needs to be a valid folder location.", call. = TRUE) }
  
  
  
  cat("Ranking drug combinations based on their synergy score.", '\n', sep = "")
  
  
  # SECTION A: Ranking Synergies ##################################################
  # Extracting the highest and lowest synergy scores for each drug combination
  synMetrics = list()
  
  for(samplename in names(bayDFS)){
    
    # Create subfolder, if it does not exist
    if(!file.exists(file.path(.saveto, "results", "synergyscores/bayesynergy", samplename))){ dir.create(file.path(.saveto, "results", "synergyscores/bayesynergy", samplename), showWarnings = FALSE, recursive = TRUE)}
    
    cat("Ranking synergies for ", samplename, ".", sep = "")
    
    
    for(drugname in names(bayDFS[[samplename]][["bayesynergy"]])){
      
      cat('\r', "Extracting synergy score for", samplename, ":", drugname, strrep(" ", 100), sep = " ")
      
      for(drugpair in names(bayDFS[[samplename]][["bayesynergy"]][[drugname]])){
        
        # Extract the relative data from the bayesynergy object
        .VUS <- bayDFS[[samplename]][["bayesynergy"]][[drugname]][[drugpair]][["VUS"]]
        .ec50 <- bayDFS[[samplename]][["bayesynergy"]][[drugname]][[drugpair]][["ec50"]]
        .summary <- bayDFS[[samplename]][["bayesynergy"]][[drugname]][[drugpair]][["summary"]]
        .quality <- bayDFS[[samplename]][["bayesynergy"]][[drugname]][[drugpair]][["quality"]]
        .bayesfactor <- bayDFS[[samplename]][["bayesynergy"]][[drugname]][[drugpair]][["bayesfactor"]]
        
        # Extract the volume under the surface (VUS) for synergy and antagonism
        .synergy <- .VUS$VUS_syn * -1
        .antagonism <- .VUS$VUS_ant
        
        
        synMetrics[["bayesynergy"]][[samplename]] <- rbind(synMetrics[["bayesynergy"]][[samplename]], data.frame(
          `Drug A` = drugname,
          `Drug B` = drugpair,
          # Extract the EC50s for each drug
          `EC50 (Drug A)` = mean(.ec50$ec50_1),
          `EC50 (Drug B)` = mean(.ec50$ec50_2),
          
          # Calculate a standardized synergy score from the mean of the VUS
          `Synergy Score` = median(.synergy) / mad(.synergy),
          
          # Calculating additional statistics and parameters
          `Median (syn)` = median(.synergy),
          `MAD (syn)` = mad(.synergy),
          `Mean (syn)` = .summary["VUS_syn", "mean"] * -1,
          `SEM (syn)` = .summary["VUS_syn", "se_mean"],
          `SD (syn)` = .summary["VUS_syn", "sd"],
          `Median/MAD (syn)` = median(.synergy) / mad(.synergy),
          `Mean/SD (syn)` =  mean(.synergy) / sd(.synergy),
          `Maximum (syn)` = max(.synergy),
          
          `Antagonism Score` = median(.antagonism) / mad(.antagonism),
          `Median (ant)` = median(.antagonism),
          `MAD (ant)` = mad(.antagonism),
          `Mean (ant)` = .summary["VUS_ant", "mean"],
          `SEM (ant)` = .summary["VUS_ant", "se_mean"],
          `SD (ant)` = .summary["VUS_ant", "sd"],
          `Median/MAD (ant)` = median(.antagonism) / mad(.antagonism),
          `Mean/SD (ant)` = mean(.antagonism) / sd(.antagonism),
          `Maximum (ant)` = max(.antagonism),
          
          # Calculate QC parameters
          `bayes factor` = .bayesfactor,
          `observation noise` = .quality[["s"]],
          `kernel length scale` = .quality[["ell"]],
          `kernel amplitude` = .quality[["sigma_f"]],
          `divergent` = .quality[["divergent"]],
          
          
          check.names = FALSE, stringsAsFactors = FALSE))
        
        rm(.VUS, .ec50, .summary, .synergy, .antagonism)
        
      }
      
    }
    
    
    # Annotate data set with additional information
    synMetrics[["bayesynergy"]][[samplename]]$`Cell Line` <- samplename
    synMetrics[["bayesynergy"]][[samplename]]$`Package` <- "bayesynergy"
    synMetrics[["bayesynergy"]][[samplename]]$`Version` <- packageVersion("bayesynergy")
    synMetrics[["bayesynergy"]][[samplename]]$`Model` <- "Bayesian"
    
    # Remove duplicate datasets and sort by maximum score for which no antagonism is found and return the top 10 entries
    synMetrics[["bayesynergy"]][[samplename]] <- synMetrics[["bayesynergy"]][[samplename]][!duplicated(t(apply(synMetrics[["bayesynergy"]][[samplename]][c("Drug A", "Drug B")], 1, sort))),]
    
    # Rank drug combinations by synergy score
    synMetrics[["bayesynergy"]][[samplename]] <- synMetrics[["bayesynergy"]][[samplename]][with(synMetrics[["bayesynergy"]][[samplename]], order(`Synergy Score`, decreasing = TRUE)),]
    synMetrics[["bayesynergy"]][[samplename]]$Rank <- 1:nrow(synMetrics[["bayesynergy"]][[samplename]])
    
    # Add a column indicating, if the score has passed or failed the qc-check
    synMetrics[["bayesynergy"]][[samplename]]$QC <- ifelse(synMetrics[["bayesynergy"]][[samplename]][["divergent"]] > 100, "FAIL", "PASS")
    
    # Split data set into pass and fail, and combine both back together to have the failed scores in the end
    synMetrics[["bayesynergy"]][[samplename]] <- with(synMetrics[["bayesynergy"]][[samplename]], split(synMetrics[["bayesynergy"]][[samplename]] , divergent > 100))
    synMetrics[["bayesynergy"]][[samplename]] <- do.call(rbind, setNames(synMetrics[["bayesynergy"]][[samplename]], NULL))
    
    synMetrics[["bayesynergy"]][[samplename]] <- synMetrics[["bayesynergy"]][[samplename]][c("Cell Line", "Package", "Version", "Model", "Rank", "QC", "Drug A", "Drug B", 
                                                                                             "Synergy Score", "Median (syn)", "MAD (syn)", "Mean (syn)", "SEM (syn)", "SD (syn)", "Median/MAD (syn)", "Mean/SD (syn)", "Maximum (syn)",
                                                                                             "Antagonism Score", "Median (ant)", "MAD (ant)", "Mean (ant)", "SEM (ant)", "SD (ant)", "Median/MAD (ant)", "Mean/SD (ant)", "Maximum (ant)",
                                                                                             "EC50 (Drug A)", "EC50 (Drug B)", "bayes factor", "observation noise", "kernel length scale", "kernel amplitude", "divergent")]
      


    
    # Save  output to file
    if(isTRUE(.saveoutput)){
      write.csv2(synMetrics[["bayesynergy"]][[samplename]], file = file.path(.saveto, "results", "synergyscores/bayesynergy", samplename, paste(samplename, "synergyscores.csv", sep = "")), row.names = FALSE, quote = FALSE)
    }
    
    
    if(samplename == utils::tail(names(bayDFS), n = 1L)){
      cat('\r', "Finished ranking synergies.", strrep(" ", 100), '\n\n', sep = "")
      rm(samplename, drugname, drugpair)
    }

  }
  
  if(isTRUE(.saveoutput)){
    message("Ranked synergy scores saved to: ", '\n', file.path(.saveto, "results", "synergyscores/bayesynergy"), strrep(" ", 100), sep = "")
  }
  
  

  # Plot synergy scores and rankings
  if(isTRUE(.plot)){
    
    cat("Plotting synergy ~ antagonism scores.", '\n', sep = "")
    
    # SECTION B: Creating plots visualizing the synergy scores for all the samples ############################
    # Plot the highest synergy and antagonism score for all drugs in a single plot
    for(samplename in names(bayDFS)){
      
      
      # Creating a function to mark top synergy scores
      in_quantile <- function(x) {
        return(x > quantile(min(x):max(x), 0.75))
        # return(x > 1.0)
      }
      
      cat('\r', "Plotting synergy scores for ", samplename, ".", strrep(" ", 100), sep = "")
      
      
      for(drugname in names(bayDFS[[samplename]][["bayesynergy"]])){
        
        cat('\r', " > Fetching synergy scores for ", samplename, ": ", drugname, strrep(" ", 100), sep = "")
        
        for(drugpair in names(bayDFS[[samplename]][["bayesynergy"]][[drugname]])){
          
          # QC-check, exclude scores with a divergent > 100
          if(bayDFS[[samplename]][["bayesynergy"]][[drugname]][[drugpair]][["quality"]][["divergent"]] > 100){ next() }
          
          
          # Check if data frame exists. Create it if not.
          if(!exists(".synscores", inherits = FALSE)){ .synscores = data.frame() }
          
          # Load and extract the volume under the surface (VUS) to
          .VUS <- bayDFS[[samplename]][["bayesynergy"]][[drugname]][[drugpair]][["VUS"]]
          
          .synergy <- .VUS$VUS_syn * -1
          .antagonism <- .VUS$VUS_ant
          
          .synscores <- rbind(.synscores, data.frame(
            Sample = samplename,
            Drug.A = drugname,
            Drug.B = drugpair,
            Synergy.Score = median(.synergy) / mad(.synergy),
            Antagonism.Score = median(.antagonism) / mad(.antagonism),
            
            check.names = FALSE, stringsAsFactors = FALSE))
          
          rm(.VUS, .synergy, .antagonism)
          
        }
        
      }
      
      # Quantile determination: label top synergy and antagonism hits
      .synscores <- transform(.synscores, Quantile = ifelse(in_quantile(as.numeric(Synergy.Score))|in_quantile(as.numeric(Antagonism.Score)),
                                                              paste(Drug.A, "+", Drug.B, sep = " "), as.character(NA)))
      
      # Remove duplicate quantile annotations. Drug combinations are repeating itself (transposed), but only the annotation for one drug pair is kept.
      .synscores <- transform(.synscores, Quantile = Quantile[!duplicated(t(apply(.synscores[c("Drug.A", "Drug.B")], 1, sort)))])
      .synscores <- .synscores[!duplicated(t(apply(.synscores[c("Drug.A", "Drug.B")], 1, sort))),]
      
      
      cat('\r', " > Plotting synergy scores for ", samplename, ".", strrep(" ", 100), sep = "")
      
      ggplot2::ggplot(.synscores, aes(x = Antagonism.Score, y = Synergy.Score)) +
        annotate(geom = "rect", xmin = -Inf, xmax = 1, ymin = -Inf, ymax = 1, fill = "#FCFCFC", size = 0.5, color="#FFFFFF", alpha = 0.45) +
        annotate(geom = "text", x = -Inf, y = 1, label = "  synergies > 1.0", color = "#0A0A0A", size = 4.5, hjust = 0, vjust = 1) + 
        annotate(geom = "text", x = 1, y = -Inf, label = "antagonsim > 1.0  ", color = "#0A0A0A", size = 4.5, hjust = 1, vjust = 1, angle = -90) + 
        geom_point(aes(colour = Drug.A, size = abs(ifelse(abs(Antagonism.Score) > abs(Synergy.Score), Antagonism.Score, Synergy.Score)), alpha = 
                         ifelse(abs(Antagonism.Score) <= 1 & abs(Synergy.Score) <= 1, 0.4, 0.8) )) +
        
        # geom_point(aes(colour = `Drug B`), size = 1) +
        scale_x_continuous(labels = as.character(floor(min(.synscores$Antagonism.Score, .synscores$Synergy.Score)):ceiling(max(.synscores$Antagonism.Score, .synscores$Synergy.Score))),
                           breaks = floor(min(.synscores$Antagonism.Score, .synscores$Synergy.Score)):ceiling(max(.synscores$Antagonism.Score, .synscores$Synergy.Score))) +
        scale_size_continuous(range = c(0.2, 8), guide = "none") +
        scale_alpha(guide = "none") +
        
        # geom_abline(intercept = 0, slope = 2.5, size = 0.5) +
        # scale_size(range = c(0.1, 10)) +
        # scale_x_reverse(name = "Antagonism") +
        # scale_y_continuous(name = "Synergy", breaks = seq(0, 100, 25), labels = c("0", "25", "50", "75", "100")) +
        # scale_x_reverse(name = "Antagonism", limits = c(0, -100)) +
        # scale_y_continuous(name = "Synergy", limits = c(0, 100)) +
        # geom_hline(aes(yintercept=0.5), lwd=0.2, colour="red") +
        # geom_text(aes(label = quantile), na.rm = TRUE, hjust = -0.2) +
        # geom_text_repel(aes(label = quantile), na.rm = TRUE, hjust = -0.2) +
        geom_text_repel(aes(label = Quantile), na.rm = TRUE, force = 1, direction = "both", point.padding = unit(1.2, "lines"), box.padding = unit(0.5, "lines"), 
                        segment.size = 0.2, segment.color = "black", nudge_x = 0.01, nudge_y = 0, hjust = 0.5, size = 5) +
        facet_grid(Sample~.) +
        xlab("Antagonism") +
        ylab("Synergy") +
        ggtitle("Drug synergy vs. antagonism scores") +
        labs(tag = paste("Package: ", synergyMethod, " (v", packageVersion("bayesynergy"), ")", sep = "")) +
        guides(colour = guide_legend(title="Drugs", nrow = 4), size = "none") +
        theme(plot.title = element_text(hjust = 0.5), legend.title = element_blank(), legend.text = element_text(size = 10),
              legend.position = "bottom", legend.justification="center", legend.direction = "horizontal", legend.box.margin=margin(20,0,0,0),
              text = element_text(size = 20), plot.caption = element_text(size = 12), plot.caption.position = "plot",
              plot.tag = element_text(hjust = 1, size = 18), plot.tag.position = c(1, 0.99),
              axis.text = element_text(angle = 45, vjust = 1, hjust = 1))
      
      
      cat('\r', " > Saving synergy plot for ", samplename, ".", strrep(" ", 100), sep = "")
      
      # Create subfolder, if it does not exist
      if(!file.exists(file.path(.saveto, "results", "graphs/synergyscores", synergyMethod, "allsynergies", samplename))){ dir.create(file.path(.saveto, "results", "graphs/synergyscores", synergyMethod, "allsynergies", samplename), showWarnings = FALSE, recursive = TRUE)}
      
      ggsave(filename = file.path(.saveto, "results", "graphs/synergyscores", synergyMethod, "allsynergies", samplename, paste(samplename, "synergyScores.png", sep = " ")), device = "png", width = 804, height = 512, units = "mm", dpi = 300, limitsize = FALSE)
      
      
      if(samplename == utils::tail(names(bayDFS), n = 1L)){
        cat('\r', "Finished plotting synergy ~ antagonism scores.", strrep(" ", 100), '\n\n', sep = "")
        rm(.synscores, samplename, drugname, drugpair, in_quantile)
      }

    }
    
    
    
    cat("Plotting synergy scores by drug.", '\n', sep = "")
    
    # SECTION C: Creating plots visualizing the synergy scores for all the samples and drugs
    # Print the highest synergy score for each drug combination by each individual drug
    
    for(samplename in names(bayDFS)){
      

      # Creating a function to mark top synergy scores
      in_quantile <- function(x) {
        return(x > quantile(min(x):max(x), 0.75))
        # return(x > 1.0)
      }
      
      cat('\r', "Plotting synergy scores for ", samplename, ".", strrep(" ", 100), sep = "")
      
      
      for(drugname in names(bayDFS[[samplename]][["bayesynergy"]])){
        
        cat('\r', " > Fetching synergy scores for ", samplename, ": ", drugname, strrep(" ", 100), sep = "")
        
        for(drugpair in names(bayDFS[[samplename]][["bayesynergy"]][[drugname]])){
          
          # QC-check, exclude scores with a divergent > 100
          if(bayDFS[[samplename]][["bayesynergy"]][[drugname]][[drugpair]][["quality"]][["divergent"]] > 100){ next() }
          
          
          # Check if data frame exists. Create it if not.
          if(!exists(".synscores", inherits = FALSE)){ .synscores = data.frame() }
          
          # Load and extract the volume under the surface (VUS) to
          .VUS <- bayDFS[[samplename]][["bayesynergy"]][[drugname]][[drugpair]][["VUS"]]
          
          .synergy <- .VUS$VUS_syn * -1
          .antagonism <- .VUS$VUS_ant
          
          .synscores <- rbind(.synscores, data.frame(
            Sample = samplename,
            Drug.A = drugname,
            Drug.B = drugpair,
            Synergy.Score = median(.synergy) / mad(.synergy),
            Antagonism.Score = median(.antagonism) / mad(.antagonism),
            
            check.names = FALSE, stringsAsFactors = FALSE))
          
          rm(.VUS, .synergy, .antagonism)
          
        }
        
      }
      
      # Quantile determination: label top synergy and antagonism hits
      .synscores <- transform(.synscores, Quantile = ifelse(in_quantile(as.numeric(Synergy.Score)), paste(Drug.A, "+", Drug.B, sep = " "), as.character(NA)))
      
      
      ggplot2::ggplot(.synscores, aes(x = factor(Drug.A), y = Synergy.Score)) +
        annotate(geom = "text", x = -Inf, y = 1, label = "  synergies > 1.0", color = "#0A0A0A", size = 4.5, hjust = 0, vjust = 1) + 
        geom_abline(intercept = 1, slope = 0, size = 1.2, linetype = "dotted", color = "#A0A0A0", alpha = 0.8) +
        geom_point(aes(colour = Drug.B, size = Synergy.Score)) +
        
        # geom_point(aes(colour = `Drug B`), size = 1) +
        scale_size_continuous(range = c(0.2, 8), guide = "none") +
        
        # geom_abline(intercept = 74, slope = 0, size = 0.5) +
        # scale_size(range = c(0.1, 10)) +
        # scale_x_reverse(name = "Antagonism") +
        scale_y_continuous(name = "Synergy") +
        # scale_x_reverse(name = "Antagonism", limits = c(0, -100)) +
        # scale_y_continuous(name = "Synergy", limits = c(0, 100)) +
        # geom_hline(aes(yintercept=0.5), lwd=0.2, colour="red") +
        # geom_text(aes(label = quantile), na.rm = TRUE, hjust = -0.2) +
        # geom_text_repel(aes(label = quantile), na.rm = TRUE, hjust = -0.2) +
        geom_text_repel(aes(label = Quantile), na.rm = TRUE, force = 1, direction = "both", point.padding = unit(1.2, "lines"), box.padding = unit(0.5, "lines"), 
                        segment.size = 0.2, segment.color = "black", nudge_x = 0.02, nudge_y = 0, hjust = 0.5, size = 5) +
        facet_grid(Sample~.) +
        xlab("") +
        ylab("Synergy") +
        ggtitle("Synergy scores") +
        labs(tag = paste("Package: ", synergyMethod, " (v", packageVersion("bayesynergy"), ")", sep = "")) +
        guides(colour = guide_legend(title="Drugs", nrow = 4), size = "none") +
        theme(plot.title = element_text(hjust = 0.5), legend.title = element_blank(), legend.text = element_text(size = 10),
              legend.position = "bottom", legend.justification="center", legend.direction = "horizontal", legend.box.margin=margin(20,0,0,0),
              text = element_text(size = 20),  plot.caption = element_text(size = 12), plot.caption.position = "plot",
              plot.tag = element_text(hjust = 1, size = 18), plot.tag.position = c(1, 0.99),
              axis.text.x = element_text(size = 10, angle = 45, vjust = 1, hjust = 1))
      
      
      cat('\r', " > Saving synergy plot for ", samplename, ".", strrep(" ", 100), sep = "")
      
      # Create subfolder, if it does not exist
      if(!file.exists(file.path(.saveto, "results", "graphs/synergyscores", synergyMethod, "synergiesbydrug", samplename))){ dir.create(file.path(.saveto, "results", "graphs/synergyscores", synergyMethod, "synergiesbydrug", samplename), showWarnings = FALSE, recursive = TRUE)}
      
      
      ggsave(filename = file.path(.saveto, "results", "graphs/synergyscores", synergyMethod, "synergiesbydrug", samplename, paste(samplename, "synergyScores-bydrug.png", sep = " ")), device = "png", width = 804, height = 512, units = "mm", dpi = 300, limitsize = FALSE)
      
      
      if(samplename == utils::tail(names(bayDFS), n = 1L)){
        cat('\r', "Finished plotting synergy scores by drug.", strrep(" ", 100), '\n\n', sep = "")
        rm(.synscores, samplename, drugname, drugpair, in_quantile)
      }
      
    }
    
    
    
    cat("Plotting synergy ~ antagonsim scores for each individual drug.", '\n', sep = "")
  
    # SECTION D: Creating plots visualizing the synergy scores for all the samples and drugs
    # Print the highest synergy and antagonism score for each individual drug
    
    for(samplename in names(bayDFS)){
      
      # Creating a function to mark top synergy scores
      in_quantile <- function(x) {
        return(x > quantile(min(x):max(x), 0.75))
        # return(x > 1.0)
      }
      
      cat('\r', "Plotting synergy scores for ", samplename, ".", strrep(" ", 100), sep = "")
      
      
      for(drugname in names(bayDFS[[samplename]][["bayesynergy"]])){
        
        cat('\r', " > Fetching synergy scores for ", samplename, ": ", drugname, strrep(" ", 100), sep = "")
        
        for(drugpair in names(bayDFS[[samplename]][["bayesynergy"]][[drugname]])){
          
          # QC-check, exclude scores with a divergent > 100
          if(bayDFS[[samplename]][["bayesynergy"]][[drugname]][[drugpair]][["quality"]][["divergent"]] > 100){ next() }
          
          
          # Check if data frame exists. Create it if not.
          if(!exists(".synscores", inherits = FALSE)){ .synscores = data.frame() }
          
          # Load and extract the volume under the surface (VUS) to
          .VUS <- bayDFS[[samplename]][["bayesynergy"]][[drugname]][[drugpair]][["VUS"]]
          
          .synergy <- .VUS$VUS_syn * -1
          .antagonism <- .VUS$VUS_ant
          
          .synscores <- rbind(.synscores, data.frame(
            Sample = samplename,
            Drug.A = drugname,
            Drug.B = drugpair,
            Synergy.Score = median(.synergy) / mad(.synergy),
            Antagonism.Score = median(.antagonism) / mad(.antagonism),
            
            check.names = FALSE, stringsAsFactors = FALSE))
          
          rm(.VUS, .synergy, .antagonism)
          
        }
        
        # Quantile determination: label top synergy and antagonism hits
        .synscores <- transform(.synscores, Quantile = ifelse(in_quantile(as.numeric(Synergy.Score))|in_quantile(as.numeric(Antagonism.Score)),
                                                              paste(Drug.A, "+", Drug.B, sep = " "), as.character(NA)))
        
        # Set the limits for the x- and y-axis
        yScaleLimit <- c(floor(min(.synscores$Synerg.Score, .synscores$Antagonism.Score)), ceiling(max(.synscores$Synergy.Score, .synscores$Antagonism.Score)))
        
        ggplot2::ggplot(.synscores, aes(x = Antagonism.Score, y = Synergy.Score)) +
          annotate(geom = "rect", xmin = -Inf, xmax = 1, ymin = -Inf, ymax = 1, fill = "#FCFCFC", size = 0.5, color="#FFFFFF", alpha = 0.25) +
          annotate(geom = "text", x = -Inf, y = 1, label = "  synergies > 1.0", color = "#0A0A0A", size = 4.5, hjust = 0, vjust = 1) + 
          annotate(geom = "text", x = 1, y = -Inf, label = "antagonsim > 1.0  ", color = "#0A0A0A", size = 4.5, hjust = 1, vjust = 1, angle = -90) + 
          geom_abline(intercept = 0, slope = 1, size = 1.2, linetype = "dotted", color = "#CCCCCC", alpha = 0.8) +
          geom_point(aes(size = abs(ifelse(abs(Antagonism.Score) > abs(Synergy.Score), Antagonism.Score, Synergy.Score)),
                         fill = ifelse(abs(Antagonism.Score) > abs(Synergy.Score), Antagonism.Score*-1, Synergy.Score),
                         colour = ifelse(abs(Antagonism.Score) > abs(Synergy.Score), Antagonism.Score*-1, Synergy.Score)), stroke = 1, shape = 21, alpha = 1.0) +
          scale_size_continuous(range = c(0.2, 12)) +
          scale_fill_gradient2(low = "#FF0000", high = "#00CC00") +
          scale_colour_gradient2(low = "#660000", mid = "#AAAAAA", high = "#006600") +
          scale_x_continuous(name = "Antagonism", limits = yScaleLimit) +
          scale_y_continuous(name = "Synergy", limits = yScaleLimit) +
          # geom_text(aes(label = Quantile), na.rm = TRUE, size = 5, hjust = -0.3) +
          geom_text_repel(aes(label = Quantile), na.rm = TRUE, force = 1, direction = "both", point.padding = unit(1.2, "lines"), box.padding = unit(0.5, "lines"),
                          segment.size = 0, segment.color = "black", nudge_x = 0.5, nudge_y = 0, hjust = 0.5, size = 5) +
          facet_grid(Sample~.) +
          xlab("Antagonism") +
          ylab("Synergy") +
          ggtitle(paste(drugname, " Synergy & Antagonism Scores", sep = "")) +
          labs(tag = paste("Package: ", synergyMethod, " (v", packageVersion("bayesynergy"), ")", sep = "")) +
          theme(plot.title = element_text(size = 20, hjust = 0.5), legend.title = element_blank(), legend.text = element_text(size = 10),
                legend.position = "none", legend.justification="center", legend.direction = "horizontal", legend.box.margin=margin(20,0,0,0),
                text = element_text(size = 20),  plot.caption = element_text(size = 12), plot.caption.position = "plot",
                plot.tag = element_text(hjust = 1, size = 18), plot.tag.position = c(1, 0.99),
                axis.text = element_text(size = 12, angle = 45, vjust = 1, hjust = 1))
        
        
        cat('\r', " > Saving synergy plot for ", samplename, ": ", drugname, ".", strrep(" ", 100), sep = "")
        
        # Create subfolder, if it does not exist
        if(!file.exists(file.path(.saveto, "results", "graphs/synergyscores", synergyMethod, "individualsynergies", samplename, gsub("/", " ", drugname)))){ dir.create(file.path(.saveto, "results", "graphs/synergyscores", synergyMethod, "individualsynergies", samplename, gsub("/", " ", drugname)), showWarnings = FALSE, recursive = TRUE)}
        
        ggsave(filename = file.path(.saveto, "results", "graphs/synergyscores", synergyMethod, "individualsynergies", samplename, gsub("/", " ", drugname), paste(samplename, " ", gsub("/", " ", drugname), ".png", sep = "")), device = "png", width = 512, height = 512, units = "mm", dpi = 300, limitsize = FALSE)
        
        rm(.synscores, yScaleLimit)
        
      }
      
      if(samplename == utils::tail(names(bayDFS), n = 1L)){
        cat('\r', "Finished plotting synergy ~ antagonsim scores for each individual drug.", strrep(" ", 100), '\n\n', sep = "")
        rm(samplename, drugname, drugpair, in_quantile)
      }
      
    }
    
    
    message("Synergy plots saved to: ", '\n', file.path(.saveto, "results", "graphs/synergyscores", synergyMethod), strrep(" ", 100), sep = "")
    
  }
  
  
  # Return object of class synMetrics
  class(synMetrics) <- "synMetrics"
  return(synMetrics)
  
}
