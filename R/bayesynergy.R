#' Calculating synergies with bayesynergy
#'
#' @description
#' An essential component of the modular library \pkg{metascreen}, and imperative for the analysis of the experimental data from a drug sensitivity screen.
#' \emph{\code{bayesynergy}} is a function that calculates synergies between two drugs in a drug combination screen based on a bayesian semi-parametric model.
#' 
#' @param processedData an object of class 'processedData'.
#' @param .saveoutput logical; if TRUE, then the output from code{bayesynergy} will be saved to file. The default is FALSE, in which the output is not saved.
#' @param .plot logical; if TRUE, then plots from code{bayesynergy} will be saved to file. The default is FALSE, in which no plots will be generated.
#' @param .saveto string; path to a folder location where the results are saved to.
#' 
#' 
#' @details The function \code{bayesynergy} is used to assess the interaction between two drugs across their dose ranges based on a bayesian semi-parametric model.
#' The analysis is based on the package \code{\link[bayesynergy]{bayesynergy}} as described in the original paper by Rønneberg, 2021, Brief Bioinform (see references).
#' 
#' 
#' Files are saved either to the specified location or the default working environment, with the corresponding folder structure: 'results/data/bayesynergy' and 'results/graphs/bayesynergy'
#' 
#' 
#' @return Returns an class S3 object of type list with the volumetric surfaces (VUS), the EC50s, summary statistics and additional quality parameters, along with a measure of synergy (bayesfactor).
#' 
#' @references Rønneberg, L., Cremaschi, A., Hanes, R., Enserink, J. M., Zucknick, M. (2021) bayesynergy: flexible Bayesian modelling of synergistic interaction effects in in vitro drug combination experiments. Brief Bioinform 22 (6), bbab251, DOI: 10.1093/bib/bbab251
#' 
#' @seealso \code{\link[bayesynergy]{bayesynergy}}
#' 
#' @examples
#' \donttest{\dontrun{
#' # Run bayesnergy and save both the output and plots
#' bayesynergy(processedData, .saveoutput = TRUE, .plot = TRUE, .saveto = "path/to/folder/")
#' 
#' # Run bayesynergy without generating plots or saving the output
#' bayesynergy(processedData)
#' }}
#'
#' @keywords drug screen analysis dose response curve matrix
#' 
#' @importFrom bayesynergy bayesynergy
#' @importFrom rstan extract get_num_divergent summary
#' 
#' @export

bayesynergy <- function(processedData, .saveoutput, .plot, .saveto){
  
  # Check, if the data has been provided as an object of class S3:processedData
  if(missing(processedData)){stop("Data missing! Please provide a data set as an object of class 'processedData' ", call. = TRUE)}
  if(class(processedData) != "processedData"){
    stop("Provided data not of class 'processedData'!", call. = TRUE)
  } else {
    # Extract the required data from the class S3 object
    analysisData <- processedData[["analysisData"]]
    synDataset <- processedData[["splitDataset"]]
    synDRM <- processedData[["doserespMatrix"]]
  }
  
  
  # Check, if bayesynergy output should be saved
  if(missing(.saveoutput)){ .saveoutput <- FALSE } 
  else if(!is.logical(.saveoutput)){
    stop("in '.saveoutput'. Argument needs to be a logical value, either TRUE or FALSE.", call. = TRUE)
  }
  
  # Check, if bayesynergy plots should be generated
  if(missing(.plot)){ .plot <- FALSE } 
  else if(!is.logical(.plot)){
    stop("in '.plot'. Argument needs to be a logical value, either TRUE or FALSE.", call. = TRUE)
  }
  
  # Check, if a folder location has been provided
  if(missing(.saveto)){ .saveto <- getwd() }
  # Create folder, if it does not exist
  if(!file.exists(file.path(.saveto))){ dir.create(file.path(.saveto), showWarnings = FALSE, recursive = TRUE) }
  if(!utils::file_test("-d", .saveto)){ stop("in '.saveto'. Argument needs to be a valid folder location.", call. = TRUE) }

  
  
  cat("Running synergy analysis with bayesynergy ", paste(packageVersion("bayesynergy")), ".", '\n\n', sep = "")
  
  
  # SECTION A: REFORMAT DATA SET TO BAYESYNERGY REQUIREMENTS  ##############################
  # Merge individual data sets before dividing data sets based on individual drugs
  
  bayDFS = list()
  
  for (samplename in names(synDataset)){

    cat("Assembling data set for ", samplename, ".", '\n', sep = "")
    
    
    # Split the synergy data for each sample to create a dose response matrix (drm)
    bayDFS[[samplename]][["drm"]] <- split(synDataset[[samplename]][["combinationData"]], synDataset[[samplename]][["combinationData"]]$Drug)
    
    
    for(drugname in names(bayDFS[[samplename]][["drm"]])){
      
      cat('\r', " > Assembling data set for ", samplename, ": ", drugname, strrep(" ", 100), sep = "")
      
      
      # Add the corresponding drug pair to each drug
      bayDFS[[samplename]][["drm"]][[drugname]] <- merge(bayDFS[[samplename]][["drm"]][[drugname]], synDataset[[samplename]][["combinationData"]], 
                                                         by = c("Sample", "Destination.Plate.Barcode", "Plate.Number", "Combination.ID",
                                                                 "Destination.Well", "Viability"), suffixes = c(".1", ".2"))
      # Remove same drug pairs and select specific columns
      bayDFS[[samplename]][["drm"]][[drugname]] <- subset(bayDFS[[samplename]][["drm"]][[drugname]], Drug.1 != Drug.2)
      bayDFS[[samplename]][["drm"]][[drugname]] <- bayDFS[[samplename]][["drm"]][[drugname]][c("Drug.1", "Drug.2", "Drug.Concentration.1", "Drug.Concentration.2", "Unit.1", "Unit.2", "Viability")]
      
      # Rename selected columns
      names(bayDFS[[samplename]][["drm"]][[drugname]])[names(bayDFS[[samplename]][["drm"]][[drugname]]) %in%  c("Drug.1", "Drug.2", "Drug.Concentration.1", "Drug.Concentration.2", "Unit.1", "Unit.2", "Viability")] <- c("drugA", "drugB", "drugAconc", "drugBconc", "drugAunit", "drugBunit", "viability")

      
      # Split the dose response matrix (drm) for each drugpair
      bayDFS[[samplename]][["drm"]][[drugname]] <- split(bayDFS[[samplename]][["drm"]][[drugname]], bayDFS[[samplename]][["drm"]][[drugname]]$drugB)
      
      
      for(drugpair in names(bayDFS[[samplename]][["drm"]][[drugname]])){
        
        # Aggregating single drug treatments into the drug pair matrix
        bayDFS[[samplename]][["drm"]][[drugname]][[drugpair]] <- do.call("rbind", list(bayDFS[[samplename]][["drm"]][[drugname]][[drugpair]], 
                                                                                        transform(subset(synDataset[[samplename]][["singleDrugResponseData"]], Drug == drugname),
                                                                                          drugA = Drug, drugB = NA, drugAconc = Drug.Concentration, drugBconc = 0, drugAunit = Unit, drugBunit = Unit, viability = Viability)[c("drugA", "drugB", "drugAconc", "drugBconc", "drugAunit", "drugBunit", "viability")],
                                                                                        transform(subset(synDataset[[samplename]][["singleDrugResponseData"]], Drug == drugpair),
                                                                                                  drugA = NA, drugB = Drug, drugAconc = 0, drugBconc = Drug.Concentration, drugAunit = Unit, drugBunit = Unit, viability = Viability)[c("drugA", "drugB", "drugAconc", "drugBconc", "drugAunit", "drugBunit", "viability")]
                                                                                        ))
        
        # Re-formatting data to bayesynergy requirements
        # Select individual columns and drop ones that are not required
        # bayDFS[[samplename]][["drm"]][[drugname]][[drugpair]] <- bayDFS[[samplename]][["drm"]][[drugname]][[drugpair]][c("drugA", "drugB", "drugAconc", "drugBconc", "viability")]
        
        # Create indices for each replicate, which is then spread from long to wide-format
        bayDFS[[samplename]][["drm"]][[drugname]][[drugpair]]$id <- with(bayDFS[[samplename]][["drm"]][[drugname]][[drugpair]], ave(seq_along(paste(drugA, drugB, drugAconc, drugBconc)), paste(drugA, drugB, drugAconc, drugBconc), FUN=seq_along))
    
        # Replace, or convert NA to character in order to be able to reshape the data from long to wide-format
        bayDFS[[samplename]][["drm"]][[drugname]][[drugpair]][is.na(bayDFS[[samplename]][["drm"]][[drugname]][[drugpair]])] <- "NA"
        bayDFS[[samplename]][["drm"]][[drugname]][[drugpair]] <- stats::reshape(bayDFS[[samplename]][["drm"]][[drugname]][[drugpair]], idvar = c("drugA", "drugB", "drugAconc", "drugBconc"), v.names = "viability", timevar = "id", direction = "wide", new.row.names = NULL, sep = "")
        bayDFS[[samplename]][["drm"]][[drugname]][[drugpair]][bayDFS[[samplename]][["drm"]][[drugname]][[drugpair]] == "NA"] <- NA
        
        
        # Reorder rows based on drug names and drug doses
        bayDFS[[samplename]][["drm"]][[drugname]][[drugpair]] <- bayDFS[[samplename]][["drm"]][[drugname]][[drugpair]][with(bayDFS[[samplename]][["drm"]][[drugname]][[drugpair]], order(drugA, drugB, drugAconc, drugBconc)),]
        rownames(bayDFS[[samplename]][["drm"]][[drugname]][[drugpair]]) <- NULL
        

        
        # Convert the data set to individual matrices, in which matrix X contains the concentrations
        # and matrix y the viability values
        
        # Creating the X matrix
        bayDFS[[samplename]][["drm"]][[drugname]][[drugpair]][["X"]] <- as.matrix(bayDFS[[samplename]][["drm"]][[drugname]][[drugpair]][c( "drugAconc", "drugBconc")])
        colnames(bayDFS[[samplename]][["drm"]][[drugname]][[drugpair]][["X"]]) <- c(drugname, drugpair)
        
        # Creating the Y matrix
        bayDFS[[samplename]][["drm"]][[drugname]][[drugpair]][["Y"]] <- as.matrix(bayDFS[[samplename]][["drm"]][[drugname]][[drugpair]][grepl("viability", names(bayDFS[[samplename]][["drm"]][[drugname]][[drugpair]]))])
        
      }
      
    }
    
    if(samplename == utils::tail(names(synDataset), n = 1L)){
      cat('\r', "Finished assembling data set for ", samplename, ".", strrep(" ", 100), '\n\n', sep = "")
      rm(samplename, drugname, drugpair)
    }
    
  }
  
  
  
  # SECTION B: RUNNING ANALYSIS WITH  BAYESYNERGY ##############################
  # Calculating the synergy for each cell line and drug combination
  
  for (samplename in names(bayDFS)){
    
    # Skip samples that have already been analyzed
    # if(exists("bayesynergy", where = bayDFS[[samplename]])){next}
    
    cat("Analyzing synergies for ", samplename, ".", sep = "")
    
    
    for(drugname in names(bayDFS[[samplename]][["drm"]])){
      
      for(drugpair in names(bayDFS[[samplename]][["drm"]][[drugname]])){
        
        cat('\r', " > Analyzing synergies for ", samplename, ":", "[", match(drugname, names(bayDFS[[samplename]][["drm"]])), "]", drugname, ":", "[", match(drugpair, names(bayDFS[[samplename]][["drm"]])), "]", drugpair, strrep(" ", 50), sep = " ")
        
        
        # Run bayesynergy for synergy estimations using type 3 (GP with Matérn kernel)
        # Rønneberg, L., Cremaschi, A., Hanes, R., Enserink, J. M., Zucknick, M. (2021) bayesynergy: flexible Bayesian modelling of 
        # synergistic interaction effects in in vitro drug combination experiments. Brief Bioinform 22 (6), bbab251, DOI: 10.1093/bib/bbab251
        # package version: bayesynergy v2.4.1
        
        .bayoutput <- tryCatch(withCallingHandlers({
          suppressMessages(suppressWarnings(
            bayesynergy::bayesynergy(bayDFS[[samplename]][["drm"]][[drugname]][[drugpair]][["Y"]], 
                        bayDFS[[samplename]][["drm"]][[drugname]][[drugpair]][["X"]], type = 3,
                        drug_names = c(gsub("/", " ", drugname), gsub("/", " ", drugpair)), experiment_ID = samplename, units = c(unique(bayDFS[[samplename]][["drm"]][[drugname]][[drugpair]]$drugAunit), unique(bayDFS[[samplename]][["drm"]][[drugname]][[drugpair]]$drugBunit)),
                        lower_asymptotes = TRUE, heteroscedastic = TRUE, bayes_factor = TRUE, method = "sampling", refresh = 0)
          )) }) )
        

        # Re-run bayesynergy, if the divergent should be above a threshold of 100, with adapt_delta=0.99
        if(exists("divergent", where = .bayoutput)){
          if(all(.bayoutput$divergent > 100, !is.na(.bayoutput$divergent))){
            
            .bayoutput <- tryCatch(withCallingHandlers({
              suppressMessages(suppressWarnings(
                bayesynergy::bayesynergy(bayDFS[[samplename]][["drm"]][[drugname]][[drugpair]][["Y"]], 
                                         bayDFS[[samplename]][["drm"]][[drugname]][[drugpair]][["X"]], type = 3,
                                         drug_names = c(gsub("/", " ", drugname), gsub("/", " ", drugpair)), experiment_ID = samplename, units = c(unique(bayDFS[[samplename]][["drm"]][[drugname]][[drugpair]]$drugAunit), unique(bayDFS[[samplename]][["drm"]][[drugname]][[drugpair]]$drugBunit)),
                                         lower_asymptotes = TRUE, heteroscedastic = TRUE, bayes_factor = TRUE, method = "sampling", control = list(adapt_delta=0.99), refresh = 0)
              )) }) )
            
          }
        }

        
        
        # Save  output to file
        if(isTRUE(.saveoutput)){
          if(!file.exists(file.path(.saveto, "results", "data/bayesynergy", samplename, gsub("/", " ", drugname)))){ dir.create(file.path(.saveto, "results", "data/bayesynergy", samplename, gsub("/", " ", drugname)), showWarnings = FALSE, recursive = TRUE)}
          save(.bayoutput, file = file.path(.saveto, "results", "data/bayesynergy", samplename, gsub("/", " ", drugname), paste(gsub("/", " ", drugname), gsub("/", " ", drugpair), sep = " + ")))
        }
        
          
        # Extract data from the bayesynergy output
        .VUS <- rstan::extract(.bayoutput$stanfit, pars=c("VUS_syn","VUS_ant","VUS_Delta"))
        .ec50 <- rstan::extract(.bayoutput$stanfit, pars=c("ec50_1","ec50_2"))
        
        # Run summary function on bayesynergy output and extract additional quality control parameters
        .summary <- rstan::summary(.bayoutput$stanfit,pars=intersect(names(rstan::extract(.bayoutput$stanfit)),c("ell","sigma_f","s","ec50_1","ec50_2","VUS_Delta","VUS_syn","VUS_ant")),probs=c(0.025,.5,0.975))$summary
        .quality <- c(lapply(rstan::extract(.bayoutput$stanfit, pars=c("ell","sigma_f","s")), mean), divergent=rstan::get_num_divergent(.bayoutput$stanfit))
        
        # Extract the bayesfactor from the bayesynergy output
        .bayesfactor <- .bayoutput$bayesfactor
        
        

        
        # Generate bayesynergy plots and save them as png/pdf
        if(isTRUE(.plot)){
          if(!file.exists(file.path(.saveto, "results", "graphs/bayesynergy", samplename, gsub("/", " ", drugname), paste(gsub("/", " ", drugname), " + ", gsub("/", " ", drugpair), sep = "")))){ dir.create(file.path(.saveto, "results", "graphs/bayesynergy", samplename, gsub("/", " ", drugname), paste(gsub("/", " ", drugname), " + ", gsub("/", " ", drugpair), sep = "")), showWarnings = FALSE, recursive = TRUE) }
          suppressMessages(
            plot(.bayoutput, plot3D = FALSE, save_plot = TRUE, path = file.path(.saveto, "results", "graphs/bayesynergy", samplename, gsub("/", " ", drugname), paste(gsub("/", " ", drugname), " + ", gsub("/", " ", drugpair), "/", sep = "")),
                 plotdevice = "png", width = 178, height = 178, units = "mm", res = 300)
          )
        }
        
        
        # Assign data to an object
        .output <- list(VUS=.VUS, ec50=.ec50, summary=.summary, quality=.quality, bayesfactor=.bayesfactor)

        # Add drug names as attributes
        attr(.output, "drugname") <- drugname
        attr(.output, "drugpair") <- drugpair
        
        bayDFS[[samplename]][["bayesynergy"]][[drugname]][[drugpair]] <- .output
        
        rm(.bayoutput, .VUS, .ec50, .summary, .quality, .bayesfactor, .output)
        
      }
      
    }
    
    if(samplename == utils::tail(names(bayDFS), n = 1L)){
      cat('\r', "Finished analyzing synergies for ", samplename, ".", strrep(" ", 100), '\n\n', sep = "")
      rm(samplename, drugname, drugpair)
    }
    
  }
  
  
  # Return object of class bayesdata
  class(bayDFS) <- "bayesdata"
  return(bayDFS)
  
}

  