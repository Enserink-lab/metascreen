#' Combine drugs for a drug combination screen
#'
#' @description
#' An complementary component of the modular library \pkg{metascreen}, invaluable for the initial set-up of a drug sensitivity screen.
#' \emph{\code{combineDrugs}} is a function that combines drugs at selected doses, with each other. These drug combinations are used for the generation of a dispensing file, 
#' which is used to dispense drugs accordingly for a high-throughput drug combination screen. Drugs are combined either at all doses from a list or at defined ranges.
#' In addition, replicates can be generated either for all drug treatments, or for single drug treatments only.
#' 
#' @param listofDoses \code{data frame}; a list of drugs and doses in long-format. The data frame needs to include a single column of drug names, doses and units.
#' @param .combineDoses \code{numeric vector}; an individual selection or a range of doses at which drugs are combined.
#' @param .noReplicates \code{numeric}; a value providing the number of replicates for selected drug combinations, representing technical replicates in a drug sensitivity screen.
#' @param .drugRepAttrib \code{character}; a predefined statement specifying, which drug treatments should be replicated. "all", both single drug treatments and drug combinations to be replicated,
#' "single", only single drug treatments to be replicated.
#'
#' @details This function is useful for combining drugs at given doses with each other. The combinations are generated from a list of drugs and doses, such as generated with \code{generateListofDoses()}. 
#' Drugs can be combined at individual doses, or at a certain dose range, in which a low index corresponds with a lower dose. i.e. a value of 1 selects for the first, and consequently, for the lowest dose among the list of doses. 
#' The range of 2:5 selects the second lowest up to the fifth dose. A negative number will selected all doses and excluded the dose at the provided index.
#' The default is to combine drugs at all doses with each other.
#' 
#' Each drug combination or drug treatment can be replicated into technical replicates. Depending on the selection, either all drug treatments are replicated or only single drug treatment.
#' If the replication number is not provided, then drugs will not be replicated and only one occurrence of each unique drug dose combination is retained.
#' If no group for replication is selected, all drug combinations will be replicated by default.
#' 
#' The function then generates a list of combinations containing a column with the name, the doses and its corresponding unit for each of the drug pairs. 
#' This list can then be used for the generation of a dispensing file with \code{generateDispensingData()}.
#' 
#' @return Returns an object of class "data.frame", with an essential column of ['drug'] names, ['doses'] and ['units'] for each drug pair.
#'
#' @references Robert Hanes et al.
#' @author Robert Hanes
#' @note Version: 2021.03
#'
#' @examples
#' \donttest{\dontrun{
#' library(metascreen)
#' combineDrugs(listofDoses, .combineDoses=c(2:5), .noReplicates = 3, .drugRepAttrib = "single")
#' }}
#'
#' @keywords drug screen list drugs combination
#'
#'
#' @importFrom stats setNames
#' @importFrom utils type.convert
#' 
#' @export

combineDrugs <- function(listofDoses, .combineDoses, .noReplicates = 1, .drugRepAttrib = c("all", "single")){

  # Check, if arguments are provided and in the proper format
  if(missing(listofDoses)){stop("Data missing! Please provide a list of doses.", call. = TRUE)}
  if(!all(sapply(c("Drug", "Dose", "Unit"), function(x) any(grepl(x, names(listofDoses), ignore.case = TRUE))))){stop("List of doses has missing data! Please provide a data frame containing a column for the drug name ['Drug'], the dose ['Dose'] and the unit ['Unit'].", call. = TRUE)}

  # If arguments are missing, set them to default
  if(missing(.combineDoses)){message("Doses to be combined not specified. Combining all doses.")}
  if(missing(.combineDoses)){.combineDoses = 1:max(unique(table(listofDoses$Drug)))}
  if(missing(.noReplicates)){message("Number of drug replicates not specified. Using default: 1")}
  if(missing(.drugRepAttrib)){message("Attribute for drug replicates not set. Applying replicates to all drugs.")}

  .drugRepAttrib <- match.arg(.drugRepAttrib)

  

  cat("\n > Generating a list of combinations. Combining drugs...")
  
  # Create a list of drugs
  listofDrugs <- unique(listofDoses[[grep("Drug", names(listofDoses), ignore.case = TRUE)]])
  
  # Create a list of doses for each drug
  doseList <- sapply(listofDrugs, function(x) list(listofDoses[which(listofDoses$Drug == x), grepl("Dose", names(listofDoses))]))
  # Select only the designated four doses for the combination treatment
  doseList <- lapply(doseList, function(x) sort(as.numeric(x))[.combineDoses])
  


  # Combine each drug~dose combination with each other
  listofCombinations <- expand.grid(Drug.1 = paste(listofDoses$Drug, listofDoses$Dose, listofDoses$Unit, sep = ":"), Drug.2 = paste(listofDoses$Drug, listofDoses$Dose, listofDoses$Unit, sep = ":"), stringsAsFactors = FALSE)
  
  # Separate drug and dose from the drug~dose combinations
  listofCombinations <- cbind(setNames(data.frame(do.call('rbind', strsplit(as.character(listofCombinations$Drug.1),':', fixed=TRUE)), stringsAsFactors = FALSE), paste(c("Drug", "Dose", "Unit"), 1, sep = ".")),
                              setNames(data.frame(do.call('rbind', strsplit(as.character(listofCombinations$Drug.2),':', fixed=TRUE)), stringsAsFactors = FALSE), paste(c("Drug", "Dose", "Unit"), 2, sep = ".")))
  
  # Transform drug doses into numeric values
  listofCombinations <- transform(listofCombinations, Dose.1 = as.numeric(Dose.1), Dose.2 = as.numeric(Dose.2))
  
  # Remove same drug combinations, unless they share the same dose
  # Same drug combinations at the same dose are considered single drug treatments
  listofCombinations <- subset(listofCombinations, (Drug.1 == Drug.2 & Dose.1 == Dose.2 | Drug.1 != Drug.2))

  # Remove drug combinations at doses, that are not meant to be combined
  # Combine only doses that have been selected, while keeping the other doses as a single drug treatment
  # listofCombinations <- listofCombinations[sapply(1:nrow(listofCombinations), function(x) (listofCombinations[x, "Drug.1"] == listofCombinations[x, "Drug.2"] | 
  #                                                                                          listofCombinations[x, "Drug.1"] != listofCombinations[x, "Drug.2"] & 
  #                                                                                            listofCombinations[x,"Dose.1"] %in% doseList[[match(listofCombinations[x,"Drug.1"], names(doseList))]] & 
  #                                                                                            listofCombinations[x,"Dose.2"] %in% doseList[[match(listofCombinations[x,"Drug.2"], names(doseList))]] )),]
  listofCombinations <- listofCombinations[
    listofCombinations[, "Drug.1"] == listofCombinations[, "Drug.2"] | 
      sapply(1:nrow(listofCombinations), function(x) {
        listofCombinations[x,"Dose.1"] %in% doseList[[listofCombinations[x,"Drug.1"]]] &
          listofCombinations[x,"Dose.2"] %in% doseList[[listofCombinations[x,"Drug.2"]]]
      }), ]
  
  # Merge name, dose and unit for each drug pair 
  listofCombinations <- data.frame(Drug.1 = apply(listofCombinations[grep("1", names(listofCombinations), value=TRUE)], MARGIN = 1, FUN = function(i) paste(i, collapse = ":")),
                                   Drug.2 = apply(listofCombinations[grep("2", names(listofCombinations), value=TRUE)], MARGIN = 1, FUN = function(i) paste(i, collapse = ":")))
  # Reset row names to default
  rownames(listofCombinations) <- NULL
  
  
  # Removing duplicate combinations:
  # by sorting the rows and eliminating duplicates from a transposed matrix
  listofCombinations <- listofCombinations[!duplicated(t(apply(listofCombinations, 1, sort))), ]
  
  
  
  # CHECK number of combinations
  # Note: combine all other drugs at all their doses with one drug and all it's doses. do that for all the drugs
  #       divide that by half avoiding transposed combinations
  #       add single drug combinations each drug at each dose combined with itself
  # [ ( no. of drugs - 1 * no. of doses ) * no. of doses * no. of drugs ] / [ 2 ] + [ no. of drugs * no. of single drug doses ]
  cat('\r', ifelse(
    ((length(listofDrugs)-1) * unique(sapply(doseList, length)) * unique(sapply(doseList, length)) * length(listofDrugs)) / 2 + (length(listofDrugs) * unique(table(listofDoses$Drug)))
    == nrow(listofCombinations), "Combinations passed check! A list of combinations has been succesfully generated.", "Combinations failed check!"), "\n")
  
  
  
  # Adding replicates for an individual treatment group
  switch(as.character(.drugRepAttrib),
         "all"={
           # Adding replicates to all drug treatments
           listofCombinations <- listofCombinations[rep(row.names(listofCombinations), .noReplicates), ]},
         "single"={
           # Adding triplicates to single drug treatments,
           # only if the occurrence of the single drug treatments is equal 1 to avoid multiple generation of triplicates
           if(table(unlist(listofCombinations[which(listofCombinations$Drug.1==listofCombinations$Drug.2),][1]))[[2]] == 1){
             listofCombinations <- listofCombinations[rep(row.names(listofCombinations), ifelse(listofCombinations$Drug.1==listofCombinations$Drug.2, .noReplicates, 1)), ]
           }},{ message("Please select a valid treatment group to add replicates to. The treatment groups supported are \"all\" or \"single\" drug treatments.") })
  
  
  
  # Separate name, dose and unit for each drug pair
  listofCombinations <- cbind(setNames(data.frame(do.call('rbind', strsplit(as.character(listofCombinations$Drug.1),':', fixed=TRUE)), stringsAsFactors = FALSE), paste(c("Drug", "Dose", "Unit"), 1, sep = ".")),
                              setNames(data.frame(do.call('rbind', strsplit(as.character(listofCombinations$Drug.2),':', fixed=TRUE)), stringsAsFactors = FALSE), paste(c("Drug", "Dose", "Unit"), 2, sep = ".")))
  listofCombinations <- utils::type.convert(listofCombinations, as.is = TRUE)
  
  return(listofCombinations)
  
}
