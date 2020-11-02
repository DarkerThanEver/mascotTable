# adjust these two source() instructions to correct locations!

source("../proteinDrawing/proteinDrawing.R")
source("../textData/textData.R")


# general functions for 'translating' some of the variables/columns in the mascot tables

# take strings in the format "0.002003006000.0" (mascot pep_var_mod_pos) and returns a data.frame with the positions
# and number of the modification (0 = no modification)
getModPositions <- function(theVarModPos = NA){
  findThem <- as.data.frame(str_locate_all(theVarModPos,"[^0\\.]")[[1]])
  if (nrow(findThem)>0){
    whichOnes <- as.numeric()
    for (counter in 1:nrow(findThem)){
      whichOnes <- append(whichOnes, strtoi(substr(theVarModPos, findThem$start[counter],findThem$start[counter])))
    }
    findThem$modification <- whichOnes
    return(findThem)
  } else {
    return(NA)
  }
}

# takes a string from the mascot export (pep_scan_title) and extracts the retention time in numeric format
# note: dependent on the number being the penultimate element when the string is split on space
# example: 'MS/MS of 632.3391418 1+ at 18.6032833333333 mins' --> this string returns 18.603
getRt <- function(scan_title){
  tempStr <- strsplit(scan_title, split = " ")
  return(unlist(lapply(tempStr,function(x){return(as.numeric(x[length(x)-1]))})))
}

# example of data.frame specifying a single fixed modification
# note: you cannot combine different amino acids, aa column should be a single letter per row
addFixedModsCM <- data.frame(aa = "C",mod = "Carbamidomethyl", stringsAsFactors = FALSE)

# generates an experiment table based on mascot export data (proteinData) that is loaded via mascotFile$readRaw
# (or read.csv); must be in data.frame format (tibbles probably also work)
# whichProtein can be a number or text
#   number is based on the order of the function: unique(...data...$prot_acc)
#   use text for better reproducibility
# order is important for routines (derived from) from drawProteins: it allows for multiple chains in one
#  graph Note: when combining separate drawProteins data.frames change before combining!
#  adddFixedMods: can be used to add fixed variables used in the search.
# Note: the fixed modifications are NOT specified anywhere in the mascot export table
#       thie addFixedModifications parameter is the only way to add them
# varModTable is used together with the var_mod_pos column from the mascot export table ("0.002003010.0")
# to determine which modifications are present
# note: if varModTable is set to NA, then the part of the table with all modifications specified per location
#       will not be created, which means that there will only be type = "COVERAGE"/"CHAIN" for that table 
# minimumExpect:allows selection of only those peptide identifications which have at least the specified Mascot Expect value
# (pep_expect column)
# Note: be careful! There may be more than one hit for a peptide, the minimumExpect in the table is only one of the pep_expect values,
#       same goes (less impact but still...) for the rt (retention tine)
# Note: when exporting a csv formatted result file from mascot, the following must be included (not standard)
#       protein sequence, peptide start and end. It is advisable to be careful when setting p-value (standard = 0.05)
mascotExperimentTable <- function(proteinData, whichProtein, order = 1, addFixedMods = NA, varModTable = NA,
                                  minimumExpect = 0.05){
  if (is.numeric(whichProtein)){
    whichProtein <- unique(proteinData$prot_acc)[whichProtein]
  }
  proteinData <- proteinData %>% filter(pep_expect <= minimumExpect)
  theData <- proteinData %>%
    filter(prot_acc == whichProtein) %>%
    distinct(across(contains(c("pep_seq", "pep_var_mod", "pep_var_mod_pos"))), .keep_all = TRUE)
  whichProtein <- gsub(".*[|.]","",whichProtein)   
  tempdf <- createExperimentTable(type = as.character("CHAIN"),
                                  description = whichProtein,
                                  begin = 1,
                                  end = nchar(theData$prot_seq[1])+1,
                                  length = nchar(theData$prot_seq[1]),
                                  order = order,
                                  sequence = theData$prot_seq[1],  # chain contains the protein sequence from the database
                                  mods = as.character(NA),
                                  expMr = as.numeric(NA),
                                  calcMr = theData$prot_mass[1],
                                  rt = as.numeric(NA),
                                  expect = as.numeric(NA),
                                  varModPos = as.character(NA))
  theData <- theData %>%
    arrange(pep_start)
  for (counter in 1:nrow(theData)){
    tempdf <- bind_rows(tempdf,
                        createExperimentTable(type = as.character("COVER"),
                                              description = whichProtein,
                                              begin = theData$pep_start[counter],
                                              end = theData$pep_end[counter]+1,
                                              length = nchar(theData$pep_seq[counter]),
                                              order = order,
                                              sequence = theData$pep_seq[counter],
                                              mods = theData$pep_var_mod[counter],
                                              expMr = theData$pep_exp_mr[counter],
                                              calcMr = theData$pep_calc_mr[counter],
                                              rt = getRt(theData$pep_scan_title[counter]),
                                              expect = theData$pep_expect[counter],
                                              varModPos = theData$pep_var_mod_pos[counter]))
  }
  if (!identical(addFixedMods,NA)){
    for (counter in 1:nrow(addFixedMods)){
      locations <- as.data.frame(str_locate_all(theData$prot_seq[1],addFixedMods$aa)[[counter]])
      if (nrow(locations)>0){
        for (counter2 in 1:nrow(locations)){
          tempdf <- bind_rows(tempdf,
                              createExperimentTable(type = as.character("MOD_RES"),
                                                    description = addFixedMods$mod[counter],
                                                    begin = locations$start[counter2],
                                                    end = locations$end[counter2],
                                                    length = locations$end[counter2]-locations$start[counter2]+1,
                                                    order = order,
                                                    sequence = as.character(NA),
                                                    mods = as.character(NA),
                                                    expMr = as.numeric(NA),
                                                    calcMr = as.numeric(NA),
                                                    rt = NA,
                                                    expect = as.numeric(NA),
                                                    varModPos = as.character(NA)))
        }
      }
    }
  }
  if (!identical(varModTable,NA)){
    theData2 <- theData %>%
      filter(nchar(pep_var_mod)>0)
    if (nrow(theData2)>0){
      modLocations <- lapply(theData2$pep_var_mod_pos,function(x){getModPositions(x)})
      for (counter in 1:length(modLocations)){
        currentLocation <- modLocations[[counter]]
        for (counter2 in 1:nrow(currentLocation)){
          if ((currentLocation$start[counter2]) %in% c(3:(nchar(theData2$pep_var_mod_pos[counter])-2))){
            tempdf <- bind_rows(tempdf,
                                createExperimentTable(type = as.character("MOD_RES"),
                                                      description = varModTable[currentLocation$modification[counter2]],
                                                      begin = theData2$pep_start[counter]+currentLocation$start[counter2]-3,
                                                      end = theData2$pep_start[counter]+currentLocation$start[counter2]-3,
                                                      length = 1,
                                                      order = order,
                                                      sequence = as.character(NA),
                                                      mods = as.character(NA),
                                                      expMr = as.numeric(NA),
                                                      calcMr = as.numeric(NA),
                                                      rt = NA,
                                                      expect = as.numeric(NA),
                                                      varModPos = as.character(NA)))
          } else {
            tempdf <- bind_rows(tempdf,
                                createExperimentTable(type = as.character("MOD_RES"),
                                                      description = varModTable[currentLocation$modification[counter2]],
                                                      begin = theData2$pep_start[counter]+currentLocation$start[counter2]-1,
                                                      end = theData2$pep_start[counter]+currentLocation$start[counter2]-1,
                                                      length = 1,
                                                      order = order,
                                                      sequence = as.character(NA),
                                                      mods = as.character(NA),
                                                      expMr = as.numeric(NA),
                                                      calcMr = as.numeric(NA),
                                                      rt = NA,
                                                      expect = as.numeric(NA),
                                                      varModPos = as.character(NA)))
          }
        }
      }
    }
  }
  startcount <- 1
  for (counter in 1:nrow(tempdf)){
    if ((tempdf$type[counter] == "MOD_RES") & (grepl(";",tempdf$description[counter]))){
      tempdf$description[counter] <- str_trim(unlist(str_split(tempdf$description[counter],";"))[startcount], side = "both")
      startcount <- startcount + 1
    } else {
      startcount <- 1
    }
  }
  tempdf <- tempdf %>% distinct(.keep_all = TRUE) # remove duplicate entries
  # note: there will be seq with different expMr but same sequence in the cover section, this is because of modified peptides 
  tempdf <- tempdf[!(is.na(tempdf$description)),] # for some reason, <NA> mod_res come in table sometimes, this is to remove
  return(tempdf)
}

mascotFile <- R6Class("mascotFile",
                       inherit = textData,
                       public = list(
                         # overrides the textData initialize --> lines should NOT be permanent in memory
                         initialize = function(fileName){
                           super$initialize(theData = fileName, fileMode = TRUE)
                           invisible(self)
                         },
                         # to find the place where the table with mascot hit info starts eturns the
                         # skip = parameter for read.csv of mascot files
                         # note: in case of Mascot export, dependent on the text
                         #       "Protein hits"  being present (ONCE) one line before the actual table
                         # indicator : under normal circumstances, the string "Protein hits" signals
                         #             the beginning of the mascot search data.
                         searchSkip = function(indicator = "Protein hits"){
                           tempLines <- readLines(self$fileName)
                           return(which(grepl(indicator,tempLines))+1)
                         },
                         # reads the data from a mascot export CSV file
                         # if skip == NA then it uses searchSkip to attempt to determine the number
                         # of lines to skip before the proteomics data starts
                         # indicator : see searchSkip function
                         readRaw = function(skip = NA, indicator = "Protein hits"){
                           if (is.na(skip)){
                             toSkip <- self$searchSkip(indicator = indicator)
                           } else {
                             toSkip <- skip
                           }
                           theData <- read.csv(private$dataFile, skip = toSkip,
                                               header = TRUE, stringsAsFactors = FALSE)
                           return(theData)
                         },
                         # gives a list of the proteins via prot_acc column in the mascot table
                         # this can be used if initially unsure about the protein(s) to select
                         proteinsInFile = function(){
                           tempList <- self$readRaw()
                           return(unique(tempList$prot_acc))
                         },
                         # to find the names of the tables present in mascot export files
                         # returns those name
                         # note: dependence on presence of "--------------------------------------------------------" in the file
                         #       any deviation of this pattern will disrupt this routine
                         getTableNames = function(){
                           tempLines <- readLines(self$fileName)
                           return(unlist(lapply(
                             unlist(lapply(strsplit(tempLines[which(grepl("--------------------------------------------------------",tempLines))],
                                                    split = ","),function(x){x[[1]][1]})),
                             function(x){substr(x,2, nchar(x)-1)})))
                         },
                         # returns the raw lines of the specified table in mascot export files
                         # similar to getTableNames: dependent on "--------------------------------------------------------" being present
                         getTableLines = function(whichTable = "Header"){
                           tempLines <- readLines(self$fileName)
                           startPosition <- which(grepl(whichTable,tempLines))[1]
                           endPosition <- which(grepl("--------------------------------------------------------",
                                                      tempLines[startPosition+2:length(tempLines)]))[1]-1
                           if (is.na(endPosition)){
                             endPosition <- length(tempLines)
                           }
                           return(tempLines[startPosition+2:endPosition])
                         },
                         # returns the table in data.frame format from mascot export files
                         # header = FALSE in most cases, however variable & fixed modification tables should TRUE
                         # note : this function does NOT work properly on 'Protein hits' tables --> use readRaw
                         #        function for this
                         getTable = function(whichTable = "Header", header = FALSE){
                           tempLines <- self$getTableLines(whichTable = whichTable)
                           tempLines <- gsub('\"','',tempLines) # remove \"  
                           if (header){
                             tempdf <- data.frame()
                             cols <- unlist(strsplit(tempLines[1], split = ","))
                             tempdf <- matrix(as.character(),ncol = length(cols))
                             for (counter in 2:length(tempLines)){
                               templ <- unlist(strsplit(tempLines[counter], split = ","))
                               if (length(cols) > length(templ)){
                                 templ <- append(templ, rep("",(length(cols)-length(templ))))
                               } else {
                                 if (length(cols < length(templ))){
                                   templ2 <- templ
                                   templ <- templ[1:length(cols)]
                                   templ[length(templ)] <- paste(templ2[length(templ):length(templ2)], collapse = ",")
                                 }
                               }
                               tempdf <- rbind(tempdf, templ)
                             }
                             tempdf <- as.data.frame(tempdf)
                             colnames(tempdf) <- cols
                             rownames(tempdf) <- 1:nrow(tempdf)
                           } else {
                             tempdf <- data.frame(variable = as.character(), value = as.character())
                             for (counter in 1:length(tempLines)){
                               templ <- unlist(strsplit(tempLines[counter],split = ","))
                               if (is.na(templ[2])){
                                 templ <- append(templ,"")
                               }
                               tempdf <- bind_rows(tempdf, data.frame(variable = templ[1], value = templ[2]))
                             }
                           }
                           return(tempdf)
                         },
                         # to get the variable modification names
                         # note: the way grepl works this may give problems with multiple "Variable modifications"
                         #       strings present in the mascot export file
                         # note: the sequence of the strings is number 1.. as used in the tables
                         getVarModTable = function(){
                           tempLines <- readLines(self$fileName)
                           thePosition <- which(grepl("Variable modifications",str_split(tempLines, '\\"')))[2]
                           return(unlist(strsplit(unlist(strsplit(tempLines[thePosition],'\"'))[2],",")))
                         },
                         # reads the data from a mascot export table (CSV format) and runs it through the function proteinFormatted
                         # proteinFormatted requires protein sequnece and peptide begin/start & end to be present in the output
                         # this function attempts to determine how many lines to skip before the protein/peptide info table starts
                         # if it doesn't work properly, the skip parameter overrides this. Similarly getVarModTable is used to
                         # determine the variable modifications that were used in the search. If it doesn't work properly a
                         # character() list can be provided with the correct modifications. This character() list must be in
                         # the correct order! the parameter addFixed can be used to specify all the fixed modifications
                         # present. It takes the same format as with the proteinFormatted function (see there for details)
                         experimentTable = function(skip = NA, addFixedMods = NA,
                                                    whichProtein = 1, order = 1, minimumExpect = 0.001,
                                                    varModTable = NA){
                           theData <- self$readRaw(skip = skip)
                           if (identical(varModTable,NA)){
                             varModTable <- self$getVarModTable()
                           }
                           theData <- mascotExperimentTable(theData,whichProtein = whichProtein,order = order,
                                                            addFixedMods = addFixedMods, varModTable = varModTable,
                                                            minimumExpect = minimumExpect)
                           return(theData)
                         }
                       )
)