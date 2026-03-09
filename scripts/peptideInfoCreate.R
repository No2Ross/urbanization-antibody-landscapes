library(dplyr)
library(ggplot2)
library(stringr)
library(ggpubr)
library(rstatix)
library(readxl)
library(reshape2)
library(resample)

library(factoextra)

#library(beer)

setwd("C:/Users/rflaidlaw/Documents/CapTan/AnalysisV2/PhIPseq/")

area_palette <- c("Rural Senegalese" = "#6592c3", "Urban Senegalese" = "#ee8d58", "Urban Dutch" = "#b0aeae")

metaData <- as.data.frame(read_xlsx("phipseq_metadata.xlsx"))
row.names(metaData) <- metaData$ID

phipSeq_counts <- read.csv("res_counts_annotated.csv", row.names = 2)
phipSeq_FC <- read.csv("res_foldchange_annotated.csv", row.names = 2)

#Get a dataframe with the peptide information 
peptideInfo <- phipSeq_counts[, c(2:7)]

#We need to create from 'uniref_func' two new columns: Tax and Protein
peptideInfo$Tax_uniref <- ""
peptideInfo$Protein_uniref <- ""
peptideInfo$intCol <- ""

peptideInfo$intCol[peptideInfo$uniref_func != ""] <- unlist(str_split(peptideInfo$uniref_func[peptideInfo$uniref_func != ""], " TaxID"))[c(TRUE, FALSE)]
peptideInfo$Tax_uniref[peptideInfo$intCol != ""] <- unlist(str_split(peptideInfo$intCol[peptideInfo$intCol != ""], " n=\\d+ Tax="))[c(FALSE, TRUE)]
peptideInfo$Protein_uniref[peptideInfo$intCol != ""] <- unlist(str_split(peptideInfo$intCol[peptideInfo$intCol != ""], " n=\\d+ Tax="))[c(TRUE, FALSE)]

peptideInfo$intCol <- NULL

peptideInfo$Pathogen <- ""
peptideInfo$Protein <- ""

####### Allergen processing ######
#There are columns which contain allergen information from different databases
#Create a column for each database so we can get the allergen name

allergenSplitString <- " \\(\\d+/\\d+\\) \\(\\d+/\\d+\\) | False | \\(\\d+/\\d+\\)\\,"

#There is one allergen entry which also includes a fummy name. Remove the fummy name part
peptideInfo$full.name[str_detect(peptideInfo$full.name, "iedb_name False fummy_name")] <- str_replace(peptideInfo$full.name[str_detect(peptideInfo$full.name, "iedb_name False fummy_name")],
                                                                                                      " fummy_name .*", "")

head(peptideInfo$full.name[str_detect(peptideInfo$full.name, "allergome_name")])

peptideInfo$allergen_name <- ""
peptideInfo$allergome_name <- ""
peptideInfo$allergenonline_name <- ""
peptideInfo$SDAP_name <- ""
peptideInfo$iedb_name <- ""

peptideInfo$allergen_name[str_detect(peptideInfo$full.name, "allergome_name")] <- unlist(str_extract_all( peptideInfo$full.name[str_detect(peptideInfo$full.name, "allergome_name")]
                                                                                                            , "(?<=allergen_name ).+(?= allergome_name)" ))
peptideInfo$allergen_name[str_detect(peptideInfo$full.name, "allergome_name")] <- str_replace_all(peptideInfo$allergen_name[str_detect(peptideInfo$full.name, "allergome_name")], 
                                                                                                  " \\(\\d+/\\d+\\) \\(\\d+/\\d+\\)",
                                                                                                  "")

peptideInfo$allergen_name[str_detect(peptideInfo$full.name, "allergome_name")] <- unlist(str_extract_all( peptideInfo$full.name[str_detect(peptideInfo$full.name, "allergome_name")]
                                                                                                          , "(?<=allergen_name ).+(?= allergome_name)" ))
peptideInfo$allergen_name[str_detect(peptideInfo$full.name, "allergome_name")] <- str_replace_all(peptideInfo$allergen_name[str_detect(peptideInfo$full.name, "allergome_name")], 
                                                                                                  " \\(\\d+/\\d+\\) \\(\\d+/\\d+\\)",
                                                                                                  "")

allergenDBNames <- c("allergen_name", "allergome_name", "allergenonline_name", "SDAP_name")

for(i in 2:length(allergenDBNames)){
  
  currentDB <- allergenDBNames[i-1]
  nextDB <- allergenDBNames[i]
  
  current <- peptideInfo$full.name[str_detect(peptideInfo$full.name, currentDB)]
  
  peptideInfo[[currentDB]][str_detect(peptideInfo$full.name, "allergen_name")] <- unlist(
                                                                              str_extract_all( current,
                                                                                             paste0("(?<=", currentDB, " ).+(?= ",nextDB,")" ) 
                                                                                             ) 
                                                                              )
  
  peptideInfo[[currentDB]][str_detect(peptideInfo$full.name, "allergen_name")] <- str_replace_all(peptideInfo[[currentDB]][str_detect(peptideInfo$full.name, "allergen_name")], 
                                                                                                    " \\(\\d+/\\d+\\) \\(\\d+/\\d+\\)",
                                                                                                    "")
  
}

#We do not format the iedb name, as the information encoded is not compatible to the other columns


#The "Pathogen" column value will be "Allergen"
peptideInfo$Pathogen[peptideInfo$allergen_name != ""] <- "Allergen"

#The allergenonline and SDAP names sometimes have this pattern following the names: ".0101". The numbers can vary, but its always got a full stop
#Remove this
sum(str_detect(peptideInfo$allergen_name, "\\."))
sum(str_detect(peptideInfo$allergome_name, "\\."))
sum(str_detect(peptideInfo$allergenonline_name, "\\."))
sum(str_detect(peptideInfo$SDAP_name, "\\."))

peptideInfo$allergenonline_name <- str_replace_all(peptideInfo$allergenonline_name, "\\..*", "")
peptideInfo$SDAP_name <- str_replace_all(peptideInfo$SDAP_name, "\\..*", "")

sum(str_detect(peptideInfo$allergenonline_name, "\\."))
sum(str_detect(peptideInfo$SDAP_name, "\\."))

#Some databases (allergome and allergenonline) have multiple results. These need to be reckoned with
#The pattern to split by is " (n/n), ". There is none for allergenome or allergenonline which have double digits in the numerator or denominator
sum(str_detect(peptideInfo$allergen_name, "\\,"))
sum(str_detect(peptideInfo$allergome_name, "\\,"))
sum(str_detect(peptideInfo$allergenonline_name, "\\,"))
sum(str_detect(peptideInfo$SDAP_name, "\\,"))

pattern <- " \\(\\d+/\\d+\\)\\, "

#Create new columns for storing the split values
#The max number of options in the column is the maximum number of commas detected + 1
for(i in 1:(max(str_count(peptideInfo$allergome_name, "\\,")) + 1) ){
  peptideInfo[[paste0("allergome_name_", i)]] <- ""
  
}

for(i in 1:(max(str_count(peptideInfo$allergenonline_name, "\\,")) + 1) ){
  peptideInfo[[paste0("allergenonline_name_", i)]] <- ""
  
}

for(i in 1:length(peptideInfo$allergome_name)){
  
  if( str_detect(peptideInfo$allergome_name[i], "\\,") ){
    
    current <- unlist(str_split(peptideInfo$allergome_name[i], pattern))
    
    #Iterate through each possible allergen, and add to a new column in peptideInfo
    for(j in 1:length(current)){
      
      peptideInfo[[paste0("allergome_name_", j)]][i] <- current[j]
      
    }
    
  }
  
  #If there is only one possible epitope, assign it to allergome_name_1 column
  else{peptideInfo[[paste0("allergome_name_1")]][i] <- peptideInfo$allergome_name[i]}
  
  
}

for(i in 1:length(peptideInfo$allergenonline_name)){
  
  if( str_detect(peptideInfo$allergenonline_name[i], "\\,") ){
    
    current <- unlist(str_split(peptideInfo$allergenonline_name[i], pattern))
    
    #Iterate through each possible allergen, and add to a new column in peptideInfo
    for(j in 1:length(current)){
      
      peptideInfo[[paste0("allergenonline_name_", j)]][i] <- current[j]
      
    }
    
  }
  
  #If there is only one possible epitope, assign it to allergenonline_name_1 column
  else{peptideInfo[[paste0("allergenonline_name_1")]][i] <- peptideInfo$allergenonline_name[i]}
  
}

#Remove the original allergenonline and allergome columns as they are now redundant
peptideInfo$allergome_name <- NULL
peptideInfo$allergenonline_name <- NULL

#Identify the allergen name for each row and assign it the "Protein" column

allergenCols <- c("allergen_name", "SDAP_name", 
                  paste0("allergome_name_", seq(1,3,1)), 
                  paste0("allergenonline_name_", seq(1,2,1))
                  )

#Iterate through each row and find the occurence of each of the values
for(i in 1:nrow(peptideInfo)){
  
  if(peptideInfo[i, "Pathogen"] == "Allergen"){
    
    current <- as.character(peptideInfo[i, allergenCols])
    
    #Remove "" and "False"
    current <- current[!(current %in% c("", "False", "Unassigned"))]
    
    currentCount <- table(current)
    
    #Find the options with the highest occurence
    maxIndex <- unname( which(currentCount == max(currentCount)) )
    
    #Give the unknown's a unique name so we can keep them seperate later
    if(length(currentCount) == 0){peptideInfo[i, "Protein"] <- paste0("Unknown_", i)}
    
    else if(length(maxIndex) == 1){peptideInfo[i, "Protein"] <- names(currentCount)[maxIndex]}
    
    #If we have multiple epitopes which are tied
    else{
      
      combinedCurrent <- paste(names(currentCount)[maxIndex], collapse = "/")
      
      peptideInfo[i, "Protein_uniref"] <- combinedCurrent
      
    }
    
  }
  
}

#Remove the created allergen columns
peptideInfo <- peptideInfo[, setdiff(colnames(peptideInfo), allergenCols)]


###### Corona processing ######

#If the peptide ID contains corona, assign the 'virus_name' to Pathogen
peptideInfo$Pathogen[str_detect(row.names(peptideInfo), "corona")] <- peptideInfo$virus_name[str_detect(row.names(peptideInfo), "corona")]

##### other processing #####
capturePattern <- "^VFG"
splitPattern <- "\\] \\["

VFG_complexNames <- unlist(
              str_split(peptideInfo[str_detect(peptideInfo$full.name, capturePattern),"full.name"], splitPattern)
              )[c(FALSE, TRUE)]

#Remove the end square bracket
VFG_complexNames <- str_replace_all(VFG_complexNames, "\\]", "")
VFG_simpleNames <- stringr::word(VFG_complexNames, 1,2, sep = " ")

peptideInfo$VFGComplex <- ""
peptideInfo$VFGSimple <- ""

peptideInfo$VFGComplex[str_detect(peptideInfo$full.name, capturePattern)] <- VFG_complexNames
peptideInfo$VFGSimple[str_detect(peptideInfo$full.name, capturePattern)] <- VFG_simpleNames

#there are others which start with "positive () 
#It seems for this one, the pathogen and the structure are seperated by a double space
capturePattern <- "^positive \\("
replacePattern <- "positive \\(.*\\): "

positiveVirus <- peptideInfo[str_detect(peptideInfo$full.name, capturePattern),"full.name"]
positiveVirus <- str_replace_all(positiveVirus, replacePattern, "")

#Human respiratory syncytial virus and Human immunodeficiency virus 1 do not follow this pattern of splitting the pathogen from the epitope
splitPattern <- "  "

#Because of this, we need to add it
positiveVirus <- str_replace_all(positiveVirus, 
                                 "Human respiratory syncytial virus ", 
                                 "Human respiratory syncytial virus  ")

positiveVirus <- str_replace_all(positiveVirus, 
                                 "Human immunodeficiency virus 1 ", 
                                 "Human immunodeficiency virus 1  ")

positiveVirusNames <- str_split(positiveVirus, splitPattern)       

positiveVirusNames <- unlist(positiveVirusNames)[c(TRUE, FALSE)]

peptideInfo$positiveVirusNames <- ""
peptideInfo$positiveVirusNames[str_detect(peptideInfo$full.name, "^positive \\(")] <- positiveVirusNames

#Others start with "NP_" and ends with "]" (this has the name of the pathogen in it, surrounded by square brackets)
capturePattern <- "^NP_.*\\]$"
replacePattern1 <- "^NP_\\d+\\.[12345] "

test <- peptideInfo[str_detect(peptideInfo$full.name, capturePattern),]

np <- peptideInfo[str_detect(peptideInfo$full.name, capturePattern),"full.name"]
np <- str_replace_all(np, replacePattern1, "")

splitPattern <- " \\["
npName <- unlist(str_split(np, splitPattern))[c(FALSE, TRUE)]
npProtein <- unlist(str_split(np, splitPattern))[c(TRUE, FALSE)]

replacePattern2 <- "\\]"
npName <- str_replace_all(npName, replacePattern2, "")

peptideInfo$npName <- ""
peptideInfo$npName[str_detect(peptideInfo$full.name, capturePattern)] <- npName

peptideInfo$npProtein <- ""
peptideInfo$npProtein[str_detect(peptideInfo$full.name, capturePattern)] <- npProtein

#Combine the other information into the Pathogen column (and protein, when applicable)
peptideInfo$Pathogen[peptideInfo$positiveVirusNames != ""] <- peptideInfo$positiveVirusNames[peptideInfo$positiveVirusNames != ""]
peptideInfo$Pathogen[peptideInfo$npName != ""] <- peptideInfo$npName[peptideInfo$npName != ""]
peptideInfo$Pathogen[peptideInfo$VFGSimple != ""] <- peptideInfo$VFGSimple[peptideInfo$VFGSimple != ""]

peptideInfo$Protein[peptideInfo$npProtein != ""] <- peptideInfo$npProtein[peptideInfo$npProtein != ""]

#If the information in Pathogen is blank, take it from the Tax_uniref column
peptideInfo$Pathogen[peptideInfo$Pathogen == ""] <- peptideInfo$Tax_uniref[peptideInfo$Pathogen == ""]

#Create new column to store more complex pathogen distinctions (i.e. subspecies) 
#Only makes a difference for the full.name values starting with VFG
peptideInfo$PathogenComplex <- peptideInfo$Pathogen
peptideInfo$PathogenComplex[peptideInfo$VFGComplex != ""] <- peptideInfo$VFGComplex[peptideInfo$VFGComplex != ""]

#Remove unneeded columns
peptideInfo$VFGComplex <- NULL
peptideInfo$VFGSimple <- NULL
peptideInfo$npName <- NULL
peptideInfo$npProtein <- NULL
peptideInfo$positiveVirusNames <- NULL

##### Save data #####

#Create two new columns: "epitopeGroup" & "epitopeUnique"
#"epitopeGroup" contains the merging of the "Protein" and the "Tax" column
#"epitopeUnique" contains the merging of the "Protein", "Tax" and row names
peptideInfo$epitopeGroup <- paste(peptideInfo$Protein, peptideInfo$Tax, sep = "_")
peptideInfo$epitopeUnique <- paste(peptideInfo$Protein, peptideInfo$Tax, row.names(peptideInfo), sep = "_")



write.csv(peptideInfo, "PhIPseq_peptideInfo_processed.csv")
