library(dplyr)
library(ggplot2)
library(stringr)
library(ggpubr)
library(rstatix)
library(readxl)
library(reshape2)
library(resample)
library(Seurat)

library(XML)
library(xml2)

#For the cochran trend test
library(DescTools)

library(wordcloud2)

library(factoextra)

#library(beer)

cleanBLASTResults <- function(currentEpitope){
  #If there are multiple options, choose one
  if(str_detect(currentEpitope, "\\>")){
    currentEpitope <- unlist(str_split(currentEpitope, " \\>"))
    currentEpitope <- currentEpitope[c(TRUE, rep(FALSE, (length(currentEpitope) -1) ))]
  }
  
  if(str_detect(currentEpitope, "RecName")){
    currentEpitope <- str_replace_all(currentEpitope, "RecName: Full=", "")
    
    #Some of these don't have a semicolon
    if(str_detect(currentEpitope, ";")){
      semiColonSplit <- unlist(str_split(currentEpitope, ";"))
      
      currentEpitope_protein <- semiColonSplit[c(TRUE, rep(FALSE, (length(semiColonSplit) -1) ))]
      currentEpitope_pathogen <- unlist(str_split(currentEpitope, " \\["))[c(FALSE, TRUE)]
      
      currentEpitope <- paste0(currentEpitope_protein, " [", currentEpitope_pathogen)
      
    }
    
    rm(currentEpitope_protein)
    rm(currentEpitope_pathogen)
    gc()
  }
  
  return(currentEpitope)
}

setwd("C:/Users/rflaidlaw/Documents/CapTan/AnalysisV2/PhIPseq/")

peptideInfo <- read.csv("PhIPseq_peptideInfo_processed.csv", row.names = 1)

#I should consider looking at all the pathogen results for a given epitope and picking the most common one

#Add a new column to show whether the information is from BLAST or not
peptideInfo$BLAST <- "No"

#What is the lengths of the epitopes
epitopeLengths <- str_length(peptideInfo$aa_seq)

hist(epitopeLengths)

#### Corona ####
coronaInfo <- peptideInfo[str_detect(row.names(peptideInfo), "corona"),]

#Extract the coronavirus epitopes which do not have a virus name
virusUnnamedInfo <- subset(coronaInfo, virus_name == "")

#Submit in the format of a fasta file
#Each sequence is two lines
#Line one is the epitope ID. It is structured: > [epitopeID]
#Line two is the epitope amino acid sequence

output <- "BLAST/input/virus.fasta"
fileConn <- file(output, open = "wt")

for (i in 1:nrow(virusUnnamedInfo)) {
  cat(paste0("> ", row.names(virusUnnamedInfo)[i]), file = output, sep = "\n", append = TRUE)
  
  #Need to remove the brackets from the end of the amino acid sequence
  aaSeq <- unlist(str_split(virusUnnamedInfo$aa_seq[i], " \\("))[c(TRUE, FALSE)]
  
  cat(paste0(aaSeq), file = output, sep = "\n", append = TRUE)
} 

close(fileConn)
rm(fileConn)
gc()

#Load in BLASTed virusEpitope data
virusBLAST <- read_xml("BLAST/output/virusBLAST.xml")

#Extract the IDs and the BLAST results
epitopeID <- xml_text(xml_find_all(virusBLAST, ".//Iteration_query-def"))
epitope_poss <- xml_text(xml_find_all(virusBLAST, ".//Hit_def"))
epitope_hitNum <- xml_integer(xml_find_all(virusBLAST, ".//Hit_num"))

epitope_firstHit <- which(epitope_hitNum == 1)

virusDF <- data.frame()

counter <- 1

#This code works on the assumption, that the most likely hit (in terms of E value) is displayed first
for(i in epitope_firstHit){
  
  #Some of the results are formatted differently. 
  #I try and clean them down in this section of code
  
  currentEpitope <- epitope_poss[i]
  
  #if(epitopeID[counter] == "corona2_12968"){stop()}
  
  if(str_detect(currentEpitope, "\\>")){
    currentEpitope <- unlist(str_split(currentEpitope, " \\>"))
    currentEpitope <- currentEpitope[c(TRUE, rep(FALSE, (length(currentEpitope) -1) ))]
  }
  
  if(str_detect(currentEpitope, "RecName")){
    currentEpitope <- str_replace_all(currentEpitope, "RecName: Full=", "")
    
    semiColonSplit <- unlist(str_split(currentEpitope, ";"))
    
    currentEpitope_protein <- semiColonSplit[c(TRUE, rep(FALSE, (length(semiColonSplit) -1) ))]
    currentEpitope_pathogen <- unlist(str_split(currentEpitope, " \\["))[c(FALSE, TRUE)]
    
    currentEpitope <- paste0(currentEpitope_protein, " [", currentEpitope_pathogen)
  }
  
  #extract the protein and pathogen information
  infoSplit <- unlist(str_split(currentEpitope, " \\["))
  
  protein <- infoSplit[c(TRUE, FALSE)]
  pathogen <- infoSplit[c(FALSE, TRUE)]
  
  #Remove the trailing ] from pathogen
  pathogen <- str_replace_all(pathogen, "\\]", "")
  
  currentDF <- data.frame("ID" = epitopeID[counter],
                          "Pathogen" = pathogen,
                          "Protein" = protein)
  
  virusDF <- rbind(virusDF, currentDF)
  
  counter <- counter + 1
  
}

#Add information to peptideInfo
peptideInfo[virusDF$ID,"Pathogen"] <- virusDF$Pathogen
peptideInfo[virusDF$ID,"PathogenComplex"] <- virusDF$Pathogen
peptideInfo[virusDF$ID,"Protein"] <- virusDF$Protein
peptideInfo[virusDF$ID,"BLAST"] <- "Yes"

rm(virusDF)
rm(currentDF)
rm(counter)
rm(protein)
rm(pathogen)
rm(infoSplit)
rm(currentEpitope)
rm(currentEpitope_pathogen)
rm(currentEpitope_protein)
rm(epitope_poss)
rm(epitopeID)
rm(virusBLAST)
rm(virusUnnamedInfo)
gc()



#### twist ####

#There are some twist epitopes which are shared across different things

twistInfo <- peptideInfo[str_detect(row.names(peptideInfo), "twist"),]

#Extract the twist epitopes which doesn't have both protein and pathogen information
twistUnnamedInfo <- subset(twistInfo, Pathogen == "" | Protein == "")

#Submit in the format of a fasta file
#Each sequence is two lines
#Line one is the epitope ID. It is structured: > [epitopeID]
#Line two is the epitope amino acid sequence

#twist has too many sequences to blast, so will need to do it in batches
#I could try BLAST+ locally, but the protein database is 394GB

#Try splitting it (roughly) 6 ways
twistList <- list()
twistList[["one"]] <- twistUnnamedInfo[1:200,]
twistList[["two"]] <- twistUnnamedInfo[201:400,]
twistList[["three"]] <- twistUnnamedInfo[401:600,]
twistList[["four"]] <- twistUnnamedInfo[601:750,]
twistList[["five"]] <- twistUnnamedInfo[751:900,]
twistList[["six"]] <- twistUnnamedInfo[901:1000,]
twistList[["seven"]] <- twistUnnamedInfo[1001:nrow(twistUnnamedInfo),]


for(j in names(twistList)){
  
  currentDF <- twistList[[j]]
  
  output <- paste0("BLAST/input/twist_",j,".fasta")
  fileConn <- file(output, open = "wt")
  
  for (i in 1:nrow(currentDF)) {
    cat(paste0("> ", row.names(currentDF)[i]), file = output, sep = "\n", append = TRUE)
    
    #Need to remove the brackets from the end of the amino acid sequence
    aaSeq <- unlist(str_split(currentDF$aa_seq[i], " \\("))[c(TRUE, FALSE)]
    
    cat(paste0(aaSeq), file = output, sep = "\n", append = TRUE)
  } 
  
  close(fileConn)
  rm(fileConn)
  gc()
  
}

twistDF <- data.frame()

for(i in list.files("BLAST/output/", pattern = "twist")){
  print(i)
  
  #Load in BLASTed twistEpitope data
  twistBLAST <- read_xml(paste0("BLAST/output/", i))
  
  #Extract the IDs and the BLAST results
  epitopeID <- xml_text(xml_find_all(twistBLAST, ".//Iteration_query-def"))
  epitope_poss <- xml_text(xml_find_all(twistBLAST, ".//Hit_def"))
  epitopeIteration_hits <- xml_find_all(twistBLAST, ".//Iteration_hits")
  
  counter <- 1
  
  #This code works on the assumption that the hits are ordered by e-value
  
  for(j in 1:length(epitopeIteration_hits)){
    
    #Some of the results are formatted differently. 
    #I try and clean them down in this section of code
    
    currentEpitope <- xml_text(xml_find_all(epitopeIteration_hits[[j]], ".//Hit_def"))
    
    #Deal with the possibility that there is no BLAST results for the epitope
    if(is.na(currentEpitope[1])){
      pathogen <- "unknownBLAST"
      protein = "unknownBLAST"
    }
    
    else{
      
      pathogenVector <- c()
      proteinVector <- c()
      
      #Extract the top 15 (if available) choices in terms of e-value
      if(length(currentEpitope) < 15){currentEpitope <- currentEpitope[1:length(currentEpitope)]}
      else{currentEpitope <- currentEpitope[1:15]}
      
      #Iterate through the different options and clean the
      for(choiceOption in 1:length(currentEpitope)){
        
        currentHit <- currentEpitope[choiceOption]
        
        currentHit <- cleanBLASTResults(currentHit)
        
        #extract the protein and pathogen information
        infoSplit <- unlist(str_split(currentHit, " \\["))
        
        #Some can have two or more square bracket sections. This can confuse things
        protein <- infoSplit[c(TRUE, FALSE, rep(FALSE, (length(infoSplit)-2) ) )]
        pathogen <- infoSplit[c(FALSE, TRUE, rep(FALSE, (length(infoSplit)-2) ) )]
        
        #Remove the trailing ] from pathogen
        pathogen <- str_replace_all(pathogen, "\\]", "")
        
        #To avoid case sensitive issues, convert the pathogen and protein to lower
        pathogen <- tolower(pathogen)
        protein <- tolower(protein)
        
        pathogenVector <- append(pathogenVector, pathogen)
        proteinVector <- append(proteinVector, protein)
        
      }
      
      pathogen <- table(pathogenVector)
      pathogen <- names(pathogen)[which(pathogen == max(pathogen))]
      
      protein <- table(proteinVector)
      protein <- names(protein)[which(protein == max(protein))]
      
      #If there are tied values, pick the first one
      if(length(protein) > 1){
        protein <- protein[1]
      }
      
      if(length(pathogen) > 1){
        pathogen <- pathogen[1]
      }
      
    }
    
    currentDF <- data.frame("ID" = epitopeID[counter],
                            "Pathogen" = pathogen,
                            "Protein" = protein)
    
    twistDF <- rbind(twistDF, currentDF)
    
    counter <- counter + 1
    
  }
  
}


#Add information to peptideInfo
peptideInfo[twistDF$ID,"Pathogen"] <- twistDF$Pathogen
peptideInfo[twistDF$ID,"PathogenComplex"] <- twistDF$Pathogen
peptideInfo[twistDF$ID,"Protein"] <- twistDF$Protein
peptideInfo[twistDF$ID,"BLAST"] <- "Yes"




#### other ####

#Lots of the uni_ref names are incorrect. For example it classified haemophilus influenza protein as "Bacteria"
#The proteins however are probably correct

#Also, some epitopes don't have any information (in full.name, corona or uniref)

#Extract the epitopes which don't have protein or pathogen information and are not corona or twist peptide list IDs.
otherInfo <- subset(peptideInfo, Pathogen == "" | Protein == "")
otherInfo <- subset(otherInfo, str_detect(row.names(otherInfo), "agilent"))

#Submit in the format of a fasta file
#Each sequence is two lines
#Line one is the epitope ID. It is structured: > [epitopeID]
#Line two is the epitope amino acid sequence

otherList <- list()

fileIndex <- {set.seed(42); sample( seq(1,round(nrow(otherInfo) / 200), 1), nrow(otherInfo), replace = TRUE )}

hist(fileIndex)

for(i in 1:round(nrow(otherInfo) / 200)){
  
  otherList[[as.character(i)]] <- otherInfo[which(fileIndex == i),]
  
}



for(j in names(otherList)){
  
  currentDF <- otherList[[j]]
  
  output <- paste0("BLAST/input/other_",j,".fasta")
  fileConn <- file(output, open = "wt")
  
  for (i in 1:nrow(currentDF)) {
    cat(paste0("> ", row.names(currentDF)[i]), file = output, sep = "\n", append = TRUE)
    
    #Need to remove the brackets from the end of the amino acid sequence
    aaSeq <- unlist(str_split(currentDF$aa_seq[i], " \\("))[c(TRUE, FALSE)]
    
    cat(paste0(aaSeq), file = output, sep = "\n", append = TRUE)
  } 
  
  close(fileConn)
  rm(fileConn)
  gc()
  
}

#The organism choice should be the one which has the most hits. This of course could be biased towards the
#organism with the most entries in the database, but I imagine the people who made the library
#did not pick incredibly obscure organisms to get the sequences from

#We take the top 15 (where available) hits in terms of E-value and pick the pathogen and protein which is the most common response
#If there are multiple of the same, pick one randomly

otherDF <- data.frame()

for(i in list.files("BLAST/output/", pattern = "other")){
  print(i)
  
  #Load in BLASTed twistEpitope data
  otherBLAST <- read_xml(paste0("BLAST/output/", i))
  
  #Extract the IDs and the BLAST results
  epitopeID <- xml_text(xml_find_all(otherBLAST, ".//Iteration_query-def"))
  epitope_poss <- xml_text(xml_find_all(otherBLAST, ".//Hit_def"))
  epitopeIteration_hits <- xml_find_all(otherBLAST, ".//Iteration_hits")
  
  counter <- 1
  
  #This code works on the assumption that the hits are ordered by e-value
  
  for(j in 1:length(epitopeIteration_hits)){
    
    #Some of the results are formatted differently. 
    #I try and clean them down in this section of code
    
    currentEpitope <- xml_text(xml_find_all(epitopeIteration_hits[[j]], ".//Hit_def"))
    
    #Deal with the possibility that there is no BLAST results for the epitope
    if(is.na(currentEpitope[1])){
      pathogen <- "unknownBLAST"
      protein = "unknownBLAST"
    }
      
    else{
      
      pathogenVector <- c()
      proteinVector <- c()
      
      #Extract the top 15 (if available) choices in terms of e-value
      if(length(currentEpitope) < 15){currentEpitope <- currentEpitope[1:length(currentEpitope)]}
      else{currentEpitope <- currentEpitope[1:15]}
      
      #Iterate through the different options and clean the
      for(choiceOption in 1:length(currentEpitope)){
        
        currentHit <- currentEpitope[choiceOption]
        
        currentHit <- cleanBLASTResults(currentHit)
        
        #extract the protein and pathogen information
        infoSplit <- unlist(str_split(currentHit, " \\["))
        
        #Some can have two or more square bracket sections. This can confuse things
        protein <- infoSplit[c(TRUE, FALSE, rep(FALSE, (length(infoSplit)-2) ) )]
        pathogen <- infoSplit[c(FALSE, TRUE, rep(FALSE, (length(infoSplit)-2) ) )]
        
        #Remove the trailing ] from pathogen
        pathogen <- str_replace_all(pathogen, "\\]", "")
        
        #To avoid case sensitive issues, convert the pathogen and protein to lower
        pathogen <- tolower(pathogen)
        protein <- tolower(protein)
        
        pathogenVector <- append(pathogenVector, pathogen)
        proteinVector <- append(proteinVector, protein)
        
      }
      
      pathogen <- table(pathogenVector)
      pathogen <- names(pathogen)[which(pathogen == max(pathogen))]
      
      protein <- table(proteinVector)
      protein <- names(protein)[which(protein == max(protein))]
      
      #If there are tied values, pick the first one
      if(length(protein) > 1){
        protein <- protein[1]
      }
      
      if(length(pathogen) > 1){
        pathogen <- pathogen[1]
      }
      
    }
    
    currentDF <- data.frame("ID" = epitopeID[counter],
                            "Pathogen" = pathogen,
                            "Protein" = protein)
    
    otherDF <- rbind(otherDF, currentDF)
    
    counter <- counter + 1
    
  }
  
}

head(sort(table(otherDF$ID), decreasing = T))

#Don't add the protein information unless it is missing
#Add information to peptideInfo
peptideInfo[otherDF$ID,"Pathogen"] <- otherDF$Pathogen
peptideInfo[otherDF$ID,"PathogenComplex"] <- otherDF$Pathogen
peptideInfo[otherDF$ID,"Protein"] <- otherDF$Protein
peptideInfo[otherDF$ID,"BLAST"] <- "Yes"

#### Format pathogen names ####

#To make comparisons easier, I need to make the pathogen names less complex so they are grouped better

#First I convert everything to lower case
peptideInfo$Pathogen <- tolower(peptideInfo$Pathogen)

#Some pathogen names still contain "[" or "]", remove them
peptideInfo$Pathogen <- str_replace_all(peptideInfo$Pathogen, "\\[|\\]", "")

unique_pathogens <- names(table(peptideInfo$Pathogen))

#Define simple groups

#Virus
peptideInfo$Pathogen[str_detect(peptideInfo$Pathogen,"^influenza a virus")] <- "influenza a virus"
peptideInfo$Pathogen[str_detect(peptideInfo$Pathogen,"^influenza b virus")] <- "influenza b virus"
peptideInfo$Pathogen[str_detect(peptideInfo$Pathogen,"^human respiratory syncytial virus")] <- "human respiratory syncytial virus"
peptideInfo$Pathogen[str_detect(peptideInfo$Pathogen,"^human papillomavirus")] <- "human papillomavirus"
peptideInfo$Pathogen[str_detect(peptideInfo$Pathogen,"^human adenovirus")] <- "human adenovirus"
peptideInfo$Pathogen[str_detect(peptideInfo$Pathogen,"^echovirus")] <- "echovirus"
peptideInfo$Pathogen[str_detect(peptideInfo$Pathogen,"^dengue virus")] <- "dengue virus"
peptideInfo$Pathogen[str_detect(peptideInfo$Pathogen,"^coxsackievirus")] <- "coxsackievirus"
peptideInfo$Pathogen[str_detect(peptideInfo$Pathogen,"^rhinovirus a")] <- "rhinovirus a"
peptideInfo$Pathogen[str_detect(peptideInfo$Pathogen,"^rhinovirus b")] <- "rhinovirus b"
peptideInfo$Pathogen[str_detect(peptideInfo$Pathogen,"^poliovirus")] <- "poliovirus"
peptideInfo$Pathogen[str_detect(peptideInfo$Pathogen,"^enterovirus")] <- "enterovirus"
peptideInfo$Pathogen[str_detect(peptideInfo$Pathogen,"^human herpesvirus 1")] <- "human herpesvirus 1"

#Human gamma herpesvirus 4 is EBV
peptideInfo$Pathogen[str_detect(peptideInfo$Pathogen,"^epstein-barr virus")] <- "epstein-barr virus"
peptideInfo$Pathogen[str_detect(peptideInfo$Pathogen,"^human gammaherpesvirus 4")] <- "epstein-barr virus"

#human herpesvirus 5/Human betaherpesvirus 5 is CMV
peptideInfo$Pathogen[str_detect(peptideInfo$Pathogen,"^human herpesvirus 5")] <- "human cytomegalovirus"
peptideInfo$Pathogen[str_detect(peptideInfo$Pathogen,"^human betaherpesvirus 5")] <- "human cytomegalovirus"

#Human alphaherpesvirus 1 is human herpesvirus 1
peptideInfo$Pathogen[str_detect(peptideInfo$Pathogen,"^human alphaherpesvirus 1")] <- "human herpesvirus 1"

#Bacteria
peptideInfo$Pathogen[str_detect(peptideInfo$Pathogen,"^mycobacterium tuberculosis")] <- "mycobacterium tuberculosis"
peptideInfo$Pathogen[str_detect(peptideInfo$Pathogen,"^listeria monocytogenes")] <- "listeria monocytogenes"
peptideInfo$Pathogen[str_detect(peptideInfo$Pathogen,"^limosilactobacillus reuteri")] <- "limosilactobacillus reuteri"
peptideInfo$Pathogen[str_detect(peptideInfo$Pathogen,"^lactobacillus delbrueckii")] <- "lactobacillus delbrueckii"
peptideInfo$Pathogen[str_detect(peptideInfo$Pathogen,"^lactobacillus gasseri")] <- "lactobacillus gasseri"
peptideInfo$Pathogen[str_detect(peptideInfo$Pathogen,"^lactobacillus acidophilus")] <- "lactobacillus acidophilus"
peptideInfo$Pathogen[str_detect(peptideInfo$Pathogen,"^lactiplantibacillus plantarum")] <- "lactiplantibacillus plantarum"
peptideInfo$Pathogen[str_detect(peptideInfo$Pathogen,"^lacticaseibacillus rhamnosus")] <- "lacticaseibacillus rhamnosus"
peptideInfo$Pathogen[str_detect(peptideInfo$Pathogen,"^lachnospiraceae bacterium")] <- "lachnospiraceae bacterium"
peptideInfo$Pathogen[str_detect(peptideInfo$Pathogen,"^haemophilus parainfluenzae")] <- "haemophilus parainfluenzae"
peptideInfo$Pathogen[str_detect(peptideInfo$Pathogen,"^francisella tularensis")] <- "francisella tularensis"
peptideInfo$Pathogen[str_detect(peptideInfo$Pathogen,"^escherichia coli")] <- "escherichia coli"
peptideInfo$Pathogen[str_detect(peptideInfo$Pathogen,"^clostridium sp")] <- "clostridium sp"
peptideInfo$Pathogen[str_detect(peptideInfo$Pathogen,"^campylobacter sp\\.")] <- "campylobacter sp."
peptideInfo$Pathogen[str_detect(peptideInfo$Pathogen,"^borreliella burgdorferi")] <- "borreliella burgdorferi"
peptideInfo$Pathogen[str_detect(peptideInfo$Pathogen,"^bacillus cereus")] <- "bacillus cereus"
peptideInfo$Pathogen[str_detect(peptideInfo$Pathogen,"^vibrio cholerae")] <- "vibrio cholerae"
peptideInfo$Pathogen[str_detect(peptideInfo$Pathogen,"^streptococcus sp\\.")] <- "streptococcus sp."
peptideInfo$Pathogen[str_detect(peptideInfo$Pathogen,"^streptococcus pneumoniae")] <- "streptococcus pneumoniae"
peptideInfo$Pathogen[str_detect(peptideInfo$Pathogen,"^staphylococcus aureus")] <- "staphylococcus aureus"
peptideInfo$Pathogen[str_detect(peptideInfo$Pathogen,"^shigella flexneri")] <- "shigella flexneri"
peptideInfo$Pathogen[str_detect(peptideInfo$Pathogen,"^shigella dysenteriae")] <- "shigella dysenteriae"
peptideInfo$Pathogen[str_detect(peptideInfo$Pathogen,"^salmonella enterica")] <- "salmonella enterica"
peptideInfo$Pathogen[str_detect(peptideInfo$Pathogen,"^prevotella sp\\.")] <- "prevotella sp."
peptideInfo$Pathogen[str_detect(peptideInfo$Pathogen,"^dorea sp\\.")] <- "dorea sp."
peptideInfo$Pathogen[str_detect(peptideInfo$Pathogen,"^alistipes shahii")] <- "alistipes shahii"
peptideInfo$Pathogen[str_detect(peptideInfo$Pathogen,"^anaerobutyricum hallii")] <- "anaerobutyricum hallii"
peptideInfo$Pathogen[str_detect(peptideInfo$Pathogen,"^bacteroides dorei")] <- "bacteroides dorei"
peptideInfo$Pathogen[str_detect(peptideInfo$Pathogen,"^bacteroides massiliensis")] <- "bacteroides massiliensis"
peptideInfo$Pathogen[str_detect(peptideInfo$Pathogen,"^barnesiella intestinihominis")] <- "barnesiella intestinihominis"
peptideInfo$Pathogen[str_detect(peptideInfo$Pathogen,"^clostridiales bacterium")] <- "clostridiales bacterium"
peptideInfo$Pathogen[str_detect(peptideInfo$Pathogen,"^clostridium perfringens")] <- "clostridium perfringens"
peptideInfo$Pathogen[str_detect(peptideInfo$Pathogen,"^coprococcus sp\\.")] <- "coprococcus sp."
peptideInfo$Pathogen[str_detect(peptideInfo$Pathogen,"^dorea formicigenerans")] <- "dorea formicigenerans"
peptideInfo$Pathogen[str_detect(peptideInfo$Pathogen,"^dialister sp\\.")] <- "dialister sp."
peptideInfo$Pathogen[str_detect(peptideInfo$Pathogen,"^eubacterium sp\\.")] <- "eubacterium sp."
peptideInfo$Pathogen[str_detect(peptideInfo$Pathogen,"^firmicutes bacterium")] <- "firmicutes bacterium"
peptideInfo$Pathogen[str_detect(peptideInfo$Pathogen,"^mediterraneibacter gnavus")] <- "mediterraneibacter gnavus"
peptideInfo$Pathogen[str_detect(peptideInfo$Pathogen,"^muribaculaceae bacterium")] <- "muribaculaceae bacterium"
peptideInfo$Pathogen[str_detect(peptideInfo$Pathogen,"^subdoligranulum sp\\.")] <- "subdoligranulum sp."

#prevotella copri is now called segatella copri
peptideInfo$Pathogen[str_detect(peptideInfo$Pathogen,"^prevotella copri")] <- "segatella copri"


#Parasite
peptideInfo$Pathogen[str_detect(peptideInfo$Pathogen,"^toxoplasma gondii")] <- "toxoplasma gondii"
peptideInfo$Pathogen[str_detect(peptideInfo$Pathogen,"^plasmodium falciparum")] <- "plasmodium falciparum"

#Only defined as plasmodium reichenowi through blast as the plasmodium falciparum scientific names have got strain information on them
peptideInfo$Pathogen[str_detect(peptideInfo$Pathogen,"^plasmodium reichenowi")] <- "plasmodium falciparum"


#Produces butyrate/good microbiome 
peptideInfo$Pathogen[str_detect(peptideInfo$Pathogen,"^faecalibacterium sp\\.")] <- "faecalibacterium sp."
peptideInfo$Pathogen[str_detect(peptideInfo$Pathogen,"^blautia sp\\.")] <- "blautia sp."
peptideInfo$Pathogen[str_detect(peptideInfo$Pathogen,"^bacteroides uniformis")] <- "bacteroides uniformis"
peptideInfo$Pathogen[str_detect(peptideInfo$Pathogen,"^bacteroides thetaiotaomicron")] <- "bacteroides thetaiotaomicron"
peptideInfo$Pathogen[str_detect(peptideInfo$Pathogen,"^bacteroides sp\\.")] <- "bacteroides sp."
peptideInfo$Pathogen[str_detect(peptideInfo$Pathogen,"^bacteroides ovatus")] <- "bacteroides ovatus"
peptideInfo$Pathogen[str_detect(peptideInfo$Pathogen,"^bacteroides fragilis")] <- "bacteroides fragilis"
peptideInfo$Pathogen[str_detect(peptideInfo$Pathogen,"^roseburia sp\\.")] <- "roseburia sp."
peptideInfo$Pathogen[str_detect(peptideInfo$Pathogen,"^roseburia intestinalis")] <- "roseburia intestinalis"
peptideInfo$Pathogen[str_detect(peptideInfo$Pathogen,"^parabacteroides sp\\.")] <- "parabacteroides sp."
peptideInfo$Pathogen[str_detect(peptideInfo$Pathogen,"^parabacteroides distasonis")] <- "parabacteroides distasonis"
peptideInfo$Pathogen[str_detect(peptideInfo$Pathogen,"^agathobacter rectalis")] <- "agathobacter rectalis"
peptideInfo$Pathogen[str_detect(peptideInfo$Pathogen,"^anaerostipes sp\\.")] <- "anaerostipes sp."
peptideInfo$Pathogen[str_detect(peptideInfo$Pathogen,"^bifidobacterium adolescentis")] <- "bifidobacterium adolescentis"
peptideInfo$Pathogen[str_detect(peptideInfo$Pathogen,"^bifidobacterium bifidum")] <- "bifidobacterium bifidum"
peptideInfo$Pathogen[str_detect(peptideInfo$Pathogen,"^blautia obeum")] <- "blautia obeum"
peptideInfo$Pathogen[str_detect(peptideInfo$Pathogen,"^clostridium scindens")] <- "clostridium scindens"
peptideInfo$Pathogen[str_detect(peptideInfo$Pathogen,"^coprococcus comes")] <- "coprococcus comes"
peptideInfo$Pathogen[str_detect(peptideInfo$Pathogen,"^lactobacillus plantarum")] <- "lactobacillus plantarum"

#Maybe good microbiome
peptideInfo$Pathogen[str_detect(peptideInfo$Pathogen,"^ruminococcus sp\\.")] <- "ruminococcus sp."
peptideInfo$Pathogen[str_detect(peptideInfo$Pathogen,"^ruminococcus torques")] <- "ruminococcus torques"
peptideInfo$Pathogen[str_detect(peptideInfo$Pathogen,"^ruminococcus lactaris")] <- "ruminococcus lactaris"
peptideInfo$Pathogen[str_detect(peptideInfo$Pathogen,"^phocaeicola dorei")] <- "phocaeicola dorei"
peptideInfo$Pathogen[str_detect(peptideInfo$Pathogen,"^akkermansia muciniphila")] <- "akkermansia muciniphila"
peptideInfo$Pathogen[str_detect(peptideInfo$Pathogen,"^alistipes putredinis")] <- "alistipes putredinis"
peptideInfo$Pathogen[str_detect(peptideInfo$Pathogen,"^oscillibacter sp\\.")] <- "oscillibacter sp."
peptideInfo$Pathogen[str_detect(peptideInfo$Pathogen,"^oscillospiraceae bacterium")] <- "oscillospiraceae bacterium"

#bacteroides vulgatus is the old name for phocaeicola vulgatus
peptideInfo$Pathogen[str_detect(peptideInfo$Pathogen,"^bacteroides vulgatus")] <- "phocaeicola vulgatus"
peptideInfo$Pathogen[str_detect(peptideInfo$Pathogen,"^phocaeicola vulgatus")] <- "phocaeicola vulgatus"

peptideInfo$Pathogen[str_detect(peptideInfo$Pathogen,"^enterococcus faecalis")] <- "enterococcus faecalis"


#Phage
peptideInfo$Pathogen[str_detect(peptideInfo$Pathogen,"^lactobacillus phage")] <- "lactobacillus phage"
peptideInfo$Pathogen[str_detect(peptideInfo$Pathogen,"^lactococcus phage")] <- "lactococcus phage"
peptideInfo$Pathogen[str_detect(peptideInfo$Pathogen,"^enterococcus phage")] <- "enterococcus phage"
peptideInfo$Pathogen[str_detect(peptideInfo$Pathogen,"^enterobacteria phage")] <- "enterobacteria phage"
peptideInfo$Pathogen[str_detect(peptideInfo$Pathogen,"^clostridium phage")] <- "clostridium phage"
peptideInfo$Pathogen[str_detect(peptideInfo$Pathogen,"^campylobacter phage")] <- "campylobacter phage"
peptideInfo$Pathogen[str_detect(peptideInfo$Pathogen,"^staphylococcus phage")] <- "staphylococcus phage"
peptideInfo$Pathogen[str_detect(peptideInfo$Pathogen,"^shigella phage")] <- "shigella phage"
peptideInfo$Pathogen[str_detect(peptideInfo$Pathogen,"^salmonella phage")] <- "salmonella phage"
peptideInfo$Pathogen[str_detect(peptideInfo$Pathogen,"^klebsiella phage")] <- "klebsiella phage"
peptideInfo$Pathogen[str_detect(peptideInfo$Pathogen,"^escherichia phage")] <- "escherichia phage"

unique_pathogensPost <- names(table(peptideInfo$Pathogen))

#If Protein_uniref values exists, make the Protein column equal to Protein_uniref
peptideInfo$Protein[which(peptideInfo$Protein_uniref != "")] <- peptideInfo$Protein_uniref[which(peptideInfo$Protein_uniref != "")]

#Convert protein names to lower case
peptideInfo$Protein <- tolower(peptideInfo$Protein)

#### Save peptideInfo ####
peptideInfo[,"epitopeGroup"] <- paste0(peptideInfo$Protein, "_", peptideInfo$Pathogen)
peptideInfo[,"epitopeUnique"] <- paste0(peptideInfo$Protein, "_", peptideInfo$Pathogen, "_", row.names(peptideInfo))

#For ease of analysis, remove all spaces and punctuation 

write.csv(peptideInfo, "PhIPseq_peptideInfo_processed_BLAST.csv")



#### outtakes ####
for(i in list.files("BLAST/output/", pattern = "other")){
  print(i)
  
  #Load in BLASTed twistEpitope data
  otherBLAST <- read_xml(paste0("BLAST/output/", i))
  
  #Extract the IDs and the BLAST results
  epitopeID <- xml_text(xml_find_all(otherBLAST, ".//Iteration_query-def"))
  epitope_poss <- xml_text(xml_find_all(otherBLAST, ".//Hit_def"))
  epitope_hitNum <- xml_integer(xml_find_all(otherBLAST, ".//Hit_num"))
  epitopeIteration_hits <- xml_text(xml_find_all(otherBLAST, ".//Iteration_hits"))
  
  epitope_firstHit <- which(epitope_hitNum == 1)
  
  counter <- 1
  
  if(max(table(epitope_hitNum)) != length(epitopeID)){stop()}
  
  #This code works on the assumption, that the most likely hit (in terms of E value) is displayed first
  
  #There are times where there isn't 100 matches. 
  #Further complicating things is when an amino acid sequence does not have a result as it shows up as a query name but doesn't have any results
  
  #If a result was not returned for an epitope, the value of epitopeIteration_hits will be "\n"
  
  for(j in epitope_firstHit){
    
    #Some of the results are formatted differently. 
    #I try and clean them down in this section of code
    
    currentEpitope <- epitope_poss[j]
    
    #If there are multiple options, choose one
    if(str_detect(currentEpitope, "\\>")){
      currentEpitope <- unlist(str_split(currentEpitope, " \\>"))
      currentEpitope <- currentEpitope[c(TRUE, rep(FALSE, (length(currentEpitope) -1) ))]
    }
    
    if(str_detect(currentEpitope, "RecName")){
      currentEpitope <- str_replace_all(currentEpitope, "RecName: Full=", "")
      
      #Some of these don't have a semicolon
      if(str_detect(currentEpitope, ";")){
        semiColonSplit <- unlist(str_split(currentEpitope, ";"))
        
        currentEpitope_protein <- semiColonSplit[c(TRUE, rep(FALSE, (length(semiColonSplit) -1) ))]
        currentEpitope_pathogen <- unlist(str_split(currentEpitope, " \\["))[c(FALSE, TRUE)]
        
        currentEpitope <- paste0(currentEpitope_protein, " [", currentEpitope_pathogen)
        
      }
      
      rm(currentEpitope_protein)
      rm(currentEpitope_pathogen)
      gc()
    }
    
    #extract the protein and pathogen information
    infoSplit <- unlist(str_split(currentEpitope, " \\["))
    
    #Some can have two or more square bracket sections. This can confuse things
    protein <- infoSplit[c(TRUE, FALSE, rep(FALSE, (length(infoSplit)-2) ) )]
    pathogen <- infoSplit[c(FALSE, TRUE, rep(FALSE, (length(infoSplit)-2) ) )]
    
    #Remove the trailing ] from pathogen
    pathogen <- str_replace_all(pathogen, "\\]", "")
    
    currentDF <- data.frame("ID" = epitopeID[counter],
                            "Pathogen" = pathogen,
                            "Protein" = protein)
    
    # if(str_detect(pathogen, "Priestia")){stop()}
    
    otherDF <- rbind(otherDF, currentDF)
    
    counter <- counter + 1
    
  }
  
}


