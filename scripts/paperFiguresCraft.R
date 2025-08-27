library(dplyr)
library(ggplot2)
library(stringr)
library(ggpubr)
library(rstatix)
library(readxl)
library(reshape2)
library(resample)
library(Seurat)
library(gridExtra)
library("FactoMineR")
library("factoextra")
library(ggtext)
library(scales) 


#For the cochran trend test
library(DescTools)

library(wordcloud2)

library(factoextra)

library(vegan)

library(emmeans)

#Correlation plots
library(PerformanceAnalytics)

#library(beer)

#######
#Perhaps there are people who are more enriched for certain type of pathogens and thus have immune frequencies attuned towards combating those?
#######

#Surely its the current location that matters? That will determine what pathogens they've been exposed to

setwd("C:/Users/rflaidlaw/Documents/CapTan/AnalysisV2/PhIPseq/")

microbiomePaperResUrbanization <- read.delim("C:/Users/rflaidlaw/Downloads/media-5.csv")

firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

area_palette <- c("Rural Senegalese" = "#e59f01", "Urban Senegalese" = "#54b3e8", "Urban Dutch" = "#009e72")
areaComplex_palette <- c("Rural Senegalese" = "#e59f01", "Semi-Urban Senegalese" = "#9a5b90" ,"Urban Senegalese" = "#54b3e8", "Urban Dutch" = "#009e72")
direction_palette <- c("RUR SEN" = "#e59f01", "URB NLD" = "#009e72")

metaData <- as.data.frame(read.csv("phipseq_metadata.csv"))
row.names(metaData) <- metaData$ID

metaData_complex <- as.data.frame(read.csv("phipseq_metadata_extended.csv", row.names = 1))
row.names(metaData_complex) <- metaData_complex$ID

#Make age groupings
#0-17, 18-29, 30-45, 46+

# phipSeq_FC <- read.csv("res_foldchange_annotated.csv", row.names = 2)
phipSeq_count <- read.csv("res_counts_annotated.csv", row.names = 2)
sampleNames <- colnames(phipSeq_count)[8:ncol(phipSeq_count)]
phipSeq_count <- t(as.matrix( phipSeq_count[, sampleNames]) )
#Change the column names to reflect the metaData
row.names(phipSeq_count) <- unlist(str_split( row.names(phipSeq_count), "_" ))[c(FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE)]

phipSeq_p <- read.csv("res_padj_annotated.csv", row.names = 2)

peptideInfo <- read.csv("PhIPseq_peptideInfo_processed_BLAST.csv", row.names = 1)

phipSeq_p <- read.csv("data/individualAll_epitope_binary.csv", row.names = 1)
phipSeq_p_freq_residence <- read.csv("data/residenceAll_epitope_percentage.csv", row.names = 1)
phipSeq_p_freq_sex <- read.csv("data/sexSenegal_epitope_percentage.csv", row.names = 1)
phipSeq_p_freq_sexResidence <- read.csv("data/sexResidenceAll_epitope_percentage.csv", row.names = 1)
phipSeq_p_pathogenTotal <- read.csv("data/pathogenEpitopeTotals.csv", row.names = 1)

senegalIndividuals <- metaData[str_detect(metaData$Residence, "enegal"),"ID"]

#Number of organisms
length(unique(peptideInfo$Pathogen))

#### Filter features old ####

phipSeq_p_df <- as.data.frame(phipSeq_p)

phipSeq_p_df$Residence <- metaData[row.names(phipSeq_p_df), "Residence"]
phipSeq_p_df$Sex <- metaData[row.names(phipSeq_p_df), "Sex"]
phipSeq_p_df$SexResidence <- paste0(phipSeq_p_df$Sex, "_", phipSeq_p_df$Residence)

epitopeByResidence <- phipSeq_p_df[, -which(colnames(phipSeq_p_df) %in% c("Sex", "SexResidence"))]  %>%
  dplyr::group_by(Residence) %>%
  summarise(across(everything(), ~ sum(.x)))

epitopeByResidence <- as.data.frame(epitopeByResidence)

row.names(epitopeByResidence) <- epitopeByResidence$Residence
epitopeByResidence$Residence <- NULL

epitopeByResidence["Rural Senegalese",] <- (epitopeByResidence["Rural Senegalese",] / nrow(subset(metaData, Residence == "Rural Senegalese"))) * 100
epitopeByResidence["Urban Senegalese",] <- (epitopeByResidence["Urban Senegalese",] / nrow(subset(metaData, Residence == "Urban Senegalese"))) * 100
epitopeByResidence["Urban Dutch",] <- (epitopeByResidence["Urban Dutch",] / nrow(subset(metaData, Residence == "Urban Dutch"))) * 100

#Create copy for only senegalese
epitopeByResidenceSenegal <- epitopeByResidence[str_detect(row.names(epitopeByResidence), "enegal"),]

#This code is quite convoluted, but it basically just finds what epitopes are present in 25% of individuals in at least one residence group
epitopeByResidence <- (epitopeByResidence + 0.000000000000000000000000000001) / 25
epitopeByResidence <- log(epitopeByResidence)
epitopeByResidence <- sign(epitopeByResidence)
epitopeByResidence <- colSums(epitopeByResidence)
epitopeByResidence <- epitopeByResidence[epitopeByResidence !=-3 ]


epitopeByResidenceSenegal <- (epitopeByResidenceSenegal + 0.000000000000000000000000000001) / 25
epitopeByResidenceSenegal <- log(epitopeByResidenceSenegal)
epitopeByResidenceSenegal <- sign(epitopeByResidenceSenegal)
epitopeByResidenceSenegal <- colSums(epitopeByResidenceSenegal)
epitopeByResidenceSenegal <- epitopeByResidenceSenegal[epitopeByResidenceSenegal != -2]

#Do the same but for Sex
epitopeBySexSenegal <- phipSeq_p_df[senegalIndividuals, -which(colnames(phipSeq_p_df) %in% c("Residence", "SexResidence"))  ] %>%
  dplyr::select(contains( c("agilent","twist","coron","Sex") )) %>%
  dplyr::group_by(Sex) %>%
  summarise(across(everything(), ~ sum(.x)))

epitopeBySexSenegal <- as.data.frame(epitopeBySexSenegal)

row.names(epitopeBySexSenegal) <- epitopeBySexSenegal$Sex
epitopeBySexSenegal$Sex <- NULL

epitopeBySexSenegal["M",] <- (epitopeBySexSenegal["M",] / nrow(subset(metaData, Residence != "Urban Dutch" & Sex == "M"))) * 100
epitopeBySexSenegal["F",] <- (epitopeBySexSenegal["F",] / nrow(subset(metaData, Residence != "Urban Dutch"& Sex == "F"))) * 100

epitopeBySexSenegal <- (epitopeBySexSenegal + 0.000000000000000000000000000001) / 25
epitopeBySexSenegal <- log(epitopeBySexSenegal)
epitopeBySexSenegal <- sign(epitopeBySexSenegal)
epitopeBySexSenegal <- colSums(epitopeBySexSenegal)
epitopeBySexSenegal <- epitopeBySexSenegal[epitopeBySexSenegal != -2]

#Finally, for combination of sex and residence
epitopeBySexResidence <- phipSeq_p_df[, -which(colnames(phipSeq_p_df) %in% c("Residence", "Sex"))  ]  %>%
  dplyr::group_by(SexResidence) %>%
  summarise(across(everything(), ~ sum(.x)))

epitopeBySexResidence <- as.data.frame(epitopeBySexResidence)

row.names(epitopeBySexResidence) <- epitopeBySexResidence$SexResidence
epitopeBySexResidence$SexResidence <- NULL

epitopeBySexResidence["F_Rural Senegalese",] <- (epitopeBySexResidence["F_Rural Senegalese",] / nrow(subset(metaData, Residence == "Rural Senegalese" & Sex == "F"))) * 100
epitopeBySexResidence["M_Rural Senegalese",] <- (epitopeBySexResidence["M_Rural Senegalese",] / nrow(subset(metaData, Residence == "Rural Senegalese" & Sex == "M"))) * 100

epitopeBySexResidence["F_Urban Senegalese",] <- (epitopeBySexResidence["F_Urban Senegalese",] / nrow(subset(metaData, Residence == "Urban Senegalese" & Sex == "F"))) * 100
epitopeBySexResidence["M_Urban Senegalese",] <- (epitopeBySexResidence["M_Urban Senegalese",] / nrow(subset(metaData, Residence == "Urban Senegalese" & Sex == "M"))) * 100

epitopeBySexResidence["F_Urban Dutch",] <- (epitopeBySexResidence["F_Urban Dutch",] / nrow(subset(metaData, Residence == "Urban Dutch" & Sex == "F"))) * 100
epitopeBySexResidence["M_Urban Dutch",] <- (epitopeBySexResidence["M_Urban Dutch",] / nrow(subset(metaData, Residence == "Urban Dutch" & Sex == "M"))) * 100

#As we are dealing with smaller groups, i'm making the percentage threshold higher
epitopeBySexResidence <- (epitopeBySexResidence + 0.000000000000000000000000000001) / 35
epitopeBySexResidence <- log(epitopeBySexResidence)
epitopeBySexResidence <- sign(epitopeBySexResidence)
epitopeBySexResidence <- colSums(epitopeBySexResidence)
epitopeBySexResidence <- epitopeBySexResidence[epitopeBySexResidence != -6]

#Find ones specific for rural and urban senegalese
epitopeBySexRuralSen <- subset(phipSeq_p_df, Residence == "Rural Senegalese")[, -which(colnames(phipSeq_p_df) %in% c("Residence", "Sex"))  ]  %>%
  dplyr::group_by(SexResidence) %>%
  summarise(across(everything(), ~ sum(.x)))

epitopeBySexRuralSen <- as.data.frame(epitopeBySexRuralSen)

row.names(epitopeBySexRuralSen) <- epitopeBySexRuralSen$SexResidence
epitopeBySexRuralSen$SexResidence <- NULL

epitopeBySexRuralSen["F_Rural Senegalese",] <- (epitopeBySexRuralSen["F_Rural Senegalese",] / nrow(subset(metaData, Residence == "Rural Senegalese" & Sex == "F"))) * 100
epitopeBySexRuralSen["M_Rural Senegalese",] <- (epitopeBySexRuralSen["M_Rural Senegalese",] / nrow(subset(metaData, Residence == "Rural Senegalese" & Sex == "M"))) * 100

#As we are dealing with smaller groups, i'm making the percentage threshold higher
epitopeBySexRuralSen <- (epitopeBySexRuralSen + 0.000000000000000000000000000001) / 35
epitopeBySexRuralSen <- log(epitopeBySexRuralSen)
epitopeBySexRuralSen <- sign(epitopeBySexRuralSen)
epitopeBySexRuralSen <- colSums(epitopeBySexRuralSen)
epitopeBySexRuralSen <- epitopeBySexRuralSen[epitopeBySexRuralSen != -2]

epitopeBySexUrbanSen <- subset(phipSeq_p_df, Residence == "Urban Senegalese")[, -which(colnames(phipSeq_p_df) %in% c("Residence", "Sex"))  ]  %>%
  dplyr::group_by(SexResidence) %>%
  summarise(across(everything(), ~ sum(.x)))

epitopeBySexUrbanSen <- as.data.frame(epitopeBySexUrbanSen)

row.names(epitopeBySexUrbanSen) <- epitopeBySexUrbanSen$SexResidence
epitopeBySexUrbanSen$SexResidence <- NULL

epitopeBySexUrbanSen["F_Urban Senegalese",] <- (epitopeBySexUrbanSen["F_Urban Senegalese",] / nrow(subset(metaData, Residence == "Urban Senegalese" & Sex == "F"))) * 100
epitopeBySexUrbanSen["M_Urban Senegalese",] <- (epitopeBySexUrbanSen["M_Urban Senegalese",] / nrow(subset(metaData, Residence == "Urban Senegalese" & Sex == "M"))) * 100

#As we are dealing with smaller groups, i'm making the percentage threshold higher
epitopeBySexUrbanSen <- (epitopeBySexUrbanSen + 0.000000000000000000000000000001) / 35
epitopeBySexUrbanSen <- log(epitopeBySexUrbanSen)
epitopeBySexUrbanSen <- sign(epitopeBySexUrbanSen)
epitopeBySexUrbanSen <- colSums(epitopeBySexUrbanSen)
epitopeBySexUrbanSen <- epitopeBySexUrbanSen[epitopeBySexUrbanSen != -2]


#### Define objects ####

phipSeq_p_raw_df <- phipSeq_p
phipSeq_p_freq_residence_raw_df <- phipSeq_p_freq_residence
phipSeq_p_freq_sex_raw_df <- phipSeq_p_freq_sex
phipSeq_p_freq_sexResidence_raw_df <- phipSeq_p_freq_sexResidence

phipSeq_p_raw <- as.matrix(phipSeq_p)
phipSeq_p_freq_residence_raw <- as.matrix(phipSeq_p_freq_residence)
phipSeq_p_freq_sex_raw <- as.matrix(phipSeq_p_freq_sex)
phipSeq_p_freq_sexResidence_raw <- as.matrix(phipSeq_p_freq_sexResidence)

phipSeq_p_df <- phipSeq_p_raw_df[,names(epitopeByResidence)]
phipSeq_p_freq_residence_df <- phipSeq_p_freq_residence_raw_df[,names(epitopeByResidence)]
phipSeq_p_freq_sex_df <- phipSeq_p_freq_sex_raw_df[,names(epitopeByResidence)]
phipSeq_p_freq_sexResidence_df <- phipSeq_p_freq_sexResidence_raw_df[,names(epitopeByResidence)]

phipSeq_p <- as.matrix(phipSeq_p_df)
phipSeq_p_freq_residence <- as.matrix(phipSeq_p_freq_residence_df)
phipSeq_p_freq_sex <- as.matrix(phipSeq_p_freq_sex_df)
phipSeq_p_freq_sexResidence <- as.matrix(phipSeq_p_freq_sexResidence_df)

#### Investigate the epitopes ####

#Look at the epitopes and categorise based on how many individuals they are present in
#Bins: 1, 2-3, 4-9, 10-24, 25-49, 50+

epitopeCount <- colSums(phipSeq_p_raw)

head(sort(epitopeCount, decreasing = T))

epitopeCountGrouping_key <- c("0" = "1", "1" = "2-3", "2" = "4-9", "3" = "10-24", "4" = "25-49", "5" = "50+")

epitopeCount_groups <- epitopeCount

epitopeCount_groups[epitopeCount <= 1] <- 0
epitopeCount_groups[epitopeCount >=2 & epitopeCount <= 3] <- 1
epitopeCount_groups[epitopeCount >=4 & epitopeCount <= 9] <- 2
epitopeCount_groups[epitopeCount >=10 & epitopeCount <= 24] <- 3
epitopeCount_groups[epitopeCount >=25 & epitopeCount <=49 ] <- 4
epitopeCount_groups[epitopeCount >=50 ] <- 5

epitopeCount_groups <- as.character(epitopeCount_groups)

epitopeCount_groups <- epitopeCountGrouping_key[epitopeCount_groups]

plotDF <- data.frame(table(epitopeCount_groups))

plotDF$epitopeCount_groups <- factor(plotDF$epitopeCount_groups,
                                     levels = c("1", "2-3","4-9", "10-24", "25-49", "50+"))

ggOut <- ggplot(plotDF, aes(x = epitopeCount_groups, y = Freq, fill = epitopeCount_groups)) +
  geom_bar(stat="identity") +
  theme_minimal() +
  xlab("No. individuals in which anti-epitope is present") +
  ylab("Number of anti-epitopes") 
ggOut
ggsave(filename = "//vf-lucid-r-o.lumcnet.prod.intern/lucid-r-o$/Projects/Senegal_phipseq/plots/epitopeCountGrouping.pdf", plot = ggOut, 
       width =12, height = 8)

plotDF <- data.frame("epitopeCount" = epitopeCount)
ggOut <- ggplot(plotDF, aes(x = epitopeCount)) +
  geom_histogram() +
  theme_minimal() +
  xlab("No. individuals in which anti-epitope is present") + 
  ylab("No. of significant anti-epitopes (log10)")+
  scale_y_log10(guide = "axis_logticks")
ggOut
ggsave(filename = "//vf-lucid-r-o.lumcnet.prod.intern/lucid-r-o$/Projects/Senegal_phipseq/plots/epitopeCount_histogram.pdf", plot = ggOut, 
       width =9, height = 8)

# #Save as stacked barplot instead (use percentages)
# plotDF$percentage <- plotDF$Freq
# plotDF$percentage <- (plotDF$percentage / sum(plotDF$Freq)) * 100
# 
# ggOut <- ggplot(plotDF, aes(x = "hi", y = percentage, fill = epitopeCount_groups)) +
#   geom_bar(stat="identity", position = "stack") +
#   theme_minimal() +
#   xlab("") +
#   ylab("Percentage ")
# ggOut
# ggsave(plot = ggOut, filename = "plots/figure_components/epitopePercentGrouping.pdf")

#Five peptides are present in 61/62 individuals (the highest amount):
#agilent_133222 - staphylococcal protein A (Ig binding protein), Staphylococcus aureus
#agilent_236864 - staphylococcal protein A (Ig binding protein), Staphylococcus aureus
#agilent_560 - Pneumococcal histidine triad protein D, Streptococcus pneumoniae
#agilent_7538 - pneumococcal histidine triad protein E & hydrolase & HIT family hydrolase, Streptococcus pneumoniae
#twist_47588 - Attachment glycoprotein, Human respiratory syncytial virus

mostCommonEpitopes <- c("agilent_133222", "agilent_236864", "agilent_560", "agilent_7538", "twist_47588")

#Look at the epitopes that are present in 50 and more individuals
mostPrevalentEpitopes <- epitopeCount[epitopeCount >=50]

#look at the rarest epitopes
mostRareEpitopes <- epitopeCount[epitopeCount ==1]


#Word cloud the most frequent and rarest epitopes
commonEpitopes <- table(peptideInfo[names(mostPrevalentEpitopes), "Pathogen"])
rareEpitopes <- table(peptideInfo[names(mostRareEpitopes), "Pathogen"])

#Normalise for how many times an epitope of that organism is present in the dataset
pathogenTotals <- table(peptideInfo$Pathogen)

#If a pathogen is only present 2 or fewer times, don't investigate it
keepPathogen <- pathogenTotals[pathogenTotals > 2]

commonEpitopes <- commonEpitopes / pathogenTotals[names(commonEpitopes)]
rareEpitopes <- rareEpitopes / pathogenTotals[names(rareEpitopes)]

commonEpitopes <- commonEpitopes[intersect(names(keepPathogen), names(commonEpitopes))]

wordcloud2(commonEpitopes, size = 0.2)


wordcloud2(rareEpitopes, size = 0.8)


rm(plotDF)
rm(epitopeCount_groups)
rm(epitopeCount)
rm(epitopeCountGrouping_key)
rm(rareEpitopes)
rm(commonEpitopes)
rm(keepPathogen)
rm(pathogenTotals)
rm(mostRareEpitopes)
gc()

#### Get normalised total epitopes for different pathogens ####

#row names are individual's ID.
#columns are total epitope count followed by the normalised epitope count for certain pathogens

normTotalEpitope_DF <- data.frame(epitopeTotals = rowSums(phipSeq_p_raw[metaData$ID,]),
                                  row.names = metaData$ID)

pathogens <- c("shigella flexneri", "rhinovirus b", "human adenovirus", "epstein-barr virus",
               "homo sapiens", "aureus", "streptococcus pneumoniae", "streptococcus pyogenes",
               "human respiratory syncytial virus", "barnesiella intestinihominis", "gammaherpesvirus 8",
               "plasmodium reichenowi", "human cytomegalovirus")

pathogens <- unique(peptideInfo$Pathogen)

for(i in pathogens){
  
  peptideInfoCurrent <- peptideInfo[peptideInfo$Pathogen == i,]
  
  currentEpitopes <- row.names(peptideInfoCurrent)
  
  
  
  #Add current epitopes bound to data frame
  if(length(currentEpitopes)==1){
    currentEpitopesTotals <- phipSeq_p_raw[,currentEpitopes]
    
  }
  
  
  else{
    currentEpitopesTotals <- rowSums(phipSeq_p_raw[,currentEpitopes])
    
  }
  
  #Sort name 
  i <- str_replace_all(i, "-| |\\/", "")
  
  normTotalEpitope_DF[[paste0(i, "TotalNorm")]] <- currentEpitopesTotals[row.names(normTotalEpitope_DF)] / normTotalEpitope_DF$epitopeTotals
  
}

normTotalEpitope_DF$Residence <- metaData[row.names(normTotalEpitope_DF), "Residence"]
normTotalEpitope_DF$Residence <- factor(normTotalEpitope_DF$Residence, levels = c("Rural Senegalese", "Urban Senegalese", "Urban Dutch"))
normTotalEpitope_DF$Donor <- row.names(normTotalEpitope_DF)

test_normTotal <- normTotalEpitope_DF
test_normTotal <- test_normTotal[,str_detect(colnames(test_normTotal), "TotalNorm")]

rowSums(test_normTotal)
#### Get normalised total epitopes for different proteins ####

#row names are individual's ID.
#columns are total epitope count followed by the normalised epitope count for certain pathogens

normTotalEpitope_DF_protein <- data.frame(epitopeTotals = rowSums(phipSeq_p_raw[metaData$ID,]),
                                  row.names = metaData$ID)

peptideInfo$protein_pathogen <- peptideInfo$epitopeGroup

peptideInfo$protein_pathogen <- str_replace_all(peptideInfo$protein_pathogen, "-| ", "")

proteins <- unique(peptideInfo$protein_pathogen)

for(i in proteins){
  
  peptideInfoCurrent <- peptideInfo[peptideInfo$protein_pathogen == i,]
  
  currentEpitopes <- row.names(peptideInfoCurrent)
  
  #Add current epitopes bound to data frame
  if(length(currentEpitopes)==1){
    currentEpitopesTotals <- phipSeq_p_raw[,currentEpitopes]
    
  }
  
  
  else{
    currentEpitopesTotals <- rowSums(phipSeq_p_raw[,currentEpitopes])
    
  }
  
  normTotalEpitope_DF_protein[[paste0(i, "TotalNorm")]] <- currentEpitopesTotals[row.names(normTotalEpitope_DF_protein)] / normTotalEpitope_DF_protein$epitopeTotals
  
}

normTotalEpitope_DF_protein$Residence <- metaData[row.names(normTotalEpitope_DF_protein), "Residence"]
normTotalEpitope_DF_protein$Residence <- factor(normTotalEpitope_DF_protein$Residence, levels = c("Rural Senegalese", "Urban Senegalese", "Urban Dutch"))
normTotalEpitope_DF_protein$Donor <- row.names(normTotalEpitope_DF_protein)

test_normTotal <- normTotalEpitope_DF_protein
test_normTotal <- test_normTotal[,str_detect(colnames(test_normTotal), "TotalNorm")]

rowSums(test_normTotal)

#### Investigate rural, urban and dutch: group level antibody detected ####

# epitopeDetectedGroup <- phipSeq_p_raw_df
# 
# epitopeDetectedGroup$Residence <- metaData[row.names(epitopeDetectedGroup), "Residence"]
# 
# #iterate three times and select different samples of rural and urban senegalese of equal size to urban dutch
# sampleSize <- sum(str_detect(epitopeDetectedGroup$Residence, "utch"))
# 
# 
# # epitopeDetectedGroup <- epitopeDetectedGroup %>%
# #                                 group_by(Residence) %>%
# #                                 summarise(across(everything(), sum))
# 
# inputList <- list()
# 
# for(i in c(42, 64, 1)){
#   
#   set.seed(i)
#   
#   ruralSenegal <- subset(epitopeDetectedGroup, Residence == "Rural Senegalese")
#   ruralSenegal <- ruralSenegal[sample(row.names(ruralSenegal), size = sampleSize),]
#   
#   urbanSenegal <- subset(epitopeDetectedGroup, Residence == "Urban Senegalese")
#   urbanSenegal <- urbanSenegal[sample(row.names(urbanSenegal), size = sampleSize),]
#   
#   urbanDutch <- subset(epitopeDetectedGroup, Residence == "Urban Dutch")
#   
#   currentCombination <- rbind(ruralSenegal, urbanSenegal)
#   currentCombination <- rbind(currentCombination, urbanDutch)
#   
#   inputList[[as.character(i)]] <- currentCombination
#   
#   currentDetected <- currentCombination %>%
#                                   group_by(Residence) %>%
#                                   summarise(across(everything(), sum))
#   currentDetected$total <- rep(sampleSize, 3)
#   
#   currentDetected <- as.data.frame(currentDetected)
#   
#   row.names(currentDetected) <- currentDetected$Residence
#   currentDetected$Residence <- NULL
#   
#   plotDF <- data.frame(Residence = row.names(currentDetected),
#                        total = ncol(currentDetected))
#   
#   epitopeCounts <- rowSums(currentDetected)
#   plotDF$count <- epitopeCounts[plotDF$Residence]
#   plotDF$percent <- (plotDF$count / plotDF$total ) * 100
#   
#   plotDF$Residence <- factor(plotDF$Residence, levels = c("Rural Senegalese", "Urban Senegalese", "Urban Dutch"))
#   
#   print(ggplot(plotDF, aes(x = Residence, y = percent, fill = Residence)) +
#           geom_bar(stat = "identity") +
#           ylim(0,100) +
#           scale_fill_manual(values = area_palette)
#   )
#   
#   
# }


#### Investigate rural, urban and dutch: PCA and cochran-armitage test (epitope presence/absence as input) ####

#Get the number of different organisms bound for each individual
diffOrganisms <- reshape2::melt(phipSeq_p_raw)
diffOrganisms$Pathogen <- peptideInfo[as.character(diffOrganisms$Var2), "Pathogen"]

diffOrganisms <- subset(diffOrganisms, value == 1)

diffOrganisms <- diffOrganisms %>%
                    group_by(Var1) %>%
                  summarise(nOrganisms = length(unique(Pathogen)),
                            nEpitope = length(Pathogen))

diffOrganisms$Residence <- metaData[as.character(diffOrganisms$Var1), "Residence"]
diffOrganisms$Residence <- factor(diffOrganisms$Residence, levels = c("Rural Senegalese", "Urban Senegalese", "Urban Dutch") )

diffOrganisms$norm_nOrganisms <- diffOrganisms$nOrganisms / diffOrganisms$nEpitope

ggplot(diffOrganisms, aes(x = Residence, y = nOrganisms, fill = Residence)) +
                          geom_boxplot() +
                          geom_jitter() +
                          scale_fill_manual(values = area_palette)

ggplot(diffOrganisms, aes(x = Residence, y = nEpitope, fill = Residence)) +
  geom_boxplot() +
  geom_jitter() +
  scale_fill_manual(values = area_palette)

my_comparisons <- list( c("Rural Senegalese", "Urban Senegalese"), c("Urban Senegalese", "Urban Dutch"), c("Rural Senegalese", "Urban Dutch") )
ggOut <- ggboxplot(diffOrganisms, x = "Residence", y = "nEpitope",
          fill = "Residence")+ 
  stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 50) +
  scale_fill_manual(values = area_palette) +
  geom_jitter()
ggsave(filename = "//vf-lucid-r-o.lumcnet.prod.intern/lucid-r-o$/Projects/Senegal_phipseq/plots/nEpitopes_residence_comparison.pdf", plot = ggOut, 
       width =12, height = 8)

ggplot(diffOrganisms, aes(x = nOrganisms, y = nEpitope, color = Residence)) +
  geom_point() +
  geom_jitter() +
  scale_color_manual(values = area_palette)

ggOut <- ggplot(diffOrganisms, aes(x = "", y = nEpitope)) +
  geom_boxplot() +
  geom_jitter() +
  theme_minimal() +
  xlab("") +
  ylab("Number of significant anti-epitopes") 
ggOut
ggsave(filename = "//vf-lucid-r-o.lumcnet.prod.intern/lucid-r-o$/Projects/Senegal_phipseq/plots/n_sigAntiEpitopes.pdf", plot = ggOut, 
       width =6, height = 8)

median(diffOrganisms$nEpitope)

#Perform PCA

#Keep RT149
PCA_phip <- phipSeq_p

#Remove RT149, as they skew the data and epitopes which everyone has
# PCA_phip <- phipSeq_p[-which(row.names(phipSeq_p) == "RT149"),]
# PCA_phip <- PCA_phip[,names(colSums(PCA_phip)[colSums(PCA_phip) < nrow(PCA_phip)])]

PCA_phip <- {set.seed(50) ;prcomp(PCA_phip, scale = TRUE, center = TRUE)}

fviz_pca_ind(PCA_phip)

embedPC <- as.data.frame(PCA_phip$x)
embedPC$Residence <- metaData[row.names(embedPC), "Residence"]
embedPC$sampleID <- row.names(embedPC)
embedPC$Sex <- metaData[row.names(embedPC), "Sex"]
embedPC$Age <- metaData[row.names(embedPC), "Age"]
embedPC$job_amalgam <- metaData_complex[row.names(embedPC), "job_amalgam"]
embedPC$Residence_complex <- metaData_complex[row.names(embedPC), "Cur_Res_G"]

#Create age categories
embedPC$Age_group[embedPC$Age < 26] <- "18-25"
embedPC$Age_group[embedPC$Age >= 26 & embedPC$Age < 36] <- "26-35"
embedPC$Age_group[embedPC$Age >= 36] <- "36+"

#What PC correlates best with residence and sex
embedPC$Residence_numeric <- embedPC$Residence
embedPC$Residence_numeric[str_detect(embedPC$Residence_numeric, "Rural")] <- "0"
embedPC$Residence_numeric[str_detect(embedPC$Residence_numeric, "Urban Senegalese")] <- "1"
embedPC$Residence_numeric[str_detect(embedPC$Residence_numeric, "Dutch")] <- "2"
embedPC$Residence_numeric <- as.numeric(embedPC$Residence_numeric)

embedPC$Sex_numeric <- embedPC$Sex
embedPC$Sex_numeric[str_detect(embedPC$Sex_numeric, "M")] <- "0"
embedPC$Sex_numeric[str_detect(embedPC$Sex_numeric, "F")] <- "1"
embedPC$Sex_numeric <- as.numeric(embedPC$Sex_numeric)

#PC2 captures urbanization
embedPC_cor <- cor(as.matrix(embedPC[, -which(colnames(embedPC) %in% c("Sex", "Residence", "sampleID", 
                                                                       "Age_group", "job_amalgam", "Residence_complex"))]))

#Get loadings of PC2 (top 10 in each direction)
PC2_loadings <- PCA_phip$rotation[, "PC2"]

PC2_loadings_pos <- PC2_loadings[PC2_loadings > 0]
PC2_loadings_pos <- PC2_loadings_pos[order(PC2_loadings_pos, decreasing = T)]
PC2_loadings_pos <- PC2_loadings_pos[1:10]

PC2_loadings_neg <- PC2_loadings[PC2_loadings < 0]
PC2_loadings_neg <- PC2_loadings_neg[order(PC2_loadings_neg, decreasing = F)]
PC2_loadings_neg <- PC2_loadings_neg[1:10]

PC2_loadings_df <- data.frame("epitope" = c(names(PC2_loadings_pos), names(PC2_loadings_neg)),
                              "PC2_loading" = as.numeric( c(PC2_loadings_pos, PC2_loadings_neg) ),
                              row.names = c(names(PC2_loadings_pos), names(PC2_loadings_neg)) )

PC2_loadings_df$Organism <- firstup(peptideInfo[row.names(PC2_loadings_df), "Pathogen"])
PC2_loadings_df$Protein <- firstup(peptideInfo[row.names(PC2_loadings_df), "Protein"])
PC2_loadings_df$Organism_protein <- paste0(PC2_loadings_df$epitope, ": " ,PC2_loadings_df$Protein, " (",PC2_loadings_df$Organism ,")")

PC2_loadings_df$direction <- "Urban Dutch"
PC2_loadings_df$direction[which(PC2_loadings_df$PC2_loading < 0)] <- "Rural Senegalese"

PC2_loadings_df$epitope <- factor(PC2_loadings_df$epitope, levels = )

PC2_loadings_df <- PC2_loadings_df[order(PC2_loadings_df$PC2_loading),]


ggOut <- ggplot(PC2_loadings_df, aes(x  = PC2_loading, y = forcats::fct_inorder(Organism_protein), fill = direction)) +
        geom_bar(stat = "identity") +
        scale_fill_manual(values = area_palette) +
        theme_minimal() +
        ylab("") +
        xlab("PC2 loadings")
ggsave(filename = "//vf-lucid-r-o.lumcnet.prod.intern/lucid-r-o$/Projects/Senegal_phipseq/plots/PC2_loadings.pdf", plot = ggOut, 
       width =15, height = 8)

#Boxplot of PC2 values split by Residence
PC2_lm <- lm(PC2 ~ Residence, data = embedPC)

PC2_emm <- emmeans(PC2_lm, "Residence")
trendOutput <- as.data.frame(contrast(PC2_emm, "poly"))

pairwiseOutput <- as.data.frame(pairs(PC2_emm))

group1 <- unlist(str_split(pairwiseOutput$contrast, " - "))[c(TRUE, FALSE)]
group2 <- unlist(str_split(pairwiseOutput$contrast, " - "))[c(FALSE, TRUE)]

ggpubr_formatted <- data.frame(".y."= "test", "group1" = group1, "group2" = group2,
                               p = pairwiseOutput$p.value)

dataMax <- max(embedPC$PC2)

ggpubr_formatted <- ggpubr_formatted %>%
  mutate(y.position = c( dataMax + (dataMax * 0.1) , dataMax + (dataMax * 0.2), dataMax + (dataMax * 0.3)))

ggpubr_formatted$p.symbol <- gtools::stars.pval(ggpubr_formatted$p)
ggpubr_formatted$p.symbol[ggpubr_formatted$p >= 0.05] <- "NS"

trendP <- subset(trendOutput, contrast == "linear")$p.value
if(trendP < 0.05){trendP <- "< 0.05"}else{trendP <- "\u2265 0.05"}

ggOut <- ggboxplot(embedPC, y = "PC2", x = "Residence", fill = "Residence", outlier.shape = NA) +
  stat_pvalue_manual(ggpubr_formatted, label = "p.symbol") +
  scale_fill_manual(values = area_palette) +
  ggtitle(label = paste0("Trend p-value ", trendP)) +
  ylab("PC2") +
  xlab("") +
  geom_jitter() +
  grids(linetype = "solid")
ggOut
ggsave(filename = "//vf-lucid-r-o.lumcnet.prod.intern/lucid-r-o$/Projects/Senegal_phipseq/plots/PC2_loadings_boxplot.pdf", plot = ggOut, 
       width =6, height = 8)

#PC2 and PC6 split sex

#Get variance explained
PC1_varExplained <- (PCA_phip$sdev[1]^2 / sum(PCA_phip$sdev^2)) * 100
PC2_varExplained <- (PCA_phip$sdev[2]^2 / sum(PCA_phip$sdev^2)) * 100
PC3_varExplained <- (PCA_phip$sdev[3]^2 / sum(PCA_phip$sdev^2)) * 100

embedPC$Residence <- factor(embedPC$Residence, levels = c("Rural Senegalese", "Urban Senegalese", "Urban Dutch"))

#Correlation of Age and PC2
cor.test(embedPC$PC2, embedPC$Age, method = "pearson")
cor.test(embedPC$PC1, embedPC$Age, method = "pearson")

#Save PCA data
write.csv(embedPC, "//vf-lucid-r-o.lumcnet.prod.intern/lucid-r-o$/Projects/Senegal_phipseq/data/processed/PC_coords.csv")

#PCA plot coloured by Area
gg <- merge(embedPC,aggregate(cbind(mean.PC1=PC1,mean.PC2=PC2)~Residence,embedPC,mean),by="Residence")
ggOut <- ggplot(gg, aes(PC1,PC2))+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0) +
  geom_point(size=5, aes(color=Residence))+
  geom_segment(aes(x=mean.PC1, y=mean.PC2, xend=PC1, yend=PC2, color = Residence)) +
  geom_point(aes(x=mean.PC1,y=mean.PC2, fill=Residence), size=7, shape = 21) +
  theme_minimal()+
  scale_fill_manual(values = area_palette) +
  scale_color_manual(values = area_palette)+
  xlab(paste0("PC1 (", round(PC1_varExplained, digits = 2), "%)"))+
  ylab(paste0("PC2 (", round(PC2_varExplained, digits = 2), "%)"))
ggOut
ggsave(filename = "plots/figure_components/PC1PC2Plot_filteredData_residenceColour.pdf", plot = ggOut, 
       width = 9, height = 8)
ggsave(filename = "//vf-lucid-r-o.lumcnet.prod.intern/lucid-r-o$/Projects/Senegal_phipseq/plots/PC1PC2Plot_filteredData_residenceColour.pdf", plot = ggOut, 
       width = 9, height = 8)

#PCA plot coloured by Sex
gg <- merge(embedPC,aggregate(cbind(mean.PC1=PC1,mean.PC2=PC2)~Sex,embedPC,mean),by="Sex")
ggOut <- ggplot(gg, aes(PC1,PC2))+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0) +
  geom_point(size=5, aes(color=Sex))+
  geom_segment(aes(x=mean.PC1, y=mean.PC2, xend=PC1, yend=PC2, color = Sex)) +
  geom_point(aes(x=mean.PC1,y=mean.PC2, fill=Sex), size=7, shape = 21) +
  theme_minimal() +
  xlab(paste0("PC1 (", round(PC1_varExplained, digits = 2), "%)"))+
  ylab(paste0("PC2 (", round(PC2_varExplained, digits = 2), "%)"))
ggOut
ggsave(filename = "plots/figure_components/PC1PC2Plot_filteredData_sexColour.pdf", plot = ggOut, 
       width = 9, height = 8)
ggsave(filename = "//vf-lucid-r-o.lumcnet.prod.intern/lucid-r-o$/Projects/Senegal_phipseq/plots/PC1PC2Plot_filteredData_sexColour.pdf", plot = ggOut, 
       width = 9, height = 8)

#PCA plot coloured by Age
ggOut <- ggplot(embedPC, aes(PC1,PC2, color = Age))+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0) +
  geom_point(size=5) +
  theme_minimal() +
  xlab(paste0("PC1 (", round(PC1_varExplained, digits = 2), "%)"))+
  ylab(paste0("PC2 (", round(PC2_varExplained, digits = 2), "%)")) +
  scale_color_viridis_b()
ggOut
ggsave(filename = "plots/figure_components/PC1PC2Plot_filteredData_ageColour.pdf", plot = ggOut, 
       width = 9, height = 8)
ggsave(filename = "//vf-lucid-r-o.lumcnet.prod.intern/lucid-r-o$/Projects/Senegal_phipseq/plots/PC1PC2Plot_filteredData_ageColour.pdf", plot = ggOut, 
       width = 9, height = 8)

#PCA plot coloured by age group
#18-25, 26-35, 36+

ageGroup_palette <- c("18-25" = "#bea0cc", "26-35" = "#ad68af", "36+" = "#562888")

gg <- merge(embedPC,aggregate(cbind(mean.PC1=PC1,mean.PC2=PC2)~job_amalgam,embedPC,mean),by="job_amalgam")
ggOut <- ggplot(gg, aes(PC1,PC2))+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0) +
  geom_point(size=3, aes(color=job_amalgam)) +
  geom_segment(aes(x=mean.PC1, y=mean.PC2, xend=PC1, yend=PC2, color = job_amalgam)) +
  geom_point(aes(x=mean.PC1,y=mean.PC2, fill=job_amalgam), size=5, shape = 21) +
  theme_minimal() + 
  xlab(paste0("PC1 (", round(PC1_varExplained, digits = 2), "%)"))+
  ylab(paste0("PC2 (", round(PC2_varExplained, digits = 2), "%)"))
ggOut
ggsave(filename = "plots/figure_components/PC1PC2Plot_filteredData_ageGroupColour.pdf", plot = ggOut, 
       width = 9, height = 8)
ggsave(filename = "//vf-lucid-r-o.lumcnet.prod.intern/lucid-r-o$/Projects/Senegal_phipseq/plots/PC1PC2Plot_filteredData_ageGroupColour.pdf", plot = ggOut, 
       width = 9, height = 8)


#PCA plot coloured by working environment
gg <- merge(embedPC,aggregate(cbind(mean.PC1=PC1,mean.PC2=PC2)~Residence,embedPC,mean),by="Residence")
ggOut <- ggplot(gg, aes(PC1,PC2))+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0) +
  geom_point(size=3, aes(color=Residence))+
  geom_segment(aes(x=mean.PC1, y=mean.PC2, xend=PC1, yend=PC2, color = Residence)) +
  geom_point(aes(x=mean.PC1,y=mean.PC2, fill=Residence), size=5, shape = 21) +
  theme_minimal() +
  xlab(paste0("PC1 (", round(PC1_varExplained, digits = 2), "%)"))+
  ylab(paste0("PC2 (", round(PC2_varExplained, digits = 2), "%)")) +
  facet_wrap(~job_amalgam)+
  scale_fill_manual(values = area_palette) +
  scale_color_manual(values = area_palette)
ggOut
ggsave(filename = "plots/figure_components/PC1PC2Plot_filteredData_residenceColour_splitWorkingEnv.pdf", plot = ggOut, 
       width = 9, height = 8)
ggsave(filename = "//vf-lucid-r-o.lumcnet.prod.intern/lucid-r-o$/Projects/Senegal_phipseq/plots/PC1PC2Plot_filteredData_residenceColour_splitWorkingEnv.pdf", plot = ggOut, 
       width = 9, height = 8)

#Does the total epitope count correlate with PC1? - Yes!!
totalEpitopeCount <- rowSums(phipSeq_p)

cor(totalEpitopeCount, embedPC[names(totalEpitopeCount),"PC1"])


#Perform pairwise permanova on the different groups
pairwise_permanova <- function(sp_matrix, group_var, dist = "bray", adj = "fdr", perm = 10000) {
  
  require(vegan)
  
  ## list contrasts
  group_var <- as.character(group_var)
  groups <- as.data.frame(t(combn(unique(group_var), m = 2)))
  
  contrasts <- data.frame(
    group1 = groups$V1, group2 = groups$V2,
    R2 = NA, F_value = NA, df1 = NA, df2 = NA, p_value = NA
  )
  
  for (i in seq(nrow(contrasts))) {
    sp_subset <- group_var == contrasts$group1[i] | group_var == contrasts$group2[i] 
    contrast_matrix <- sp_matrix[sp_subset,]
    
    ## fit contrast using adonis
    fit <- vegan::adonis2(
      contrast_matrix ~ group_var[sp_subset],
      method = dist, 
      perm = perm
    )
    
    contrasts$R2[i] <- round(fit$R2[1], digits = 3)
    contrasts$F_value[i] <- round(fit[["F"]][1], digits = 3)
    contrasts$df1[i] <- fit$Df[1]
    contrasts$df2[i] <- fit$Df[2]
    contrasts$p_value[i] <- fit$`Pr(>F)`[1]
  }
  
  ## adjust p-values for multiple comparisons
  contrasts$p_value <- round(p.adjust(contrasts$p_value, method = adj), digits = 3)
  
  return(list(
    contrasts = contrasts, 
    "p-value adjustment" = adj, 
    permutations = perm
  ))
}

#Residence
pairwise_permanova(sp_matrix = embedPC[, c("PC1", "PC2")], group_var = embedPC$Residence)
pairwise_permanova(sp_matrix = embedPC[, c("PC1", "PC2")], group_var = embedPC$Sex)
pairwise_permanova(sp_matrix = embedPC[, c("PC1", "PC2")], group_var = embedPC$Age_group)


#Perform Cochran Armitage test to test for trend across Urbanization. Adjust p-value with bonferronni

phipSeq_p_df$Residence <- metaData[row.names(phipSeq_p_df),"Residence"]

convertBinary <- reshape2::melt(phipSeq_p_df)
convertBinary$variable <- as.character(convertBinary$variable)
convertBinary$Residence <- factor(convertBinary$Residence, levels =  c("Rural Senegalese" , "Urban Senegalese" , "Urban Dutch" ))

testDF <- data.frame()

for(i in unique(convertBinary$variable)){
  current <- subset(convertBinary, variable == i)
  
  testResult <- CochranArmitageTest( table(current$Residence, current$value) )
  
  testDF <- rbind(testDF, data.frame("ID" = i,
                                     "p.value" = testResult$p.value,
                                     "Z" = testResult$statistic))
  
}

#Add meta data
testDF$Pathogen <- peptideInfo[testDF$ID, "Pathogen"]
testDF$PathogenComplex <- peptideInfo[testDF$ID, "PathogenComplex"]
testDF$Protein <- peptideInfo[testDF$ID, "Protein"]
testDF$protein_pathogen <- peptideInfo[testDF$ID, "protein_pathogen"]
testDF$aa_seq <- peptideInfo[testDF$ID, "aa_seq"]

testDF$p.adjust <- p.adjust(testDF$p.value, method = "BH")

testSigDF <- subset(testDF, p.adjust < 0.05)

testSigDF$direction <- "increaseUrbanization"
testSigDF$direction[testSigDF$Z > 0] <- "decreaseUrbanization"

ggplot() + geom_point(aes(x= testDF$Z,y= testDF$p.value))


#Save significance results
write.csv(testDF, "//vf-lucid-r-o.lumcnet.prod.intern/lucid-r-o$/Projects/Senegal_phipseq/data/processed/cochranArmitage_urbanization.csv")

#Save amino acid sequence of sig epitopes
testDF$aa_seq <- peptideInfo[testDF$ID ,"aa_seq"]
write.csv(testDF, "//vf-lucid-r-o.lumcnet.prod.intern/lucid-r-o$/Projects/Senegal_phipseq/data/processed/cochranArmitage_urbanization_aaSeq.csv")

### OLD SECTION start ###
# 
# #Investigate organisms
# decreaseWithUrbanization <- table(subset(testSigDF, Z > 0)$Pathogen)
# increaseWithUrbanization <- table(subset(testSigDF, Z < 0)$Pathogen)
# 
# #Normalise for how many times an epitope of that organism is present in the dataset
# pathogenTotals <- table(peptideInfo$Pathogen)
# 
# decreaseWithUrbanization <- decreaseWithUrbanization / pathogenTotals[names(decreaseWithUrbanization)]
# increaseWithUrbanization <- increaseWithUrbanization / pathogenTotals[names(increaseWithUrbanization)]
# 
# wordcloud2(decreaseWithUrbanization, size = 0.4)
# wordcloud2(increaseWithUrbanization, size = 0.8)
# 
# #Plot the enrichment of the organisms 
# plotDF <- data.frame("enrichment" = c(as.numeric(decreaseWithUrbanization), as.numeric(increaseWithUrbanization) * -1 ),
#                      "organism" = c(as.character(names(decreaseWithUrbanization)), as.character(names(increaseWithUrbanization)) ),
#                      "direction" = c( rep("Rural Senegalese", length(decreaseWithUrbanization)), rep("Urban Dutch", length(increaseWithUrbanization)) ))
# 
# #Order by enrichment score
# plotDF <- plotDF[order(plotDF$enrichment, decreasing = F), ]
# 
# plotDF$organism <- factor(plotDF$organism, unique(plotDF$organism))

### OLD SECTION end ###

#Average trend z value across direction and organism

plotDF <- testSigDF %>%
                    group_by(Pathogen, direction) %>%
                    summarise(medianZ = median(Z))

plotDF$direction[str_detect(plotDF$direction, "increase")] <- "Urban Dutch"
plotDF$direction[str_detect(plotDF$direction, "decrease")] <- "Rural Senegalese"

plotDF <- plotDF[order(plotDF$medianZ, decreasing = F), ]
plotDF$Pathogen <- factor(plotDF$Pathogen, unique(plotDF$Pathogen))

#Get the most common protein from the different pathogens. If there are equal numbers of protein, choose one randomly
mostCommon <- testSigDF %>%
                group_by(Pathogen) %>%
                summarise(commonEpitope = paste0( names(which.max(table(Protein))), " (n=", table(Protein)[which.max(table(Protein))]  ,")" ) )
mostCommon <- as.data.frame(mostCommon)
row.names(mostCommon) <- mostCommon$Pathogen

#Several proteins are uncharacterized, add the organism to them to distinguish
mostCommon$commonEpitope[str_detect(mostCommon$commonEpitope, "uncharacterized")] <- paste0(mostCommon$Pathogen[str_detect(mostCommon$commonEpitope, "uncharacterized")], " ", mostCommon$commonEpitope[str_detect(mostCommon$commonEpitope, "uncharacterized")])

#Add to the plotDF
plotDF$commonEpitope <- mostCommon[as.character(plotDF$Pathogen),"commonEpitope"]
plotDF$commonEpitope <- firstup(plotDF$commonEpitope)
plotDF$commonEpitope <- factor(plotDF$commonEpitope, unique(plotDF$commonEpitope))

plotDF$Pathogen_commonEpitope <- paste0(plotDF$Pathogen, ": ", plotDF$commonEpitope)
plotDF$Pathogen_commonEpitope<- firstup(plotDF$Pathogen_commonEpitope)

plotDF <- plotDF[order(plotDF$medianZ, decreasing = T), ]
plotDF$Pathogen_commonEpitope <- factor(plotDF$Pathogen_commonEpitope, unique(plotDF$Pathogen_commonEpitope))

plotDF$medianZ <- plotDF$medianZ * -1

ggOut <- ggplot(plotDF, aes(x = medianZ, y = Pathogen_commonEpitope, fill = direction)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  xlab("Median Trend Z value") +
  scale_fill_manual(values = area_palette)
ggOut 
ggsave(filename = paste0("plots/figure_components/commonProtein_cochranArmitage_medianTrendZ_barplot.pdf"), plot = ggOut)
ggsave(filename = "//vf-lucid-r-o.lumcnet.prod.intern/lucid-r-o$/Projects/Senegal_phipseq/plots/commonProtein_cochranArmitage_medianTrendZ_barplot.pdf", plot = ggOut, 
       width = 16, height = 8)

#NLR relationship with shigella epitopes
normTotalEpitope_DF$NLR <- metaData_complex[row.names(normTotalEpitope_DF),"NLR"]

ggplot(normTotalEpitope_DF, aes(x = NLR, y  = shigellaflexneriTotalNorm, color = Residence)) +
        geom_point()

ggplot(normTotalEpitope_DF, aes(x = NLR, y  = shigellasonneiTotalNorm, color = Residence)) +
  geom_point()

#IL-8 circulating levels and shigella epitopes
luminexData <- readRDS("urban-rural-study-manuscript/data/data_luminex.rds")
row.names(luminexData) <- luminexData$donor

normTotalEpitope_DF$shigellaTotalNorm <- rowSums(normTotalEpitope_DF[, str_detect(colnames(normTotalEpitope_DF), "shigelladys|shigellaflex|shigellason|shigellaboy|shigellasen")])

normTotalEpitope_DF$IL8 <- luminexData[row.names(normTotalEpitope_DF),"il8_cxcl8"]

ggplot(normTotalEpitope_DF, aes(x = IL8, y  = shigellaflexneriTotalNorm, color = Residence)) +
  geom_point() +
  stat_cor()
  

ggplot(normTotalEpitope_DF, aes(x = shigellaflexneriTotalNorm, y  = shigellaphageTotalNorm, color = Residence)) +
  geom_point()

#Look at the significant results
phipSeq_p_freq_residence_raw_df$Residence <- row.names(phipSeq_p_freq_residence_raw_df)
plotDF <- melt(phipSeq_p_freq_residence_raw_df)

plotDF$Residence <- factor(plotDF$Residence, c("Rural Senegalese", "Urban Senegalese", "Urban Dutch") )

#What pathogens have different epitopes present in different groups

#Increasing with urbanization
ggOut <- ggplot(subset(plotDF, variable %in% subset(testSigDF, Z < 0)$ID), aes(x = Residence, y = value, fill = Residence)) +
  geom_bar(stat = "identity") +
  facet_wrap(~variable) +
  scale_fill_manual(values = area_palette) +
  ylim(0, 100)+ 
  ylab("Ab against epitope present (%)")+
  theme_minimal()+
  rotate_x_text()
ggOut
ggsave(filename = paste0("plots/figure_components/sup_figures/cochranArmitage_increaseWithUrbanization_epitopePercentPresent.pdf"), plot = ggOut, width = 20, height = 20)
ggsave(filename = "//vf-lucid-r-o.lumcnet.prod.intern/lucid-r-o$/Projects/Senegal_phipseq/plots/cochranArmitage_increaseWithUrbanization_epitopePercentPresent.pdf", plot = ggOut, 
       width = 14, height = 12)

#Decreasing with urbanization
ggOut <- ggplot(subset(plotDF, variable %in% subset(testSigDF, Z > 0)$ID), aes(x = Residence, y = value, fill = Residence)) +
  geom_bar(stat = "identity") +
  facet_wrap(~variable) +
  scale_fill_manual(values = area_palette) +
  ylim(0, 100) + 
  ylab("Ab against epitope present (%)")+
  theme_minimal() +
  rotate_x_text()
ggOut
ggsave(filename = paste0("plots/figure_components/sup_figures/cochranArmitage_decreaseWithUrbanization_epitopePercentPresent.pdf"), plot = ggOut, width = 20, height = 20)
ggsave(filename = "//vf-lucid-r-o.lumcnet.prod.intern/lucid-r-o$/Projects/Senegal_phipseq/plots/cochranArmitage_decreaseWithUrbanization_epitopePercentPresent.pdf", plot = ggOut, 
       width = 14, height = 12)

#Look at staph aureus exported protein
aureusEMP<- subset(peptideInfo,len_seq == "315.0" & Protein == "exported protein" )

ggOut <- ggplot(subset(plotDF, variable %in% row.names(aureusEMP)), aes(x = Residence, y = value, fill = Residence)) +
  geom_bar(stat = "identity") +
  facet_wrap(~variable) +
  scale_fill_manual(values = area_palette) +
  ylim(0, 100) +
  theme_minimal() +
  rotate_x_text()
ggOut

#Plot the top three enriched pathogens as the normalised frequencies
plotDF <- normTotalEpitope_DF[,-which(colnames(normTotalEpitope_DF) == "epitopeTotals")]

ggplot(plotDF, aes(x = staphylococcusaureusTotalNorm, y = Residence, fill = Residence)) +
        geom_boxplot() +
        geom_jitter() +
        theme_minimal()  +
        scale_fill_manual(values = area_palette) +
        xlab("Proportion of Ab-OME against organism epitopes") +
        ylab("")

ggplot(plotDF, aes(x = streptococcuspyogenesTotalNorm, y = Residence, fill = Residence)) +
  geom_boxplot() +
  geom_jitter() +
  theme_minimal()  +
  scale_fill_manual(values = area_palette) +
  xlab("Proportion of Ab-OME against organism epitopes") +
  ylab("")

ggplot(plotDF, aes(x = homosapiensTotalNorm, y = Residence, fill = Residence)) +
  geom_boxplot() +
  geom_jitter() +
  theme_minimal()  +
  scale_fill_manual(values = area_palette) +
  xlab("Proportion of Ab-OME against organism epitopes") +
  ylab("")


top4 <- c("humanherpesvirus1TotalNorm", "shigellasonneiTotalNorm", "shigellaflexneriTotalNorm", "humangammaherpesvirus8TotalNorm",
          "haemophilusinfluenzaeTotalNorm", "epsteinbarrvirusTotalNorm", "escherichiacoliTotalNorm", "streptococcusdysgalactiaeTotalNorm")

#Calculate significance
featuresTest <- normTotalEpitope_DF[, c(top4, "Residence")]
unique(featuresTest$Residence)

plotList <- list()

for(i in top4){
  plotName <- str_replace_all(i, "TotalNorm", "")
  
  featuresTest_lm <- lm(eval(parse(text = i)) ~ Residence, data = featuresTest)
  
  featureTest_emm <- emmeans(featuresTest_lm, "Residence")
  trendOutput <- as.data.frame(contrast(featureTest_emm, "poly"))
  
  pairwiseOutput <- as.data.frame(pairs(featureTest_emm))
  
  group1 <- unlist(str_split(pairwiseOutput$contrast, " - "))[c(TRUE, FALSE)]
  group2 <- unlist(str_split(pairwiseOutput$contrast, " - "))[c(FALSE, TRUE)]
  
  ggpubr_formatted <- data.frame(".y."= "test", "group1" = group1, "group2" = group2,
                                 p = pairwiseOutput$p.value)
  
  dataMax <- max(featuresTest[, i])
  
  ggpubr_formatted <- ggpubr_formatted %>%
    mutate(y.position = c( dataMax + (dataMax * 0.1) , dataMax + (dataMax * 0.2), dataMax + (dataMax * 0.3)))
  
  ggpubr_formatted$p.symbol <- gtools::stars.pval(ggpubr_formatted$p)
  ggpubr_formatted$p.symbol[ggpubr_formatted$p >= 0.05] <- "NS"
  
  trendP <- subset(trendOutput, contrast == "linear")$p.value
  if(trendP < 0.05){trendP <- "< 0.05"}else{trendP <- "\u2265 0.05"}
  
  ggOut <- ggboxplot(featuresTest, y = i, x = "Residence", fill = "Residence", outlier.shape = NA) +
    stat_pvalue_manual(ggpubr_formatted, label = "p.symbol") +
    scale_fill_manual(values = area_palette) +
    ggtitle(label = paste0("Trend p-value ", trendP)) +
    ylab(plotName) +
    xlab("") +
    geom_jitter() +
    grids(linetype = "solid") +
    rotate() + 
    guides(fill="none")
  
  ggsave(filename = paste0("//vf-lucid-r-o.lumcnet.prod.intern/lucid-r-o$/Projects/Senegal_phipseq/plots/normPathogenPlots/",i,"_enrichment_urbanization_normPathogen.pdf"), plot = ggOut, 
         width = 10, height = 4)

  plotList[[i]] <- ggOut
}


ggsave(filename = "//vf-lucid-r-o.lumcnet.prod.intern/lucid-r-o$/Projects/Senegal_phipseq/plots/organism_enrichment_urbanization_normPathogen.pdf", 
       plot = ggarrange(plotlist = plotList, ncol = 6, nrow = 6), 
       width = 28, height = 12)


#Test all organisms in the significant epitopes list
sigOrganisms <- c(paste0(unique(str_replace_all(testSigDF$Pathogen, "-| |\\/", "")), "TotalNorm"))

featuresTest <- normTotalEpitope_DF[, c(sigOrganisms, "Residence")]
unique(featuresTest$Residence)

plotList <- list()
pvalDF <- data.frame()
ggpubrFormatList <- list()

for(i in sigOrganisms){
  
  featuresTest_lm <- lm(eval(parse(text = i)) ~ Residence, data = featuresTest)
  
  featureTest_emm <- emmeans(featuresTest_lm, "Residence")
  trendOutput <- as.data.frame(contrast(featureTest_emm, "poly"))
  
  pairwiseOutput <- as.data.frame(pairs(featureTest_emm))
  
  group1 <- unlist(str_split(pairwiseOutput$contrast, " - "))[c(TRUE, FALSE)]
  group2 <- unlist(str_split(pairwiseOutput$contrast, " - "))[c(FALSE, TRUE)]
  
  ggpubr_formatted <- data.frame(".y."= "test", "group1" = group1, "group2" = group2,
                                 p = pairwiseOutput$p.value)
  
  dataMax <- max(featuresTest[, i])
  
  ggpubr_formatted <- ggpubr_formatted %>%
    mutate(y.position = c( dataMax + (dataMax * 0.1) , dataMax + (dataMax * 0.2), dataMax + (dataMax * 0.3)))
  
  ggpubr_formatted$p.symbol <- gtools::stars.pval(ggpubr_formatted$p)
  ggpubr_formatted$p.symbol[ggpubr_formatted$p >= 0.05] <- "NS"
  
  trendP <- subset(trendOutput, contrast == "linear")$p.value
  #if(trendP < 0.05){trendP <- "< 0.05"}else{trendP <- "\u2265 0.05"}
  
  pvalDF <- rbind(pvalDF, data.frame(organism = i,
                                     pval = trendP))
  
  ggpubrFormatList[[i]] <- ggpubr_formatted
  
}

pvalDF$p_adj <- p.adjust(pvalDF$pval, method = "BH")

#Make plots
for(i in sigOrganisms){

  currentTrendP <- subset(pvalDF, organism == i)$p_adj
  current_ggpubr_formatted <- ggpubrFormatList[[i]]
  plotName <- str_replace_all(i, "TotalNorm", "")
  
  ggOut <- ggboxplot(featuresTest, y = i, x = "Residence", fill = "Residence", outlier.shape = NA) +
    stat_pvalue_manual(current_ggpubr_formatted, label = "p.symbol") +
    scale_fill_manual(values = area_palette) +
    ggtitle(label = paste0("Trend p-value ", currentTrendP)) +
    ylab(plotName) +
    xlab("") +
    geom_jitter() +
    grids(linetype = "solid") +
    rotate()+ 
    guides(fill="none")

  ggsave(filename = paste0("//vf-lucid-r-o.lumcnet.prod.intern/lucid-r-o$/Projects/Senegal_phipseq/plots/normPathogenPlots/",i,"_enrichment_urbanization_normPathogen.pdf"), plot = ggOut, 
       width = 5.5, height = 4)
  
  plotList[[i]] <- ggOut
  

}



ggsave(filename = "//vf-lucid-r-o.lumcnet.prod.intern/lucid-r-o$/Projects/Senegal_phipseq/plots/organism_enrichment_urbanization_normPathogen.pdf", 
       plot = ggarrange(plotlist = plotList, ncol = 6, nrow = 6), 
       width = 28, height = 22)


#Get the individuals who have antibody against gammaherpes virus 8
individualsGH8 <- row.names(subset(normTotalEpitope_DF, humangammaherpesvirus8TotalNorm > 0))
individualsGH8_meta <- metaData_complex[individualsGH8,]

individualsCMV <- row.names(subset(normTotalEpitope_DF, humancytomegalovirusTotalNorm > 0))
individualsCMV_meta <- metaData_complex[individualsCMV,]

embedPC$GH8Pos <- "negative"
embedPC[individualsGH8,"GH8Pos"] <- "positive"

embedPC$CMVPos <- "negative"
embedPC[individualsCMV,"CMVPos"] <- "positive"


ggplot(embedPC, aes(x = PC1, y = PC2, color = GH8Pos)) +
          geom_point() +
          facet_wrap(~Residence)

ggplot(embedPC, aes(x = PC1, y = PC2, color = CMVPos)) +
  geom_point()

#Linear trend for all pathogens
testTotalNorm_df <- data.frame()

colnames(normTotalEpitope_DF) <- str_replace_all(colnames(normTotalEpitope_DF), "_|:|\\(|\\)|\\.|\\=|\\/|\\||\\,|\\&", "")

for(i in colnames(normTotalEpitope_DF)[str_detect(colnames(normTotalEpitope_DF), "TotalNorm")]){
  plotName <- str_replace_all(i, "TotalNorm", "")
  
  featuresTest_lm <- lm(eval(parse(text = i)) ~ Residence, data = normTotalEpitope_DF)
  
  featureTest_emm <- emmeans(featuresTest_lm, "Residence")
  trendOutput <- as.data.frame(contrast(featureTest_emm, "poly"))
  
  pairwiseOutput <- as.data.frame(pairs(featureTest_emm))
  
  testCurrent_df <- data.frame("pathogen" = i, 
                               "trendLinear" = subset(trendOutput, contrast == "linear")$p.value,
                               "trendQuadratic" = subset(trendOutput, contrast == "quadratic")$p.value,
                               "estimateLinear" = subset(trendOutput, contrast == "quadratic")$estimate)
  testTotalNorm_df <- rbind(testTotalNorm_df, testCurrent_df)
  
}

testTotalNorm_df$trendLinear <- p.adjust(testTotalNorm_df$trendLinear, method = "BH")
testTotalNorm_df$trendQuadratic <- p.adjust(testTotalNorm_df$trendQuadratic, method = "BH")

#Heatmap of presence/absence of significant epitopes
phipSeq_p_raw_df_sig <- phipSeq_p_raw_df[, testSigDF$ID]
phipSeq_p_raw_df_sig$Residence <- metaData[row.names(phipSeq_p_raw_df_sig), "Residence"]
phipSeq_p_raw_df_sig$Donor <- row.names(phipSeq_p_raw_df_sig)

phipSeq_p_raw_df_sig_melt <- reshape2::melt(phipSeq_p_raw_df_sig)

phipSeq_p_raw_df_sig_melt$value <- as.character(phipSeq_p_raw_df_sig_melt$value)

ggplot(phipSeq_p_raw_df_sig_melt, aes(x = Donor, y = variable, fill = value)) +
                      geom_tile()

heatmapMtx <- t(as.matrix(phipSeq_p_raw_df_sig[, -which(colnames(phipSeq_p_raw_df_sig) %in% c("Residence", "Donor"))]))

heatmapMtx <- heatmapMtx[, rev(colnames(heatmapMtx))]

library(circlize)
col_fun <- colorRamp2(c(0, 1), c("black", "red"))

column_ha <- ComplexHeatmap::HeatmapAnnotation(Residence = metaData[colnames(heatmapMtx), "Residence"],
                               col = list(Residence = c("Rural Senegalese" = "#e59f01", "Urban Senegalese" = "#54b3e8", "Urban Dutch" = "#009e72")))

#Get new row names: [proteinName] ([organismName])
row.names(heatmapMtx) <- paste0(row.names(peptideInfo[row.names(heatmapMtx),]),": ",firstup(peptideInfo[row.names(heatmapMtx),"Protein"]), " (", firstup(peptideInfo[row.names(heatmapMtx),"Pathogen"]) ,")")

dendrogramRes <- {set.seed(42); hclust(dist(heatmapMtx))}
dendrogramCut <- {set.seed(42); cutree(dendrogramRes, k = 4)}
dendrogramCut <- factor(dendrogramCut, levels = c(4,2,1,3))

#Get selected row names
selectedRowNames <- c("twist_58790: Genome polyprotein  (Rhinovirus b)", "agilent_12017: Envelope glycoprotein d  (Human herpesvirus 1)",
                      "agilent_12023: Envelope glycoprotein m  (Human cytomegalovirus)", "agilent_141247: Invasin ipac (Shigella sonnei)",
                      "agilent_230549: Invasin ipac (Shigella flexneri)", "agilent_219967: Grab (Streptococcus pyogenes)",
                      "agilent_241538: Invasin ipaa (Shigella flexneri)", "twist_41031: K8.1 (Human gammaherpesvirus 8)",
                      "twist_54039: Interspersed repeat antigen (Plasmodium falciparum)", "agilent_237918: Adhesin (Haemophilus influenzae)",
                      "agilent_8991: Epstein-barr nuclear antigen 6 (Epstein-barr virus)", "agilent_2709: Syncytin-1 (Homo sapiens)",
                      "agilent_223729: Iga-specific serine endopeptidase autotransporter (Neisseria meningitidis)")
intersect(row.names(heatmapMtx), selectedRowNames)
setdiff(selectedRowNames, row.names(heatmapMtx))

epitopeIndex_DF <- data.frame("index" = seq(1, nrow(heatmapMtx), 1),
                              row.names = row.names(heatmapMtx))

epitopeAnnotate = ComplexHeatmap::rowAnnotation(foo = ComplexHeatmap::anno_mark(at = epitopeIndex_DF[selectedRowNames,"index"], 
                                   labels = selectedRowNames))

pdf(paste0("//vf-lucid-r-o.lumcnet.prod.intern/lucid-r-o$/Projects/Senegal_phipseq/plots/peptideSig_heatmap_selected.pdf"), height = 18, width = 26)
heatmapPlot <- ComplexHeatmap::Heatmap(heatmapMtx, col = col_fun, 
                        cluster_rows = TRUE,  row_split = dendrogramCut, cluster_row_slices = FALSE, right_annotation = epitopeAnnotate,
                        cluster_columns = FALSE, top_annotation = column_ha,
                        height = unit(40, "cm"), width = unit(30, "cm"))
ComplexHeatmap::draw(heatmapPlot)
dev.off()

pdf(paste0("//vf-lucid-r-o.lumcnet.prod.intern/lucid-r-o$/Projects/Senegal_phipseq/plots/peptideSig_heatmap.pdf"), height = 18, width = 26)
heatmapPlot <- ComplexHeatmap::Heatmap(heatmapMtx, col = col_fun, 
                                       cluster_rows = TRUE,  row_split = dendrogramCut, cluster_row_slices = FALSE,
                                       cluster_columns = FALSE, top_annotation = column_ha,
                                       height = unit(40, "cm"), width = unit(30, "cm"))
ComplexHeatmap::draw(heatmapPlot)
dev.off()

#Save epitope info and epitope significance
exportEpitopeStats <- testSigDF

exportEpitopeStats$aa_seq <- NULL
exportEpitopeStats$PathogenComplex <- NULL
exportEpitopeStats$protein_pathogen <- NULL

#Replace Pathogen with organism

colnames(exportEpitopeStats) <- str_replace_all(colnames(exportEpitopeStats), "Pathogen", "Organism")

write.csv(x = exportEpitopeStats, 
          file = "//vf-lucid-r-o.lumcnet.prod.intern/lucid-r-o$/Projects/Senegal_phipseq/data/cochranArmitage_sig_epitopes.csv")

exportEpitopeInfo <- peptideInfo
exportEpitopeInfo$protein_pathogen <- NULL
exportEpitopeInfo$Protein_uniref <- NULL
exportEpitopeInfo$Tax_uniref <- NULL
exportEpitopeInfo$iedb_name <- NULL
exportEpitopeInfo$twistName <- NULL
exportEpitopeInfo$twistProtein <- NULL
exportEpitopeInfo$PathogenComplex <- NULL
exportEpitopeInfo$epitopeGroup <- NULL
exportEpitopeInfo$epitopeUnique <- NULL

#Calculate the number of individuals the epitopes are expressed in
nIndividuals <- colSums(phipSeq_p_raw)

exportEpitopeInfo$nIndividuals_detected <- nIndividuals[row.names(exportEpitopeInfo)]

colnames(exportEpitopeInfo) <- str_replace_all(colnames(exportEpitopeInfo), "Pathogen", "Organism")

write.csv(x = exportEpitopeInfo, 
          file = "//vf-lucid-r-o.lumcnet.prod.intern/lucid-r-o$/Projects/Senegal_phipseq/data/epitope_information.csv")


#### Investigate urban, rural and dutch: cochran-armitage test (epitope presence/absence as input) ####

#Perform Cochran Armitage test to test for trend across Urbanization. Adjust p-value with bonferronni

phipSeq_p_df$Residence <- metaData[row.names(phipSeq_p_df),"Residence"]

convertBinary <- reshape2::melt(phipSeq_p_df)
convertBinary$variable <- as.character(convertBinary$variable)
convertBinary$Residence <- factor(convertBinary$Residence, levels =  c("Urban Senegalese" , "Rural Senegalese" , "Urban Dutch" ))

testDF <- data.frame()

for(i in unique(convertBinary$variable)){
  current <- subset(convertBinary, variable == i)
  
  testResult <- CochranArmitageTest( table(current$Residence, current$value) )
  
  testDF <- rbind(testDF, data.frame("ID" = i,
                                     "p.value" = testResult$p.value,
                                     "Z" = testResult$statistic))
  
}

#Add meta data
testDF$Pathogen <- peptideInfo[testDF$ID, "Pathogen"]
testDF$PathogenComplex <- peptideInfo[testDF$ID, "PathogenComplex"]
testDF$Protein <- peptideInfo[testDF$ID, "Protein"]
testDF$protein_pathogen <- peptideInfo[testDF$ID, "protein_pathogen"]
testDF$aa_seq <- peptideInfo[testDF$ID, "aa_seq"]

testDF$p.adjust <- p.adjust(testDF$p.value, method = "BH")

testSigDF <- subset(testDF, p.adjust < 0.05)

testSigDF$direction <- "increaseUrbanization"
testSigDF$direction[testSigDF$Z > 0] <- "decreaseUrbanization"

ggplot() + geom_point(aes(x= testDF$Z,y= testDF$p.value))

#write.csv(testDF, "//vf-lucid-r-o.lumcnet.prod.intern/lucid-r-o$/Projects/Senegal_phipseq/data/processed/cochranArmitage_urbanization.csv")

#Investigate organisms
decreaseWithUrbanization <- table(subset(testSigDF, Z > 0)$Pathogen)
increaseWithUrbanization <- table(subset(testSigDF, Z < 0)$Pathogen)

#Normalise for how many times an epitope of that organism is present in the dataset
pathogenTotals <- table(peptideInfo$Pathogen)

decreaseWithUrbanization <- decreaseWithUrbanization / pathogenTotals[names(decreaseWithUrbanization)]
increaseWithUrbanization <- increaseWithUrbanization / pathogenTotals[names(increaseWithUrbanization)]

wordcloud2(decreaseWithUrbanization, size = 0.4)
wordcloud2(increaseWithUrbanization, size = 0.8)

#Plot the enrichment of the organisms 
plotDF <- data.frame("enrichment" = c(as.numeric(decreaseWithUrbanization), as.numeric(increaseWithUrbanization) * -1 ),
                     "organism" = c(as.character(names(decreaseWithUrbanization)), as.character(names(increaseWithUrbanization)) ),
                     "direction" = c( rep("Rural Senegalese", length(decreaseWithUrbanization)), rep("Urban Dutch", length(increaseWithUrbanization)) ))

#Order by enrichment score
plotDF <- plotDF[order(plotDF$enrichment, decreasing = F), ]

plotDF$organism <- factor(plotDF$organism, unique(plotDF$organism))

ggOut <- ggplot(plotDF, aes(x = enrichment, y = organism, fill = direction)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  xlab("Epitope enrichment Score") +
  scale_fill_manual(values = area_palette)
ggOut
#ggsave(filename = paste0("plots/figure_components/cochranArmitage_enrichment_barplot.pdf"), plot = ggOut)

#Look at the significant results
phipSeq_p_freq_residence_raw_df$Residence <- row.names(phipSeq_p_freq_residence_raw_df)
plotDF <- melt(phipSeq_p_freq_residence_raw_df)

plotDF$Residence <- factor(plotDF$Residence, c("Rural Senegalese", "Urban Senegalese", "Urban Dutch") )

#What pathogens have different epitopes present in different groups

#Increasing with urbanization
ggOut <- ggplot(subset(plotDF, variable %in% subset(testSigDF, Z < 0)$ID), aes(x = Residence, y = value, fill = Residence)) +
  geom_bar(stat = "identity") +
  facet_wrap(~variable) +
  scale_fill_manual(values = area_palette) +
  ylim(0, 100)+
  theme_minimal()
ggOut
#ggsave(filename = paste0("plots/figure_components/sup_figures/cochranArmitage_increaseWithUrbanization_epitopePercentPresent.pdf"), plot = ggOut, width = 20, height = 20)

#Decreasing with urbanization
ggOut <- ggplot(subset(plotDF, variable %in% subset(testSigDF, Z > 0)$ID), aes(x = Residence, y = value, fill = Residence)) +
  geom_bar(stat = "identity") +
  facet_wrap(~variable) +
  scale_fill_manual(values = area_palette) +
  ylim(0, 100) +
  theme_minimal()
ggOut
#ggsave(filename = paste0("plots/figure_components/sup_figures/cochranArmitage_decreaseWithUrbanization_epitopePercentPresent.pdf"), plot = ggOut, width = 20, height = 20)


#Plot the top three enriched pathogens as the normalised frequencies
plotDF <- normTotalEpitope_DF[,-which(colnames(normTotalEpitope_DF) == "epitopeTotals")]

top3 <- c("humanherpesvirus1TotalNorm", "shigellasonneiTotalNorm", "shigellaflexneriTotalNorm",
          "haemophilusinfluenzaeTotalNorm", "epsteinbarrvirusTotalNorm", "rhinovirusbTotalNorm")

plotDF <- melt(plotDF[, c(top3, "Residence")])

plotDF$Residence <- factor(plotDF$Residence, levels = c("Urban Dutch", "Urban Senegalese", "Rural Senegalese"))

ggOut <- ggplot(plotDF, aes(x = value, y = Residence, fill = Residence)) +
  geom_boxplot() +
  geom_jitter() +
  theme_minimal() +
  facet_wrap(~variable, scales = "free") +
  scale_fill_manual(values = area_palette) +
  xlab("Proportion of Ab-OME against organism epitopes") +
  ylab("")
ggOut
#ggsave(filename = paste0("plots/figure_components/top3_enrichment_urbanization_normPathogen.pdf"), plot = ggOut, width = 14)


#### Investigate rural, urban and dutch: Specific protein binding differences ####

phipSeq_p_freq_residence_raw_df$Residence <- row.names(phipSeq_p_freq_residence_raw_df)
plotDF <- melt(phipSeq_p_freq_residence_raw_df)

plotDF$Residence <- factor(plotDF$Residence, rev(c("Rural Senegalese", "Urban Senegalese", "Urban Dutch")) )
normTotalEpitope_DF$Residence <- factor(normTotalEpitope_DF$Residence, rev(c("Rural Senegalese", "Urban Senegalese", "Urban Dutch")) )

#Plot 6 interesting epitopes that are specific for the different residences
epitopeIDs <- c("agilent_219967", "agilent_204192", "agilent_104875","agilent_198635",
                "agilent_214358", "agilent_150541", "twist_37139", "agilent_227745")

epitopeNames <- c("GRAB (<i>S. pyogenes</i>)", "Uncharacterized Protein (<i>P. dorei</i>)", "Bifunctional metallophosphatase/5'-nucleotidase (<i>P. vulgatus</i>)", "Exotoxin 5 (<i>S. aureus</i>)",
                  "Exported protein (<i>S. aureus</i>)", "Port family protein (<i>P. vulgatus</i>)", "Albumin-binding GA domain-containing protein, partial (<i>S. dysgalactiae</i>)", "Adhesin (<i>H. influenzae</i>)")

names(epitopeNames) <- epitopeIDs

selectedEpitope_plotDF <- subset(plotDF, variable %in% epitopeIDs)

selectedEpitope_plotDF$variable <- as.character(selectedEpitope_plotDF$variable)

selectedEpitope_plotDF$epitopeName <- epitopeNames[selectedEpitope_plotDF$variable]

selectedEpitope_plotDF$epitopeName <- factor(selectedEpitope_plotDF$epitopeName, levels = rev( as.character(epitopeNames) ) )

ggplot(subset(plotDF, variable == "agilent_2709"), aes(x = Residence, y = value, fill = Residence)) +
  geom_bar(stat = "identity")  +
  scale_fill_manual(values = area_palette) +
  ylim(0, 100) +
  theme_minimal()

ggplot(subset(plotDF, variable == "twist_42197"), aes(x = Residence, y = value, fill = Residence)) +
  geom_bar(stat = "identity")  +
  scale_fill_manual(values = area_palette) +
  ylim(0, 100) +
  theme_minimal()


ggplot(selectedEpitope_plotDF, aes(x = Residence, y = value, fill = Residence)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = area_palette) + 
  facet_wrap(~epitopeName) +
  ylim(0, 100) +
  theme_minimal()

ggOut <- ggplot(selectedEpitope_plotDF, aes(x = value, y = epitopeName, fill = Residence, group = Residence)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = area_palette)  +
  xlim(0, 100) +
  theme_minimal() +
  xlab("% of individuals where Ab detected") +
  ylab("") +
  theme(axis.text.y = ggtext::element_markdown())
ggOut
ggsave(filename = "plots/figure_components/epitopes_across_urbanization_barplot.pdf", plot = ggOut, 
       width = 11, height = 6)

#Both Strep pyogenes and staph aureus have epitopes that are either increasing or decreasing with urbanization



ggplot(normTotalEpitope_DF, aes(x= streptococcuspyogenesTotalNorm ,y= Residence, fill = Residence )) +
          geom_boxplot() +
          geom_jitter() +
          scale_fill_manual(values = area_palette) +
          theme_minimal()+
  ylab("") +
  xlab("Proportion of Ab against Streptococcus pyogenes epitopes")

ggplot(normTotalEpitope_DF, aes(x= staphylococcusaureusTotalNorm ,y= Residence, fill = Residence )) +
  geom_boxplot() +
  geom_jitter() +
  scale_fill_manual(values = area_palette) +
  theme_minimal() +
  ylab("")+
  xlab("Proportion of Ab against Staphylococcus aureus epitopes")



#### Investigate rural, urban and dutch: Shannon index ####

#Get the occurrence of each organism, grouped by Residence
phipSeq_p_raw_df$Residence <- metaData[row.names(phipSeq_p_raw_df), "Residence"]

organismCounts <- phipSeq_p_raw_df

organismCounts$Donor <- row.names(organismCounts)

organismCounts_melt <- reshape2::melt(organismCounts)

organismCounts_melt$variable <- as.character(organismCounts_melt$variable)

#Add pathogen name
organismCounts_melt$Pathogen <- peptideInfo[organismCounts_melt$variable,"Pathogen"]

#Calculate pathogen sums for each Residence
organismCounts_Donor <-  organismCounts_melt %>%
                    group_by(Donor, Pathogen) %>%
                    summarise(organismCount_Donor = sum(value))
organismCounts_Donor <- as.data.frame(organismCounts_Donor)

organismCounts_melt <- organismCounts_melt %>%
                        group_by(Donor) %>%
                        summarise(organismCount = sum(value))
organismCounts_melt <- as.data.frame(organismCounts_melt)


row.names(organismCounts_melt) <- organismCounts_melt$Donor

organismCounts_Donor$organismCount <- organismCounts_melt[organismCounts_Donor$Donor, "organismCount"]

#Remove organisms which are present in less than 3 individuals
organismCounts_Donor <- subset(organismCounts_Donor, organismCount_Donor > 0)

organismCounts_Donor$organismProportions <- organismCounts_Donor$organismCount_Donor / organismCounts_Donor$organismCount 

#Add small pseudocount
# organismCounts_Donor$organismProportions <- organismCounts_Donor$organismProportions + 0.0000000000000001

shannonDiversity <- organismCounts_Donor %>%
                          dplyr::group_by(Donor) %>%
                          summarise(shannonDiv = -sum(organismProportions * log(organismProportions) ))

shannonDiversity$Residence <- metaData[shannonDiversity$Donor, "Residence"]
shannonDiversity$Residence <- factor(shannonDiversity$Residence, levels = c("Rural Senegalese", "Urban Senegalese", "Urban Dutch"))


ggplot(shannonDiversity, aes(x = Residence, y = shannonDiv, fill = Residence)) +
            geom_boxplot() +
            geom_jitter() +
            theme_minimal() +
            scale_fill_manual(values = area_palette)

#Relationship with total epitopes
shannonDiversity$totalEpitopes <- rowSums(phipSeq_p_raw)[shannonDiversity$Donor]

ggplot(shannonDiversity, aes(x = totalEpitopes, y = shannonDiv, color = Residence)) +
  geom_point() +
  theme_minimal() +
  scale_color_manual(values = area_palette)

ggplot(shannonDiversity, aes(x = Residence, y = totalEpitopes, fill = Residence)) +
  geom_boxplot() +
  geom_jitter() +
  theme_minimal() +
  scale_fill_manual(values = area_palette)

#### Investigate rural and urban and dutch: microbiome ####

#Bacteriodes and prevotella

normTotalEpitope_microbiome_DF <- data.frame(epitopeTotals = rowSums(phipSeq_p_raw[metaData$ID,]),
                                  row.names = metaData$ID)

pathogens <- c("^bacteroides|phocaeicola vulgatus", "segatella copri")
pathogenNames <- c("bacteroides", "prevotella")

for(i in 1:length(pathogens)){
  
  peptideInfoCurrent <- peptideInfo[str_detect(peptideInfo$Pathogen, pathogens[i]),]
  
  currentEpitopes <- row.names(peptideInfoCurrent)
  
  
  
  #Add current epitopes bound to data frame
  if(length(currentEpitopes)==1){
    currentEpitopesTotals <- phipSeq_p_raw[,currentEpitopes]
    
  }
  
  
  else{
    currentEpitopesTotals <- rowSums(phipSeq_p_raw[,currentEpitopes])
    
  }

  
  normTotalEpitope_microbiome_DF[[paste0(pathogenNames[i], "TotalNorm")]] <- currentEpitopesTotals[row.names(normTotalEpitope_microbiome_DF)] / normTotalEpitope_microbiome_DF$epitopeTotals
  
}

normTotalEpitope_microbiome_DF$Residence <- metaData[row.names(normTotalEpitope_microbiome_DF), "Residence"]
normTotalEpitope_microbiome_DF$Residence <- factor(normTotalEpitope_microbiome_DF$Residence, levels = c("Rural Senegalese", "Urban Senegalese", "Urban Dutch"))
normTotalEpitope_microbiome_DF$Donor <- row.names(normTotalEpitope_microbiome_DF)
normTotalEpitope_microbiome_DF$Sex <- metaData[row.names(normTotalEpitope_microbiome_DF), "Sex"]

ggplot(normTotalEpitope_microbiome_DF, aes(x = bacteroidesTotalNorm, y = Residence, fill = Residence)) +
          geom_boxplot() +
          geom_jitter() + 
          theme_minimal() +
          scale_fill_manual(values = area_palette)

ggplot(normTotalEpitope_microbiome_DF, aes(x = prevotellaTotalNorm, y = Residence, fill = Residence)) +
  geom_boxplot() +
  geom_jitter() + 
  theme_minimal() +
  scale_fill_manual(values = area_palette)

test_normTotal <- normTotalEpitope_DF
test_normTotal <- test_normTotal[,str_detect(colnames(test_normTotal), "TotalNorm")]

rowSums(test_normTotal)

#### Investigate Rural and urban senegalese differences ####

metaData_complex_senegal <- subset(metaData_complex, Residence != "Urban Dutch")

table(metaData_complex_senegal$Residence, metaData_complex_senegal$job_amalgam)
table(metaData_complex_senegal$Sex, metaData_complex_senegal$job_amalgam)

phipSeq_p_senegal <- phipSeq_p_raw[senegalIndividuals,names(epitopeByResidenceSenegal) ]

#Remove columns where everyone has the antibody
phipSeq_p_senegal <- phipSeq_p_senegal[,colSums(phipSeq_p_senegal) != nrow(phipSeq_p_senegal)]

# #Remove people whose working/job environment is not indoors or outdoors
# phipSeq_p_senegal <- phipSeq_p_senegal[subset(metaData_complex_senegal, job_amalgam %in% c("indoor", "outdoor") )$ID ,]

PCA_phip <- prcomp(phipSeq_p_senegal, scale = TRUE, center = TRUE)

fviz_pca_ind(PCA_phip)

embedPC <- as.data.frame(PCA_phip$x)
embedPC$Residence <- as.character(metaData[row.names(embedPC), "Residence"])
embedPC$sampleID <- row.names(embedPC)
embedPC$Sex <- metaData[row.names(embedPC), "Sex"]
embedPC$Age <- metaData[row.names(embedPC), "Age"]
embedPC$job_amalgam <- metaData_complex[row.names(embedPC), "job_amalgam"]
embedPC$Residence_complex <- metaData_complex[row.names(embedPC), "Cur_Res_G"]


ggplot(embedPC, aes(x = PC1, y = PC2, color = Residence)) +
  geom_point(size = 2) +
  theme_minimal() +
  scale_color_manual(values = area_palette)

#Get variance explained
PC1_varExplained <- (PCA_phip$sdev[1]^2 / sum(PCA_phip$sdev^2)) * 100
PC2_varExplained <- (PCA_phip$sdev[2]^2 / sum(PCA_phip$sdev^2)) * 100
PC10_varExplained <- (PCA_phip$sdev[10]^2 / sum(PCA_phip$sdev^2)) * 100

embedPC$Residence <- factor(embedPC$Residence, levels = c("Rural Senegalese", "Urban Senegalese"))

#PCA plot coloured by Area
gg <- merge(embedPC,aggregate(cbind(mean.PC1=PC1,mean.PC2=PC2)~Residence,embedPC,mean),by="Residence")
ggOut <- ggplot(gg, aes(PC1,PC2))+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0) +
  geom_point(size=3, aes(color=Residence))+
  geom_segment(aes(x=mean.PC1, y=mean.PC2, xend=PC1, yend=PC2, color = Residence)) +
  geom_point(aes(x=mean.PC1,y=mean.PC2, fill=Residence), size=5, shape = 21) +
  theme_minimal()+
  scale_fill_manual(values = area_palette) +
  scale_color_manual(values = area_palette)+
  xlab(paste0("PC1 (", round(PC1_varExplained, digits = 2), "%)"))+
  ylab(paste0("PC2 (", round(PC2_varExplained, digits = 2), "%)"))
ggOut
ggsave(filename = "plots/figure_components/PC1PC2Plot_senegal_filteredData_residenceColour.pdf", plot = ggOut, 
       width = 9, height = 8)

#PCA plot coloured by job_amalgam
gg <- merge(embedPC,aggregate(cbind(mean.PC1=PC1,mean.PC2=PC2)~job_amalgam,embedPC,mean),by="job_amalgam")
ggOut <- ggplot(gg, aes(PC1,PC2))+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0) +
  geom_point(size=3, aes(color=job_amalgam))+
  geom_segment(aes(x=mean.PC1, y=mean.PC2, xend=PC1, yend=PC2, color = job_amalgam)) +
  geom_point(aes(x=mean.PC1,y=mean.PC2, fill=job_amalgam), size=5, shape = 21) +
  theme_minimal() +
  xlab(paste0("PC1 (", round(PC1_varExplained, digits = 2), "%)"))+
  ylab(paste0("PC2 (", round(PC2_varExplained, digits = 2), "%)"))
ggOut

#PCA plot coloured by Sex
gg <- merge(embedPC,aggregate(cbind(mean.PC1=PC1,mean.PC2=PC2)~Sex,embedPC,mean),by="Sex")
ggOut <- ggplot(gg, aes(PC1,PC2))+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0) +
  geom_point(size=3, aes(color=Sex))+
  geom_segment(aes(x=mean.PC1, y=mean.PC2, xend=PC1, yend=PC2, color = Sex)) +
  geom_point(aes(x=mean.PC1,y=mean.PC2, fill=Sex), size=5, shape = 21) +
  theme_minimal() +
  xlab(paste0("PC1 (", round(PC1_varExplained, digits = 2), "%)"))+
  ylab(paste0("PC2 (", round(PC2_varExplained, digits = 2), "%)"))
ggOut

#PCA plot coloured by residence, split by job_amalgam
gg <- merge(embedPC,aggregate(cbind(mean.PC1=PC1,mean.PC2=PC2)~Residence,embedPC,mean),by="Residence")
ggOut <- ggplot(gg, aes(PC1,PC2))+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0) +
  geom_point(size=3, aes(color=Residence))+
  geom_segment(aes(x=mean.PC1, y=mean.PC2, xend=PC1, yend=PC2, color = Residence)) +
  geom_point(aes(x=mean.PC1,y=mean.PC2, fill=Residence), size=5, shape = 21) +
  theme_minimal() +
  xlab(paste0("PC1 (", round(PC1_varExplained, digits = 2), "%)"))+
  ylab(paste0("PC2 (", round(PC2_varExplained, digits = 2), "%)")) +
  facet_wrap(~job_amalgam)+
  scale_fill_manual(values = area_palette) +
  scale_color_manual(values = area_palette)
ggOut

gg <- merge(embedPC,aggregate(cbind(mean.PC1=PC1,mean.PC2=PC2)~Residence,embedPC,mean),by="Residence")
ggOut <- ggplot(gg, aes(PC1,PC2))+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0) +
  geom_point(size=3, aes(color=Sex))+
  geom_segment(aes(x=mean.PC1, y=mean.PC2, xend=PC1, yend=PC2, color = Residence)) +
  geom_point(aes(x=mean.PC1,y=mean.PC2, fill=Residence), size=5, shape = 21) +
  theme_minimal() +
  xlab(paste0("PC1 (", round(PC1_varExplained, digits = 2), "%)"))+
  ylab(paste0("PC2 (", round(PC2_varExplained, digits = 2), "%)")) +
  facet_wrap(~job_amalgam)+
  scale_fill_manual(values = area_palette) +
  scale_color_manual(values = append(area_palette, c("M" = "red", "F" = "blue")))
ggOut

gg <- merge(embedPC,aggregate(cbind(mean.PC1=PC1,mean.PC2=PC2)~Residence_complex,embedPC,mean),by="Residence_complex")
ggOut <- ggplot(gg, aes(PC1,PC2))+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0) +
  geom_point(size=3, aes(color=Sex))+
  geom_segment(aes(x=mean.PC1, y=mean.PC2, xend=PC1, yend=PC2, color = Residence_complex)) +
  geom_point(aes(x=mean.PC1,y=mean.PC2, fill=Residence_complex), size=5, shape = 21) +
  theme_minimal() +
  xlab(paste0("PC1 (", round(PC1_varExplained, digits = 2), "%)"))+
  ylab(paste0("PC2 (", round(PC2_varExplained, digits = 2), "%)")) +
  facet_wrap(~job_amalgam)+
  scale_fill_manual(values = areaComplex_palette) +
  scale_color_manual(values = append(areaComplex_palette, c("M" = "red", "F" = "blue")))
ggOut

#Perform chi-sq test to compare epitopes between rural and urban senegalese

phipSeq_p_senegal_df <- as.data.frame(phipSeq_p_senegal)

#Remove columns where everyone has the antibody

phipSeq_p_senegal_df$Residence <- metaData[row.names(phipSeq_p_senegal_df),"Residence"]

convertBinary <- reshape2::melt(phipSeq_p_senegal_df)
convertBinary$variable <- as.character(convertBinary$variable)
convertBinary$Residence <- factor(convertBinary$Residence, levels =  c("Rural Senegalese" , "Urban Senegalese"))

fishersDF <- data.frame()

for(i in unique(convertBinary$variable)){
  current <- subset(convertBinary, variable == i)
  
  testResult <- fisher.test( table(current$Residence, current$value) )
  
  fishersDF <- rbind(fishersDF, data.frame("ID" = i,
                                           "p.value" = testResult$p.value,
                                           "Z" = testResult$estimate))
  
}

#Add meta data
fishersDF$Pathogen <- peptideInfo[fishersDF$ID, "Pathogen"]
fishersDF$PathogenComplex <- peptideInfo[fishersDF$ID, "PathogenComplex"]
fishersDF$Protein <- peptideInfo[fishersDF$ID, "Protein"]

fishersDF$p.adjust <- p.adjust(fishersDF$p.value, method = "BH")

fishersSigDF <- subset(fishersDF, p.adjust < 0.05)
fishersSigDF <- subset(fishersDF, p.value < 0.05)

ggplot() + geom_point(aes(x= fishersDF$Z,y= fishersDF$p.value)) +
  geom_vline(xintercept = median(fishersDF$Z))

freqPercent <- as.data.frame(t(phipSeq_p_freq_residence_df[c("Rural Senegalese" , "Urban Senegalese"),fishersSigDF$ID]))
freqPercent$ResidenceDiff <- freqPercent$`Rural Senegalese` - freqPercent$`Urban Senegalese`

ruralEpitopes <- table(peptideInfo[row.names(subset(freqPercent, ResidenceDiff > 0)), "Pathogen"])
urbanEpitopes <- table(peptideInfo[row.names(subset(freqPercent, ResidenceDiff < 0)), "Pathogen"])

#Normalise for how many times an epitope of that organism is present in the dataset
pathogenTotals <- table(peptideInfo$Pathogen)

ruralEpitopes <- ruralEpitopes / pathogenTotals[names(ruralEpitopes)]
urbanEpitopes <- urbanEpitopes / pathogenTotals[names(urbanEpitopes)]

wordcloud2(ruralEpitopes, size = 0.6)
wordcloud2(urbanEpitopes, size = 0.8)

#Look at the significant results
phipSeq_p_freq_residence_df$Residence <- row.names(phipSeq_p_freq_residence_df)
plotDF <- melt(phipSeq_p_freq_residence_df[c("Rural Senegalese" , "Urban Senegalese"),c(row.names(subset(freqPercent, ResidenceDiff > 0)), "Residence")] )
plotDF$Residence <- factor(plotDF$Residence, c("Rural Senegalese", "Urban Senegalese") )

#What pathogens have different epitopes present in different groups
ggplot(plotDF, aes(x = Residence, y = value, fill = Residence)) +
  geom_bar(stat = "identity") +
  facet_wrap(~variable) +
  scale_fill_manual(values = area_palette) +
  ylim(0, 100)

plotDF <- melt(phipSeq_p_freq_residence_df[c("Rural Senegalese" , "Urban Senegalese"),c(row.names(subset(freqPercent, ResidenceDiff < 0)), "Residence")] )
plotDF$Residence <- factor(plotDF$Residence, c("Rural Senegalese", "Urban Senegalese") )

#What pathogens have different epitopes present in different groups
ggplot(plotDF, aes(x = Residence, y = value, fill = Residence)) +
  geom_bar(stat = "identity") +
  facet_wrap(~variable) +
  scale_fill_manual(values = area_palette) +
  ylim(0, 100)

#### Predict Rural vs urban senegalese ####







library(randomForest)
library(caret)

x <- phipSeq_p_senegal_df

#Add current residence area
x$currentResidence <- metaData_complex_senegal[, ""]

x$Residence <- factor(x$Residence)

{set.seed(42); trainTest <- sample(2, nrow(x), replace = TRUE, prob = c(0.7, 0.3)) }

train <- x[trainTest==1,]
test <- x[trainTest==2,]

#Random forest
rfResults <- randomForest(Residence~., data=train)

print(rfResults)

predictions_test_rf <- as.data.frame(predict(rfResults, test, type = "prob"))

predictions_train_rf <- as.data.frame(predict(rfResults, type = "prob"))

rf.roc<-pROC::roc(test$Residence,predictions_test[,2])
plot(rf.roc)

#Logisitic regression

lrResults <- glm(Residence ~.,family=binomial(link='logit'),data=train)

predictions_test_lr <- as.data.frame(predict(lrResults, test,type='response'))
predictions_test_lr <- ifelse(predictions_test_lr > 0.5,1,0)

misClasificError <- mean(predictions_test_lr != test$Residence)

#Decision tree

library(partykit)

dtResults <- ctree(Residence ~., 
            data=train)
plot(dtResults)

table(predict(dtResults), train$Residence)

predictions_test_dt <- predict(dtResults, newdata = test, type = "prob")



