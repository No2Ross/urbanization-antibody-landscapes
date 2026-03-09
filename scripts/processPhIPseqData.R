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

#For the cochran trend test
library(DescTools)

library(wordcloud2)

library(factoextra)

#library(beer)

setwd("C:/Users/rflaidlaw/Documents/CapTan/AnalysisV2/PhIPseq/")

area_palette <- c("Rural Senegalese" = "#e59f01", "Urban Senegalese" = "#54b3e8", "Urban Dutch" = "#009e72")

metaData <- as.data.frame(read_xlsx("phipseq_metadata.xlsx"))
row.names(metaData) <- metaData$ID

#Make age groupings
#0-17, 18-29, 30-45, 46+

# phipSeq_FC <- read.csv("res_foldchange_annotated.csv", row.names = 2)
phipSeq_p <- read.csv("res_padj_annotated.csv", row.names = 2)

peptideInfo <- read.csv("PhIPseq_peptideInfo_processed_BLAST.csv", row.names = 1)

senegalIndividuals <- metaData[str_detect(metaData$Residence, "enegal"), "ID"]

#### Transform and subset data ####

sampleNames <- colnames(phipSeq_p)[8:ncol(phipSeq_p)]

#Prepare data
phipSeq_p <- t(as.matrix( phipSeq_p[, sampleNames]) )

#Change the row names to reflect the metaData
row.names(phipSeq_p) <- unlist(str_split( row.names(phipSeq_p), "_" ))[c(FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE)]

#Convert NA into 0 and anything else into 1
phipSeq_p[is.na(phipSeq_p)] <- 0
phipSeq_p <- sign(phipSeq_p)

#Remove epitopes which are not present in any individuals
ncol(phipSeq_p)
phipSeq_p <- phipSeq_p[, colSums(phipSeq_p) > 0]
ncol(phipSeq_p)

#Plot peptide variation and keep only certain peptides
plot( colMeans( phipSeq_p ), colVars( phipSeq_p ))

#Keep peptides which are found in 25% of at least one of the three residences
phipSeq_p_df <- as.data.frame(phipSeq_p)
phipSeq_p_df$Residence <- metaData[row.names(phipSeq_p_df), "Residence"]
phipSeq_p_df$Sex <- metaData[row.names(phipSeq_p_df), "Sex"]


#Calculate frequencies for residence
phipSeq_p_count <- phipSeq_p_df[, -which(colnames(phipSeq_p_df) == "Sex")] %>%
  dplyr::group_by(Residence) %>%
  summarise(across(everything(),sum)) %>%
  as.data.frame()

row.names(phipSeq_p_count) <- phipSeq_p_count$Residence
phipSeq_p_count$Residence <- NULL

residenceTotals <- table(metaData$Residence)
phipSeq_p_freq_residence <- (phipSeq_p_count / residenceTotals[row.names(phipSeq_p_count)]) * 100

#Calculate frequencies for sex (only for senegal)
phipSeq_p_count <- phipSeq_p_df[senegalIndividuals, -which(colnames(phipSeq_p_df) == "Residence")] %>%
  dplyr::group_by(Sex) %>%
  summarise(across(everything(),sum)) %>%
  as.data.frame()

row.names(phipSeq_p_count) <- phipSeq_p_count$Sex
phipSeq_p_count$Sex <- NULL

SexTotals <- table(metaData$Sex)
phipSeq_p_freq_sex <- (phipSeq_p_count / SexTotals[row.names(phipSeq_p_count)]) * 100

#Calculate frequencies for combination of sex and residence
phipSeq_p_count <- phipSeq_p_df %>%
  dplyr::group_by(Sex, Residence) %>%
  summarise(across(everything(),sum)) %>%
  as.data.frame()

row.names(phipSeq_p_count) <- paste0(phipSeq_p_count$Residence,"_",phipSeq_p_count$Sex)
phipSeq_p_count$Residence <- NULL
phipSeq_p_count$Sex <- NULL

SexResidenceTotals <- table(paste0(metaData$Residence,"_",metaData$Sex))
phipSeq_p_freq_sexResidence <- (phipSeq_p_count / SexResidenceTotals[row.names(phipSeq_p_count)]) * 100

#### Create normalised total epitopes bound for each pathogen ####

phipSeq_p_pathogenTotals <- reshape2::melt(phipSeq_p)
phipSeq_p_pathogenTotals$Var2 <- as.character(phipSeq_p_pathogenTotals$Var2)
phipSeq_p_pathogenTotals$Var1 <- as.character(phipSeq_p_pathogenTotals$Var1)

#Add pathogen information
phipSeq_p_pathogenTotals$Pathogen <- peptideInfo[phipSeq_p_pathogenTotals$Var2, "Pathogen"]

phipSeq_p_pathogenTotals <- phipSeq_p_pathogenTotals %>%
        dplyr::group_by(Var1, Pathogen) %>%
        summarise(totalEpitopes = sum(value))

phipSeq_p_pathogenTotalsMtx <- reshape2::dcast(phipSeq_p_pathogenTotals, formula = Var1 ~ Pathogen, value.var = "totalEpitopes")
row.names(phipSeq_p_pathogenTotalsMtx) <- phipSeq_p_pathogenTotalsMtx$Var1
phipSeq_p_pathogenTotalsMtx$Var1 <- NULL

#### Save out datasets ####

#rows are individuals, cols are epitope IDs. Values are whether an individual has antibodies against that epitope
#rows are residences, cols are epitope IDs. Values are the percentage of individuals in that residence that have antibodies against the epitope
#rows are sex (no dutch), cols are epitope IDs. Values are the percentage of individuals for that Sex that have antibodies against the epitope
#rows are combination of residence and Sex, cols are epitope IDs. Values are the percentage of individuals for that Sex-residence combination that have antibodies against the epitope

write.csv(phipSeq_p,"data/individualAll_epitope_binary.csv")
write.csv(phipSeq_p_freq_residence,"data/residenceAll_epitope_percentage.csv")
write.csv(phipSeq_p_freq_sex,"data/sexSenegal_epitope_percentage.csv")
write.csv(phipSeq_p_freq_sexResidence,"data/sexResidenceAll_epitope_percentage.csv")
write.csv(phipSeq_p_pathogenTotalsMtx,"data/pathogenEpitopeTotals.csv")
