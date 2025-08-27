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

dbInfo <- as.data.frame(read_xlsx("Clean GDIR Database 08102020.xlsx"))
#Some rows have the same ID

dbInfo <- dbInfo[which(dbInfo$ID %in% row.names(metaData)),]
row.names(dbInfo) <- dbInfo$ID
table(dbInfo$ID)

#Remove columns that aren't of interest
dbInfo$Date <- NULL
dbInfo$Sex_G <- NULL
dbInfo$Date_birth <- NULL

#Remove columns where every value is na
keepCols <- colSums(is.na(dbInfo))
keepCols <- keepCols[keepCols < nrow(dbInfo)]

dbInfo <- dbInfo[, names(keepCols)]

#The past residence 2, 3 and 4 columns are only relevant for the dutch, so remove them
dbInfo$Past_Ressid2 <- NULL
dbInfo$Past_Ressid3 <- NULL
dbInfo$Past_Ressid4 <- NULL
dbInfo$Dur_past_Ressid2 <- NULL
dbInfo$Dur_past_Ressid3 <- NULL
dbInfo$Dur_past_Ressid4 <- NULL

#Convert all string columns to lowercase
dbInfo <- dbInfo %>%
          mutate(across(where(is.character), ~tolower(.x)))

#I want ID and Sex to remain uppercase
dbInfo$ID <- toupper(dbInfo$ID)
dbInfo$Sex <- toupper(dbInfo$Sex)


#Change the place of birth, growing up and residence from numeric to character encoding
dbInfo$Place_birth_G <- as.character(dbInfo$Place_birth_G)
dbInfo$Place_birth_G <- c("1" = "Rural Senegalese", "2" = "Semi-Urban Senegalese", "3" = "Urban Senegalese", "4" = "Urban Dutch")[dbInfo$Place_birth_G]

dbInfo$Place_grow_G <- as.character(dbInfo$Place_grow_G)
dbInfo$Place_grow_G <- c("1" = "Rural Senegalese", "2" = "Semi-Urban Senegalese", "3" = "Urban Senegalese", "4" = "Urban Dutch")[dbInfo$Place_grow_G]

dbInfo$Cur_Res_G <- as.character(dbInfo$Cur_Res_G)
dbInfo$Cur_Res_G <- c("1" = "Rural Senegalese", "2" = "Semi-Urban Senegalese", "3" = "Urban Senegalese", "4" = "Urban Dutch")[dbInfo$Cur_Res_G]

#Add in Mikhaels Residence designation column
dbInfo$Residence <- metaData[dbInfo$ID,"Residence"]

#Where Residence is "urban Dutch", change Place_birth_G to "Not Senegal"
dbInfo$Place_birth_G[which(dbInfo$Residence == "Urban Dutch")] <- "Not Senegal"

#Calculate neutrophil to lymphocyte ratio
dbInfo$NLR <- dbInfo$Neutro_ab / dbInfo$Lymp_ab

metaData$NLR <- dbInfo[metaData$ID, "NLR"]

test <- dbInfo[, c("Residence", "Place_birth_G", "Place_grow_G")]

#Working conditions: "outdoor" is sometimes spelled "outdor"
dbInfo$Working_Env[str_detect(dbInfo$Working_Env, "outdor")] <- "outdoor"

#Create a new column which combines the job/working environment/jobless columns
#If job or jobless type is housewife, the column value becomes "housewife"
#If they don't have a working environment column, set them to unknown

dbInfo$job_amalgam <- dbInfo$Working_Env
dbInfo$job_amalgam[str_detect(dbInfo$Job_type, "housewife")] <- "housewife"
dbInfo$job_amalgam[str_detect(dbInfo$Occup_if_jobless, "housewife")] <- "housewife"
dbInfo$job_amalgam[is.na(dbInfo$job_amalgam) == TRUE] <- "unknown"

dbInfo$job_amalgam[str_detect(dbInfo$job_amalgam, "outside")] <- "outdoor"
dbInfo$job_amalgam[str_detect(dbInfo$job_amalgam, "outdoor")] <- "outdoor"

dbInfo$job_amalgam[str_detect(dbInfo$job_amalgam, "indoor")] <- "indoor"


write.csv(metaData, "phipseq_metadata.csv")
write.csv(dbInfo, "phipseq_metadata_extended.csv")
writexl::write_xlsx(dbInfo, "phipseq_metadata_extended.xlsx")
