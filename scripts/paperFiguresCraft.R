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
library(grid)
library(ggridges)

library(ape)
library(prabclus)

#For the cochran trend test
library(DescTools)

library(wordcloud2)

library(factoextra)

library(vegan)

library(emmeans)

#Correlation plots
library(PerformanceAnalytics)

#library(beer)

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


setwd("C:/Users/rflaidlaw/Documents/CapTan/AnalysisV2/PhIPseq/")

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


# phipSeq_FC <- read.csv("res_foldchange_annotated.csv", row.names = 2)
phipSeq_count <- read.csv("res_counts_annotated.csv", row.names = 2)
sampleNames <- colnames(phipSeq_count)[8:ncol(phipSeq_count)]
phipSeq_count <- t(as.matrix( phipSeq_count[, sampleNames]) )
#Change the column names to reflect the metaData
row.names(phipSeq_count) <- unlist(str_split( row.names(phipSeq_count), "_" ))[c(FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE)]

peptideInfo <- read.csv("PhIPseq_peptideInfo_processed_BLAST.csv", row.names = 1)

phipSeq_p <- read.csv("data/individualAll_epitope_binary.csv", row.names = 1)
phipSeq_p_freq_residence <- read.csv("data/residenceAll_epitope_percentage.csv", row.names = 1)
phipSeq_p_freq_sex <- read.csv("data/sexSenegal_epitope_percentage.csv", row.names = 1)
phipSeq_p_freq_sexResidence <- read.csv("data/sexResidenceAll_epitope_percentage.csv", row.names = 1)
phipSeq_p_pathogenTotal <- read.csv("data/pathogenEpitopeTotals.csv", row.names = 1)

#### Filter features ####

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

# #This code is quite convoluted, but it basically just finds what epitopes are present in 25% of individuals in at least one residence group
# epitopeByResidence["Rural Senegalese",] <- (epitopeByResidence["Rural Senegalese",] / nrow(subset(metaData, Residence == "Rural Senegalese"))) * 100
# epitopeByResidence["Urban Senegalese",] <- (epitopeByResidence["Urban Senegalese",] / nrow(subset(metaData, Residence == "Urban Senegalese"))) * 100
# epitopeByResidence["Urban Dutch",] <- (epitopeByResidence["Urban Dutch",] / nrow(subset(metaData, Residence == "Urban Dutch"))) * 100
# 
# percentageCutoff <- 25
# 
# epitopeByResidence <- (epitopeByResidence + 0.000000000000000000000000000001) / percentageCutoff
# epitopeByResidence <- log(epitopeByResidence)
# epitopeByResidence <- sign(epitopeByResidence)
# epitopeByResidence <- colSums(epitopeByResidence)
# epitopeByResidence <- epitopeByResidence[epitopeByResidence !=-3 ]
# 
# subsetPeptideInfo <- peptideInfo[names(epitopeByResidence),]
# 
# dim(subsetPeptideInfo)

#This code keeps epitopes that are found in 4 or more individuals
epitopeByResidence <- epitopeByResidence / 4
epitopeByResidence <- log(epitopeByResidence)

#This section rests on log(0) = -Inf
epitopeByResidence <- sign(epitopeByResidence)

epitopeByResidence <- colSums(epitopeByResidence)
epitopeByResidence <- epitopeByResidence[epitopeByResidence !=-3 ]

subsetPeptideInfo <- peptideInfo[names(epitopeByResidence),]

dim(subsetPeptideInfo)

#### Define objects ####

#Remove the unknown BLAST epitopes
known_epitopes <- row.names(subset(peptideInfo, Pathogen != "unknownblast"))

phipSeq_p_raw_df <- phipSeq_p[,known_epitopes]
phipSeq_p_freq_residence_raw_df <- phipSeq_p_freq_residence
phipSeq_p_freq_sex_raw_df <- phipSeq_p_freq_sex
phipSeq_p_freq_sexResidence_raw_df <- phipSeq_p_freq_sexResidence

phipSeq_p_raw <- as.matrix(phipSeq_p[,known_epitopes])
phipSeq_p_freq_residence_raw <- as.matrix(phipSeq_p_freq_residence)
phipSeq_p_freq_sex_raw <- as.matrix(phipSeq_p_freq_sex)
phipSeq_p_freq_sexResidence_raw <- as.matrix(phipSeq_p_freq_sexResidence)

epitopeByResidence <- epitopeByResidence[intersect(known_epitopes, names(epitopeByResidence))]

phipSeq_p_df <- phipSeq_p_raw_df[,names(epitopeByResidence)]
phipSeq_p_freq_residence_df <- phipSeq_p_freq_residence_raw_df[,names(epitopeByResidence)]
phipSeq_p_freq_sex_df <- phipSeq_p_freq_sex_raw_df[,names(epitopeByResidence)]
phipSeq_p_freq_sexResidence_df <- phipSeq_p_freq_sexResidence_raw_df[,names(epitopeByResidence)]

phipSeq_p <- as.matrix(phipSeq_p_df)
phipSeq_p_freq_residence <- as.matrix(phipSeq_p_freq_residence_df)
phipSeq_p_freq_sex <- as.matrix(phipSeq_p_freq_sex_df)
phipSeq_p_freq_sexResidence <- as.matrix(phipSeq_p_freq_sexResidence_df)

peptideInfo <- subset(peptideInfo, Pathogen != "unknownblast")

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
ggsave(filename = "//vf-lucid-r-o.lumcnet.prod.intern/lucid-r-o$/Projects/Senegal_phipseq/plotsCell/epitopeCountGrouping.pdf", plot = ggOut, 
       width =12, height = 8)

plotDF <- data.frame("epitopeCount" = epitopeCount)
ggOut <- ggplot(plotDF, aes(x = epitopeCount)) +
  geom_histogram() +
  theme_minimal() +
  xlab("No. individuals in which anti-epitope is present") + 
  ylab("No. of significant anti-epitopes (log10)")+
  scale_y_log10(guide = "axis_logticks")
ggOut
ggsave(filename = "//vf-lucid-r-o.lumcnet.prod.intern/lucid-r-o$/Projects/Senegal_phipseq/plotsCell/epitopeCount_histogram.pdf", plot = ggOut, 
       width =9, height = 11)

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


# #Word cloud the most frequent and rarest epitopes
# commonEpitopes <- table(peptideInfo[names(mostPrevalentEpitopes), "Pathogen"])
# rareEpitopes <- table(peptideInfo[names(mostRareEpitopes), "Pathogen"])
# 
# #Normalise for how many times an epitope of that organism is present in the dataset
# pathogenTotals <- table(peptideInfo$Pathogen)
# 
# #If a pathogen is only present 2 or fewer times, don't investigate it
# keepPathogen <- pathogenTotals[pathogenTotals > 2]
# 
# commonEpitopes <- commonEpitopes / pathogenTotals[names(commonEpitopes)]
# rareEpitopes <- rareEpitopes / pathogenTotals[names(rareEpitopes)]
# 
# commonEpitopes <- commonEpitopes[intersect(names(keepPathogen), names(commonEpitopes))]
# 
# wordcloud2(commonEpitopes, size = 0.2)
# 
# 
# wordcloud2(rareEpitopes, size = 0.8)


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

rowSums(normTotalEpitope_DF[,which(colnames(normTotalEpitope_DF) != "epitopeTotals")])

#Remove organisms which are only detected in 4 or fewer individuals
organismSum <- colSums(normTotalEpitope_DF > 0)
organismSum <- organismSum[organismSum >= 4]

normTotalEpitope_DF <- normTotalEpitope_DF[, names(organismSum)]

normTotalEpitope_DF$Residence <- metaData[row.names(normTotalEpitope_DF), "Residence"]
normTotalEpitope_DF$Residence <- factor(normTotalEpitope_DF$Residence, levels = c("Rural Senegalese", "Urban Senegalese", "Urban Dutch"))
normTotalEpitope_DF$Donor <- row.names(normTotalEpitope_DF)

test_normTotal <- normTotalEpitope_DF
test_normTotal <- test_normTotal[,str_detect(colnames(test_normTotal), "TotalNorm")]

rowSums(test_normTotal)

#### Look at the total number of epitopes and species found across the groups and individuals ####
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
ggsave(filename = "//vf-lucid-r-o.lumcnet.prod.intern/lucid-r-o$/Projects/Senegal_phipseq/plotsCell/nEpitopes_residence_comparison.pdf", plot = ggOut, 
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
ggsave(filename = "//vf-lucid-r-o.lumcnet.prod.intern/lucid-r-o$/Projects/Senegal_phipseq/plotsCell/n_sigAntiEpitopes.pdf", plot = ggOut, 
       width =4, height = 8)


#### Investigate rural, urban and dutch: PCA and cochran-armitage test (epitope presence/absence as input) ####

PCoA_phip <- pcoa(jaccard(t( as.matrix(phipSeq_p) )))

embedPCoord <- as.data.frame(PCoA_phip$vectors)
embedPCoord$Residence <- metaData[row.names(embedPCoord), "Residence"]
embedPCoord$sampleID <- row.names(embedPCoord)
embedPCoord$Sex <- metaData[row.names(embedPCoord), "Sex"]
embedPCoord$Age <- metaData[row.names(embedPCoord), "Age"]

embedPCoord$Residence <- factor(embedPCoord$Residence, levels = c("Rural Senegalese", "Urban Senegalese", "Urban Dutch"))

ggplot(embedPCoord, aes(x= Axis.1, y = Axis.2, color = Residence)) +
        geom_point()

#Boxplot of PCoord2 values split by Residence
PC2_lm <- lm(Axis.2 ~ Residence, data = embedPCoord)

PC2_emm <- emmeans(PC2_lm, "Residence")
trendOutput <- as.data.frame(contrast(PC2_emm, "poly"))

pairwiseOutput <- as.data.frame(pairs(PC2_emm))

group1 <- unlist(str_split(pairwiseOutput$contrast, " - "))[c(TRUE, FALSE)]
group2 <- unlist(str_split(pairwiseOutput$contrast, " - "))[c(FALSE, TRUE)]

ggpubr_formatted <- data.frame(".y."= "test", "group1" = group1, "group2" = group2,
                               p = pairwiseOutput$p.value)

dataMax <- max(embedPCoord$Axis.2)

ggpubr_formatted <- ggpubr_formatted %>%
  mutate(y.position = c( dataMax + (dataMax * 0.1) , dataMax + (dataMax * 0.2), dataMax + (dataMax * 0.3)))

ggpubr_formatted$p.symbol <- gtools::stars.pval(ggpubr_formatted$p)
ggpubr_formatted$p.symbol[ggpubr_formatted$p >= 0.05] <- "NS"

trendP <- subset(trendOutput, contrast == "linear")$p.value
if(trendP < 0.05){trendP <- "< 0.05"}else{trendP <- "\u2265 0.05"}

ggOut <- ggboxplot(embedPCoord, y = "Axis.2", x = "Residence", fill = "Residence", outlier.shape = NA) +
  stat_pvalue_manual(ggpubr_formatted, label = "p.symbol") +
  scale_fill_manual(values = area_palette) +
  ggtitle(label = paste0("Trend p-value ", trendP)) +
  ylab("MDS2") +
  xlab("") +
  geom_jitter() +
  grids(linetype = "solid")
ggOut
ggsave(filename = "//vf-lucid-r-o.lumcnet.prod.intern/lucid-r-o$/Projects/Senegal_phipseq/plotsCell/PC2_boxplot.pdf", plot = ggOut, 
       width =6, height = 8)

#Save PCoA data
write.csv(embedPCoord, "//vf-lucid-r-o.lumcnet.prod.intern/lucid-r-o$/Projects/Senegal_phipseq/dataCell/processed/PCoA_coords.csv")

PCoA_phip$vectors

PC1_varExplained <- (PCoA_phip[["values"]][["Relative_eig"]][1] / sum(PCoA_phip[["values"]][["Relative_eig"]])) * 100
PC2_varExplained <- (PCoA_phip[["values"]][["Relative_eig"]][2] / sum(PCoA_phip[["values"]][["Relative_eig"]])) * 100

##PCoA plot coloured by Area

permanovaOut <- adonis2(phipSeq_p ~ Residence, data = metaData, method = "jaccard")
permanovaOut <- permanovaOut$`Pr(>F)`[1]
if(permanovaOut < 0.05){permanovaOut <- "PERMANOVA < 0.05"} else{permanovaOut <- "PERMANOVA <= 0.05"}

gg <- merge(embedPCoord,aggregate(cbind(mean.Axis.1=Axis.1,mean.Axis.2=Axis.2)~Residence,embedPCoord,mean),by="Residence")
ggOut <- ggplot(gg, aes(Axis.1,Axis.2))+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0) +
  geom_point(size=5, aes(color=Residence))+
  geom_segment(aes(x=mean.Axis.1, y=mean.Axis.2, xend=Axis.1, yend=Axis.2, color = Residence)) +
  geom_point(aes(x=mean.Axis.1,y=mean.Axis.2, fill=Residence), size=7, shape = 21) +
  theme_minimal()+
  scale_fill_manual(values = area_palette) +
  scale_color_manual(values = area_palette)+
  xlab(paste0("MDS1 (", round(PC1_varExplained, digits = 2), "%)"))+
  ylab(paste0("MDS2 (", round(PC2_varExplained, digits = 2), "%)")) +
  ggtitle(permanovaOut)
ggOut
ggsave(filename = "//vf-lucid-r-o.lumcnet.prod.intern/lucid-r-o$/Projects/Senegal_phipseq/plotsCell/MDSPlot_filteredData_residenceColour.pdf", plot = ggOut, 
       width = 9, height = 8)

##PCoA plot coloured by Sex
permanovaOut <- adonis2(phipSeq_p ~ Sex, data = metaData, method = "jaccard")
permanovaOut <- permanovaOut$`Pr(>F)`[1]
if(permanovaOut < 0.05){permanovaOut <- "PERMANOVA < 0.05"} else{permanovaOut <- "PERMANOVA <= 0.05"}

#PERMANOVA for sex when dutch removed
adonis2(phipSeq_p[row.names(subset(metaData, Residence != "Urban Dutch")),] ~ Sex, data = subset(metaData, Residence != "Urban Dutch"), method = "jaccard")

gg <- merge(embedPCoord,aggregate(cbind(mean.Axis.1=Axis.1,mean.Axis.2=Axis.2)~Sex,embedPCoord,mean),by="Sex")
ggOut <- ggplot(gg, aes(Axis.1,Axis.2))+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0) +
  geom_point(size=5, aes(color=Sex))+
  geom_segment(aes(x=mean.Axis.1, y=mean.Axis.2, xend=Axis.1, yend=Axis.2, color = Sex)) +
  geom_point(aes(x=mean.Axis.1,y=mean.Axis.2, fill=Sex), size=7, shape = 21) +
  theme_minimal()+
  xlab(paste0("MDS1 (", round(PC1_varExplained, digits = 2), "%)"))+
  ylab(paste0("MDS2 (", round(PC2_varExplained, digits = 2), "%)")) +
  ggtitle(permanovaOut)
ggOut
ggsave(filename = "//vf-lucid-r-o.lumcnet.prod.intern/lucid-r-o$/Projects/Senegal_phipseq/plotsCell/MDSPlot_filteredData_sexColour.pdf", plot = ggOut, 
       width = 9, height = 8)

##PCoA plot coloured by Age
axis1Res <- cor.test(embedPCoord$Axis.1, embedPCoord$Age)
axis1Res <- if(axis1Res$p.value < 0.05){axis1Res <- paste0("MDS1 vs Age: p < 0.05 cor = ", axis1Res$estimate)}else{axis1Res <- paste0("MDS1 vs Age: p >= 0.05 cor = ", axis1Res$estimate)}

axis2Res <- cor.test(embedPCoord$Axis.2, embedPCoord$Age)
axis2Res <- if(axis2Res$p.value < 0.05){axis2Res <- paste0("MDS2 vs Age: p < 0.05 cor = ", axis2Res$estimate)}else{axis2Res <- paste0("MDS2 vs Age: p >= 0.05 cor = ", axis2Res$estimate)}

ggOut <- ggplot(embedPCoord, aes(Axis.1,Axis.2, color = Age))+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0) +
  geom_point(size=5) +
  theme_minimal() +
  xlab(paste0("MDS1 (", round(PC1_varExplained, digits = 2), "%)"))+
  ylab(paste0("MDS2 (", round(PC2_varExplained, digits = 2), "%)")) +
  scale_color_viridis_b() +
  ggtitle(paste0(axis1Res, "\n", axis2Res)) 
ggOut
ggsave(filename = "//vf-lucid-r-o.lumcnet.prod.intern/lucid-r-o$/Projects/Senegal_phipseq/plotsCell/MDSPlot_filteredData_ageColour.pdf", plot = ggOut, 
       width = 9, height = 8)


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
write.csv(testDF, "//vf-lucid-r-o.lumcnet.prod.intern/lucid-r-o$/Projects/Senegal_phipseq/dataCell//processed/cochranArmitage_urbanization.csv")

#Save amino acid sequence of sig epitopes
testDF$aa_seq <- peptideInfo[testDF$ID ,"aa_seq"]
write.csv(testDF, "//vf-lucid-r-o.lumcnet.prod.intern/lucid-r-o$/Projects/Senegal_phipseq/dataCell//processed/cochranArmitage_urbanization_aaSeq.csv")

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
ggsave(filename = "//vf-lucid-r-o.lumcnet.prod.intern/lucid-r-o$/Projects/Senegal_phipseq/plotsCell/commonProtein_cochranArmitage_medianTrendZ_barplot.pdf", plot = ggOut, 
       width = 16, height = 8)

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
ggsave(filename = "//vf-lucid-r-o.lumcnet.prod.intern/lucid-r-o$/Projects/Senegal_phipseq/plotsCell/cochranArmitage_increaseWithUrbanization_epitopePercentPresent.pdf", plot = ggOut, 
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
ggsave(filename = "//vf-lucid-r-o.lumcnet.prod.intern/lucid-r-o$/Projects/Senegal_phipseq/plotsCell/cochranArmitage_decreaseWithUrbanization_epitopePercentPresent.pdf", plot = ggOut, 
       width = 14, height = 12)

#### Carry out trend test on normalised pathogen totals ####

plotDF <- normTotalEpitope_DF

plotDF$Residence <- factor(plotDF$Residence, levels =  c("Rural Senegalese" , "Urban Senegalese" , "Urban Dutch" ))

#Remove _ from the colnames
colnames(plotDF) <- str_replace_all(colnames(plotDF), "_|\\:|\\.|\\(|\\)|\\=|\\,|\\||\\&", "")

plotDF <- plotDF %>%
                dplyr::select(matches("TotalNorm|Residence"))

plotList <- list()
pvalDF <- data.frame()
ggpubrFormatList <- list()

pathogenNames <- colnames(plotDF)
pathogenNames <- pathogenNames[-which(pathogenNames == "Residence")]

for(i in unique(pathogenNames)){
  
  featuresTest_lm <- lm(eval(parse(text = i)) ~ Residence, data = plotDF)
  
  featureTest_emm <- emmeans(featuresTest_lm, "Residence")
  trendOutput <- as.data.frame(contrast(featureTest_emm, "poly"))
  
  pairwiseOutput <- as.data.frame(pairs(featureTest_emm))
  
  group1 <- unlist(str_split(pairwiseOutput$contrast, " - "))[c(TRUE, FALSE)]
  group2 <- unlist(str_split(pairwiseOutput$contrast, " - "))[c(FALSE, TRUE)]
  
  ggpubr_formatted <- data.frame(".y."= "test", "group1" = group1, "group2" = group2,
                                 p = pairwiseOutput$p.value)
  
  dataMax <- max(plotDF[, i])
  
  ggpubr_formatted <- ggpubr_formatted %>%
    mutate(y.position = c( dataMax + (dataMax * 0.1) , dataMax + (dataMax * 0.2), dataMax + (dataMax * 0.3)))
  
  ggpubr_formatted$p.symbol <- gtools::stars.pval(ggpubr_formatted$p)
  ggpubr_formatted$p.symbol[ggpubr_formatted$p >= 0.05] <- "NS"
  
  trendP <- subset(trendOutput, contrast == "linear")$p.value
  trendEstimate <- subset(trendOutput, contrast == "linear")$estimate
  #if(trendP < 0.05){trendP <- "< 0.05"}else{trendP <- "\u2265 0.05"}
  
  pvalDF <- rbind(pvalDF, data.frame(measurement = i,
                                     pval = trendP,
                                     estimate = trendEstimate))
  
  ggpubrFormatList[[i]] <- ggpubr_formatted
  

}

pvalDF$p_adj <- p.adjust(pvalDF$pval, method = "BH")

ggOut <- ggboxplot(plotDF, y = "humangammaherpesvirus8TotalNorm", x = "Residence", fill = "Residence", outlier.shape = NA) + 
  scale_fill_manual(values = area_palette) +
  xlab("") +
  geom_jitter() +
  grids(linetype = "solid") +
  coord_flip() + 
  guides(fill="none")
ggOut

#Plot the significant pathogens
featuresTest <- normTotalEpitope_DF

#Remove _ from the colnames
colnames(featuresTest) <- str_replace_all(colnames(featuresTest), "_|\\:|\\.|\\(|\\)|\\=|\\,|\\||\\&", "")

#Seperate the data into the different trend directions, then take the top 4 in terms of adjusted p-value
posDirection <- subset(pvalDF, estimate > 0)
negDirection <- subset(pvalDF, estimate < 0)

top4Pos <- posDirection$p_adj
names(top4Pos) <- posDirection$measurement
top4Pos <- names(top4Pos[order(top4Pos, decreasing = FALSE)])[1:4]

top4Neg <- negDirection$p_adj
names(top4Neg) <- negDirection$measurement
top4Neg <- names(top4Neg[order(top4Neg, decreasing = FALSE)])[1:4]

featuresTest <- featuresTest[,c(top4Pos, top4Neg, "Residence")]


for(i in c(top4Pos, top4Neg) ){
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
  
  trendP <- pvalDF[pvalDF$measurement == i,"p_adj"]
  #if(trendP < 0.05){trendP <- "< 0.05"}else{trendP <- "\u2265 0.05"}
  trendP <- signif(trendP, digits = 3)
  
  ggOut <- ggboxplot(featuresTest, y = i, x = "Residence", fill = "Residence", outlier.shape = NA) +
    stat_pvalue_manual(ggpubr_formatted, label = "p.symbol") +
    scale_fill_manual(values = area_palette) +
    ggtitle(label = paste0("Trend p-value ", trendP)) +
    ylab(plotName) +
    xlab("") +
    geom_jitter() +
    grids(linetype = "solid") +
    coord_flip() + 
    guides(fill="none")
  
  ggsave(filename = paste0("//vf-lucid-r-o.lumcnet.prod.intern/lucid-r-o$/Projects/Senegal_phipseq/plotsCell/normPathogenPlots_figure/",i,"_enrichment_urbanization_normPathogen.pdf"), plot = ggOut, 
         width = 6, height = 4)
  
  plotList[[i]] <- ggOut
}


#Calculate significance
featuresTest <- normTotalEpitope_DF
unique(featuresTest$Residence)

#Remove _ from the colnames
colnames(featuresTest) <- str_replace_all(colnames(featuresTest), "_|\\:|\\.|\\(|\\)|\\=|\\,|\\||\\&", "")

plotList <- list()


for(i in c(pvalDF$measurement) ){
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
  
  trendP <- pvalDF[pvalDF$measurement == i,"p_adj"]
  #if(trendP < 0.05){trendP <- "< 0.05"}else{trendP <- "\u2265 0.05"}
  trendP <- signif(trendP, digits = 3)
  
  ggOut <- ggboxplot(featuresTest, y = i, x = "Residence", fill = "Residence", outlier.shape = NA) +
    stat_pvalue_manual(ggpubr_formatted, label = "p.symbol") +
    scale_fill_manual(values = area_palette) +
    ggtitle(label = paste0("Trend p-value ", trendP)) +
    ylab(plotName) +
    xlab("") +
    geom_jitter() +
    grids(linetype = "solid") +
    coord_flip() + 
    guides(fill="none")
  
  ggsave(filename = paste0("//vf-lucid-r-o.lumcnet.prod.intern/lucid-r-o$/Projects/Senegal_phipseq/plotsCell/normPathogenPlots/",i,"_enrichment_urbanization_normPathogen.pdf"), plot = ggOut, 
         width = 6, height = 4)
  
  plotList[[i]] <- ggOut
}

#### Calculate diversity metrics at the organism level ####

#The results seem to show
#The richness is not different between the groups, nor is shannon diversity
#PERMANOVA shows differences
#Betadispersion isn't different, but the senegalese ones have greater variability for distance from centroid

library("ape")
library("vegan")
library("dplyr")
library("hagis")
library("ggplot2")

data(BCI)
data(dune)
data(dune.env)
data(sipoo)
data(varespec)


#We calculate diversity not on the epitopes, but on the organisms

totalEpitope_DF <- data.frame(epitopeTotals = rowSums(phipSeq_p_raw[metaData$ID,]),
                                  row.names = metaData$ID)

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
  
  totalEpitope_DF[[paste0(i, "TotalNorm")]] <- currentEpitopesTotals[row.names(totalEpitope_DF)] 
  
}

allData <- totalEpitope_DF

allData$Residence <- metaData[row.names(allData),"Residence"]

urbanPhipSeq <- subset(allData, Residence == "Urban Senegalese")
ruralPhipSeq <- subset(allData, Residence == "Rural Senegalese")
dutchPhipSeq <- subset(allData, Residence == "Urban Dutch")

urbanPhipSeq <- urbanPhipSeq %>% dplyr::select(!c(Residence, epitopeTotals))
ruralPhipSeq <- ruralPhipSeq %>% dplyr::select(!c(Residence, epitopeTotals))
dutchPhipSeq <- dutchPhipSeq %>% dplyr::select(!c(Residence, epitopeTotals))

#New section
organismCounts <- rbind(rbind(ruralPhipSeq, urbanPhipSeq), dutchPhipSeq)
organismCounts_binary <- as.matrix(organismCounts/organismCounts)
organismCounts_binary[is.na(organismCounts_binary)] <- 0
epitopeCounts <- phipSeq_p_raw_df

sppr <- specnumber(organismCounts_binary)

#Species richness is not different across the sites
sppr_aov <- aov(sppr ~ Residence, data = metaData[row.names(organismCounts_binary),])
summary(sppr_aov)

shannonOut <- diversity(organismCounts_binary)
shannonPlot <- data.frame("shannon" = shannonOut,
                          "Residence" = factor(metaData[names(shannonOut),"Residence"], levels = c("Rural Senegalese", "Urban Senegalese", "Urban Dutch")))

shannon_lm <- lm(shannon ~ Residence, data = shannonPlot)

shannon_emm <- emmeans(shannon_lm, "Residence")
trendOutput <- as.data.frame(contrast(shannon_emm, "poly"))

pairwiseOutput <- as.data.frame(pairs(shannon_emm))

group1 <- unlist(str_split(pairwiseOutput$contrast, " - "))[c(TRUE, FALSE)]
group2 <- unlist(str_split(pairwiseOutput$contrast, " - "))[c(FALSE, TRUE)]

ggpubr_formatted <- data.frame(".y."= "test", "group1" = group1, "group2" = group2,
                               p = pairwiseOutput$p.value)

dataMax <- max(shannonPlot$shannon)

ggpubr_formatted <- ggpubr_formatted %>%
  mutate(y.position = c( dataMax + (dataMax * 0.1) , dataMax + (dataMax * 0.2), dataMax + (dataMax * 0.3)))

ggpubr_formatted$p.symbol <- gtools::stars.pval(ggpubr_formatted$p)
ggpubr_formatted$p.symbol[ggpubr_formatted$p >= 0.05] <- "NS"

trendP <- subset(trendOutput, contrast == "linear")$p.value
if(trendP < 0.05){trendP <- "< 0.05"}else{trendP <- "\u2265 0.05"}

ggOut <- ggboxplot(shannonPlot, y = "shannon", x = "Residence", fill = "Residence", outlier.shape = NA) +
  stat_pvalue_manual(ggpubr_formatted, label = "p.symbol") +
  scale_fill_manual(values = area_palette) +
  ggtitle(label = paste0("Trend p-value ", trendP)) +
  ylab("Shannon entropy") +
  xlab("") +
  geom_jitter() +
  grids(linetype = "solid")
ggOut

bird_perm <- adonis2(organismCounts ~ Residence, data = metaData[row.names(organismCounts),], method = "bray")
bird_perm

x <- betadisper(vegdist(organismCounts, method = "bray"), group = metaData$Residence, type = "median")
anova(x)
TukeyHSD(x)
boxplot(x)

#Because of the class imbalance, bootstrap 8 individuals per group and then perform dispersion tests
selectedN <- 8

distMedianAvg <- data.frame("sum" = rep(0, nrow(metaData)), 
                            "i" = rep(0, nrow(metaData)),
                            "Residence" = metaData$Residence,
                            row.names = row.names(metaData))

for(i in 1:50){
  
  urbDutch_sample <- subset(metaData, Residence == "Urban Dutch")
  urbDutch_sample <- row.names(urbDutch_sample[sample(row.names(urbDutch_sample), size = selectedN, replace = FALSE),])
  
  urbSenegal_sample <- subset(metaData, Residence == "Urban Senegalese")
  urbSenegal_sample <- row.names(urbSenegal_sample[sample(row.names(urbSenegal_sample), size = selectedN, replace = FALSE),])
  
  rurSenegal_sample <- subset(metaData, Residence == "Rural Senegalese")
  rurSenegal_sample <- row.names(rurSenegal_sample[sample(row.names(rurSenegal_sample), size = selectedN, replace = FALSE),])
  
  selectedSamples <- c(urbDutch_sample, urbSenegal_sample, rurSenegal_sample)
  
  currentData <- organismCounts[selectedSamples,]
  currentMetaData <- metaData[selectedSamples,]
  
  betaDisp_res <- betadisper(vegdist(currentData, method = "bray"), group = currentMetaData$Residence, type = "median")
  
  distMedianAvg[selectedSamples,"i"]  <- distMedianAvg[selectedSamples,"i"] + 1
  distMedianAvg[selectedSamples,"sum"] <- distMedianAvg[selectedSamples,"sum"] + betaDisp_res$distances[selectedSamples]
  
}

distMedianAvg$avg <- distMedianAvg$sum / distMedianAvg$i

ggplot(distMedianAvg, aes(x= Residence, y = avg)) +
        geom_boxplot()

mod3 <- betadisper(vegdist(organismCounts, method = "bray"), metaData$Residence, type = "centroid")
mod3
permutest(mod3, permutations = 99)
anova(mod3)
plot(mod3)
boxplot(mod3)
plot(TukeyHSD(mod3))

#Jaccard on epitope presence
epitopeCounts.jaccard <-
  vegdist(epitopeCounts, "jaccard", na.rm = TRUE)
princoor.epitopeCounts <- pcoa(epitopeCounts.jaccard) 

princoor.epitopeCounts.data <-
  data.frame(
    Sample = rownames(princoor.epitopeCounts$vectors),
    X = princoor.epitopeCounts$vectors[, 1],
    Y = princoor.epitopeCounts$vectors[, 2]
  )

princoor.epitopeCounts.data$Residence <- metaData[princoor.epitopeCounts.data$Sample,"Residence"]

ggplot(data = princoor.epitopeCounts.data, aes(x = X, y = Y)) +
  geom_point(aes(colour = Residence))  +
  theme_bw() +
  theme(
    axis.title.x = element_text(face = "bold", size = 15),
    axis.title.y = element_text(face = "bold", size = 15),
    axis.text = element_text(face = "bold", size = 10),
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(face = "bold", size = 10),
    legend.key.size = unit(1, 'lines')
  ) +
  stat_ellipse(data = princoor.epitopeCounts.data, aes(x = X, y = Y),
               level = 0.95) +
  ggtitle("Pathotype Jaccard Distances PCOA") +
  scale_color_manual(values = area_palette)

epitope.disp <-
  betadisper(epitopeCounts.jaccard, princoor.epitopeCounts.data$Residence)
anova(epitope.disp)
TukeyHSD(epitope.disp)
plot(epitope.disp, hull = FALSE, ellipse = TRUE)

adonisOut <- adonis2(epitopeCounts.jaccard ~  princoor.epitopeCounts.data$Residence)
adonisOut

pairwise_permanova(sp_matrix = epitopeCounts, group_var = embedPCoord$Residence, dist = "jaccard")

#Jaccard on organism presence
organismCounts.jaccard <-
  vegdist(organismCounts_binary, "jaccard", na.rm = TRUE)
princoor.organismCounts <- pcoa(organismCounts.jaccard) 

princoor.organismCounts.data <-
  data.frame(
    Sample = rownames(princoor.organismCounts$vectors),
    X = princoor.organismCounts$vectors[, 1],
    Y = princoor.organismCounts$vectors[, 2]
  )

princoor.organismCounts.data$Residence <- metaData[princoor.organismCounts.data$Sample,"Residence"]

ggplot(data = princoor.organismCounts.data, aes(x = X, y = Y)) +
  geom_point(aes(colour = Residence))  +
  theme_bw() +
  theme(
    axis.title.x = element_text(face = "bold", size = 15),
    axis.title.y = element_text(face = "bold", size = 15),
    axis.text = element_text(face = "bold", size = 10),
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(face = "bold", size = 10),
    legend.key.size = unit(1, 'lines')
  ) +
  stat_ellipse(data = princoor.organismCounts.data, aes(x = X, y = Y),
               level = 0.95) +
  ggtitle("Pathotype Jaccard Distances PCOA") +
  scale_color_manual(values = area_palette)

organism.disp <-
  betadisper(organismCounts.jaccard, princoor.organismCounts.data$Residence)
anova(organism.disp)
TukeyHSD(organism.disp)
plot(organism.disp, hull = FALSE, ellipse = TRUE)

adonisOut <- adonis2(organismCounts.jaccard ~  princoor.organismCounts.data$Residence)
adonisOut

head(BCI)

epitopeDiv <- diversity(epitopeCounts)
organismDiv <- diversity(organismCounts)

epitope_rare <- rarefy(epitopeCounts, min(rowSums(epitopeCounts)))
organism_rare <- rarefy(organismCounts, min(rowSums(organismCounts)))

organismRare_DF <- data.frame("rare_rich" = organism_rare,
                              "Residence" = factor(metaData[names(organism_rare),"Residence"], levels = c("Rural Senegalese", "Urban Senegalese", "Urban Dutch") ))

stat.test <- organismRare_DF %>%
                wilcox_test(rare_rich ~ Residence) %>%
                add_significance() %>%
                adjust_pvalue(method = "BH") %>%
                add_xy_position("Residence")

ggplot(organismRare_DF, aes(x = Residence, y = rare_rich, fill = Residence)) +
          geom_boxplot(outlier.shape = NA) +
          geom_point(position = position_jitterdodge())+
          scale_fill_manual(values = area_palette) +
          stat_compare_means()

z <- betadiver(organismCounts, "z")
mod <- with(metaData, betadisper(z, Residence))
mod

outputDF <- rbind(rbind(ruralPhipSeq, urbanPhipSeq), dutchPhipSeq)

simpsonRes <- diversity(outputDF,index = "simpson")

simpsonPlot <- data.frame(simpsonRes)
simpsonPlot$Residence <- metaData[row.names(simpsonPlot), "Residence"]

simpsonPlot$Residence <- factor(simpsonPlot$Residence, levels = c("Rural Senegalese", "Urban Senegalese", "Urban Dutch"))

simpsonTest <- simpsonPlot %>%
  wilcox_test(simpsonRes ~ Residence) %>%
  add_significance() %>%
  adjust_pvalue(method = "BH") %>%
  add_xy_position("Residence")

ggOut <- ggboxplot(simpsonPlot, x= "Residence", y = "simpsonRes", fill = "Residence", outlier.shape = NA)  +
  geom_point(position = position_jitterdodge(), aes(fill = Residence)) +
  stat_pvalue_manual(simpsonTest)+
  scale_fill_manual(values= area_palette) +
  stat_compare_means()+
  ylab("Simpson's Index") +
  xlab("")

ggOut
ggsave(filename = paste0("//vf-lucid-r-o.lumcnet.prod.intern/lucid-r-o$/Projects/Senegal_phipseq/plotsCell/simpsonDiversity_organism.pdf"), plot = ggOut)

shannonRes <- diversity(outputDF,index = "shannon")

shannonPlot <- data.frame(shannonRes)
shannonPlot$Residence <- metaData[row.names(shannonPlot), "Residence"]

shannonPlot$Residence <- factor(shannonPlot$Residence, levels = c("Rural Senegalese", "Urban Senegalese", "Urban Dutch"))

shannonTest <- shannonPlot %>%
                wilcox_test(shannonRes ~ Residence) %>%
                add_significance() %>%
                adjust_pvalue(method = "BH") %>%
                add_xy_position("Residence")

ggOut <- ggboxplot(shannonPlot, x= "Residence", y = "shannonRes", fill = "Residence", outlier.shape = NA)  +
  geom_point(position = position_jitterdodge(), aes(fill = Residence)) +
  stat_pvalue_manual(shannonTest)+
  scale_fill_manual(values= area_palette) +
  stat_compare_means()+
  ylab("Simpson's Index") +
  xlab("")

ggsave(filename = paste0("//vf-lucid-r-o.lumcnet.prod.intern/lucid-r-o$/Projects/Senegal_phipseq/plotsCell/shannonDiversity_organism.pdf"), plot = ggOut)

ggplot(simpsonPlot, aes(x = simpsonRes, y = Residence, fill = Residence)) + 
    geom_density_ridges() +
    geom_point() +
    scale_fill_manual(values = area_palette) +
    theme_classic() +
    xlab("Simpson's Index") +
    ylab("")

#### Representative urban, rural and dutch proportion barplot ####

phipSeq_p_residence_raw_df <- phipSeq_p_raw_df

phipSeq_p_residence_raw_df$Residence <- metaData[row.names(phipSeq_p_residence_raw_df), "Residence"]

phipSeq_p_residence_raw_df <- phipSeq_p_residence_raw_df %>%
                                group_by(Residence) %>%
                                summarise(across(everything(), ~ sum(.x)))
phipSeq_p_residence_raw_df <- as.data.frame(phipSeq_p_residence_raw_df)
row.names(phipSeq_p_residence_raw_df) <- phipSeq_p_residence_raw_df$Residence
phipSeq_p_residence_raw_df$Residence <- NULL

normTotalEpitope_residence_DF <- data.frame(epitopeTotals = rowSums(phipSeq_p_residence_raw_df),
                                            row.names = row.names(phipSeq_p_residence_raw_df))


pathogens <- unique(peptideInfo$Pathogen)

for(i in pathogens){
  
  peptideInfoCurrent <- peptideInfo[peptideInfo$Pathogen == i,]
  
  currentEpitopes <- row.names(peptideInfoCurrent)
  
  #Add current epitopes bound to data frame
  if(length(currentEpitopes)==1){
    currentEpitopesTotals <- phipSeq_p_residence_raw_df[,currentEpitopes]
    
  }
  
  
  else{
    currentEpitopesTotals <- rowSums(phipSeq_p_residence_raw_df[,currentEpitopes])
    
  }
  
  #Sort name 
  i <- str_replace_all(i, "-| |\\/", "")
  
  normTotalEpitope_residence_DF[[paste0(i, "TotalNorm")]] <- (currentEpitopesTotals[row.names(normTotalEpitope_residence_DF)] / normTotalEpitope_residence_DF$epitopeTotals) * 100
  
}




#### Compare correlations in epitope-ome proportions for senegalese ####

senegaleseData <- normTotalEpitope_DF[str_detect(row.names(normTotalEpitope_DF), "LD", negate = TRUE),]
senegaleseData <- senegaleseData[,-which(colnames(senegaleseData) == "epitopeTotals")]
senegaleseData$MDS2 <- embedPCoord[row.names(senegaleseData), "Axis.2"]

corSenegal <- cor( as.matrix( senegaleseData[,-which(colnames(senegaleseData) %in% c("Residence", "Donor") )] ))

flexneri_negCor <- names(sort(corSenegal[, str_detect(colnames(corSenegal), "flexneri")], decreasing = FALSE)[1:10])
#sort(corSenegal[, str_detect(colnames(corSenegal), "flexneri")], decreasing = TRUE)[1:10]

corRural <- cor( as.matrix( senegaleseData[ senegaleseData$Residence == "Rural Senegalese" ,-which(colnames(senegaleseData) %in% c("Residence", "Donor") )] ))
corUrban <- cor( as.matrix( senegaleseData[ senegaleseData$Residence == "Urban Senegalese" ,-which(colnames(senegaleseData) %in% c("Residence", "Donor") )] ))

corRural <- corRural["shigellaflexneriTotalNorm",flexneri_negCor]
corUrban <- corUrban["shigellaflexneriTotalNorm",flexneri_negCor]

plotCor <- rbind(corRural, corUrban)

ComplexHeatmap::Heatmap(plotCor, cluster_rows = FALSE, cluster_columns = FALSE)

ggplot(senegaleseData, aes(x = shigellaflexneriTotalNorm, y = `bacteroidesfragilisTotalNorm`, color = Residence)) +
  geom_point() +
  stat_cor()

ggplot(senegaleseData, aes(x = shigellaflexneriTotalNorm, y = `agathobacterrectalisTotalNorm`, color = Residence)) +
  geom_point()+
  stat_cor()

ggplot(senegaleseData, aes(x = shigellaflexneriTotalNorm, y = `phocaeicolavulgatusTotalNorm`, color = Residence)) +
  geom_point()+
  stat_cor()

ggplot(senegaleseData, aes(x = shigellaflexneriTotalNorm, y = `dorealongicatenaTotalNorm`, color = Residence)) +
  geom_point()+
  stat_cor()

ggplot(senegaleseData, aes(x = shigellaflexneriTotalNorm, y = shigellasonneiTotalNorm, color = Residence)) +
  geom_point()+
  stat_cor(method = "spearman")

ggplot(senegaleseData, aes(x = shigellaflexneriTotalNorm, y = coxsackievirusTotalNorm, color = Residence)) +
  geom_point()+
  stat_cor(method = "spearman")

#### Heatmap of presence/absence of significant epitopes ####
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
col_fun <- colorRamp2(c(0, 1), c("#afd8e5", "#ee6550"))

column_ha <- ComplexHeatmap::HeatmapAnnotation(Residence = metaData[colnames(heatmapMtx), "Residence"],
                               col = list(Residence = c("Rural Senegalese" = "#e59f01", "Urban Senegalese" = "#54b3e8", "Urban Dutch" = "#009e72")))

#Get new row names: [proteinName] ([organismName])
row.names(heatmapMtx) <- paste0(row.names(peptideInfo[row.names(heatmapMtx),]),": ",firstup(peptideInfo[row.names(heatmapMtx),"Protein"]), " (", firstup(peptideInfo[row.names(heatmapMtx),"Pathogen"]) ,")")

dendrogramRes <- {set.seed(42); hclust(dist(heatmapMtx))}
dendrogramCut <- {set.seed(42); cutree(dendrogramRes, k = 4)}
dendrogramCut <- factor(dendrogramCut, levels = c(4,2,1,3))

#Get selected row names
selectedRowNames <- c("twist_58790: Genome polyprotein  (Rhinovirus b)", "agilent_12017: Envelope glycoprotein d  (Human herpesvirus 1)", 
                      "agilent_141247: Invasin ipac (Shigella sonnei)",
                      "agilent_230549: Invasin ipac (Shigella flexneri)", "agilent_219967: Grab (Streptococcus pyogenes)",
                      "twist_41031: K8.1 (Human gammaherpesvirus 8)",
                      "twist_54039: Interspersed repeat antigen (Plasmodium falciparum)", "agilent_237918: Adhesin (Haemophilus influenzae)",
                      "agilent_8991: Epstein-barr nuclear antigen 6 (Epstein-barr virus)", "agilent_2709: Syncytin-1 (Homo sapiens)",
                      "agilent_223729: Iga-specific serine endopeptidase autotransporter (Neisseria meningitidis)")
intersect(row.names(heatmapMtx), selectedRowNames)
setdiff(selectedRowNames, row.names(heatmapMtx))

epitopeIndex_DF <- data.frame("index" = seq(1, nrow(heatmapMtx), 1),
                              row.names = row.names(heatmapMtx))

########### Add Z scores to the left side of the heatmap, that means I can get rid of the barplot
testDF

rowAnnotation_data <- testDF
row.names(rowAnnotation_data) <- rowAnnotation_data$ID

rowAnnotation_data <- rowAnnotation_data[str_split_fixed(row.names(heatmapMtx), pattern = "\\: ", n = 2 )[,1],]

#This code plots barplots of Z score
ha <- ComplexHeatmap::rowAnnotation(numeric = ComplexHeatmap::anno_numeric(rowAnnotation_data$Z,  
                                          bg_gp = gpar(fill = c("#069e73", "#e5a124"))), 
                   annotation_name_rot = 0)

#This code plots colours 
ha <- ComplexHeatmap::rowAnnotation(bar2 = ComplexHeatmap::anno_barplot(rowAnnotation_data$Z * -1))

epitopeAnnotate = ComplexHeatmap::rowAnnotation(foo = ComplexHeatmap::anno_mark(at = epitopeIndex_DF[selectedRowNames,"index"], 
                                   labels = selectedRowNames))

pdf(paste0("//vf-lucid-r-o.lumcnet.prod.intern/lucid-r-o$/Projects/Senegal_phipseq/plotsCell/peptideSig_heatmap_selected.pdf"), height = 18, width = 26)
heatmapPlot <- ComplexHeatmap::Heatmap(heatmapMtx, col = col_fun, 
                        cluster_rows = TRUE,  row_split = dendrogramCut, cluster_row_slices = FALSE, right_annotation = epitopeAnnotate, left_annotation = ha, 
                        cluster_columns = FALSE, top_annotation = column_ha,
                        show_row_dend = FALSE)
ComplexHeatmap::draw(heatmapPlot)
dev.off()

pdf(paste0("//vf-lucid-r-o.lumcnet.prod.intern/lucid-r-o$/Projects/Senegal_phipseq/plotsCell/peptideSig_heatmap.pdf"), height = 18, width = 26)
heatmapPlot <- ComplexHeatmap::Heatmap(heatmapMtx, col = col_fun, 
                                       cluster_rows = TRUE,  row_split = dendrogramCut, cluster_row_slices = FALSE, left_annotation = ha,
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
          file = "//vf-lucid-r-o.lumcnet.prod.intern/lucid-r-o$/Projects/Senegal_phipseq/dataCell//cochranArmitage_sig_epitopes.csv")

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
          file = "//vf-lucid-r-o.lumcnet.prod.intern/lucid-r-o$/Projects/Senegal_phipseq/dataCell//epitope_information.csv")

#### Compare different organism results ####
ruralNormTotal <- subset(normTotalEpitope_DF, Residence == "Rural Senegalese")
ruralNormTotal <- subset(normTotalEpitope_DF, Residence != "Urban Dutch")

head(ruralNormTotal)

#Remove non-useful Info from ruralNormTotal
ruralNormTotal <- ruralNormTotal[-which(colnames(ruralNormTotal) %in% c("epitopeTotals", "Residence", "Donor", "Sex"))]

#Remove columns which have a variance of 0
ruralNormTotal <- ruralNormTotal[, colVars(ruralNormTotal) != 0]

#Remove columns where the pathogens are not in the selected epitope list
chosenPathogens <- paste(str_replace_all(subsetPeptideInfo$Pathogen, "-| ", ""), collapse = "|")
ruralNormTotal <- ruralNormTotal[,str_detect(colnames(ruralNormTotal), chosenPathogens)]

combData <- ruralNormTotal

head(combData)

combData$MDS2 <- embedPCoord[row.names(combData),"Axis.2"]

combData <- as.matrix(combData)

#Correlation with shigella flexneri 
corRes <- cor(combData[,str_detect(colnames(combData), "TotalNorm")])

sort(corRes[, str_detect(colnames(corRes), "shigellaflexneriTotalNorm")], decreasing = FALSE)[1:5]
sort(corRes[, str_detect(colnames(corRes), "shigellaflexneriTotalNorm")], decreasing = TRUE)[1:6]

ggplot(combData, aes(x = shigellaflexneriTotalNorm, y = `phocaeicolavulgatusTotalNorm`)) +
  geom_point()

