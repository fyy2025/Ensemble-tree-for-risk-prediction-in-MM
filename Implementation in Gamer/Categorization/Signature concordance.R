library(readr)
setwd("/Users/fanyiyang/Desktop/Research/summer/Ensemble-tree-for-risk-prediction-in-MM/Implementation\ in\ Gamer/Categorization")

# filter for high risk group
Gamer <- read.csv("./merged_Gamer_scores_with_patients.csv")
GamerHR <- subset(Gamer, group=="Group1")
binary <- GamerHR
for (i in 2:(ncol(binary)-5)){
  threshold = quantile(binary[,i],0.8)
  binary[,i] = ifelse(binary[,i]>threshold,1,0)
}

GPI_threshold = quantile(GamerHR[,4],0.92)
binary[,4] = ifelse(GamerHR[,4]>GPI_threshold,1,0)

EMC <- subset(binary, EMC==1)$patient
UAMS80 <- subset(binary, UAMS80==1)$patient
UAMS70 <- subset(binary, UAMS70==1)$patient
UAMS17 <- subset(binary, UAMS17==1)$patient
HM19 <- subset(binary, HM19==1)$patient
IFM15 <- subset(binary, IFM15==1)$patient
GPI <- subset(binary, GPI==1)$patient


all_patients <- unique(c(EMC, UAMS70, UAMS17, UAMS80, IFM15, HM19, GPI))

# Create a binary membership matrix
membership_matrix <- data.frame(
  EMC = as.integer(all_patients %in% EMC),
  UAMS17 = as.integer(all_patients %in% UAMS17),
  UAMS70 = as.integer(all_patients %in% UAMS70),
  UAMS80 = as.integer(all_patients %in% UAMS80),
  IFM15 = as.integer(all_patients %in% IFM15),
  HM19 = as.integer(all_patients %in% HM19),
  GPI = as.integer(all_patients %in% GPI)
)

rownames(membership_matrix)=all_patients

library(UpSetR)
library(grid)
upset(
  membership_matrix,
  sets = c("EMC", "UAMS70", "UAMS17", "UAMS80", "IFM15", "HM19", "GPI"),
  order.by = "freq"
)

grid.text("UpSet Plot for high risk patients in GAMER HR",x = 0.65, y=0.95, gp=gpar(fontsize=15))
