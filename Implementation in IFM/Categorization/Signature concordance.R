library(readr)
setwd("/Users/fanyiyang/Desktop/Research/summer/Ensemble-tree-for-risk-prediction-in-MM/Implementation\ in\ IFM/Categorization")
EMC <- read.table("./IFM-EMC92.txt")$V1
UAMS70 <- read.table("./IFM-UAMS70.txt")$V1
UAMS17 <- read.table("./IFM-UAMS17.txt")$V1
UAMS80 <- read.table("./IFM-UAMS80.txt")$V1
IFM15 <- read.table("./IFM-IFM15.txt")$V1
HM19 <- read.table("./IFM-HM19.txt")$V1
GPI <- read.table("./IFM-GPI.txt")$V1

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

grid.text("UpSet Plot for high risk patients in IFM",x = 0.65, y=0.95, gp=gpar(fontsize=15))
