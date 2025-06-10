#GSVA on UAMS
library(readr)
library(tidyr)
library(dplyr)
library(janitor)
library(openxlsx)
library(GSVA)
library(biomaRt) 
library(rio)
UAMS <- read_csv("./UAMS_expression.csv")
patients <- UAMS["...1"]
gene_set <- read.xlsx("./DGE\ common\ genes.xlsx",colNames = FALSE)
gene_set$up_gene <- gsub("\\..*","",gene_set$X1)
gene_set$down_gene <- gsub("\\..*","",gene_set$X2)

ensembl = useEnsembl("ensembl",dataset="hsapiens_gene_ensembl",mirror="www") 
affy_ensembl= c("affy_hg_u133_plus_2", "ensembl_gene_id","hgnc_symbol")
dict=getBM(attributes=affy_ensembl, filters="hgnc_symbol",values = c(gene_set$up_gene,gene_set$down_gene), mart = ensembl)
dict <- subset(dict, affy_hg_u133_plus_2 %in% colnames(UAMS))

up_symbol_matrix <- matrix(0, nrow = nrow(UAMS), ncol = nrow(gene_set)-2)
for (i in 1:(nrow(gene_set)-2)){
  gene <- gene_set$up_gene[i]
  sub_dict <- subset(dict, hgnc_symbol == gene)
  sub_expression <- UAMS[names(UAMS) %in% sub_dict$affy_hg_u133_plus_2]
  
  up_symbol_matrix[,i] <- rowMeans(sub_expression)
}

colnames(up_symbol_matrix) <- na.omit(gene_set$up_gene)

down_symbol_matrix <- matrix(0, nrow = nrow(UAMS), ncol = nrow(gene_set))
for (i in 1:nrow(gene_set)){
  gene <- gene_set$down_gene[i]
  sub_dict <- subset(dict, hgnc_symbol == gene)
  sub_expression <- UAMS[names(UAMS) %in% sub_dict$affy_hg_u133_plus_2]
  
  down_symbol_matrix[,i] <- rowMeans(sub_expression)
}

colnames(down_symbol_matrix) <- gene_set$down_gene

expression_matrix <- cbind(up_symbol_matrix,down_symbol_matrix)

GSVA::geneSets(gene_set)
gbmPar <- gsvaParam(t(expression_matrix), list(gene_set[1:238,3],gene_set[,4]), maxDiff=FALSE)

gsva_results <- gsva(
  gbmPar
)

colnames(gsva_results) <- UAMS$...1
rownames(gsva_results) <- c("Up-regulated", "Down-regulated")

write.csv(t(gsva_results), "./gsva_result.csv")
