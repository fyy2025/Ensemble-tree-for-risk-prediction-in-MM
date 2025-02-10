library(geneClassifiers)
library(biomaRt) 
library(readr)
library(tidyr)
library(dplyr)
library(janitor)

ensembl <- useEnsembl("ensembl",dataset="hsapiens_gene_ensembl",mirror="www") 
affy_ensembl <- c("affy_hg_u133_plus_2", "ensembl_gene_id")
dict <- getBM(attributes=affy_ensembl, filters="affy_hg_u133_plus_2",values = probeset_id, mart = ensembl)

MMRF_TPM <- read_tsv("./MMRF_TPM.tsv")
TPM_first_visit <- MMRF_TPM %>% dplyr::select("Gene",ends_with("_1_BM_CD138pos"))
TPM_first_EI  <-  subset(TPM_first_visit,Gene %in% dict$ensembl_gene_id)
tidydata <- t(TPM_first_EI)%>%
  row_to_names(row_number=1)
patient <- row.names(tidydata)
tidydata <- as.data.frame(tidydata)
tidydata <- as.data.frame(sapply(tidydata,as.numeric ))
scaled_data <- scale(tidydata)

HGU133Plus2_Hs_ENSG_mapping <- read_delim("./HGU133Plus2_Hs_ENSG_mapping.txt",delim = "\t", escape_double = FALSE,  trim_ws = TRUE)
GPI_gene_list <- read.csv("./Proliferation\ Index\ Gene\ Symbols.csv")$Gene.Symbol
GPI_dict <- subset(annotLookup,hgnc_symbol %in% gene_list)

gene_signature_score = function(signature,scaled_data){
  if (!(signature %in% c("GPI,EI"))){
    Classifier <- getClassifier(signature)
    weights <- getWeights(Classifier)
    probeset_id <- names(weights)
    
    HGU133Plus2_Hs_ENSG_mapping2 <- HGU133Plus2_Hs_ENSG_mapping[HGU133Plus2_Hs_ENSG_mapping$`Affy Probe Set Name` %in% probeset_id,c(1,7)]
    HGU133Plus2_Hs_ENSG_mapping2 <- HGU133Plus2_Hs_ENSG_mapping2[!duplicated(HGU133Plus2_Hs_ENSG_mapping2),]
    brainarray_out <- setdiff(probeset_id,HGU133Plus2_Hs_ENSG_mapping2$`Affy Probe Set Name`) # not in brainarray, try bioMart
    HGU133Plus2_Hs_ENSG_mapping2$`Probe Set Name` <- gsub(pattern = "_at",replacement = "",x = HGU133Plus2_Hs_ENSG_mapping2$`Probe Set Name`)
    # brainarray_one=intersect(colnames(scaled_data),HGU133Plus2_Hs_ENSG_mapping2$`Probe Set Name`) 
    brainarray_dict <- subset(HGU133Plus2_Hs_ENSG_mapping2,`Probe Set Name` %in% sub('\\.[0-9]*$', '', colnames(scaled_data)))# in TPM and have one ensembl to one probe set in brainarray
    
    count <- rep(0,length(weights))
    for (i in 1:length(probeset_id)){
      filtered_dict <- subset(dict,affy_hg_u133_plus_2==probeset_id[i] & ensembl_gene_id %in% colnames(scaled_data))
      count[i] <- nrow(filtered_dict)
    }
   
    first <- 0
    second <- 0
    revised_weight <- c()
    revised_ensembl <- c()
    revised_probeset <- c()
    for (i in 1:length(weights)){
      probeset <- probeset_id[i]
      if (probeset %in% brainarray_dict$`Affy Probe Set Name`){
        revised_weight <- c(revised_weight,weights[i])
        revised_ensembl <- c(revised_ensembl,brainarray_dict$`Probe Set Name`[brainarray_dict$`Affy Probe Set Name`==probeset])
        revised_probeset <- c(revised_probeset,probeset_id[i])
        first <- first+1
      }
      else if (count[i]<5 & count[i]>0 ){
        filtered_dict <- subset(dict,affy_hg_u133_plus_2==probeset_id[i] & ensembl_gene_id %in% colnames(scaled_data))
        revised_weight <- c(revised_weight,rep(weights[i]/count[i],count[i])) # take the average of multiple ensembl genes mapping to the same probe set
        revised_ensembl <- c(revised_ensembl,filtered_dict$ensembl_gene_id)
        revised_probeset <- c(revised_probeset,rep(probeset_id[i],count[i]))
        second <- second+1
      }
    }
    
    score <- rep(0,nrow(scaled_data))
    for (i in 1:length(revised_weight)){
      score <- score+revised_weight[i]*scaled_data[,revised_ensembl[i]]
    }
    
    if (signature == "GPI"){
      
    }
    
    return(score)
  }
}

gene_signature_score("UAMS70",scaled_data)



