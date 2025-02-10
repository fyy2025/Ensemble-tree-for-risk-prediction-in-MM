library(geneClassifiers)
library(biomaRt) 
library(readr)
library(tidyr)
library(dplyr)
library(janitor)
library(viridis)
library(ggplot2)
library(rawr)
library(survival)
library(survminer)

# read in the TPM data, first column name is Gene, containing gene ensembl ids, and the rest are patient names, each row is one gene
setwd("/Users/fanyiyang/Desktop/Research/summer/Ensemble-tree-for-risk-prediction-in-MM/Plotter")
MMRF_TPM <- read_tsv("./MMRF_TPM.tsv")
TPM_first_visit <- MMRF_TPM %>% dplyr::select("Gene",ends_with("_1_BM_CD138pos"))

# read in the brainarray dataset for matching gene ensembl id with probeset id
HGU133Plus2_Hs_ENSG_mapping <- read_delim("./HGU133Plus2_Hs_ENSG_mapping.txt",delim = "\t", escape_double = FALSE,  trim_ws = TRUE)

# read in the bioMart dataset for matching gene ensembl id with gene symbol id
mart <- useMart('ENSEMBL_MART_ENSEMBL')
mart <- useDataset('hsapiens_gene_ensembl', mart)
annotLookup <- getBM(
  mart = mart,
  attributes = c(
    'hgnc_symbol',
    'ensembl_gene_id',
    'gene_biotype'),
  uniqueRows = TRUE)

gene_signature_score = function(signature, TPM_first_visit, surv="OS"){
  if (!(signature %in% c("GPI","EI"))){
    Classifier <- getClassifier(signature)
    weights <- getWeights(Classifier)
    probeset_id <- names(weights)
    
    ensembl <- useEnsembl("ensembl",dataset="hsapiens_gene_ensembl",mirror="www") 
    affy_ensembl <- c("affy_hg_u133_plus_2", "ensembl_gene_id")
    dict <- getBM(attributes=affy_ensembl, filters="affy_hg_u133_plus_2",values = probeset_id, mart = ensembl)
    
    TPM_first <- subset(TPM_first_visit,Gene %in% dict$ensembl_gene_id)
    tidydata <- t(TPM_first)%>%
      row_to_names(row_number=1)
    patient <- row.names(tidydata)
    tidydata <- as.data.frame(tidydata)
    tidydata <- as.data.frame(sapply(tidydata,as.numeric ))
    scaled_data <- scale(tidydata) # scale before computing risk scores
    
    # try to match probeset id using the brainarray dataset
    HGU133Plus2_Hs_ENSG_mapping2 <- HGU133Plus2_Hs_ENSG_mapping[HGU133Plus2_Hs_ENSG_mapping$`Affy Probe Set Name` %in% probeset_id,c(1,7)]
    HGU133Plus2_Hs_ENSG_mapping2 <- HGU133Plus2_Hs_ENSG_mapping2[!duplicated(HGU133Plus2_Hs_ENSG_mapping2),]
    
    # if not in brainarray, try bioMart
    brainarray_out <- setdiff(probeset_id,HGU133Plus2_Hs_ENSG_mapping2$`Affy Probe Set Name`) 
    HGU133Plus2_Hs_ENSG_mapping2$`Probe Set Name` <- gsub(pattern = "_at",replacement = "",x = HGU133Plus2_Hs_ENSG_mapping2$`Probe Set Name`)
    brainarray_dict <- subset(HGU133Plus2_Hs_ENSG_mapping2,`Probe Set Name` %in% sub('\\.[0-9]*$', '', colnames(scaled_data)))# in TPM and have one ensembl to one probe set in brainarray
    
    count <- rep(0,length(weights))
    for (i in 1:length(probeset_id)){
      filtered_dict <- subset(dict,affy_hg_u133_plus_2==probeset_id[i] & ensembl_gene_id %in% colnames(scaled_data))
      count[i] <- nrow(filtered_dict)
    }
   
    first <- 0 # number of probesets in brainarray
    second <- 0 # number of probesets in bioMart
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
    
  }
  if (signature=="GPI"){
    # read in genes for proliferation index
    GPI_gene_list <- read.csv("./Proliferation\ Index\ Gene\ Symbols.csv")$Gene.Symbol
    GPI_dict <- subset(annotLookup,hgnc_symbol %in% GPI_gene_list)
    
    TPM_first_GPI <- subset(TPM_first_visit,Gene %in% GPI_dict$ensembl_gene_id)
    for (i in 1:5){
      TPM_first_GPI$Gene[i] <- GPI_dict$hgnc_symbol[which(GPI_dict$ensembl_gene_id==TPM_first_GPI$Gene[i])]
    }
    tidydata <- t(TPM_first_GPI)%>%
      row_to_names(row_number=1)
    patient <- row.names(tidydata)
    patient <- gsub(pattern = "_1_BM_CD138pos",replacement = "",x = patient)
    tidydata <- as.data.frame(tidydata)
    tidydata <- as.data.frame(sapply(tidydata,as.numeric ))
    scaled_data <- scale(tidydata)
    
    score <- rep(0,nrow(scaled_data))
    for (i in 1:ncol(scaled_data)){
      score <- score+scaled_data[,i]
    }
    
    return(data.frame(score,patient))
    
  }
  if (signature=="EI"){
    # read in blood parameters for calculating EI scores, should include ldh and beta2 microglobulin level
    blood <- read_tsv("./MMRF_CoMMpass_IA21_PER_PATIENT_VISIT.tsv")
    bloods <- drop_na(subset(blood,endsWith(SPECTRUM_SEQ,"_1")==1)[c("D_LAB_chem_ldh","D_LAB_serum_beta2_microglobulin","PUBLIC_ID","SPECTRUM_SEQ")])
    
    EI_gene_list <- c("APOBEC2","APOBEC3B","IL11","TGFB1","TGFB3")
    EI_dict <- subset(annotLookup,hgnc_symbol %in% EI_gene_list) # matching gene ensembl id with gene symbol id
    
    TPM_first_EI  <-  subset(TPM_first_visit,Gene %in% EI_dict$ensembl_gene_id)
    for (i in 1:5){
      TPM_first_EI$Gene[i] <- EI_dict$hgnc_symbol[which(EI_dict$ensembl_gene_id==TPM_first_EI$Gene[i])]
    }
    tidydata <- t(TPM_first_EI)%>%
      row_to_names(row_number=1)
    patient <- row.names(tidydata)
    tidydata <- as.data.frame(tidydata)
    tidydata <- as.data.frame(sapply(tidydata,as.numeric ))
    public_id <- gsub(pattern = "_1_BM_CD138pos",replacement = "",x = patient)
    tidydata$PUBLIC_ID <- public_id
    
    merged_data <- merge(bloods,tidydata,by="PUBLIC_ID")
    
    OS_score=rep(0,nrow(merged_data))
    OS_score=OS_score+ifelse(merged_data$D_LAB_chem_ldh>3.18,2,0)
    OS_score=OS_score+ifelse(merged_data$D_LAB_serum_beta2_microglobulin>4.22,4,0)
    OS_score=OS_score+ifelse(merged_data$APOBEC2>0.2,3,0)
    OS_score=OS_score+ifelse(merged_data$APOBEC3B>3.79,3,0)
    OS_score=OS_score+ifelse(merged_data$IL11>0.48,2.5,0)
    OS_score=OS_score+ifelse(merged_data$TGFB1>0.11,1,0)
    OS_score=OS_score+ifelse(merged_data$TGFB3>0.1,2,0)
    
    PFS_score=rep(0,nrow(merged_data))
    PFS_score=PFS_score+ifelse(merged_data$D_LAB_chem_ldh>3.08,2,0)
    PFS_score=PFS_score+ifelse(merged_data$D_LAB_serum_beta2_microglobulin>3.33,4,0)
    PFS_score=PFS_score+ifelse(merged_data$APOBEC2>0.19,4,0)
    PFS_score=PFS_score+ifelse(merged_data$APOBEC3B>6.91,3,0)
    PFS_score=PFS_score+ifelse(merged_data$IL11>0.56,2,0)
    PFS_score=PFS_score+ifelse(merged_data$TGFB1>0.43,1,0)
    PFS_score=PFS_score+ifelse(merged_data$TGFB3>0.03,1,0)
    
    if (surv == "OS"){
      return(data.frame(EI=OS_score,patient=merged_data$PUBLIC_ID))
    }
    else if (surv == "PFS"){
      return(data.frame(EI=PFS_score,patient=merged_data$PUBLIC_ID))
    }
  }
    
  return(score)
}

EMC <- gene_signature_score("EMC92",TPM_first_visit)
UAMS70 <- gene_signature_score("UAMS70",TPM_first_visit)
UAMS17 <- gene_signature_score("UAMS17",TPM_first_visit)
UAMS80 <- gene_signature_score("UAMS80",TPM_first_visit)
IFM15 <- gene_signature_score("IFM15",TPM_first_visit)
HM19 <- gene_signature_score("HM19",TPM_first_visit)
EI <- gene_signature_score("EI",TPM_first_visit)
EI_PFS <- gene_signature_score("EI",TPM_first_visit,surv="PFS")
GPI <- gene_signature_score("GPI",TPM_first_visit)

merged <- data.frame(EMC,UAMS70,GPI,UAMS17,UAMS80,HM19,IFM15)
merged2 <- merge(merged, EI, by="patient")

# Use 80 percentile as cutoff, count the number of signatures catagorizing each patient as high risk
binary <- merged2
count <- rep(0,nrow(binary))
for (i in 2:ncol(binary)){
  threshold <- quantile(binary[,i],0.8)
  binary[,i] <- ifelse(binary[,i]>threshold,1,0)
  count <- count+binary[,i]
}
binary$count <- count

# Produce a barplot and heatmap for those patients identified as high risk by more than zero signatures
positive_count <- subset(binary,count>0)
positive_count <- positive_count[order(positive_count$count, decreasing = TRUE), ]
positive_count$patient <- factor(positive_count$patient, levels = positive_count$patient)
ggplot(positive_count,aes(x=patient,y=count))+
  geom_bar(stat="identity",width=0.6)+
  labs(title = "Bar Plot ordered by the count of signatures indicating high risk", x = "Patient", y = "Count") +
  theme(axis.text.x=element_blank())


long_df <- positive_count[,-10] |> pivot_longer(cols = colnames(positive_count)[-c(1,10)], names_to = "X", values_to = "Z")
long_df$patient <- factor(long_df$patient,levels=unique(long_df$patient))
long_df |> ggplot(aes(patient, X, fill = as.factor(Z))) + geom_tile()+
  scale_fill_viridis_d()+
  labs(fill = "High Risk",y="Signature",x="Patient") +
  theme(axis.text.x=element_blank())

# read in survival data
survival=read_tsv("./MMRF_CoMMpass_IA21_STAND_ALONE_SURVIVAL.tsv")%>%
  dplyr::select("ttcos","censos","ttcpfs","censpfs","PUBLIC_ID")%>%
  rename(patient=PUBLIC_ID,time=ttcos,status=censos)

full_data <- merge(binary,survival,by="patient")

### This is a data frame with 3 columns,
### time: should be PFS or OS time depending your KM plot. Please convert your time to months if it is coded as days or years
### event: should be your censor information, by default 1 is event observed, 0 ins censor
### group: should be you sample groupings like mutation present vs. WT vs. deletion etc.
full_data$risk <- ifelse(full_data$count>5,1,ifelse(full_data$count>2,2,ifelse(full_data$count>0,3,4)))
survData <- full_data[,c(13,14,15)]
colnames(survData) <- c("Time","Event","Class")
survData$Time <- survData$Time/30.42
table(survData$Class)

survData <- na.omit(survData)

### This function will create survival curves and create survfit object required for plotting below.
### It simply use the data frame generated above
surv.fit <- survfit (Surv(Time,Event)~ Class, data=survData)
paircomp_OS <- coxph_pairs(Surv(Time,Event)~ Class, data=survData)

### ggsurvplot funtion creates KM plot and risk table in JCO format. 
### Please do not change the palette option
### It takes survival curves and data table generated above
### If your data has a short follow up you may consider adjusting break.time.by = 6 to add a risk table entry for every 6 months instead of 12 months 
gg_os <- ggsurvplot(
  title = "Overall Survival for MMRF dataset",
  surv.fit,                     # survfit object with calculated statistics.
  data = survData,             # data used to fit survival   curves.
  risk.table = TRUE,       # show risk table.
  tables.theme = clean_theme(),
  pval = TRUE,             # show p-value of log-rank test.
  pval.coord = c(-0.75,0.20),
  #  conf.int = TRUE,         # show confidence intervals for 
  # point estimates of survival curves.
  palette = "jco",
  # survival estimates.
  xlab = "Time in months",   # customize X axis label. 
  ylab = " OS probability",
  break.time.by = 12,     # break X axis in time intervals by 500.
  ggtheme = theme_light(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T,# color risk table text annotations.
  risk.table.height = 0.25, # the height of the risk table
  risk.table.y.text = FALSE,# show bars instead of names in text annotations
  # in legend of risk table.
  conf.int.style = "step",  # customize style of confidence intervals
  surv.median.line = "hv",  # add the median survival pointer.
  legend.labs = c(">5 signatures","3-5 signatures","1-2 signatures","0 signatures"),    # change legend labels, please double check your edited labels here and your input data!!
  #palette = c("#CD534CFF","#0073C2FF","#EFC000FF","#6a994e")
)

gg_os # KM curves and risk table for overall survival



merged2 <- merge(merged, EI_PFS, by="patient")

# Use 80 percentile as cutoff, count the number of signatures catagorizing each patient as high risk
binary <- merged2
count <- rep(0,nrow(binary))
for (i in 2:ncol(binary)){
  threshold <- quantile(binary[,i],0.8)
  binary[,i] <- ifelse(binary[,i]>threshold,1,0)
  count <- count+binary[,i]
}
binary$count <- count

# Produce a barplot and heatmap for those patients identified as high risk by more than zero signatures
positive_count <- subset(binary,count>0)
positive_count <- positive_count[order(positive_count$count, decreasing = TRUE), ]
positive_count$patient <- factor(positive_count$patient, levels = positive_count$patient)
ggplot(positive_count,aes(x=patient,y=count))+
  geom_bar(stat="identity",width=0.6)+
  labs(title = "Bar Plot ordered by the count of signatures indicating high risk", x = "Patient", y = "Count") +
  theme(axis.text.x=element_blank())


long_df <- positive_count[,-10] |> pivot_longer(cols = colnames(positive_count)[-c(1,10)], names_to = "X", values_to = "Z")
long_df$patient <- factor(long_df$patient,levels=unique(long_df$patient))
long_df |> ggplot(aes(patient, X, fill = as.factor(Z))) + geom_tile()+
  scale_fill_viridis_d()+
  labs(fill = "High Risk",y="Signature",x="Patient") +
  theme(axis.text.x=element_blank())

# read in survival data
survival=read_tsv("./MMRF_CoMMpass_IA21_STAND_ALONE_SURVIVAL.tsv")%>%
  dplyr::select("ttcos","censos","ttcpfs","censpfs","PUBLIC_ID")%>%
  rename(patient=PUBLIC_ID,time=ttcos,status=censos)

full_data <- merge(binary,survival,by="patient")

### This is a data frame with 3 columns,
### time: should be PFS or OS time depending your KM plot. Please convert your time to months if it is coded as days or years
### event: should be your censor information, by default 1 is event observed, 0 ins censor
### group: should be you sample groupings like mutation present vs. WT vs. deletion etc.
full_data$risk <- ifelse(full_data$count>5,1,ifelse(full_data$count>2,2,ifelse(full_data$count>0,3,4)))
survData <- full_data[,c(11,12,15)]
colnames(survData) <- c("Time","Event","Class")
survData$Time <- survData$Time/30.42
table(survData$Class)

survData <- na.omit(survData)

### This function will create survival curves and create survfit object required for plotting below.
### It simply use the data frame generated above
surv.fit <- survfit (Surv(Time,Event)~ Class, data=survData)
paircomp_OS <- coxph_pairs(Surv(Time,Event)~ Class, data=survData)

### ggsurvplot funtion creates KM plot and risk table in JCO format. 
### Please do not change the palette option
### It takes survival curves and data table generated above
### If your data has a short follow up you may consider adjusting break.time.by = 6 to add a risk table entry for every 6 months instead of 12 months 
gg_pfs <- ggsurvplot(
  title = "Overall Survival for MMRF dataset",
  surv.fit,                     # survfit object with calculated statistics.
  data = survData,             # data used to fit survival   curves.
  risk.table = TRUE,       # show risk table.
  tables.theme = clean_theme(),
  pval = TRUE,             # show p-value of log-rank test.
  pval.coord = c(-0.75,0.20),
  #  conf.int = TRUE,         # show confidence intervals for 
  # point estimates of survival curves.
  palette = "jco",
  # survival estimates.
  xlab = "Time in months",   # customize X axis label. 
  ylab = " OS probability",
  break.time.by = 12,     # break X axis in time intervals by 500.
  ggtheme = theme_light(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T,# color risk table text annotations.
  risk.table.height = 0.25, # the height of the risk table
  risk.table.y.text = FALSE,# show bars instead of names in text annotations
  # in legend of risk table.
  conf.int.style = "step",  # customize style of confidence intervals
  surv.median.line = "hv",  # add the median survival pointer.
  legend.labs = c(">5 signatures","3-5 signatures","1-2 signatures","0 signatures"),    # change legend labels, please double check your edited labels here and your input data!!
  #palette = c("#CD534CFF","#0073C2FF","#EFC000FF","#6a994e")
)

gg_pfs # KM curves and risk table for overall survival
