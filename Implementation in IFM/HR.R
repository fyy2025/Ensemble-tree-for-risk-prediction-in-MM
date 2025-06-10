library(readr)
setwd("/Users/fanyiyang/Desktop/Research/summer/Ensemble-tree-for-risk-prediction-in-MM/Implementation\ in\ IFM")
IFM <- read.csv("./merged_IFM_scores.csv")
binary <- IFM
for (i in 1:(ncol(binary)-4)){
  threshold = quantile(binary[,i],0.8)
  binary[,i] = ifelse(binary[,i]>threshold,1,0)
}

GPI_threshold = quantile(IFM[,4],0.92)
binary[,4] = ifelse(IFM[,4]>GPI_threshold,1,0)

extract_summary <- function(model, model_label) {
  summary_model <- summary(model)
  data.frame(
    Variable = rownames(summary_model$coefficients),
    HR = summary_model$coefficients[, "exp(coef)"],
    Lower = summary_model$conf.int[, "lower .95"],
    Upper = summary_model$conf.int[, "upper .95"],
    Model = model_label
  )
}

EMC_model <- coxph(Surv(TTRelpase, Relapse) ~ EMC,data = binary)
UAMS70_model <- coxph(Surv(TTRelpase, Relapse) ~ UAMS70,data = binary)
UAMS80_model <- coxph(Surv(TTRelpase, Relapse) ~ UAMS80,data = binary)
UAMS17_model <- coxph(Surv(TTRelpase, Relapse) ~ UAMS17,data = binary)
IFM15_model <- coxph(Surv(TTRelpase, Relapse) ~ IFM15,data = binary)
HM19_model <- coxph(Surv(TTRelpase, Relapse) ~ HM19,data = binary)
GPI_model <- coxph(Surv(TTRelpase, Relapse) ~ GPI,data = binary)

df1 <- extract_summary(EMC_model, "EMC92")
df2 <- extract_summary(UAMS70_model, "UAMS70")
df3 <- extract_summary(UAMS17_model, "UAMS17")
df4 <- extract_summary(UAMS80_model, "UAMS80")
df5 <- extract_summary(IFM15_model, "IFM15")
df6 <- extract_summary(HM19_model, "HM19")
df7 <- extract_summary(GPI_model, "GPI")


combined_df <- rbind(df1, df2, df3, df4, df5, df6, df7)

library(ggplot2)

ggplot(combined_df, aes(x = Variable, y = HR, ymin = Lower, ymax = Upper, color = Model)) +
  geom_pointrange(position = position_dodge(width = 0.7)) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  coord_flip() +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(title = "Forest Plot of individual signature models in IFM PFS",
       y = "Hazard Ratio (HR)",
       x = "")

# OS

EMC_model <- coxph(Surv(time,status) ~ EMC,data = binary)
UAMS70_model <- coxph(Surv(time,status) ~ UAMS70,data = binary)
UAMS80_model <- coxph(Surv(time,status) ~ UAMS80,data = binary)
UAMS17_model <- coxph(Surv(time,status) ~ UAMS17,data = binary)
IFM15_model <- coxph(Surv(time,status) ~ IFM15,data = binary)
HM19_model <- coxph(Surv(time,status) ~ HM19,data = binary)
GPI_model <- coxph(Surv(time,status) ~ GPI,data = binary)

df1 <- extract_summary(EMC_model, "EMC92")
df2 <- extract_summary(UAMS70_model, "UAMS70")
df3 <- extract_summary(UAMS17_model, "UAMS17")
df4 <- extract_summary(UAMS80_model, "UAMS80")
df5 <- extract_summary(IFM15_model, "IFM15")
df6 <- extract_summary(HM19_model, "HM19")
df7 <- extract_summary(GPI_model, "GPI")

combined_df <- rbind(df1, df2, df3, df4, df5, df6, df7)

ggplot(combined_df, aes(x = Variable, y = HR, ymin = Lower, ymax = Upper, color = Model)) +
  geom_pointrange(position = position_dodge(width = 0.7)) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  coord_flip() +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(title = "Forest Plot of individual signature models in IFM OS",
       y = "Hazard Ratio (HR)",
       x = "")