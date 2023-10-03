library(survival)
library(survminer)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(lubridate)
library(grid)
library(gridExtra)

clinical <- read.csv2("TCGA_BRCA_Survival_myc.csv", as.is = T, check.names = F, row.names = 1) #TCGA BRCA Clinical data with cluster information

clinical$DSS_Year <- clinical$DSS_Days/365
clinical <- clinical %>% filter(Subtype == "Basal")
clinical <- clinical %>% filter(Cluster == "Cluster 1" | Cluster == "Cluster 5")

survival <- Surv(time = clinical$DSS_Year, event = clinical$DSS)
survival

survival_fit <- survfit(formula = survival ~ Cluster, data = clinical)
survival_fit

Basal <- ggsurvplot(fit = survival_fit, pval = TRUE, legend = "right",title = "Basal-like Patients",
           xlab = "Years", ylab = "DSS Probability", ylim=c(0.6,1), xlim=c(0,5), 
            pval.coord=c(0,0.62), conf.int = T, risk.table = T, legend.labs=c("Cluster 1", "Cluster 5"),
           legend.title="Patient Cluster", palette = c("#d28168","#608244"), risk.table.y.text = TRUE, risk.table.title="",
           risk.table.height=0.15, risk.table.fontsize=3.5, tables.theme = theme(axis.text.x = element_blank(), axis.ticks.x=element_blank(), 
                                                       axis.title.x = element_blank(), axis.title.y = element_blank(),
                                                       axis.text.y = element_text(size = 10)),
           font.x=10, font.y=11, font.tickslab=10)
Basal

pdf("All Basal DSS_TCGA.pdf", width = 8, height = 7, onefile = F)
print(Basal)
dev.off()

#PFI
survival <- Surv(time = clinical$PFI_Year, event = clinical$PFI)
survival
survival_fit <- survfit(formula = survival ~ Cluster, data = clinical)
survival_fit
Basal <- ggsurvplot(fit = survival_fit, pval = TRUE, legend = "right",title = "Basal-like Patients",
                    xlab = "Years", ylab = "PFI Probability", ylim=c(0.5,1), xlim=c(0,5), 
                    pval.coord=c(0,0.52), conf.int = T, risk.table = T, legend.labs=c("Cluster 1", "Cluster 5"),
                    legend.title="Patient Cluster", palette = c("#d28168","#608244"), risk.table.y.text = TRUE, risk.table.title="",
                    risk.table.height=0.15, risk.table.fontsize=3.5, tables.theme = theme(axis.text.x = element_blank(), axis.ticks.x=element_blank(), 
                                                                                          axis.title.x = element_blank(), axis.title.y = element_blank(),
                                                                                          axis.text.y = element_text(size = 10)),
                    font.x=10, font.y=11, font.tickslab=10)
Basal

pdf("All Basal PFI_TCGA.pdf", width = 8, height = 7, onefile = F)
print(Basal)
dev.off()


