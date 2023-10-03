library(survival)
library(ggplot2)
library(dplyr)
library(survminer)
library(readxl)

clinical <- read.csv2("Final_clinical_data_types.csv",as.is = T,check.names = F, row.names = 1)

clinical$DSS_year <- as.numeric(clinical$DSS_year)
clinical$PFI_year <- as.numeric(clinical$PFI_Year)


Basal <- subset(clinical,`Pam50 + Claudin-low subtype`%in% c("Basal", "claudin-low"))
Basal <- subset(Basal, Basal$Cluster == "Cluster 1" | Basal$Cluster == "Cluster 3")

survival <- Surv(time = Basal$DSS_year, event = Basal$DSS)
survival_fit <- survfit(formula = survival ~ Cluster, data = Basal)

Basal_plot <- ggsurvplot(fit = survival_fit, pval = TRUE, title="Basal in METABRIC", 
                         xlab = "Years", ylab = "DSS Probability", ylim=c(0.4,1), 
                         xlim=c(0,5), break.x.by=1, 
                         pval.coord=c(0,0.42), 
                         legend="right",
                         legend.labs=c("Cluster 1","Cluster 3"), conf.int = T,risk.table = T, legend.title="Patient Cluster", conf.int.alpha = c(0.2),
                         palette = c("#d28168","#608244"), risk.table.y.text = TRUE, risk.table.title="",
                         risk.table.height=0.15, risk.table.fontsize=3.5, tables.theme = theme(axis.text.x = element_blank(), axis.ticks.x=element_blank(), 
                                                                                               axis.title.x = element_blank(), axis.title.y = element_blank(),
                                                                                               axis.text.y = element_text(size = 10)),
                         font.x=10, font.y=11, font.tickslab=10)
Basal_plot

pdf("All Basal DSS_METABRIC.pdf", width = 8, height = 7, onefile = F)
print(Basal_plot)
dev.off()

survival <- Surv(time = Basal$PFI_Year, event = Basal$`Relapse Free Status`)
survival_fit <- survfit(formula = survival ~ Cluster, data = Basal)

Basal_plot <- ggsurvplot(fit = survival_fit, pval = TRUE, title="Basal", 
                         xlab = "PFI Time (Years)", ylab = "PFI Probability", ylim=c(0.5,1), xlim=c(0,5), break.x.by=1, pval.coord=c(0,0.52), legend="right",
                         legend.labs=c("Cluster 1","Cluster 3"), conf.int = T,risk.table = T, legend.title="Patient Cluster", conf.int.alpha = c(0.2))
Basal_plot

pdf("All Basal PFI_METABRIC.pdf", width = 8, height = 7, onefile = F)
print(Basal_plot)
dev.off()
