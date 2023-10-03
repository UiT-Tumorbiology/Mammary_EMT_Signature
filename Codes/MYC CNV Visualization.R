library(dplyr)
library(ggpubr)
library(tibble)
library(tidyr)

#CNV profiles were downloaded from Xena browser
cnv <- read.table("TCGA-BRCA_CopyNumber_Gistic2_all_thresholded.by_genes", sep = "\t", header = TRUE, row.names = 1, check.names = F)
clusters <- read.csv2("TCGA_Patient_Cluster.csv", as.is = T, check.names = F)
subtype <- read.csv2("TCGA_Patient_Info.csv", as.is = T, check.names = F)
clusters <- merge(clusters, subtype, by = "Patient_ID")

#Adding "-01" to patient ids of the metadata object to match the naming of CNV
clusters$Patient_ID <- paste(clusters$Patient_ID, "-01",sep = "")

#Taking out cluster 5 patients
cl5 <- subset(clusters, (TCGA_Cluster == "Cluster 5"))
cl5 <- subset(cl5, (Subtype == "Basal-like"))
cl5_cnv <- select(cnv, cl5$Patient_ID)
cl5_cnv <- cl5_cnv["MYC",]

cl5_cnv <- tibble::add_column(cl5_cnv, GeneID="MYC", .before=1) #Adding a column

cl5_cnv_long <- pivot_longer(cl5_cnv, cols = 2:128, names_to = "PatientID", values_to = "CNV")

cl5_cnv_long$Cluster <- "Cluster 5"

#Taking out cluster 1 patients
cl1 <- subset(clusters, (TCGA_Cluster == "Cluster 1"))
cl1 <- subset(cl1, (Subtype == "Basal-like"))
cl1_cnv <- select(cnv, cl1$Patient_ID)
cl1_cnv <- cl1_cnv["MYC",]
cl1_cnv <- add_column(cl1_cnv, GeneID="MYC", .before=1)
cl1_cnv_long <- pivot_longer(cl1_cnv, cols = 2:50, names_to = "PatientID", values_to = "CNV")
cl1_cnv_long$Cluster <- "Cluster 1"

combined <- rbind(cl1_cnv_long, cl5_cnv_long)

#Writing a function to add column explaining the CNV status
map_cnv_to_comment <- function(value) {
  if (value < 0 ) {
    return("Del")
  } else if (value > 0) {
    return("Amp")
  } else {
    return("Norm")
  }
}

combined$Status <- sapply(combined$CNV, map_cnv_to_comment)

percent_data<-combined %>%
  select(`Cluster`,`Status`)%>%
  table()%>%
  as.data.frame()%>%
  group_by(`Cluster`)%>%
  mutate(percent = paste0(round(Freq/sum(Freq)*100,1)))%>%
  as.data.frame()
percent_data$percent <- as.numeric(percent_data$percent)
percent_data$Cluster <-factor(percent_data$Cluster, levels = c("Cluster 1", "Cluster 5"))

#Changing the order so that higher values stay on top of the lower ones
status_order <- percent_data %>%
  group_by(Status) %>%
  summarise(median_percent = median(percent)) %>%
  arrange(desc(median_percent)) %>%
  pull(Status)


percent_data$Status <- factor(percent_data$Status, levels = status_order)

p <- ggbarplot(percent_data, x="Cluster", y="percent", fill = "Status", palette = c("#608244","#1e90ff","#d28168"),
               position = position_stack(0.4), 
               xlab = "Cluster", ylab = "Percentage",
               width = 0.3) + 
  #scale_x_discrete(labels = my_labels) + 
  #scale_y_continuous(labels = function(x) paste0(x, "%")) +
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size=10), axis.ticks.x=element_blank(), 
        axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 11),
        axis.text.y = element_text(size = 10)) + ggtitle("MYC Copy-number")
pdf("MYC CNV Number.pdf", width = 4, height = 4)
print(p)
dev.off()
