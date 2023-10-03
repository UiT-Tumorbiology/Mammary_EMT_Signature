library(dplyr)
library(forcats)
library(biomaRt)
library(ggpubr)
library(psych)


#TCGA
emt <- read.csv2("265_EMT_raw_tcga_id.csv", as.is = T, check.names = F, row.names = 1) #Taking EMT expression data from TCGA in linear form

emt1 <- as.data.frame(t(emt)) #emt is the expression matrix
###
cor_results <- data.frame(Gene = character(0), Correlation = numeric(0))

# Calculating correlations between "ZEB1" and all other genes
for (gene_name in colnames(emt1)) {
  if (gene_name != "ENSG00000148516") {  # Exclude self-correlation
    corr_result <- corr.test(x = emt1$ENSG00000148516, y = emt1[, gene_name], adjust = "holm")
    cor_results <- rbind(cor_results, data.frame(Gene = gene_name, Correlation = corr_result$r))
  }
}

cor_results <- cor_results[order(-cor_results$Correlation), ]
info <- read.csv2("Supplementary_table_6.csv", as.is = T, check.names = F, row.names = 1)

zeb1 <- merge(cor_results, info, by.x = "Gene", by.y = "Identifier")

zeb1$State <- ifelse(zeb1$`Gene cluster, TCGA breast cancer patients` == "Cluster A", "EMT-Down",
                      ifelse(zeb1$`Gene cluster, TCGA breast cancer patients` == "Cluster C", "EMT-Up", "Partial EMT"))

zeb1$State <- factor(zeb1$State, levels = c("EMT-Down", "Partial EMT", "EMT-Up"))
library(ggpubr)
set.seed(12)

p <- ggboxplot(data = zeb1, x="State", y="Correlation", color = "State", palette = c("#1e90ff", "darkseagreen4", "#ff0000"), outlier.shape = NA,
               ylab = "Pearson Correlation Coefficient") + 
    geom_jitter(aes(group=State, color=State, stroke=0.5), position = position_jitter(width = 0.2, height = 0.2)) + 
  stat_compare_means(label = "p.signif") + 
  coord_cartesian(ylim = c(-1,1)) + labs(title = "ZEB1")+ 
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, face ="italic"),
        axis.text.x = element_text(size=10), axis.ticks.x=element_blank(), 
        axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 11),
        axis.text.y = element_text(size = 10))

p

pdf("ZEB1 Correlation_TCGA_new.pdf", width = 5, height = 5, onefile = F)
print(p)
dev.off()



###SNAI1
cor_results <- data.frame(Gene = character(0), Correlation = numeric(0))

# Calculate correlations between "SNAI1" and all other genes
for (gene_name in colnames(emt1)) {
  if (gene_name != "ENSG00000124216") {  # Exclude self-correlation
    corr_result <- corr.test(x = emt1$ENSG00000124216, y = emt1[, gene_name], adjust = "holm")
    cor_results <- rbind(cor_results, data.frame(Gene = gene_name, Correlation = corr_result$r))
  }
}

cor_results <- cor_results[order(-cor_results$Correlation), ]

SNAI1 <- merge(cor_results, info, by.x = "Gene", by.y = "Identifier")
SNAI1$State <- ifelse(SNAI1$`Gene cluster, TCGA breast cancer patients` == "Cluster A", "EMT-Down",
                      ifelse(SNAI1$`Gene cluster, TCGA breast cancer patients` == "Cluster C", "EMT-Up", "Partial EMT"))

SNAI1$State <- factor(SNAI1$State, levels = c("EMT-Down", "Partial EMT", "EMT-Up"))
library(ggpubr)
set.seed(12)
p <- ggboxplot(data = SNAI1, x="State", y="Correlation", color = "State", palette = c("#1e90ff", "darkseagreen4", "#ff0000"), outlier.shape = NA,
               ylab = "Pearson Correlation Coefficient") + 
  geom_jitter(aes(group=State, color=State, stroke=0.5), position = position_jitter(width = 0.2, height = 0.2)) + 
  stat_compare_means(label = "p.signif") + 
  coord_cartesian(ylim = c(-1,1)) + labs(title = "SNAI1")+ 
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, face ="italic"),
        axis.text.x = element_text(size=10), axis.ticks.x=element_blank(), 
        axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 11),
        axis.text.y = element_text(size = 10))
p

pdf("SNAI1 Correlation_TCGA_new.pdf", width = 5, height = 5, onefile = F)
print(p)
dev.off()

##METABRIC
#METABRIC
emt <- read.csv2("265_EMT_raw_METABRIC.csv", as.is = T, check.names = F, row.names = 1) #Taking EMT expression data from METABRIC in linear form

emt1 <- as.data.frame(t(emt)) #emt is the expression matrix
###
cor_results <- data.frame(Gene = character(0), Correlation = numeric(0))

# Calculating correlations between "ZEB1" and all other genes
for (gene_name in colnames(emt1)) {
  if (gene_name != "ZEB1") {  # Exclude self-correlation
    corr_result <- corr.test(x = emt1$ZEB1, y = emt1[, gene_name], adjust = "holm")
    cor_results <- rbind(cor_results, data.frame(Gene = gene_name, Correlation = corr_result$r))
  }
}

cor_results <- cor_results[order(-cor_results$Correlation), ]
info <- read.csv2("Supplementary_table_6.csv", as.is = T, check.names = F, row.names = 1)

zeb1 <- merge(cor_results, info, by.x = "Gene", by.y = "Identifier")

zeb1$State <- ifelse(zeb1$`Gene cluster, METABRIC breast cancer patients` == "Cluster A1" | zeb1$`Gene cluster, METABRIC breast cancer patients` == "Cluster A2", "EMT-Down",
                     ifelse(zeb1$`Gene cluster, METABRIC breast cancer patients` == "Cluster C", "EMT-Up", "Partial EMT"))

zeb1$State <- factor(zeb1$State, levels = c("EMT-Down", "Partial EMT", "EMT-Up"))
library(ggpubr)
set.seed(12)

p <- ggboxplot(data = zeb1, x="State", y="Correlation", color = "State", palette = c("#1e90ff", "darkseagreen4", "#ff0000"), outlier.shape = NA,
               ylab = "Pearson Correlation Coefficient") + 
  geom_jitter(aes(group=State, color=State, stroke=0.5), position = position_jitter(width = 0.2, height = 0.2)) + 
  stat_compare_means(label = "p.signif") + 
  coord_cartesian(ylim = c(-1,1)) + labs(title = "ZEB1")+ 
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, face ="italic"),
        axis.text.x = element_text(size=10), axis.ticks.x=element_blank(), 
        axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 11),
        axis.text.y = element_text(size = 10))

p

pdf("ZEB1 Correlation_METABRIC_new.pdf", width = 5, height = 5, onefile = F)
print(p)
dev.off()



###SNAI1
cor_results <- data.frame(Gene = character(0), Correlation = numeric(0))

# Calculate correlations between "SNAI1" and all other genes
for (gene_name in colnames(emt1)) {
  if (gene_name != "SNAI1") {  # Exclude self-correlation
    corr_result <- corr.test(x = emt1$SNAI1, y = emt1[, gene_name], adjust = "holm")
    cor_results <- rbind(cor_results, data.frame(Gene = gene_name, Correlation = corr_result$r))
  }
}

cor_results <- cor_results[order(-cor_results$Correlation), ]

SNAI1 <- merge(cor_results, info, by.x = "Gene", by.y = "Identifier")
SNAI1$State <- ifelse(SNAI1$`Gene cluster, METABRIC breast cancer patients` == "Cluster A1" | SNAI1$`Gene cluster, METABRIC breast cancer patients` == "Cluster A2","EMT-Down",
                      ifelse(SNAI1$`Gene cluster, METABRIC breast cancer patients` == "Cluster C", "EMT-Up", "Partial EMT"))

SNAI1$State <- factor(SNAI1$State, levels = c("EMT-Down", "Partial EMT", "EMT-Up"))
library(ggpubr)
set.seed(12)
p <- ggboxplot(data = SNAI1, x="State", y="Correlation", color = "State", palette = c("#1e90ff", "darkseagreen4", "#ff0000"), outlier.shape = NA,
               ylab = "Pearson Correlation Coefficient") + 
  geom_jitter(aes(group=State, color=State, stroke=0.5), position = position_jitter(width = 0.2, height = 0.2)) + 
  stat_compare_means(label = "p.signif") + 
  coord_cartesian(ylim = c(-1,1)) + labs(title = "SNAI1")+ 
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, face ="italic"),
        axis.text.x = element_text(size=10), axis.ticks.x=element_blank(), 
        axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 11),
        axis.text.y = element_text(size = 10))
p

pdf("SNAI1 Correlation_METABRIC_new.pdf", width = 5, height = 5, onefile = F)
print(p)
dev.off()

