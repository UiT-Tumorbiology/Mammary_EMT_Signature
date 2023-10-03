library(dplyr)
library(ggplot2)
library(clValid)
library(readxl)

#TCGA
expr <- read_xlsx("TCGA - 265 EMT genes - z-normalized.xlsx") #Scaled count matrix of 265 Genes in TCGA
expr <- expr[c(1:265),c(1,3:1043)]
expr <- as.data.frame(expr)
rownames(expr) <- expr[,1]
expr <- expr[,-1]

internal <- clValid::clValid(t(expr), method = "complete", metric = "correlation", clMethods = "hierarchical", nClust = 2:10, validation = "internal")

##Visualization
pdf("TCGA_Stability.pdf", width = 5, height = 5, fonts = "serif")
plot(internal, legend = FALSE, main = "")
dev.off()


#Metabric
expr <- read.csv2("EMT_scaled_METABRIC.csv", sep = ";", as.is = T,check.names = F)
rownames(expr) <- expr[,1]
expr <- expr[,-c(1:2)]

internal <- clValid::clValid(t(expr), method = "complete", metric = "correlation", clMethods = "hierarchical", nClust = 2:10, validation = "internal")


pdf("METABRIC_Stability.pdf", width = 5, height = 5, fonts = "serif")
plot(internal, legend = FALSE, main = "")
dev.off()

#CCLE
expr <- read.csv2("CCLE_EMT_Expression_Scaled.csv", as.is = T, check.names = F)
rownames(expr) <- expr$Name
expr$Name <- NULL

internal <- clValid::clValid(t(expr), method = "complete", metric = "correlation", clMethods = "hierarchical", nClust = 2:10, validation = "internal")

pdf("CCLE_Stability.pdf", width = 5, height = 5, fonts = "serif")
plot(internal, legend = FALSE, main = "")
dev.off()
