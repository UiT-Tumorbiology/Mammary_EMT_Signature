#Import the libraries
library(dplyr)
library(forcats)
library(biomaRt)
library(readxl)
#library(tidyverse)
library(dendsort)
library(ComplexHeatmap)
library(circlize)
expr <- read.csv2("metabric_mrna.csv", as.is = T, check.names = F)
#Using entrez id to get the ensembl id
meta_entrez <- expr$Entrez_Gene_Id

#Using mart to get the ensembl ids for the entrez ids
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
G_list_ensembl <- getBM(filters= "entrezgene_id", attributes= c("ensembl_gene_id",
                                                                "entrezgene_id"),
                        values=meta_entrez,mart= mart)
emt <- read.csv2("265_EMT.csv", as.is = T, check.names = F) #265 EMT Genes
ensembl <- emt$Identifier

overlap <- subset(G_list_ensembl, (ensembl_gene_id %in% ensembl))

emt_metabric <- subset(expr,(Entrez_Gene_Id %in% overlap$entrezgene_id))

emt_metabric <- distinct(emt_metabric, Entrez_Gene_Id, .keep_all = T)

rownames(emt_metabric) <- emt_metabric$Hugo_Symbol
emt_metabric <- emt_metabric[,-c(1:2)]
emt_metabric <- as.data.frame(t(scale(t(emt_metabric))))
#Getting the clinical data of these patients
clinical <- as.data.frame(read_xlsx("brca_metabric_clinical_data.tsv.xlsx"))
clinical <- subset(clinical, (`Patient ID` %in% colnames(emt_metabric)))
#Only taking the subtypes not classified as NC
clinical <- clinical %>% filter(!(`Pam50 + Claudin-low subtype` == "NC" | `Pam50 + Claudin-low subtype` == "NA"))
emt_metabric <- dplyr::select(emt_metabric, clinical$`Patient ID`)
all(colnames(emt_metabric)==clinical$`Patient ID`)

annotation2 <- read.csv2("Supplementary_table_6.csv", as.is = T, check.names = F)

all(annotation2$Name == rownames(emt_metabric))

#Complexheatmap needs a matrix, so converting the dataframe to matrix
h <- as.matrix(emt_metabric)


#Reposition the breaks in quantile positions
quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}

#Make clusters based on their pearson correlation
row_dist <- as.dist(1-cor(t(h), method = "pearson"))
col_dist <- as.dist(1-cor(h, method = "pearson"))
col_hc <- hclust(col_dist, method = "complete")
row_hc <- hclust(row_dist, method = "complete")

#We can use the dendsort package and reorder the clustering
Rowv=dendsort(as.dendrogram(row_hc), isRevers=TRUE, type = "average") 
Colv=dendsort(as.dendrogram(col_hc), type = "average")

col_fun <- colorRamp2(quantile_breaks(h, n = 11), c("#1e90ff", "#0479f5", "#0467d2", "#0356af", "#02458c","#023369","#012246", "#130202", "#6b0002","#b20004","#ff0000"))

#Setting up the annotation files to be used in the heatmap
#Column Annotation
ha <- HeatmapAnnotation(Subtype=clinical$`Pam50 + Claudin-low subtype`,
                 col=list(Subtype=c("Basal"="#F22233",
                                           "Normal"="#8C51A6",
                                           "Her2"="#3C5FA6",
                                           "LumA"="sienna1",
                                           "LumB"="darkseagreen4",
                                           "claudin-low"="#00FFCC",
                                           "NC"="black")),
                        simple_anno_size = unit(0.30, "cm"),
                        border = T,
                        annotation_name_gp= gpar(fontsize = 8.5,fontface="bold"),
                        show_legend = T,
                        annotation_legend_param = list(title_gp = gpar(fontsize=10, fontface="bold"), labels_gp=gpar(fontsize=8)))

#Row annotation
hr <- rowAnnotation("Up/Down"=annotation2$`Up/Down regulated in cell lines`,
                    "Gene Cluster - TCGA" = annotation2$`Gene cluster patients`,
                    col=list("Gene Cluster - TCGA"=c("cluster1"="cadetblue2",
                                                     "cluster2"="burlywood2",
                                                     "cluster3"="#ff3300"),
                             "Up/Down"=c("Down" = "#1e90ff",
                                                    "Up" = "#ff0000")),
                    simple_anno_size = unit(0.30, "cm"),
                    border = T,
                    annotation_name_gp= gpar(fontsize = 8.5,fontface="bold"),
                    show_legend = T,
                    annotation_legend_param = list(title_gp = gpar(fontsize=10, fontface="bold"), labels_gp=gpar(fontsize=8)))
#Final Heatmap
ht <- Heatmap(h,col = col_fun,
              cluster_columns = Colv,
              name = "Expression Values",
              show_heatmap_legend = T ,
              top_annotation = ha,
              left_annotation = hr,
              show_column_names = F,
              show_row_names = F,cluster_rows = Rowv,
              row_title_gp = gpar(fontsize=12),
              column_title_gp = gpar(fontsize=12),
              height = unit(12, "cm"),
              width = unit(10.55, "cm"),
              column_dend_reorder = T,
              row_dend_reorder = T,
              show_row_dend = T,
              border = T,
              column_dend_height = unit(2, "cm"),
              column_names_rot = 90,
              column_names_gp = grid::gpar(fontsize = 7,fontface="bold"),
              heatmap_legend_param = list(title="Scaled Expression", legend_height=unit(3, "cm"),title_gp=gpar(fontsize=10, fontface="bold"),labels_gp = gpar(fontsize=8)),
              row_split = 4, column_split = 4,
              column_gap = unit(c(0.5,0.5,0.5,0.5), "mm"),
              row_gap = unit(c(0.5,0.5,0.5, 0.5), "mm"))
           
#Exporting to a pdf file
pdf("Metabric_4_heatmap.pdf",  width=13,height=16)

draw(ht, merge_legend=TRUE, padding = unit(c(2, 2, 2, 2), "mm"))

dev.off()

#To find the gene and patient orders, draw the heatmap object first so that it will not change with every run
ht <- draw(ht)

#For getting the column orders
for (i in 1:length(column_order(ht))){   if (i == 1) {
  clu <- t(t(colnames(emt_metabric[,column_order(ht)[[i]]])))
  out <- cbind(clu, paste("cluster", i, sep=""))
  colnames(out) <- c("Patient_ID", "Cluster")   } else {
    clu <- t(t(colnames(emt_metabric[,column_order(ht)[[i]]])))
    clu <- cbind(clu, paste("cluster", i, sep=""))
    out <- rbind(out, clu)   } 
}

write.csv2(out, "Patient_orders_Metabric.csv")

#For getting the gene orders
for (i in 1:length(row_order(ht))){   if (i == 1) {
  clu <- t(t(row.names(emt_metabric[row_order(ht)[[i]],])))
  out1 <- cbind(clu, paste("cluster", i, sep=""))
  colnames(out1) <- c("GeneID", "Cluster")   } else {
    clu <- t(t(row.names(emt_metabric[row_order(ht)[[i]],])))
    clu <- cbind(clu, paste("cluster", i, sep=""))
    out1 <- rbind(out1, clu)   } 
}
out1

write.csv(out1, file = "Gene_order_Metabric.csv")