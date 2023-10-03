#Import the libraries
library(ComplexHeatmap)
library(RColorBrewer)
library(gplots)
library(circlize)
library(dendsort)
library(dplyr)
library(ggplot2)
library(readxl)
#Import the data
expr <- read.csv2("TCGA_log.csv", as.is = T, check.names = F) #TCGA BRCA Count Matrix
annotation1 <- read.csv2("TCGA_Patient_info.csv", as.is = T, check.names = F) 
annotation2 <- read.csv2("Supplementary_table_6.csv", as.is = T, check.names = F)
#Subset the 265 genes based on ensembl id
heatmap_data <- subset(expr, (id %in% annotation2$Identifier))

#Transform the expression values into linear form
rownames(heatmap_data) <- heatmap_data$id
heatmap_data <- heatmap_data[,-c(1,2)]
heatmap_data <- 2^(heatmap_data - 1)
heatmap_data <- as.data.frame(t(scale(t(heatmap_data))))

#Setting expression matrix and information files in same order
rownames(annotation1) <- annotation1$Patients
rownames(annotation2) <- annotation2$Identifier
annotation1 <- annotation1[colnames(heatmap_data),]
annotation2 <- annotation2[rownames(heatmap_data),]
all(rownames(annotation2) == rownames(heatmap_data))
all(rownames(annotation1) == colnames(heatmap_data))

write.csv2(heatmap_data, "TCGA - 265 EMT genes_Scaled.csv")
#Complexheatmap needs a matrix, so convert the dataframe to matrix
h <- as.matrix(heatmap_data)

#Reposition the breaks in quantile positions
quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}

#Hierarchical Clustering
row_dist <- as.dist(1-cor(t(h), method = "pearson"))
col_dist <- as.dist(1-cor(h, method = "pearson"))
col_hc <- hclust(col_dist, method = "complete")
row_hc <- hclust(row_dist, method = "complete")

#Sorting the dendogram
Rowv=dendsort(as.dendrogram(row_hc), isRevers=TRUE, type = "average") 
Colv=dendsort(as.dendrogram(col_hc), type = "average")
#Choose the color gradient
col_fun <- colorRamp2(quantile_breaks(h, n = 11), c("#1e90ff", "#0479f5", "#0467d2", "#0356af", "#02458c","#023369","#012246", "#130202", "#6b0002","#b20004","#ff0000"))
ha <- HeatmapAnnotation(Subtype=annotation1$Subtype,
                        col=list(Subtype=c("Basal-like"="#F22233",
                                           "Normal-like"="#8C51A6",
                                           "Her2"="#3C5FA6",
                                           "Luminal A"="sienna1",
                                           "Luminal B"="darkseagreen4")),
                        simple_anno_size = unit(0.30, "cm"),
                        border = T,
                        annotation_name_gp= gpar(fontsize = 8.5,fontface="bold"),
                        show_legend = T,
                        annotation_legend_param = list(title_gp = gpar(fontsize=10, fontface="bold"), labels_gp=gpar(fontsize=8)))

hr <- rowAnnotation("Up/Down"=annotation2$`Up/down regulated in EMT models`,
                    "Gene Cluster - Cell Lines" = annotation2$`Gene cluster, breast cancer cell lines`,
                    col=list("Gene Cluster - Cell Lines"=c("Cluster A"="cadetblue2",
                                                           "Cluster B"="burlywood2",
                                                           "Cluster C"="chocolate2",
                                                           "#N/A" = "#CCCCCC"),
                             "Up/Down"=c("DOWN" = "#1e90ff",
                                         "UP" = "#ff0000")),
                    simple_anno_size = unit(0.30, "cm"),
                    border = T,
                    annotation_name_gp= gpar(fontsize = 8.5,fontface="bold"),
                    show_legend = T,
                    annotation_legend_param = list(title_gp = gpar(fontsize=10, fontface="bold"), labels_gp=gpar(fontsize=8)))

#Lets make the heatmap now
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
              heatmap_legend_param = list(title="Expression", legend_height=unit(3, "cm"),title_gp=gpar(fontsize=10, fontface="bold"),labels_gp = gpar(fontsize=8)),
              row_split = 3, column_split = 6,
              column_gap = unit(c(0.5,0.5,0.5,0.5,0.5), "mm"),
              row_gap = unit(c(0.5,0.5,0.5), "mm"))

pdf("TCGA_added_final.pdf",  width=13,height=16)
draw(ht, merge_legend=TRUE, padding = unit(c(2, 2, 2, 2), "mm"))
dev.off()
draw(ht, merge_legend=TRUE)
#Print the heatmap
ht <- draw(ht)



for (i in 1:length(column_order(ht))){   if (i == 1) {
  clu <- t(t(colnames(heatmap_data[,column_order(ht)[[i]]])))
  out <- cbind(clu, paste("Cluster", i, sep=" "))
  colnames(out) <- c("Patient_ID", "TCGA Cluster")   } else {
    clu <- t(t(colnames(heatmap_data[,column_order(ht)[[i]]])))
    clu <- cbind(clu, paste("Cluster", i, sep=" "))
    out <- rbind(out, clu)   } 
}
out

write.csv2(out, "TCGA_Patient_Cluster.csv")

#Gene cluster order
for (i in 1:length(row_order(ht))){   if (i == 1) {
  clu <- t(t(row.names(heatmap_data[row_order(ht)[[i]],])))
  out <- cbind(clu, paste("Cluster", i, sep=" "))
  colnames(out) <- c("GeneID", "TCGA Cluster")   } else {
    clu <- t(t(row.names(heatmap_data[row_order(ht)[[i]],])))
    clu <- cbind(clu, paste("Cluster", i, sep=" "))
    out <- rbind(out, clu)   } 
}
out

write.csv2(out, file = "TCGA_Gene_Cluster.csv")

