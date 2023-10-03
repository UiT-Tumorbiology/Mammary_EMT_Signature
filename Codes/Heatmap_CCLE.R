#Import the libraries
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(ggplot2)
library(dendsort)
library(readxl)
#Import the data
expr <- read.csv2("CCLE_Scaled_Data.csv", as.is = T, check.names = F, row.names = 1) #Count matrix
annotation1 <- read.csv2("Sample_info.csv", as.is = T, check.names = F, row.names = 1) #Info about the cell state
annotation2 <- read.csv2("Supplementary_table_6.csv", as.is = T, check.names = F, row.names = 1) #Classification of genes (Up or Down)

#Subsetting the 265 genes from CCLE data
expr <- subset(expr, (rownames(expr) %in% rownames(annotation2)))
annotation2 <- annotation2[rownames(expr),]

annotation1 <- annotation1[colnames(expr),]

#Ensuring the orders are exactly the same
all(rownames(annotation2) == rownames(expr))
all(rownames(annotation1) == colnames(expr))



#Complexheatmap needs a matrix, so convert the dataframe to matrix
h <- as.matrix(expr)


#Reposition the breaks in quantile positions
quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}

h_breaks <- quantile_breaks(h, n = 11)
h_breaks



#Hierarchical Clustering of Genes and Cell lines
row_dist <- as.dist(1-cor(t(h), method = "pearson"))
col_dist <- as.dist(1-cor(h, method = "pearson"))
col_hc <- hclust(col_dist, method = "complete")
row_hc <- hclust(row_dist, method = "complete")

#Sorting the dendograms
Rowv=dendsort(as.dendrogram(row_hc), isRevers=TRUE, type = "average") 
Colv=dendsort(as.dendrogram(col_hc), type = "average")

#Defining the colors
col_fun <- colorRamp2(quantile_breaks(h, n = 11), c("#1e90ff", "#0479f5", "#0467d2", "#0356af", "#02458c","#023369","#012246", "#130202", "#6b0002","#b20004","#ff0000"))


#Annotation for the cell lines
ha <- HeatmapAnnotation(Cell_State = annotation1$Cell_State,
                        Subtype=annotation1$Subtype,
                        col=list(Subtype=c("BasalA"="#F22233",
                                           "BasalB"="#8C51A6",
                                           "HER2"="#3C5FA6",
                                           "Luminal"="sienna1"),
                                 Cell_State = c("Epithelial" = "lightslateblue",
                                                "Plastic" = "darkseagreen4",
                                                "Mesenchymal" = "darkred")),
                        simple_anno_size = unit(0.30, "cm"),
                        border = T,
                        annotation_name_gp= gpar(fontsize = 8.5,fontface="bold"),
                        show_legend = c(T,T),
                        annotation_legend_param = list(title_gp = gpar(fontsize=10, fontface="bold"), labels_gp=gpar(fontsize=8)))

#Annotation for the genes
hr <- rowAnnotation("Up/Down"=annotation2$`Up/Down regulated in cell lines`,
                    col=list("Up/Down"=c("Down" = "#1e90ff",
                                         "Up" = "#ff0000")),
                    simple_anno_size = unit(0.30, "cm"),
                    border = T,
                    annotation_name_gp= gpar(fontsize = 8.5,fontface="bold"),
                    show_legend = T,
                    annotation_legend_param = list(title_gp = gpar(fontsize=10, fontface="bold"), labels_gp=gpar(fontsize=8)))



ht <- Heatmap(h,col = col_fun,
              cluster_columns = Colv,
              name = "Expression Values",
              show_heatmap_legend = T ,
              top_annotation = ha,
              left_annotation = hr,
              show_column_names = T,
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
              row_split = 3, column_split = 3,
              column_gap = unit(c(0.5,0.5), "mm"),
              row_gap = unit(c(0.5,0.5), "mm"))

pdf("CCLE EMT Genes.pdf",  width=13,height=16)
draw(ht,  merge_legend=TRUE, padding = unit(c(2, 2, 2, 2), "mm"))
dev.off()


#Getting the genes from each column
ht <- draw(ht)

for (i in 1:length(row_order(ht))){   if (i == 1) {
  clu <- t(t(row.names(expr[row_order(ht)[[i]],])))
  out <- cbind(clu, paste("Cluster", i, sep=" "))
  colnames(out) <- c("GeneID", "CCLE_Cluster")   } else {
    clu <- t(t(row.names(expr[row_order(ht)[[i]],])))
    clu <- cbind(clu, paste("Cluster", i, sep=" "))
    out <- rbind(out, clu)   } 
}
out
#Merge with previous gene metadata file
x <- merge(out, annotation2, by.x = "GeneID", by.y = 0)
write.csv2(x, "Supplementary_table_6.csv")
