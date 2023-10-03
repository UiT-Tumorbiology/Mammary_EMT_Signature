library(dplyr)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(dendsort)

expr <- read.csv2("EMT_matrix.csv", sep = ";", as.is = T, check.names = F, row.names = 1)
expr <- expr[,c(13,16,19)]

annotation <- read.csv("Supplementary_table_6.csv", sep=";", as.is = T, check.names = F)

up_an <- annotation[annotation$`Up/down regulated in EMT models` == "UP",]

expr <- subset(expr, (row.names(expr) %in% up_an$Name))

expr <- as.data.frame(t(scale(t(expr))))


h <- as.matrix(expr)


#Reposition the breaks in quantile positions
quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}

#Set the order of the dendogram
row_dist <- as.dist(1-cor(t(h), method = "pearson"))
col_dist <- as.dist(1-cor(h, method = "pearson"))
col_hc <- hclust(col_dist, method = "complete")
row_hc <- hclust(row_dist, method = "complete")
#Average
Rowv=dendsort(as.dendrogram(row_hc), isRevers=TRUE, type = "average") 
Colv=dendsort(as.dendrogram(col_hc), type = "average")

col_fun <- colorRamp2(quantile_breaks(h, n = 11), c("#1e90ff", "#0479f5", "#0467d2", "#0356af", "#02458c","#023369","#012246", "#130202", "#6b0002","#b20004","#ff0000"))


#Annotation
ht <- Heatmap(h,col = col_fun,
              cluster_columns = F,
              name = "Expression Values",
              show_heatmap_legend = T ,
              show_column_names = T,
              show_row_names = F,cluster_rows = Rowv,
              row_title_gp = gpar(fontsize=12),
              column_title_gp = gpar(fontsize=12),
              height = unit(12, "cm"),
              width = unit(10.55, "cm"),
              column_dend_reorder = T,
              row_dend_reorder = F,
              show_row_dend = T,
              border = T,
              column_dend_height = unit(2, "cm"),
              column_names_rot = 90,
              
              row_split = 3, 
              
              column_names_gp = grid::gpar(fontsize = 7,fontface="bold"),
              row_names_gp = grid::gpar(fontsize = 6),
              heatmap_legend_param = list(title="Expression", legend_height=unit(3, "cm"),title_gp=gpar(fontsize=10, fontface="bold"),labels_gp = gpar(fontsize=8)))

pdf("Heatmap_HMLER_Up.pdf", width=13,height=16)

draw(ht, merge_legend=T, padding = unit(c(2, 2, 2, 2), "mm"))
dev.off()
