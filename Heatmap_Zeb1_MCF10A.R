library(dplyr)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(dendsort)

expr <- read.csv2("MCF10A_replicates_with_EMT_expressions.csv",as.is = T, check.names = F, row.names = 1)
expr <- expr[complete.cases(expr),]
expr <- as.data.frame(t(scale(t(expr))))

annotation <- read.csv2("Supplementary_table_6.csv", as.is = T, check.names = F)
annotation <- subset(annotation, (Name %in% rownames(expr)))

expr <- expr[annotation$Name,]
all(rownames(expr) == annotation$Name)


h <- as.matrix(expr)
mean(h)

h_breaks <- seq(min(h), max(h), length.out = 10)
h_breaks

#Reposition the breaks in quantile positions
quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}


#Set the order of the dendogram
library(dendsort)
row_dist <- as.dist(1-cor(t(h), method = "pearson"))
col_dist <- as.dist(1-cor(h, method = "pearson"))
col_hc <- hclust(col_dist, method = "complete")
row_hc <- hclust(row_dist, method = "complete")
#Average
Rowv=dendsort(as.dendrogram(row_hc), isRevers=TRUE, type = "average") 
Colv=dendsort(as.dendrogram(col_hc), type = "average")

col_fun <- colorRamp2(quantile_breaks(h, n = 11), c("#1e90ff", "#0479f5", "#0467d2", "#0356af", "#02458c","#023369","#012246", "#130202", "#6b0002","#b20004","#ff0000"))


#Annotation

hr <- rowAnnotation("Up/Down"=annotation$`Up/down regulated in EMT models`,
                    "Gene Cluster - Cell Lines" = annotation$`Gene cluster, breast cancer cell lines`,
                    "Gene Cluster - Patients" = annotation$`Gene cluster, TCGA breast cancer patients`,
                    col=list("Gene Cluster - Cell Lines"=c("Cluster A"="cadetblue2",
                                                          "Cluster B"="burlywood2",
                                                           "Cluster C"="chocolate2",
                                                          "#N/A" = "#CCCCCC"),
                             "Gene Cluster - Patients"=c("Cluster A"="#66C2A5",
                                                           "Cluster B"="#FC8D62",
                                                           "Cluster C"="#8DA0CB"),
                             "Up/Down"=c("DOWN" = "#1e90ff",
                                         "UP" = "#ff0000")),
                    simple_anno_size = unit(0.30, "cm"),
                    border = T,
                    annotation_name_gp= gpar(fontsize = 8.5,fontface="bold"),
                    show_legend = T,
                    annotation_legend_param = list(title_gp = gpar(fontsize=10, fontface="bold"), labels_gp=gpar(fontsize=8)))

col.subsections <- c(4,4,4,4) #As we want the heatmap to be clustered into four sample clusters
ht <- Heatmap(h,col = col_fun,
              cluster_columns = F,
              name = "Expression Values",
              show_heatmap_legend = T ,
              #top_annotation = ha,
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
              #legend_height = unit(4, "cm"),
              row_split = 5, 
              column_split = data.frame(rep(c("A", "B", "C", "D"), col.subsections)),
              column_names_gp = grid::gpar(fontsize = 7,fontface="bold"),
              row_names_gp = grid::gpar(fontsize = 6),
              heatmap_legend_param = list(title="Expression", legend_height=unit(3, "cm"),title_gp=gpar(fontsize=10, fontface="bold"),labels_gp = gpar(fontsize=8)))

pdf("Heatmap_ZEB1KO.pdf", width=13,height=16)
draw(ht, merge_legend=T, padding = unit(c(2, 2, 2, 2), "mm"))
dev.off()
ht <- draw(ht)


for (i in 1:length(row_order(ht))){   if (i == 1) {
  clu <- t(t(row.names(heatmap_data[row_order(ht)[[i]],])))
  out <- cbind(clu, paste("cluster", i, sep=""))
  colnames(out) <- c("GeneID", "Cluster")   } else {
    clu <- t(t(row.names(heatmap_data[row_order(ht)[[i]],])))
    clu <- cbind(clu, paste("cluster", i, sep=""))
    out <- rbind(out, clu)   } 
}
out

write.csv2(out, file = "Genes_cluster_MCF10A.csv")


