#Import the libraries
library(ComplexHeatmap)
library(RColorBrewer)
library(gplots)
library(circlize)
library(dplyr)
library(dendsort)

#Changing the work directory
setwd("~/Github")
#Import the data
expr <- read.csv("Genes_in_model.csv", sep = ";", as.is = T, row.names = 1) #Expression matrix
annotation <- read.csv("EMT_Model_Samples.csv", sep = ";", as.is = T, row.names = 1)


#Check if the matrix and sample files are in the same order
all(colnames(expr)==rownames(annotation))

#Log2 Normalize the data
expr <- log2(1+expr)


#Complexheatmap needs a matrix, so convert the dataframe to matrix
h <- as.matrix(expr)


#Determine the quantile breaks for using as color scale
quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}

h_breaks <- quantile_breaks(h, n = 11)
h_breaks

column_dend = as.dendrogram(hclust(dist(t(h)))) #To find the sample order and if they are separating based on epithelial and mesenchymal

plot(column_dend) #Visualize the dendogram

rm(col_dend)


#Choosing the color gradient
col_fun <- colorRamp2(quantile_breaks(h, n = 11), c("#1e90ff", "#0479f5", "#0467d2", "#0356af", "#02458c","#023369","#012246", "#130202", "#6b0002","#b20004","#ff0000"))


#Setting up the row annotation
ha <- HeatmapAnnotation("Cell State" = annotation$Cell_State,
                          col=list("Cell State" = c("Epithelial" = "lightslateblue",
                                                    "Mesenchymal" = "darkred")),
                        simple_anno_size = unit(0.30, "cm"),
                        border = T,
                        annotation_name_gp= gpar(fontsize = 8.5,fontface="bold"),
                        show_legend = T,
                        annotation_legend_param = list(title_gp = gpar(fontsize=10, fontface="bold"), labels_gp=gpar(fontsize=8)))

#Creating the heatmap object
ht <- Heatmap(h,col = col_fun,
              cluster_columns = FALSE,
              name = "Expression Values",
              show_heatmap_legend = T ,
              top_annotation = ha,
              show_column_names = T,
              show_row_names = F,cluster_rows = T, #Rows will be clustered
              row_title_gp = gpar(fontsize=12),
              column_title_gp = gpar(fontsize=12),
              height = unit(12, "cm"),
              width = unit(10.55, "cm"),
              show_row_dend = T,
              border = T,
              column_dend_height = unit(2, "cm"),
              column_names_rot = 45,
              row_title = NULL,
              column_names_gp = grid::gpar(fontsize = 7,fontface="bold"),
              heatmap_legend_param = list(title="Expression", legend_height=unit(3, "cm"),title_gp=gpar(fontsize=10, fontface="bold"),labels_gp = gpar(fontsize=8)),
              row_split = 5, column_split = 2,
              column_gap = unit(c(0.5,0.5), "mm"),
              row_gap = unit(c(0.5,0.5,0.5,0.5), "mm")) #setting the gap for each clusters

#Export the figure to a pdf file
pdf("EMT_Model.pdf",  width=12,height=16)
draw(ht,  merge_legend=TRUE, padding = unit(c(2, 2, 2, 2), "mm"))
dev.off()

