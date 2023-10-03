library(ComplexHeatmap)
library(dplyr)
library(dendextend)
library(circlize)
library(RColorBrewer)
library(readxl)
markers <- c("CDH1", "OCLN", "VIM", "FN1", "CDH2", "SNAI1",  "ZEB1", "ZEB2", "TWIST1", "GAPDH", "ACTB", "RPLP0", "B2M", "PPIA", "UBC",
         "TBP", "GUSB")
raw <- as.data.frame(read_xlsx("C:/Users/ssa214/UiT Office 365/O365-PhD_Saikat - General/bulk_EMT/Combined EMT Models/Heatmap on 3 cell lines/HMLE, D492, and MCF10a GRCh38.104 October 2021.xlsx"))

expr1 <- subset(raw,(Name %in% markers))
rownames(expr1) <- expr1$Name

expr1 <- expr1[match(markers, expr1$Name),]
expr1 <- expr1[,c(1,9,12,15)]

expr1$Name <- NULL
names(expr1) <- c("HMLE", "D492", "MCF10A")

#Taking logFC
heatmap_data <- log2(abs(expr1)) *sign(expr1) 

h <- as.matrix(heatmap_data)

#Reposition the breaks in quantile positions
quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}


col_fun <- colorRamp2(quantile_breaks(h, n = 11), c("#1e90ff", "#0356af", "#000000", "#000000", "#000000", 
                                                    "#000000", "#000000", "#b20004", "#b20004", "#F22233", "#ff0000"))

ht <- Heatmap(h,col = col_fun,
              cluster_columns = F,
              name = "Expression Values",
              show_heatmap_legend = T ,
              #top_annotation = ha,
              #left_annotation = hr,
              show_column_names = T,
              show_row_names = T,cluster_rows = F,
              row_title_gp = gpar(fontsize=11),
              column_title_gp = gpar(fontsize=11),
              height = unit(12, "cm"),
              width = unit(10.55, "cm"),
              #column_dend_reorder = T,
              #row_dend_reorder = T,
              #show_row_dend = T,
              border = T,
              #column_dend_height = unit(2, "cm"),
              column_names_rot = 0,
              column_dend_reorder = F,
              row_dend_reorder = F,
              #legend_height = unit(4, "cm"),
              row_names_gp = grid::gpar(fontsize = 10,fontface="italic"),
              column_names_gp = grid::gpar(fontsize = 11),
              heatmap_legend_param = list(title="Log2 Fold Change", legend_height=unit(3, "cm"),title_gp=gpar(fontsize=10, fontface="bold"),labels_gp = gpar(fontsize=8)),
              #column_split = 2,
              )
ht
pdf("EMT Model with markers.csv",  width=10,height=13)
draw(ht,  merge_legend=TRUE, padding = unit(c(2, 2, 2, 2), "mm"))
dev.off()
