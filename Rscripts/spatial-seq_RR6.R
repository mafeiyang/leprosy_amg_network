library(Seurat)
library(ggplot2)
library(gplots)

# Load in the data
load("RR6_spatial-seq_seurat_object.Rda")

# Emulate ggplot default colors
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

# Set the colors to use
my_palette = gg_color_hue(6)[c(2,5,1,6,3)]
my_palette = c(my_palette, "grey60")

# Plot Figure 6b
SpatialPlot(spatial_rr6, image.alpha=0, pt.size.factor=2, stroke=0, cols=my_palette) +
  scale_y_reverse() + scale_x_reverse() + labs(fill=NULL) +
  guides(fill=guide_legend(override.aes=list(size=5)))

# Plot Figure 6c
SpatialFeaturePlot(spatial_rr6, features="CYBB", alpha=c(1, 1), pt.size.factor=1.8, min.cutoff="q70", 
                   max.cutoff="q90") + scale_fill_gradient(low="grey90", high=gg_color_hue(6)[6])
SpatialFeaturePlot(spatial_rr6, features="CCL5", alpha=c(1, 1), pt.size.factor=1.8, min.cutoff="q70", 
                   max.cutoff="q90") + scale_fill_gradient(low="grey90", high=gg_color_hue(6)[5])
SpatialFeaturePlot(spatial_rr6, features="MMP2", alpha=c(1, 1), pt.size.factor=1.8, min.cutoff="q70", 
                   max.cutoff="q90") + scale_fill_gradient(low="grey90", high=gg_color_hue(6)[3])
SpatialFeaturePlot(spatial_rr6, features="KLK5", alpha=c(1, 1), pt.size.factor=1.8, min.cutoff="q70", 
                   max.cutoff="q90") + scale_fill_gradient(low="grey90", high=gg_color_hue(6)[1])

# Plot Figure 6d
# Average expression by cluster for RR6
ave_expr_rr6 = read.csv("RR6_ave_expr_by_cluster.csv", row.names=1)
colnames(ave_expr_rr6) = gsub("X", "", colnames(ave_expr_rr6))
colnames(ave_expr_rr6) = paste(colnames(ave_expr_rr6), "RR6", sep="-")
ave_expr_rr6 = ave_expr_rr6[,1:5]
ave_expr_rr6 = t(scale(t(ave_expr_rr6)))
# Average expression by cluster for RR7
ave_expr_rr7 = read.csv("RR7_ave_expr_by_cluster.csv", row.names=1)
colnames(ave_expr_rr7) = gsub("X", "", colnames(ave_expr_rr7))
colnames(ave_expr_rr7) = paste(colnames(ave_expr_rr7), "RR7", sep="-")
ave_expr_rr7 = ave_expr_rr7[,c(2,4,1,5,3)]
ave_expr_rr7 = t(scale(t(ave_expr_rr7)))
# Average expression by cluster for RR8
ave_expr_rr8 = read.csv("RR8_ave_expr_by_cluster.csv", row.names=1)
colnames(ave_expr_rr8) = gsub("X", "", colnames(ave_expr_rr8))
colnames(ave_expr_rr8) = paste(colnames(ave_expr_rr8), "RR8", sep="-")
ave_expr_rr8 = ave_expr_rr8[,c(3,5,1,6,2)]
ave_expr_rr8 = t(scale(t(ave_expr_rr8)))
# Average expression by cluster for T-lep1
ave_expr_tlep1 = read.csv("Tlep1_ave_expr_by_cluster.csv", row.names=1)
colnames(ave_expr_tlep1) = gsub("X", "", colnames(ave_expr_tlep1))
colnames(ave_expr_tlep1) = paste(colnames(ave_expr_tlep1), "T-lep1", sep="-")
ave_expr_tlep1 = ave_expr_tlep1[,1:5]
ave_expr_tlep1 = t(scale(t(ave_expr_tlep1)))
# Combine the four datasets
genes = c("APP","CCL14","CCL21","CD40LG","CXCL13","IL32","CCL5","GZMB","PRF1","GNLY","DEFB1",
          "KRT6A","FLG2","KLK5","RNASE7","S100A7","CYBB","IDO1","MMP9","RNASE6","TXN",
          "CCL18","CYP27B1","CCL3","LYZ","ADM","CXCL12","MMP2","CXCL2","CXCL3","PLA2G2A")
ave_expr_rr6 = ave_expr_rr6[genes,]
ave_expr_rr7 = ave_expr_rr7[genes,]
ave_expr_rr8 = ave_expr_rr8[genes,]
ave_expr_tlep1 = ave_expr_tlep1[genes,]
ave_expr = cbind(ave_expr_rr6, ave_expr_rr7, ave_expr_rr8, ave_expr_tlep1)
ave_expr = ave_expr[,c(1,6,11,16,2,7,12,17,3,8,13,18,4,9,14,19,5,10,15,20)]
colnames(ave_expr) = paste(rep(c("EC","TC","KC","ML","FB"), each=4), rep(c("RR6","RR7","RR8","T-lep1"), 5), sep=" ")
# Plot the heatmap
my_palette = colorRampPalette(c("blue","white","red"))(n=100)
heatmap.2(as.matrix(t(ave_expr)), Rowv=FALSE, Colv=FALSE, dendrogram="none", scale="none",
          trace="none", key=TRUE, margins=c(10,10), adjCol=c(0.8,0.5),
          col=my_palette, key.title=NA, key.ylab=NA, cexCol=1, density.info="none", sepcolor="grey60",
          sepwidth=c(0.01,0.01), colsep=1:ncol(ave_expr), rowsep=1:nrow(ave_expr))
