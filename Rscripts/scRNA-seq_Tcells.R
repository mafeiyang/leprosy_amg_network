library(Seurat)
library(scales)
library(ggplot2)
library(dplyr)

# Load the data
load("tcell_scRNA-seq_seurat_object.Rda")

# Random order of the cells in the UMAP plot
set.seed(1)
random_order = sample(colnames(tcells), ncol(tcells))

# Colors for each origin
my_color = brewer_pal(palette="RdBu")(10)

# Make Figure 2a
DimPlot(tcells, label=TRUE, reduction="umap", group.by="ident", 
        pt.size=0.8, cells=random_order, label.size=4.3)

# Select 100 cells for each subtype to plot the marker genes in the heatmap
bars = c()
set.seed(1)
for (i in names(table(Idents(tcells)))){
  new_bars = colnames(tcells)[Idents(tcells) == i]
  new_bars = new_bars[order(-tcells$nCount_RNA[new_bars])][1:100]
  new_bars = sample(new_bars)
  bars = c(bars, new_bars)
}

# Select the top 100 marker genes for each cluster
cluster_markers = read.csv("tcell_scRNA-seq_cluster_markers.csv", row.names=1, stringsAsFactors=FALSE)
top = cluster_markers %>% group_by(cluster) %>% top_n(100, avg_logFC)
genes = unique(as.character(top$gene))

# Plot Figure 2b
tcells = ScaleData(object=tcells, features=genes)
DoHeatmap(object=tcells[,bars], features=genes, angle=0, size=4, raster=FALSE, disp.min=0, group.bar.height=0.01) +
  scale_fill_gradientn(colors=c("#FF00FF","000000","#FFFF00"), na.value='white') +
  guides(color=FALSE, fill=FALSE) + theme(axis.text.y=element_blank())

# Plot Figure 2c
FeaturePlot(object=tcells, features=c("IFNG"), cols=c("grey90", "red3"), pt.size=0.5, order=TRUE,
            max.cutoff="q80", min.cutoff="q20") + theme(plot.title=element_text(size=20))
FeaturePlot(object=tcells, features=c("GZMB"), cols=c("grey90", "red3"), pt.size=0.5, order=TRUE,
            max.cutoff="q80", min.cutoff="q20") + theme(plot.title=element_text(size=20))
FeaturePlot(object=tcells, features=c("RORC"), cols=c("grey90", "red3"), pt.size=0.5, order=TRUE,
            max.cutoff="q80", min.cutoff="q20") + theme(plot.title=element_text(size=20))
FeaturePlot(object=tcells, features=c("PRF1"), cols=c("grey90", "red3"), pt.size=0.5, order=TRUE,
            max.cutoff="q80", min.cutoff="q20") + theme(plot.title=element_text(size=20))
FeaturePlot(object=tcells, features=c("FOXP3"), cols=c("grey90", "red3"), pt.size=0.5, order=TRUE,
            max.cutoff="q80", min.cutoff="q20") + theme(plot.title=element_text(size=20))
FeaturePlot(object=tcells, features=c("GNLY"), cols=c("grey90", "red3"), pt.size=0.5, order=TRUE,
            max.cutoff="q80", min.cutoff="q20") + theme(plot.title=element_text(size=20))

# Count number of cells for cluster for each origin
orig = tcells$orig.ident
u_orig = sort(unique(orig))
cluster = as.character(Idents(tcells))
uct = sort(unique(cluster))
z = matrix(nrow=length(uct), ncol=length(u_orig))
rownames(z) = uct
colnames(z) = u_orig
for (i in 1:nrow(z)){
  for (j in 1:ncol(z)){
    z[i,j] = sum(orig[cluster==uct[i]] == u_orig[j])
  }
}

# Normalize by the total number of cells for each sample and then scale by each cluster
z = t(t(z) / rowSums(t(z)))
z = z / rowSums(z) * 100

# Plot Figure 2d
z1 = data.frame(ncells=as.vector(z), type=rep(as.integer(rownames(z)), ncol(z)), orig=rep(colnames(z), each=nrow(z)))
my_color = brewer_pal(palette="RdBu")(10)
z1$orig = gsub("[.]", "-", z1$orig)
ggplot(z1, aes(x=type, y=ncells, fill=orig, group=orig)) +
  geom_bar(stat="identity", width=0.7) +
  scale_fill_manual(values=my_color) + 
  theme_classic() +
  labs(x="Cluster", y="% Component of Samples", fill=NULL)+
  scale_x_continuous(breaks=0:max(z1$type))
