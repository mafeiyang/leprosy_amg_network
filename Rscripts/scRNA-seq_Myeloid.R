library(Seurat)
library(scales)
library(ggplot2)
library(dplyr)

# Load the data
load("myeloid_scRNA-seq_seurat_object.Rda")

# Random order of the cells in the UMAP plot
set.seed(1)
random_order = sample(colnames(myeloid), ncol(myeloid))

# Colors for each origin
my_color = brewer_pal(palette="RdBu")(10)

# Plot Figure 3a
DimPlot(myeloid, label=TRUE, reduction="umap", group.by="ident", 
        pt.size=1, cells=random_order, label.size=4.3)

# Plot Figure 3b
DimPlot(myeloid, label=FALSE, reduction="umap", group.by="CF",
        pt.size=1, cells=random_order, cols=my_color[c(3,8)])

# Select 100 cells for each subtype to plot the marker genes in the heatmap
bars = c()
set.seed(1)
for (i in names(table(Idents(myeloid)))){
  new_bars = colnames(myeloid)[Idents(myeloid) == i]
  new_bars = new_bars[order(-myeloid$nCount_RNA[new_bars])][1:100]
  new_bars = sample(new_bars)
  bars = c(bars, new_bars)
}

# Select the top 100 marker genes for each cluster
cluster_markers = read.csv("myeloid_scRNA-seq_cluster_markers.csv", row.names=1, stringsAsFactors=FALSE)
top = cluster_markers %>% group_by(cluster) %>% top_n(100, avg_logFC)
genes = unique(as.character(top$gene))

# Plot Figure 3c
myeloid = ScaleData(object=myeloid, features=genes)
DoHeatmap(object=myeloid[,bars], features=genes, angle=0, size=4, raster=FALSE, disp.min=0, group.bar.height=0.01) +
  scale_fill_gradientn(colors=c("#FF00FF","000000","#FFFF00"), na.value='white') +
  guides(color=FALSE, fill=FALSE) + theme(axis.text.y=element_blank())

# Count number of cells for cluster for each origin
orig = myeloid$orig.ident
u_orig = sort(unique(orig))
cluster = as.character(Idents(myeloid))
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

# Plot Figure 3d
z1 = data.frame(ncells=as.vector(z),
                type=rep(as.integer(rownames(z)), ncol(z)),
                orig=rep(colnames(z), each=nrow(z)))
my_color = brewer_pal(palette="RdBu")(10)
my_color = my_color[c(1:6,8,10)]
z1$orig = gsub("[.]", "-", z1$orig)
ggplot(z1, aes(x=type, y=ncells, fill=orig, group=orig)) + geom_bar(stat="identity", width=0.7) + 
  scale_fill_manual(values=my_color) + theme_classic() +
  labs(x="Cluster", y="% Component of Samples", fill=NULL)+
  scale_x_continuous(breaks=0:max(z1$type))

# Plot Figure 3e
VlnPlot(myeloid, features=c("TREM2"), pt.size=0.1) + guides(fill=FALSE) +
  theme(axis.text.x=element_text(angle=0, hjust=0.5)) + labs(x=NULL)
VlnPlot(myeloid, features=c("APOE"), pt.size=0.1) + guides(fill=FALSE) +
  theme(axis.text.x=element_text(angle=0, hjust=0.5)) + labs(x=NULL)
genes = c("APOE","CTSB","TREM2","CD68","GPNMB","LGALS3","LPL","SPP1","TYROBP")
genes = intersect(genes, rownames(myeloid))
myeloid = AddModuleScore(myeloid, features=list(genes), name="TREM2_Score")
VlnPlot(myeloid, features=c("TREM2_Score1"), pt.size=0.1) + guides(fill=FALSE) +
  theme(axis.text.x=element_text(angle=0, hjust=0.5)) + 
  labs(x=NULL, y="Module Score", title="TREM2 Module Score")
