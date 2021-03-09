library(Seurat)
library(ggplot2)
library(scales)

# Load the data
load("all_cells_scRNA-seq_seurat_object.Rda")

# Random order of the cells in the UMAP plot
set.seed(1)
random_order = sample(colnames(allcells), ncol(allcells))

# Colors for each origin
my_color = brewer_pal(palette="RdBu")(10)

# Plot Figure 1a
DimPlot(allcells, label=TRUE, reduction="umap", group.by="celltype", 
        pt.size=0.5, cells=random_order, label.size=4.3)

# Plot Figure 1b
DimPlot(allcells, label=FALSE, reduction="umap", group.by="CF",
        pt.size=0.5, cells=random_order, cols=my_color[c(3,8)])

# Select 200 cells for each cell type to plot the marker genes in the heatmap
table(allcells$celltype)
bars = c()
set.seed(1)
for (i in as.character(unique(allcells$celltype))){
  new_bars = colnames(allcells)[allcells$celltype == i]
  new_bars = new_bars[order(-allcells$nCount_RNA[new_bars])][1:200]
  new_bars = sample(new_bars)
  bars = c(bars, new_bars)
}

# Three signature genes for each cell type
genes = c("CD3D","CD3E","TRBC2","MS4A1","BANK1","CD79A","IGHG3","IGHG1","IGHG4",
          "LYZ","C1QC","CXCL8","CCL22","CD1A","CD207","TPSB2","CPA3","CTSG",
          "KRT10","KRT2","KRT1","C3","DCN","COL1A1","MYH11","ACTA2","TAGLN",
          "PECAM1","CDH5","AQP1","DCD","MUCL1","SCGB2A2","TYRP1","DCT","PMEL")

# Plot Figure 1c
heatmap_object = subset(allcells, cells=bars)
Idents(heatmap_object) = heatmap_object$celltype
heatmap_object = ScaleData(object=heatmap_object, features=rownames(heatmap_object))
DoHeatmap(object=heatmap_object, features=genes, angle=45, size=4, raster=FALSE, 
          cells=bars, disp.min=0, group.bar.height=0.01) +
  scale_fill_gradientn(colors=c("#FF00FF","white","#FFFF00"), na.value='white') + 
  guides(color=FALSE)

# Count the number of cells for each cell types by each sample
orig = allcells$orig.ident
u_orig = sort(unique(orig))
celltype = as.character(allcells$celltype)
uct = sort(unique(celltype))
z = matrix(nrow=length(uct), ncol=length(u_orig))
rownames(z) = uct
colnames(z) = u_orig
for (i in 1:nrow(z)){
  for (j in 1:ncol(z)){
    z[i,j] = sum(orig[celltype==uct[i]] == u_orig[j])
  }
}

# Normalize by the total number of cells for each sample and then scale by each cell type
z = t(t(z) / rowSums(t(z)))
z = z / rowSums(z) * 100

# Plot Figure 1d
z1 = data.frame(ncells=as.vector(z), type=rep(rownames(z), ncol(z)), orig=rep(colnames(z), each=nrow(z)))
celltype_order = c("T cell","B cell","Plasma cell","Myeloid","Langerhans","Mast cell","Keratinocyte",
                   "Fibroblast","Smooth muscle","Endothelial","Eccrine gland","Melanocyte")
z1$type = factor(z1$type, levels=celltype_order)
z1$orig = gsub("[.]", "-", z1$orig)
my_color = brewer_pal(palette="RdBu")(10)
ggplot(z1, aes(x=type, y=ncells, fill=orig, group=orig)) +
  theme_classic() + geom_bar(stat="identity", width=0.7) +
  scale_fill_manual(values=my_color) + 
  labs(x="", y="% Component of allcellss", fill=NULL) +
  theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1, size=12))

# Plot Figure 1e
VlnPlot(allcells, features=c("MLEP_nUMIs"), group.by="orig.ident", pt.size=0.1) +
  labs(y="log10 ( Mlep nUMIs )", title="MLEP nUMIs") + scale_fill_manual(values=my_color) +
  labs(x=NULL) + guides(fill=FALSE)
