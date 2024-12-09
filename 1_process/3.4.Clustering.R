# loading packages
library(Seurat)
library(harmony)
library(scales)
library(clustree)
library(clusterProfiler)
library(org.Mm.eg.db)
library(tidyverse)
library(dplyr)
library(Matrix)
library(glue)
library(ggplot2)
library(RColorBrewer)
library(grDevices)
library(ggthemes)
library(gridExtra)
library(patchwork)
library(cowplot)
library(ComplexHeatmap)
rm(list = ls());gc()

source("~/project/lung_10xST/src/utilities.R")

# Set work directory
projectDir = '~/project/lung_10xST'
setwd(projectDir)
stDir = file.path(projectDir, "results","E125_P0_merge8")
scDir = file.path(projectDir, "data/reference_sc/GSE165063_Development_2021/processing")
tgDir = file.path(projectDir, "doc", "Tangram_mapping")
outDir <- file.path(tgDir,"Clustering_celltype89_Deconvolution")

# load ST seurat object for 8 slices
st_ob <- readRDS(file.path(stDir, "seurat_Lung_E125_P0s1_8sections_merge_Log.rds"))
st_ob <- st_ob[!grepl("Gm42418", rownames(st_ob)),]
st_ob

# load the decomposition resluts with 89 celltypes from SC reference for origQC ST data
tg89 <- read.csv(file.path(tgDir, "Clustering_celltype89_Deconvolution/E125_P0s1_origQC_tangram_89CT_modeCell_trainGene_deconvolution_mtx.csv"), row.names = 1, check.names = FALSE)
tg89[1:3,1:5]

# create TangramNew89 assays
st_ob <- subset(st_ob, cells = colnames(tg89))                 
tg_assay89 <- CreateAssayObject(data = as.matrix(tg89))
st_ob[['TangramNew89']] <- tg_assay89
DefaultAssay(st_ob) <- "TangramNew89"

# Clustering
st_ob <- FindVariableFeatures(st_ob, assay = "TangramNew89", nfeatures = nrow(st_ob))
st_ob <- ScaleData(st_ob, features = rownames(st_ob))
st_ob <- RunPCA(st_ob, assay = "TangramNew89", npcs = 50, verbose = FALSE)
st_ob <- FindNeighbors(st_ob, reduction = 'pca', k.param = 50, dims = 1:20, verbose = FALSE)
resolutions <- c(seq(.5,2.0,.1))
st_ob <- FindClusters(st_ob, resolution = resolutions, algorithm=4, verbose = FALSE)
st_ob <- RunUMAP(st_ob, reduction = 'pca', 
                umap.method = 'umap-learn', metric = 'correlation', 
                n.neighbors = 30, 
                min.dist = 0.2, 
                dims = 1:20, 
                verbose = FALSE)
saveRDS(st_ob, file.path(outDir,"seurat_E125_P0s1_origQC_tangram_89CT_modeCell_trainGene_deconvolution.rds"))

# visualizaton
resolution = 1.1
image = Images(st_ob)
prefix = "TangramNew89_snn_res."
Idents(st_ob) = paste0(prefix, resolution)
dims = 20
cols2use = CustomCol(1:length(unique(Idents(st_ob))))
cols2use

groupby = paste0(prefix, resolution)
p1 <- DimPlot(st_ob, reduction = "umap", label = TRUE, group.by =  groupby, cols = cols2use, label.size = 6) + 
        labs (title = paste0('UMAP_res', resolution))  +
        theme(plot.title = element_text(hjust = 0.5))
p11 <- DimPlot(st_ob, reduction = "umap", group.by="Stage",label = TRUE,label.size = 4, repel = TRUE)
p1 | p11

groupby = "orig.ident"
propby = paste0(prefix, resolution)
p12 <- PlotAbundances(st_ob,  groupby, propby) + scale_fill_manual(values = cols2use)
p12

# SpatialDimPlot
slides = unique(st_ob$orig.ident)
spot.size = c(7,4,3,3,2.5,2,1.5,1.5)
names(spot.size) = slides
get_color <- function(slice_query) {
    a <- FetchData(st_ob, var = c(propby, 'orig.ident')) %>% 
        as.data.frame() %>% 
        filter(orig.ident == slice_query)
    b <- sort(unique(a[[propby]]))
    col_pat <- cols2use[b]
    return(col_pat)
}

gsp_lst = list()
for(slide in slides){
    gsp = SpatialDimPlot(st_ob, label = FALSE, images = slide,
                         pt.size.factor = as.numeric(spot.size[[slide]]), 
                         group.by = propby, repel = TRUE, combine = FALSE)   
    gsp = gsp[[1]] + scale_fill_manual(values = get_color(slide)) + 
          ggtitle(slide) + 
          theme(plot.title = element_text(size = 15, hjust = 0.5)) + NoLegend()     
    gsp_lst = c(gsp_lst, list(gsp)) 
}

gsp.all <- do.call(ggpubr::ggarrange, c(gsp_lst, list(ncol = 4, nrow = 2)))

cluster_cols <- data.frame(clusters = levels(Idents(st_ob)), cols2use) 
options(repr.plot.width = 18, repr.plot.height = 6)
lgd3 <- Legend(
    labels = cluster_cols$clusters,
    type = "points",
    size = unit(5, "mm"), 
    legend_gp = gpar(col = cluster_cols$cols2use), 
    background = "white",              
    labels_gp = gpar(fontsize = 15),
    grid_height = unit(6, "mm"), grid_width = unit(6, "mm"),
    gap = unit(1.6, "cm"),
    nrow = length(unique(Idents(st_ob))),
    row_gap = unit(1, "mm"),
    title = "Clusters",
    title_gp = gpar(fontsize = 16),
    title_gap = unit(5, "mm"),
    legend_height = unit(7, "cm")
)

options(repr.plot.width = 14, repr.plot.height= 8)
grid.newpage()
vp.sp <- viewport(x=0, y=0.5, width = 0.9, height = 0.9, just = c("left","center"))
pushViewport(vp.sp)
print(gsp.all, newpage = F)
upViewport()
vp.lgd3 <- viewport(x=0.9, y=0.9, width=0.1, height=0.9, just=c("left","center"))
pushViewport(vp.lgd3)
draw(lgd3, x = unit(0.5, "cm"), y = unit(0.09, "npc"), just = c("left","center"))
upViewport()

# Find marker genes for 10 Domains
Idents(st_ob) <- paste0(prefix, resolution)
DefaultAssay(st_ob) <- "Spatial"
# marker gene analysis
global_DEGs <- FindAllMarkers(st_ob, only.pos = T, test.use = "wilcox", logfc.threshold = 0.25, min.pct = 0.1)
dim(global_DEGs)

global_DEGs <- global_DEGs %>% mutate( gene_diff = round(global_DEGs$pct.1 / global_DEGs$pct.2, 3)) %>% dplyr::select( gene, everything())
filter_DEGs <- global_DEGs %>% filter( p_val < 0.05 & p_val_adj < 0.05 ) %>% dplyr::select( gene, everything())
top10_avg_logFC_markers <- filter_DEGs %>% group_by(cluster) %>% arrange(p_val,desc(avg_log2FC),desc(gene_diff)) %>% top_n(10,avg_log2FC) %>% arrange(cluster)
write.csv(global_DEGs, file.path(outDir, paste0("E125_P0s1_Tangram_celltype89_dim20_res1.2_cluster_wilcox_markers.csv")), row.names = FALSE)
write.csv(filter_DEGs, file.path(outDir, paste0("E125_P0s1_Tangram_celltype89_dim20_res1.2_cluster_wilcox_filtered_markers.csv")), row.names = FALSE)
write.csv(top10_avg_logFC_markers, file.path(outDir, paste0("E125_P0s1_Tangram_celltype89_dim20_res1.2_cluster_wilcox_filtered_top10_markers.csv")), row.names = FALSE)
# heatmap
topn_markers <- top10_avg_logFC_markers %>% group_by(cluster) %>% top_n(10, gene_diff)
topn_markers2vis <- unique(as.vector(topn_markers$gene))
colors2use <- CustomCol(1:length(unique(st_ob@meta.data[,propby])))
ggheatmap <- DoHeatmap( object = st_ob,
                           features = topn_markers2vis,
                           group.colors = colors2use,
                           group.by = levels(propby), 
                           group.bar = T, label = T, 
                           size = 4, hjust = 0.5, angle = 0) + 
                theme(axis.text.y = element_text(size = 10 * 0.9, face = "bold"),
                      legend.text = element_text(size = 6),
                      legend.title = element_text(size = 7),
                      legend.key.size = unit(0.8, "lines")) + 
                guides(fill = guide_colorbar( title.position = "top")) + 
                       scale_fill_viridis(option = "D")
ggheatmap
ggsave(file.path(outDir, paste0("E125_P0s1_Tangram_celltype89_dim20_res1.2_cluster_wilcox_filtered_top10_markers_heatmap.pdf")), plot = ggheatmap, width = 10, height = 12)