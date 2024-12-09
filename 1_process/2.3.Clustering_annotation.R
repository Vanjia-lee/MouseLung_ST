library(Seurat)
library(future)
library(future.apply)
library(stringr)
library(ggplot2)
library(RColorBrewer)

projectDir <- "~/project/lung_10xST/"
sctdir <- file.path(projectDir, "results/4.Development_lineage_construction/single_stage_SCT_cluster")
outdir <- file.path(sctdir, "TOME_result/harmony_SCT")
slides = c("E125_S2","E135_S3","E145_S4","E155_S4","E165_S2","E175_S2","E185_S3","P0_S1")
time_point <-  c(paste0("E",seq(125,185,10)), "P0")
resols <- c(1.2, 1.2, 1.2, 1.3, 1.9, 1.8, 1.2, 1.9) 
names(resols) <- slides
spot.size <- c(7,4,3,3,2.5,2,2,2)
names(spot.size) <- slides

#8
fib.group <- c("Fibro_1", "Fibro_2", "Fibro_3", "Fibro_1&2", "Fibro_4", "Fibro_5", "Fibro_6")
#10
mes.group <- c("Mes_1", "Mes_2", "Mes_3", "Mes_4", "Mes_5", "Mes_6", "Mes_7", "Mes_8", "Mes_9")
#12
dist.group <- c("Distal_1","Alv_1","Alv_2","Alv_3","Alv_4","Alv_5","Alv_6")
#6
prox.group <- c("Proximal_1", "Proximal_2","Proximal_3","Proximal_4","Proximal_5","Proximal_6")
#1
immune.group <- c("Immune")

ct.group <- c(fib.group,mes.group,dist.group,prox.group,immune.group)
length(ct.group)

#cols1 <- RColorBrewer::brewer.pal(n = 9, name = "Set1")
cols2 <- RColorBrewer::brewer.pal(n = 8, name = "Dark2")
cols3 <- RColorBrewer::brewer.pal(n = 12, name = "Paired")
cols4 <- RColorBrewer::brewer.pal(n = 12, name = "Set3")
cols2use <- sample(c(cols4[-c(2,9)],cols3,cols2[-8],"#8CA1CC"))[1:length(ct.group)]
names(cols2use) = ct.group
cols2use

plotDir <- file.path(sctdir, "Dimplot_new_annotation")
if(!dir.exists(plotDir)) dir.create(plotDir, recursive = T)

# E125
slide1 <- slides[1]
res1 <- resols[1]
anno1 <- readRDS(paste0(sctdir, "/", slide1, "/seurat_", slide1, "_SCT_res", res1, "_dim30_cluster.rds"))
anno1$cell_state2 = "other"
anno1$cell_state2[anno1$SCT_snn_res.1.2 == 1] = "Fibro_1" 
anno1$cell_state2[anno1$SCT_snn_res.1.2 == 2] = "Distal_1" 
anno1$cell_state2[anno1$SCT_snn_res.1.2 == 3] = "Mes_1" 
anno1$cell_state2[anno1$SCT_snn_res.1.2 == 4] = "Proximal_1" 
table(anno1$cell_state2)
saveRDS(anno1,paste0(plotDir, "/seurat_", slide1, "_SCT_res", res1, "_dim30_cluster.rds"))

cols.slice = cols2use[unique(anno1$cell_state2)]
dp = DimPlot(anno1, reduction = "umap", group.by = "cell_state2",label = FALSE, cols = cols.slice, label.size = 4) +
         labs (title = paste0(slide1,'_UMAP_res', res1)) + theme(plot.title=element_text(hjust=0.5))
gsp = SpatialDimPlot(anno1, images = slide1, group.by = "cell_state2", label = FALSE, label.size = 6,  image.alpha = 0.5,
                     pt.size.factor = spot.size[slide1],repel = TRUE) + scale_fill_manual(values = cols.slice) + NoLegend()
p1 = dp + gsp
ggsave(file.path(plotDir, paste0(slide1, "_DimPlot.pdf")), p1, width = 10, height = 5)

# E135
slide1 <- slides[2]
res1 <- resols[2]
anno1 <- readRDS(paste0(sctdir, "/", slide1, "/seurat_", slide1, "_SCT_res", res1, "_dim30_cluster.rds"))
anno1$cell_state2 = "other"
anno1$cell_state2[anno1$SCT_snn_res.1.2 == 1] = "Mes_1"
anno1$cell_state2[anno1$SCT_snn_res.1.2 == 2] = "Fibro_1"  
anno1$cell_state2[anno1$SCT_snn_res.1.2 == 3] = "Distal_1"
anno1$cell_state2[anno1$SCT_snn_res.1.2 == 4] = "Fibro_2" 
anno1$cell_state2[anno1$SCT_snn_res.1.2 == 5] = "Proximal_1"
table(anno1$cell_state2)
saveRDS(anno1,paste0(plotDir, "/seurat_", slide1, "_SCT_res", res1, "_dim30_cluster.rds"))
cols.slice = cols2use[unique(anno1$cell_state2)]
dp = DimPlot(anno1, reduction = "umap", group.by = "cell_state2",label = FALSE, cols = cols.slice, label.size = 4) +
         labs (title = paste0(slide1,'_UMAP_res', res1)) + theme(plot.title=element_text(hjust=0.5))
gsp = SpatialDimPlot(anno1, images = slide1, group.by = "cell_state2", label = FALSE, label.size = 6,  image.alpha = 0.5,
                     pt.size.factor = spot.size[slide1],repel = TRUE) + scale_fill_manual(values = cols.slice) + NoLegend()
p1 = dp + gsp
options(repr.plot.width = 10, repr.plot.height = 4)
p1

ggsave(file.path(plotDir, paste0(slide1, "_DimPlot.pdf")), p1, width = 10, height = 5)

# E145
slide1 <- slides[3]
res1 <- resols[3]
anno1 <- readRDS(paste0(sctdir, "/", slide1, "/seurat_", slide1, "_SCT_res", res1, "_dim30_cluster.rds"))
anno1$cell_state2 = "other"
anno1$cell_state2[anno1$SCT_snn_res.1.2 == 1] = "Mes_3" 
anno1$cell_state2[anno1$SCT_snn_res.1.2 == 2] = "Distal_1" 
anno1$cell_state2[anno1$SCT_snn_res.1.2 == 3] = "Mes_2"
anno1$cell_state2[anno1$SCT_snn_res.1.2 == 4] = "Proximal_2"
anno1$cell_state2[anno1$SCT_snn_res.1.2 == 5] = "Fibro_2"
anno1$cell_state2[anno1$SCT_snn_res.1.2 == 6] = "Fibro_3"
anno1$cell_state2[anno1$SCT_snn_res.1.2 == 7] = "Fibro_1"
table(anno1$cell_state2)
saveRDS(anno1,paste0(plotDir, "/seurat_", slide1, "_SCT_res", res1, "_dim30_cluster.rds"))

cols.slice = cols2use[unique(anno1$cell_state2)]
dp = DimPlot(anno1, reduction = "umap", group.by = "cell_state2",label = FALSE, cols = cols.slice, label.size = 4) +
         labs (title = paste0(slide1,'_UMAP_res', res1)) + theme(plot.title=element_text(hjust=0.5))
gsp = SpatialDimPlot(anno1, images = slide1, group.by = "cell_state2", label = FALSE, label.size = 6,  image.alpha = 0.5,
                     pt.size.factor = spot.size[slide1],repel = TRUE) + scale_fill_manual(values = cols.slice) + NoLegend()
p1 = dp + gsp
options(repr.plot.width = 10, repr.plot.height = 4)
p1

ggsave(file.path(plotDir, paste0(slide1, "_DimPlot.pdf")), p1, width = 10, height = 5)

# E155
slide1 <- slides[4]
res1 <- resols[4]
anno1 <- readRDS(paste0(sctdir, "/", slide1, "/seurat_", slide1, "_SCT_res", res1, "_dim30_cluster.rds"))
table(anno1$SCT_snn_res.1.3)
anno1$cell_state2 = "other"
anno1$cell_state2[anno1$SCT_snn_res.1.3 == 1] = "Proximal_2" 
anno1$cell_state2[anno1$SCT_snn_res.1.3 == 2] = "Mes_3" 
anno1$cell_state2[anno1$SCT_snn_res.1.3 == 3] = "Fibro_1&2"
anno1$cell_state2[anno1$SCT_snn_res.1.3 == 4] = "Distal_1" 
anno1$cell_state2[anno1$SCT_snn_res.1.3 == 5] = "Mes_2"
anno1$cell_state2[anno1$SCT_snn_res.1.3 == 6] = "Fibro_3"
anno1$cell_state2[anno1$SCT_snn_res.1.3 == 7] = "Distal_1"
table(anno1$cell_state2)
saveRDS(anno1,paste0(plotDir, "/seurat_", slide1, "_SCT_res", res1, "_dim30_cluster.rds"))

cols.slice = cols2use[unique(anno1$cell_state2)]
dp = DimPlot(anno1, reduction = "umap", group.by = "cell_state2",label = FALSE, cols = cols.slice, label.size = 4) +
         labs (title = paste0(slide1,'_UMAP_res', res1)) + theme(plot.title=element_text(hjust=0.5))
gsp = SpatialDimPlot(anno1, images = slide1, group.by = "cell_state2", label = FALSE, label.size = 6,  image.alpha = 0.5,
                     pt.size.factor = spot.size[slide1],repel = TRUE) + scale_fill_manual(values = cols.slice) + NoLegend()
p1 = dp + gsp
options(repr.plot.width = 10, repr.plot.height = 4)
p1

ggsave(file.path(plotDir, paste0(slide1, "_DimPlot.pdf")), p1, width = 10, height = 5)

# E165
slide1 <- slides[5]
res1 <- resols[5]
anno1 <- readRDS(paste0(sctdir, "/", slide1, "/seurat_", slide1, "_SCT_res", res1, "_dim30_cluster.rds"))
table(anno1$SCT_snn_res.1.9)
anno1$cell_state2 = "other"
anno1$cell_state2[anno1$SCT_snn_res.1.9 == 1] = "Alv_1" 
anno1$cell_state2[anno1$SCT_snn_res.1.9 == 2] = "Alv_3" 
anno1$cell_state2[anno1$SCT_snn_res.1.9 == 3] = "Alv_2"
anno1$cell_state2[anno1$SCT_snn_res.1.9 == 4] = "Mes_2" 
anno1$cell_state2[anno1$SCT_snn_res.1.9 == 5] = "Mes_4"
anno1$cell_state2[anno1$SCT_snn_res.1.9 == 6] = "Proximal_2"
anno1$cell_state2[anno1$SCT_snn_res.1.9 == 7] = "Fibro_2"
anno1$cell_state2[anno1$SCT_snn_res.1.9 == 8] = "Mes_5"
anno1$cell_state2[anno1$SCT_snn_res.1.9 == 9] = "Fibro_1"
table(anno1$cell_state2)
saveRDS(anno1,paste0(plotDir, "/seurat_", slide1, "_SCT_res", 1.9, "_dim30_cluster.rds"))

cols.slice = cols2use[unique(anno1$cell_state2)]
dp = DimPlot(anno1, reduction = "umap", group.by = "cell_state2",label = FALSE, cols = cols.slice, label.size = 4) +
         labs (title = paste0(slide1,'_UMAP_res', res1)) + theme(plot.title=element_text(hjust=0.5))
gsp = SpatialDimPlot(anno1, images = slide1, group.by = "cell_state2", label = FALSE, label.size = 6,  image.alpha = 0.5,
                     pt.size.factor = spot.size[slide1],repel = TRUE) + scale_fill_manual(values = cols.slice) + NoLegend()
p1 = dp + gsp
options(repr.plot.width = 10, repr.plot.height = 4)
p1

ggsave(file.path(plotDir, paste0(slide1, "_DimPlot.pdf")), p1, width = 10, height = 5)

# E175
slide1 <- slides[6]
res1 <- resols[6]
anno1 <- readRDS(paste0(sctdir, "/", slide1, "/seurat_", slide1, "_SCT_res", res1, "_dim30_cluster.rds"))
table(anno1$SCT_snn_res.1.8)
anno1$cell_state2 = "other"
anno1$cell_state2[anno1$SCT_snn_res.1.8 == 1] = "Proximal_3" 
anno1$cell_state2[anno1$SCT_snn_res.1.8 == 2] = "Mes_4" 
anno1$cell_state2[anno1$SCT_snn_res.1.8 == 3] = "Alv_3"
anno1$cell_state2[anno1$SCT_snn_res.1.8 == 4] = "Proximal_4" 
anno1$cell_state2[anno1$SCT_snn_res.1.8 == 5] = "Alv_4"
anno1$cell_state2[anno1$SCT_snn_res.1.8 == 6] = "Mes_6"
anno1$cell_state2[anno1$SCT_snn_res.1.8 == 7] = "Mes_7"
anno1$cell_state2[anno1$SCT_snn_res.1.8 == 8] = "Fibro_5"
anno1$cell_state2[anno1$SCT_snn_res.1.8 == 9] = "Immune"
anno1$cell_state2[anno1$SCT_snn_res.1.8 == 10] = "Fibro_4"
table(anno1$cell_state2)
saveRDS(anno1,paste0(plotDir, "/seurat_", slide1, "_SCT_res", res1, "_dim30_cluster.rds"))

cols.slice = cols2use[unique(anno1$cell_state2)]
dp = DimPlot(anno1, reduction = "umap", group.by = "cell_state2",label = FALSE, cols = cols.slice, label.size = 4) +
         labs (title = paste0(slide1,'_UMAP_res', res1)) + theme(plot.title=element_text(hjust=0.5))
gsp = SpatialDimPlot(anno1, images = slide1, group.by = "cell_state2", label = FALSE, label.size = 6,  image.alpha = 0.5,
                     pt.size.factor = spot.size[slide1],repel = TRUE) + scale_fill_manual(values = cols.slice) + NoLegend()
p1 = dp + gsp
options(repr.plot.width = 10, repr.plot.height = 4)
p1

ggsave(file.path(plotDir, paste0(slide1, "_DimPlot.pdf")), p1, width = 10, height = 5)

# E185
slide1 <- slides[7]
res1 <- resols[7]
anno1 <- readRDS(paste0(sctdir, "/", slide1, "/seurat_", slide1, "_SCT_res", res1, "_dim30_cluster.rds"))
table(anno1$SCT_snn_res.1.9)
anno1$cell_state2 = "other"
anno1$cell_state2[anno1$SCT_snn_res.1.2 == 1] = "Alv_3" 
anno1$cell_state2[anno1$SCT_snn_res.1.2 == 2] = "Alv_4" 
anno1$cell_state2[anno1$SCT_snn_res.1.2 == 3] = "Fibro_1"
anno1$cell_state2[anno1$SCT_snn_res.1.2 == 4] = "Mes_4" 
anno1$cell_state2[anno1$SCT_snn_res.1.2 == 5] = "Proximal_3"
anno1$cell_state2[anno1$SCT_snn_res.1.2 == 6] = "Mes_9"
anno1$cell_state2[anno1$SCT_snn_res.1.2 == 7] = "Proximal_4"
anno1$cell_state2[anno1$SCT_snn_res.1.2 == 8] = "Mes_8"
anno1$cell_state2[anno1$SCT_snn_res.1.2 == 9] = "Fibro_6"
anno1$cell_state2[anno1$SCT_snn_res.1.2 == 10] = "Fibro_6"
anno1$cell_state2[anno1$SCT_snn_res.1.2 == 11] = "Proximal_5"
anno1$cell_state2[anno1$SCT_snn_res.1.2 == 12] = "Immune"
table(anno1$cell_state2)
saveRDS(anno1,paste0(plotDir,"/seurat_", slide1, "_SCT_res", res1, "_dim30_cluster.rds"))

cols.slice = cols2use[unique(anno1$cell_state2)]
dp = DimPlot(anno1, reduction = "umap", group.by = "cell_state2",label = FALSE, cols = cols.slice, label.size = 4) +
         labs (title = paste0(slide1,'_UMAP_res', res1)) + theme(plot.title=element_text(hjust=0.5))
gsp = SpatialDimPlot(anno1, images = slide1, group.by = "cell_state2", label = FALSE, label.size = 6,  image.alpha = 0.5,
                     pt.size.factor = spot.size[slide1],repel = TRUE) + scale_fill_manual(values = cols.slice) + NoLegend()
p1 = dp + gsp
options(repr.plot.width = 10, repr.plot.height = 4)
p1

ggsave(file.path(plotDir, paste0(slide1, "_DimPlot.pdf")), p1, width = 10, height = 5)

# P0
slide1 <- slides[8]
res1 <- resols[8]
anno1 <- readRDS(paste0(sctdir, "/", slide1, "/seurat_", slide1, "_SCT_res", res1, "_dim30_cluster.rds"))
table(anno1$SCT_snn_res.1.2)
anno1$cell_state2 = "other"
anno1$cell_state2[anno1$SCT_snn_res.1.9 == 1] = "Alv_4" 
anno1$cell_state2[anno1$SCT_snn_res.1.9 == 2] = "Alv_5" 
anno1$cell_state2[anno1$SCT_snn_res.1.9 == 3] = "Alv_6"
anno1$cell_state2[anno1$SCT_snn_res.1.9 == 4] = "Alv_3" 
anno1$cell_state2[anno1$SCT_snn_res.1.9 == 5] = "Proximal_3"
anno1$cell_state2[anno1$SCT_snn_res.1.9 == 6] = "Immune"
anno1$cell_state2[anno1$SCT_snn_res.1.9 == 7] = "Fibro_2"
anno1$cell_state2[anno1$SCT_snn_res.1.9 == 8] = "Proximal_6"
anno1$cell_state2[anno1$SCT_snn_res.1.9 == 9] = "Proximal_4"
anno1$cell_state2[anno1$SCT_snn_res.1.9 == 10] = "Fibro_1"
anno1$cell_state2[anno1$SCT_snn_res.1.9 == 11] = "Fibro_6"
anno1$cell_state2[anno1$SCT_snn_res.1.9 == 12] = "Mes_9"
table(anno1$cell_state2)
saveRDS(anno1,paste0(plotDir,  "/seurat_", slide1, "_SCT_res", res1, "_dim30_cluster.rds"))

cols.slice = cols2use[unique(anno1$cell_state2)]
dp = DimPlot(anno1, reduction = "umap", group.by = "cell_state2",label = FALSE, cols = cols.slice, label.size = 4) +
         labs (title = paste0(slide1,'_UMAP_res', res1)) + theme(plot.title=element_text(hjust=0.5))
gsp = SpatialDimPlot(anno1, images = slide1, group.by = "cell_state2", label = FALSE, label.size = 6,  image.alpha = 0.5,
                     pt.size.factor = spot.size[slide1],repel = TRUE) + scale_fill_manual(values = cols.slice) + NoLegend()
p1 = dp + gsp
options(repr.plot.width = 10, repr.plot.height = 4)
p1

ggsave(file.path(plotDir, paste0(slide1, "_DimPlot.pdf")), p1, width = 10, height = 5)