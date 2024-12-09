library(Seurat)
library(tidyverse)
library(patchwork)
library(openxlsx)
library(pheatmap)
rm(list=ls())

source("~/project/lung_10xST/src/0.utilities/utilities.R")

# set work directory
projectDir = '~/project/lung_10xST'
tgDir = file.path(projectDir, "results", "3.2.Cluster_deconvolution_Tangram","tangram_origQC_training_genes_mode_cells/Clustering_celltype89_Deconvolution")
miaDir <- file.path(outDir, "MIA_results")
if(!dir.exists(miaDir)) dir.create(miaDir, recursive = T)

# load ST seurat object
st_ob <- readRDS(file.path(tgDir,"seurat_E125_P0s1_origQC_tangram_89CT_modeCell_trainGene_deconvolution_res1.1.rds"))
Idents(st_ob) = 'TangramNew89_snn_res.1.1'
st27 <- subset(st_ob, subset = TangramNew89_snn_res.1.1 %in% c(2,7))

# load SC data
sc_ob <- readRDS(file.path(projectDir, "data/reference_sc/GSE165063_Development_2021/processing/SC_ref_E12_P0_new_subcelltype.rds"))
dim(sc_ob)
sc.degs <- read.xlsx(file.path(projectDir, "data/reference_sc/GSE165063_Development_2021/processing/DEG_clusters.xlsx"), sheet = "Epi - AT2", rowNames = TRUE)
sc_ob@meta.data$celltype[sc_ob@meta.data$celltype == "Ciliated"] = "Cilliated"

epi <- c("Early Epithelium", "Transitional", "AT1", "AT2", "Mki67+ AT2", "Cilliated", "Secretory", "Neuroendocrine")
endo <- c("Arterial maEC","Venous maEC","Prolif. gCap","gCap", "aCap", "Lymphatic")
meso <- c( "Prolif. Wnt2+ FB", "Wnt2+ FB", "Prolif. Myo FB", "Myofibroblast", "Adventitial FB", "Pericyte", "Mesothelium", "Smooth Muscle", "Cardiomyocyte")
degs_lst <- list()
for(ct in unique(sc_ob@meta.data$celltype)){
    prefix = ifelse(ct %in% epi, "Epi", ifelse(ct %in% endo, "Endo", "Meso"))
    sc.degs <- read.xlsx(file.path(projectDir, "data/reference_sc/GSE165063_Development_2021/processing/DEG_clusters.xlsx"), sheet = paste(prefix, "-", ct), rowNames = TRUE)
    sc.degs$Celltype = ct
    sc.degs$gene = rownames(sc.degs)
    degs_lst <- c(degs_lst, list(sc.degs))
}
names(degs_lst) <- unique(sc_ob@meta.data$celltype)
sc.degs <- do.call(rbind,degs_lst)
colnames(sc.degs)[6] = "cluster"
sc.degs.filter <- sc.degs %>% filter(p_val_adj <= 0.05) %>% group_by(cluster) %>% top_n(n = 250, wt = avg_log2FC) %>% filter(avg_log2FC > 0.25) 
length(sc.degs.filter$gene)

# load ST cluster markers
st.marker <- read.csv(file.path(tgDir,"E125_P0s1_Tangram_celltype89_dim20_res1.1_cluster_wilcox_filtered_markers.csv"))
head(st.marker)
st.marker.filter <- st.marker %>% filter(p_val_adj <= 0.05) %>% 
                    group_by(cluster) %>% arrange(desc(avg_log2FC)) %>%
                    top_n(300, wt = avg_log2FC) %>% filter(avg_log2FC > 0.25) 
st.genes <- unique(rownames(st_ob@assays$Spatial@counts))
sc.genes <- unique(rownames(sc_ob@assays$RNA@counts))
all.genes.scrna_and_spt <- unique(intersect(sc.genes, st.genes))
# Run MIA Results for all spots
MIA_all <- MIA(total_genes = length(all.genes.scrna_and_spt),
               single_cell.markers = sc.degs.filter,
               spatial.markers = st.marker.filter
              )
# Order MIA results & Filter NAs
E.all <- MIA_all %>% column_to_rownames("cluster")
E.all <- E.all[, order(colnames(E.all))]
is.na(E.all) <- do.call(cbind, lapply(E.all, is.infinite))

p <- pheatmap(E.all, cluster_cols = FALSE, cluster_rows = FALSE, fontsize=15, color = col.pal)
pdf(file.path(miaDir, "MIA_all_domains.pdf"), width = 8, height = 5)
print(p)
dev.off()

# MIA for D27
st.marker27.filter <- st.marker %>% filter(cluster %in% c(2,7)) %>% 
                        filter(p_val_adj <= 0.05) %>% 
                        group_by(cluster) %>% 
                        arrange(desc(avg_log2FC)) %>%
                        filter(avg_log2FC > 0.1)
st27.genes <- unique(rownames(st27@assays$Spatial@counts))
sc.genes <- unique(rownames(sc_ob@assays$RNA@counts))
D27.genes.scrna_and_spt <- unique(intersect(sc.genes, st27.genes))

# Run MIA Results
MIA27 <- MIA(total_genes = length(D27.genes.scrna_and_spt),
               single_cell.markers = sc.degs.filter,
               spatial.markers = st.marker27.filter
              )

# Order MIA results & Filter NAs
E.27 <- MIA27 %>% column_to_rownames("cluster")
E.27 <- E.27[, order(colnames(E.27))]
is.na(E.27) <- do.call(cbind, lapply(E.27, is.infinite))

# Plot MIA Cell Type Enrichment Heat Map
options(repr.plot.width = 6, repr.plot.height= 4)
pheatmap(E.27, cluster_cols = FALSE, cluster_rows = FALSE, fontsize=15, color = col.pal)

p <- pheatmap(E.27, cluster_cols = FALSE, cluster_rows = FALSE, fontsize=15, color = col.pal)
pdf(file.path(miaDir, "MIA_D27.pdf"), width = 6, height = 4)
print(p)
dev.off()