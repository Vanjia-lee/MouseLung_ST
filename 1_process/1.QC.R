library(Seurat)
library(scales)
library(Matrix)
library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(ggthemes)
library(patchwork)
library(cowplot)

projectDir = './'

#  Create Seurat Object for each slide (eg:E155)
outDir <- file.path(projectDir, "results", "SpaceRanger_RDS_rmImgBG")
ifelse(!dir.exists(outDir)) dir.create(outDir, recursive = T)
ymd <- format(Sys.time(), format="%Y%m%d" )
stage <- 'E155'

st_lst <- list()
options(repr.plot.width = 18, repr.plot.height = 10)
for(id in c("2A","2B","2C","2D")){
    slide = paste0(stage,"_", id)
    lung_dir = file.path(projectDir,'data/SpaceRanger',slide)
    print(lung_dir)
    
    slice = paste0(stage,"_S", ifelse(id == "2A", 1, ifelse(id == "2B", 2, ifelse(id == "2C", 3, 4))))
    slice_ob = Load10X_Spatial(lung_dir, slice = slice)
    slice_ob[["orig.ident"]] <- slice 
    slice_ob@project.name <- slice
    levels(slice_ob@active.ident) <- slice
    print(dim(slice_ob))
    
    # compute the percentage of mitochondria genes
    slice_ob <- PercentageFeatureSet(slice_ob, pattern = '^mt-', col.name = 'percent_mito')
    # compute the percentage of hemoglobin genes
    slice_ob <- PercentageFeatureSet(slice_ob, pattern = "^Hb[^(eps)]", col.name = "percent_hb")
    # compute the percentage of ribosome genes
    slice_ob <- PercentageFeatureSet(slice_ob, pattern = '^Rp[sl]', col.name = 'percent_ribo')

    saveRDS(slide_ob, file.path(outDir, paste0("seurat_", section, "_afterQC_rmHbMtRp_", ymd, ".rds")))
    
    # visualization before QC
    vln <- VlnPlot(slice_ob, features = c("nCount_Spatial","nFeature_Spatial","percent_mito","percent_hb","percent_ribo"),
               pt.size = 0.1, ncol = 5) + NoLegend() 
    st_feature <- SpatialFeaturePlot(slice_ob, features = c("nCount_Spatial","nFeature_Spatial","percent_mito","percent_hb","percent_ribo"),
                                pt.size.factor = pt.size, alpha = 1, stroke = 0, ncol = 5) 
    plot <- wrap_plots(vln, st_feature, ncol=1)   
    print(plot)
    
    st_lst <- c(st_lst, list(slice_ob))
}

st_ob <- merge(st_lst[[1]],st_lst[2:length(st_lst)])

# visualization before QC
options(repr.plot.width = 6, repr.plot.height = 4)
beforeQC_vln <- VlnPlot(st_ob, group.by = "orig.ident", 
                       features = c("nCount_Spatial", "nFeature_Spatial","percent_mito"), 
                       assay = "Spatial", cols = brewer.pal(11,"RdYlGn")[8:11],
                       pt.size = 0.1,  fill.by = 'ident',
                       stack = T, flip = T) + NoLegend() +
                 theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
beforeQC_vln <- beforeQC_vln + geom_boxplot(width=.2,col="black",fill="white")
beforeQC_vln

# QC 
genenum <- 1000
for(section in slides){
    slide_ob <- readRDS(file.path(rmImgBG_dir,paste0("seurat_",section,"_spaceranger_rmImgBG_",ymd,".rds")))
    # QC by gene filter
    slide_ob <- subset(slide_ob, subset = nFeature_Spatial >= genenum & percent_mito <= 15)
    print(dim(slide_ob))
    # Filter Hemoglobin gene
    slide_ob <- slide_ob[!grepl("^Hb[^(eps)]", rownames(slide_ob)),]
    # Filter mito genes
    slide_ob <- slide_ob[!grepl("^mt-", rownames(slide_ob)),]
    # Filter ribosome genes
    slide_ob <- slide_ob[!grepl("^Rp[sl]", rownames(slide_ob)),]

    saveRDS(slide_ob, file.path(outrDir, paste0("seurat_", section, "_afterQC_rmHbMtRp_", ymd, ".rds")))
}

sub_ob = subset(st_ob, subset = nFeature_Spatial >= genenum & percent_mito <= 15)
sub_ob = sub_ob[!grepl("^Hb[^(eps)]", rownames(sub_ob)),]
sub_ob = sub_ob[!grepl("^Hb[^(eps)]", rownames(sub_ob)),]
sub_ob = sub_ob[!grepl("^mt-", rownames(sub_ob)),]
sub_ob = sub_ob[!grepl("^Rp[sl]", rownames(sub_ob)),]
dim(sub_ob)
saveRDS(sub_ob, file.path(outDdir, paste0("seurat_",stage,"_afterQC_rmHbMtRp_",ymd,".rds")))

# visualization afterQC
afterQC_vln <- VlnPlot(sub_ob, group.by = "orig.ident", 
                       features = c("nCount_Spatial", "nFeature_Spatial","percent_mito"), 
                       assay = "Spatial", cols = brewer.pal(11,"RdYlGn")[8:11],
                       pt.size = 0.1,  fill.by = 'ident',
                       stack = T, flip = T) + NoLegend() +
                 theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
afterQC_vln <- afterQC_vln + geom_boxplot(width=.2,col="black",fill="white")
afterQC_vln