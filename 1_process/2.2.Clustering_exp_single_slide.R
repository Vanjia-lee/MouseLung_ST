# loading packages
library(Seurat)
library(harmony)
library(scales)
library(clustree)
library(clusterProfiler)
library(org.Mm.eg.db)
library(preprocessCore)
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

source("~/project/lung_10xST/src/0.utilities/utilities.R")
projectDir = "~/project/lung_10xST"
slides = c("E125_S2","E135_S3","E145_S4","E155_S4","E165_S2","E175_S2","E185_S3","P0_S1")
stages = c("E125","E135","E145","E155","E165","E175","E185","P0")

options(repr.plot.width = 8, repr.plot.height = 8)
ob_lst = list()
for(id in seq(length(slides))){
    slide = slides[id]
    stage = stages[id]
    
    # Set Directory
    svgDir = file.path(projectDir,"results",stage,"2.svg_giotto")
    rds_inputdir = file.path(projectDir, "results", stage, "0.rds")
    
    # load st object
    print(paste0("Now is loading the spatial data of slide: ", slide))
    st_file = file.path(rds_inputdir,  paste0("seurat_",slide,"_afterQC_rmHbMtRp.rds"))
    st_ob = readRDS(st_file)
    st_ob@meta.data$Stage = stage    
    print(dim(st_ob))
    
    # rm Gm42418
    st_ob = st_ob[rownames(st_ob)!= "Gm42418",]
    
    # SCT normalization
    print('Now is normalization analysis with SCTransform')
    vars.to.regress = c('nCount_Spatial')
    st_ob = SCTransform(st_ob, assay = "Spatial", vars.to.regress = vars.to.regress, return.only.var.genes = FALSE, verbose = FALSE)
    
    # Load svgs
    load(file.path(projectDir,"results",stage,"2.svg_giotto",paste0(stage,"_4sections_spatial_genes.RData")))
    svg_slide = svg_lst[[slide]]
    print("Now is merge the svgs from Giotto_kmeans and Giotto_rank.")
    if(length(svg_slide[[1]]$genes) < 2500 & length(svg_slide[[2]]$genes) < 2500){
        print("the number of svgs is less than 2500 ")
        svgs = intersect(svg_slide[[1]]$genes[1:length(svg_slide[[1]]$genes)],svg_slide[[2]]$genes[1:length(svg_slide[[2]]$genes)])
        print(length(svgs))
    }else{
        svgs = intersect(svg_slide[[1]]$genes[1:2500],svg_slide[[2]]$genes[1:2500])
        print(length(svgs))
    }
    #svgs <- svgs[-c(grep("^mt-", svgs), grep("^Rp[sl]",svgs))]
    hvgs = VariableFeatures(st_ob)[1:2000]
    shvgs = union(svgs,hvgs)
    print(length(shvgs))
    
    st_ob = RunPCA(st_ob, features = shvgs, assay = "SCT", npcs = 50, verbose = FALSE)
    elbow.plot = ElbowPlot(st_ob, ndims = 50, reduction = "pca")
    print(elbow.plot)
    
    st_ob = FindNeighbors(st_ob, reduction = 'pca', k.param = 50, dims = 1:30, verbose = FALSE)
    resolutions <- c(seq(.5,2.0,.1))
    st_ob = FindClusters(st_ob, resolution = resolutions, algorithm=4, verbose = FALSE)
    
    # 借助clustree查看聚类状况    
    p0 = clustree(st_ob@meta.data, prefix = "SCT_snn_res.")
    print(p0)
    
    ob_lst = c(ob_lst, list(st_ob))
    
}

spot.size = c(7,4,3,3,2.5,2,2,2)
names(spot.size) = slides

resols <- c(1.2, 1.2, 1, 1.3, 0.8, 1.0, 1.0, 1.6)
names(resols) <- slides

for(id in seq(length(slides))){
    slide = slides[id]
    st_ob = ob_lst[[id]]
    
    # set outDir
    outDir = file.path(projectDir, "results", "single_stage_SCT_cluster",slide)
    if(!dir.exists(outDir)) dir.create(outDir, recursive = T)
    
    # set parameters
    #res = ifelse(slide == "E125_S2", 1.1, ifelse(slide == "E165_S2", 0.8, ifelse(slide == "P0_S1", 1.5, ifelse(slide == "E155_S4", 1.2, 1.0))))
    res = resols[slide]
    image = Images(st_ob)
    prefix = "SCT_snn_res."
    Idents(st_ob) = paste0(prefix, res)
    n = length(unique(Idents(st_ob)))
    cols2use = CustomCol(1:n)
    groupby = paste0(prefix, res)
    dims = 30
    slide_dim = paste0(slide, "_rmHbMtRp_SCT_hvg_svg_dim", dims,"_snn_res", res)
    
    # UMAP
    st_ob = RunUMAP(st_ob, reduction = 'pca', umap.method = 'umap-learn', metric = 'correlation', n.neighbors = 30, min.dist = 0.3, dims = 1:30, verbose = FALSE)
    saveRDS(st_ob, file.path(outDir, paste0("seurat_",slide,"_SCT_res",res,"_dim30_cluster.rds")))
    
    # cluster plot
    dp = DimPlot(st_ob, reduction = "umap", label = TRUE, cols = cols2use, label.size = 6) +
         labs (title = paste0('UMAP_res', res)) + theme(plot.title=element_text(hjust=0.5))
    gsp = SpatialDimPlot(st_ob, images = image, group.by = groupby,label = FALSE, label.size = 6, pt.size.factor = spot.size[id],repel = TRUE) + scale_fill_manual(values = cols2use) + NoLegend()
    p1 = wrap_plots(dp, gsp)
    ggsave(file.path(outDir, paste0(slide_dim,"_cluster_dimplot.pdf")), plot = p1, width = 7, height = 4)
    ggsave(file.path(outDir, paste0(slide_dim,"_cluster_dimplot.png")), plot = p1, width = 7, height = 4, dpi = 300, units = "in", device = "png")
    
    # nCount & nGenes for each cluster
    p2 = VlnPlot(st_ob, group.by = paste0(prefix, res), 
                 features = c("nCount_Spatial", "nFeature_Spatial"), 
                 assay = "Spatial", cols = cols2use,
                 pt.size = 0, fill.by = 'ident',stack = T, flip = T) + NoLegend() + 
         theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5)) + 
         geom_boxplot(width=.2,col="black",fill="white")
    print(p2)
    
    # marker gene analysis
    global_DEGs <- FindAllMarkers(st_ob, only.pos = T, test.use = "wilcox", logfc.threshold = 0.25, min.pct = 0.1)
    global_DEGs = global_DEGs %>% mutate( gene_diff = round(global_DEGs$pct.1 / global_DEGs$pct.2, 3)) %>% dplyr::select( gene, everything())
    filter_DEGs = global_DEGs %>% filter( p_val < 0.05 & p_val_adj < 0.05 ) %>% dplyr::select( gene, everything())
    top10_avg_logFC_markers <- filter_DEGs %>% group_by(cluster) %>% arrange(p_val,desc(avg_log2FC)) %>% top_n(10,avg_log2FC) %>% arrange(cluster)
    write.csv(filter_DEGs, file.path(outDir, paste0(slide_dim,"_cluster_wilcox_markers.csv")), row.names = F)
    write.csv(top10_avg_logFC_markers, file.path(outDir, paste0(slide_dim,"_cluster_wilcox_top10_avg_logFC_markers.csv")), row.names = F)
    
    # heatmap_top10_avg_logFC_gene
    topn_markers = top10_avg_logFC_markers %>% group_by(cluster) %>% top_n(10, gene_diff)
    topn_markers2vis = unique(as.vector(topn_markers$gene))
    colors2use = CustomCol(1:length(unique(st_ob@meta.data[,groupby])))
    ggheatmap = DoHeatmap( object = st_ob,
                           features = topn_markers2vis,
                           group.colors = colors2use,
                           group.by = levels(groupby), 
                           group.bar = T, label = T, 
                           size = 3, hjust = 0.5, angle = 0) + 
                theme(axis.text.y = element_text(size = 10 * 0.6, colour = "black"),
                      legend.text = element_text(size = 6),
                      legend.title = element_text(size = 7),
                      legend.key.size = unit(0.8, "lines")) + 
                guides(fill = guide_colorbar( title.position = "top")) + 
                       scale_fill_viridis(option = "D")
    ggsave(file.path(outDir, paste0(slide_dim,"_cluster_wilcox_top1_avg_logFC_gene_diff_markers_heatmap.pdf")), plot = ggheatmap, width = n * 0.7, height = n * 0.9)
    ggsave(file.path(outDir, paste0(slide_dim,"_cluster_wilcox_top1_avg_logFC_gene_diff_markers_heatmap.png")), plot = ggheatmap, width = n * 0.7, height = n * 0.9, dpi = 300, units = "in", device = "png")
    
   # BP enrichement
    marker_list <- split(topn_markers$gene, topn_markers$cluster)
    enrich_BP <- compareCluster(marker_list,
                                OrgDb = org.Mm.eg.db,
                                fun = "enrichGO",
                                keyType = "SYMBOL",
                                ont = "BP",
                                pAdjustMethod = "BH",
                                pvalueCutoff = 0.05,
                                qvalueCutoff = 0.05 )
    enrichplot <- dotplot(enrich_BP, showCategory = 10, label_format = 100) + 
                  theme(panel.border = element_rect(size = 2), 
                        plot.title = element_text(size = 18, hjust = 0.5),
                        axis.text.y = element_text(size = 13, margin = margin(0,5,0,0)),
                        axis.text.x = element_text(size = 16, margin = margin(5,0,0,0)),
                        axis.ticks = element_line(size = 1),
                        axis.ticks.length = unit(.2, "cm"),
                        legend.title = element_text(size = 16),
                        legend.text = element_text(size = 13) )
    ggsave(file.path(outDir, paste0(slide_dim,"_cluster_wilcox_top1_avg_logFC_gene_diff_markers_GO_BP_dotplot.pdf")), plot = enrichplot, width = n * 2, height = n * 2)
    ggsave(file.path(outDir, paste0(slide_dim,"_cluster_wilcox_top1_avg_logFC_gene_diff_markers_GO_BP_dotplot.png")), plot = enrichplot, width = n * 2, height = n * 2, dpi = 300, units = "in", device = "png") 
    
}