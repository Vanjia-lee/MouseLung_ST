library(Seurat)
library(tidyverse)
library(future)
library(pheatmap)
library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(cowplot)

rm(list=ls())
options(stringsAsFactors = F)

source("~/project/lung_10xST/src/0.utilities/utilities.R")
source("~/project/lung_10xST/src/0.utilities/multiSpatialPlot.R"))

# set work directory
projectDir = '~/project/lung_10xST'
outDir <- file.path(projectDir, "results/5.comparison_D4distal_D6proximal/smoothGBA_similar_PD_pattern_genes")
tgDir = file.path(projectDir, "results", "3.2.Cluster_deconvolution_Tangram","tangram_origQC_training_genes_mode_cells/Clustering_celltype89_Deconvolution")
st_ob <- readRDS(file.path(tgDir,"seurat_E125_P0s1_origQC_tangram_89CT_modeCell_trainGene_deconvolution_res1.1.rds"))
st_ob

Idents(st_ob) = 'TangramNew89_snn_res.1.1'
n = length(unique(Idents(st_ob)))
cols2use = CustomCol(1:n)

# define function
smooth_exp <- function(object, slide = NULL, sigma = 0.5){
    st <- subset(object, orig.ident == slide)
    # calculate W matrix accroding to the coordinate of spots
    coords <- GetTissueCoordinates(st, image = slide, scale = NULL, cols = c("imagerow","imagecol"))
    d.all <- proxy::dist(coords)
    d.min <- min(d.all)
    d.all <- as.matrix(d.all/d.min)
    d.w <- exp(-(d.all ** 2)/ (2 * sigma ** 2))
    
    # smooth exprssion
    st.exp <- GetAssayData(st, assay = "Spatial", slot = "data")
    exp.smooth <- future.apply::future_apply(st.exp, 1, function(x) x %*% d.w)
    exp.smooth <- t(exp.smooth)
    colnames(exp.smooth) <- colnames(st.exp)  
    return(exp.smooth)                                         
}
                                             
fisherT_exp <- function(smooth.exp, meta.data, group.by = NULL, 
                        gscore = c("norm","mean","q3"), 
                        fscore = c("norm","mean","q3"),
                        alternative = c("greater","less","two.sided")
                       ){
    #gene.exp <- smooth.exp[gene,]
    # z-score of gene expression
    #gene.score <- (gene.exp - mean(gene.exp))/sd(gene.exp)
    gene.score <- (smooth.exp - mean(smooth.exp))/sd(smooth.exp)
    
    # binary gene sorce
    Gscore <- ifelse(gscore == "norm", 0, ifelse(gscore == "mean", mean(gene.score), quantile(gene.score, 0.75)))
    Fscore <- ifelse(fscore == "norm", 0, ifelse(fscore == "mean", mean(meta.data[[group.by]]), quantile(meta.data[[group.by]],0.75)))
  
    gene.binary <- ifelse(gene.score > Gscore, 1, 0)
    Fsore.binary <- ifelse(meta.data[[group.by]] > Fscore, 1, 0)
    p.value <- fisher.test(as.matrix(table(Fsore.binary, gene.binary))[2:1,2:1], alternative = alternative)$p.value
    return(p.value)
}
                                             
# smooth expression
load(file.path(outDir, "ST_expression_guassian_smoothed_sigma0.5.RData"))
# add smooth_sigma0.5 assay
st_ob[["Smooth0.5"]] <- CreateAssayObject(data = exp.allsmooth2)

# genes with proximal and distal pattern
p.genes <- c('Sox2','Foxj1','Krt15','Muc5b','Ascl1')
d.genes <- c('Etv5','Sftpc','Wnt2') 
                                             
# Fisher test for proximal genes                                           
options(future.globals.maxSize = 15000 * 1024^2)
plan(multisession, workers = 10)    
DefaultAssay(st_ob) <- "Smooth0.5"
st_ob@meta.data <- st_ob@meta.data[,c("orig.ident","Stage","TangramNew89_snn_res.1.1","proximal","distal")]
# binaration with z-score > q3
st.meta <- st_ob@meta.data
pval5p.q3 <- future.apply::future_apply(exp.allsmooth2[apply(exp.allsmooth2,1,sum)>0,], 1, function(x) 
    fisherT_exp(x, meta.data = st.meta, group.by = "proximal", gscore = "mean", fscore = "q3", alternative = "two.sided"))
# correct P-value
adjp5p.q3 <- p.adjust(pval5p.q3,"fdr",n=length(pval5p.q3))  
proximal5.q3 <- data.frame(pval = pval5p.q3, adj.pval = adjp5p.q3)
fproximal5.q3 <- proximal5.q3 %>% filter(adj.pval < 0.001) %>% arrange(adj.pval) 
fproximal5.q3 <- mutate(fproximal5.q3, rank = seq(nrow(fproximal5.q3)))                                        
dim(fproximal5.q3)
head(fproximal5.q3)
write.csv(fproximal5.q3, file.path(outDir, "proximal5_pattern_genes_smooth0.5_filter_q3_adjP0.001.csv"))

# Fisher test for distal genes     
DefaultAssay(st_ob) <- "Smooth0.5"
# binaration with z-score > q3
st.meta <- st_ob@meta.data
pval5d.q3 <- future.apply::future_apply(exp.allsmooth2[apply(exp.allsmooth2,1,sum)>0,], 1, function(x) 
    fisherT_exp(x, meta.data = st.meta, group.by = "distal", gscore = "mean", fscore = "q3", alternative = "two.sided"))
# correct P-value
adjp5d.q3 <- p.adjust(pval5d.q3,"fdr",n=length(pval5d.q3))                                       
distal5.q3 <- data.frame(pval = pval5d.q3, adj.pval = adjp5d.q3)
fdistal5.q3 <- distal5.q3 %>% filter(adj.pval < 0.001) %>% arrange(adj.pval)
fdistal5.q3 <- mutate(fdistal5.q3, rank = seq(nrow(fdistal5.q3)))
dim(fdistal5.q3)
head(fdistal5.q3) 
write.csv(fdistal5.q3, file.path(outDir, "distal4_pattern_genes_smooth0.5_filter_q3_mean_adjP0.001.csv"))                                        