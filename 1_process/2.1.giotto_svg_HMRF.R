library(Seurat)
library(Giotto)
library(smfishHmrf)
library(tidyverse)
library(Matrix)
library(data.table)

# E155
stage <- "E155"
projectDir <- '~project/lung_10xST'
stage_dir <- file.path(projectDir, "results", stage)
svg_dir = file.path(projectDir,"results","2.SVGs_giotto")
if(!dir.exists(svg_dir)) dir.create(svg_dir, recursive = TRUE)

# Find SVGs
source(file.path(projectDir,"src","step2.svg_giotto.R"))
ymd_qc = 'XXXXXX'
slides = paste0(stage,"_S", 1:4)
svg_lst = list()
for(id in seq(slides)){
    slide_svg = svg_gio(stage = stage, slide = id, ymd = ymd_qc, projDir = stage_dir)
    svg_lst = c(svg_lst, list(slide_svg))
}

save(svg_lst, file = file.path(svg_dir,paste0(stage,"_4sections_spatial_genes.RData")))

svg_s1 = intersect(svg_lst[[slides[1]]][[1]]$gene[1:630], svg_lst[[slides[1]]][[2]]$gene[1:630])
length(svg_s1)
svg_s2 = intersect(svg_lst[[slides[2]]][[1]]$gene[1:580], svg_lst[[slides[2]]][[2]]$gene[1:580])
length(svg_s2)
svg_s3 = intersect(svg_lst[[slides[3]]][[1]]$gene[1:750], svg_lst[[slides[3]]][[2]]$gene[1:750])
length(svg_s3)
svg_s4 = intersect(svg_lst[[slides[4]]][[1]]$gene[1:600], svg_lst[[slides[4]]][[2]]$gene[1:600])
length(svg_s4)

selected_svg_lst = list(svg_s1, svg_s2, svg_s3, svg_s4)
names(selected_svg_lst) = slides
save(selected_svg_lst, file = file.path(svg_dir,paste0(stage,"_4sections_spatial_genes_selected.RData")))

# HMRF
source(file.path(projectDir,"src","step3.HMRF_giotto.R"))
svg_dir = file.path(stage_dir,"2.svg_giotto")
hmrf_lst = list()
for(slide in slides){
    hmrf_dir = file.path(stage_dir, "3.Clustering", paste0("3.2.cluster_HMRF_",slide))
    if(!dir.exists(hmrf_dir)) dir.create(hmrf_dir, recursive = TRUE)
    
    slide_gio = readRDS(file.path(svg_dir, paste0("giotto_",slide,"_svg.rds") ))
    svgs = selected_svg_lst[[slide]]
    K = ifelse(slide == "E155_S4",7,6)
    hmrf_svg = hmrf(gobject = slide_gio, slide = slide, k = K, spatial.genes = svgs, nSGs = length(svgs), spatial.network = "kNN_network", point.size = 2, hmrf_dir = hmrf_dir)
    hmrf_lst = c(hmrf_lst, list(hmrf_svg))
}