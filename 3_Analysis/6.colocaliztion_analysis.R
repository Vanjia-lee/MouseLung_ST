library(Seurat)
library(proxy)
library(data.table)
library(tidyverse)
library(foreach)
library(doParallel)
library(patchwork)

rm(list=ls())

# set work directory
projectDir <- '~/project/lung_10xST'
stDir <- file.path(projectDir, "results", "3.2.Cluster_deconvolution_Tangram","tangram_origQC_training_genes_mode_cells/Clustering_celltype89_Deconvolution")
mapDir <- file.path(projectDir, "results/cytoSPACE_origQC")
outDir <- file.path(projectDir, "results/6.CCI_cellchat2_stComm/colocalization")
if(!dir.exists(outDir)) dir.create(outDir, recursive = T)

# load ST object
st_ob <- readRDS(file.path(stDir,"seurat_E125_P0s1_origQC_tangram_89CT_modeCell_trainGene_deconvolution_res1.1.rds"))
Idents(st_ob) = 'TangramNew89_snn_res.1.1'
slides <- unique(st_ob$orig.ident)
names(slides) <- unique(st_ob$Stage)

# load functions
source(file.path(projectDir, "src/colocalization/coloc_permutate.R"))
source(file.path(projectDir, "src/colocalization/visualization.R"))

# colocalization analysis for cell types (SC) mapping to ST domain 2 and 7 st E185 (eg:E185)

stage <- "E185"
domains <- c("2","7")
slide <- slides[stage]

for( dm in domains){
    st_dm <- subset(st_ob, subset = Stage == stage, idents = dm)
    coords <- GetTissueCoordinates(st_dm, image = slide, scale = NULL) 
    # load cytospace results
    map.all <- read.csv(file.path(mapDir, paste0(slide, "_cytospace"),"assigned_locations.csv"), row.names = 1)
    map.all$SpotID <- paste0(map.all$SpotID, "_",which(slide == slides))
    map_dm <-  map.all[map.all$SpotID %in% Cells(st_dm),]
    # add original coordinate infomation for the full resolution image
    map_dm$imagerow <- 0
    map_dm$imagecol <- 0
    for(i in rownames(coords)){
        map_dm$imagerow[which(i == map_dm$SpotID)] = coords[i,"imagerow"]
        map_dm$imagecol[which(i == map_dm$SpotID)] = coords[i,"imagecol"] 
    }
    map_dm$Domain <- st_dm$`TangramNew89_snn_res.1.1`[match(map_dm$SpotID,Cells(st_dm))]
    # calculate the dist
    d.all <- getD(map_dm, distant = 0)
    d.all$reorder.pair <- apply(d.all, 1, FUN = sort_repair)
    d.all$spotID <- map_dm$SpotID[match(d.all$UCID, rownames(map_dm))]
    d.unique <- d.all %>% group_by(spotID) %>% dplyr::distinct(reorder.pair, .keep_all = TRUE)

    celltypes <- sort(unique(map_dm$CellType))
    dist.lst <- list()
    for(celltype in celltypes){
        dist.cell = d.unique %>% filter(CellType == celltype)
        if(nrow(dist.cell)>0){
            df.cell = as.data.frame(table(dist.cell$adjCellType))
            colnames(df.cell) = c("adjCellType","Freq")
            df.cell$perc = df.cell$Freq/sum(df.cell$Freq) * 100
            dist.lst <- c(dist.lst, list(df.cell))
            names(dist.lst)[length(dist.lst)] <- celltype
        }
    }
    
    # permutation test 1000 times
    cts.coloc <- permutateD(map_dm, dist_lst = dist.lst, threads = 10, 
                            distant = 0, outdir = outDir, slide = slide, 
                            domain = dm)
    print(dim(cts.coloc))
    
    # visualization
    p1 <- circle_coloc(cts.coloc, slide = slide, domain = dm, outdir = outDir)
    p2 <- Dotplot_coloc(cts.coloc, slide = slide, domain = dm, outdir = outDir)
    print(p1 | p2)
}
