#!/bin/Rscript
#===============================================================================
# This script is used to get the cell-cell communication with spaTalk.
# Created by Wenjia Li
# Created time: 2024.03.12

#========================= import packages =====================================
rm(list=ls())
suppressWarnings({
    suppressPackageStartupMessages(library("optparse"))
    suppressPackageStartupMessages(library("Seurat"))
    suppressPackageStartupMessages(library("CellChat"))
    suppressPackageStartupMessages(library("tidyverse"))
})

#========================= command line parameters setting ======================

option_list = list(
    make_option(c("--slide","-s"), type = "character", help = "The id of ST slide."),
    make_option(c("--slide2","-S"), type = "character", help = "The origin names of ST slide"),
    make_option(c("--threads","-c"), type = "integer", default = 4, help = "Number of CPU cores to use")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#================================================================================
# parse the command line parameters
#================================================================================
if (is.null(opt$slide)){
    print("Slide is unavaliable!")
}else{
    slide = opt$slide
}

if (is.null(opt$slide2)){
    print("Slide2 is unavaliable!")
}else{
    slide2 = opt$slide2
}

if (is.null(opt$threads)){
    threads = 4
}else{
    threads = opt$threads
}

################## set work directory ###########################################
projectDir <- '~/project/lung_10xST'
stDir <- file.path(projectDir, "results", "3.2.Cluster_deconvolution_Tangram","tangram_origQC_training_genes_mode_cells/Clustering_celltype89_Deconvolution")
cccDir <- file.path(projectDir, "results/6.CCI_cellchat2_stComm")
if(!dir.exists(cccDir)) dir.create(cccDir, recursive = T)
outDir <- file.path(cccDir, "CellChat2_ST")
if(!dir.exists(outDir)) dir.create(outDir, recursive = T)

st_ob <- readRDS(file.path(stDir,"seurat_E125_P0s1_origQC_tangram_89CT_modeCell_trainGene_deconvolution_res1.1.rds"))
sub.st <- subset(st_ob, subset = orig.ident == slide)

# Prepare input data for CelChat analysis
data.input <- GetAssayData(sub.st, slot = "data", assay = "Spatial")
meta <- data.frame(labels = paste0("D",sub.st$TangramNew89_snn_res.1.1), samples = slide, stages = sub.st$Stage, row.names = Cells(sub.st))
meta$labels <- factor(meta$labels, levels = paste0("D",levels(sub.st$TangramNew89_snn_res.1.1)))
spatial.locs <- GetTissueCoordinates(sub.st, image = slide, scale = NULL, cols = c("imagerow", "imagecol")) 
# Spatial factors of spatial coordinates
scalefactors <- jsonlite::fromJSON(txt = file.path(projectDir, "data/SpaceRanger", slide2, "spatial/scalefactors_json.json"))
spot.size <- 65 # the theoretical spot size (um) in 10X Visium
conversion.factor <- spot.size/scalefactors$spot_diameter_fullres
spatial.factors <- data.frame(ratio = conversion.factor, tol = spot.size/2)
d.spatial <- computeCellDistance(coordinates = spatial.locs, ratio = spatial.factors$ratio, tol = spatial.factors$tol)
print(min(d.spatial[d.spatial!=0])) # this value should approximately equal 100um for 10X Visium data

# 1.Create a CellChat object
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels",
                           datatype = "spatial", coordinates = spatial.locs, spatial.factors = spatial.factors)
# 2.Set the ligand-receptor interaction database
CellChatDB <- CellChatDB.mouse
# use a subset of CellChatDB for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB)
# set the used database in the object
cellchat@DB <- CellChatDB.use

# 3.Preprocessing the expression data for cell-cell communication analysis
cellchat <- subsetData(cellchat)
future::plan("multisession", workers = threads) 
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat, variable.both = F)
cellchat <- computeCommunProb(cellchat, 
                              type = "truncatedMean", 
                              trim = 0.1,
                              distance.use = TRUE, 
                              interaction.range = 250, 
                              scale.distance = 0.01,
                              contact.dependent = TRUE, 
                              contact.range = 100)
cellchat <- filterCommunication(cellchat, min.cells = 10)
# Infer the cell-cell communication at a signaling pathway level
cellchat <- computeCommunProbPathway(cellchat)
print("Now have inferred ccc at singaling pathway level.")
# Calculate the aggregated cell-cell communication network
print("Now is aggregating.")
cellchat <- aggregateNet(cellchat)

# save cellchat obj
print("Now is saving the cellchat object.")
saveRDS(cellchat, file = file.path(outDir, paste0(slide, "_CellChat.rds")))
