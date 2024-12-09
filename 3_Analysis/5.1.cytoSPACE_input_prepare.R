
library(Seurat)
library(data.table)
library(tidyr)
library(Matrix)

rm(list = ls());gc()

source("~/biosoft/cytospace/cytospace/Prepare_input_files/generate_cytospace_from_seurat_object.R")

# set dir
projectDir <- '~/project/lung_10xST'
refDir = file.path(projectDir, "data/reference_sc/GSE165063_Development_2021/processing/")
ref_file = file.path(refDir, "SC_ref_E12_P0.rds")
outdir <- file.path(projectDir, "data/cytoSPACE_input")
if(!dir.exists(outdir)) dir.create(outdir, recursive = T)

# SC data prepare
sc_ob = readRDS(ref_file)
Idents(sc_ob) = "celltype"
generate_cytospace_from_scRNA_seurat_object(sc_ob, dir_out=outdir, fout_prefix='', write_sparse=FALSE, rna_assay='RAW_COUNTS')

# ST data prepare
slides = c("E125_S2","E135_S3","E145_S4","E155_S4","E165_S2","E175_S2","E185_S3","P0_S1")
names(ymd) <- slides
for(slide in slides){
    print(paste0("Now extracting the STdata for ", slide))
    outdir_st <- file.path(outdir, paste0(slide,"_ST"))
    if(!dir.exists(outdir_st)) dir.create(outdir_st, recursive = T)
    
    st_file = file.path(projectDir, "data/SpaceRanger_RDS_rmImgBG",paste0("seurat_", slide, "_spaceranger_rmImgBG.rds"))
    st_ob = readRDS(st_file)
    generate_cytospace_from_ST_seurat_object(st_ob, dir_out=outdir_st, fout_prefix='', write_sparse=FALSE, slice=slide)
}
