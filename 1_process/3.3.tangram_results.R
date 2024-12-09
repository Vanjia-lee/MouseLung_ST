library(Seurat)
library(data.table)
library(future)

rm(list = ls());gc()

# Set work directory
projectDir = '~/project/lung_10xST'
scDir = file.path(projectDir, "data/reference_sc/GSE165063_Development_2021/processing")
tgDir = file.path(projectDir, "results", "3.2.Cluster_deconvolution_Tangram")
sc_meta = read.csv(file.path(scDir,"SC_ref_E12_P0_timepoint_celltype_new_subcelltype_metadata.csv"), row.names = 1)
head(sc_meta)

slides <- c("E125_S2","E135_S3","E145_S4","E155_S4","E165_S2","E175_S2","E185_S3","P0_S1")
cts_all2 <- c(paste0("Epithelial_",1:18),paste0("Endothelial_",1:28),paste0("Mesenchymal_",1:43))

# Tidy the results of Tangram mapping with the mode of cells
tg_lst2 = list()
for(slide in slides){
    tg.df = read.csv(file.path(tgDir, "20230823_tangram_origQC_training_genes_mode_cells",paste0(slide, "_tangram_mode_cells_training_genes_ct89.csv")),check.names = F, row.names = 1)
    tt = apply(tg.df,1, function(row) ifelse(row < max(row), 0, 1)) 
    tg = data.frame(index = rownames(tt))
    for(ct in cts_all2){
        cells = rownames(sc_meta)[sc_meta$new_subcelltype == ct]
        tg[,ct] = apply(tt[,cells], 1, sum)
    }
    rownames(tg) = paste0(tg[,1],"_",which(slides == slide))
    tg = tg[,-1]
    #tg = apply(tg,1,function(row) row/sum(row))  
    print(table(apply(tg,1,function(row) length(which(row>0)))))   
    print(table(apply(tg,1,function(row) length(which(row==1)))))                  
               
    tg_lst2 = c(tg_lst2, list(tg))           
}
                      
for(id in seq(length(tg_lst2))){
    tg = tg_lst2[[id]]
    print(rownames(tg)[apply(tg,1,sum)==0])                     
}
                      
tg.all89 = do.call(rbind, tg_lst2)
tg.all89 = tg.all89[rownames(tg.all89) != "TCGTATAGTGCAATTA-1_8",]

# Count the proportion of each cell type in each spot
tg.all89 = apply(tg.all89, 1, function(x) x/sum(x))
outDir = file.path(tgDir, "tangram_origQC_training_genes_mode_cells","Clustering_celltype89_Deconvolution")
if(!dir.exists(outDir)) dir.create(outDir, recursive = T)
write.csv(tg.all89, file.path(outDir,"E125_P0s1_origQC_tangram_89CT_modeCell_trainGene_deconvolution_mtx.csv"), row.names = TRUE)                 
                 