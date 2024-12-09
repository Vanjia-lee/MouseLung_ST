library(Seurat)
library(dior)

projectDir = '~/project/lung_10xST'
outDir = file.path(getwd(), "ST_h5")
if(!dir.exists(outDir)) dir.create(outDir)

slides = c("E125_S2", "E135_S3", "E145_S4", "E155_S4","E165_S2", "E175_S2", "E185_S3", "P0_S1")

for(slide in rev(slides)){
    stage = strsplit(slide,"_")[[1]][1]
    seurat_ob = readRDS(file.path(projectDir, "results", "0.rds", paste0("seurat_",slide,"_afterQC_rmHbMtRp.rds")))
    seurat_ob@meta.data$Stage = stage
    
    coords = GetTissueCoordinates(seurat_ob, image = slide)
    seurat_ob@meta.data$x = coords$imagecol
    seurat_ob@meta.data$y = coords$imagerow
    
    dior::write_h5(seurat_ob, 
                   file= file.path(outDir, paste0(slide,"_afterQC_rmHbMtRp.h5")),
                   object.type = 'seurat', assay.name = 'Spatial',
                   save.graphs = TRUE, # determing whether to save the graph(cell-cell similarity network)
                   save.scale = FALSE # determint whether to save the scale.data(dense matrix)
                  )
}

# save the common genes  between svgs and hvgs  
shvgs_dir = file.path(getwd(), "ST_shvgs")
if(!dir.create(shvgs_dir)) dir.create(shvgs_dir)
for(slide in slides){
    stage = strsplit(slide,"_")[[1]][1]
    seurat_ob = readRDS(file.path(projectDir, "results", stage, "0.rds", paste0("seurat_",slide,"_afterQC_rmHbMtRp.rds")))
    seurat_ob@meta.data$Stage = stage
    
    # shvgs
    load(file.path(projectDir,"results",stage,"2.svg_giotto",paste0(stage,"_4sections_spatial_genes.RData")))
    svg_slide = svg_lst[[slide]]
    print(length(svg_slide[[1]]$genes))
    print(length(svg_slide[[2]]$genes))
    if(length(svg_slide[[1]]$genes) < 2500 & length(svg_slide[[2]]$genes) < 2500){
        print("the number of svgs is less than 2500 ")
        svgs = intersect(svg_slide[[1]]$genes[1:length(svg_slide[[1]]$genes)],svg_slide[[2]]$genes[1:length(svg_slide[[2]]$genes)])
        print(length(svgs))
    }else{
        svgs = intersect(svg_slide[[1]]$genes[1:2500],svg_slide[[2]]$genes[1:2500])
        print(length(svgs))
    }
    
    shvgs = union(svgs,VariableFeatures(seurat_ob)[1:2000])
    print(length(shvgs))
    write.csv(shvgs, file.path(shvgs_dir,paste0(slide,"_shvgs.csv")))
}