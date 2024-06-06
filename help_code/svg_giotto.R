#'@ stage the development of embryo, for example: E12.5, E13.5, E14.5, E15.5, E16.5, E18.5 and PO
#'@ slide slide id, eg: 1,2,3,4
#'@ point.size spot size when plot
#'@ k number of nearest neighbors based on physical distance when creating spatial network
#'@ ymd the date of QC

svg_gio = function(stage = stage, 
                   slide = slide, 
                   point.size = 1.2, 
                   k = 5, 
                   ymd = ymd,
                   projDir = projDir
                  ){
    section = paste0(stage,"_S",slide)
    rds_dir = file.path(projDir, "0.rds")
    slide_ob = readRDS(file.path(rds_dir, paste0("seurat_", section, "_afterQC_rmHbMtRp_",ymd,".rds")))
    
    spatial_count = as.matrix(slide_ob@assays$Spatial@counts)
    spatial_location = GetTissueCoordinates(slide_ob, image = section)
    spatial_location$x_round = round(spatial_location$imagecol)
    spatial_location$y_round = round(spatial_location$imagerow) * -1
    
    # set dir to save giotto
    svg_dir = file.path(projDir,"2.svg_giotto")
    if(!dir.exists(svg_dir)){
        dir.create(svg_dir, recursive = TRUE)
    }
    
    # 1. Create specific instructions for our Giotto analysis workflow
    instrs = createGiottoInstructions(save_plot = TRUE,
                                      show_plot = FALSE,
                                      return_plot = TRUE,
                                      save_dir = svg_dir
                                     )
    # 2. Create Giotto object
    gio_ob = createGiottoObject(raw_exprs = spatial_count,
                                spatial_locs = spatial_location[,c('x_round','y_round')],
                                instructions = instrs,
                                cell_metadata = spatial_location
                               )
    spatPlot(gobject = gio_ob,
             point_size = point.size,
             save_param = c(save_name = paste0("1.",section,"_spatplot"), 
                            base_width = 3.6, base_height = 3)
            )
    # 3. Filter genes and cells
    #@expression_threshold: threshold to consider a gene expressed
    gio_ob = filterGiotto(gobject = gio_ob, 
                          expression_threshold = 1,
                          gene_det_in_min_cells = 10, 
                          min_det_genes_per_cell = 500,
                          expression_values = c('raw'),
                          verbose = T)
    ## normalize
    gio_ob <- normalizeGiotto(gobject = gio_ob,
                              norm_methods = "standard",
                              scalefactor = 6000, 
                              verbose = T)
    ## add gene & cell statistics
    gio_ob <- addStatistics(gobject = gio_ob)
    spatPlot2D(gobject = gio_ob, 
               cell_color = 'nr_genes', 
               color_as_factor = F, 
               point_size = point.size,
               show_image = F,
               save_param = list(save_name = paste0("2.",section,"_nr_genes"),
                                 base_width = 4, base_height = 3))
    # 4. Create spatial network based on cell centroid physical distances
    ## 4.1.Spatial network_kNN
    ## k: number of nearest neighbors based on physical distance, we recommend approximately 5-10 neighbors per cell.
    ## maximum_distance_knn: distance cuttof for nearest neighbors to consider for kNN network
    ## minimum_k: minimum nearest neigbhours if maximum_distance != NULL
    print("Creating spatial network with KNN now.")
    gio_ob = createSpatialNetwork(gobject = gio_ob,
                                  method = "kNN",
                                  k = k,
                                  maximum_distance_knn = 400,
                                  minimum_k = 1,
                                  name = "kNN_network"
                                ) 
    spatPlot(gio_ob, show_network = T, 
             point_shape = "no_border",
             network_color = "black",
             spatial_network_name = "kNN_network",
             point_size = 0.5,
             save_param = list(save_name = paste0("3.",section,"_knn_network_k",k),
                               base_width = 3.6, base_height = 3
                            )
            )
    
    ## 4.2.Spatial network_Delaunay
    ## maximum_distance_delaunay: distance cuttof for nearest neighbors to consider for Delaunay network
    ## minimum_k: minimum nearest neigbhours if maximum_distance != NULL
    print("Creating spatial network with Delauany now.")
    gio_ob = createSpatialNetwork(gobject = gio_ob,
                                  method = "Delaunay",
                                  minimum_k = 0,
                                  maximum_distance_delaunay = "auto",
                                  name = "Delaunay_network"
                                )
    spatPlot(gio_ob, show_network = T,
             point_shape = "no_border",
             network_color = "black",
             spatial_network_name = "Delaunay_network",
             point_size = 0.5,
             save_param = list(save_name = paste0("3.",section,"_delaunay_network"),
                               base_width = 3.6, base_height = 3
                            )
            )
    
    
    # 5. Spatial genes
    ## Calculate spatial genes based on kmeans binarization or rank binarization
    ## bin_method(): method to binarize gene expression
    ## binarize: Each gene is binarized (0 or 1) in each cell with kmeans (k = 2) or based on rank percentile
    print("Calculate spatial genes now.")
    svg_lst = c()
    for(net in c("kNN_network","Delaunay_network")){
        bin_lst = list()
        for(method in c("kmeans", "rank")){
            svgs = binSpect(gio_ob, 
                            bin_method = method,
                            calc_hub = T, 
                            hub_min_int = 5,
                            spatial_network_name = net
                            )
            spatGenePlot(gio_ob, expression_values = "scaled",
                         genes = svgs$genes[1:6], 
                         cow_n_col = 2, 
                         point_size = point.size,
                         cell_color_gradient = c("white", "white", "red"),
                         point_shape = 'border', point_border_stroke = 0.1,
                         show_network = TRUE, network_color = 'lightgrey',
                         save_param = list(save_name = paste0("4.",section,"_",net,"_spatial_genes_",method),
                                           base_width = 7, base_height = 8)
                         )
            svgs = svgs[svgs$adj.p.value < 0.05,]
            bin_lst = c(bin_lst,list(svgs))
        }
        svg_lst = c(svg_lst,bin_lst)
    } 
    
    
    
    saveRDS(gio_ob, file.path(svg_dir,paste0("giotto_",section,"_svg.rds")))
    
    return(svg_lst)
}