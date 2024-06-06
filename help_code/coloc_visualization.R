circle_coloc <- function(coloc.all, slide = NULL, domain = NULL, outdir = NULL){
    require(igraph)
    require(ggraph)
    require(tidygraph)
    
    # set color code for all cell types
    epi <- c("Early Epithelium", "Transitional", "AT1", "AT2", "Mki67 AT2", "Ciliated", "Secretory", "Neuroendocrine")
    endo <- c("Arterial maEC","Venous maEC","Prolif gCap","gCap", "aCap", "Lymphatic")
    fibro <- c("Prolif Wnt2 FB", "Wnt2 FB", "Prolif Myo FB", "Myofibroblast", "Adventitial FB", "Pericyte", "Mesothelium", "Smooth Muscle", "Cardiomyocyte")   
    celltypes.all <- c(epi, endo, fibro)
    color.code <- c(rep("#A0D3FF",8), rep("#FF8AD0",6), rep("#A886F6",9))
    names(color.code) <- celltypes.all
    
    # filter the repetitive pairs of celltype
    sort_repair <- function(row1) {
        cells <- sort(row1[1:2])
        repair <- paste(cells[1], cells[2], sep = "_")
        return(repair)
    }
    coloc.all$reorder_pair <- apply(coloc.all, 1, FUN = sort_repair)
    co.locs <- coloc.all %>% arrange(reorder_pair,desc(perc)) %>% dplyr::distinct(reorder_pair, .keep_all = TRUE)
    
    # create node df
    nodes <- data.frame(id = seq(length(celltypes.all)), label = factor(celltypes.all, levels = celltypes.all))
    edges <- co.locs %>% left_join(nodes, by = c("CellType" = "label")) %>% 
             rename(from = id) %>% left_join(nodes, by = c("adjCellType" = "label")) %>% 
             rename(to = id) %>% select(from, to, perc, pvalue)
    g <- tbl_graph(nodes = nodes, edges = edges)
    # set features of point and edge
    V(g)$color <- color.code[match(V(g)$label, names(color.code))] 
    E(g)$weight <- -log(edges$pvalue + 0.000001)
    
    angle <- 360 * (c(1:length(celltypes.all)) - 0.5)/length(celltypes.all)
    hjust <- ifelse(angle > 180, 1, 0)
    angle <- ifelse(angle > 180, 90 - angle + 180, 90 - angle)
    p <- ggraph(g, layout = 'linear', circular = TRUE) +
         geom_edge_arc2(aes(edge_width = weight), color="grey") +
         geom_node_point(aes(color = label),size = 10, show.legend = FALSE) +
         geom_node_text(aes(x = x * 1.1, y = y * 1.1, label = label), size = 5,  angle = angle, hjust = hjust) +
         expand_limits(x = c(-1.6, 1.6), y = c(-1.6, 1.6)) +
         scale_color_manual(values = color.code) +
         scale_edge_width(range = c(0.5,2))+
         theme_graph() +
         theme(legend.position = "right", legend.key.size = unit(0.4, "cm"),
               text = element_text(family = "sans"),
               plot.title = element_text(hjust = 0.5, family = "sans", color = "#404040"),
               legend.title = element_text(family = "sans", face = "bold", size = 13, color = "#404040"),
               legend.text = element_text(family = "sans", size = 12, color = "#404040")
              ) +
         labs(title = paste0(stage,"_D",domain), color = "Cell Types", edge_width = "Significance") 
    ggsave(file.path(outdir, paste0(slide,"_D",domain,"_Spots_colocalization_celltype_pairs_significant_CirclePlot.pdf")), p, width = 10, height = 10)    
    return(p)
}

Dotplot_coloc <- function(coloc.all, slide = NULL, domain = NULL, outdir = NULL){
    require(ggplot2)
    
    # set color code for all cell types
    epi <- c("Early Epithelium", "Transitional", "AT1", "AT2", "Mki67 AT2", "Ciliated", "Secretory", "Neuroendocrine")
    endo <- c("Arterial maEC","Venous maEC","Prolif gCap","gCap", "aCap", "Lymphatic")
    fibro <- c("Prolif Wnt2 FB", "Wnt2 FB", "Prolif Myo FB", "Myofibroblast", "Adventitial FB", "Pericyte", "Mesothelium", "Smooth Muscle", "Cardiomyocyte")   
    celltypes.all <- c(epi, endo, fibro)
    color.code <- c(rep("#A0D3FF",8), rep("#FF8AD0",6), rep("#A886F6",9))
    names(color.code) <- celltypes.all
    
    coloc.all$CellType <- factor(coloc.all$CellType, levels = intersect(celltypes.all, coloc.all$CellType))
    coloc.all$adjCellType <- factor(coloc.all$adjCellType, levels = intersect(celltypes.all, coloc.all$adjCellType))
    
    xmin.epi <- head(intersect(epi,cts.coloc$CellType),1) 
    xmax.epi <- tail(intersect(epi,cts.coloc$CellType),1) 
    ymin.epi <- head(intersect(epi,cts.coloc$adjCellType),1)
    ymax.epi <- tail(intersect(epi,cts.coloc$adjCellType),1)  
    xmin.endo <- head(intersect(endo,cts.coloc$CellType),1) 
    xmax.endo <- tail(intersect(endo,cts.coloc$CellType),1) 
    ymin.endo <- head(intersect(endo,cts.coloc$adjCellType),1) 
    ymax.endo <- tail(intersect(endo,cts.coloc$adjCellType),1) 
    xmin.fibro <- head(intersect(fibro,cts.coloc$CellType),1) 
    xmax.fibro <- tail(intersect(fibro,cts.coloc$CellType),1) 
    ymin.fibro <- head(intersect(fibro,cts.coloc$adjCellType),1) 
    ymax.fibro <- tail(intersect(fibro,cts.coloc$adjCellType),1) 
    p <- ggplot(coloc.all, aes(x = CellType, y = adjCellType, color = -log(pvalue + 0.000001))) +
         geom_point(aes(size = perc), alpha = 0.6) +
         scale_color_gradientn(colours = c('skyblue', 'seagreen','gold','tomato'), name = "Significance") +
         labs(title = "", x = "CellType", y = "adjCellType", size = "Percentage") +
         theme_minimal() +
         theme(axis.text = element_text(size = 12, color = "black"),
               axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
               axis.title = element_text(size = 12, color = "black"),
               legend.position = "right") +
         annotate("rect",xmin = xmin.epi, xmax = xmax.epi, ymin = ymin.epi,  ymax = ymax.epi, alpha = .2,fill = "#A0D3FF") +
         annotate("rect",xmin = xmin.endo, xmax = xmax.endo, ymin = ymin.endo,  ymax = ymax.endo, alpha = .2,fill = "#FF8AD0") + 
         annotate("rect",xmin = xmin.fibro, xmax = xmax.fibro, ymin = ymin.fibro,  ymax = ymax.fibro, alpha = .2,fill = "#A886F6") 
      
    ggsave(file.path(outdir, paste0(slide,"_D",domain,"_Spots_colocalization_celltype_pairs_significant_DotPlot.pdf")), p, width = 10, height = 10)
    return(p)
}