# spatialfeatureplot for multi slides with common legend

MultiSpatialFeaturePlot <- function(st_ob = st_ob, 
                                    gene = gene, 
                                    slides = NULL,
                                    plotTitle = NULL, 
                                    slot = "data", 
                                    image.alpha = 0.2,
                                    alpha = c(0.1, 1),
                                    stroke = 0,
                                    min.cutoff = NULL, 
                                    max.cutoff = NULL,
                                    #breaks = NULL,
                                    cols = colorRampPalette(colors = rev(x = brewer.pal(n = 11, name = "Spectral")))(n = 100),
                                    nrow = 1,
                                    addsubtitle = FALSE, subtitle.size = 12,
                                    common_legend = T){
    if (!is.null(min.cutoff)) {
        min_scale <- min.cutoff
    } else {
        min_scale <- min(GetAssayData(st_ob, assay = DefaultAssay(st_ob), slot = slot)[gene,]) 
    }
    if (!is.null(max.cutoff)) {
        max_scale <- max.cutoff
    } else {
        max_scale <- max(GetAssayData(st_ob, assay = DefaultAssay(st_ob), slot = slot)[gene,])
    }
    
    slides.all = c("E125_S2","E135_S3","E145_S4","E155_S4","E165_S2","E175_S2","E185_S3","P0_S1")
    spot.size = c(7.5,4.3,3,3,2.5,2,1.9,2)
    names(spot.size) = slides.all
    if(is.null(slides)){
        slides = slides.all
    }else{
        slides = slides
    }     
    gsp_lst = list()
    for(slide in slides){
        gsp <- suppressMessages(SpatialFeaturePlot(st_ob, features = gene, images = slide, image.alpha = image.alpha,
                                                   pt.size.factor = as.numeric(spot.size[[slide]]),
                                                   slot = slot, alpha = alpha, stroke = stroke, combine = FALSE
                                                  ))
        
        gsp <- suppressMessages(gsp[[1]] + 
                                scale_fill_gradientn(limits = c(min_scale,max_scale), colours = cols) +
                                scale_alpha(range = alpha) + guides(alpha = "none") +
                                theme(plot.background = element_rect(color = "white"),
                                    panel.border = element_blank())
                                )
        if (addsubtitle) {
            gsp <- suppressMessages(gsp + ggtitle(slide) + theme(plot.title = element_text(size = subtitle.size, hjust = 0.5)))
        } else {
            gsp <- suppressMessages(gsp + theme(plot.title = element_blank(), axis.title.x = element_blank(),axis.title.y = element_blank()))
        }
        gsp_lst = c(gsp_lst, list(gsp))
    }

    gsp.all <- do.call(ggpubr::ggarrange, c(gsp_lst, list(nrow = nrow, common.legend = common_legend, legend = "right", align = "hv")))                 
    plot <- gsp.all + ggtitle(plotTitle) + theme(plot.title = element_text(hjust = 0.5, size = 35, margin = margin(20,0,20,0)),panel.border = element_blank())
    #return(gsp_lst)
    return(plot)
}