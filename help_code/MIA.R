## download from https://github.com/IkjotSidhu/Spatial-HP-Skin/
## MULTI-MODAL INTERSECTION ANALYSIS

MIA_ENRICH <- function(stlist,sclist,total){
  overlap <- length(intersect(stlist$gene,sclist$gene))
  C <- length(sclist$gene)
  D <- length(stlist$gene)
  
  # ENRICHMENT CALCULATION
  e <- -log10(phyper(overlap,C, total - C, D, lower.tail = FALSE))
  return(e)
}
MIA_DEPLETE <- function(stlist,sclist){
  overlap <- length(intersect(stlist$gene,sclist$gene))
  C <- length(sclist$gene)
  D <- length(stlist$gene)
  
  # ENRICHMENT CALCULATION
  d <- -log10(1- (phyper(overlap,C, total - C, D, lower.tail = FALSE)))
  return(d)
}
MIA <- function(total_genes,single_cell.markers,spatial.markers)
{
  #D.SCORES <- c()
  # Perform this operation for every cell type
  #single_cell.markers <- FindAllMarkers(single_cell,assay = assay_use,logfc.threshold = 0.25,return.thresh = p_val_adj < 0.1)
  #spatial.markers <- FindAllMarkers(spatial_data,assay = assay_use,logfc.threshold = 0.25,return.thresh = p_val_adj < 0.1)
  cell.types <- single_cell.markers %>% dplyr::select(cluster) %>% unique() %>% as.list()
  spatial.regions <- spatial.markers %>% dplyr::select(cluster) %>% unique() %>% as.list()
  E.SCORES <- data.frame(spatial.regions)
  for(i in cell.types){
    for (x in i){
      e_list <- c()
    #list.append(e_list,i)
      for(y in spatial.regions){
        for(z in y){
          single_cell <- single_cell.markers %>% filter(cluster==x)
          spatial_data <- spatial.markers %>% filter(cluster==z)
          e <- MIA_ENRICH(single_cell,spatial_data,total = total_genes)
          #d <- MIA_DEPLETE(single_cell,spatial_data,total = total_genes)
          e_list <- c(e_list,e)
        }
        #D.SCORES <- append(D.SCORES,d)
      }
      E.SCORES[paste(x)] <- e_list
    #E.SCORES <- append(E.SCORES,e)
    }
  }
  #e.data <- data.frame("GA"=E.SCORES[1],"ER"=E.SCORES[2],"C0L17A1+"=E.SCORES[3])
  #d.data <- data.frame("GA"=D.SCORES[1],"ER"=D.SCORES[2],"C0L17A1+"=D.SCORES[3])
  #res <- list(E.SCORES,e.data)
  return(E.SCORES)
}
