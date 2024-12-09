library(Seurat)
library(tidyverse)
library(proxy)
library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(cowplot)

rm(list=ls())

# set work directory
projectDir <- '~/project/lung_10xST'
outDir <- file.path(projectDir, "results/adjacent_domains_enrichment")
if(!dir.exists(outDir)) dir.create(outDir, recursive = T)

# load ST object
stDir <- file.path(projectDir, "results", "3.2.Cluster_deconvolution_Tangram","tangram_origQC_training_genes_mode_cells/Clustering_celltype89_Deconvolution")
st_ob <- readRDS(file.path(stDir,"seurat_E125_P0s1_origQC_tangram_89CT_modeCell_trainGene_deconvolution_res1.1.rds"))
st_ob$seurat_clusters <- st_ob$TangramNew89_snn_res.1.1
st_ob$Domain <- paste0('D', st_ob$`TangramNew89_snn_res.1.1`)
st_ob$Domain <- factor(st_ob$Domain, levels =  paste0('D', levels(st_ob$`TangramNew89_snn_res.1.1`)))
st_ob@meta.data <- st_ob@meta.data %>% select(orig.ident, Stage, TangramNew89_snn_res.1.1, Domain)

col2domain <- c("#A5CFE3","#1C79B3","#FF8000","#FFBD6F",'#31A229','#B4DF8C',"#E6407A","#F89B99","#683E9A","#DBB9EC")
names(col2domain) <- levels(st_ob$Domain)
slides <- unique(st_ob$orig.ident)
spot.size = c(7.5,4.3,3,3,2.5,2,1.9,2)
names(spot.size) = slides

# define functions
getD <- function(obj, slide = NULL, domain = NULL, group.by = 'Clusters', min.d = 1){
    coords <- GetTissueCoordinates(obj, image = slide)
    cellsx <- Cells(obj)[obj@meta.data[[group.by]] == domain & obj$orig.ident == slide]
    if(length(cellsx) > 0){
        coordsx <- coords[cellsx,]
        #coordsy <- coords[!rownames(coords) %in% cellsx,]
        #计算Domain 点与所有点的距离
        cell.locs <- as.matrix(proxy::dist(coordsx, coords))
        # 筛选出周围一圈的spots
        dist.scale <-  round(cell.locs / min(cell.locs[cell.locs > 0]), digits = 0) #计算获得几个距离
        adj.cells <- unique(as.character(unlist(apply(dist.scale, 1, function(x) names(which(x==min.d))))))
        # 保存周围一圈spots的domain信息                                 
        adj.domain <- obj@meta.data[adj.cells, group.by] 
        #adj.cd <- list(adj.cells, adj.domain) 
        domain.spots <- data.frame(spots = adj.cells, domain = adj.domain)                                 
        return(domain.spots)
    }
}

for(dom in levels(st_ob1$Domain)){
    st_ob <- st_ob1

    # Real adjacent domains percentage
    dm.lst <- list()
    # Calculate the adjacent cells around target dom with 1 standard distance
    for(slide in slides){
        domain.stage <- getD(st_ob, slide = slide, domain = dom, group.by = 'Domain', min.d = 1)
        dm.lst <- c(dm.lst, list(domain.stage))
    }
    save(dm.lst, file = file.path(outDir, paste0(dom, "_adjacent_domain_enrichment_dist2_cellbarcode.RData")))
    dm.df <- do.call(rbind, dm.lst)
    real.adj.dm <- data.frame(table(dm.df$domain))
    colnames(real.adj.dm) <- c("Domain","No.spots")
    real.adj.dm$Real.perc <- real.adj.dm$No.spots/sum(real.adj.dm$No.spots) * 100
    print(real.adj.dm)
    # barplot shows domain percentage
    sadj.dm <- real.adj.dm %>% arrange(desc(Real.perc))
    sadj.dm$Domain <- factor(sadj.dm$Domain, levels = sadj.dm$Domain)
    p <- ggplot(sadj.dm, aes(Domain,weight = Real.perc,fill = Domain)) + 
                geom_bar(position = 'dodge') + 
         scale_fill_manual(values = col2domain) +
         theme_bw() +
         theme(axis.text = element_text(size = 12, colour = "black"), 
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank()
              ) +
         geom_text(aes(x = Domain, y = Real.perc,label = round(Real.perc,digits = 2)), vjust =  -0.2, hjust = 0.5) +
         ylab("Frequency")
    print(p)

    # 置换检验1000次
    # 每张切片的spot随机打乱，然后统计D1周围的spot属于哪些Domain
    rep.lst <- list()
    permutation_times <- 1000
    for(i in 1:permutation_times){
        dm.lst <- list()
        for(slide in slides){
            st_ob <- st_ob1
            coords <- GetTissueCoordinates(st_ob, image = slide)
            neworder <- sample(rownames(coords))
            rownames(st_ob[[slide]]@coordinates) <- neworder
            domain.stage <- getD(st_ob, slide = slide, domain = dom, group.by = 'Domain', min.d = 1)
            dm.lst <- c(dm.lst, list(as.character(domain.stage[[2]])))
        }
        adj.dmp <- data.frame(table(unlist(dm.lst)))
        colnames(adj.dmp) <- c("Domain","No.spots")
        adj.dmp$perc <- adj.dmp$No.spots/sum(adj.dmp$No.spots) * 100
        adj.dmp <- adj.dmp[,c(1,3)]
        rep.lst <- c(rep.lst, list(adj.dmp))
    }
    # extract domain percentage after permutating 1000 times
    pt <- Reduce(cbind, rep.lst) %>% tibble::column_to_rownames('Domain') 
    pt[,grep('Domain', colnames(pt))] <- NULL
    colnames(pt) <- paste0("perc",1:1000)
    pt <- pt[levels(real.adj.dm$Domain),]
    pt.all <- cbind(real.adj.dm$Real.perc, pt)
    #rownames(pt.all) <- real.adj.dm$Domain
    colnames(pt.all)[1] <- "Real.perc"
    # calculate the pvalue
    real.adj.dm$permutation.Pval <- 0
    for(domain in rownames(pt.all)){
        t1 <- sort(as.numeric(pt.all[rownames(pt.all) == domain,]), decreasing = T)
        pval <- which(t1 == pt.all$Real.perc[rownames(pt.all) == domain])/length(t1)
        print(pval)
        real.adj.dm$permutation.Pval[real.adj.dm$Domain == domain] <- pval
    }
    real.adj.dm$adj.Pval <- p.adjust(real.adj.dm$permutation.Pval, method = 'BH')
    real.adj.dm$Domain <- factor(real.adj.dm$Domain, levels = real.adj.dm$Domain[order(-real.adj.dm$Real.perc)])
    real.adj.dm <- real.adj.dm %>% arrange(desc(Real.perc), permutation.Pval)
    write.csv(real.adj.dm, file.path(outDir, paste0(dom, "_adjacent_domain_enrichment_permutation1000_pval.csv")))

    p <- ggplot(real.adj.dm, aes(x = Domain, y = Real.perc, fill = Domain)) + 
         geom_col(width = 0.9, color = 'white', alpha = 0.9) +
         geom_text(data = real.adj.dm %>% filter(Real.perc > 5), aes(x = Domain, y = Real.perc, label = round(Real.perc,digits = 2)), nudge_y = 2.2,
                   color = '#585858', size = 4.5) +
         scale_y_continuous(limits = c(0, round(real.adj.dm$Real.perc,0) + 4)) +
         coord_polar(start = 0) +
         scale_fill_manual(values = col2domain) +
         theme_minimal() +
         theme(axis.title = element_blank(), axis.text.y = element_blank(),
               axis.text.x = element_text(size = 11, face = 'bold', color = ifelse(real.adj.dm$adj.Pval < 0.05, 'red', '#454545')),
               #  legend.position = 'none'
              )
    print(p)
    ggsave(file.path(outDir, paste0(dom, '_adjacent_domain_enrichment_permutation1000_Dist2_pval_signif.pdf')), p, width = 10, height = 8)
    
}                                                      