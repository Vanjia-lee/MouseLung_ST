#'This function is used to calculated the distant based on the mapped single cells, 
#'and to filtered the adjacent cell types 
#'cytospace: The cytospace result with assigned locations, 
#'which include the single-cell UniqueCID, CellType, SpotID, imagerow, imagecol.
#'distant: The distant selcted for co-localized.eg: 0 for coloc, 1 for one minimal distant unit, 2 for two minimal distant unit. The minimal distant unit represents the minimal distance between two adjacent spots. 
#'
getD <- function(cytospace, distant = c("coloc","adjac")){  
    require(proxy)
    #distance <- ifelse(distant == "coloc", 0, 2)
    cell.locs <- cytospace %>% select(imagerow, imagecol) %>% proxy::dist() %>% as.matrix() 
    dist.scale <-  round(cell.locs / min(cell.locs[cell.locs > 0]), digits = 0) #计算获得几个距离                      
    dist.all <- data.frame(dist.scale) %>% 
             mutate(UCID = rownames(.)) %>% 
             gather(adjUCID, Distance, -UCID) %>%
             filter(Distance <= distant, UCID != adjUCID) %>% arrange(UCID) %>%
             mutate(CellType = cytospace$CellType[match(.$UCID, rownames(cytospace))],
                    adjCellType = cytospace$CellType[match(.$adjUCID, rownames(cytospace))]) %>%
             select(UCID, CellType, adjUCID, adjCellType, Distance)
    return(dist.all)
}

#' The function is used to calculate the percentage of adjacent celltype for each cell type
#'dist.all: The results of getD, which include five columns information: 
#'UCID, CellType, adjUCID, adjCellType, distance
#'celltypes: All cell types included witch was mapped to the spatial domains
#'
cell.stat <- function(dist.all){
    require(data.table)
    dist.cell <- data.table(dist.all)[, .(Freq = .N), by = .(CellType, adjCellType)]
    dist.cell <- dist.cell[, .(adjCellType, Freq, perc = Freq/sum(Freq) * 100), by = CellType]
    dist.cell <- dist.cell %>% dcast(adjCellType ~ CellType, value.var = "perc", fill = 0)
    return(dist.cell)
}

sort_repair <- function(row1) {
            cells <- sort(row1[c(2,4)])
            repair <- paste(cells[1], cells[2], sep = "_")
            return(repair)
        }
#' The function is used to do permutation analysis 1000 times to filter the significant co-localized #'cell type pairs with test pvaule < 0.05.
permutateD <- function(cytospace, dist_lst = NULL, threads = 10, distant = 0,
                       outdir = NULL, slide = NULL, domain = NULL){ 
    require(foreach)
    require(doParallel)
    if(!is.null(threads)){
        threads = threads
    }else{
        threads = 10
    }
    if(!dir.exists(outdir)){
        dir.create(outdir, recursive = T)
    }
    setwd(outdir)
    # permutation analysis
    cl <- makeCluster(threads)
    registerDoParallel(cl)
    permutation_times <- 1000
    bind_fun<-function(a,b){
        dat <- full_join(a,b,by="adjCellType",suffix = c("1","2")) 
        return(dat)
    }
    rep.df <- foreach(i = 1:permutation_times, .combine = bind_fun, 
                      .packages = c('Seurat','tidyverse','proxy','data.table'),
                      .export = c('getD', 'cell.stat','sort_repair')
                     ) %dopar% { 
        map_test <- cytospace
        rownames(map_test) <- sample(rownames(map_test)) # 将细胞随机打乱
        map_test$CellType <- cytospace$CellType[match(rownames(map_test),rownames(cytospace))]        
        d.test <- getD(map_test, distant = distant)
        d.test$reorder.pair <- apply(d.test, 1, FUN = sort_repair)
        d.test$spotID <- map_test$SpotID[match(d.test$UCID, rownames(map_test))]
        dt.unique <- d.test %>% group_by(spotID) %>% dplyr::distinct(reorder.pair, .keep_all = TRUE)
        cts.test <- cell.stat(dt.unique)
        return(cts.test)
    }
    # close clusters resource after permutation
    stopImplicitCluster()
    stopCluster(cl)
    # calculate the pvalue of percentage for adjCellType of each CellType according to the permutation results
    rep.df <- rep.df %>% tibble::column_to_rownames("adjCellType")
    rep.df <- rep.df[,!grepl("adjCellType", colnames(rep.df))] 
    rep.df[is.na(rep.df)] <- 0
    for(cts in names(dist_lst)){
        pos <- paste0("^", cts)
        ctp <- rep.df[, grep(pos, colnames(rep.df))]
        if(!is.null(dist_lst[[cts]])){
            n <- ncol(ctp)
            ctp$Real.perc <- dist_lst[[cts]]$perc[match(rownames(ctp), dist_lst[[cts]]$adjCellType)]
            ctp$pvalue <- apply(ctp, 1, function(x) sum(x[1:n] > x[["Real.perc"]])/1000) # permutation pvalue
            dist_lst[[cts]]$pvalue <- ctp$pvalue[match(dist_lst[[cts]]$adjCellType, rownames(ctp))] 
        }                 
    }
    save(dist_lst, file = paste0(slide,"_D",domain,"_Spots_colocalization_celltype_pairs_results.RData"))                       

    # filter the significant results with pvalue < 0.05                        
    pcts.lst <- lapply(dist_lst, function(x) pvalue <- x[x$pvalue < 0.05,])
    pcts.lst <- imap(pcts.lst, ~ .x %>% mutate(CellType = .y) %>% select(CellType, adjCellType, Freq, perc, pvalue)) 
    co.locs <- bind_rows(pcts.lst)                   
    write.csv(co.locs, paste0(slide,"_D",domain,"_Spots_colocalization_celltype_pairs_permutation1000_results_significant.csv"))                   
    
    for(id in names(dist_lst)){
        colnames(dist_lst[[id]])[2:4] <- paste0(id, ".", colnames(dist_lst[[id]])[2:4])
    }
    dist.df <- reduce(dist_lst,  full_join, by = join_by(adjCellType))
    dist.df[is.na(dist.df)] <- 0
    write.csv(dist.df, paste0(slide,"_D",domain,"_Spots_colocalization_celltype_pairs_permutation1000_results.csv"))
                       
    return(co.locs)                        
}