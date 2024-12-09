library(Seurat)
library(dplyr)
library(ggplot2)
library(Matrix)
library(htmlwidgets)
library(plotly)
library(FNN)
plan("multicore", workers = 4)
options(future.globals.maxSize = 100000 * 1024^2)

rm(list=ls());gc()

projectDir <- "~/project/lung_10xST/"
sctdir <- file.path(projectDir, "results/4.Development_lineage_construction/single_stage_SCT_cluster")
outdir <- file.path(sctdir, "TOME_result/harmony_SCT")
if(!dir.exists(outdir)) dir.create(outdir, recursive = T)

slides = c("E125_S2","E135_S3","E145_S4","E155_S4","E165_S2","E175_S2","E185_S3","P0_S1")
time_point <-  c(paste0("E",seq(125,185,10)), "P0")
resols <- c(1.2, 1.2, 1.2, 1.3, 1.9, 1.8, 1.2, 1.9) 
names(resols) <- slides
spot.size <- c(7,4,3,3,2.5,2,2,2)
names(spot.size) <- slides

source("~/data/tome_code/tome_code/help_code/help_code.R")

# running Knn to find ancestor
for(id in seq(length(time_point)-1)){
    slide1 <- slides[id]; slide2 <- slides[id + 1]
    res1 <- resols[slide1]; res2 <- resols[slide2]
    time_1 <- time_point[id]; time_2 <- time_point[id + 1]
    print(time_1); print(time_2)
    
    emb = readRDS(paste0(harmonyDir,"/", time_1, "_", time_2, "_umap3.rds"))
    emb = data.frame(emb)
    print(dim(emb))
    
    anno1 <- readRDS(paste0(sctdir, "/", slide1, "/seurat_", slide1, "_SCT_res", res1, "_dim30_cluster.rds"))
    anno1$Anno = anno1$cell_state
    anno1 = anno1[[]][,c("Stage", "Anno")]
    anno1$day = "pre"
    rownames(anno1) = paste0(rownames(anno1), "_1")
    print(head(anno1, 3))
    
    anno2 <- readRDS(paste0(sctdir, "/", slide2, "/seurat_", slide2, "_SCT_res", res2, "_dim30_cluster.rds"))
    anno2$Anno = anno2$cell_state
    anno2 = anno2[[]][,c("Stage", "Anno")]
    anno2$day = "nex"
    rownames(anno2) = paste0(rownames(anno2), "_2")
    print(head(anno2, 3))
    
    anno = rbind(anno1, anno2)
    if(nrow(emb) != nrow(anno)){
        print("Error!")
    }
    anno = anno[rownames(emb),]
    res = createLineage_Knn(emb, anno,  k_neigh = 20)
    saveRDS(res, paste0(outdir, "/", time_1, "_", time_2, "_Knn_umap.rds"))
}

# performing a permutation (shuffling the Anno labels for each connection) by 1,000 times
for(id in seq(length(time_point)-1)){
    slide1 <- slides[id]; slide2 <- slides[id + 1]
    res1 <- resols[slide1]; res2 <- resols[slide2]
    time_1 <- time_point[id]; time_2 <- time_point[id + 1]
    print(time_1); print(time_2)
    
    emb = readRDS(paste0(harmonyDir,"/", time_1, "_", time_2, "_umap3.rds"))
    emb = data.frame(emb)
    print(dim(emb))
    
    anno1 <- readRDS(paste0(sctdir, "/", slide1, "/seurat_", slide1, "_SCT_res", res1, "_dim30_cluster.rds"))
    anno1$Anno = anno1$cell_state
    anno1 = anno1[[]][,c("Stage", "Anno")]
    anno1$day = "pre"
    rownames(anno1) = paste0(rownames(anno1), "_1")
    print(head(anno1, 3))
    anno2 <- readRDS(paste0(sctdir, "/", slide2, "/seurat_", slide2, "_SCT_res", res2, "_dim30_cluster.rds"))
    anno2$Anno = anno2$cell_state
    anno2 = anno2[[]][,c("Stage", "Anno")]
    anno2$day = "nex"
    rownames(anno2) = paste0(rownames(anno2), "_2")
    print(head(anno2, 3))
    
    permutation_times = 1000
    k_neigh = 20
    
    res = list()
    for(rep_i in 1:permutation_times){
    
        anno1$state = anno1$Anno[sample(1:nrow(anno1))]
        anno2$state = anno2$Anno[sample(1:nrow(anno2))]
    
        anno = rbind(anno1, anno2)
        if(nrow(emb) != nrow(anno)){
            print("Error!")
            print(xxx)
        }
        pd = anno[rownames(emb),]
    
        emb_sub = emb
        pd_sub = pd
    
        irlba_pca_res_1 <- emb_sub[as.vector(pd_sub$day)=="pre",]
        irlba_pca_res_2 <- emb_sub[as.vector(pd_sub$day)=="nex",]
        pd_sub1 <- pd_sub[pd_sub$day == "pre",]
        pd_sub2 <- pd_sub[pd_sub$day == "nex",]
    
        pre_state_min = min(table(as.vector(pd_sub1$state)))
    
        if (pre_state_min < k_neigh & pre_state_min >= 3){
            k_neigh = pre_state_min
            print(k_neigh)
        }
    
        if (pre_state_min < 3){
            next
        }
    
        neighbors <- get.knnx(irlba_pca_res_1, irlba_pca_res_2, k = k_neigh)$nn.index
    
        tmp1 <- matrix(NA,nrow(neighbors),ncol(neighbors))
        for(i in 1:k_neigh){
            tmp1[,i] <- as.vector(pd_sub1$state)[neighbors[,i]]
        }
        state1 <- names(table(as.vector(pd_sub1$state)))
        state2 <- names(table(as.vector(pd_sub2$state)))
    
        tmp2 <- matrix(NA,length(state2),length(state1))
        for(i in 1:length(state2)){
            x <- c(tmp1[as.vector(pd_sub2$state)==state2[i],])
            for(j in 1:length(state1)){
                tmp2[i,j] <- sum(x==state1[j])
            }
        }
        tmp2 <- tmp2/apply(tmp2,1,sum)
        tmp2 <- data.frame(tmp2)
        row.names(tmp2) = state2
        names(tmp2) = state1
    
        res[[rep_i]] = tmp2
    
    }
    saveRDS(res, paste0(outdir, "/", time_1, "_", time_2, "_Knn_umap_permutation.rds"))
}

x = list()
z = NULL
for(cnt in 1:(length(time_point)-1)){
    time_i = time_point[cnt]
    time_j = time_point[cnt+1]
    
    dat = readRDS(paste0(outdir, "/", time_i, "_", time_j, "_Knn_umap_permutation.rds"))
    
    permutation_times = 1000
    y = NULL
    
    for(i in 1:permutation_times){
        y = c(y, as.vector(as.matrix(dat[[i]])))
    }
    
    x[[cnt]] = y
    z = c(z, y)
    
    print(paste0(time_i, ":", sum(y > 0.2)/length(y)))
}

print(sum(z >= 0.2)/length(z))

dat = data.frame(edge_weights_by_permutation = z)
p <- ggplot(dat, aes(x=edge_weights_by_permutation)) +
    geom_histogram(position="identity", alpha=0.5, binwidth=0.01) + 
    geom_vline(xintercept = 0.2, colour = "red") +
    theme_classic(base_size = 15)
p

# summarize results
library(reshape2)
library(scales)
library(ggplot2)
library(dplyr)
library(gplots)
library(viridis)

time_i = 6
dat = readRDS(paste0(outdir, "/",time_point[time_i],"_",time_point[time_i+1],"_Knn_umap.rds"))

replication_times=500
res_median_umap = list()

for(time_i in 1:(length(time_point)-1)){
  print(paste0(time_point[time_i], ":", time_point[time_i+1]))  
  dat = readRDS(paste0(outdir, "/",time_point[time_i],"_",time_point[time_i+1],"_Knn_umap.rds"))
  state_1 = row.names(dat[[1]])
  state_1 = paste0(time_point[time_i+1], ":", gsub(paste0(time_point[time_i+1], ":"), "", state_1))
  state_2 = names(dat[[1]])
  state_2 = paste0(time_point[time_i], ":", gsub(paste0(time_point[time_i], ":"), "", state_2))
  tmp_1 = matrix(NA,nrow(dat[[1]]),ncol(dat[[1]]))
  for(i in 1:nrow(dat[[1]])){
    for(j in 1:ncol(dat[[1]])){
      xx = NULL
      for(k in 1:replication_times){
        xx = c(xx, dat[[k]][i,j])
      }
      tmp_1[i,j] = median(xx[!is.na(xx)])
    }
  }
  tmp_1 = data.frame(tmp_1)
  row.names(tmp_1) = state_1
  names(tmp_1) = state_2
  res_median_umap[[time_i]] = tmp_1
}

dat = NULL
for(i in 1:length(res_median_umap)){
  print(time_point[i])
  dat = rbind(dat, melt(as.matrix(res_median_umap[[i]])))
}

dat = data.frame(dat)
names(dat) = c("nex", "pre", "prob")
head(dat)

dat$pre_time = unlist(lapply(as.vector(dat$pre), function(x) strsplit(x,"[:]")[[1]][1]))
dat$pre_cell = unlist(lapply(as.vector(dat$pre), function(x) strsplit(x,"[:]")[[1]][2]))
dat$nex_time = unlist(lapply(as.vector(dat$nex), function(x) strsplit(x,"[:]")[[1]][1]))
dat$nex_cell = unlist(lapply(as.vector(dat$nex), function(x) strsplit(x,"[:]")[[1]][2]))

saveRDS(dat, paste0(outdir, "/edge_all1.rds"))

### Here we use “cell state” to mean an annotated cluster at a given stage. 

print(paste0("how many edges: ", nrow(dat)))
print(paste0("how many edges (> 0): ", nrow(dat[dat$prob>0,])))
print(paste0("how many edges (> 0.1): ", nrow(dat[dat$prob>=0.1,])))
print(paste0("how many edges (> 0.15): ", nrow(dat[dat$prob>=0.15,])))
print(paste0("how many edges (> 0.2): ", nrow(dat[dat$prob>=0.2,])))
print(paste0("how many edges (> 0.2): ", nrow(dat[dat$prob>=0.3,])))
print(paste0("how many edges (> 0.5): ", nrow(dat[dat$prob>=0.5,])))
print(paste0("how many edges (> 0.7): ", nrow(dat[dat$prob>=0.7,])))
print(paste0("how many edges (> 0.8): ", nrow(dat[dat$prob>=0.8,])))
print(paste0("how many nodes: ", length(unique(c(as.vector(dat$pre), as.vector(dat$nex))))))
print(paste0("how many cell types: ", length(unique(c(as.vector(dat$pre_cell), as.vector(dat$nex_cell))))))

res = dat[dat$prob > 0,]
res = res[,c(2,1,3:7)]
source.max <- res %>% group_by(pre) %>% slice_max(prob, n = 1)
target.max <- res %>% group_by(nex) %>% slice_max(prob, n = 1)
res.max <- source.max %>% bind_rows(target.max) %>% arrange(pre) %>% unique()
head(res.max)

write.table(res.max[,1:3], file.path(outdir, paste0("edge_prob_max_source_target.txt")), row.names = F, quote = F, sep = "\t")          res.max <- read.csv(file.path(outdir, paste0("edge_prob_max_source_target.txt")),sep = "\t")
head(res.max) 

fib.group <- c("E125:Fibroblasts","E135:Fibroblasts","E135:Myocardium",
               "E145:Fibroblasts","E145:Trachea","E145:Myocardium",
               "E155:Trachea","E155:Myocardium & Fibroblasts",
               "E165:Fibroblasts","E165:Myocardium",
               "E175:Pulmonary artery","E175:Pulmonary vessel (Vwf+)",
               "E185:Fibroblasts","E185:Pulmonary vein",
               "P0:Fibroblasts","P0:Trachea","P0:Pulmonary vein")
mes.group <- c("E125:Pulmonary mesenchyme","E135:Pulmonary mesenchyme",
               "E145:Pulmonary mesenchyme (Cdh5+Aplnr+)","E145:Pulmonary mesenchyme (Tbx2+Wnt2+)",
               "E155:Pulmonary mesenchyme (Aplnr+)","E155:Pulmonary mesenchyme (Tbx2+Wnt2+)",
               "E165:Pulmonary mesenchyme (Aplnr+)","E165:Pulmonary mesenchyme (Upk3b+)",
               "E165:Pulmonary mesenchyme (Tgfbi+)",
               "E175:Pulmonary mesenchyme (Upk3b+)","E175:Capillary","E175:Matrix fibroblasts",
               "E185:Mesothelium","E185:Fibroblasts & SM","E185:Pulmonary mesenchyme",
               "P0:Mesothelium")
dist.group <- c("E125:Distal airways","E135:Distal airways","E145:Distal airways","E155:Distal airways",
                "E165:Alveolar saccules (Hopx+Ager+)","E165:Alveolar saccules (Mfap4+Sftpa1+)",
                "E165:Alveolar saccules (Actb+Sftpc+Upk3b+)",
                "E175:Alveolar saccules (Actb+Akap5+Tgfbi+)",
                "E175:Alveolar saccules (Ager+Hopx+)",
                "E185:Alveolar saccules (Tmsb10+Hopx+Ager+)",
                "E185:Alveolar saccules (Sftpc+Cxcl15+)",
                "P0:Alveolar saccules (Malat1+Retnla+)",
                "P0:Alveolar saccules (Sftpc+Ager+)",
                "P0:Alveolar saccules (Sftpc+Sftpa1+)",
                "P0:Alveolar saccules (Cdh5+Sftpc+)"
                )
prox.group <- c("E125:Proximal airways","E135:Proximal airways",
                "E145:Pulmonary airways","E155:Pulmonary airways","E165:Pulmonary airways",
                "E175:Small airways","E175:Large airways",
                "E185:Small airways","E185:Large airways",
                "P0:Small airways","P0:Medium airways","P0:Large airways",
                "E185:Bronchi")
immune.group <- c("E175:Immune","E185:Immune","P0:Immune")

annotations <- c(paste0("E125:",c("Fibroblasts", 
                                  "Pulmonary mesenchyme",
                                  "Distal airways",
                                  "Proximal airways")),
              paste0("E135:",c("Fibroblasts", 
                               "Myocardium",
                               "Pulmonary mesenchyme", 
                               "Distal airways", 
                               "Proximal airways")),
              paste0("E145:",c("Trachea",
                               "Fibroblasts",
                               "Myocardium",
                               "Pulmonary mesenchyme (Cdh5+Aplnr+)",
                               "Pulmonary mesenchyme (Tbx2+Wnt2+)",
                               "Distal airways",
                               "Pulmonary airways")),
              paste0("E155:",c("Trachea",
                               "Myocardium & Fibroblasts",
                               "Pulmonary mesenchyme (Aplnr+)",
                               "Pulmonary mesenchyme (Tbx2+Wnt2+)",
                               "Distal airways",
                               "Pulmonary airways")),
              paste0("E165:",c("Fibroblasts",
                               "Myocardium",
                               "Pulmonary mesenchyme (Aplnr+)",
                               "Pulmonary mesenchyme (Upk3b+)",
                               "Pulmonary mesenchyme (Tgfbi+)",
                               "Alveolar saccules (Hopx+Ager+)",
                               "Alveolar saccules (Mfap4+Sftpa1+)",
                               "Alveolar saccules (Actb+Sftpc+Upk3b+)",
                                "Pulmonary airways")),
              paste0("E175:",c("Pulmonary artery",
                               "Pulmonary vessel (Vwf+)",
                               "Pulmonary mesenchyme (Upk3b+)", 
                               "Capillary",
                                "Alveolar saccules (Actb+Akap5+Tgfbi+)",
                               "Alveolar saccules (Ager+Hopx+)",
                               "Immune",
                               "Matrix fibroblasts",
                               "Small airways",
                               "Large airways"
                              )),
              paste0("E185:",c("Fibroblasts",
                               "Pulmonary vein",
                               "Fibroblasts & SM",
                               "Mesothelium",
                               "Pulmonary mesenchyme",
                               "Alveolar saccules (Tmsb10+Hopx+Ager+)",
                               "Alveolar saccules (Sftpc+Cxcl15+)",
                               "Bronchi",
                               "Immune",
                               "Small airways",
                               "Large airways")),
              paste0("P0:",c("Fibroblasts",
                             "Pulmonary vein",
                             "Trachea",
                             "Mesothelium",
                             "Alveolar saccules (Malat1+Retnla+)",
                             "Alveolar saccules (Sftpc+Ager+)",
                             "Alveolar saccules (Sftpc+Sftpa1+)",
                             "Alveolar saccules (Cdh5+Sftpc+)",
                             "Immune",
                             "Small airways",
                             "Large airways",
                             "Medium airways"))         
             )
group <- ifelse(annotations %in% fib.group, 1, 
                ifelse(annotations %in% mes.group, 2,
                       ifelse(annotations %in% dist.group, 3,
                              ifelse(annotations %in% prox.group, 4, 5))))
df <- data.frame(annotations, group)
write.table(df,file.path(outdir, "celltype_group_max.txt"), row.names = F, sep = "\t",quote =F)

colors = c("#1B9B52","#BDAED4","#1DBEC2","#FDC087","#87A8D8")
color_map = data.frame(group = 1:5, colors = colors)

d.anno = df
d.anno$id = seq(nrow(d.anno))
d.anno$colors = color_map$colors[match(d.anno$group, color_map$group)]
d.anno = d.anno[order(d.anno$id),]
head(d.anno)

edge_color <- data.frame(stage1=d.anno$annotations, colr.l=d.anno$colors)
snak <- data.frame(stage1=res.max$pre,stage2=res.max$nex,prob=res.max$prob)
snak$IDsource <- match(snak$stage1, d.anno$annotations) - 1 
snak$IDtarget <- match(snak$stage2, d.anno$annotations) - 1
snak$label <- as.character(round(snak$prob, 3))
snak$color <- d.anno$colors[match(snak$stage1, d.anno$annotations)]
snak$weight <- 100^(snak$prob)
str(snak)
head(snak)

p <- plot_ly(
  type = 'sankey', orientation = 'h', arrangement = "freeform",
  
  #node
  node = list(
      label = d.anno$annotations, 
      color = d.anno$colors,
      # x = d.anno$pos.x,
      pad = 50, 
      thickness = 20,
      line = list(color = 'white', width = 0.5)
  ),
  
  #id for link
  link = list(
      source = snak$IDsource, 
      target = snak$IDtarget,
      value = snak$prob, 
      #label = snak$label,
      color = snak$color
  )
)

p <- p%>% layout(
  title = "Lung development trajectory weighted by prop",
  font = list(
      size = 15
  )
)

htmlwidgets::saveWidget(as_widget(p), paste0(outdir,"/snakey.KNN.pca.prob.max.nodePosition",".html"))

p <- plot_ly(
  type = 'sankey', orientation = 'h', arrangement = "freeform",
  
  #node
  node = list(
      #label = d.anno$annotations, 
      color = d.anno$colors,
      # x = d.anno$pos.x,
      pad = 50, 
      thickness = 20,
      line = list(color = 'white', width = 0.5)
  ),
  
  #id for link
  link = list(
      source = snak$IDsource, 
      target = snak$IDtarget,
      value = snak$prob, 
      #label = snak$label,
      color = snak$color
  )
)

p <- p%>% layout(
  title = "Lung development trajectory weighted by prop",
  font = list(
      size = 15
  )
)

htmlwidgets::saveWidget(as_widget(p), paste0(outdir,"/snakey.KNN.pca.prob.max.nodePosition.nolabel",".html"))

dat <- readRDS(paste0(outdir, "/edge_all1.rds"))
res = dat[dat$prob > 0,]
res = res[,c(2,1,3:7)]
head(res)
res[res$pre == "E165:Alveolar saccules (Mfap4+Sftpa1+)",]
res.max[res.max$pre == "E165:Alveolar saccules (Mfap4+Sftpa1+)",]
res.max1 <- res.max
res.max1[res.max1$pre == "E165:Alveolar saccules (Mfap4+Sftpa1+)","nex"] = "E175:Capillary"
res.max1[res.max1$pre == "E165:Alveolar saccules (Mfap4+Sftpa1+)","prob"] = res[res$pre == "E165:Alveolar saccules (Mfap4+Sftpa1+)" & res$nex == "E175:Capillary","prob"]
res.max1[res.max1$pre == "E165:Alveolar saccules (Mfap4+Sftpa1+)",] 
d.anno1 = df
d.anno1$id = seq(nrow(d.anno1))
d.anno1$colors = color_map$colors[match(d.anno1$group,color_map$group)]
d.anno1 = d.anno1[order(d.anno1$id),]
edge_color1 <- data.frame(stage1=d.anno1$annotations, colr.l=d.anno1$colors)
snak1 <- data.frame(stage1 = res.max1$pre, stage2 = res.max1$nex, prob = res.max1$prob)
snak1$IDsource <- match(snak1$stage1, d.anno1$annotations) - 1 
snak1$IDtarget <- match(snak1$stage2, d.anno1$annotations) - 1
snak1$label <- as.character(round(snak1$prob, 3))
snak1$color <- d.anno1$colors[match(snak1$stage1, d.anno1$annotations)]
snak1$weight <- 100^(snak1$prob)
p <- plot_ly(
  type = 'sankey', orientation = 'h', arrangement = "freeform",
  
  #node
  node = list(
      label = d.anno1$annotations, 
      color = d.anno1$colors,
      # x = d.anno$pos.x,
      pad = 50, 
      thickness = 20,
      line = list(color = 'white', width = 0.5)
  ),
  
  #id for link
  link = list(
      source = snak1$IDsource, 
      target = snak1$IDtarget,
      value = snak1$prob, 
      label = snak1$label,
      color = snak1$color
  )
)

p <- p%>% layout(
  title = "Lung development trajectory weighted by prop",
  font = list(
      size = 15
  )
)

htmlwidgets::saveWidget(as_widget(p), paste0(outdir,"/snakey.KNN.pca.prob.all.nodePosition.Immune",".html"))
p <- plot_ly(
  type = 'sankey', orientation = 'h', arrangement = "freeform",
  
  #node
  node = list(
      #label = d.anno1$annotations, 
      color = d.anno1$colors,
      # x = d.anno$pos.x,
      pad = 50, 
      thickness = 20,
      line = list(color = 'white', width = 0.5)
  ),
  
  #id for link
  link = list(
      source = snak1$IDsource, 
      target = snak1$IDtarget,
      value = snak1$prob, 
      label = snak1$label,
      color = snak1$color
  )
)

p <- p%>% layout(
  title = "Lung development trajectory weighted by prop",
  font = list(
      size = 15
  )
)

htmlwidgets::saveWidget(as_widget(p), paste0(outdir,"/snakey.KNN.pca.prob.all.nodePosition.Immune.nolabel",".html"))                             
                             