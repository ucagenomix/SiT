library(Seurat)
library(dplyr)
library(pheatmap)
library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(knitr)
library(xtable)
library(ggrepel)
library(cowplot)
library(metap)
library(formattable)
library(STutility)
library(Nebulosa)
library(DoubletFinder)
library(slingshot)
library(formattable)
library(EnhancedVolcano)
library(ReactomeGSA)
library(iDEA)
library(AUCell)
library(GSEABase)
library(fgsea)
library(clusterProfiler)
library(enrichplot)
library(reshape2)
library(dittoSeq)
library(STutility)
library(gameofthrones)
library(SPOTlight)
library(patchwork)

custom.pal <- c("grey","#FF9B00","#EC380B")

ggplotColours <- function(n = 6, h = c(0, 360) + 15){
  if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}
col=ggplotColours(n = 12)


# set scale.data slot from data slot
scale_my_data <- function(seurat, assay){
  data_in <- as.matrix(seurat@assays[[assay]]@data)
  mat <- sapply(strsplit(rownames(data_in), "\\.\\."), `[`, 1)
  
  all.genes <- unique(mat)
  all.genes <- as.data.frame(all.genes)
  
  print(paste("Unique features =", dim(all.genes)[1], "out of", dim(data_in)[1], "total features", sep=" "))
  
  colnames(all.genes)<-c("geneId")
  rownames(all.genes) <- all.genes$geneId
  
  output <- data_in
  
  # for all genes
  for(i in 1:dim(all.genes)[1]){
    
    #print(paste("i",i,sep="="))
    x <- which(mat == all.genes[i,])
    m <- mean(unlist(data_in[x,]))
    sd <- sd(unlist(data_in[x,]))
    
    # for all gene isoforms
    for(j in 1:length(x)){
      #print(paste("j",j,length(x),sep="="))
      
      # fo all cells
      for(k in 1:dim(data_in)[2]){
        #print(paste("k",k,sep="="))
        output[x[j],k] <- (data_in[x[j],k] - m)/sd
      }
    }
  }
  
  #seurat@assays[[assay]]@scale.data <- output
  seurat <- SetAssayData(seurat, slot = "scale.data", output, assay = assay)
  return(seurat)
}

## remove the x-axis text and tick
## plot.margin to adjust the white space between each plot.
## ... pass any arguments to VlnPlot in Seurat
modify_vlnplot<- function(obj, 
                          feature, 
                          pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  p<- VlnPlot(obj, features = feature, pt.size = pt.size, ... )  + 
    xlab("") + ylab(feature) + ggtitle("") + 
    theme(legend.position = "none", 
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          axis.title.y = element_text(size = rel(1), angle = 0), 
          axis.text.y = element_text(size = rel(1)), 
          plot.margin = plot.margin ) 
  return(p)
}

## extract the max value of the y axis
extract_max<- function(p){
  ymax<- max(ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)
  return(ceiling(ymax))
}


## main function
StackedVlnPlot<- function(obj, features,
                          pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  
  plot_list<- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))
  
  # Add back x-axis title to bottom plot. patchwork is going to support this?
  plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
    theme(axis.text.x=element_text(), axis.ticks.x = element_line())
  
  # change the y-axis tick to only max value 
  ymaxs<- purrr::map_dbl(plot_list, extract_max)
  plot_list<- purrr::map2(plot_list, ymaxs, function(x,y) x + 
                            scale_y_continuous(breaks = c(y)) + 
                            expand_limits(y = y))
  
  p<- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
  return(p)
}

