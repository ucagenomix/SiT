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

custom.pal <- c("grey","#FF9B00","#EC380B")

markers <- c("Gad1","Gad2","Maf","Sst","Foxp1","Isl1","Htr3a","Dlx6os1","Dlx2","Nrxn3","Synpr","Slc17a6","Fabp7","Reln","Eomes","Vim","Top2a","Ube2c","Mki67","Ccnd2","Olig1","Pecam1","Igfbp7","Id2","Neurod6","Neurod2")
markers2 <- c("Olig1","Gad1","Gad2","Slc17a6","Slc17a7","Snap25","Gfap","Id2","Neurod1","Neurod2","Neurod6","Cd63","Rora","Calb2","Nr2f2","Htr3a","Dcx","Pax6","Eomes","Vim","S100b","Rbfox1","Ccnd2","Synpr","Mef2c","Pde1c")
precursors <- c("Mef2c","Erbb4","Plcxd3","Tspan7","Satb1","Synpr","Reln","Mpped1","Id2","Top2a","Cenpf","Ube2c","Olig1","Calb2","Pecam1","Cldn5","Igfbp7","Gad1","Gad2","Nrgn","Rora","Unc5c","Mdk","Neurod2","Slc17a6","Dlx1","Pbx3","Htr3a","Ckb","Cd63","Cd9","Slc6a5")
sub.markers <- c("Npy","Calb2","Nr2f2","Th","Cck","Reln","Pvalb","Nos1","Htr3a","Ccnd2","Id2","Synpr","Mef2c","Reln","Neurod6")

ggplotColours <- function(n = 6, h = c(0, 360) + 15){
  if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}
col=ggplotColours(n = 12)

setwd("/data/10x_data/10x_visium/")

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
