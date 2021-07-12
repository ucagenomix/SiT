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

imprinted <- c('Znf264','Gpr1','Zdbf2','Mcts2','Mir296','Mir298','Sfmbt2','Gatm','H13','Blcap','Nnat','Nespas','Gnas','Lin28a','Magi2','Il6','Mkrn1-ps1','Mir335','Peg10','Ppp1r9a','Asb4','Sgce','Tfpi2','Calcr','Mest','Copg2','Klf14','Nap1l5','Dhcr7','LOC101055709','Peg3os','Ano1','Zim1','Usp29','Peg3','Zim3','Zfp264','Magel2','Peg12','Ndn','Ube3a','Pwcr1','Zfp127as','Mkrn3','Snrpn','Ampd3','H19','Igf2as','Igf2','Ins2','Th','Ascl2','Kcnq1ot1','Tssc4','Cd81','Kcnq1','Phlda2','Slc22a18','AF313042','Cdkn1c','Nap1l4','Tnfrsf23','Tnfrsf22','AK155734','Tnfrsf26','Nctc1','Zim2','Snurf','Inpp5f','Gab1','Snx14','Mir184','Ntm','Musd2','Rasgrf1','Hymai','Plagl1','Dcn','Ccdc40','Ddc','Grb10','Commd1','Zrsr1','U2af1-rs1','AF357359','Mir410','Mir154','Mir370','Mir376b','Mir136','Mir134','AF357355','B830012L14Rik','Mir127','Mir380','Mir337','Mir411','Mir431','Gtl2','Dlk1','AF357341','AF357426','AF357425','AF357428','Rian','Dio3','Mirg','Begain','Rtl1','Htr2a','Peg13','Kcnk9','Slc38a4','Slc22a3','Slc22a2','Igf2r','Air','Impact','Rhox5','Tsix','Jpx','Ftx','Zcchc13','Xist')

epigenes <- c( "MKI67","TOP2A",                               # Cycling basals
                "VIM",'CAV1',"FN1",'KIT', "DLK2", "KRT5", "TP63",  # Basals
                "KRT6B","FABP5",                               # Suprabasals
                "SCGB1A1","KRT13","SPRR3",                     # Secretory
                "MUC5AC", "MUC5B", "BPIFA1",'BPIFB1',                   # Goblets
                "FOXJ1","TPPP3",                               # Multiciliated
                "CDC20B", "DEUP1", "FOXN4",'CCNO','MCIDAS',                  # Deuterosomal
                "CFTR", "ASCL3"                      # ionocytes
)

immgenes <- c("CD79A","CD3E","CD3D","NKG7","KLRB1","GZMA","CD4","S100A8","S100A9","BCL2A1","G0S2","FCN1","LYZ","MAF","IRF8","CLEC9A","JCHAIN","IGKC","IRF7","TPSB2","PF4","PPBP","GP9")

cc.genes <- readLines(con = "/data/10x_data/cell_cycle_vignette_files/regev_lab_cell_cycle_genes.txt")
# We can segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes[1:43]
g2m.genes <- cc.genes[44:97]

interferon <- c("IFNB1","IFNL1","IFNL2","IFNL3","IFNLR1","IFNLR2","IFI27","IFITM3","IFI6","IFIT1","MX1","ISG15","CCL2","CXCL9","CXCL10","CXCL11","CXCL16","IL1A","IL1B","IL1RN","IL6","IL10","TNF")


library(patchwork)

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

