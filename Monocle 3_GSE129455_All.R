#### Reference ####
  ## https://cole-trapnell-lab.github.io/monocle3/
  ## http://cole-trapnell-lab.github.io/monocle-release/docs/#installing-monocle

##### Presetting ######
  rm(list = ls()) # Clean variable
  memory.limit(150000)

#### Version information ####
  # platform       x86_64-w64-mingw32          
  # arch           x86_64                      
  # os             mingw32                     
  # system         x86_64, mingw32             
  # status                                     
  # major          3                           
  # minor          6.3                         
  # year           2020                        
  # month          02                          
  # day            29                          
  # svn rev        77875                       
  # language       R                           
  # version.string R version 3.6.3 (2020-02-29)
  # nickname       Holding the Windsock      

#### Load Packages ####
  library(monocle3)
  library(ggplot2)
  library(magrittr) # need to run every time you start R and want to use %>%
  library(dplyr)
  
  library(reticulate)
  library(Matrix)
  # reticulate::py_install("louvain")
  # py_install("louvain")
  library(AnnotationDbi)
  library(org.Mm.eg.db)
  library(Hmisc)
  

#### Load files ####
  setwd(getwd()) ## Set current working directory
  PathName <- getwd() ## Set output directroy
  FilesName = paste0(PathName,"/")

  RVersion = "20210307AllV1"
  dir.create(paste0(FilesName,RVersion))

  # Build # expression_matrix
  data <- read.csv(paste0(FilesName,"GSE129455_All_Viable_expression.csv"),   
                   header=T,          
                   sep=",")           
  
  row.names(data) <- data[,1]
  
  expression_matrix2 <- data[1:length(data[,1]), 2:length(data[1,])]
  ession_matrix3 <- as(as.matrix(expression_matrix2), "dgCMatrix")

# Bu#### Build monocle obj ####  
  # Build # gene_metadata
  dataG <- as.data.frame(data[,1])
  row.names(dataG) <- data[,1]
  colnames(dataG) <- c("gene_short_name")
  
  # Build # cell_metadata
  dataC2 <- as.data.frame(colnames(expression_matrix2))
  row.names(dataC2) <- dataC2[,1]
  colnames(dataC2) <- c("samples")
  dataC2$id  <- 1:nrow(dataC2)
  
  ## Make the CDS object
  cell_metadata2 <- dataC2
  gene_annotation2 <- dataG
  
  cds <- new_cell_data_set(expression_matrix3,
                           cell_metadata = cell_metadata2,
                       gene_metadata = gene_annotation2)

# cds@rowRanges@elementMetadata@listData$gene_short_name <- toupper(cds@rowRanges@elementMetadata@listData$gene_short_name)

#################### Gene name conversion #################### 
libreNameEnsembl <- as.matrix(cds@rowRanges@elementMetadata@listData$gene_short_name)
GeneNameSymbol <- AnnotationDbi::select(org.Mm.eg.db, keys=GeneNameEnsembl, columns='SYMBOL', keytype='ENSEMBL')

## Clear the repeat
# GeneNameSymbol2=as.matrix(GeneNameSymbol[unique(GeneNameSymbol[,1]),]) #error
GeneNameSymbol2=GeneNameSymbol[!duplicated(GeneNameSymbol$ENSEMBL),]

## Number for null
GNS2_NA <- which(is.na(GeneNameSymbol2[,2])) 
GeneNameSymbol2[GNS2_NA,2] <- GeneNameSymbol2[GNS2_NA,1]

cds@rowRanges@elementMetadata@listData$gene_ENSEMBL <- cds@rowRanges@elementMetadata@listData$gene_short_name
cds@rowRanges@elementMetadata@listData$gene_short_name <- GeneNameSymbol2[,2]
# Main = c("TOP2A")
# Main <- capitalize(tolower(Main))
# Main <- AnnotationDbi::select(org.Mm.eg.db, keys=Main, columns='ENSEMBL', keytype='SYMBOL')


####################
MainCD <- c("Cd4",
          "Cd8a","Cd8b1")

MainMeth <- c("Rbpjl","Nr5a2","Ptf1a","Amy2a5","Amy2a3","Cela2b","Cpa1","Sycn","Ctrc","Pdx1","Sox9","Hnf1b","Krt19")

Main <- c("Cd248",
          "Ptk2",
          "Ptk7",
          "Des",
          "Wnt5a",
          "Il1b",
          "Cxcl15", # IL-8 ENSMUSG00000029375 https://www.ensembl.org/Mus_musculus/Gene/Summary?db=core;g=ENSMUSG00000029375;r=5:90942393-90950926;t=ENSMUST00000031322
          "Mki67",
          "Top2a")


CD248 = c("ENSMUSG00000056481")
CD248 = c("Cd248")

NSUN_Main = c("Nsun1","Nop2","Nsun2","Nsun3","Nsun4","Nsun5","Nsun6","Nsun7")

####################

All_Macrophages_Main = c("Apoe","Saa3","C1qc","C1qa","C1qb","Lyz2")
All_Myeloid_Main = c("Retnla","Ear2","Lyz1","ENSMUSG00000098178","Crip1","Lpl")
All_BCells_Main = c("Cd79a","Ly6d","Ms4a1","Cd79b","Ebf1","H2-DMb2")
All_Neutrophils_Main = c("S100a8","S100a9","G0s2","Il1r2","Retnlg","ENSMUSG00000079597")
All_Fibroblasts = c("Dcn","Col3a1","Col1a1","Sparc","Col1a2","Mgp")
All_TandNKCells = c("Ccl5","Cd3g","Ltb","Cd3d","Cd3e","Gzma")
All_Ductal_Main = c("Clu","Tff2","Tff1","Krt18","Gkn3","Npy")
All_Endothelial_Main = c("Igfbp7","Plvap","Sparc","Igfbp3","Fabp4","Aqp1")
All_Dendritic_Main = c("Ccl5","Fscn1","Ccr7","Ccl22","Tmem123","Ccl17")
All_Acinar_Main = c("Ctrb1","Clps","Prss2","Try5","Reg3b","Reg1")
All_Perivascular_Main = c("Rgs5","Igfbp7","Acta2","Sparcl1","Sparc","Serpine2")
All_EMTlike_Main = c("Cdkn2a","S100a6","Timp1","Tm4sf1","Igfbp4","Lgals1")

####################


### Clustering and classifying your cells
## Pre-process the data

cds = preprocess_cds(cds, num_dim = 100)

png(paste0(FilesName,RVersion,"/",RVersion,"_","PCA.png"), width = 640, height = 360) # ?]?w???X????
plot_pc_variance_explained(cds)

dev.off() # ???????X????



######################## tSNE ###################################

cds <- reduce_dimension(cds,reduction_method = c( "tSNE") , perplexity = 60)

png(paste0(FilesName,RVersion,"/",RVersion,"_","tSNE.png"), width = 640, height = 360) # ?]?w???X????
plot_cells(cds,reduction_method = c( "tSNE"))
dev.off() # ???????X????

# plot_cells(cds, genes=Main)
plot_cells(cds, genes=iCAF_Main,reduction_method = c( "tSNE"))
plot_cells(cds, genes=apCAF_Main,reduction_method = c( "tSNE"))
plot_cells(cds, genes=myCAF_Main,reduction_method = c( "tSNE"))
plot_cells(cds, genes=CD248,reduction_method = c( "tSNE"),cell_size = 1)

mypalette <- colorRampPalette(c("white" , "red"))
plot_cells(cds, genes=CD248,reduction_method = c( "tSNE"),cell_size = 1,label_cell_groups = FALSE)+ theme(legend.position="top")

library(ggplot2)
# http://r-statistics.co/Complete-Ggplot2-Tutorial-Part2-Customizing-Theme-With-R-Code.html
# https://stackoverflow.com/questions/43359050/error-continuous-value-supplied-to-discrete-scale-in-default-data-set-example
# http://www.sthda.com/english/wiki/ggplot2-axis-ticks-a-guide-to-customize-tick-marks-and-labels
# http://www.sthda.com/english/wiki/ggplot2-legend-easy-steps-to-change-the-position-and-the-appearance-of-a-graph-legend-in-r-software
# http://www.sthda.com/english/wiki/ggplot2-title-main-axis-and-legend-titles
# https://bookdown.org/asmundhreinn/r4ds-master/graphics-for-communication.html
# https://gist.github.com/jrnold/79d76e927386298688c9

plot_cells(cds, genes=CD248,reduction_method = c( "tSNE"),cell_size = 1,label_cell_groups = FALSE)+
 theme(axis.text.x = element_text(face="bold", color="#993333", 
                                   size=14, angle=45),
        axis.text.y = element_text(face="bold", color="#993333", 
                                   size=14, angle=45))
#?u gend ?L?tps://www.r-graph-gallery.com/275-add-text-labels-with-ggplot2.html

# Cd248
plot_cells(cds, genes=CD248,reduction_method = c( "tSNE"),cell_size = 1,label_cell_groups = FALSE)+
  theme(axis.text.x = element_text(face="bold",  size=14),
        axis.text.y = element_text(face="bold",size=14),
        axis.line = element_line(colour = "darkblue", size = 2, linetype = "solid"),
        axis.title = element_text(size = rel(1.5),face="bold"),
        plot.title = element_text(color="black", size=20, face="bold.italic",hjust = 0.5,vjust = - 6),
#        panel.background = element_rect(fill = "lightblue"),
        legend.title = element_text(size=12, color = "black", face="bold",inherit.blank=TRUE),
        legend.text = element_text(colour="black", size=12,face="bold"),
        aspect.ratio=1  #square plot
        )+ labs(fill = "log10(Exp)")  # +labs(title="CD248")

#+ theme(legend.position="top")

# MainCD
plot_cells(cds, genes=MainCD,reduction_method = c( "tSNE"),cell_size = 1,label_cell_groups = FALSE)+
  theme(axis.text.x = element_text(face="bold",  size=12),
        axis.text.y = element_text(face="bold",size=12),
        axis.line = element_line(colour = "darkblue", size = 2, linetype = "solid"),
        axis.title = element_text(size = rel(1.5),face="bold"),
        plot.title = element_text(color="black", size=17, face="bold.italic",hjust = 0.5,vjust =-1),
        legend.title = element_text(size=12, color = "black", face="bold"),
        legend.text = element_text(colour="black", size=12,face="bold"),
        aspect.ratio=1  #square plot
  )+labs(title="MainCD")+ labs(fill = "log10(Exp)")

# MainMeth
plot_cells(cds, genes=MainMeth,reduction_method = c( "tSNE"),cell_size = 1,label_cell_groups = FALSE)+
  theme(axis.text.x = element_text(face="bold",  size=12),
        axis.text.y = element_text(face="bold",size=12),
        axis.line = element_line(colour = "darkblue", size = 2, linetype = "solid"),
        axis.title = element_text(size = rel(1.5),face="bold"),
        plot.title = element_text(color="black", size=17, face="bold.italic",hjust = 0.5,vjust =-1),
        legend.title = element_text(size=12, color = "black", face="bold"),
        legend.text = element_text(colour="black", size=12,face="bold"),
        aspect.ratio=1  #square plot
  )+labs(title="MainMeth")+ labs(fill = "log10(Exp)")


# Main
plot_cells(cds, genes=Main,reduction_method = c( "tSNE"),cell_size = 1,label_cell_groups = FALSE)+
  theme(axis.text.x = element_text(face="bold",  size=12),
        axis.text.y = element_text(face="bold",size=12),
        axis.line = element_line(colour = "darkblue", size = 2, linetype = "solid"),
        axis.title = element_text(size = rel(1.5),face="bold"),
        plot.title = element_text(color="black", size=17, face="bold.italic",hjust = 0.5,vjust =-1),
        legend.title = element_text(size=12, color = "black", face="bold"),
        legend.text = element_text(colour="black", size=12,face="bold"),
        aspect.ratio=1  #square plot
  )+labs(title="Main")+ labs(fill = "log10(Exp)")

# Macrophages
plot_cells(cds, genes=All_Macrophages_Main,reduction_method = c( "tSNE"),cell_size = 1,label_cell_groups = FALSE)+
  theme(axis.text.x = element_text(face="bold",  size=14),
        axis.text.y = element_text(face="bold",size=14),
        axis.line = element_line(colour = "darkblue", size = 2, linetype = "solid"),
        axis.title = element_text(size = rel(1.5),face="bold"),
        plot.title = element_text(color="black", size=20, face="bold.italic",hjust = 0.5,vjust =-1),
        legend.title = element_text(size=12, color = "black", face="bold"),
        legend.text = element_text(colour="black", size=12,face="bold"),
        aspect.ratio=1  #square plot
  )+labs(title="Macrophages")+ labs(fill = "log10(Exp)")

# Myeloid
plot_cells(cds, genes=All_Myeloid_Main,reduction_method = c( "tSNE"),cell_size = 1,label_cell_groups = FALSE)+
  theme(axis.text.x = element_text(face="bold",  size=14),
        axis.text.y = element_text(face="bold",size=14),
        axis.line = element_line(colour = "darkblue", size = 2, linetype = "solid"),
        axis.title = element_text(size = rel(1.5),face="bold"),
        plot.title = element_text(color="black", size=20, face="bold.italic",hjust = 0.5,vjust =-1),
        legend.title = element_text(size=12, color = "black", face="bold"),
        legend.text = element_text(colour="black", size=12,face="bold"),
        aspect.ratio=1  #square plot
  )+labs(title="Myeloid")+ labs(fill = "log10(Exp)")

# B Cells
plot_cells(cds, genes=All_BCells_Main,reduction_method = c( "tSNE"),cell_size = 1,label_cell_groups = FALSE)+
  theme(axis.text.x = element_text(face="bold",  size=14),
        axis.text.y = element_text(face="bold",size=14),
        axis.line = element_line(colour = "darkblue", size = 2, linetype = "solid"),
        axis.title = element_text(size = rel(1.5),face="bold"),
        plot.title = element_text(color="black", size=20, face="bold.italic",hjust = 0.5,vjust =-1),
        legend.title = element_text(size=12, color = "black", face="bold"),
        legend.text = element_text(colour="black", size=12,face="bold"),
        aspect.ratio=1  #square plot
  )+labs(title="B Cells")+ labs(fill = "log10(Exp)")

# Neutrophils
plot_cells(cds, genes=All_Neutrophils_Main,reduction_method = c( "tSNE"),cell_size = 1,label_cell_groups = FALSE)+
  theme(axis.text.x = element_text(face="bold",  size=12),
        axis.text.y = element_text(face="bold",size=12),
        axis.line = element_line(colour = "darkblue", size = 2, linetype = "solid"),
        axis.title = element_text(size = rel(1.5),face="bold"),
        plot.title = element_text(color="black", size=17, face="bold.italic",hjust = 0.5,vjust =-1),
        legend.title = element_text(size=12, color = "black", face="bold"),
        legend.text = element_text(colour="black", size=12,face="bold"),
        aspect.ratio=1  #square plot
  )+labs(title="Neutrophils")+ labs(fill = "log10(Exp)")

# Fibroblasts
plot_cells(cds, genes=All_Fibroblasts,reduction_method = c( "tSNE"),cell_size = 1,label_cell_groups = FALSE)+
  theme(axis.text.x = element_text(face="bold",  size=12),
        axis.text.y = element_text(face="bold",size=12),
        axis.line = element_line(colour = "darkblue", size = 2, linetype = "solid"),
        axis.title = element_text(size = rel(1.5),face="bold"),
        plot.title = element_text(color="black", size=17, face="bold.italic",hjust = 0.5,vjust =-1),
        legend.title = element_text(size=12, color = "black", face="bold"),
        legend.text = element_text(colour="black", size=12,face="bold"),
        aspect.ratio=1  #square plot
  )+labs(title="Fibroblasts")+ labs(fill = "log10(Exp)")

# T and NK Cells
plot_cells(cds, genes= All_TandNKCells,reduction_method = c( "tSNE"),cell_size = 1,label_cell_groups = FALSE)+
  theme(axis.text.x = element_text(face="bold",  size=12),
        axis.text.y = element_text(face="bold",size=12),
        axis.line = element_line(colour = "darkblue", size = 2, linetype = "solid"),
        axis.title = element_text(size = rel(1.5),face="bold"),
        plot.title = element_text(color="black", size=17, face="bold.italic",hjust = 0.5,vjust =-1),
        legend.title = element_text(size=12, color = "black", face="bold"),
        legend.text = element_text(colour="black", size=12,face="bold"),
        aspect.ratio=1  #square plot
  )+labs(title="T and NK Cells")+ labs(fill = "log10(Exp)")

# Ductal Cells 
plot_cells(cds, genes= All_Ductal_Main,reduction_method = c( "tSNE"),cell_size = 1,label_cell_groups = FALSE)+
  theme(axis.text.x = element_text(face="bold",  size=12),
        axis.text.y = element_text(face="bold",size=12),
        axis.line = element_line(colour = "darkblue", size = 2, linetype = "solid"),
        axis.title = element_text(size = rel(1.5),face="bold"),
        plot.title = element_text(color="black", size=17, face="bold.italic",hjust = 0.5,vjust =-1),
        legend.title = element_text(size=12, color = "black", face="bold"),
        legend.text = element_text(colour="black", size=12,face="bold"),
        aspect.ratio=1  #square plot
  )+labs(title="Ductal Cells ")+ labs(fill = "log10(Exp)")


# Endothelial Cells
plot_cells(cds, genes= All_Endothelial_Main,reduction_method = c( "tSNE"),cell_size = 1,label_cell_groups = FALSE)+
  theme(axis.text.x = element_text(face="bold",  size=12),
        axis.text.y = element_text(face="bold",size=12),
        axis.line = element_line(colour = "darkblue", size = 2, linetype = "solid"),
        axis.title = element_text(size = rel(1.5),face="bold"),
        plot.title = element_text(color="black", size=17, face="bold.italic",hjust = 0.5,vjust =-1),
        legend.title = element_text(size=12, color = "black", face="bold"),
        legend.text = element_text(colour="black", size=12,face="bold"),
        aspect.ratio=1  #square plot
  )+labs(title="Endothelial Cells")+ labs(fill = "log10(Exp)")

# Dendritic Cells
plot_cells(cds, genes= All_Dendritic_Main,reduction_method = c( "tSNE"),cell_size = 1,label_cell_groups = FALSE)+
  theme(axis.text.x = element_text(face="bold",  size=12),
        axis.text.y = element_text(face="bold",size=12),
        axis.line = element_line(colour = "darkblue", size = 2, linetype = "solid"),
        axis.title = element_text(size = rel(1.5),face="bold"),
        plot.title = element_text(color="black", size=17, face="bold.italic",hjust = 0.5,vjust =-1),
        legend.title = element_text(size=12, color = "black", face="bold"),
        legend.text = element_text(colour="black", size=12,face="bold"),
        aspect.ratio=1  #square plot
  )+labs(title="Dendritic Cells")+ labs(fill = "log10(Exp)")

# Acinar Cells
plot_cells(cds, genes= All_Acinar_Main,reduction_method = c( "tSNE"),cell_size = 1,label_cell_groups = FALSE)+
  theme(axis.text.x = element_text(face="bold",  size=12),
        axis.text.y = element_text(face="bold",size=12),
        axis.line = element_line(colour = "darkblue", size = 2, linetype = "solid"),
        axis.title = element_text(size = rel(1.5),face="bold"),
        plot.title = element_text(color="black", size=17, face="bold.italic",hjust = 0.5,vjust =-1),
        legend.title = element_text(size=12, color = "black", face="bold"),
        legend.text = element_text(colour="black", size=12,face="bold"),
        aspect.ratio=1  #square plot
  )+labs(title="Acinar Cells")+ labs(fill = "log10(Exp)")

# Perivascular Cells
plot_cells(cds, genes= All_Perivascular_Main,reduction_method = c( "tSNE"),cell_size = 1,label_cell_groups = FALSE)+
  theme(axis.text.x = element_text(face="bold",  size=12),
        axis.text.y = element_text(face="bold",size=12),
        axis.line = element_line(colour = "darkblue", size = 2, linetype = "solid"),
        axis.title = element_text(size = rel(1.5),face="bold"),
        plot.title = element_text(color="black", size=17, face="bold.italic",hjust = 0.5,vjust =-1),
        legend.title = element_text(size=12, color = "black", face="bold"),
        legend.text = element_text(colour="black", size=12,face="bold"),
        aspect.ratio=1  #square plot
  )+labs(title="Perivascular Cells")+ labs(fill = "log10(Exp)")

# EMT-like Cells
plot_cells(cds, genes= All_EMTlike_Main,reduction_method = c( "tSNE"),cell_size = 1,label_cell_groups = FALSE)+
  theme(axis.text.x = element_text(face="bold",  size=12),
        axis.text.y = element_text(face="bold",size=12),
        axis.line = element_line(colour = "darkblue", size = 2, linetype = "solid"),
        axis.title = element_text(size = rel(1.5),face="bold"),
        plot.title = element_text(color="black", size=17, face="bold.italic",hjust = 0.5,vjust =-1),
        legend.title = element_text(size=12, color = "black", face="bold"),
        legend.text = element_text(colour="black", size=12,face="bold"),
        aspect.ratio=1  #square plot
  )+labs(title="EMT-like Cells")+ labs(fill = "log10(Exp)")

# NSUN_Main
plot_cells(cds, genes= NSUN_Main,reduction_method = c( "tSNE"),cell_size = 1,label_cell_groups = FALSE)+
  theme(axis.text.x = element_text(face="bold",  size=12),
        axis.text.y = element_text(face="bold",size=12),
        axis.line = element_line(colour = "darkblue", size = 2, linetype = "solid"),
        axis.title = element_text(size = rel(1.5),face="bold"),
        plot.title = element_text(color="black", size=17, face="bold.italic",hjust = 0.5,vjust =-1),
        legend.title = element_text(size=12, color = "black", face="bold"),
        legend.text = element_text(colour="black", size=12,face="bold"),
        aspect.ratio=1  #square plot
  )+labs(title="NSUN")+ labs(fill = "log10(Exp)")



# png(paste0(FilesName,RVersion,"/",RVersion,"_","tSNE_iCAF.png"), width = 640, height = 360) # ?]?w???X????
# plot_cells(cds, genes=iCAF_Main,reduction_method = c( "tSNE"))
# dev.off() # ???????X????
# 
# png(paste0(FilesName,RVersion,"/",RVersion,"_","tSNE_apCAF.png"), width = 640, height = 360) # ?]?w???X????
# plot_cells(cds, genes=apCAF_Main,reduction_method = c( "tSNE"))
# dev.off() # ???????X????
# 
# png(paste0(FilesName,RVersion,"/",RVersion,"_","tSNE_myCAF.png"), width = 640, height = 360) # ?]?w???X????
# plot_cells(cds, genes=myCAF_Main,reduction_method = c( "tSNE"))
# dev.off() # ???????X????
# 
# png(paste0(FilesName,RVersion,"/",RVersion,"_","myCAF_Cd248.png"), width = 640, height = 360) # ?]?w???X????
# plot_cells(cds, genes=CD248,reduction_method = c( "tSNE"),cell_size = 1.5)
# dev.off() # ???????X????

# Group cells into clusters
cds = cluster_cells(cds, resolution=1e-5,reduction_method = c( "tSNE"),k=25)

png(paste0(FilesName,RVersion,"/",RVersion,"_","myCAF_clusterK5.png"), width = 640, height = 360) # ?]?w???X????
plot_cells(cds,reduction_method = c( "tSNE"),label_cell_groups = FALSE,cell_size = 1)
dev.off() # ???????X????

plot_cells(cds,reduction_method = c( "tSNE"),cell_size = 2,group_label_size = 5)
plot_cells(cds,reduction_method = c( "tSNE"),group_cells_by="partition",cell_size = 2,group_label_size = 5)
plot_cells(cds,reduction_method = c( "tSNE"),group_cells_by="cluster",cell_size = 2,group_label_size = 5)

######################## UMAP ###################################

cds <- reduce_dimension(cds)

png(paste0(FilesName,RVersion,"/",RVersion,"_","UMAP.png"), width = 640, height = 360) # ?]?w???X????
plot_cells(cds)
dev.off() # ???????X????

# plot_cells(cds, genes=Main)
plot_cells(cds, genes=iCAF_Main)
plot_cells(cds, genes=apCAF_Main)
plot_cells(cds, genes=myCAF_Main)
plot_cells(cds, genes=CD248,cell_size = 1)

png(paste0(FilesName,RVersion,"/",RVersion,"_","UMAP_iCAF.png"), width = 640, height = 360) # ?]?w???X????
plot_cells(cds, genes=iCAF_Main)
dev.off() # ???????X????

png(paste0(FilesName,RVersion,"/",RVersion,"_","UMAP_apCAF.png"), width = 640, height = 360) # ?]?w???X????
plot_cells(cds, genes=apCAF_Main)
dev.off() # ???????X????

png(paste0(FilesName,RVersion,"/",RVersion,"_","UMAP_myCAF.png"), width = 640, height = 360) # ?]?w???X????
plot_cells(cds, genes=myCAF_Main)
dev.off() # ???????X????

png(paste0(FilesName,RVersion,"/",RVersion,"_","UMAP_Cd248.png"), width = 640, height = 360) # ?]?w???X????
plot_cells(cds, genes=CD248,cell_size = 1)
dev.off() # ???????X????

# Group cells into clusters
cds = cluster_cells(cds, resolution=1e-5,k=25)

png(paste0(FilesName,RVersion,"/",RVersion,"_","UMAP_clusterK6.png"), width = 640, height = 360) # ?]?w???X????
plot_cells(cds,label_cell_groups = FALSE)
dev.off() # ???????X????

plot_cells(cds,cell_size = 1.5,group_label_size = 4)

plot_cells(cds, genes=All_Macrophages_Main)+labs(title="Macrophages")
plot_cells(cds, genes=All_Myeloid_Main)+labs(title="Myeloid")
plot_cells(cds, genes=All_BCells_Main)+labs(title="B Cells")
plot_cells(cds, genes=All_Neutrophils_Main)+labs(title="Neutrophils")
plot_cells(cds, genes=All_Fibroblasts)+labs(title="Fibroblasts")
plot_cells(cds, genes=All_TandNKCells)+labs(title="T and NK Cells")
plot_cells(cds, genes=All_Ductal_Main)+labs(title="Ductal Cells")
plot_cells(cds, genes=All_Endothelial_Main)+labs(title="Endothelial Cells")
plot_cells(cds, genes=All_Dendritic_Main)+labs(title="Dendritic Cells")
plot_cells(cds, genes=All_Acinar_Main)+labs(title="Acinar Cells")
plot_cells(cds, genes=All_Perivascular_Main)+labs(title="Perivascular Cells")
plot_cells(cds, genes=All_EMTlike_Main)+labs(title="EMTlike Cells")


######################## Find marker genes expressed by each cluster ###################################

plot_cells(cds,cell_size = 2,group_label_size = 5)

cds2 <- cds
cds2 = cluster_cells(cds2, reduction_method = c( "tSNE"), resolution=1e-5,k=25)

plot_cells(cds2, reduction_method = c( "tSNE"),group_cells_by="partition",cell_size = 2,group_label_size = 5,label_cell_groups = FALSE)
plot_cells(cds2, reduction_method = c( "tSNE"),group_cells_by="partition",cell_size = 2,group_label_size = 5)

plot_cells(cds2, reduction_method = c( "tSNE"),group_cells_by="cluster",cell_size = 2,group_label_size = 5,label_cell_groups = FALSE)+
theme(axis.text.x = element_text(face="bold",  size=12),
      axis.text.y = element_text(face="bold",size=12),
      axis.line = element_line(colour = "darkblue", size = 2, linetype = "solid"),
      axis.title = element_text(size = rel(1.5),face="bold"),
      plot.title = element_text(color="black", size=17, face="bold.italic",hjust = 0.5,vjust =-1),
      legend.title = element_text(size=12, color = "black", face="bold"),
      legend.text = element_text(colour="black", size=12,face="bold"),
      aspect.ratio=1  #square plot
)+labs(title="PDAC (fibroblast-enriched)")+ labs(fill = "log10(Exp)")

plot_cells(cds2, reduction_method = c( "tSNE"),group_cells_by="cluster",cell_size = 2,group_label_size = 5)+
theme(axis.text.x = element_text(face="bold",  size=12),
      axis.text.y = element_text(face="bold",size=12),
      axis.line = element_line(colour = "darkblue", size = 2, linetype = "solid"),
      axis.title = element_text(size = rel(1.5),face="bold"),
      plot.title = element_text(color="black", size=17, face="bold.italic",hjust = 0.5,vjust =-1),
      legend.title = element_text(size=12, color = "black", face="bold"),
      legend.text = element_text(colour="black", size=12,face="bold"),
      aspect.ratio=1  #square plot
)+labs(title="PDAC (fibroblast-enriched)")+ labs(fill = "log10(Exp)")

marker_test_res2 <- top_markers(cds2, group_cells_by="cluster", reduction_method = c( "tSNE"), 
                               reference_cells=1000, cores=8)

top_specific_markers2 <- marker_test_res2 %>%
  filter(fraction_expressing >= 0.10) %>%
  group_by(cell_group) %>%
  top_n(1, pseudo_R2)

top_specific_marker_ids2 <- unique(top_specific_markers2 %>% pull(gene_id))

plot_genes_by_group(cds2, reduction_method = c( "tSNE"),
                    top_specific_marker_ids2,
                    group_cells_by="cluster",
                    ordering_type="maximal_on_diag",
                    max.size=3)

top_specific_markers2 <- marker_test_res2 %>%
  filter(fraction_expressing >= 0.10) %>%
  group_by(cell_group) %>%
  top_n(10, pseudo_R2)

top_specific_marker_ids2 <- unique(top_specific_markers2 %>% pull(gene_id))

plot_genes_by_group(cds2, reduction_method = c( "tSNE"),
                    top_specific_marker_ids2,
                    group_cells_by="cluster",
                    ordering_type="cluster_row_col",
                    max.size=3)

plot_cells(cds2, reduction_method = c( "tSNE"),group_cells_by="cluster",cell_size = 2,group_label_size = 5,label_cell_groups = FALSE)
plot_cells(cds2, reduction_method = c( "tSNE"),group_cells_by="cluster",cell_size = 2,group_label_size = 5)

plot_cells(cds2, color_cells_by="cluster",cell_size = 2,group_label_size = 5)
colData(cds2)$assigned_cell_type <- as.character(clusters(cds2, reduction_method = c( "tSNE"))[colnames(cds2)])
colData(cds2)$assigned_cell_type <- dplyr::recode(colData(cds2, reduction_method = c( "tSNE"))$assigned_cell_type,
                                                        "1"="01.Macrophages",
                                                        "11"="01.Macrophages", 
                                                  
                                                        "2"="02.Myeloid Cells",
                                                        "12"="02.Myeloid Cells",                                                  
                                                        "16"="02.Myeloid Cells",
                                                  
                                                        "8"="03.B cells",
                                                  
                                                        "3"="04.Neutrophils",
                                                  
                                                        "5"="05.Fibroblasts",
  
                                                        "6"="06.T and NK Cells",
                                                  
                                                        "4"="07.Ductal cell 1",
                                                  
                                                        "10"="08.Endothelial Cells",
                                                  
                                                        "9"="09.Dendritic Cells",
                                                        "12"="09.Dendritic Cells",
                                                  
#                                                        "15"="10.Acinar Cells",

                                                        "13"="11.Perivascular Cells",
                                                                                                         
                                                        "7"="12.EMT-like Cells"
)


cds2@colData@listData[["PDACAll"]] <- cds2@colData@listData[["assigned_cell_type"]]

plot_cells(cds2, group_cells_by="cluster", color_cells_by="PDACAll",
           reduction_method = c( "tSNE"),cell_size = 2,group_label_size = 5,label_cell_groups = FALSE)+
  theme(axis.text.x = element_text(face="bold",  size=10),
        axis.text.y = element_text(face="bold",size=10),
        axis.line = element_line(colour = "darkblue", size = 2, linetype = "solid"),
        axis.title = element_text(size = rel(1.5),face="bold"),
        plot.title = element_text(color="black", size=17, face="bold.italic",hjust = 0.5,vjust =-1),
        legend.title = element_text(size=12, color = "black", face="bold"),
        legend.text = element_text(colour="black", size=12,face="bold"),
        aspect.ratio=1  #square plot
  )+labs(title="PDAC All")+ labs(fill = "log10(Exp)")

plot_cells(cds2, group_cells_by="cluster", color_cells_by="PDACAll",
           reduction_method = c( "tSNE"),cell_size = 2,group_label_size = 5)
  

######################## Isolate_subset ######################## 
cds_subsetV2Ori <- choose_cells(cds2,reduction_method = c( "tSNE"),
                            clear_cds = FALSE,
                            return_list = FALSE)
cds_subsetV2 <- cds_subsetV2Ori

cds_subsetV2 = cluster_cells(cds_subsetV2, resolution=1e-2,K=2)
# plot_cells(cds_subsetV2, color_cells_by="cluster",reduction_method = c( "tSNE"),cell_size = 2,group_label_size = 5)

plot_cells(cds_subsetV2, color_cells_by="cluster",reduction_method = c( "tSNE"),cell_size = 2,label_cell_groups = FALSE)
plot_cells(cds_subsetV2, color_cells_by="cluster",reduction_method = c( "tSNE"),cell_size = 2,group_label_size = 5)

plot_cells(cds_subsetV2, color_cells_by="cluster",cell_size = 2,group_label_size = 5,label_cell_groups = FALSE)
plot_cells(cds_subsetV2, color_cells_by="cluster",cell_size = 2,group_label_size = 5)
colData(cds_subsetV2)$assigned_cell_type <- as.character(clusters(cds_subsetV2)[colnames(cds_subsetV2)])
colData(cds_subsetV2)$assigned_cell_type <- dplyr::recode(colData(cds_subsetV2)$assigned_cell_type,
                                                        "1"="10.Acinar Cells"
)

cds_subsetV2@colData@listData[["Acinar"]] <- cds_subsetV2@colData@listData[["assigned_cell_type"]]

plot_cells(cds_subsetV2, group_cells_by="cluster", color_cells_by="assigned_cell_type",label_cell_groups = FALSE,
           reduction_method = c( "tSNE"),cell_size = 2,group_label_size = 5)

######################## Isolate_subset ######################## 

#######Now we can transfer the annotations from the cds_subset object back to the full dataset. We'll also filter out low-quality cells at this stage ###
colData(cds2)[colnames(cds_subsetV2),]$assigned_cell_type <- colData(cds_subsetV2)$assigned_cell_type
plot_cells(cds2, group_cells_by="cluster", color_cells_by="assigned_cell_type",
           reduction_method = c( "tSNE"),cell_size = 2,group_label_size = 5,label_cell_groups = FALSE)+
  theme(axis.text.x = element_text(face="bold",  size=10),
        axis.text.y = element_text(face="bold",size=10),
        axis.line = element_line(colour = "darkblue", size = 2, linetype = "solid"),
        axis.title = element_text(size = rel(1.5),face="bold"),
        plot.title = element_text(color="black", size=17, face="bold.italic",hjust = 0.5,vjust =-1),
        legend.title = element_text(size=12, color = "black", face="bold"),
        legend.text = element_text(colour="black", size=12,face="bold"),
        aspect.ratio=1  #square plot
  )+labs(title="PDAC All")+ labs(fill = "log10(Exp)")


plot_cells(cds2, group_cells_by="cluster", color_cells_by="assigned_cell_type",
           reduction_method = c( "tSNE"),cell_size = 2,group_label_size = 5)


### Main Violin
cds2_Violin <- cds2
cds2_Violin <- cds2_Violin[rowData(cds2_Violin)$gene_short_name %in% Main,]

plot_genes_violin(cds2_Violin, group_cells_by="assigned_cell_type",  log_scale = FALSE, ncol=3) +
  theme(axis.text.x=element_text(angle=65, hjust=1))+
  theme(axis.text.x = element_text(face="bold",  size=10),
        axis.text.y = element_text(face="bold",size=12),
        axis.line = element_line(colour = "darkblue", size = 2, linetype = "solid"),
        axis.title = element_text(size = rel(1.5),face="bold"),
        plot.title = element_text(color="black", size=17, face="bold.italic",hjust = 0.5,vjust =-1),
        legend.title = element_text(size=12, color = "black", face="bold"),
        legend.text = element_text(colour="black", size=12,face="bold"),
        aspect.ratio=1  #square plot
  )+labs(title="Main")+ labs(fill = "log10(Exp)")


### MainCD Violin
cds2_ViolinCD <- cds2
cds2_ViolinCD <- cds2_ViolinCD[rowData(cds2_ViolinCD)$gene_short_name %in% MainCD,]

plot_genes_violin(cds2_ViolinCD, group_cells_by="assigned_cell_type",  log_scale = FALSE, ncol=3) +
  theme(axis.text.x=element_text(angle=65, hjust=1))+
  theme(axis.text.x = element_text(face="bold",  size=10),
        axis.text.y = element_text(face="bold",size=12),
        axis.line = element_line(colour = "darkblue", size = 2, linetype = "solid"),
        axis.title = element_text(size = rel(1.5),face="bold"),
        plot.title = element_text(color="black", size=17, face="bold.italic",hjust = 0.5,vjust =-1),
        legend.title = element_text(size=12, color = "black", face="bold"),
        legend.text = element_text(colour="black", size=12,face="bold"),
        aspect.ratio=1  #square plot
  )+labs(title="MainCD")+ labs(fill = "log10(Exp)")

### MainMeth Violin
cds2_ViolinMeth <- cds2
cds2_ViolinMeth <- cds2_ViolinMeth[rowData(cds2_ViolinMeth)$gene_short_name %in% MainMeth,]

plot_genes_violin(cds2_ViolinMeth, group_cells_by="assigned_cell_type",  log_scale = FALSE, ncol=3) +
  theme(axis.text.x=element_text(angle=65, hjust=1))+
  theme(axis.text.x = element_text(face="bold",  size=10),
        axis.text.y = element_text(face="bold",size=12),
        axis.line = element_line(colour = "darkblue", size = 2, linetype = "solid"),
        axis.title = element_text(size = rel(1.5),face="bold"),
        plot.title = element_text(color="black", size=17, face="bold.italic",hjust = 0.5,vjust =-1),
        legend.title = element_text(size=12, color = "black", face="bold"),
        legend.text = element_text(colour="black", size=12,face="bold"),
        aspect.ratio=1  #square plot
  )+labs(title="MainMeth")+ labs(fill = "log10(Exp)")


# NSUN_Main
cds2_Violin_NSUN <- cds2
cds2_Violin_NSUN <- cds2_Violin_NSUN[rowData(cds2_Violin_NSUN)$gene_short_name %in% NSUN_Main,]

plot_genes_violin(cds2_Violin_NSUN, group_cells_by="assigned_cell_type",  log_scale = FALSE, ncol=3) +
  theme(axis.text.x=element_text(angle=65, hjust=1))+
  theme(axis.text.x = element_text(face="bold",  size=10),
        axis.text.y = element_text(face="bold",size=12),
        axis.line = element_line(colour = "darkblue", size = 2, linetype = "solid"),
        axis.title = element_text(size = rel(1.5),face="bold"),
        plot.title = element_text(color="black", size=17, face="bold.italic",hjust = 0.5,vjust =-1),
        legend.title = element_text(size=12, color = "black", face="bold"),
        legend.text = element_text(colour="black", size=12,face="bold"),
        aspect.ratio=1  #square plot
  )+labs(title="NSUN")+ labs(fill = "log10(Exp)")

