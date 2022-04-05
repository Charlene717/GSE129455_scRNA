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
  library(AnnotationDbi)
  library(org.Mm.eg.db)
  library(Hmisc)
  
  
#### Load files ####
  setwd(getwd()) ## Set current working directory
  PathName <- getwd() ## Set output directroy
  FilesName = paste0(PathName,"/")
  
  RVersion = "20210304V1"
  dir.create(paste0(FilesName,RVersion))

  data <- read.csv(paste0(FilesName,"GSE129455_Fibroblast-enriched.csv"),  
                   header=T,         
                   sep=",")          
  
  row.names(data) <- data[,1]
  
  expression_matrix2 <- data[1:length(data[,1]), 2:length(data[1,])]
  expression_matrix3 <- as(as.matrix(expression_matrix2), "dgCMatrix")

  
#### Build monocle obj ####   
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
  
#### Gene name conversion #### 
  GeneNameEnsembl <- as.matrix(cds@rowRanges@elementMetadata@listData$gene_short_name)
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

#### Main gene #### 
  iCAF_Main = c("ENSMUSG00000025784","ENSMUSG00000022371","ENSMUSG00000003665","ENSMUSG00000025746")
  # Clec3b,Col14a1,Has1,Il6
  iCAF_Main = c("Clec3b","Col14a1","Has1","Il6")
  
  apCAF_Main = c("ENSMUSG00000073421","ENSMUSG00000024610","ENSMUSG00000040026","ENSMUSG00000017002")
  # H2???Ab1, Cd74, Saa3, Slpi
  apCAF_Main = c("ENSMUSG00000073421","Cd74","Saa3","Slpi")
  apCAF_Main = c("H2-Ab1","Cd74","Saa3","Slpi")
  
  myCAF_Main = c("ENSMUSG00000032085","ENSMUSG00000032011","ENSMUSG00000032332","ENSMUSG00000023885")
  #Tagln Thy1 Col12a1 Thbs2
  myCAF_Main = c("Tagln","Thy1","Col12a1","Thbs2")
  
  # Lipofibroblasts_Main = c("Car3","Lrg1","Apoe","Tmem176a")
  Lipofibroblasts_Main = c("Car3","Lrg1","Apoe","Tmem176a","Mgst1","Tmem176b")
  
  # Perivascular_cell_Main = c("Rgs5","Acta2","Sparc","Vim","Igfbp7")
  Perivascular_cell_Main = c("Rgs5","Ndufa4l2","Mfge8","Sparcl1","Des","Esam")
  
  # Neutrophils_Main = c("S100a8","S100a9","G0s2")
  Neutrophils_Main = c("S100a9","S100a8","Tyrobp","Fcer1g","Srgn","Il1r2")
  
  # Macrophages_Main = c("Apoe","Saa3","C1qc","Fcerlg","Fcgr3","Laptm5")
  Macrophages_Main = c("Fcer1g","Fcgr3","Laptm5","Tyrobp","Ctss","Cfp")
  
  Ductal1_Main = c("Npy","Tgfb1","Cd9","Plaur","Cadm1","Nt5e")
  Ductal2_Main = c("Perp","Ctsd","Cdkn1a","Phlda3","Cd81","Cd9")
  Ductal3_Main = c("Cdkn2a","Pkm","Msn","Stmn2","Klra4","Upp1")
  Ductal4_Main = c("Cystm1","Clu","Wfdc2","Ctse","2200002D01Rik","Spint2")
  
  NSUN_Main = c("Nsun1","Nop2","Nsun2","Nsun3","Nsun4","Nsun5","Nsun6","Nsun7")
  NSUN2 = "Nsun2"


#### Pre-process the data ####
  cds = preprocess_cds(cds, num_dim = 100)
  
  png(paste0(FilesName,RVersion,"/",RVersion,"_","PCA.png"), width = 640, height = 360) # ?]?w???X????
  plot_pc_variance_explained(cds)
  
  dev.off() 



#### tSNE ####
  cds <- reduce_dimension(cds,reduction_method = c( "tSNE") , perplexity = 60)
  
  png(paste0(FilesName,RVersion,"/",RVersion,"_","tSNE.png"), width = 640, height = 360) # ?]?w???X????
  plot_cells(cds,reduction_method = c( "tSNE"))
  dev.off() 
  
  # plot_cells(cds, genes=Main)
  plot_cells(cds, genes=iCAF_Main,reduction_method = c( "tSNE"))
  plot_cells(cds, genes=apCAF_Main,reduction_method = c( "tSNE"))
  plot_cells(cds, genes=myCAF_Main,reduction_method = c( "tSNE"))
  plot_cells(cds, genes=NSUN2,reduction_method = c( "tSNE"),cell_size = 1)
  
  mypalette <- colorRampPalette(c("white" , "red"))
  plot_cells(cds, genes=NSUN2,reduction_method = c( "tSNE"),cell_size = 1,label_cell_groups = FALSE)+ theme(legend.position="top")
  

  plot_cells(cds, genes=NSUN2,reduction_method = c( "tSNE"),cell_size = 1,label_cell_groups = FALSE)+
   theme(axis.text.x = element_text(face="bold", color="#993333", 
                                     size=14, angle=45),
          axis.text.y = element_text(face="bold", color="#993333", 
                                     size=14, angle=45))

  # NSUN2
  plot_cells(cds, genes=NSUN2,reduction_method = c( "tSNE"),cell_size = 1,label_cell_groups = FALSE)+
    theme(axis.text.x = element_text(face="bold",  size=14),
          axis.text.y = element_text(face="bold",size=14),
          axis.line = element_line(colour = "darkblue", size = 2, linetype = "solid"),
          axis.title = element_text(size = rel(1.5),face="bold"),
          plot.title = element_text(color="black", size=20, face="bold.italic",hjust = 0.5,vjust = - 6),
  #        panel.background = element_rect(fill = "lightblue"),
          legend.title = element_text(size=12, color = "black", face="bold",inherit.blank=TRUE),
          legend.text = element_text(colour="black", size=12,face="bold"),
          aspect.ratio=1  #square plot
          )+ labs(fill = "log10(Exp)")  # +labs(title="Nsun2")
  
  #+ theme(legend.position="top")
  
  # iCAF
  plot_cells(cds, genes=iCAF_Main,reduction_method = c( "tSNE"),cell_size = 1,label_cell_groups = FALSE)+
    theme(axis.text.x = element_text(face="bold",  size=14),
          axis.text.y = element_text(face="bold",size=14),
          axis.line = element_line(colour = "darkblue", size = 2, linetype = "solid"),
          axis.title = element_text(size = rel(1.5),face="bold"),
          plot.title = element_text(color="black", size=20, face="bold.italic",hjust = 0.5,vjust =-1),
          legend.title = element_text(size=12, color = "black", face="bold"),
          legend.text = element_text(colour="black", size=12,face="bold"),
          aspect.ratio=1  #square plot
    )+labs(title="iCAF")+ labs(fill = "log10(Exp)")
  
  # myCAF
  plot_cells(cds, genes=myCAF_Main,reduction_method = c( "tSNE"),cell_size = 1,label_cell_groups = FALSE)+
    theme(axis.text.x = element_text(face="bold",  size=14),
          axis.text.y = element_text(face="bold",size=14),
          axis.line = element_line(colour = "darkblue", size = 2, linetype = "solid"),
          axis.title = element_text(size = rel(1.5),face="bold"),
          plot.title = element_text(color="black", size=20, face="bold.italic",hjust = 0.5,vjust =-1),
          legend.title = element_text(size=12, color = "black", face="bold"),
          legend.text = element_text(colour="black", size=12,face="bold"),
          aspect.ratio=1  #square plot
    )+labs(title="myCAF")+ labs(fill = "log10(Exp)")
  
  # apCAF
  plot_cells(cds, genes=apCAF_Main,reduction_method = c( "tSNE"),cell_size = 1,label_cell_groups = FALSE)+
    theme(axis.text.x = element_text(face="bold",  size=14),
          axis.text.y = element_text(face="bold",size=14),
          axis.line = element_line(colour = "darkblue", size = 2, linetype = "solid"),
          axis.title = element_text(size = rel(1.5),face="bold"),
          plot.title = element_text(color="black", size=20, face="bold.italic",hjust = 0.5,vjust =-1),
          legend.title = element_text(size=12, color = "black", face="bold"),
          legend.text = element_text(colour="black", size=12,face="bold"),
          aspect.ratio=1  #square plot
    )+labs(title="apCAF")+ labs(fill = "log10(Exp)")
  
  # Perivascular_cell
  plot_cells(cds, genes=Perivascular_cell_Main,reduction_method = c( "tSNE"),cell_size = 1,label_cell_groups = FALSE)+
    theme(axis.text.x = element_text(face="bold",  size=12),
          axis.text.y = element_text(face="bold",size=12),
          axis.line = element_line(colour = "darkblue", size = 2, linetype = "solid"),
          axis.title = element_text(size = rel(1.5),face="bold"),
          plot.title = element_text(color="black", size=17, face="bold.italic",hjust = 0.5,vjust =-1),
          legend.title = element_text(size=12, color = "black", face="bold"),
          legend.text = element_text(colour="black", size=12,face="bold"),
          aspect.ratio=1  #square plot
    )+labs(title="Perivascular cell")+ labs(fill = "log10(Exp)")
  
  # Neutrophils_Main
  plot_cells(cds, genes=Neutrophils_Main,reduction_method = c( "tSNE"),cell_size = 1,label_cell_groups = FALSE)+
    theme(axis.text.x = element_text(face="bold",  size=12),
          axis.text.y = element_text(face="bold",size=12),
          axis.line = element_line(colour = "darkblue", size = 2, linetype = "solid"),
          axis.title = element_text(size = rel(1.5),face="bold"),
          plot.title = element_text(color="black", size=17, face="bold.italic",hjust = 0.5,vjust =-1),
          legend.title = element_text(size=12, color = "black", face="bold"),
          legend.text = element_text(colour="black", size=12,face="bold"),
          aspect.ratio=1  #square plot
    )+labs(title="Neutrophils")+ labs(fill = "log10(Exp)")
  
  # Macrophages_Main
  plot_cells(cds, genes= Macrophages_Main,reduction_method = c( "tSNE"),cell_size = 1,label_cell_groups = FALSE)+
    theme(axis.text.x = element_text(face="bold",  size=12),
          axis.text.y = element_text(face="bold",size=12),
          axis.line = element_line(colour = "darkblue", size = 2, linetype = "solid"),
          axis.title = element_text(size = rel(1.5),face="bold"),
          plot.title = element_text(color="black", size=17, face="bold.italic",hjust = 0.5,vjust =-1),
          legend.title = element_text(size=12, color = "black", face="bold"),
          legend.text = element_text(colour="black", size=12,face="bold"),
          aspect.ratio=1  #square plot
    )+labs(title="Macrophages")+ labs(fill = "log10(Exp)")
  
  # Lipofibroblasts_Main 
  plot_cells(cds, genes= Lipofibroblasts_Main,reduction_method = c( "tSNE"),cell_size = 1,label_cell_groups = FALSE)+
    theme(axis.text.x = element_text(face="bold",  size=12),
          axis.text.y = element_text(face="bold",size=12),
          axis.line = element_line(colour = "darkblue", size = 2, linetype = "solid"),
          axis.title = element_text(size = rel(1.5),face="bold"),
          plot.title = element_text(color="black", size=17, face="bold.italic",hjust = 0.5,vjust =-1),
          legend.title = element_text(size=12, color = "black", face="bold"),
          legend.text = element_text(colour="black", size=12,face="bold"),
          aspect.ratio=1  #square plot
    )+labs(title="Lipofibroblasts")+ labs(fill = "log10(Exp)")
  
  
  # Ductal1_Main
  plot_cells(cds, genes= Ductal1_Main,reduction_method = c( "tSNE"),cell_size = 1,label_cell_groups = FALSE)+
    theme(axis.text.x = element_text(face="bold",  size=12),
          axis.text.y = element_text(face="bold",size=12),
          axis.line = element_line(colour = "darkblue", size = 2, linetype = "solid"),
          axis.title = element_text(size = rel(1.5),face="bold"),
          plot.title = element_text(color="black", size=17, face="bold.italic",hjust = 0.5,vjust =-1),
          legend.title = element_text(size=12, color = "black", face="bold"),
          legend.text = element_text(colour="black", size=12,face="bold"),
          aspect.ratio=1  #square plot
    )+labs(title="Ductal1")+ labs(fill = "log10(Exp)")
  
  # Ductal2_Main
  plot_cells(cds, genes= Ductal2_Main,reduction_method = c( "tSNE"),cell_size = 1,label_cell_groups = FALSE)+
    theme(axis.text.x = element_text(face="bold",  size=12),
          axis.text.y = element_text(face="bold",size=12),
          axis.line = element_line(colour = "darkblue", size = 2, linetype = "solid"),
          axis.title = element_text(size = rel(1.5),face="bold"),
          plot.title = element_text(color="black", size=17, face="bold.italic",hjust = 0.5,vjust =-1),
          legend.title = element_text(size=12, color = "black", face="bold"),
          legend.text = element_text(colour="black", size=12,face="bold"),
          aspect.ratio=1  #square plot
    )+labs(title="Ductal2")+ labs(fill = "log10(Exp)")
  
  # Ductal3_Main
  plot_cells(cds, genes= Ductal3_Main,reduction_method = c( "tSNE"),cell_size = 1,label_cell_groups = FALSE)+
    theme(axis.text.x = element_text(face="bold",  size=12),
          axis.text.y = element_text(face="bold",size=12),
          axis.line = element_line(colour = "darkblue", size = 2, linetype = "solid"),
          axis.title = element_text(size = rel(1.5),face="bold"),
          plot.title = element_text(color="black", size=17, face="bold.italic",hjust = 0.5,vjust =-1),
          legend.title = element_text(size=12, color = "black", face="bold"),
          legend.text = element_text(colour="black", size=12,face="bold"),
          aspect.ratio=1  #square plot
    )+labs(title="Ductal3")+ labs(fill = "log10(Exp)")
  
  # Ductal4_Main
  plot_cells(cds, genes= Ductal4_Main,reduction_method = c( "tSNE"),cell_size = 1,label_cell_groups = FALSE)+
    theme(axis.text.x = element_text(face="bold",  size=12),
          axis.text.y = element_text(face="bold",size=12),
          axis.line = element_line(colour = "darkblue", size = 2, linetype = "solid"),
          axis.title = element_text(size = rel(1.5),face="bold"),
          plot.title = element_text(color="black", size=17, face="bold.italic",hjust = 0.5,vjust =-1),
          legend.title = element_text(size=12, color = "black", face="bold"),
          legend.text = element_text(colour="black", size=12,face="bold"),
          aspect.ratio=1  #square plot
    )+labs(title="Ductal4")+ labs(fill = "log10(Exp)")
  
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
  


  cds = cluster_cells(cds, resolution=1e-5,reduction_method = c( "tSNE"),k=5)
  
  png(paste0(FilesName,RVersion,"/",RVersion,"_","myCAF_clusterK5.png"), width = 640, height = 360) # ?]?w???X????
  plot_cells(cds,reduction_method = c( "tSNE"),label_cell_groups = FALSE,cell_size = 1)
  dev.off() # ???????X????
  
  plot_cells(cds,reduction_method = c( "tSNE"),cell_size = 2,group_label_size = 5)
  plot_cells(cds,reduction_method = c( "tSNE"),group_cells_by="partition",cell_size = 2,group_label_size = 5)
  plot_cells(cds,reduction_method = c( "tSNE"),group_cells_by="cluster",cell_size = 2,group_label_size = 5)

##### UMAP #####
  cds <- reduce_dimension(cds)
  
  png(paste0(FilesName,RVersion,"/",RVersion,"_","UMAP.png"), width = 640, height = 360) # ?]?w???X????
  plot_cells(cds)
  dev.off() 
  
  plot_cells(cds, genes=iCAF_Main)
  plot_cells(cds, genes=apCAF_Main)
  plot_cells(cds, genes=myCAF_Main)
  plot_cells(cds, genes=NSUN2,cell_size = 1)
  
  png(paste0(FilesName,RVersion,"/",RVersion,"_","UMAP_iCAF.png"), width = 640, height = 360) # ?]?w???X????
  plot_cells(cds, genes=iCAF_Main)
  dev.off()
  
  png(paste0(FilesName,RVersion,"/",RVersion,"_","UMAP_apCAF.png"), width = 640, height = 360) # ?]?w???X????
  plot_cells(cds, genes=apCAF_Main)
  dev.off()
  
  png(paste0(FilesName,RVersion,"/",RVersion,"_","UMAP_myCAF.png"), width = 640, height = 360) # ?]?w???X????
  plot_cells(cds, genes=myCAF_Main)
  dev.off()
  
  png(paste0(FilesName,RVersion,"/",RVersion,"_","UMAP_NSUN2.png"), width = 640, height = 360) # ?]?w???X????
  plot_cells(cds, genes=NSUN2,cell_size = 1)
  dev.off() 
  
  # Group cells into clusters
  cds = cluster_cells(cds, resolution=1e-5,k=6)
  
  png(paste0(FilesName,RVersion,"/",RVersion,"_","UMAP_clusterK6.png"), width = 640, height = 360) # ?]?w???X????
  plot_cells(cds,label_cell_groups = FALSE)
  dev.off() 


#### Isolate_subset ####
  cds_subset2 <- choose_cells(cds,reduction_method = c( "tSNE"),
                              clear_cds = FALSE,
                              return_list = FALSE)
  cds_subset <- cds_subset2
  
  cds_subset = cluster_cells(cds_subset, resolution=1e-2)
  # plot_cells(cds_subset, color_cells_by="cluster",reduction_method = c( "tSNE"),cell_size = 2,group_label_size = 5)
  
  plot_cells(cds_subset, color_cells_by="cluster",reduction_method = c( "tSNE"),cell_size = 2,label_cell_groups = FALSE)
  plot_cells(cds_subset, color_cells_by="cluster",reduction_method = c( "tSNE"),cell_size = 2,group_label_size = 5)

  plot_cells(cds_subset, color_cells_by="cluster",cell_size = 2,group_label_size = 5)
  colData(cds_subset)$assigned_cell_type <- as.character(clusters(cds_subset)[colnames(cds_subset)])
  colData(cds_subset)$assigned_cell_type <- dplyr::recode(colData(cds_subset)$assigned_cell_type,
                                                          "4"="apCAF",
                                                          "5"="apCAF",
                                                          "26"="apCAF",
                                                          "29"="apCAF",
                                                          "6"="apCAF",
                                                          "10"="apCAF",
                                                          
                                                          "22"="myCAF",
                                                          "18"="myCAF",
                                                          "34"="myCAF",
                                                          "1"="myCAF",
                                                          "25"="myCAF",
                                                          "27"="myCAF",
                                                          "35"="myCAF",
                                                          "2"="myCAF",
                                                          "11"="myCAF",
                                                          "28"="myCAF",
                                                          "14"="myCAF",
                                                          "19"="myCAF",
                                                          "30"="myCAF",
                                                          "31"="myCAF",
                                                          
                                                          "24"="iCAF",
                                                          "15"="iCAF",
                                                          "16"="iCAF",
                                                          "3"="iCAF",
                                                          "21"="iCAF",
                                                          "7"="iCAF",
                                                          "23"="iCAF",
                                                          "13"="iCAF",
                                                          "9"="iCAF",
                                                          "8"="iCAF",
                                                          "12"="iCAF",
                                                          "20"="iCAF",
                                                          "32"="iCAF",
                                                          "17"="iCAF",
                                                          "33"="iCAF",
                                                          "36"="iCAF"
  )
  
  
  cds_subset@colData@listData[["CAF"]] <- cds_subset@colData@listData[["assigned_cell_type"]]
  
  #
  plot_cells(cds_subset, group_cells_by="cluster", color_cells_by="assigned_cell_type",
             reduction_method = c( "tSNE"),cell_size = 2,group_label_size = 5,label_cell_groups = FALSE)+
    theme(axis.text.x = element_text(face="bold",  size=10),
          axis.text.y = element_text(face="bold",size=10),
          axis.line = element_line(colour = "darkblue", size = 2, linetype = "solid"),
          axis.title = element_text(size = rel(1.5),face="bold"),
          plot.title = element_text(color="black", size=17, face="bold.italic",hjust = 0.5,vjust =-1),
          legend.title = element_text(size=12, color = "black", face="bold"),
          legend.text = element_text(colour="black", size=12,face="bold"),
          aspect.ratio=1  #square plot
    )+labs(title="CAF")+ labs(fill = "log10(Exp)")
  
  #
  plot_cells(cds_subset, group_cells_by="cluster", color_cells_by="assigned_cell_type",
             reduction_method = c( "tSNE"),cell_size = 2,group_label_size = 5)
  
  ###
  cds_subset_Violin <- cds_subset
  
  ciliated_genes <- c("Nsun2")
  ### Main
  cds_subset_Violin <- cds_subset_Violin[rowData(cds_subset_Violin)$gene_short_name %in% ciliated_genes,]
  
  plot_genes_violin(cds_subset_Violin, group_cells_by="CAF",  log_scale = FALSE, ncol=3) +
    theme(axis.text.x=element_text(angle=45, hjust=1))+
    theme(axis.text.x = element_text(face="bold",  size=12),
          axis.text.y = element_text(face="bold",size=12),
          axis.line = element_line(colour = "darkblue", size = 2, linetype = "solid"),
          axis.title = element_text(size = rel(1.5),face="bold"),
          plot.title = element_text(color="black", size=17, face="bold.italic",hjust = 0.5,vjust =-1),
          legend.title = element_text(size=12, color = "black", face="bold"),
          legend.text = element_text(colour="black", size=12,face="bold"),
          aspect.ratio=1  #square plot
    )+labs(title="Main")+ labs(fill = "log10(Exp)")
  
  ### iCAF
  cds_subset_Violin_iCAF <- cds_subset
  
  cds_subset_Violin_iCAF <- cds_subset_Violin_iCAF[rowData(cds_subset_Violin_iCAF)$gene_short_name %in% iCAF_Main,]
  
  plot_genes_violin(cds_subset_Violin_iCAF, group_cells_by="CAF",  log_scale = FALSE, ncol=2) +
    theme(axis.text.x=element_text(angle=45, hjust=1))+
    theme(axis.text.x = element_text(face="bold",  size=10),
          axis.text.y = element_text(face="bold",size=10),
          axis.line = element_line(colour = "darkblue", size = 2, linetype = "solid"),
          axis.title = element_text(size = rel(1.5),face="bold"),
          plot.title = element_text(color="black", size=17, face="bold.italic",hjust = 0.5,vjust =-1),
          legend.title = element_text(size=12, color = "black", face="bold"),
          legend.text = element_text(colour="black", size=12,face="bold"),
          aspect.ratio=1  #square plot
    )+labs(title="iCAF markers")+ labs(fill = "log10(Exp)")
  
  ### myCAF
  cds_subset_Violin_myCAF <- cds_subset
  
  cds_subset_Violin_myCAF <- cds_subset_Violin_myCAF[rowData(cds_subset_Violin_myCAF)$gene_short_name %in% myCAF_Main,]
  
  plot_genes_violin(cds_subset_Violin_myCAF, group_cells_by="CAF",  log_scale = FALSE, ncol=2) +
    theme(axis.text.x=element_text(angle=45, hjust=1))+
    theme(axis.text.x = element_text(face="bold",  size=10),
          axis.text.y = element_text(face="bold",size=10),
          axis.line = element_line(colour = "darkblue", size = 2, linetype = "solid"),
          axis.title = element_text(size = rel(1.5),face="bold"),
          plot.title = element_text(color="black", size=17, face="bold.italic",hjust = 0.5,vjust =-1),
          legend.title = element_text(size=12, color = "black", face="bold"),
          legend.text = element_text(colour="black", size=12,face="bold"),
          aspect.ratio=1  #square plot
    )+labs(title="myCAF markers")+ labs(fill = "log10(Exp)")
  
  
  ### apCAF
  cds_subset_Violin_apCAF <- cds_subset
  
  cds_subset_Violin_apCAF <- cds_subset_Violin_apCAF[rowData(cds_subset_Violin_apCAF)$gene_short_name %in% apCAF_Main,]
  
  plot_genes_violin(cds_subset_Violin_apCAF, group_cells_by="CAF",  log_scale = FALSE, ncol=2) +
    theme(axis.text.x=element_text(angle=45, hjust=1))+
    theme(axis.text.x = element_text(face="bold",  size=10),
          axis.text.y = element_text(face="bold",size=10),
          axis.line = element_line(colour = "darkblue", size = 2, linetype = "solid"),
          axis.title = element_text(size = rel(1.5),face="bold"),
          plot.title = element_text(color="black", size=17, face="bold.italic",hjust = 0.5,vjust =-1),
          legend.title = element_text(size=12, color = "black", face="bold"),
          legend.text = element_text(colour="black", size=12,face="bold"),
          aspect.ratio=1  #square plot
    )+labs(title="apCAF markers")+ labs(fill = "log10(Exp)")




#### Find marker genes expressed by each cluster ####
  plot_cells(cds,cell_size = 2,group_label_size = 5)
  
  cds2 <- cds
  cds2 = cluster_cells(cds2, reduction_method = c( "tSNE"), resolution=1e-5,k=12)
  
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
  colData(cds2)$assigned_cell_type <- as.character(clusters(cds2)[colnames(cds2)])
  colData(cds2)$assigned_cell_type <- dplyr::recode(colData(cds2)$assigned_cell_type,
                                                          "7"="07.Perivascular cell",
                                                          
                                                          "1"="09.iCAF",                                                       
                                                          "5"="11.myCAF",
                                                          "6"="10.apCAF",
                                                          "19"="02.Ductal cell 2",
    
                                                          "11"="01.Ductal cell 1",
                                                          "12"="01.Ductal cell 1",
                                                          "13"="01.Ductal cell 1",
                                                    
                                                          "3"="02.Ductal cell 2",
                                                          "4"="02.Ductal cell 2",
                                                          "16"="02.Ductal cell 2",
                                                          "18"="02.Ductal cell 2",
  
                                                    
                                                          "15"="05.Macrophages",
                                                          "17"="06.Neutrophils",
                                                                                                           
                                                          "8"="03.Ductal cell 3",
                                                          "9"="03.Ductal cell 3",
                                                          "14"="03.Ductal cell 3",   
                                                    
                                                          "2"="04.Ductal cell 4",
                                                          "10"="04.Ductal cell 4"
  
  )
  
  
  cds2@colData@listData[["PDACFE"]] <- cds2@colData@listData[["assigned_cell_type"]]
  
  plot_cells(cds2, group_cells_by="cluster", color_cells_by="PDACFE",
             reduction_method = c( "tSNE"),cell_size = 2,group_label_size = 5,label_cell_groups = FALSE)+
    theme(axis.text.x = element_text(face="bold",  size=10),
          axis.text.y = element_text(face="bold",size=10),
          axis.line = element_line(colour = "darkblue", size = 2, linetype = "solid"),
          axis.title = element_text(size = rel(1.5),face="bold"),
          plot.title = element_text(color="black", size=17, face="bold.italic",hjust = 0.5,vjust =-1),
          legend.title = element_text(size=12, color = "black", face="bold"),
          legend.text = element_text(colour="black", size=12,face="bold"),
          aspect.ratio=1  #square plot
    )+labs(title="PDAC (fibroblast-enriched)")+ labs(fill = "log10(Exp)")
  
  plot_cells(cds2, group_cells_by="cluster", color_cells_by="PDACFE",
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
                                                          "4"="10.apCAF",
                                                          "7"="10.apCAF",
                                                          "24"="10.apCAF",
                                                          "28"="10.apCAF",                                                        
                                                          "29"="10.apCAF",
                                                          "6"="10.apCAF",
                                                          "11"="10.apCAF",
                                                          
                                                          "2"="11.myCAF",
                                                          "27"="11.myCAF",
                                                          "5"="11.myCAF",
                                                          "9"="11.myCAF",
                                                          "16"="11.myCAF",
                                                          "20"="11.myCAF",
                                                          "18"="11.myCAF",
                                                          "25"="11.myCAF",
                                                          "26"="11.myCAF",
                                                          
                                                          "23"="09.iCAF",
                                                          "15"="09.iCAF",
                                                          "17"="09.iCAF",
                                                          "3"="09.iCAF",
                                                          "8"="09.iCAF",
                                                          "19"="09.iCAF",
                                                          "13"="09.iCAF",
                                                          "22"="09.iCAF",
                                                          "12"="09.iCAF",
                                                          "10"="09.iCAF",
                                                          "14"="09.iCAF",
                                                          "1"="09.iCAF",
                                                          "21"="09.iCAF",
                                                          
                                                          "30"="08.Lipofibroblasts"
  )
  
  cds_subsetV2@colData@listData[["Fib"]] <- cds_subsetV2@colData@listData[["assigned_cell_type"]]
  
  #
  plot_cells(cds_subsetV2, group_cells_by="cluster", color_cells_by="assigned_cell_type",
             reduction_method = c( "tSNE"),cell_size = 2,group_label_size = 5,label_cell_groups = FALSE)+
    theme(axis.text.x = element_text(face="bold",  size=10),
          axis.text.y = element_text(face="bold",size=10),
          axis.line = element_line(colour = "darkblue", size = 2, linetype = "solid"),
          axis.title = element_text(size = rel(1.5),face="bold"),
          plot.title = element_text(color="black", size=17, face="bold.italic",hjust = 0.5,vjust =-1),
          legend.title = element_text(size=12, color = "black", face="bold"),
          legend.text = element_text(colour="black", size=12,face="bold"),
          aspect.ratio=1  #square plot
    )+labs(title="CAF")+ labs(fill = "log10(Exp)")
  
  #
  plot_cells(cds_subsetV2, group_cells_by="cluster", color_cells_by="assigned_cell_type",
             reduction_method = c( "tSNE"),cell_size = 2,group_label_size = 5)
  
#### Isolate_subset #### 
  
  ### Now we can transfer the annotations from the cds_subset object back to the full dataset. We'll also filter out low-quality cells at this stage ###
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
    )+labs(title="PDAC (fibroblast-enriched)")+ labs(fill = "log10(Exp)")
  
  
  plot_cells(cds2, group_cells_by="cluster", color_cells_by="assigned_cell_type",
             reduction_method = c( "tSNE"),cell_size = 2,group_label_size = 5)
  
  
  ### Main
  cds2_Violin <- cds2
  cds2_Violin <- cds2_Violin[rowData(cds2_Violin)$gene_short_name %in% ciliated_genes,]
  
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
  
#### Isolate_CAF ####  
  cds_subset_CAFOri <- choose_cells(cds2,reduction_method = c( "tSNE"),
                                  clear_cds = FALSE,
                                  return_list = FALSE)
  cds_subset_CAF <- cds_subset_CAFOri
  
  cds_subset_CAF = cluster_cells(cds_subset_CAF, resolution=1e-2,K=2)
  # plot_cells(cds_subset_CAF, color_cells_by="cluster",reduction_method = c( "tSNE"),cell_size = 2,group_label_size = 5)
  
  plot_cells(cds_subset_CAF, color_cells_by="cluster",reduction_method = c( "tSNE"),cell_size = 2,label_cell_groups = FALSE)
  plot_cells(cds_subset_CAF, color_cells_by="cluster",reduction_method = c( "tSNE"),cell_size = 2,group_label_size = 5)
  
  plot_cells(cds_subset_CAF, color_cells_by="cluster",cell_size = 2,group_label_size = 5,label_cell_groups = FALSE)
  plot_cells(cds_subset_CAF, color_cells_by="cluster",cell_size = 2,group_label_size = 5)
  colData(cds_subset_CAF)$assigned_cell_type <- as.character(clusters(cds_subset_CAF)[colnames(cds_subset_CAF)])
  colData(cds_subset_CAF)$assigned_cell_type <- dplyr::recode(colData(cds_subset_CAF)$assigned_cell_type,
                                                            "3"="10.apCAF",
                                                            "9"="10.apCAF",
                                                            "23"="10.apCAF",
                                                            "25"="10.apCAF",                                                        
                                                            "5"="10.apCAF",
                                                            "11"="10.apCAF",
                                                            
                                                            "16"="11.myCAF",
                                                            "22"="11.myCAF",
                                                            "29"="11.myCAF",
                                                            "1"="11.myCAF",
                                                            "4"="11.myCAF",
                                                            "10"="11.myCAF",
                                                            "26"="11.myCAF",
                                                            "18"="11.myCAF",
                                                            "24"="11.myCAF",
  
                                                            "27"="09.iCAF",                                                          
                                                            "13"="09.iCAF",
                                                            "14"="09.iCAF",
                                                            "2"="09.iCAF",
                                                            "19"="09.iCAF",
                                                            "12"="09.iCAF",
                                                            "21"="09.iCAF",
                                                            "6"="09.iCAF",
                                                            "7"="09.iCAF",
                                                            "20"="09.iCAF",
                                                            "8"="09.iCAF",
                                                            "30"="09.iCAF",
                                                            "15"="09.iCAF",
                                                            "17"="09.iCAF",
                                                            "28"="09.iCAF"
  )
  
  
  cds_subset_CAF@colData@listData[["CAF"]] <- cds_subset_CAF@colData@listData[["assigned_cell_type"]]
  
  #
  plot_cells(cds_subset_CAF, group_cells_by="cluster", color_cells_by="assigned_cell_type",
             reduction_method = c( "tSNE"),cell_size = 2,group_label_size = 5,label_cell_groups = FALSE)+
    theme(axis.text.x = element_text(face="bold",  size=10),
          axis.text.y = element_text(face="bold",size=10),
          axis.line = element_line(colour = "darkblue", size = 2, linetype = "solid"),
          axis.title = element_text(size = rel(1.5),face="bold"),
          plot.title = element_text(color="black", size=17, face="bold.italic",hjust = 0.5,vjust =-1),
          legend.title = element_text(size=12, color = "black", face="bold"),
          legend.text = element_text(colour="black", size=12,face="bold"),
          aspect.ratio=1  #square plot
    )+labs(title="CAF")+ labs(fill = "log10(Exp)")
  
  #
  plot_cells(cds_subset_CAF, group_cells_by="cluster", color_cells_by="assigned_cell_type",
             reduction_method = c( "tSNE"),cell_size = 2,group_label_size = 5)
  
  ###
  cds_subset_Violin_V2 <- cds_subset_CAF
  
  ### Main
  cds_subset_Violin_V2 <- cds_subset_Violin_V2[rowData(cds_subset_Violin_V2)$gene_short_name %in% ciliated_genes,]
  
  plot_genes_violin(cds_subset_Violin_V2, group_cells_by="CAF",  log_scale = FALSE, ncol=3) +
    theme(axis.text.x=element_text(angle=45, hjust=1))+
    theme(axis.text.x = element_text(face="bold",  size=12),
          axis.text.y = element_text(face="bold",size=12),
          axis.line = element_line(colour = "darkblue", size = 2, linetype = "solid"),
          axis.title = element_text(size = rel(1.5),face="bold"),
          plot.title = element_text(color="black", size=17, face="bold.italic",hjust = 0.5,vjust =-1),
          legend.title = element_text(size=12, color = "black", face="bold"),
          legend.text = element_text(colour="black", size=12,face="bold"),
          aspect.ratio=1  #square plot
    )+labs(title="Main")+ labs(fill = "log10(Exp)")
  
  ### iCAF
  cds_subset_Violin_iCAF_V2 <- cds_subset_CAF
  
  cds_subset_Violin_iCAF_V2 <- cds_subset_Violin_iCAF_V2[rowData(cds_subset_Violin_iCAF_V2)$gene_short_name %in% iCAF_Main,]
  
  plot_genes_violin(cds_subset_Violin_iCAF_V2, group_cells_by="CAF",  log_scale = FALSE, ncol=2) +
    theme(axis.text.x=element_text(angle=45, hjust=1))+
    theme(axis.text.x = element_text(face="bold",  size=10),
          axis.text.y = element_text(face="bold",size=10),
          axis.line = element_line(colour = "darkblue", size = 2, linetype = "solid"),
          axis.title = element_text(size = rel(1.5),face="bold"),
          plot.title = element_text(color="black", size=17, face="bold.italic",hjust = 0.5,vjust =-1),
          legend.title = element_text(size=12, color = "black", face="bold"),
          legend.text = element_text(colour="black", size=12,face="bold"),
          aspect.ratio=1  #square plot
    )+labs(title="iCAF markers")+ labs(fill = "log10(Exp)")
  
  ### myCAF
  cds_subset_Violin_myCAF_V2 <- cds_subset_CAF
  
  cds_subset_Violin_myCAF_V2 <- cds_subset_Violin_myCAF_V2[rowData(cds_subset_Violin_myCAF_V2)$gene_short_name %in% myCAF_Main,]
  
  plot_genes_violin(cds_subset_Violin_myCAF_V2, group_cells_by="CAF",  log_scale = FALSE, ncol=2) +
    theme(axis.text.x=element_text(angle=45, hjust=1))+
    theme(axis.text.x = element_text(face="bold",  size=10),
          axis.text.y = element_text(face="bold",size=10),
          axis.line = element_line(colour = "darkblue", size = 2, linetype = "solid"),
          axis.title = element_text(size = rel(1.5),face="bold"),
          plot.title = element_text(color="black", size=17, face="bold.italic",hjust = 0.5,vjust =-1),
          legend.title = element_text(size=12, color = "black", face="bold"),
          legend.text = element_text(colour="black", size=12,face="bold"),
          aspect.ratio=1  #square plot
    )+labs(title="myCAF markers")+ labs(fill = "log10(Exp)")
  
  
  ### apCAF
  cds_subset_Violin_apCAF_V2 <- cds_subset_CAF
  
  cds_subset_Violin_apCAF_V2 <- cds_subset_Violin_apCAF_V2[rowData(cds_subset_Violin_apCAF_V2)$gene_short_name %in% apCAF_Main,]
  
  plot_genes_violin(cds_subset_Violin_apCAF_V2, group_cells_by="CAF",  log_scale = FALSE, ncol=2) +
    theme(axis.text.x=element_text(angle=45, hjust=1))+
    theme(axis.text.x = element_text(face="bold",  size=10),
          axis.text.y = element_text(face="bold",size=10),
          axis.line = element_line(colour = "darkblue", size = 2, linetype = "solid"),
          axis.title = element_text(size = rel(1.5),face="bold"),
          plot.title = element_text(color="black", size=17, face="bold.italic",hjust = 0.5,vjust =-1),
          legend.title = element_text(size=12, color = "black", face="bold"),
          legend.text = element_text(colour="black", size=12,face="bold"),
          aspect.ratio=1  #square plot
    )+labs(title="apCAF markers")+ labs(fill = "log10(Exp)")
