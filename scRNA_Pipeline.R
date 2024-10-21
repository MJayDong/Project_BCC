
setwd('/GPFS/Magda_lab_temp/dongmingjie/Project_MJ/Project_1/scRNA')

scRNA <- readRDS('scRNA_T.rds')

scRNA;table(scRNA$ident);table(scRNA$preorpost);table(scRNA$patient);table(scRNA$celltype_citation);table(scRNA$cellcondition)

# celltype_citation
MYCOLOR          <- c(
                      "#91D1C24C","#DC00004C","#fcec7c","#ff92e0","#E64B357F","#4DBBD57F","#00A0877F","#3C54887F",
                      "#F39B7F7F","#8491B47F","#91D1C27F","#DC00007F","#eae8c1","#ef4db1","#E64B35B2","#4DBBD5B2",
                      "#00A087B2","#3C5488B2","#F39B7FB2","#8491B4B2","#91D1C2B2","#DC0000B2","#edd600","#c64ea4",
                      "#E64B35FF","#4DBBD5FF","#00A087FF","#3C5488FF","#F39B7FFF","#8491B4FF","#91D1C2FF","#DC0000FF",
                      "#d1cc78","#b2387e","#E64B354C","#4DBBD54C","#00A0874C","#3C54884C","#F39B7F4C","#8491B44C")

# celltype
MYCOLOR_1        <- c('#D51F26','#272E6A','#208A42','#89288F','#F47D2B','#FEE500','#8A9FD1','#C06CAB','#D8A767')
names(MYCOLOR_1) <- c('CD8_act','CD8_eff','CD8_ex','CD8_ex_act','CD8_mem','Naive','Tfh','Th17','Tregs')
# patient
MYCOLOR_2        <- c("#A1A9D0","#F0988C","#B883D4","#CFEAF1","#C4A5DE","#F6CAE5","#96CCCB","#FFBE7A","#82B0D2","#BEB8DC","#E7DAD2")
names(MYCOLOR_2) <- c('su001','su002','su003','su004','su005','su006','su007','su008','su009','su010','su012')
# preorpost
MYCOLOR_3        <- c("#F8766D","#00BFC4")
names(MYCOLOR_3) <- c('pre','post')

MYCOLOR_1;MYCOLOR_2;MYCOLOR_3
# 0.Preprocessing ---------------------------------------------------------

#library(harmony)

library(dplyr)
library(Seurat)
library(SingleR)
library(patchwork)
library(tidyverse)
library(data.table)
library(sctransform) 
library(SeuratObject)

p1   <- VlnPlot(scRNA,group.by = 'ident',
               features = c("nFeature_RNA", "nCount_RNA", "Percent_MT"), 
               cols = MYCOLOR, pt.size = 0, ncol = 1)
ggsave(p1, file='Figure/1.QC.pdf', width=8, height=12)


# 1.DEAnalysis ------------------------------------------------------------

library(Seurat)
library(ggplot2)
library(ggrepel)
library(stringr)   
library(dplyr)
library(pheatmap)
library(ComplexHeatmap)
library(circlize)
library(plyr)
library(tidyverse)
library(AnnotationDbi)
library(org.Hs.eg.db)#基因注释包
library(clusterProfiler)#富集包
library(msigdbr)
library(BoutrosLab.plotting.general)
library(gdata)
library(VennDiagram)
library(dendextend)
library(igraph)
#library(qusage)
library(readxl)
library(lattice)
library(ragg)
library(reshape2)

library(Seurat)
library(ggplot2)
library(cowplot)
library(dplyr)
#library(ggunchull)
#library(tidydr)
library(ggsci)
library(Cairo)

library(Seurat)
library(reshape2)
library(tidyverse)
library(ggplot2)
library(ggsignif)
library(dplyr)
library(viridis)
library(ComplexHeatmap)
library(genefilter)


#celltype
p2 <- DimPlot(scRNA,group.by='celltype_citation',reduction = "umap",label = T,label.size = 4,repel = T,cols = MYCOLOR_1,label.box = T,pt.size = 1)
#样本来源
p3 <- DimPlot(scRNA,group.by='patient'          ,reduction = "umap",label = T,label.size = 4,repel = T,cols = MYCOLOR_2,label.box = T,pt.size = 1)
#治疗前后来源
p4 <- DimPlot(scRNA,group.by='preorpost'        ,reduction = "umap",label = T,label.size = 4,repel = T,cols = MYCOLOR_3,label.box = T,pt.size = 1)

ggsave(p2, file='Figure/2.umap_celltype_citation.pdf', width=8, height=8)
ggsave(p3, file='Figure/2.umap_patient.pdf'          , width=8, height=8)
ggsave(p4, file='Figure/2.umap_preorpost.pdf'        , width=8, height=8)

markerGenes <- c(
  'CD3G', 'CD3D', 'CD3E',			# T cells
  
  'CD4',											# CD4+ T cells
  'IL2RA','FOXP3','CTLA4',			  # Treg cells
  'TBX21','IFNG',									# Th1 cells
  'IL4','GATA3',									# Th2 cells
  'IL9','SPI1',									  # Th9 cells
  'TGFB1','RORC','IL17A','IL17F','IL26','CCR6','KLRB1','CTSH',	# Th17
  'IL22','AHR','IL13',						# Th22
  'CXCR5','IL21','BCL6','CD200','IL26','CD200','PTPN13','BTLA',	# Tfh
  
  'CD8A', 'GZMA',									# CD8+ T cells
  'CD69','IL2RA','TNF','IFNG','FOS','JUN',	# Activated
  'GZMB','PRF1',								  # Cytotoxic
  'PTPRC','SELL','CCR7',					# Naive+++/Effector、Memory+--
  'EOMES','GZMK','CXCR3',					# Memory
  'KLRD1',									      # Effect Memory
  'IL7R'  ,'CCR7',							  # Naive CD4+ T
  'IL7R'  ,'S100A4',							# Memory CD4+ T
  'PDCD1','HAVCR2',							  # T cell exhaustion genes
  'IFNG','TNF' ,								  # T cell activation genes
  'ENTPD1', 'ITGAE',							# 'CD39','CD103' Exhausted CD8+ T cells
  'TCF7','PDCD1','HAVCR2',				# Progenitor Exhausted Tcell
  'TOX','TOX2','PDCD1','TIGIT',		# Terminally Exhausted Tcell
  

  'EPCAM','BCAM','TP63',				# BCC (malignant cells)
  
  # 5 resident markers 
  'RUNX3', 'NR4A1', 'CD69', 'CXCR6','NR4A3',
  # 7 cytotoxicity associated genes
  'IFNG', 'GNLY', 'NKG7', 'GZMB', 'GZMA', 'CST7','TNFSF10',
  # 5 exhausted markers 
  'CTLA4', 'HAVCR2', 'LAG3', 'PDCD1', 'TIGIT',
  # 6 costimulatory molecular genes
  'ICOS', 'CD226', 'TNFRSF14', 'TNFRSF25', 'TNFRSF9', 'CD28'
  
  #'CD14'  ,'LYZ',				  # CD14+ Mono
  #'FCGR3A','MS4A7',			  # FCGR3A+ Mono
  #'FAP', 'PDPN',				    # cancer-associated fibroblasts
  #'FCGR2A', 'CSF1R',			  # macrophages
  #'FLT3','FCER1A','CST3' ,	# dendritic cells
  #'PPBP',						      # Platelet
  #'CD34',						      # Early Progenitor
  #'GATA1'						      # Erythrocyte
)


inferno_mod <- viridis(20)[3:20]

p5 <- DotPlot(scRNA,  group.by='celltype_citation', features = rev(unique(markerGenes[1:60])),  
              assay = "RNA", dot.scale = 4.5) + 
              scale_color_gradient2(low = inferno_mod[3], 
                                    mid = "#F1ED6FFF", high = "#D64B40FF", midpoint = 0.25) +
              theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1, size = 12),
                    axis.line = element_line(size = 0.3), axis.ticks = element_line(size = 0.3)) + 
              coord_flip()


p6 <- DotPlot(scRNA,  group.by='celltype_citation', features = rev(unique(markerGenes)),  
              assay = "RNA", dot.scale = 4.5) + 
              scale_color_gradient2(low = inferno_mod[3], 
                                    mid = "#F1ED6FFF", high = "#D64B40FF", midpoint = 0.25) +
              theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1, size = 12),
                    axis.line = element_line(size = 0.3), axis.ticks = element_line(size = 0.3)) + 
              coord_flip()

ggsave(p5, file='Figure/3.MarkerGenes_scRNA.pdf', width=8, height=10)
ggsave(p6, file='Figure/3.AllGenes_scRNA.pdf'   , width=8, height=14)


library(Seurat)
library(ggplot2)
p7 <- FeaturePlot(scRNA,features = unique(markerGenes) , cols = c("lightgrey", "red"),min.cutoff = 0, max.cutoff = 3)
ggsave(p7, file='Figure/7.AllGenes_FeaturePlot.pdf'   , width=5, height=5)
# 2.Pseudotime ------------------------------------------------------------


# 3.CellphoneDB -------------------------------------------------------------


# 4.SCENIC ----------------------------------------------------------------


# 5.Addition --------------------------------------------------------------


