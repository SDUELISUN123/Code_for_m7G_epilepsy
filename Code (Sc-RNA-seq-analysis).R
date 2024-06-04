#>Tittle:
#>Integrative Role of RNA N7-methylguanosine in Epilepsy: 
#>Regulation of Neuronal Oxidative Phosphorylation, 
#>Programmed Death and Immune Microenvironment
#>
#>
#>Author:
#>Qingyuan Sun (孙清源)
#>
#>
#>Institute:
#>Shandong University
#>
#>
#>Date:
#>June 4th, 2024




library(Seurat)
library(dplyr)
library(future)
library(future.apply)
library(dplyr)
library(msigdbr)
#library(clusterProfiler)
library(tidyverse)
library(patchwork)
library(monocle)
library(ggpubr)

setwd("")
cors<-ggsci::pal_igv()(51)

nbrOfWorkers()
plan("multicore", workers = 3) ###set the compute core
options(future.globals.maxSize = 14000 * 1024^2)


#### 1 data reading ####
N1<-readRDS("Normal1.rds")
N2<-readRDS("Normal2.rds")
N3<-readRDS("Normal3.rds")
N4<-readRDS("Normal4.rds")

T1<-readRDS("TLE1.rds")
T2<-readRDS("TLE2.rds")
T3<-readRDS("TLE3.rds")
T4<-readRDS("TLE4.rds")

sqy<-cbind(N1,N2,N3,N4,T1,T2,T3,T4)
rm(N1,N2,N3,N4,T1,T2,T3,T4)
gc()



#### 2 processing data ####

sqy=CreateSeuratObject(counts = sqy,project = "seurat", min.cells=3, min.features=50, names.delim = "_")
gc()
dir.create("3 seurat")
setwd("./3 seurat/")
saveRDS(sqy,"initial.rds")
beepr::beep(1)
gc()


library(devtools)
#install_github("immunogenomics/harmony")
library(harmony)
gc()

#sqy@meta.data[1:5,]

#使用PercentageFeatureSet函数计算线粒体基因的百分比
sqy[["percent.mt"]] <- PercentageFeatureSet(object = sqy, pattern = "^MT-")#鼠源的换成mt

#质控
VlnPlot(object = sqy, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,group.by = "orig.ident",pt.size = 0,cols = cors)  #ncol是指每行放几张图
sqy=subset(x = sqy, subset = nFeature_RNA > 50 & percent.mt <0.1)    #对数据进行过滤
gc()

#测序深度的相关性图
plot1 <- FeatureScatter(object = sqy, feature1 = "nCount_RNA", feature2 = "percent.mt",pt.size=1.5,cols = cors)
plot2 <- FeatureScatter(object = sqy, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",pt.size=1.5,cols = cors)
CombinePlots(plots = list(plot1, plot2))
rm(plot1,plot2)

#看细胞属于哪一期，并加到矩阵里
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
sqy <- CellCycleScoring(sqy, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
sqy@meta.data[1:5,]
######标准化
sqy<-NormalizeData(sqy,verbose = T)   #标准化

####PCA
sqy<-FindVariableFeatures(sqy,selection.method = "vst", nfeatures = 2000)   #找前2000个差异显著的基因
sqy<-ScaleData(sqy,vars.to.regress = c("percent.mt","S.Score","G2M.Score"),verbose = T) #去除线粒体基因和分裂期的影响
sqy<-RunPCA(sqy,verbose = T,npcs = 70)  #pca降维
ElbowPlot(sqy,ndims = 70)  #看拐点
#选定pc30

#矫正前的降维图和vln
p1 <- DimPlot(object = sqy, reduction = "pca", pt.size = .1, group.by = "orig.ident")
p2 <- VlnPlot(object = sqy, features = "PC_1", group.by = "orig.ident", pt.size = .1)
p1|p2

#矫正
sqy<-RunHarmony(sqy,group.by.vars = c("orig.ident"), plot_convergence = TRUE)
harmony_embeddings <- Embeddings(sqy, 'harmony')
#dim(harmony_embeddings)



#矫正后的降维图和vln
p3 <- DimPlot(object = sqy, reduction = "harmony", pt.size = .1, group.by = "orig.ident")
p4 <- VlnPlot(object = sqy, features = "harmony_1", group.by = "orig.ident", pt.size = .1)
p3|p4


p1+p3#去除批次前后的降维图

p2|p4#去批次前后的vlnplot


rm(p1,p2,p3,p4,harmony_embeddings,g2m.genes,s.genes,plot1,plot2)


###2.5 计算通路
#文件输入
pathways<-read.table("gene.txt",sep = "\t",header = T)
pathways<-as.list(pathways)
sqy<-AddModuleScore(sqy,features = pathways,ctrl = 100,name = "m7G_Score")
rm(pathways,i,gene)

#序号输入
library("KEGGREST") 
listDatabases()  
gs<-keggGet('hsa00190')
#获取通路中gene信息 
gs[[1]]$GENE 
#查找所有基因 
genes<-unlist(lapply(gs[[1]]$GENE,function(x) strsplit(x,';'))) 
pathways <- genes[1:length(genes)%%3 ==2] 
pathways <- data.frame(pathways)  
rm(gs,genes)
gene<-as.list(pathways)
sqy<-AddModuleScore(object = sqy,features = gene,name = "Oxidative_phosphorylation",ctrl = 100)#更改path名字

rm(pathways,i,gene)

#批量处理
pathways<-read.table("celldeathmode.txt",header = T,sep = "\t")
for (i in 1:ncol(pathways)) {
  genes<-as.data.frame(pathways[,i])
  colnames(genes)<-colnames(pathways)[i]
  genes<-subset(genes,genes[,1]!="")
  genes<-as.list(genes)
  path_name<-colnames(pathways)[i]
  sqy<-AddModuleScore(sqy,features = genes,ctrl = 100,name = path_name)
}
rm(genes,my.comparasion,pathways,i,path_name)

##
table(sqy$orig.ident)
sqy@active.ident<-sqy$orig.ident
sqy<-RenameIdents(sqy,
                  'Normal1'='Normal',
                  'Normal2'='Normal',
                  'Normal3'='Normal',
                  'Normal4'='Normal',
                  'TLE1'='TLE',
                  'TLE2'='TLE',
                  'TLE3'='TLE',
                  'TLE4'='TLE')
sqy$group<-sqy@active.ident


#### 3 dimention reduction and annotation ####
#UMAP/tsne
sqy <- FindNeighbors(object = sqy, dims = 1:30)       #计算邻接距离
sqy <- FindClusters(object = sqy, resolution = 1)         #对细胞分组,对细胞标准模块化,resolution在0.4-1.2，越大分的群越多
sqy <- RunTSNE(object = sqy, dims = 1:30) 

pdf(file = "TSNE_cluster.pdf",width=6,height = 5)
TSNEPlot(object = sqy, label = TRUE, cols=cors)  #整体的umap
dev.off()

Excitatory_neurons=c('SYN3','RBFOX3','CAMK2A')
Inhibitory_neurons=c('ERBB4','NXPH1','GAD2')
Astrocytes=c('GFAP','FGFR3','GJA1','AQP4','ALDH1L1')
OPCs=c('PCDH15','LHFPL3','VCAN','PDGFRA','CSPG4')
Oligodendrocytes=c('MOBP','PLP1','MBP')
Microglias=c('DOCK8','SFMBT2','CX3CR1','C1QC')
Macrophages=c('CD68','IL1B','CD14')
ECs=c('FLT1','VWF','PECAM1','CDH5')

DotPlot(sqy,features = c(Excitatory_neurons,
                         Inhibitory_neurons,
                         Astrocytes,
                         OPCs,
                         Oligodendrocytes,
                         Microglias,
                         Macrophages,
                         ECs),
        cols = c("grey","red"),
        group.by = "seurat_clusters") + NoLegend()  + RotatedAxis()
sqy@active.ident<-sqy$seurat_clusters
sqy<-RenameIdents(sqy,
                  '0'='Oligodendrocytes',
                  '1'='OPCs',
                  '2'='Inhibitory_neurons',
                  '3'='Excitatory_neurons',
                  '6'='Excitatory_neurons',
                  '7'='Excitatory_neurons',
                  '10'='Excitatory_neurons',
                  '13'='Excitatory_neurons',
                  '17'='Excitatory_neurons',
                  '19'='Excitatory_neurons',
                  '20'='Excitatory_neurons',
                  '21'='Excitatory_neurons',
                  '22'='Excitatory_neurons',
                  '24'='Excitatory_neurons',
                  '25'='Excitatory_neurons',
                  '28'='Excitatory_neurons',
                  '29'='Excitatory_neurons',
                  '31'='Excitatory_neurons',
                  '4'='Inhibitory_neurons',
                  '9'='Inhibitory_neurons',
                  '12'='Inhibitory_neurons',
                  '14'='Inhibitory_neurons',
                  '34'='Inhibitory_neurons',
                  '30'='Astrocytes',
                  '23'='Astrocytes',
                  '18'='Astrocytes',
                  '15'='Astrocytes',
                  '11'='Astrocytes',
                  '8'='Astrocytes',
                  '26'='OPCs',
                  '35'='OPCs',
                  '16'='Oligodendrocytes',
                  '27'='Macrophages/Monocytes',
                  '5'='Microglias',
                  '35'='Microglias',
                  '33'='ECs',
                  '32'='Astrocytes')
sqy$celltype<-sqy@active.ident
table(sqy$celltype)
DotPlot(sqy,features = c(ECs,Microglias,Macrophages,Oligodendrocytes,OPCs,Astrocytes,Inhibitory_neurons,Excitatory_neurons),
        cols = c("grey","red"),
        group.by = "celltype") + RotatedAxis()

pdf(file = "TSNE_cluster_celltype.pdf",width=7,height = 5)
TSNEPlot(object = sqy, label = TRUE, cols=cors,group.by="celltype",repel=T)  #整体的umap
dev.off()



### 4 m7G Activity ####
pdf('featureplot_m7G得分.pdf',width = 5,height = 5)
sqy@active.ident<-sqy$celltype
FeaturePlot(sqy,features = "m7G_Score1",label = T,cols = c('grey','red'),repel = T)
dev.off()

pdf('vln_m7G得分.pdf',width = 4,height = 5)
VlnPlot(sqy,features = 'm7G_Score1',sort = T,cols = cors,pt.size = 0,group.by = "celltype")+
  NoLegend()+
  stat_compare_means(label.x = 2)+
  geom_boxplot()
dev.off()

pdf('vln_m7G得分_各位患者.pdf',width = 4,height = 4)
VlnPlot(sqy[,sqy$group=="TLE"],
        features = 'm7G_Score1',
        cols = cors,
        pt.size = 0,
        sort = T,
        group.by = "orig.ident")+
  stat_compare_means()+
  geom_boxplot()+
  NoLegend()
dev.off()

pdf('vln_m7G得分_各位患者_神经元.pdf',width = 4,height = 4)
VlnPlot(sqy[,sqy$group=="TLE"&sqy$celltype%in%c("Inhibitory_neurons","Excitatory_neurons")],
        features = 'm7G_Score1',
        cols = cors,
        pt.size = 0,
        sort = T,
        group.by = "orig.ident")+
  stat_compare_means()+
  geom_boxplot()+
  NoLegend()
dev.off()

sqy@active.ident<-sqy$orig.ident
sqy<-RenameIdents(sqy,
                  'Normal1'='Normal',
                  'Normal2'='Normal',
                  'Normal3'='Normal',
                  'Normal4'='Normal',
                  'TLE1'='high_m7G_TLE',
                  'TLE2'='low_m7G_TLE',
                  'TLE3'='low_m7G_TLE',
                  'TLE4'='low_m7G_TLE')
sqy$m7G_group<-sqy@active.ident


m7Ggenes<-read.table("gene.txt",header = T)[,1]


### 5 m7G regulators' expression pattern ####
pdf('dotplot_m7G基因在各类细胞表达.pdf',width = 10.5,height = 3.5)
DotPlot(sqy,
        features = c(m7Ggenes),
        group.by = "celltype",
        cols = c('grey','red'))+
  RotatedAxis()
dev.off()



#### 6 Different functions between high and low m7G TLE ####
#细胞死亡
celldeathmode<-read.table("celldeathmode.txt",header = T,sep = "\t")
celldeath_path<-paste0(colnames(celldeathmode),"1")
rm(celldeathmode)

dir.create("m7G高低患者之间细胞死亡")
setwd("./m7G高低患者之间细胞死亡/")

for (i in 1:12) {
  filename<-paste0("vln_",celldeath_path[i],"_m7G高低之间",".pdf")
  pdf(file = filename,width = 3,height = 5)
  VlnPlot(sqy[,sqy$group=="TLE"],
          features = celldeath_path[i],
          pt.size = 0,
          sort = F,
          group.by = "m7G_group",
          cols = c('pink','lightblue'))+
    geom_boxplot()+
    stat_compare_means()+
    NoLegend()
  dev.off()
}

celldeath<-
  DotPlot(sqy[,sqy$group=="TLE"&sqy$celltype%in%c("Inhibitory_neurons","Excitatory_neurons")],
          features = c('Autophagy1','Ferroptosis1','Necroptosis1','Pyroptosis1'),
          group.by = "orig.ident")+
  RotatedAxis()

celldeath_vln_autophagy_neurons<-
  VlnPlot(sqy[,sqy$group=="TLE"&sqy$celltype%in%c("Inhibitory_neurons","Excitatory_neurons")],
          features = c('Autophagy1'),
          pt.size = 0,
          cols=c("pink","lightblue"),
          group.by = 'm7G_group')+stat_compare_means(method = "wilcox")+geom_boxplot()+
  NoLegend()

celldeath_vln_Ferroptosis_neurons<-
  VlnPlot(sqy[,sqy$group=="TLE"&sqy$celltype%in%c("Inhibitory_neurons","Excitatory_neurons")],
          features = c('Ferroptosis1'),
          pt.size = 0,
          cols=c("pink","lightblue"),
          group.by = 'm7G_group')+stat_compare_means(method = "wilcox")+geom_boxplot()+
  NoLegend()

celldeath_vln_Necroptosis_neurons<-
  VlnPlot(sqy[,sqy$group=="TLE"&sqy$celltype%in%c("Inhibitory_neurons","Excitatory_neurons")],
          features = c('Necroptosis1'),
          pt.size = 0,
          cols=c("pink","lightblue"),
          group.by = 'm7G_group')+stat_compare_means(method = "wilcox")+geom_boxplot()+
  NoLegend()

celldeath_vln_Pyroptosis_neurons<-
  VlnPlot(sqy[,sqy$group=="TLE"&sqy$celltype%in%c("Inhibitory_neurons","Excitatory_neurons")],
          features = c('Pyroptosis1'),
          pt.size = 0,
          cols=c("pink","lightblue"),
          group.by = 'm7G_group')+stat_compare_means(method = "wilcox")+geom_boxplot()+
  NoLegend()

pdf("celldeath.vln.celldeath4.neurons.pdf",width = 6,height = 10)
(celldeath_vln_autophagy_neurons|celldeath_vln_Ferroptosis_neurons)/
  (celldeath_vln_Necroptosis_neurons|celldeath_vln_Pyroptosis_neurons)
dev.off()

#op(全部细胞)#

pdf(file = "vln_op.pdf",width = 3,height = 5)
VlnPlot(sqy[,sqy$group=="TLE"],
        features = "Oxidative_phosphorylation1",
        pt.size = 0,
        sort = F,
        group.by = "m7G_group",
        cols = c('pink','lightblue'))+
  geom_boxplot()+
  stat_compare_means()+
  NoLegend()
dev.off()

library("KEGGREST") 
listDatabases()  
gs<-keggGet('hsa00190')
gs[[1]]$GENE 
genes<-unlist(lapply(gs[[1]]$GENE,function(x) strsplit(x,';'))) 
pathways <- genes[1:length(genes)%%3 ==2] 
pathways <- data.frame(pathways)  

VlnPlot(sqy[,sqy$group=="TLE"],
        features = "NDUFS4",
        pt.size = 0,
        sort = F,
        group.by = "m7G_group",
        cols = c('pink','lightblue'))+
  geom_boxplot()+
  stat_compare_means()+
  NoLegend()

NDUFS_family<-c("NDUFS1","NDUFS2","NDUFS3","NDUFS4","NDUFS5",
                "NDUFS6","NDUFS7","NDUFS8")
NDUFV_family<-c('NDUFV1','NDUFV2','NDUFV3')
NDUFA_family<-c('NDUFA1','NDUFA2','NDUFA3','NDUFA4','NDUFA5','NDUFA6','NDUFA7','NDUFA8','NDUFA9','NDUFA10','NDUFAB1','NDUFA11','NDUFA12','NDUFA13')
NUDFB_family<-c('NDUFB1','NDUFB2','NDUFB3','NDUFB4','NDUFB5','NDUFB6','NDUFB7','NDUFB8','NDUFB9','NDUFB10','NDUFB11')
NUDFC_family<-c('NDUFC1','NDUFC2')
NADH_dehydrogenase<-
  DotPlot(sqy[,sqy$group=="TLE"],
          features = c(NDUFS_family,NDUFV_family,NDUFA_family,NUDFB_family,NUDFC_family),
          group.by = "orig.ident")+
  RotatedAxis()#

Cyto_C_reduc_gene<-c('UQCRFS1','UQCRC1','UQCRC2','UQCRH','UQCRB','UQCRQ','UQCR10','UQCR11')
Cyto_C_reduc<-
  DotPlot(sqy[,sqy$group=="TLE"],
          features = c(Cyto_C_reduc_gene),
          group.by = "orig.ident")+
  RotatedAxis()

Cyto_C_oxida_gene<-c('COX4I1','COX5A','COX5B','COX6A1','COX6B1','COX6C','COX7A1','COX7A2','COX7A2L','COX7B','COX7C','COX8A','COX11','COX15','COX17','COX10')
Cyto_C_oxida<-
  DotPlot(sqy[,sqy$group=="TLE"],
          features = c(Cyto_C_oxida_gene),
          group.by = "orig.ident")+
  RotatedAxis()

Cyto_C_CYC<-
  DotPlot(sqy[,sqy$group=="TLE"],
          features = c("CYCS"),
          group.by = "orig.ident")+
  RotatedAxis()

ATPase_gene<-c('ATP6V1A','ATP6V1B2','ATP6V1C1','ATP6V1D','ATP6V1E1','ATP6V1E2','ATP6V1F','ATP6V1G2','ATP6V1G1','ATP6V0A1','ATP6V0A2','ATP6V0B','ATP6V0D1','ATP6V0E1','ATP6V0E2')
ATPase<-
  DotPlot(sqy[,sqy$group=="TLE"],
          features = c(ATPase_gene),
          group.by = "orig.ident")+
  RotatedAxis()

succinate_dehydrogenase_gene<-c('SDHA','SDHB','SDHC','SDHD')
succinate_dehydrogenase<-
  DotPlot(sqy[,sqy$group=="TLE"],
          features = c(succinate_dehydrogenase_gene),
          group.by = "orig.ident")+
  RotatedAxis()

OP<-
  DotPlot(sqy[,sqy$group=="TLE"],
          features = c("Oxidative_phosphorylation1"),
          group.by = "orig.ident")+
  RotatedAxis()

length(c(NDUFS_family,NDUFV_family,NDUFA_family,NUDFB_family,NUDFC_family))
length(succinate_dehydrogenase_gene)
length(Cyto_C_reduc_gene)
length(Cyto_C_oxida_gene)
length(ATPase_gene)
NADH_dehydrogenase/(succinate_dehydrogenase|Cyto_C_reduc|Cyto_C_oxida|Cyto_C_CYC|ATPase)



#op(神经元)


NADH_dehydrogenase_neuron<-
  DotPlot(sqy[,sqy$group=='TLE'&sqy$celltype%in%c('Inhibitory_neurons','Excitatory_neurons')],
          features = c(NDUFS_family,NDUFV_family,NDUFA_family,NUDFB_family,NUDFC_family),
          group.by = "orig.ident")+
  RotatedAxis()+
  NoLegend()


Cyto_C_reduc_neuron<-
  DotPlot(sqy[,sqy$group=='TLE'&sqy$celltype%in%c('Inhibitory_neurons','Excitatory_neurons')],        features = c(Cyto_C_reduc_gene),
          group.by = "orig.ident")+
  RotatedAxis()+
  NoLegend()


Cyto_C_oxida_neuron<-
  DotPlot(sqy[,sqy$group=='TLE'&sqy$celltype%in%c('Inhibitory_neurons','Excitatory_neurons')],          features = c(Cyto_C_oxida_gene),
          group.by = "orig.ident")+
  RotatedAxis()+
  NoLegend()


Cyto_C_CYC_neuron<-
  DotPlot(sqy[,sqy$group=='TLE'&sqy$celltype%in%c('Inhibitory_neurons','Excitatory_neurons')],          features = c("CYCS"),
          group.by = "orig.ident")+
  RotatedAxis()+
  NoLegend()


ATPase_neuron<-
  DotPlot(sqy[,sqy$group=='TLE'&sqy$celltype%in%c('Inhibitory_neurons','Excitatory_neurons')],          features = c(ATPase_gene),
          group.by = "orig.ident")+
  RotatedAxis()+
  NoLegend()


succinate_dehydrogenase_neuron<-
  DotPlot(sqy[,sqy$group=='TLE'&sqy$celltype%in%c('Inhibitory_neurons','Excitatory_neurons')],          features = c(succinate_dehydrogenase_gene),
          group.by = "orig.ident")+
  RotatedAxis()+
  NoLegend()

pdf("dot_op_neurons.pdf",width = 10,height = 8)
(NADH_dehydrogenase_neuron)/
  (Cyto_C_reduc_neuron|Cyto_C_oxida_neuron)/
  ((succinate_dehydrogenase_neuron+Cyto_C_CYC_neuron)|(ATPase_neuron))
dev.off()


save.image("24.2.4.RData")

####### Codes for bulk-RNA-seq analysis were in another file

##>Last version: June 4th, 2024
##>By: Qingyuan Sun
##>Shandong University