#There's no need to run this part yourself, since the raw data from both strains is preserved in the integrated object. 
#This is just to document the steps taken to integrate our two datasets.


setwd("~/Box Sync/Evan/Zhao lab notebook/Misc data/181105_Dmel_sample2/")
library(Seurat)
library(viridis)
library(MASS)
library(plyr)
library(dplyr)
library(sctransform)
library(hues)
library(ggpubr)
library(ggsignif)
library(tidyverse)
library(rstatix)
library(cowplot)



#integrate 2 strains###########
wild  <-Read10X(data.dir = "Dmel_merged/" )
wild <-CreateSeuratObject(wild, min.cells=3, min.features =200 )
wild <- NormalizeData(object = wild)

wild<- ScaleData(object = wild)
wild<-FindVariableFeatures(wild)

wild<-RunPCA(wild)
wild<-NormalizeData(wild)
wild<-ScaleData(wild)
R517 <-readRDS("~/Box Sync/Evan/Zhao lab notebook/Misc data/181217_Dmel_li_denovo/181218_testis_5000_li_denovo_nomagic.rds")
R517<-UpdateSeuratObject(R517)
avgexp<-AverageExpression(R517)
write.table(avgexp,"~/Box Sync/Evan/Zhao lab notebook/Misc data/181217_Dmel_li_denovo/avgexp.txt" )

R517<-NormalizeData(R517)
R517<-ScaleData(R517)


R517$strain<-"R517"
wild$strain<-"wild"



lapply(paste('package:',names(sessionInfo()$otherPkgs),sep=""),detach,character.only=TRUE,unload=TRUE) # unload all
library("Seurat")

#try SCtransform
R517<-SCTransform(R517)
wild<-SCTransform(wild)
testis.list<-list(R517, wild)

testis.features<- SelectIntegrationFeatures(object.list = testis.list, nfeatures = 3000)
testis.list<-PrepSCTIntegration(object.list = testis.list, anchor.features = testis.features)

testis.anchors <- FindIntegrationAnchors(object.list = testis.list, normalization.method = "SCT", 
                                       anchor.features = testis.features)
testis.integrated <- IntegrateData(anchorset = testis.anchors, normalization.method = "SCT", verbose = FALSE)
testis.integrated<- ScaleData(testis.integrated, verbose = FALSE)
testis.integrated<-FindVariableFeatures(testis.integrated)
testis.integrated <- RunPCA(testis.integrated,verbose = FALSE, npcs = 100, assay = "integrated")
testis.integrated<-JackStraw(testis.integrated, reduction='pca', assay = "integrated", dims = 100)
testis.integrated<-ScoreJackStraw(testis.integrated, dims = 1:100)
JackStrawPlot(testis.integrated, dims=1:100)
ggsave("~/Box Sync/Evan/Zhao lab notebook/Misc data/190730_R517_wild_testis_integrated_jackstraw.pdf", width=12, height=8)
PCHeatmap(testis.integrated, assay="integrated", dims=1:20)

testis.integrated<- RunUMAP(testis.integrated, dims=c(1:5, 8, 16, 18, 20), repulsion.strength = 5)
testis.integrated <- FindNeighbors(testis.integrated, reduction = "pca", dims = c(1:5, 8, 16, 18, 20))
testis.integrated <- FindClusters(testis.integrated, reduction = "umap", dims =c(1:5, 8, 16, 18, 20), resolution=2)
testis.integrated<-RunTSNE(testis.integrated,dims=c(1:5, 8, 16, 18, 20))
