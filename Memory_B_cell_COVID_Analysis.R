# This script contains a big part of the downstream analysis for the scRNAseq data in the MBC project.
# It is divided into 5 parts (see below).
# All parts use as input the preprocessed "Bmem.rds" object.
# Some chunks also use additional files as input. In these cases, I have added some comments which explain the use of those files.

# Loading packages
library(rgl)
library(ggforce)
library(Seurat)
library(scRepertoire)
library(tidyverse)
library(Biostrings)
library(patchwork)
library(RColorBrewer)
library(gprofiler2)
library(umap)
library(ComplexHeatmap)
library(paletteer)
library(SeuratData)
library(circlize)
library(cowplot)
library(dplyr)
library(data.table)
library(bruceR)
library(ggridges)
library(MetBrewer)
library(openxlsx)
library(fgsea)
library(msigdbr)
library(viridis)
library(pheatmap)
library(grid)
library(gridExtra)
library(ggpubr)
library(waffle)
library(rstatix)
library(monocle3)
library(dittoSeq)
library(airr)

# Clearing the environment
rm(list = ls())

# Creating useful functions
center.title <- function(){theme(plot.title = element_text(hjust = 0.5))} # For plotting
black.axis.text <- function(){theme(axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"))} # For plotting
myscaling <- function(x) {
  functionframe <- x
  i <- 0
  for (i in 1:length(functionframe)) {
    functionframe[i] <- scaler(functionframe[i])
  }
  return(functionframe)
  
} # Used for scaling in LIBRA score calculation
Highlightscatter <- function(x, clonesequence){
  df <- x@meta.data
  df2 <- which(df == clonesequence, arr.ind = TRUE)
  cells_to_highlight <- rownames(df2)
  x@meta.data$highlight <- "no"
  x@meta.data$highlight[rownames(df)%in%cells_to_highlight] <- "yes"
  FeatureScatter(x,  feature1 = "cor_wt_Spike", feature2 = "cor_RBD", pt.size = 1, group.by = "highlight")+
    ggtitle("Wt Spike and RBD scores")+xlab("wt-Spike scores")+ylab("RBD scores")+
    
    FeatureScatter(x,  feature1 = "cor_wt_Spike", feature2 = "cor_B.1.351", pt.size = 1, group.by = "highlight")+
    ggtitle("Wt Spike and Beta variant scores")+xlab("wt-Spike scores")+ylab("Beta variant scores")+
    
    FeatureScatter(x,  feature1 = "cor_wt_Spike", feature2 = "cor_B.1.617.2", pt.size = 1, group.by = "highlight")+
    ggtitle("Wt Spike and Delta variant scores")+xlab("wt-Spike scores")+ylab("Delta variant scores")+
    
    FeatureScatter(x,  feature1 = "cor_B.1.351", feature2 = "cor_B.1.617.2", pt.size = 1, group.by = "highlight")+
    ggtitle("Beta and Delta variant scores")+xlab("Beta scores")+ylab("Delta variant scores")
  
} # Highlight cells belonging to a clone of interest in the FeatureScatter plot
Highlightscatter_cells <- function(x, cellbarcodes){
  df <- x@meta.data
  cells_to_highlight <- rownames(df)[which(rownames(df) %in% cellbarcodes)]
  x@meta.data$highlight <- "no"
  x@meta.data$highlight[rownames(df)%in%cells_to_highlight] <- "yes"
  FeatureScatter(x,  feature1 = "cor_wt_Spike", feature2 = "cor_RBD", pt.size = 1, group.by = "highlight")+
    ggtitle("Wt Spike and RBD scores")+xlab("wt-Spike scores")+ylab("RBD scores")+
    
    FeatureScatter(x,  feature1 = "cor_wt_Spike", feature2 = "cor_B.1.351", pt.size = 1, group.by = "highlight")+
    ggtitle("Wt Spike and Beta variant scores")+xlab("wt-Spike scores")+ylab("Beta variant scores")+
    
    FeatureScatter(x,  feature1 = "cor_wt_Spike", feature2 = "cor_B.1.617.2", pt.size = 1, group.by = "highlight")+
    ggtitle("Wt Spike and Delta variant scores")+xlab("wt-Spike scores")+ylab("Delta variant scores")+
    
    FeatureScatter(x,  feature1 = "cor_B.1.351", feature2 = "cor_B.1.617.2", pt.size = 1, group.by = "highlight")+
    ggtitle("Beta and Delta variant scores")+xlab("Beta scores")+ylab("Delta variant scores")
  
} # Highlight chose cells in the FeatureScatter plot
Highlightscatter_restrictive <- function(x, clonesequence){
  df <- x@meta.data
  df2 <- which(df == clonesequence, arr.ind = TRUE)
  cells_to_highlight <- rownames(df2)
  x@meta.data$highlight <- "no"
  x@meta.data$highlight[rownames(df)%in%cells_to_highlight] <- "yes"
  x <- subset(x, subset=highlight=="yes")
  FeatureScatter(x,  feature1 = "cor_wt_Spike", feature2 = "cor_RBD", pt.size = 1, group.by = "highlight",jitter = F)+
    ggtitle("Wt Spike and RBD scores")+xlab("wt-Spike scores")+ylab("RBD scores")+xlim(0,1)+ylim(0,1)+
    
    FeatureScatter(x,  feature1 = "cor_wt_Spike", feature2 = "cor_B.1.351", pt.size = 1, group.by = "highlight",jitter = F)+
    ggtitle("Wt Spike and Beta variant scores")+xlab("wt-Spike scores")+ylab("Beta variant scores")+xlim(0,1)+ylim(0,1)+
    
    FeatureScatter(x,  feature1 = "cor_wt_Spike", feature2 = "cor_B.1.617.2", pt.size = 1, group.by = "highlight",jitter = F)+
    ggtitle("Wt Spike and Delta variant scores")+xlab("wt-Spike scores")+ylab("Delta variant scores")+xlim(0,1)+ylim(0,1)+
    
    FeatureScatter(x,  feature1 = "cor_B.1.351", feature2 = "cor_B.1.617.2", pt.size = 1, group.by = "highlight",jitter = F)+
    ggtitle("Beta and Delta variant scores")+xlab("Beta scores")+ylab("Delta variant scores")+xlim(0,1)+ylim(0,1)
  
} # Like Highlightscatter but only chosen cells are shown in the plot
asterics.function <- function(prop.test.result) {
  
  if(prop.test.result$p.value>0.05){
    return <- "ns"
  }
  if(prop.test.result$p.value<=0.05){
    retun <- "*"
  }
  if(prop.test.result$p.value<=0.01){
    return <- "**"
  }
  if(prop.test.result$p.value<=0.001){
    return <- "***"
  }
  if(prop.test.result$p.value<=0.0001){
    return <- "****"
  }
  return(return)
} # Used in Part 15 to plot significance asterics for waffle plot.
Stats <- function(x){
  print(paste0("The object has ", nrow(x@meta.data), " cells."))
  print(paste0("Of these, ",nrow(x@meta.data[x@meta.data$bait.positive=="yes",])," are Bait positive"))
  print(paste0("There are ",nrow(x@meta.data[x@meta.data$full.BCR.known=="yes",]), " cells with full BCR info"))
  print(paste0("There are ",nrow(x@meta.data[x@meta.data$common.clone =="yes",]), " cells annotated as common clone"))
  print("Sample group stats:")
  print(table(x@meta.data$sample_group))
  print("Which percentage of Spike wt binding BCRs are binders of Beta and Delta (per timepoint)")
  Early_binders <- x@meta.data[x@meta.data$Timepoint=="6 months after infection" &
                                 x@meta.data$cor_wt_Spike >0,]
  Late_binders <- x@meta.data[x@meta.data$Timepoint=="12 months after infection" &
                                x@meta.data$cor_wt_Spike >0 ,]
  
  #Beta binders
  print(paste0(100*(sum(Early_binders$B.1.351_positive=="yes")/nrow(Early_binders)), "% of early binders are Beta binders"))
  print(paste0(100*(sum(Late_binders$B.1.351_positive=="yes")/nrow(Late_binders)), "% of late binders are Beta binders"))
  
  # Delta binders
  print(paste0(100*(sum(Early_binders$B.1.617.2_positive=="yes")/nrow(Early_binders)), "% of early binders are Delta binders"))
  print(paste0(100*(sum(Late_binders$B.1.617.2_positive=="yes")/nrow(Late_binders)), "% of late binders are Delta binders"))
  
  # All binders
  print(paste0(100*(sum(Early_binders$all_positive =="yes")/nrow(Early_binders)), "% of early binders are universal binders"))
  print(paste0(100*(sum(Late_binders$all_positive=="yes")/nrow(Late_binders)), "% of late binders are universal binders"))
  
  remove(Early_binders,Late_binders)
  print("")
  print("Which percentage of RBD wt binding BCRs are binders of Beta and Delta (per timepoint)")
  Early_binders <- x@meta.data[x@meta.data$Timepoint=="6 months after infection" &
                                 x@meta.data$cor_RBD >0,]
  Late_binders <- x@meta.data[x@meta.data$Timepoint=="12 months after infection" &
                                x@meta.data$cor_RBD >0 ,]
  sixmo_binders <- x@meta.data[x@meta.data$sample_group=="6 months" &
                                 x@meta.data$cor_RBD>0 ,]
  tvelvemo_nonvax_binders <- x@meta.data[x@meta.data$sample_group=="12 months not vaccinated" &
                                           x@meta.data$cor_RBD>0 ,]
  tvelvemo_vax_binders <- x@meta.data[x@meta.data$sample_group=="12 months vaccinated" &
                                        x@meta.data$cor_RBD>0 ,]
  
  #Beta binders
  print("Beta binders:")
  print(paste0(100*(sum(Early_binders$B.1.351_positive=="yes")/nrow(Early_binders)), "% of early binders are Beta binders"))
  print(paste0(100*(sum(Late_binders$B.1.351_positive=="yes")/nrow(Late_binders)), "% of late binders are Beta binders"))
  print("")
  print(paste0(100*(sum(sixmo_binders$B.1.351_positive=="yes")/nrow(sixmo_binders)), "% of sixmo_binders are Beta binders"))
  print(paste0(100*(sum(tvelvemo_nonvax_binders$B.1.351_positive=="yes")/nrow(tvelvemo_nonvax_binders)), "% of tvelvemo_nonvax_binders are Beta binders"))
  print(paste0(100*(sum(tvelvemo_vax_binders$B.1.351_positive=="yes")/nrow(tvelvemo_vax_binders)), "% of tvelvemo_vax_binders are Beta binders"))
  print("")
  
  # Delta binders
  print("Delta binders:")
  print(paste0(100*(sum(Early_binders$B.1.617.2_positive=="yes")/nrow(Early_binders)), "% of early binders are Delta binders"))
  print(paste0(100*(sum(Late_binders$B.1.617.2_positive=="yes")/nrow(Late_binders)), "% of late binders are Delta binders"))
  print("")
  print(paste0(100*(sum(sixmo_binders$B.1.617.2_positive=="yes")/nrow(sixmo_binders)), "% of sixmo_binders are Delta binders"))
  print(paste0(100*(sum(tvelvemo_nonvax_binders$B.1.617.2_positive=="yes")/nrow(tvelvemo_nonvax_binders)), "% of tvelvemo_nonvax_binders are Delta binders"))
  print(paste0(100*(sum(tvelvemo_vax_binders$B.1.617.2_positive=="yes")/nrow(tvelvemo_vax_binders)), "% of tvelvemo_vax_binders are Delta binders"))
  print("")
  
  # All binders (among RBD binders)
  print("Universal binders:")
  print(paste0(100*(sum(Early_binders$all_positive =="yes")/nrow(Early_binders)), "% of early binders are universal binders"))
  print(paste0(100*(sum(Late_binders$all_positive=="yes")/nrow(Late_binders)), "% of late binders are universal binders"))
  print("")
  print(paste0(100*(sum(sixmo_binders$all_positive=="yes")/nrow(sixmo_binders)), "% of sixmo_binders are universal binders"))
  print(paste0(100*(sum(tvelvemo_nonvax_binders$all_positive=="yes")/nrow(tvelvemo_nonvax_binders)), "% of tvelvemo_nonvax_binders are universal binders"))
  print(paste0(100*(sum(tvelvemo_vax_binders$all_positive=="yes")/nrow(tvelvemo_vax_binders)), "% of tvelvemo_vax_binders are universal binders"))
  
  
  print("")
  print("Which percentage of full length Spike binders are RBD binders?")
  print(paste(100*nrow(Bmem@meta.data[Bmem@meta.data$spike_RBD_positive=="yes",])/nrow(Bmem@meta.data[Bmem@meta.data$cor_wt_Spike>0,]),"%"))
  
  
  
  remove(Early_binders,Late_binders)
}
`%notin%` <- Negate(`%in%`) # Creating a useful operator

######################################################################################################################
# Part 1: Exploration of the integrated dataset and plotting the results from some brief analyses.
######################################################################################################################
# Loading the data into the session.
setwd("~/NAS/Jan/Experiments and Data/210801 - Exp 035 - Profiling of COVID vaccinated patient samples/Memory B cell Fourth Experiment 220323/multi_exp035_4")  
Bmem <- readRDS("Bmem.rds")

# For comparing mutational counts of groups of cells to those of naive cells, I am adding the following objects to the session.
# They are generated as part of the preprocessing and are included in the zenodo dataset repository.
Bmem_fourth <- readRDS("Bmem_fourth.IgD.rds")
Naives <- (read.table("Naives.txt"))
Naives <- Naives$V1

# Fig.5.A.2
# Bait positive cells per timepoint
sixmo_positives <- rownames(Bmem@meta.data[Bmem@meta.data$bait.positive=="yes"&
                                             Bmem@meta.data$sample_group=="6 months",])
twelvemo_positives <- rownames(Bmem@meta.data[Bmem@meta.data$bait.positive=="yes"&
                                             Bmem@meta.data$sample_group=="12 months not vaccinated",])
twelvemo_vax_positives <- rownames(Bmem@meta.data[Bmem@meta.data$bait.positive=="yes"&
                                                Bmem@meta.data$sample_group=="12 months vaccinated",])
plot1 <- DimPlot(Bmem, reduction = 'wnn.umap', label = F, repel = TRUE, label.size = 2.5, cells.highlight = sixmo_positives,sizes.highlight = 0.5)+
  NoLegend() +NoAxes() + ggtitle("6 months")+center.title()
plot2 <- DimPlot(Bmem, reduction = 'wnn.umap', label = F, repel = TRUE, label.size = 2.5, cells.highlight = twelvemo_positives,sizes.highlight = 0.5)+
  NoLegend() +NoAxes() +ggtitle("12 months")+center.title()
plot3 <- DimPlot(Bmem, reduction = 'wnn.umap', label = F, repel = TRUE, label.size = 2.5, cells.highlight = twelvemo_vax_positives,sizes.highlight = 0.5)+
  NoLegend() +NoAxes() +ggtitle("12 months \n+ vaccinated")+center.title()

grid.arrange(plot1, plot2, plot3, ncol=3) # Export 3x12

remove(plot1,plot2,plot3, sixmo_positives, twelvemo_positives, twelvemo_vax_positives)

# Supp.Fig.5.B.2
DefaultAssay(Bmem) <- "Protein"
FeaturePlot(Bmem, features = c("FcRL5", "CD21", "CD27.1", "CXCR5.1", "CD71"), reduction = "wnn.umap",ncol = 3)& NoAxes() # Export 5x12


# Featurescatter on the integrated dataset. Save 6x9 or 5x7
# Supp.Fig.2.A
FeatureScatter(Bmem,  feature1 = "cor_wt_Spike", feature2 = "cor_RBD", pt.size = 0.6, group.by = "orig.ident", cols = "black")+
  ggtitle("")+xlab("Wuhan-Hu-1 Spike LIBRA scores")+ylab("Wuhan-Hu-1 RBD LIBRA scores")+
  theme_classic()+theme(legend.position = "none",text = element_text(size = 12))+center.title()+black.axis.text()

# One density plot for Supp.Fig.2.A
d_cor_wt_Spike <- density(Bmem@meta.data$cor_wt_Spike)
plot(d_cor_wt_Spike, main="", axes=FALSE, frame.plot=F, xaxt='n', ann=FALSE)
grid(nx = NULL, ny = NULL,
     lty = 1,      # Grid line type
     col = "gray", # Grid line color
     lwd = 1)      # Grid line width
polygon(d_cor_wt_Spike, col="black")
axis(side=1, labels=T, at=c(0,0.25,0.5,0.75,1)) # Export 3x8

# Second density plot for Supp.Fig.2.A
d_cor_RBD <- density(Bmem@meta.data$cor_RBD)
plot(d_cor_RBD, main="",axes=F,frame.plot=T,xaxt="n",ann=F,xlim=c(1,-0.07))
grid(nx = NULL, ny = NULL,
     lty = 1,      # Grid line type
     col = "gray", # Grid line color
     lwd = 1)      # Grid line width
polygon(d_cor_RBD, col="black")
Axis(side=1, labels=T, at=c(0,0.25,0.5,0.75,1)) # Export 3x5.5

remove(d_cor_RBD,d_cor_wt_Spike)

# Supp.Fig.2.C
FeatureScatter(Bmem,  feature1 = "cor_wt_Spike", feature2 = "cor_RBD", pt.size = 0.6, group.by = "Vaccinated",cols = c("chocolate3","deepskyblue3"))+
  ggtitle("Wt Spike and RBD binding scores")+xlab("Wuhan-Hu-1 Spike binding scores")+ylab("Wuhan-Hu-1 RBD binding scores")+labs(color = "Vaccinated")+
  
  FeatureScatter(Bmem,  feature1 = "cor_wt_Spike", feature2 = "cor_B.1.351", pt.size = 0.6, group.by = "Vaccinated",cols = c("chocolate3","deepskyblue3"))+
  ggtitle("Wt Spike and Beta variant scores")+xlab("wt-Spike scores")+ylab("Beta variant scores")+labs(color = "Vaccinated")+
  
  FeatureScatter(Bmem,  feature1 = "cor_wt_Spike", feature2 = "cor_B.1.617.2", pt.size = 0.6, group.by = "Vaccinated",cols = c("chocolate3","deepskyblue3"))+
  ggtitle("Wt Spike and Delta variant scores")+xlab("wt-Spike scores")+ylab("Delta variant scores")+labs(color = "Vaccinated")+
  
  FeatureScatter(Bmem,  feature1 = "cor_B.1.351", feature2 = "cor_B.1.617.2", pt.size = 0.6, group.by = "Vaccinated",cols = c("chocolate3","deepskyblue3"))+
  ggtitle("Beta and Delta variant scores")+xlab("Beta variant scores")+ylab("Delta variant scores")+labs(color = "Vaccinated") # Export 7x10

# Supp.Fig.5.C
# Percentage plot showing the percent of each patient per WNN cluster
Patients.df <- as.data.frame(with(Bmem@meta.data, table(Patient, seurat_clusters)))
Patients.df <- Patients.df %>% group_by(seurat_clusters) %>% add_tally(Freq, name = "Cells per cluster")
Patients.df$Patient <- as.character(Patients.df$Patient)
Patients <- as.character(unique(Patients.df$Patient))
charactervector <- c("A","B","C","D","E","F","G","H","I")
for (i in 1:length(Patients)) {
  Patients.df$Patient[Patients.df$Patient==Patients[i]] <- charactervector[i]
}
Patients.df$pat.perc.clust <- 100*round(Patients.df$Freq/Patients.df$`Cells per cluster`,3)
Patients.df$seurat_clusters <- as.factor(Patients.df$seurat_clusters)
Patients.df$Patient <- as.character(Patients.df$Patient)
ggplot(Patients.df, 
       aes(fill=factor(Patient), y=seurat_clusters, x=pat.perc.clust)) + 
  geom_bar(position="fill", stat="identity",colour="black")+theme_classic()+ theme(text = element_text(size = 15))+
  labs(x = "Percentages of patients per cluster", y="WNN clusters")+scale_x_continuous(labels=scales::percent)+labs(fill = "Patient")+
  scale_fill_manual(values=met.brewer("Austria", 9))+ theme(text = element_text(size = 17))+ggtitle("Patients per WNN cluster \n All cells")+
  center.title()+black.axis.text() # Export 5x8

remove(Patients.df, Patients,charactervector)

# Percentage plot showing the percent of each patient per WNN cluster - Bait positive cells only!
Bmem.subset <- subset(Bmem, subset = bait.positive =="yes")
Patients.df <- as.data.frame(with(Bmem.subset@meta.data, table(Patient, seurat_clusters)))
Patients.df$Patient <- as.character(Patients.df$Patient)
Patients <- as.character(unique(Patients.df$Patient))
for (i in 1:length(Patients)) {
  Patients.df$Patient[Patients.df$Patient==Patients[i]] <- i
}
Patients.df <- Patients.df %>% group_by(seurat_clusters) %>% add_tally(Freq, name = "Cells per cluster")
Patients.df$pat.perc.clust <- 100*round(Patients.df$Freq/Patients.df$`Cells per cluster`,3)
Patients.df$seurat_clusters <- as.factor(Patients.df$seurat_clusters)

ggplot(Patients.df, 
       aes(fill=factor(Patient), y=seurat_clusters, x=pat.perc.clust)) + 
  geom_bar(position="fill", stat="identity",colour="black")+theme_classic()+ theme(text = element_text(size = 15))+
  labs(x = "Percentages of patients per cluster", y="WNN clusters")+scale_x_continuous(labels=scales::percent)+labs(fill = "Patient")+
  scale_fill_manual(values=met.brewer("Austria", 9))+ theme(text = element_text(size = 17))+ggtitle("Patients per cluster \n COVID specific cells")+
  center.title()+black.axis.text()

remove(Patients.df,Bmem.subset,Patients)


# Percentage plot showing the percent of each patient per cell subset
Patients.df <- as.data.frame(with(Bmem@meta.data, table(Patient, named.clusters)))
Patients.df <- Patients.df %>% group_by(named.clusters) %>% add_tally(Freq, name = "Cells per Subset")
Patients.df$pat.perc.clust <- 100*round(Patients.df$Freq/Patients.df$`Cells per Subset`,3)
Patients.df$named.clusters <- factor(Patients.df$named.clusters, levels=rev(c("Unswitched", "CD27low RM","CD27high RM","Activated","Atypical")))

ggplot(Patients.df, 
       aes(fill=factor(Patient), y=named.clusters, x=pat.perc.clust)) + 
  geom_bar(position="fill", stat="identity",colour="black")+theme_classic()+ theme(text = element_text(size = 15))+
  labs(x = "Percentages of patients per Subset", y="Subsets")+scale_x_continuous(labels=scales::percent)+labs(fill = "Patient")+
  scale_fill_manual(values=met.brewer("Austria", 9))+ theme(text = element_text(size = 17))+ggtitle("Patients per Subset \n All cells")+
  center.title()

remove(Patients.df)


# Percentage plot showing the percent of each patient per  cell subset - Bait positive cells only!
Bmem.subset <- subset(Bmem, subset = bait.positive =="yes")
Patients.df <- as.data.frame(with(Bmem.subset@meta.data, table(Patient, named.clusters)))
Patients.df <- Patients.df %>% group_by(named.clusters) %>% add_tally(Freq, name = "Cells per Subset")
Patients.df$pat.perc.clust <- 100*round(Patients.df$Freq/Patients.df$`Cells per Subset`,3)
Patients.df$named.clusters <- factor(Patients.df$named.clusters, levels=rev(c("Unswitched", "CD27low RM","CD27high RM","Activated","Atypical")))

ggplot(Patients.df, 
       aes(fill=factor(Patient), y=named.clusters, x=pat.perc.clust)) + 
  geom_bar(position="fill", stat="identity",colour="black")+theme_classic()+ theme(text = element_text(size = 15))+
  labs(x = "Percentages of patients per Subset", y="Subsets")+scale_x_continuous(labels=scales::percent)+labs(fill = "Patient")+
  scale_fill_manual(values=met.brewer("Austria", 9))+ theme(text = element_text(size = 17))+ggtitle("Patients per Subset \n COVID spec. cells")+
  center.title()

remove(Patients.df,Bmem.subset)

# Supp.Fig.2.B
# Pie chart showing the percentage of RBD binders among spike binders:
Percentage.of.RBD.binders <- round(100*nrow(Bmem@meta.data[Bmem@meta.data$spike_RBD_positive=="yes",])/nrow(Bmem@meta.data[Bmem@meta.data$cor_wt_Spike>0,]),1)
Percentage.of.non.RBD.binders <- 100-Percentage.of.RBD.binders
data <- data.frame(
  Epitope_Specificity=c("RBD Binders","Non RBD Binders"),
  value=c(Percentage.of.RBD.binders,Percentage.of.non.RBD.binders))
data$Epitope_Specificity <- factor(data$Epitope_Specificity, levels = c("Non RBD Binders","RBD Binders"))
colnames(data) <- c("Specificity","value")

ggplot(data, aes(x="", y=value, fill=Specificity)) +
  geom_bar(stat="identity", width=1,colour="black") +
  coord_polar("y", start=0)+theme_void()+
  scale_fill_manual(values=met.brewer("Austria", 2))+ggtitle("% RBD binders of full length Spike binders")+
  geom_text(x=1.1, y=15, label=paste0(Percentage.of.RBD.binders,"%"), colour="white",size=6)+center.title() #Export 5x6

remove(Percentage.of.non.RBD.binders,Percentage.of.RBD.binders,data)

# Supp.Fig.5.D
# Pie charts showing the percentage of antigen specific cells per cell subset
df <- Bmem@meta.data
df <- df[,c("named.clusters","bait.positive")]
df <- as.data.frame(table(df$named.clusters,df$bait.positive))
colnames(df) <- c("Var1","COVID.specific","Freq")
Pie_Unswitched <- df[df$Var1=="Unswitched",]
Pie_27low <- df[df$Var1=="CD27low RM",]
Pie_27high <- df[df$Var1=="CD27high RM",]
Pie_Activated <- df[df$Var1=="Activated",]
Pie_Atypical <- df[df$Var1=="Atypical",]

p1 <- ggplot(Pie_Unswitched, aes(x="", y=Freq, fill=COVID.specific)) +
  geom_bar(stat="identity", width=1,colour="black") +
  coord_polar("y", start=0)+theme_void()+
  scale_fill_manual(values=met.brewer("Egypt", 2))+ggtitle("Unswitched cells")+
  geom_text(x=1.3, y=110, label=paste0(round(100*(Pie_Unswitched$Freq[2]/(Pie_Unswitched$Freq[2]+Pie_Unswitched$Freq[1])),2),"%"), colour="white",size=6)+center.title()+ theme(legend.position = "none")

p2 <- ggplot(Pie_27low, aes(x="", y=Freq, fill=COVID.specific)) +
  geom_bar(stat="identity", width=1,colour="black") +
  coord_polar("y", start=0)+theme_void()+
  scale_fill_manual(values=met.brewer("Egypt", 2))+ggtitle("CD27low RM cells")+
  geom_text(x=1.1, y=200, label=paste0(round(100*(Pie_27low$Freq[2]/(Pie_27low$Freq[2]+Pie_27low$Freq[1])),2),"%"), colour="white",size=6)+center.title()+ theme(legend.position = "none")

p3 <- ggplot(Pie_27high, aes(x="", y=Freq, fill=COVID.specific)) +
  geom_bar(stat="identity", width=1,colour="black") +
  coord_polar("y", start=0)+theme_void()+
  scale_fill_manual(values=met.brewer("Egypt", 2))+ggtitle("CD27high RM cells")+
  geom_text(x=1.1, y=600, label=paste0(round(100*(Pie_27high$Freq[2]/(Pie_27high$Freq[2]+Pie_27high$Freq[1])),2),"%"), colour="white",size=6)+center.title()+ theme(legend.position = "none")

p4 <- ggplot(Pie_Activated, aes(x="", y=Freq, fill=COVID.specific)) +
  geom_bar(stat="identity", width=1,colour="black") +
  coord_polar("y", start=0)+theme_void()+
  scale_fill_manual(values=met.brewer("Egypt", 2))+ggtitle("Activated cells")+
  geom_text(x=1.1, y=300, label=paste0(round(100*(Pie_Activated$Freq[2]/(Pie_Activated$Freq[2]+Pie_Activated$Freq[1])),2),"%"), colour="white",size=6)+center.title()+ theme(legend.position = "none")

p5 <- ggplot(Pie_Atypical, aes(x="", y=Freq, fill=COVID.specific)) +
  geom_bar(stat="identity", width=1,colour="black") +
  coord_polar("y", start=0)+theme_void()+labs(fill="SARS-CoV-2 specific")+
  scale_fill_manual(values=met.brewer("Egypt", 2))+ggtitle("Atypical cells")+
  geom_text(x=1, y=300, label=paste0(round(100*(Pie_Atypical$Freq[2]/(Pie_Atypical$Freq[2]+Pie_Atypical$Freq[1])),2),"%"), colour="white",size=6)+center.title()


grid.arrange(p1,p2,p3,p4,p5, nrow = 3,top="Antigen binding cells among cell subsets") # Export 8x8
remove(p1,p2,p3,p4,p5,Pie_27high,Pie_27low,Pie_Atypical,Pie_Activated,Pie_Unswitched)

# Percentage plot showing the number of cells in each subset in different sample groups. - Only bait positive cells # Export 5x9 - Fig.5.B
phenotype.contributions.df <- as.data.frame(with(Bmem@meta.data[Bmem@meta.data$bait.positive=="yes",], table(named.clusters, sample_group)))
phenotype.contributions.df <- phenotype.contributions.df %>% group_by(sample_group) %>% add_tally(Freq,name = "total.cells.persamplegroup") %>% ungroup()
phenotype.contributions.df <- phenotype.contributions.df %>% mutate(percentages= 100*(Freq/total.cells.persamplegroup))
phenotype.contributions.df$percentages <- round(phenotype.contributions.df$percentages,digits = 1)
phenotype.contributions.df$percentages <- as.numeric(phenotype.contributions.df$percentages)

phenotype.contributions.df$sample_group <- factor(phenotype.contributions.df$sample_group, levels = rev(c("6 months","12 months not vaccinated","12 months vaccinated")))

ggplot(phenotype.contributions.df, 
       aes(fill=factor(named.clusters, levels=rev(c("Unswitched",
                                                    "CD27low RM","CD27high RM",
                                                    "Activated","Atypical"))), y=sample_group, x=percentages)) + 
  geom_bar(position="fill", stat="identity",colour="black")+theme_classic()+ theme(text = element_text(size = 15),plot.title = element_text(hjust = 0.5),
                                                                                   plot.subtitle = element_text(hjust = 0.5))+
  labs(x = "", y="")+scale_x_continuous(labels=scales::percent)+labs(fill = "Cell subset")+scale_fill_manual(values=rev(c("#DBC35E","#228833","#4477AA","#EE6677","#AA3377")))+ 
  theme(text = element_text(size = 17))+ ggtitle("SARS-CoV-2+ MBC")+
  black.axis.text()

remove(phenotype.contributions.df) # Export 4x8


# How similar do cells of the same clone (same CTaa) behave? Supp.Fig.2.E
Bmem.intermediate <- subset(Bmem, subset = full.BCR.known=="yes")
df <- Bmem.intermediate@meta.data
df <- df[,c(1,26,37)]
df <- df %>% group_by(CTaa) %>% add_tally(name="number.of.CTaas") %>% ungroup()
df2 <- df[df$bait.positive =="yes",2]
df2 <- df2 %>% group_by(CTaa) %>% add_tally(name="number.of.positive.CTaas") %>% ungroup()
df2 <- unique(df2)
df3 <- merge(df,df2,by = "CTaa",all.x = T)
df3$number.of.positive.CTaas[is.na(df3$number.of.positive.CTaas)] <- 0
df3 <- merge(df3,aggregate(number.of.positive.CTaas~CTaa,df3,FUN=max),all.x=TRUE,by="CTaa")
df3 <- df3[,1:5]
colnames(df3) <- c("CTaa","Row.names","bait.positive","number.of.CTaas","number.of.positive.CTaas")
df3$percentage.positive <- 100*round(df3$number.of.positive.CTaas/df3$number.of.CTaas,2)
df3$percent.with.values <- paste(100*round(df3$number.of.positive.CTaas/df3$number.of.CTaas,2)," (",df3$number.of.positive.CTaas,"of",df3$number.of.CTaas,")")

# Subset to only expanded clones
df4 <- df3[df3$number.of.CTaas>1, ]
df4 <- df4[,c(1,6)]
df4 <- unique(df4)
n <- nrow(df4)
df4$percentage.positive <- as.factor(df4$percentage.positive)
df4 <- df4 %>% group_by(percentage.positive) %>% add_tally(name="number.per.percentage") %>% ungroup()
df4.sub <- unique(df4[,2:3])
ggplot(df4.sub, aes(x=percentage.positive, y=number.per.percentage, fill=percentage.positive))+geom_bar(stat = "identity",colour="black")+
  labs(x="Percent of clone classified bait positive",y="Number of expanded clones")+theme_classic()+
  scale_y_continuous(expand = c(0,0),limits = c(0,80))+black.axis.text()+ labs(fill = "Percent of clone \nclassified bait positive")+
  scale_fill_manual(values=met.brewer("Egypt", 5))

# Pie chart
df4$percent.of.total <- 100*round(df4$number.per.percentage/nrow(df4),2)
df4 <- df4[,c(2,4)]
df4 <- unique(df4)
df4
df4$percentage.positive <- c("100%", "0%","50%","67%","33%")
df4$percentage.positive <- factor(df4$percentage.positive, levels = c("33%", "67%", "50%", "0%","100%"))
df5 <- data.frame(Specificity=c("Matching specificity","Mixed specificity"),Percent=c(90,10))
ggplot(df5, aes(x="", y=Percent, fill=Specificity)) +
  geom_bar(stat="identity", width=1,colour="black") +
  coord_polar("y", start=0)+theme_void()+labs(fill = "Specificity")+
  scale_fill_manual(values=met.brewer("Egypt", 2))+ggtitle(paste0("Percentage of clones with strictly matching specificity \n (Expanded clones only)",
                                                                  "\n n=",n))+ center.title()+
  geom_text(x=1, y=34, label=paste0(df5[1,2],"%"), colour="white",size=7) # Export 5x6

remove(Bmem.intermediate, df, df2, df3, df4,df4.sub,df5,n)


# Fig.2.H
# Binding breadth plot - Using the two-proportions z-test to compare proportions.
# Use this if you want to look at cells that bind the Wuhan Hu 1 full length spike antigen.
Bmem.subset <- subset(Bmem, subset = cor_wt_Spike > 0)

# Use this if you want to subset to patients who are not vaccinated - 3 patients only then.
#Bmem.subset <- subset(Bmem.subset, subset = Patient=="PFCL1-196878-1691" | 
#                        Patient == "PFCL1-LIM-674950" |
#                        Patient == "PFCL1-UST-190762")

Bmem.subset@meta.data$Beta_and_Delta_positive <- "no"
Bmem.subset@meta.data$Beta_and_Delta_positive[Bmem.subset@meta.data$B.1.351_positive=="yes" & 
                                                Bmem.subset@meta.data$B.1.617.2_positive=="yes"] <- "yes"

# The following two short paragraphs only matter for cells from vaccinated samples. You can either include all of them, or only those from specific times after vac.
# If you only want to look at short after vaccination, use these lines
#Bmem.subset <- subset(Bmem.subset, subset = time_after_vac != "85d" &
#                        time_after_vac != "87d" & 
#                        time_after_vac != "108d")

# If you only want to look at long after vaccination (vaccination short after infection), use these lines
#Bmem.subset <- subset(Bmem.subset, subset = time_after_vac != "23d" &
#                      time_after_vac != "15d")


# Use this if you want to look at cells that bind the Wuhan Hu 1 full length spike and the RBD antigen.
#Bmem.subset <- subset(Bmem, subset = spike_RBD_positive =="yes")

# Making a data frame out of the subsetted observations. - Specify here which variant.
df <- data.frame(table(Bmem.subset@meta.data$sample_group,Bmem.subset@meta.data$B.1.351_positive))
df2 <- data.frame(table(Bmem.subset@meta.data$sample_group,Bmem.subset@meta.data$B.1.617.2_positive))
df3 <- data.frame(table(Bmem.subset@meta.data$sample_group,Bmem.subset@meta.data$Beta_and_Delta_positive))
colnames(df) <- c("Timepoint","Binding","Freq")
colnames(df2) <- c("Timepoint","Binding","Freq")
colnames(df3) <- c("Timepoint","Binding","Freq")

# For six months vs. 12 months not vaccinated - Beta
test1.Beta <- prop.test(x = c(df[df$Timepoint=="6 months"& df$Binding=="yes",3], 
                df[df$Timepoint=="12 months not vaccinated"& df$Binding=="yes",3]), 
                n = c(df[df$Timepoint=="6 months"& df$Binding=="yes",3]+
                        df[df$Timepoint=="6 months"& df$Binding=="no",3], 
                      df[df$Timepoint=="12 months not vaccinated"& df$Binding=="yes",3]+
                        df[df$Timepoint=="12 months not vaccinated"& df$Binding=="no",3]))
test1.Beta.stats <- asterics.function(test1.Beta)

# For 6 months vs. 12 months vaccinated - Beta
test2.Beta <- prop.test(x = c(df[df$Timepoint=="6 months"& df$Binding=="yes",3], 
                df[df$Timepoint=="12 months vaccinated"& df$Binding=="yes",3]), 
          n = c(df[df$Timepoint=="6 months"& df$Binding=="yes",3]+
                  df[df$Timepoint=="6 months"& df$Binding=="no",3], 
                df[df$Timepoint=="12 months vaccinated"& df$Binding=="yes",3]+
                  df[df$Timepoint=="12 months vaccinated"& df$Binding=="no",3]))
test2.Beta.stats <- asterics.function(test2.Beta)

# For 12 months not vaccinated vs. 12 months vaccinated - Beta
test3.Beta <- prop.test(x = c(df[df$Timepoint=="12 months not vaccinated"& df$Binding=="yes",3], 
                df[df$Timepoint=="12 months vaccinated"& df$Binding=="yes",3]), 
          n = c(df[df$Timepoint=="12 months not vaccinated"& df$Binding=="yes",3]+
                  df[df$Timepoint=="12 months not vaccinated"& df$Binding=="no",3], 
                df[df$Timepoint=="12 months vaccinated"& df$Binding=="yes",3]+
                  df[df$Timepoint=="12 months vaccinated"& df$Binding=="no",3]))
test3.Beta.stats <- asterics.function(test3.Beta)

# For six months vs. 12 months not vaccinated - Delta
test1.Delta <- prop.test(x = c(df2[df2$Timepoint=="6 months"& df2$Binding=="yes",3], 
                         df2[df2$Timepoint=="12 months not vaccinated"& df2$Binding=="yes",3]), 
                   n = c(df2[df2$Timepoint=="6 months"& df2$Binding=="yes",3]+
                           df2[df2$Timepoint=="6 months"& df2$Binding=="no",3], 
                         df2[df2$Timepoint=="12 months not vaccinated"& df2$Binding=="yes",3]+
                           df2[df2$Timepoint=="12 months not vaccinated"& df2$Binding=="no",3]))
test1.Delta.stats <- asterics.function(test1.Delta)

# For 6 months vs. 12 months vaccinated - Delta
test2.Delta <- prop.test(x = c(df2[df2$Timepoint=="6 months"& df2$Binding=="yes",3], 
                         df2[df2$Timepoint=="12 months vaccinated"& df2$Binding=="yes",3]), 
                   n = c(df2[df2$Timepoint=="6 months"& df2$Binding=="yes",3]+
                           df2[df2$Timepoint=="6 months"& df2$Binding=="no",3], 
                         df2[df2$Timepoint=="12 months vaccinated"& df2$Binding=="yes",3]+
                           df2[df2$Timepoint=="12 months vaccinated"& df2$Binding=="no",3]))
test2.Delta.stats <- asterics.function(test2.Delta)

# For 12 months not vaccinated vs. 12 months vaccinated - Delta
test3.Delta <- prop.test(x = c(df2[df2$Timepoint=="12 months not vaccinated"& df2$Binding=="yes",3], 
                         df2[df2$Timepoint=="12 months vaccinated"& df2$Binding=="yes",3]), 
                   n = c(df2[df2$Timepoint=="12 months not vaccinated"& df2$Binding=="yes",3]+
                           df2[df2$Timepoint=="12 months not vaccinated"& df2$Binding=="no",3], 
                         df2[df2$Timepoint=="12 months vaccinated"& df2$Binding=="yes",3]+
                           df2[df2$Timepoint=="12 months vaccinated"& df2$Binding=="no",3]))
test3.Delta.stats <- asterics.function(test3.Delta)

# For six months vs. 12 months not vaccinated - Beta and Delta
test1.Both <- prop.test(x = c(df3[df3$Timepoint=="6 months"& df3$Binding=="yes",3], 
                         df3[df3$Timepoint=="12 months not vaccinated"& df3$Binding=="yes",3]), 
                   n = c(df3[df3$Timepoint=="6 months"& df3$Binding=="yes",3]+
                           df3[df3$Timepoint=="6 months"& df3$Binding=="no",3], 
                         df3[df3$Timepoint=="12 months not vaccinated"& df3$Binding=="yes",3]+
                           df3[df3$Timepoint=="12 months not vaccinated"& df3$Binding=="no",3]))
test1.Both.stats <- asterics.function(test1.Both)

# For 6 months vs. 12 months vaccinated - Beta and Delta
test2.Both <- prop.test(x = c(df3[df3$Timepoint=="6 months"& df3$Binding=="yes",3], 
                         df3[df3$Timepoint=="12 months vaccinated"& df3$Binding=="yes",3]), 
                   n = c(df3[df3$Timepoint=="6 months"& df3$Binding=="yes",3]+
                           df3[df3$Timepoint=="6 months"& df3$Binding=="no",3], 
                         df3[df3$Timepoint=="12 months vaccinated"& df3$Binding=="yes",3]+
                           df3[df3$Timepoint=="12 months vaccinated"& df3$Binding=="no",3]))
test2.Both.stats <- asterics.function(test2.Both)

# For 12 months not vaccinated vs. 12 months vaccinated - Beta and Delta
test3.Both <- prop.test(x = c(df3[df3$Timepoint=="12 months not vaccinated"& df3$Binding=="yes",3], 
                         df3[df3$Timepoint=="12 months vaccinated"& df3$Binding=="yes",3]), 
                   n = c(df3[df3$Timepoint=="12 months not vaccinated"& df3$Binding=="yes",3]+
                           df3[df3$Timepoint=="12 months not vaccinated"& df3$Binding=="no",3], 
                         df3[df3$Timepoint=="12 months vaccinated"& df3$Binding=="yes",3]+
                           df3[df3$Timepoint=="12 months vaccinated"& df3$Binding=="no",3]))
test3.Both.stats <- asterics.function(test3.Both)

# Making a waffle chart for showing increased binding breath. 
# Only Wuhan Hu 1 positive cells will be regarded.
# Bait binding will be shown in 4 categories:
# 1) Cells which only bind Wuhan Hu 1 spike
# 2) Beta binders
# 3) Delta binders
# 4) All binders
df <- Bmem.subset@meta.data[,c("cor_wt_Spike","B.1.351_positive","B.1.617.2_positive","all_positive","sample_group")]
df_6mo <- df[df$sample_group=="6 months",]
df_12mo <- df[df$sample_group=="12 months not vaccinated",]
df_12mo_vax <- df[df$sample_group=="12 months vaccinated",]

# For 6 months time point
No_variant <- round(100*(nrow(df_6mo[df_6mo$B.1.351_positive=="no"& df_6mo$B.1.617.2_positive=="no",])/nrow(df_6mo)),0)
Beta <- round(100*(nrow(df_6mo[df_6mo$B.1.351_positive=="yes"& df_6mo$B.1.617.2_positive=="no",])/nrow(df_6mo)),0)
Delta <- round(100*(nrow(df_6mo[df_6mo$B.1.351_positive=="no"& df_6mo$B.1.617.2_positive=="yes",])/nrow(df_6mo)),0)
Both <- round(100*(nrow(df_6mo[df_6mo$B.1.351_positive=="yes"& df_6mo$B.1.617.2_positive=="yes",])/nrow(df_6mo)),0)

# For 12 months time point
No_variant12 <- round(100*(nrow(df_12mo[df_12mo$B.1.351_positive=="no"& df_12mo$B.1.617.2_positive=="no",])/nrow(df_12mo)),0)
Beta12 <- round(100*(nrow(df_12mo[df_12mo$B.1.351_positive=="yes"& df_12mo$B.1.617.2_positive=="no",])/nrow(df_12mo)),0)
Delta12 <- round(100*(nrow(df_12mo[df_12mo$B.1.351_positive=="no"& df_12mo$B.1.617.2_positive=="yes",])/nrow(df_12mo)),0)
Both12 <- round(100*(nrow(df_12mo[df_12mo$B.1.351_positive=="yes"& df_12mo$B.1.617.2_positive=="yes",])/nrow(df_12mo)),0)-1

# For 12 months vax time point
No_variant12vax <- round(100*(nrow(df_12mo_vax[df_12mo_vax$B.1.351_positive=="no"& df_12mo_vax$B.1.617.2_positive=="no",])/nrow(df_12mo_vax)),0)
Beta12vax <- round(100*(nrow(df_12mo_vax[df_12mo_vax$B.1.351_positive=="yes"& df_12mo_vax$B.1.617.2_positive=="no",])/nrow(df_12mo_vax)),0)
Delta12vax <- round(100*(nrow(df_12mo_vax[df_12mo_vax$B.1.351_positive=="no"& df_12mo_vax$B.1.617.2_positive=="yes",])/nrow(df_12mo_vax)),0)
Both12vax <- round(100*(nrow(df_12mo_vax[df_12mo_vax$B.1.351_positive=="yes"& df_12mo_vax$B.1.617.2_positive=="yes",])/nrow(df_12mo_vax)),0)

# Plotting 
p1 <- waffle(c("No Variant"=No_variant,"Only Beta"=Beta,"Only Delta"=Delta,"Beta + Delta"=Both),rows = 10,flip = F,
       reverse = F,colors = met.brewer("Hokusai3", 4),title = "6 Months",legend_pos = "none",size=1)+center.title()
p2 <- waffle(c("No Variant"=No_variant12,"Only Beta"=Beta12,"Only Delta"=Delta12,"Beta + Delta"=Both12),rows = 10,flip = F,
             reverse = F,colors = met.brewer("Hokusai3", 4),title = "12 Months",size=1)+center.title()+guides(fill=guide_legend(title="Binding specificity"))
p2.1 <- waffle(c("No Variant"=No_variant12,"Only Beta"=Beta12,"Only Delta"=Delta12,"Beta + Delta"=Both12),rows = 10,flip = F,
             reverse = F,colors = met.brewer("Hokusai3", 4),title = "12 Months",legend_pos = "none",size=1)+center.title()
p3 <- waffle(c("No Variant"=No_variant12vax,"Only Beta"=Beta12vax,"Only Delta"=Delta12vax,"Beta + Delta"=Both12vax),rows = 10,flip = F,
             reverse = F,colors = met.brewer("Hokusai3", 4),title = "12 Months vac.",size=1)+center.title()+guides(fill=guide_legend(title="Binding specificity"))

p1+p2#Export 4x8
p1+p2.1+p3

# Custom significance plot
plot.new()
plot.window(xlim=c(0.5,2), ylim=c(0,11))
segments(0.5,8,2,8, lwd=2, col="black")
options(digits = 5)
text(0.65,8.5,"Antigen")+text(1,8.5,"6 months")+text(1.4,8.5,"12 months")+text(1.8,8.5,"p val")
text(0.65,6.5,"Beta")+text(0.65,4.5,"Delta")+text(0.65,2.5,"Both")
text(1,6.5,paste0(Beta+Both,"%"))+text(1,4.5,paste0(Delta+Both,"%"))+text(1,2.5,paste0(Both,"%"))
text(1.4,6.5,paste0(Beta12+Both12,"%"))+text(1.4,4.5,paste0(Delta12+Both12,"%"))+text(1.4,2.5,paste0(Both12,"%"))
text(1.8,6.5,test1.Beta.stats)+text(1.8,4.5,test1.Delta.stats)+text(1.8,2.5,test1.Both.stats)
text(1.25,11,"Antigen Binding Breadth")
stats.plot <- recordPlot()
# Export 4x5

remove(Bmem.subset, df, df_12mo,df_12mo_vax,df_6mo,df2,df3,p1,p2,p2.1,p3,stats.plot,test1.Beta,test1.Delta,test1.Both,test2.Beta,test2.Delta,
       test2.Both,test3.Beta,test3.Delta,test3.Both,Beta,Beta12,Beta12vax,Both,Both12,Both12vax,Delta,Delta12,Delta12vax,No_variant,No_variant12,
       No_variant12vax,test1.Beta.stats,test1.Delta.stats,test1.Both.stats,test2.Beta.stats,test2.Delta.stats,test2.Both.stats,test3.Beta.stats,
       test3.Delta.stats,test3.Both.stats)


# Make a heatmap out of the chain pairing
# Specify here in the first row of the section which cells should be included!
df <- Bmem@meta.data[Bmem@meta.data$full.BCR.known=="yes" ,c("CTgene","Full.Row.names")]

df$HV <- sub("\\..*", "", df$CTgene)
df$LV <- sub(".*_", "", df$CTgene)
df$LV <- sub("\\..*", "", df$LV)

df <- df[,3:4]
df <- table(df) %>%  as.data.frame()
df_relative <- df
df_relative$Freq <- 100*(df_relative$Freq/sum(df_relative$Freq))

ggplot(df_relative,aes(x=HV, y=LV, fill=Freq)) + geom_tile()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_fill_gradient(low = "black", high = "yellow")+labs(fill = "Frequency (%)")+
  ggtitle("BCR Heavy and light chain pairings among all cells")+center.title() # Export 8x8 or 10x10



######################################################################################################################
# Part 2: Differential gene expression analysis between timepoints and sample groups
######################################################################################################################
# Should we kick out cells from Dataset 2 for transcriptional profiling?
Investigation.set <- subset(Bmem, subset = Timepoint=="12 months after infection" &
                              Vaccinated=="no" & seurat_clusters=="1")
features <- c("BACH2","CXCR4","CD19","CD55","RHOB","S100A6","CD83","KLF2","NFKBIA")
DefaultAssay(Investigation.set) <- "RNA"
VlnPlot(Investigation.set, features = features, split.by = "Dataset")

df <- as.data.frame.matrix(with(Bmem@meta.data, table(seurat_clusters, Dataset)))
df$sum <- rowSums(df)
df$Dataset.2.proportion <- df$Second/df$sum

# Differential gene expression between positive cells 6mo vs 12 mo not vaccinated and cells with a specific phenotype only
Bait.positives <- subset(Bmem, subset = bait.positive =="yes" & Dataset!= "Second" & named.clusters=="CD27high RM"|
                           bait.positive =="yes" & Dataset!= "Second" & named.clusters=="CD27low RM")
Bait.positives <- subset(Bait.positives, subset = Vaccinated =="no")

Bait.positives@meta.data$Timepoint[Bait.positives@meta.data$Timepoint=="12 months after infection"] <- "12 months"
Bait.positives@meta.data$Timepoint[Bait.positives@meta.data$Timepoint=="6 months after infection"] <- "6 months"
Idents(Bait.positives) <- "Timepoint"
Bait.positive.markers.overtime <- FindAllMarkers(Bait.positives, assay = "RNA", logfc.threshold = 0.1,min.pct = 0.1)
Bait.positive.markers.overtime <- Bait.positive.markers.overtime[Bait.positive.markers.overtime$p_val_adj<0.05,]
Bait.positive.markers.overtime <- Bait.positive.markers.overtime[Bait.positive.markers.overtime$avg_log2FC>0,]

DotPlot(Bait.positives, features = c("CXCR4","JUND","CD83","NFKBIA","NFKBID","GPR183","JUNB","NR4A2","MAP3K8","PTPN6","MS4A1","KLF2",
                                     "SELL","CCR6","CD79B","IL2RG","CD24"),assay = "RNA") + 
  RotatedAxis()+labs(x="",y="")+theme_classic()+theme(legend.position = "right",text = element_text(size = 20),axis.text.x = element_text(angle = 45, hjust=1))+black.axis.text() # Export 5x11
write.csv2(Bait.positive.markers.overtime, "DGEA_6mo_vs_12mo.csv")
remove(Bait.positives, Bait.positive.markers.overtime)


# Differential gene expression between positive cells 6mo vs 12 mo not vaccinated. All cells - Fig.2.I
Bait.positives <- subset(Bmem, subset = bait.positive =="yes" & Dataset!= "Second")
Bait.positives <- subset(Bait.positives, subset = Vaccinated =="no")
Bait.positives@meta.data$Timepoint[Bait.positives@meta.data$Timepoint=="12 months after infection"] <- "12 months"
Bait.positives@meta.data$Timepoint[Bait.positives@meta.data$Timepoint=="6 months after infection"] <- "6 months"
Idents(Bait.positives) <- "Timepoint"
Bait.positive.markers.overtime <- FindAllMarkers(Bait.positives, assay = "RNA", logfc.threshold = 0.1,min.pct = 0.1)
Bait.positive.markers.overtime <- Bait.positive.markers.overtime[Bait.positive.markers.overtime$p_val_adj<0.05,]
Bait.positive.markers.overtime <- Bait.positive.markers.overtime[Bait.positive.markers.overtime$avg_log2FC>0,]


# Dotplot
DotPlot(Bait.positives, features = c("CXCR4","CD83","NFKBIA","NFKBID","NR4A2","MAP3K8","JUND","BCL10","MS4A1","PTPN6",
                                     "IL2RG","TLR10","CD53"),assay = "RNA") + RotatedAxis()+labs(x="",y="")+theme_classic()+
  theme(legend.position = "right",text = element_text(size = 20),axis.text.x = element_text(angle = 45, hjust=1))+
  black.axis.text() # Export 5x11
write.xlsx(Bait.positive.markers.overtime, "DGEA_6mo_vs_12mo.xlsx")


# 5: Differential expression between named clusters - Fig.5.C
Bmem.subset <- subset(Bmem, subset = Dataset != "Second" & bait.positive=="yes")
Bmem.subset@meta.data$named.clusters <- factor(Bmem.subset@meta.data$named.clusters, levels=c("Unswitched","CD27low RM","CD27high RM",
                                                                                              "Activated","Atypical"))
Idents(Bmem.subset) <- "named.clusters"
Phenotype_Markers <- FindAllMarkers(Bmem.subset, assay = "RNA") 
Phenotype_Markers <- Phenotype_Markers[Phenotype_Markers$p_val_adj<0.05,]
write.xlsx(Phenotype_Markers,"Fig_5_C_DGEA.xlsx")											 
Phenotype_Markers.slice <- Phenotype_Markers %>% group_by(cluster) %>% slice_max(n = 20, order_by = avg_log2FC)
Features <- c("TCL1A", "IGHD", "PCDH9", "YBX3","BTG1","CCR7","CXCR4","CXCR5","BACH2","IGHM","IL4R","CD83","NFKBIA","JUNB","JUND","JCHAIN","FCER2","CR1",
              "CD40","CD24","CD44","CD69","SELL","CD27","CD70","CXCR3","COCH","RIN3","TBX21","FCRL5","MS4A1","FGR","SIGLEC6", "SIGLEC10","ITGB7","ITGB2","HLA-DPA1","HLA-DPB1","HLA-DRB1",
              "HLA-DRB5","CD74","CD86","TNFRSF1B","NEAT1","PLEK","LILRB1","LILRB2","CD79A","ITGAX","CD19","IFNGR1","NR4A2","IKZF3","ZBTB32","FCGR2B",
              "CD22","CD72","FCGR2B","LAIR1","FCRL2","FCRL3","ZEB2","DAPP1","PTPRC","LAPTM5","CSK","PTPN2","NKG7","ZAP70","GRAP2")
DoHeatmap(Bmem.subset, features = Features,group.colors = c("#DBC35E","#228833","#4477AA","#EE6677","#AA3377"))   # Export 9.5x13

# Making a better heatmap where the heatmap rows are split and clustered - using ComplexHeatmap
# Define categories
BCR_signaling <-c('TCL1A', 'IGHD', 'IGHM', 'MS4A1', 'CD79A', 'CD19', 'CD22', 'FCGR2B', 'FCRL5', 'FCRL2', 'FCRL3', 
                  'LAPTM5', 'GRAP2', 'DAPP1', 'SYK', 'MAP3K8')
Antigen_presentation <- c('HLA-DPA1', 'HLA-DPB1', 'HLA-DRB1', 'HLA-DRB5', 'CD86', 'IFI30', 'LAIR1', 'LILRB1', 'LILRB2',
                          'CD83', 'CD40', 'CD74')
Chemokine_signalling <- c('CCR7','CXCR4','CXCR5','SELL','CXCR3','FGR','ITGB7','ITGB2','ITGAX')
Cytokine_signaling <- c('IL4R', 'IFNGR1', 'TNFRSF1B', 'LTB', 'IL10RA', 'IL2RG', 'IL13RA1')
Transcription_factors <- c('BACH2','TBX21','IKZF3','ZBTB32','ZEB2','IRF1','NR4A2','NFKBIA','JUNB','JUND','JUN','TOX')
Surface_molecules <- c('FCER2','CR1','CD24','CD44','CD69','CD27','CD72','SIGLEC6','SIGLEC10')

heatmap.matrix <- as.matrix(Bmem.subset@assays$RNA@scale.data)
heatmap.matrix <- heatmap.matrix[rownames(heatmap.matrix) %in% c(BCR_signaling,Antigen_presentation,
                                                                 Chemokine_signalling,Cytokine_signaling,
                                                                 Transcription_factors, Surface_molecules),]
Unswitched.Cells <- rownames(Bmem.subset@meta.data[Bmem.subset@meta.data$named.clusters=="Unswitched",])
CD27low.Cells <- rownames(Bmem.subset@meta.data[Bmem.subset@meta.data$named.clusters=="CD27low RM",])
CD27high.Cells <- rownames(Bmem.subset@meta.data[Bmem.subset@meta.data$named.clusters=="CD27high RM",])
Activated.Cells <- rownames(Bmem.subset@meta.data[Bmem.subset@meta.data$named.clusters=="Activated",])
Atypical.Cells <- rownames(Bmem.subset@meta.data[Bmem.subset@meta.data$named.clusters=="Atypical",])

heatmap.matrix <- heatmap.matrix[,c(Unswitched.Cells,CD27low.Cells,CD27high.Cells,Activated.Cells,Atypical.Cells)]

heatmap.matrix <- heatmap.matrix[c(BCR_signaling,Antigen_presentation,Chemokine_signalling,
                                   Cytokine_signaling,Transcription_factors,
                                   Surface_molecules),]

column_split = rep("Unswitched", ncol(heatmap.matrix))
column_split[148:461] = "CD27low RM"
column_split[462:1450] = "CD27high RM"
column_split[1451:2301] = "Activated"
column_split[2302:2839] = "Atypical"
column_split <- factor(column_split, levels = c("Unswitched","CD27low RM","CD27high RM",
                                                "Activated","Atypical"))

row_split = rep("BCR Signaling", nrow(heatmap.matrix))
row_split[17:28] = "Antigen Presentation"
row_split[29:37] = "Chemokine Signaling"
row_split[38:44] = "Cytokine Signaling"
row_split[45:56] = "Transcription factors"
row_split[57:65] = "Surface Molecules"
row_split <- factor(row_split, levels = c("Surface Molecules","Antigen Presentation","BCR Signaling",
                                          "Chemokine Signaling","Cytokine Signaling",
                                          "Transcription factors"))

quantile(heatmap.matrix, c(0.1, 0.95))
col_fun = circlize::colorRamp2(c(-4, 0, 4), c("#5584B0", "#EACF65", "#E13E43"))
col_fun = circlize::colorRamp2(c(-4, 0, 4), c("#FF00FF", "#000000", "#FFFF00"))


map <- Heatmap(heatmap.matrix, cluster_columns = F, show_column_names = F,
        column_split = column_split,row_split = row_split,border = TRUE,
        column_gap = unit(1, "mm"),col=col_fun,name = "Expression",
        cluster_row_slices = F,use_raster = F,column_title_rot = 90,row_title_rot = 0,
        row_names_gp = gpar(fontsize = 8),column_names_gp = gpar(fontsize=8)) # Export 10x8

draw(map, padding = unit(c(5, 2, 8, 2), "mm")) #bottom, left, top, right paddings



# Generating Pseudo Bulk Heatmap for 8 patients (because one is only in Dataset 2) - before running this subset to Bmem.subset and run factorization of the named.clusters column from above.
Average_Expression <- AverageExpression(Bmem.subset, assays = "RNA",slot = "scale.data",
                                        features = c(BCR_signaling,Antigen_presentation,Chemokine_signalling,
                                                     Cytokine_signaling,Transcription_factors,
                                                     Surface_molecules),
                                        group.by = c("Patient","named.clusters"))
Average_Expression_df <- as.matrix(Average_Expression[[1]])
Average_Expression_df <- Average_Expression_df[,c(1,6,11,16,21,26,31,36,
                                                 2,7,12,17,22,27,32,37,
                                                 3,8,13,18,23,28,33,38,
                                                 4,9,14,19,24,29,34,39,
                                                 5,10,15,20,25,30,35,40)]

column_split = rep("Unswitched", ncol(Average_Expression_df))
column_split[9:16] = "CD27low RM"
column_split[17:24] = "CD27high RM"
column_split[25:32] = "Activated"
column_split[33:40] = "Atypical"
column_split <- factor(column_split, levels = c("Unswitched","CD27low RM","CD27high RM",
                                                "Activated","Atypical"))


# Use this section if you want to keep only vaccinated Patients for the heatmap
kicking.columns <- c(grep("PFCL1-196878-1691",colnames(Average_Expression_df)),grep("PFCL1-LIM-674950",colnames(Average_Expression_df)),grep("PFCL1-UST-190762",colnames(Average_Expression_df)))
Average_Expression_df <- Average_Expression_df[,-kicking.columns]
column_split = rep("Unswitched", ncol(Average_Expression_df))
column_split[6:10] = "CD27low RM"
column_split[11:15] = "CD27high RM"
column_split[16:20] = "Activated"
column_split[21:25] = "Atypical"
column_split <- factor(column_split, levels = c("Unswitched","CD27low RM","CD27high RM",
                                                "Activated","Atypical"))

# For looking at vaccinated patients only
colnames(Average_Expression_df)[grep("PFCL1-179308-1988",colnames(Average_Expression_df))] <- "Patient 1"
colnames(Average_Expression_df)[grep("PFCL1-199474-1994",colnames(Average_Expression_df))] <- "Patient 2"
colnames(Average_Expression_df)[grep("PFCL1-LIM-137402",colnames(Average_Expression_df))] <- "Patient 3"
colnames(Average_Expression_df)[grep("PFCL1-LIM-313000",colnames(Average_Expression_df))] <- "Patient 4"
colnames(Average_Expression_df)[grep("PFCL1-LIM828246",colnames(Average_Expression_df))] <- "Patient 5"
colnames(Average_Expression_df)[grep("PFCL1-154583-1943",colnames(Average_Expression_df))] <- "Patient 6"



# Renaming the colnames of the Average_Expression_df Matrix into Patients 1-8 (skip this if working with vaccinated patients only)
colnames(Average_Expression_df)[grep("PFCL1-179308-1988",colnames(Average_Expression_df))] <- "Patient 1"
colnames(Average_Expression_df)[grep("PFCL1-196878-1691",colnames(Average_Expression_df))] <- "Patient 2"
colnames(Average_Expression_df)[grep("PFCL1-199474-1994",colnames(Average_Expression_df))] <- "Patient 3"
colnames(Average_Expression_df)[grep("PFCL1-LIM-137402",colnames(Average_Expression_df))] <- "Patient 4"
colnames(Average_Expression_df)[grep("PFCL1-LIM-313000",colnames(Average_Expression_df))] <- "Patient 5"
colnames(Average_Expression_df)[grep("PFCL1-LIM-674950",colnames(Average_Expression_df))] <- "Patient 6"
colnames(Average_Expression_df)[grep("PFCL1-LIM828246",colnames(Average_Expression_df))] <- "Patient 7"
colnames(Average_Expression_df)[grep("PFCL1-UST-190762",colnames(Average_Expression_df))] <- "Patient 8"
colnames(Average_Expression_df)[grep("PFCL1-154583-1943",colnames(Average_Expression_df))] <- "Patient 9"


quantile(Average_Expression_df, c(0.1, 0.95))
Seurat::PurpleAndYellow()
col_fun = circlize::colorRamp2(c(-2, 0, 2), c("#5584B0", "#EACF65", "#E13E43"))
col_fun = circlize::colorRamp2(c(-2, 0, 2), c("#FF00FF", "#000000", "#FFFF00"))

map <- Heatmap(Average_Expression_df, cluster_columns = F, show_column_names = T,
        column_split = column_split,row_split = row_split,border = TRUE,
        column_gap = unit(1, "mm"),col = col_fun,
        cluster_row_slices = T,column_title_rot = 90,
        row_title_rot = 0,name = "Expression",
        row_names_gp = gpar(fontsize = 8),column_names_gp = gpar(fontsize=8)) # Export 10x7

draw(map, padding = unit(c(5, 0, 5, 2), "mm")) #bottom, left, top, right paddings




# Differential gene expression of bait positive cells between the two Resting memory compartments and the Atypicals:
Bmem.subset <- subset(Bmem, subset = Dataset != "Second" &
                          bait.positive =="yes" &
                        named.clusters != "Activated"&
                        named.clusters != "Unswitched")
RM <- rownames(Bmem.subset@meta.data[Bmem.subset@meta.data$named.clusters=="CD27low RM"|
                                       Bmem.subset@meta.data$named.clusters=="CD27high RM",])
Atypical <- rownames(Bmem.subset@meta.data[Bmem.subset@meta.data$named.clusters=="Atypical",])

Markers <- FindMarkers(Bmem, assay = "RNA", ident.1 = Atypical, ident.2 = RM)
Markers <- Markers[Markers$p_val_adj<0.05,]
Markers$gene <- rownames(Markers)

remove(Bmem.subset, RM, Atypical)


# Differential gene expression of bait positive cells between the two Resting memory compartments and the activated cells:
Bmem.subset <- subset(Bmem, subset = Dataset != "Second" &
                        bait.positive =="yes" &
                        named.clusters != "Atypical"&
                        named.clusters != "Unswitched")
RM <- rownames(Bmem.subset@meta.data[Bmem.subset@meta.data$named.clusters=="CD27low RM"|
                                       Bmem.subset@meta.data$named.clusters=="CD27high RM",])
Activated <- rownames(Bmem.subset@meta.data[Bmem.subset@meta.data$named.clusters=="Activated",])

Markers <- FindMarkers(Bmem, assay = "RNA", ident.1 = Activated, ident.2 = RM)
Markers <- Markers[Markers$p_val_adj<0.05,]
Markers$gene <- rownames(Markers)

remove(Bmem.subset, RM, Activated)


# 10: Differential gene expression of bait positive cells between Activated and Atypical compartments:
Bmem.subset <- subset(Bmem, subset = Dataset != "Second" &
                        bait.positive =="yes" &
                        named.clusters != "CD27high RM"&
                        named.clusters != "Unswitched"&
                        named.clusters != "CD27low RM")
Activated <- rownames(Bmem.subset@meta.data[Bmem.subset@meta.data$named.clusters=="Activated",])
Atypical <- rownames(Bmem.subset@meta.data[Bmem.subset@meta.data$named.clusters=="Atypical",])

Markers <- FindMarkers(Bmem, assay = "RNA", ident.1 = Activated, ident.2 = Atypical)
Markers <- Markers[Markers$p_val_adj<0.05,]
Markers$gene <- rownames(Markers)
Idents(Bmem.subset) <- "named.clusters"
DoHeatmap(Bmem.subset,assay = "RNA",features = candidates)
remove(Bmem.subset, Activated,Atypical)


######################################################################################################################
# Part 3: GSEA and GSVA input file preparation - Fig.5.F and Fig.5.G
######################################################################################################################
# Creating a function for GSEA
my_gsea <- function(category,category2,DGEAresult,number_of_pathways_to_plot,vector_of_sets_to_plot,title,Gene.sets.plotenrichment){
  # Loading the gene sets from the database
  if(missing(category2)){
  m_df<- msigdbr(species = "Homo sapiens", category = category)
  fgsea_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
  # Adding a gene set of interest manually
  Jenks_set <- read.xlsx("elife-41641-supp2-v2.xlsx") # This geneset was obtained from Zumaquero et al., 2019
  Jenks_set <- Jenks_set$Sanz_DN2_Up
  Jenks_set <- list(Jenks_set)
  fgsea_sets[length(fgsea_sets)+1] <- Jenks_set 
  names(fgsea_sets) <- c(names(fgsea_sets[])[1:length(fgsea_sets)-1],"HBO_Jenks_et_al_2018")
  # Take the DGEA result (that has to be given as input) and extract the list of differentially expressed genes from it. Ordered by -log(p_val)*sign(avg_log2FC)
  Markers <- DGEAresult %>% 
    mutate(ranking_metric = -log(p_val)*sign(avg_log2FC) ) %>%
    arrange(desc(ranking_metric)) %>% rownames_to_column(var = "feature")
  preranked_list <- Markers %>% arrange(desc(ranking_metric)) %>% 
    dplyr::select(feature, ranking_metric) %>% deframe
  set.seed(42)
  # Exectute the actual gene set enrichment using the fgsea function from the fgsea package. This takes the previously prepared list of gene sets
  # and the previously perpared list of differentially expressed genes as inputs.
  fgseaResults <- fgsea(fgsea_sets, stats = preranked_list) %>%
    as_tibble() %>%
    arrange(desc(NES))
  
  fgseaResults <- fgseaResults[fgseaResults$padj<0.05,]
  
  # The code that follows here, is just for plotting the result.
  if(missing(vector_of_sets_to_plot)){
    fgseaResults_to_plot <- fgseaResults[1:number_of_pathways_to_plot,]
    i <-1
    for (i in 1:nrow(fgseaResults_to_plot)) {
      fgseaResults_to_plot$pathway[i] <- substring(fgseaResults_to_plot$pathway[i], regexpr("_", fgseaResults_to_plot$pathway[i]) + 1, nchar(fgseaResults_to_plot$pathway[i]))
      fgseaResults_to_plot$pathway[i] <- str_replace_all(fgseaResults_to_plot$pathway[i],pattern = "_",replacement = " ")
      }
  plot <- ggplot(fgseaResults_to_plot %>% filter(padj < 0.05), aes(reorder(pathway, NES), NES)) +
    geom_col(fill="cadetblue3") +
    coord_flip() +
    labs(x="", y="Normalized Enrichment Score",
         title=title) + theme_cowplot()+center.title()
  plot(plot)
  } else{
    fgseaResults_to_plot <- fgseaResults[fgseaResults$pathway %in% vector_of_sets_to_plot,]
    i <-1
    for (i in 1:nrow(fgseaResults_to_plot)) {
      fgseaResults_to_plot$pathway[i] <- substring(fgseaResults_to_plot$pathway[i], regexpr("_", fgseaResults_to_plot$pathway[i]) + 1, nchar(fgseaResults_to_plot$pathway[i]))
      fgseaResults_to_plot$pathway[i] <- str_replace_all(fgseaResults_to_plot$pathway[i],pattern = "_",replacement = " ")
      }
    plot <- ggplot(fgseaResults_to_plot %>% filter(padj < 0.05), aes(reorder(pathway, NES), NES)) +
      geom_col(fill="cadetblue3") + theme_classic()+
      coord_flip() +
      labs(x="", y="Normalized Enrichment Score",
           title=title)+center.title()
    plot(plot)
    i <- 1
    for (i in 1:length(Gene.sets.plotenrichment)) {
      this.title <- substring(Gene.sets.plotenrichment[i], regexpr("_", Gene.sets.plotenrichment[i]) + 1, nchar(Gene.sets.plotenrichment[i]))
      this.title <- str_replace_all(this.title,pattern = "_",replacement = " ")
      plot2 <-plotEnrichment(fgsea_sets[[Gene.sets.plotenrichment[i]]],
               preranked_list) + labs(title=this.title)+labs(x="Rank",y="Enrichment score")+theme_cowplot()+black.axis.text()+center.title()
      plot(plot2)
    }
    
  }
  }else {
    m_df<- msigdbr(species = "Homo sapiens", category = category)
    m_df2<- msigdbr(species = "Homo sapiens", category = category2)
    fgsea_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
    fgsea_sets2<- m_df2 %>% split(x = .$gene_symbol, f = .$gs_name)
    fgsea_sets <- c(fgsea_sets,fgsea_sets2)
    # Adding a gene set of interest manually
    Jenks_set <- read.xlsx("elife-41641-supp2-v2.xlsx") # This geneset was obtained from Zumaquero et al., 2019
    Jenks_set <- Jenks_set$Sanz_DN2_Up
    Jenks_set <- list(Jenks_set)
    fgsea_sets[length(fgsea_sets)+1] <- Jenks_set 
    names(fgsea_sets) <- c(names(fgsea_sets[])[1:length(fgsea_sets)-1],"HBO_Jenks_et_al_2018")
    # Take the DGEA result (that has to be given as input) and extract the list of differentially expressed genes from it. Ordered by -log(p_val)*sign(avg_log2FC)
    Markers <- DGEAresult %>% 
      mutate(ranking_metric = -log(p_val)*sign(avg_log2FC) ) %>%
      arrange(desc(ranking_metric)) %>% rownames_to_column(var = "feature")
    preranked_list <- Markers %>% arrange(desc(ranking_metric)) %>% 
      dplyr::select(feature, ranking_metric) %>% deframe
    set.seed(42)
    # Exectute the actual gene set enrichment using the fgsea function from the fgsea package. This takes the previously prepared list of gene sets
    # and the previously perpared list of differentially expressed genes as inputs.
    fgseaResults <- fgsea(fgsea_sets, stats = preranked_list) %>%
      as_tibble() %>%
      arrange(desc(NES))
    
    fgseaResults <- fgseaResults[fgseaResults$padj<0.05,]
    
    # The code that follows here, is just for plotting the result.
    if(missing(vector_of_sets_to_plot)){
      fgseaResults_to_plot <- fgseaResults[1:number_of_pathways_to_plot,]
      i <-1
      for (i in 1:nrow(fgseaResults_to_plot)) {
        fgseaResults_to_plot$pathway[i] <- substring(fgseaResults_to_plot$pathway[i], regexpr("_", fgseaResults_to_plot$pathway[i]) + 1, nchar(fgseaResults_to_plot$pathway[i]))
        fgseaResults_to_plot$pathway[i] <- str_replace_all(fgseaResults_to_plot$pathway[i],pattern = "_",replacement = " ")
      }
      plot <- ggplot(fgseaResults_to_plot %>% filter(padj < 0.05), aes(reorder(pathway, NES), NES)) +
        geom_col(fill="cadetblue3") +
        coord_flip() +
        labs(x="", y="Normalized Enrichment Score",
             title=title) + theme_cowplot()+center.title()
      plot(plot)
    } else{
      fgseaResults_to_plot <- fgseaResults[fgseaResults$pathway %in% vector_of_sets_to_plot,]
      i <-1
      for (i in 1:nrow(fgseaResults_to_plot)) {
        fgseaResults_to_plot$pathway[i] <- substring(fgseaResults_to_plot$pathway[i], regexpr("_", fgseaResults_to_plot$pathway[i]) + 1, nchar(fgseaResults_to_plot$pathway[i]))
        fgseaResults_to_plot$pathway[i] <- str_replace_all(fgseaResults_to_plot$pathway[i],pattern = "_",replacement = " ")
      }
      plot <- ggplot(fgseaResults_to_plot %>% filter(padj < 0.05), aes(reorder(pathway, NES), NES)) +
        geom_col(fill="cadetblue3") + theme_classic()+
        coord_flip() +
        labs(x="", y="Normalized Enrichment Score",
             title=title)+center.title()
      plot(plot)
      i <- 1
      for (i in 1:length(Gene.sets.plotenrichment)) {
        this.title <- substring(Gene.sets.plotenrichment[i], regexpr("_", Gene.sets.plotenrichment[i]) + 1, nchar(Gene.sets.plotenrichment[i]))
        this.title <- str_replace_all(this.title,pattern = "_",replacement = " ")
        plot2 <-plotEnrichment(fgsea_sets[[Gene.sets.plotenrichment[i]]],
                               preranked_list) + labs(title=this.title)+labs(x="Rank",y="Enrichment score")+theme_cowplot()+black.axis.text()+center.title()
        plot(plot2)
      }
      
    }
  }

  return(fgseaResults)
}


Gene.sets <- c("HBO_Jenks_et_al_2018","GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_PEPTIDE_ANTIGEN","GOMF_MOLECULAR_TRANSDUCER_ACTIVITY",
               "GOCC_MHC_PROTEIN_COMPLEX","GOBP_REGULATION_OF_CELL_ACTIVATION","GOBP_RESPONSE_TO_INTERFERON_GAMMA",
               "GOBP_INTERFERON_GAMMA_MEDIATED_SIGNALING_PATHWAY","GOCC_COPII_COATED_ER_TO_GOLGI_TRANSPORT_VESICLE",
               "GOBP_INTEGRIN_MEDIATED_SIGNALING_PATHWAY","GOBP_CELL_ADHESION_MEDIATED_BY_INTEGRIN")

Gene.sets.plotenrichment <- c("HALLMARK_INTERFERON_GAMMA_RESPONSE","HBO_Jenks_et_al_2018","GOBP_INTERFERON_GAMMA_MEDIATED_SIGNALING_PATHWAY","GOBP_B_CELL_RECEPTOR_SIGNALING_PATHWAY",
                              "GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_PEPTIDE_ANTIGEN","GOBP_INTEGRIN_MEDIATED_SIGNALING_PATHWAY","GOBP_RESPONSE_TO_INTERFERON_GAMMA")

GSEA <- my_gsea(category = "H", category2 = "C5" ,DGEAresult = Markers ,number_of_pathways_to_plot = 10 ,
                vector_of_sets_to_plot = Gene.sets,title = "Hallmark Sets", Gene.sets.plotenrichment = Gene.sets.plotenrichment)

remove(GSEA, Markers, Gene.sets)

# Making the GSEA output exportable
GSEA <- apply(GSEA,2,as.character)
write.csv2(GSEA[,c(1,3,5,8)], "GO Sets enriched in Activated.vs.Atypical.csv")


# GSVEA - RM vs Activated vs Atypical
# Prepare a matrix with rows corresponding to genes and columns to samples
RM <- Bmem@meta.data[Bmem@meta.data$named.clusters=="CD27high RM" & Bmem@meta.data$bait.positive=="yes"|
                       Bmem@meta.data$named.clusters=="CD27low RM"& Bmem@meta.data$bait.positive=="yes","Full.Row.names"]
Activated <- Bmem@meta.data[Bmem@meta.data$named.clusters=="Activated" & Bmem@meta.data$bait.positive=="yes","Full.Row.names"]
Atypical <- Bmem@meta.data[Bmem@meta.data$named.clusters=="Atypical" & Bmem@meta.data$bait.positive=="yes","Full.Row.names"]

Bmem.subset <- subset(Bmem, subset = Full.Row.names %in% RM |
                        Full.Row.names %in% Activated |
                        Full.Row.names %in% Atypical)

saveRDS(Bmem.subset, "Bmem.subset.rds") # -> Use this as input in GSVA script.


######################################################################################################################
# Part 4: Clonal and mutational analysis and LIBRA score analysis
######################################################################################################################
# Fig.2.G
# A: How does the mutational load of COVID specific cells change over time? - Looking only at non vaccinated cells.
Bmem_fourth <- Bmem_fourth@meta.data
Bmem_fourth <- Bmem_fourth[rownames(Bmem_fourth)%in%Naives,] 

# Read in the mutational load file
# The pre processed Bmem Seurat object contains the mutational counts of all the cells which are inside it already. 
# However, the cells from Dataset 4, which are needed here because they contain the naive cells, don't have mutational counts assinged yet.
# Therefore I am loading a file which we have prepared with help of the Immcanation pipeline, which contains the naive cells from Dataset 4 - and their 
# mutational counts. This file is included in the Zenodo upload of the dataset.
Mutations <- read.table(file = 'mutational_loads.tsv', header = TRUE)

# Prepare the file so that it can be merged with the Seurat object
i <-1
for (i in 1:nrow(Mutations)) {
  if(substr(Mutations$cell_id[i],1,2)=="CC"){
    Mutations$cell_id[i] <- paste0("Dataset_4_",substr(Mutations$cell_id[i],3,nchar(Mutations$cell_id[i])))
  }
  if(substr(Mutations$cell_id[i],1,2)=="AA"){
    Mutations$cell_id[i] <- paste0("Dataset_1_",substr(Mutations$cell_id[i],3,nchar(Mutations$cell_id[i])))
  }
  if(substr(Mutations$cell_id[i],1,2)=="GG"){
    Mutations$cell_id[i] <- paste0("Dataset_2_",substr(Mutations$cell_id[i],3,nchar(Mutations$cell_id[i])))
  }
  if(substr(Mutations$cell_id[i],1,2)=="TT"){
    Mutations$cell_id[i] <- paste0("Dataset_3_",substr(Mutations$cell_id[i],3,nchar(Mutations$cell_id[i])))
  }
}
rownames(Mutations) <- Mutations$cell_id
Mutations <- Mutations[,c(49,63,64)]

# Do the merge
Bmem_fourth <- merge(Bmem_fourth, Mutations, by = 0, all.x = T)
Bmem_fourth <- Bmem_fourth[,2:ncol(Bmem_fourth)]
rownames(Bmem_fourth) <- Bmem_fourth$Full.Row.names

remove(Mutations)

Bmem.subset <- subset(Bmem, subset = mu_freq >0 |
                        mu_freq ==0)
Bmem.subset <- subset(Bmem.subset, subset = bait.positive =="yes" &
                        sample_group != "12 months vaccinated")

# Use this if you want to look only at the three not vaccinated patients.
#Bmem.subset <- subset(Bmem.subset, subset = Patient == "PFCL1-196878-1691" |
#                        Patient == "PFCL1-LIM-674950" | Patient == "PFCL1-UST-190762")

df <- Bmem.subset@meta.data
df$sample_group[df$sample_group=="12 months not vaccinated"] <- "12 months"
df <- df[,c(22,46,47)]
Bmem_fourth$sample_group <- "Naive"
Bmem_fourth <- Bmem_fourth[,c(22, 48,49)]
df <- rbind(df,Bmem_fourth)
df$sample_group <- factor(df$sample_group, levels=c("Naive","6 months","12 months"))
df <- drop_na(df)

# Visualize: Specify the comparisons you want
stat.test <- df %>% rstatix::wilcox_test(mu_count~sample_group)

# Selecting the pvalues we want to show
stat.test <- stat.test[c(3),]

ggboxplot(df, x = "sample_group", y = "mu_count",
          fill = "sample_group")+ 
  stat_pvalue_manual(stat.test, label = "p.adj", y.position = c(62))+
  scale_fill_manual(values=c("black","#004488","#997700"))+ ggtitle("")+
  ylab("SHM count")+xlab("")+ theme(legend.position="right")+labs(fill="Sample group")+
  theme(axis.text.x = element_text(angle = 45, hjust=1)) # Export 5x6

remove(Bmem.subset, df, my_comparisons,stat.test)

# Fig.4.E
# B: How does the mutational load of COVID specific cells change over time? - Looking at non vaccinated cells and cells after vaccination.
Bmem_fourth <- Bmem_fourth@meta.data
Bmem_fourth <- Bmem_fourth[rownames(Bmem_fourth)%in%Naives,] 

# Read in the mutational load file
Mutations <- read.table(file = 'mutational_loads.tsv', header = TRUE) # For inormation about this file see above in part 4.

# Prepare the file so that it can be merged with the Seurat object
i <-1
for (i in 1:nrow(Mutations)) {
  if(substr(Mutations$cell_id[i],1,2)=="CC"){
    Mutations$cell_id[i] <- paste0("Dataset_4_",substr(Mutations$cell_id[i],3,nchar(Mutations$cell_id[i])))
  }
  if(substr(Mutations$cell_id[i],1,2)=="AA"){
    Mutations$cell_id[i] <- paste0("Dataset_1_",substr(Mutations$cell_id[i],3,nchar(Mutations$cell_id[i])))
  }
  if(substr(Mutations$cell_id[i],1,2)=="GG"){
    Mutations$cell_id[i] <- paste0("Dataset_2_",substr(Mutations$cell_id[i],3,nchar(Mutations$cell_id[i])))
  }
  if(substr(Mutations$cell_id[i],1,2)=="TT"){
    Mutations$cell_id[i] <- paste0("Dataset_3_",substr(Mutations$cell_id[i],3,nchar(Mutations$cell_id[i])))
  }
}
rownames(Mutations) <- Mutations$cell_id
Mutations <- Mutations[,c(49,63,64)]
# Do the merge
Bmem_fourth <- merge(Bmem_fourth, Mutations, by = 0, all.x = T)
Bmem_fourth <- Bmem_fourth[,2:ncol(Bmem_fourth)]
rownames(Bmem_fourth) <- Bmem_fourth$Full.Row.names

remove(Mutations)

Bmem.subset <- subset(Bmem, subset = mu_freq >0 |
                        mu_freq ==0)
Bmem.subset <- subset(Bmem.subset, subset = bait.positive =="yes")

df <- Bmem.subset@meta.data
df$sample_group[df$sample_group=="12 months vaccinated" & df$time_after_vac=="15d"|
                  df$sample_group=="12 months vaccinated" & df$time_after_vac=="23d"] <- "12 months vaccinated <23d"
df$sample_group[df$sample_group=="12 months vaccinated" & df$time_after_vac=="85d"|
                  df$sample_group=="12 months vaccinated" & df$time_after_vac=="87d"|
                  df$sample_group=="12 months vaccinated" & df$time_after_vac=="108d"] <- "12 months vaccinated >85d"
df <- df[,c(22,46,47)]
Bmem_fourth$sample_group <- "Naive"
Bmem_fourth <- Bmem_fourth[,c(22, 48,49)]
df <- rbind(df,Bmem_fourth)
df$sample_group <- factor(df$sample_group, levels=c("Naive","6 months","12 months not vaccinated","12 months vaccinated <23d","12 months vaccinated >85d"))
df <- drop_na(df)

# Stat_compare means fails to show the adj. p value, therefore we use this more custom way to bring the p values to the boxplot:
# https://github.com/kassambara/rstatix
stat.test <- df %>% rstatix::wilcox_test(mu_count~sample_group)

# Selecting the pvalues we want to show
stat.test <- stat.test[c(5,6,9,10),]

ggboxplot(df, x = "sample_group", y = "mu_count",
          fill = "sample_group")+
  stat_pvalue_manual(stat.test, label = "p.adj", 
                     y.position = c(62,67,72,62))+
  scale_fill_manual(values=c("black","#004488","#997700","#d5fa43","#35e690"))+ ggtitle("")+
  ylab("SHM count")+xlab("")+ theme(legend.position="right")+labs(fill="Sample group")+
  theme(axis.text.x = element_text(angle = 45, hjust=1))

remove(Bmem.subset, df, my_comparisons) # Export 5.5x8


# Fig.6.F
# Is there a difference in SHMs between spike specific cells of different named clusters? # Export 6x7
Bmem.subset <- subset(Bmem, subset = mu_freq >0 |
                        mu_freq ==0)
Bmem.subset <- subset(Bmem.subset, subset = bait.positive =="yes"&
                        sample_group=="12 months vaccinated")
df <- Bmem.subset@meta.data
df$named.clusters <- factor(df$named.clusters, levels=c("Unswitched","CD27low RM",
                                              "CD27high RM","Activated","Atypical"))
# Run test
stat.test <- df %>% rstatix::wilcox_test(mu_count~named.clusters)

# Selecting the pvalues we want to show - this command is needed if pvalues of only some comparisons should be shown. Otherwise skip.
#stat.test <- stat.test[c(1,5,8,10),]

ggboxplot(df, x = "named.clusters", y = "mu_count",
          fill = "named.clusters")+
  stat_pvalue_manual(stat.test, label = "p.adj", 
                     y.position = c(62,67,72,77,82,87,92,57,62,82))+
  scale_fill_manual(values=c("#DBC35E","#228833","#4477AA","#EE6677","#AA3377"))+ ggtitle("")+
  ylab("SHM count")+xlab("")+ theme(legend.position="right")+labs(fill="B Cell state")+
  theme(axis.text.x = element_text(angle = 75, hjust=1))

remove(Bmem.subset, df)



# Mapping the receptor sequences of our bait specific cells to the CoV-AbDab
CoVAbDab <- read.csv("~/NAS/Jan/Experiments and Data/210801 - Exp 035 - Profiling of COVID vaccinated patient samples/CoV-AbDab_260722.csv")
setwd("~/NAS/Jan/Experiments and Data/210801 - Exp 035 - Profiling of COVID vaccinated patient samples/Memory B cell Fourth Experiment 220323/multi_exp035_4")  

All_cells <- read_airr("all_cells_0.20_heavy_germ-pass.tsv")
i <-1
for (i in 1:nrow(All_cells)) {
  if(substr(All_cells$cell_id[i],1,2)=="CC"){
    All_cells$cell_id[i] <- paste0("Dataset_4_",substr(All_cells$cell_id[i],3,nchar(All_cells$cell_id[i])))
  }
  if(substr(All_cells$cell_id[i],1,2)=="AA"){
    All_cells$cell_id[i] <- paste0("Dataset_1_",substr(All_cells$cell_id[i],3,nchar(All_cells$cell_id[i])))
  }
  if(substr(All_cells$cell_id[i],1,2)=="GG"){
    All_cells$cell_id[i] <- paste0("Dataset_2_",substr(All_cells$cell_id[i],3,nchar(All_cells$cell_id[i])))
  }
  if(substr(All_cells$cell_id[i],1,2)=="TT"){
    All_cells$cell_id[i] <- paste0("Dataset_3_",substr(All_cells$cell_id[i],3,nchar(All_cells$cell_id[i])))
  }
}

All_cells <- All_cells[All_cells$cell_id %in% rownames(Bmem@meta.data[Bmem@meta.data$bait.positive=="yes",]),c(5,6,7,48,49)]
All_cells_cdr3 <- translate(DNAStringSet(All_cells$cdr3))
All_cells_cdr3 <- as.data.frame(All_cells_cdr3)
All_cells$cdr3_as <- as.vector(All_cells_cdr3[,1])
length(intersect(All_cells$cdr3_as,CoVAbDab$CDRH3))
length(unique(CoVAbDab$CDRH3))


# Including heavy and light chains
# Creating first a file for all the heavy chain information
setwd("~/NAS/Jan/Experiments and Data/210801 - Exp 035 - Profiling of COVID vaccinated patient samples/Immcantation Analysis/Immcantation Outputs/diversity")
PFCL1_154583_1943_vacc_early <- read_airr("./PFCL1_154583_1943_vacc_early/outputs/PFCL1_154583_1943_vacc_early_0.20/PFCL1_154583_1943_vacc_early_0.20_heavy_germ-pass.tsv")
PFCL1_179308_1988_vacc_late <- read_airr("./PFCL1_179308_1988_vacc_late/outputs/PFCL1_179308_1988_vacc_late_0.20/PFCL1_179308_1988_vacc_late_0.20_heavy_germ-pass.tsv")
PFCL1_199474_1994_vacc_late <- read_airr("./PFCL1_199474_1994_vacc_late/outputs/PFCL1_199474_1994_vacc_late_0.20/PFCL1_199474_1994_vacc_late_0.20_heavy_germ-pass.tsv")
PFCL1_LIM_137402_vacc_early <- read_airr("./PFCL1_LIM_137402_vacc_early/outputs/PFCL1_LIM_137402_vacc_early_0.20/PFCL1_LIM_137402_vacc_early_0.20_heavy_germ-pass.tsv")
PFCL1_LIM_313000_vacc_early <- read_airr("./PFCL1_LIM_313000_vacc_early/outputs/PFCL1_LIM_313000_vacc_early_0.20/PFCL1_LIM_313000_vacc_early_0.20_heavy_germ-pass.tsv")
PFCL1_LIM828246_vacc_late <- read_airr("./PFCL1_LIM828246_vacc_late/outputs/PFCL1_LIM828246_vacc_late_0.20/PFCL1_LIM828246_vacc_late_0.20_heavy_germ-pass.tsv")

PFCL1_154583_1943_6mo <- read_airr("./PFCL1_154583_1943_6mo/outputs/PFCL1_154583_1943_6mo_0.20/PFCL1_154583_1943_6mo_0.20_heavy_germ-pass.tsv")
PFCL1_179308_1988_6mo <- read_airr("./PFCL1_179308_1988_6mo/outputs/PFCL1_179308_1988_6mo_0.20/PFCL1_179308_1988_6mo_0.20_heavy_germ-pass.tsv")
PFCL1_199474_1994_6mo <- read_airr("./PFCL1_199474_1994_6mo/outputs/PFCL1_199474_1994_6mo_0.20/PFCL1_199474_1994_6mo_0.20_heavy_germ-pass.tsv")
PFCL1_LIM_137402_6mo <- read_airr("./PFCL1_LIM_137402_6mo/outputs/PFCL1_LIM_137402_6mo_0.20/PFCL1_LIM_137402_6mo_0.20_heavy_germ-pass.tsv")
PFCL1_LIM_313000_6mo <- read_airr("./PFCL1_LIM_313000_6mo/outputs/PFCL1_LIM_313000_6mo_0.20/PFCL1_LIM_313000_6mo_0.20_heavy_germ-pass.tsv")
PFCL1_LIM828246_6mo <- read_airr("./PFCL1_LIM828246_6mo/outputs/PFCL1_LIM828246_6mo_0.20/PFCL1_LIM828246_6mo_0.20_heavy_germ-pass.tsv")

PFCL1_196878_1691_6mo <- read_airr("./PFCL1_196878_1691_6mo/outputs/PFCL1_196878_1691_6mo_0.20/PFCL1_196878_1691_6mo_0.20_heavy_germ-pass.tsv")
PFCL1_LIM_674950_6mo <- read_airr("./PFCL1_LIM_674950_6mo/outputs/PFCL1_LIM_674950_6mo_0.20/PFCL1_LIM_674950_6mo_0.20_heavy_germ-pass.tsv")
PFCL1_UST_190762_6mo <- read_airr("./PFCL1_UST_190762_6mo/outputs/PFCL1_UST_190762_6mo_0.20/PFCL1_UST_190762_6mo_0.20_heavy_germ-pass.tsv")
PFCL1_196878_1691_12mo <- read_airr("./PFCL1_196878_1691_12mo/outputs/PFCL1_196878_1691_12mo_0.20/PFCL1_196878_1691_12mo_0.20_heavy_germ-pass.tsv")
PFCL1_LIM_674950_12mo <- read_airr("./PFCL1_LIM_674950_12mo/outputs/PFCL1_LIM_674950_12mo_0.20/PFCL1_LIM_674950_12mo_0.20_heavy_germ-pass.tsv")
PFCL1_UST_190762_12mo <- read_airr("./PFCL1_UST_190762_12mo/outputs/PFCL1_UST_190762_12mo_0.20/PFCL1_UST_190762_12mo_0.20_heavy_germ-pass.tsv")

sixmonths <- rbind(PFCL1_154583_1943_6mo,PFCL1_179308_1988_6mo,PFCL1_199474_1994_6mo,PFCL1_LIM_137402_6mo,PFCL1_LIM_313000_6mo,PFCL1_LIM828246_6mo)
twelvemonths <- rbind(PFCL1_154583_1943_vacc_early,PFCL1_179308_1988_vacc_late,PFCL1_199474_1994_vacc_late,
                      PFCL1_LIM_137402_vacc_early,PFCL1_LIM_313000_vacc_early,PFCL1_LIM828246_vacc_late)

# Merging the files
All <- rbind(sixmonths,twelvemonths)

# Then for the light chains
PFCL1_154583_1943_vacc_early <- read_airr("./PFCL1_154583_1943_vacc_early/outputs/PFCL1_154583_1943_vacc_early_0.20/PFCL1_154583_1943_vacc_early_0.20_light_productive-T.tsv")
PFCL1_179308_1988_vacc_late <- read_airr("./PFCL1_179308_1988_vacc_late/outputs/PFCL1_179308_1988_vacc_late_0.20/PFCL1_179308_1988_vacc_late_0.20_light_productive-T.tsv")
PFCL1_199474_1994_vacc_late <- read_airr("./PFCL1_199474_1994_vacc_late/outputs/PFCL1_199474_1994_vacc_late_0.20/PFCL1_199474_1994_vacc_late_0.20_light_productive-T.tsv")
PFCL1_LIM_137402_vacc_early <- read_airr("./PFCL1_LIM_137402_vacc_early/outputs/PFCL1_LIM_137402_vacc_early_0.20/PFCL1_LIM_137402_vacc_early_0.20_light_productive-T.tsv")
PFCL1_LIM_313000_vacc_early <- read_airr("./PFCL1_LIM_313000_vacc_early/outputs/PFCL1_LIM_313000_vacc_early_0.20/PFCL1_LIM_313000_vacc_early_0.20_light_productive-T.tsv")
PFCL1_LIM828246_vacc_late <- read_airr("./PFCL1_LIM828246_vacc_late/outputs/PFCL1_LIM828246_vacc_late_0.20/PFCL1_LIM828246_vacc_late_0.20_light_productive-T.tsv")

PFCL1_154583_1943_6mo <- read_airr("./PFCL1_154583_1943_6mo/outputs/PFCL1_154583_1943_6mo_0.20/PFCL1_154583_1943_6mo_0.20_light_productive-T.tsv")
PFCL1_179308_1988_6mo <- read_airr("./PFCL1_179308_1988_6mo/outputs/PFCL1_179308_1988_6mo_0.20/PFCL1_179308_1988_6mo_0.20_light_productive-T.tsv")
PFCL1_199474_1994_6mo <- read_airr("./PFCL1_199474_1994_6mo/outputs/PFCL1_199474_1994_6mo_0.20/PFCL1_199474_1994_6mo_0.20_light_productive-T.tsv")
PFCL1_LIM_137402_6mo <- read_airr("./PFCL1_LIM_137402_6mo/outputs/PFCL1_LIM_137402_6mo_0.20/PFCL1_LIM_137402_6mo_0.20_light_productive-T.tsv")
PFCL1_LIM_313000_6mo <- read_airr("./PFCL1_LIM_313000_6mo/outputs/PFCL1_LIM_313000_6mo_0.20/PFCL1_LIM_313000_6mo_0.20_light_productive-T.tsv")
PFCL1_LIM828246_6mo <- read_airr("./PFCL1_LIM828246_6mo/outputs/PFCL1_LIM828246_6mo_0.20/PFCL1_LIM828246_6mo_0.20_light_productive-T.tsv")

PFCL1_196878_1691_6mo <- read_airr("./PFCL1_196878_1691_6mo/outputs/PFCL1_196878_1691_6mo_0.20/PFCL1_196878_1691_6mo_0.20_light_productive-T.tsv")
PFCL1_LIM_674950_6mo <- read_airr("./PFCL1_LIM_674950_6mo/outputs/PFCL1_LIM_674950_6mo_0.20/PFCL1_LIM_674950_6mo_0.20_light_productive-T.tsv")
PFCL1_UST_190762_6mo <- read_airr("./PFCL1_UST_190762_6mo/outputs/PFCL1_UST_190762_6mo_0.20/PFCL1_UST_190762_6mo_0.20_light_productive-T.tsv")
PFCL1_196878_1691_12mo <- read_airr("./PFCL1_196878_1691_12mo/outputs/PFCL1_196878_1691_12mo_0.20/PFCL1_196878_1691_12mo_0.20_light_productive-T.tsv")
PFCL1_LIM_674950_12mo <- read_airr("./PFCL1_LIM_674950_12mo/outputs/PFCL1_LIM_674950_12mo_0.20/PFCL1_LIM_674950_12mo_0.20_light_productive-T.tsv")
PFCL1_UST_190762_12mo <- read_airr("./PFCL1_UST_190762_12mo/outputs/PFCL1_UST_190762_12mo_0.20/PFCL1_UST_190762_12mo_0.20_light_productive-T.tsv")

sixmonths <- rbind(PFCL1_154583_1943_6mo,PFCL1_179308_1988_6mo,PFCL1_199474_1994_6mo,PFCL1_LIM_137402_6mo,PFCL1_LIM_313000_6mo,PFCL1_LIM828246_6mo)
twelvemonths <- rbind(PFCL1_154583_1943_vacc_early,PFCL1_179308_1988_vacc_late,PFCL1_199474_1994_vacc_late,
                      PFCL1_LIM_137402_vacc_early,PFCL1_LIM_313000_vacc_early,PFCL1_LIM828246_vacc_late)


# Merging the files
All_light <- rbind(sixmonths,twelvemonths)

# Merging everything by cell_id
All <- merge(All,All_light, by="cell_id",all.x=T)
All <- All[,c(1,49,110)]

All_cells_cdr3_heavy <- translate(DNAStringSet(All$cdr3.x))
All_cells_cdr3_light <- translate(DNAStringSet(All$cdr3.y))

All_cells_cdr3_heavy <- as.data.frame(All_cells_cdr3_heavy)
All_cells_cdr3_light <- as.data.frame(All_cells_cdr3_light)

All$cdr3_heavy_as <- as.vector(All_cells_cdr3_heavy[,1])
All$cdr3_light_as <- as.vector(All_cells_cdr3_light[,1])

intersect(All$cdr3_heavy_as,CoVAbDab$CDRH3)
intersect(All$cdr3_light_as,CoVAbDab$CDRL3)

All$combiend.cdr3 <- paste0(All$cdr3_heavy_as,All$cdr3_light_as)
CoVAbDab$CDR3_complete <- paste0(CoVAbDab$CDRH3,CoVAbDab$CDRL3)

intersect(All$combiend.cdr3,CoVAbDab$CDR3_complete)
CoVAbDab[CoVAbDab$CDR3_complete %in% intersect(All$combiend.cdr3,CoVAbDab$CDR3_complete),]


######################################################################################################################
# Part 5: Monocle 3 Analysis
######################################################################################################################
# Manually constructing a Monocle object without transferring a Seurat Object to a Monocle object.
# 1 Create one expression Matrix with the counts data for each dataset (Datasets 1, 3, 4) - Dataset 2 is excluded.

Dataset1.cells <- rownames(Bmem@meta.data[Bmem@meta.data$Dataset=="First",])
Dataset3.cells <- rownames(Bmem@meta.data[Bmem@meta.data$Dataset=="Third",])
Dataset4.cells <- rownames(Bmem@meta.data[Bmem@meta.data$Dataset=="Fourth",])

Count_matrix_Dataset_1 <- Bmem@assays$RNA@counts[,colnames(Bmem@assays$RNA@counts)%in% Dataset1.cells]
Count_matrix_Dataset_3 <- Bmem@assays$RNA@counts[,colnames(Bmem@assays$RNA@counts)%in% Dataset3.cells]
Count_matrix_Dataset_4 <- Bmem@assays$RNA@counts[,colnames(Bmem@assays$RNA@counts)%in% Dataset4.cells]

gene_metadata <- data.frame("gene_short_name" = rownames(Bmem@assays$RNA@counts),"Species"="Human")
rownames(gene_metadata) <- gene_metadata$gene_short_name

cell_metadata_Dataset_1 <- Bmem@meta.data[Bmem@meta.data$Full.Row.names %in% Dataset1.cells, c(18,19,20,21,22,37,42,53)]
cell_metadata_Dataset_3 <- Bmem@meta.data[Bmem@meta.data$Full.Row.names %in% Dataset3.cells, c(18,19,20,21,22,37,42,53)]
cell_metadata_Dataset_4 <- Bmem@meta.data[Bmem@meta.data$Full.Row.names %in% Dataset4.cells, c(18,19,20,21,22,37,42,53)]

# Renaiming colnames to make the plot legend better
colnames(cell_metadata_Dataset_1)[8] <- "Subset"
colnames(cell_metadata_Dataset_3)[8] <- "Subset"
colnames(cell_metadata_Dataset_4)[8] <- "Subset"


# Creating Monocle Datasets
cds_Dataset_1 <- new_cell_data_set(Count_matrix_Dataset_1,
                         cell_metadata = cell_metadata_Dataset_1,
                         gene_metadata = gene_metadata)
cds_Dataset_3 <- new_cell_data_set(Count_matrix_Dataset_3,
                                   cell_metadata = cell_metadata_Dataset_3,
                                   gene_metadata = gene_metadata)
cds_Dataset_4 <- new_cell_data_set(Count_matrix_Dataset_4,
                                   cell_metadata = cell_metadata_Dataset_4,
                                   gene_metadata = gene_metadata)

cds_list <- list(cds_Dataset_1,cds_Dataset_3,cds_Dataset_4)

# Monocle pipeline
Bmem_cds <- combine_cds(cds_list = cds_list, cell_names_unique = T)
Bmem_cds <- preprocess_cds(Bmem_cds, num_dim = 20)
plot_pc_variance_explained(Bmem_cds)
Bmem_cds <- align_cds(Bmem_cds, alignment_group = "sample")
Bmem_cds <- reduce_dimension(Bmem_cds)

# Cluster the cells
Bmem_cds <- cluster_cells(Bmem_cds)

# When learning trajectories, each partition will eventually become a separate trajectory.
# Exception: Call use_partition = F
Bmem_cds <- learn_graph(Bmem_cds,use_partition = F)
# Supp.Fig.6.C
plot_cells(Bmem_cds,
           color_cells_by = "cluster",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,group_label_size = 5,label_cell_groups=F,
           label_principal_points = F,cell_size = 0.4)+ggtitle("Monocle Clustering and Trajectory")+ black.axis.text()+
  theme_void()+
  theme(plot.margin = margin(t = 2,  # Top margin
                             r = 15,  # Right margin
                             b = 2,  # Bottom margin
                             l = 2))+ # Left margin
  center.title() # Export 5x7

# Fig.6.G.1
plot_cells(Bmem_cds,
           color_cells_by = "Subset",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE, 
           label_branch_points=FALSE,group_label_size = 5,label_cell_groups=F,
           label_principal_points = F,cell_size = 0.4)+ scale_color_manual(values = rev(c("#DBC35E","#228833","#4477AA","#AA3377","#EE6677")))+
           ggtitle("Monocle Clustering and Trajectory\nColored by Subsets")+ black.axis.text()+
    theme(plot.margin = margin(t = 2,  # Top margin
                             r = 15,  # Right margin
                             b = 2,  # Bottom margin
                             l = 2))+ # Left margin
  center.title() # Export 5x7

# Creating plots for a subset of the cds object that only contains bait positive cells
Bait.positives <- colnames(Bmem_cds)[colData(Bmem_cds)$bait.positive=="yes"]
Bmem_cds_subset <- Bmem_cds[,Bait.positives]
# Fig.6.G.2
plot_cells(Bmem_cds_subset,
           color_cells_by = "Subset",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE, 
           label_branch_points=FALSE,group_label_size = 5,label_cell_groups=F,
           label_principal_points = F,cell_size = 0.4)+ scale_color_manual(values = rev(c("#DBC35E","#228833","#4477AA","#AA3377","#EE6677")))+
  ggtitle("Monocle Clustering and Trajectory\nColored by Subsets, Bait positive cells only")+ black.axis.text()+
  theme(plot.margin = margin(t = 2,  # Top margin
                             r = 15,  # Right margin
                             b = 2,  # Bottom margin
                             l = 2))+ # Left margin
  center.title() # Export 5x7


# Ordering the cells
Bmem_cds <- order_cells(Bmem_cds, root_pr_nodes="Y_141")
# Supp.Fig.6.E
plot_cells(Bmem_cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)+theme_void()+ggtitle("Monocle cells colored by Pseudotime") + center.title() # Export 5x7


# Analyzing branches in single-cell trajectories
# First I am subsetting the Moncle object to only those cells which are "late" in Pseudotime
cells_to_subset <- names(Bmem_cds@principal_graph_aux@listData$UMAP$pseudotime[Bmem_cds@principal_graph_aux@listData$UMAP$pseudotime>22])
Bmem_cds_subset <- Bmem_cds[,cells_to_subset]

# Subsetting to only bait specific cells
Bait.positives <- colnames(Bmem_cds_subset)[colData(Bmem_cds_subset)$bait.positive=="yes"]
Bmem_cds_subset <- Bmem_cds_subset[,Bait.positives]


plot_cells(Bmem_cds_subset,
           color_cells_by = "Subset",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,group_label_size = 5,label_cell_groups=F,
           label_principal_points = F,cell_size = 0.4)+ggtitle("Subsetted Monocle Clustering and Trajectory")+ black.axis.text()+ theme_void()+
           theme(plot.margin = margin(t = 2,r = 15,b = 2,l = 2))+
           scale_color_manual(values = rev(c("#DBC35E","#228833","#4477AA","#AA3377","#EE6677")))+
           center.title() # Export 3x5

# Then, I run the "graph_test" and the "find_gene_modules" function to identify modules of genes which are relevant through pseudotime.
# The find_gene_modules function seems to give different results every time the analysis is run.
# However, the results are similar to each other. Every time, I see the gene module that is high in Atypicals and in CD27high RM, but lower
# in the other subsets.
Bmem_cds_subset_graph_test <- graph_test(Bmem_cds_subset, neighbor_graph="principal_graph", cores=4)
Bmem_cds_subset_graph_test <- Bmem_cds_subset_graph_test[Bmem_cds_subset_graph_test$q_value < 0.00000000000000000000005,]
Bmem_cds_subset_graph_test_genes <- row.names(subset(Bmem_cds_subset_graph_test, q_value < 0.00000000000000000000005))
gene_module_df <- find_gene_modules(Bmem_cds_subset[Bmem_cds_subset_graph_test_genes,], 
                                    resolution=0.05) # The resolution parameter defines how many gene modules will be identified. Smaller values yield less modules.

plot_cells(Bmem_cds_subset,
           genes=gene_module_df,
           label_cell_groups=FALSE,
           show_trajectory_graph=T,label_branch_points = F, label_leaves = F)

cell_group_df <- tibble::tibble(cell=row.names(colData(Bmem_cds_subset)), 
                                cell_group=colData(Bmem_cds_subset)$Subset)
agg_mat <- aggregate_gene_expression(Bmem_cds_subset, gene_module_df, cell_group_df)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))

# Supp.Fig.6.F
pheatmap::pheatmap(agg_mat, scale="column", clustering_method="ward.D2",treeheight_col=0, main = "Monocle Gene Module Scores",
                   col = inferno(100),treeheight_row=15,) # Export 4x4


# Focusing on interesting modules:
plot_cells(Bmem_cds_subset,
          genes=gene_module_df %>% filter(module %in% c(6,10,9,11,12)),
          label_cell_groups=FALSE,
          show_trajectory_graph=T)
plot_cells(Bmem_cds_subset,
           genes=gene_module_df %>% filter(module %in% c(11,12)),
           label_cell_groups=FALSE,
           show_trajectory_graph=T,label_branch_points = F,label_leaves = F)+
  ggtitle("Monocle cells colored by Gene module expression")+center.title()+
  theme_void()# Export 3x5 Fig.6.H

# Gprofiler analysis
module_11 <- gene_module_df$id[gene_module_df$module=="11"]
module_8 <- gene_module_df$id[gene_module_df$module=="8"]
module_12 <- gene_module_df$id[gene_module_df$module=="12"]
module_5 <- gene_module_df$id[gene_module_df$module=="5"]
module_17 <- gene_module_df$id[gene_module_df$module=="17"]
module_6 <- gene_module_df$id[gene_module_df$module=="6"]

gostres_module_11 <- gost(query = module_11 ,ordered_query = F, significant = T, correction_method = "gSCS")
gostres_module_8 <- gost(query = module_8 ,ordered_query = F, significant = T, correction_method = "gSCS")
gostres_module_12 <- gost(query = module_12 ,ordered_query = F, significant = T, correction_method = "gSCS")
gostres_module_5 <- gost(query = module_5 ,ordered_query = F, significant = T, correction_method = "gSCS")
gostres_module_17 <- gost(query = module_17 ,ordered_query = F, significant = T, correction_method = "gSCS")
gostres_module_6 <- gost(query = module_6 ,ordered_query = F, significant = T, correction_method = "gSCS")

gostres_module_11_df <- gostres_module_11$result
gostres_module_8_df <- gostres_module_8$result
gostres_module_12_df <- gostres_module_12$result
gostres_module_5_df <- gostres_module_5$result
gostres_module_17_df <- gostres_module_17$result
gostres_module_6_df <- gostres_module_6$result

# Making a dot plot out of the g:profiler results:
gostres_module_11_df <- gostres_module_11_df[gostres_module_11_df$p_value < 0.005,c(3,11)]
gostres_module_8_df <- gostres_module_8_df[gostres_module_8_df$p_value < 0.05,c(3,11)]
gostres_module_12_df <- gostres_module_12_df[gostres_module_12_df$p_value < 0.05,c(3,11)]
gostres_module_5_df <- gostres_module_5_df[gostres_module_5_df$p_value < 0.0005,c(3,11)]
gostres_module_17_df <- gostres_module_17_df[gostres_module_17_df$p_value < 0.01,c(3,11)]
gostres_module_6_df <- gostres_module_6_df[gostres_module_6_df$p_value < 0.01,c(3,11)]

gostres_module_11_df <- gostres_module_11_df[gostres_module_11_df$term_name %in% c("actin cytoskeleton","focal adhesion",
                                                                                   "cell-substrate junction","spleen; cells in white pulp[High]",
                                                                                   "Cell-extracellular matrix interactions","antigen processing and presentation of exogenous peptide antigen via MHC class Ib"),]
gostres_module_8_df <- gostres_module_8_df[gostres_module_8_df$term_name %in% c("lymphocyte proliferation","cell activation",
                                                                                "B cell receptor signaling pathway"),]
gostres_module_12_df <- gostres_module_12_df[gostres_module_12_df$term_name %in% c("immune system process","regulation of lymphocyte mediated immunity",
                                                                                   "leukocyte proliferation","secretory granule"),]
gostres_module_5_df <- gostres_module_5_df[gostres_module_5_df$term_name %in% c("leukocyte activation","cellular response to stimulus",
                                                                                "positive regulation of immune system process","cell communication"),]
gostres_module_17_df <- gostres_module_17_df[gostres_module_17_df$term_name %in% c("cell activation","antigen processing and presentation of endogenous peptide antigen via MHC class I via ER pathway",
                                                                                   "MHC class I protein complex","plasma membrane"),]

gostres_module_11_df$p_value <- -log10(gostres_module_11_df$p_value)
gostres_module_8_df$p_value <- -log10(gostres_module_8_df$p_value)
gostres_module_12_df$p_value <- -log10(gostres_module_12_df$p_value)
gostres_module_5_df$p_value <- -log10(gostres_module_5_df$p_value)
gostres_module_17_df$p_value <- -log10(gostres_module_17_df$p_value)

i <- 1
for (i in 1:nrow(gostres_module_11_df)) {
  gostres_module_11_df$term_name[i] <- paste(toupper(substr(gostres_module_11_df$term_name[i], 1, 1)), substr(gostres_module_11_df$term_name[i], 2, nchar(gostres_module_11_df$term_name[i])), sep="")
  i <- i+1
}
i <- 1
for (i in 1:nrow(gostres_module_8_df)) {
  gostres_module_8_df$term_name[i] <- paste(toupper(substr(gostres_module_8_df$term_name[i], 1, 1)), substr(gostres_module_8_df$term_name[i], 2, nchar(gostres_module_8_df$term_name[i])), sep="")
  i <- i+1
}
i <- 1
for (i in 1:nrow(gostres_module_12_df)) {
  gostres_module_12_df$term_name[i] <- paste(toupper(substr(gostres_module_12_df$term_name[i], 1, 1)), substr(gostres_module_12_df$term_name[i], 2, nchar(gostres_module_12_df$term_name[i])), sep="")
  i <- i+1
}
i <- 1
for (i in 1:nrow(gostres_module_5_df)) {
  gostres_module_5_df$term_name[i] <- paste(toupper(substr(gostres_module_5_df$term_name[i], 1, 1)), substr(gostres_module_5_df$term_name[i], 2, nchar(gostres_module_5_df$term_name[i])), sep="")
  i <- i+1
}
i <- 1
for (i in 1:nrow(gostres_module_17_df)) {
  gostres_module_17_df$term_name[i] <- paste(toupper(substr(gostres_module_17_df$term_name[i], 1, 1)), substr(gostres_module_17_df$term_name[i], 2, nchar(gostres_module_17_df$term_name[i])), sep="")
  i <- i+1
}

labs = c("Actin cytoskeleton","Antigen processing and presentation of \n exogenous peptide antigen via MHC class Ib",
         "Cell-extracellular matrix interactions",
         "Cell-substrate junction",
         "Focal adhesion","Spleen; Cells in white pulp[High]" )
# Fig.6.I and Supp.Fig.6.G
ggplot(gostres_module_11_df, mapping = aes(x=p_value, y=term_name)) + geom_point(size=2) + xlim(0,12) + black.axis.text() + xlab(expression(-log[10]~adjusted~italic("P")~value))+
  ylab("")+theme_linedraw()+center.title()+ 
  theme(text = element_text(size = 15),strip.background =element_rect(fill="grey"),strip.text = element_text(colour = 'black'),
        panel.grid.minor = element_blank(),panel.grid.major = element_line(colour = "azure4",size = 0.3))+
  facet_grid(. ~ "GO enrichments Module 11")+ geom_vline(xintercept=c(-log10(0.05),-log10(0.05)), linetype="dashed")+
  scale_x_continuous(limits = c(-0, 12),breaks=c(0,2,4,6,8,10,12))+scale_y_discrete(labels = labs) # Export 3.5x6

ggplot(gostres_module_8_df, mapping = aes(x=p_value, y=term_name)) + geom_point(size=2) + xlim(0,12) + black.axis.text() + xlab(expression(-log[10]~adjusted~italic("P")~value))+
  ylab("")+theme_linedraw()+center.title()+ 
  theme(text = element_text(size = 15),strip.background =element_rect(fill="grey"),strip.text = element_text(colour = 'black'),
        panel.grid.minor = element_blank(),panel.grid.major = element_line(colour = "azure4",size = 0.3))+
  facet_grid(. ~ "GO enrichments Module 8")+ geom_vline(xintercept=c(-log10(0.05),-log10(0.05)), linetype="dashed")+
  scale_x_continuous(limits = c(-0, 12),breaks=c(0,2,4,6,8,10,12)) # Export 3.5x6

ggplot(gostres_module_12_df, mapping = aes(x=p_value, y=term_name)) + geom_point(size=2) + xlim(0,12) + black.axis.text() + xlab(expression(-log[10]~adjusted~italic("P")~value))+
  ylab("")+theme_linedraw()+center.title()+ 
  theme(text = element_text(size = 15),strip.background =element_rect(fill="grey"),strip.text = element_text(colour = 'black'),
        panel.grid.minor = element_blank(),panel.grid.major = element_line(colour = "azure4",size = 0.3))+
  facet_grid(. ~ "GO enrichments Module 12")+ geom_vline(xintercept=c(-log10(0.05),-log10(0.05)), linetype="dashed")+
  scale_x_continuous(limits = c(-0, 12),breaks=c(0,2,4,6,8,10,12)) # Export 3.5x6

ggplot(gostres_module_5_df, mapping = aes(x=p_value, y=term_name)) + geom_point(size=2) + xlim(0,12) + black.axis.text() + xlab(expression(-log[10]~adjusted~italic("P")~value))+
  ylab("")+theme_linedraw()+center.title()+ 
  theme(text = element_text(size = 15),strip.background =element_rect(fill="grey"),strip.text = element_text(colour = 'black'),
        panel.grid.minor = element_blank(),panel.grid.major = element_line(colour = "azure4",size = 0.3))+
  facet_grid(. ~ "GO enrichments Module 5")+ geom_vline(xintercept=c(-log10(0.05),-log10(0.05)), linetype="dashed")+
  scale_x_continuous(limits = c(-0, 12),breaks=c(0,2,4,6,8,10,12)) # Export 3.5x6

labs = c("Antigen processing and presentation of endogenous\npeptide antigen via MHC class I via ER pathway",
         "Cell activation",
         "MHC class I protein complex",
         "Plasma membrane" )

ggplot(gostres_module_17_df, mapping = aes(x=p_value, y=term_name)) + geom_point(size=2) + xlim(0,12) + black.axis.text() + xlab(expression(-log[10]~adjusted~italic("P")~value))+
  ylab("")+theme_linedraw()+center.title()+ 
  theme(text = element_text(size = 15),strip.background =element_rect(fill="grey"),strip.text = element_text(colour = 'black'),
        panel.grid.minor = element_blank(),panel.grid.major = element_line(colour = "azure4",size = 0.3))+
  facet_grid(. ~ "GO enrichments Module 17")+ geom_vline(xintercept=c(-log10(0.05),-log10(0.05)), linetype="dashed")+
  scale_x_continuous(limits = c(-0, 12),breaks=c(0,2,4,6,8,10,12))+scale_y_discrete(labels = labs) # Export 3.5x6


# Supp.Fig.6.D
# Checking B cell subset contributions to Monocle clusters:
# Look at which named clusters make up which Monocle cluster

# Optional: Take only bait positive cells
#Bait.positives <- colnames(Bmem_cds_subset)[colData(Bmem_cds_subset)$bait.positive=="yes"]
#Bmem_cds <- Bmem_cds[,Bait.positives]

Cluster_3_IDs <- names(Bmem_cds@clusters$UMAP$clusters[Bmem_cds@clusters$UMAP$clusters==3])
Cluster_6_IDs <- names(Bmem_cds@clusters$UMAP$clusters[Bmem_cds@clusters$UMAP$clusters==6])
Cluster_1_IDs <- names(Bmem_cds@clusters$UMAP$clusters[Bmem_cds@clusters$UMAP$clusters==1])
Cluster_2_IDs <- names(Bmem_cds@clusters$UMAP$clusters[Bmem_cds@clusters$UMAP$clusters==2])
Cluster_4_IDs <- names(Bmem_cds@clusters$UMAP$clusters[Bmem_cds@clusters$UMAP$clusters==4])
Cluster_5_IDs <- names(Bmem_cds@clusters$UMAP$clusters[Bmem_cds@clusters$UMAP$clusters==5])
Cluster_7_IDs <- names(Bmem_cds@clusters$UMAP$clusters[Bmem_cds@clusters$UMAP$clusters==7])
Cluster_8_IDs <- names(Bmem_cds@clusters$UMAP$clusters[Bmem_cds@clusters$UMAP$clusters==8])

Bmem.subset_Cluster_3 <- subset(Bmem, subset= Full.Row.names %in% Cluster_3_IDs)
Bmem.subset_Cluster_6 <- subset(Bmem, subset= Full.Row.names %in% Cluster_6_IDs)
Bmem.subset_Cluster_1 <- subset(Bmem, subset= Full.Row.names %in% Cluster_1_IDs)
Bmem.subset_Cluster_2 <- subset(Bmem, subset= Full.Row.names %in% Cluster_2_IDs)
Bmem.subset_Cluster_4 <- subset(Bmem, subset= Full.Row.names %in% Cluster_4_IDs)
Bmem.subset_Cluster_5 <- subset(Bmem, subset= Full.Row.names %in% Cluster_5_IDs)
Bmem.subset_Cluster_7 <- subset(Bmem, subset= Full.Row.names %in% Cluster_7_IDs)
Bmem.subset_Cluster_8 <- subset(Bmem, subset= Full.Row.names %in% Cluster_8_IDs)

df_Cluster_3 <- data.frame(table(Bmem.subset_Cluster_3@meta.data$named.clusters))
df_Cluster_6 <- data.frame(table(Bmem.subset_Cluster_6@meta.data$named.clusters))
df_Cluster_1 <- data.frame(table(Bmem.subset_Cluster_1@meta.data$named.clusters))
df_Cluster_2 <- data.frame(table(Bmem.subset_Cluster_2@meta.data$named.clusters))
df_Cluster_4 <- data.frame(table(Bmem.subset_Cluster_4@meta.data$named.clusters))
df_Cluster_5 <- data.frame(table(Bmem.subset_Cluster_5@meta.data$named.clusters))
df_Cluster_7 <- data.frame(table(Bmem.subset_Cluster_7@meta.data$named.clusters))
df_Cluster_8 <- data.frame(table(Bmem.subset_Cluster_8@meta.data$named.clusters))

df_Cluster_3$total <- sum(df_Cluster_3$Freq)
df_Cluster_6$total <- sum(df_Cluster_6$Freq)
df_Cluster_1$total <- sum(df_Cluster_1$Freq)
df_Cluster_2$total <- sum(df_Cluster_2$Freq)
df_Cluster_4$total <- sum(df_Cluster_4$Freq)
df_Cluster_5$total <- sum(df_Cluster_5$Freq)
df_Cluster_7$total <- sum(df_Cluster_7$Freq)
df_Cluster_8$total <- sum(df_Cluster_8$Freq)

df_Cluster_3$percent <- round(100*df_Cluster_3$Freq/df_Cluster_3$total,2)
df_Cluster_6$percent <- round(100*df_Cluster_6$Freq/df_Cluster_6$total,2)
df_Cluster_1$percent <- round(100*df_Cluster_1$Freq/df_Cluster_1$total,2)
df_Cluster_2$percent <- round(100*df_Cluster_2$Freq/df_Cluster_2$total,2)
df_Cluster_4$percent <- round(100*df_Cluster_4$Freq/df_Cluster_4$total,2)
df_Cluster_5$percent <- round(100*df_Cluster_5$Freq/df_Cluster_5$total,2)
df_Cluster_7$percent <- round(100*df_Cluster_7$Freq/df_Cluster_7$total,2)
df_Cluster_8$percent <- round(100*df_Cluster_8$Freq/df_Cluster_8$total,2)

df_Cluster_3$Cluster <- "3"
df_Cluster_6$Cluster <- "6"
df_Cluster_1$Cluster <- "1"
df_Cluster_2$Cluster <- "2"
df_Cluster_4$Cluster <- "4"
df_Cluster_5$Cluster <- "5"
df_Cluster_7$Cluster <- "7"
df_Cluster_8$Cluster <- "8"

df <- rbind(df_Cluster_3,df_Cluster_6,df_Cluster_1,df_Cluster_2,df_Cluster_4,df_Cluster_5,df_Cluster_7,df_Cluster_8)
df$Var1 <- factor(df$Var1, levels = c("Unswitched", "CD27low RM","CD27high RM",
                                      "Activated","Atypical"))
df$Cluster <- factor(df$Cluster, levels = rev(c("1","2","3","4","5","6","7","8")))

ggplot(df, 
       aes(fill=Var1, y=Cluster, x=percent)) + 
  geom_bar(position="fill", stat="identity",colour="black")+theme_classic()+ theme(text = element_text(size = 15),plot.title = element_text(hjust = 0.5),
                                                                                   plot.subtitle = element_text(hjust = 0.5))+
  labs(x = "Subset contribution", y="Monocle Clusters")+scale_x_continuous(labels=scales::percent)+labs(fill = "Subset")+
  scale_fill_manual(values=c("#DBC35E","#228833","#4477AA","#AA3377","#EE6677"))+ theme(text = element_text(size = 17))+ ggtitle("Subset contributions to Monocle clusters")



# Checking Patient contributions to Monocle clusters:
# Look at which Patients make up which Monocle cluster
Cluster_3_IDs <- names(Bmem_cds@clusters$UMAP$clusters[Bmem_cds@clusters$UMAP$clusters==3])
Cluster_6_IDs <- names(Bmem_cds@clusters$UMAP$clusters[Bmem_cds@clusters$UMAP$clusters==6])
Cluster_1_IDs <- names(Bmem_cds@clusters$UMAP$clusters[Bmem_cds@clusters$UMAP$clusters==1])
Cluster_2_IDs <- names(Bmem_cds@clusters$UMAP$clusters[Bmem_cds@clusters$UMAP$clusters==2])
Cluster_4_IDs <- names(Bmem_cds@clusters$UMAP$clusters[Bmem_cds@clusters$UMAP$clusters==4])
Cluster_5_IDs <- names(Bmem_cds@clusters$UMAP$clusters[Bmem_cds@clusters$UMAP$clusters==5])
Cluster_7_IDs <- names(Bmem_cds@clusters$UMAP$clusters[Bmem_cds@clusters$UMAP$clusters==7])
Cluster_8_IDs <- names(Bmem_cds@clusters$UMAP$clusters[Bmem_cds@clusters$UMAP$clusters==8])

Bmem.subset_Cluster_3 <- subset(Bmem, subset= Full.Row.names %in% Cluster_3_IDs)
Bmem.subset_Cluster_6 <- subset(Bmem, subset= Full.Row.names %in% Cluster_6_IDs)
Bmem.subset_Cluster_1 <- subset(Bmem, subset= Full.Row.names %in% Cluster_1_IDs)
Bmem.subset_Cluster_2 <- subset(Bmem, subset= Full.Row.names %in% Cluster_2_IDs)
Bmem.subset_Cluster_4 <- subset(Bmem, subset= Full.Row.names %in% Cluster_4_IDs)
Bmem.subset_Cluster_5 <- subset(Bmem, subset= Full.Row.names %in% Cluster_5_IDs)
Bmem.subset_Cluster_7 <- subset(Bmem, subset= Full.Row.names %in% Cluster_7_IDs)
Bmem.subset_Cluster_8 <- subset(Bmem, subset= Full.Row.names %in% Cluster_8_IDs)

df_Cluster_3 <- data.frame(table(Bmem.subset_Cluster_3@meta.data$Patient))
df_Cluster_6 <- data.frame(table(Bmem.subset_Cluster_6@meta.data$Patient))
df_Cluster_1 <- data.frame(table(Bmem.subset_Cluster_1@meta.data$Patient))
df_Cluster_2 <- data.frame(table(Bmem.subset_Cluster_2@meta.data$Patient))
df_Cluster_4 <- data.frame(table(Bmem.subset_Cluster_4@meta.data$Patient))
df_Cluster_5 <- data.frame(table(Bmem.subset_Cluster_5@meta.data$Patient))
df_Cluster_7 <- data.frame(table(Bmem.subset_Cluster_7@meta.data$Patient))
df_Cluster_8 <- data.frame(table(Bmem.subset_Cluster_8@meta.data$Patient))

df_Cluster_3$total <- sum(df_Cluster_3$Freq)
df_Cluster_6$total <- sum(df_Cluster_6$Freq)
df_Cluster_1$total <- sum(df_Cluster_1$Freq)
df_Cluster_2$total <- sum(df_Cluster_2$Freq)
df_Cluster_4$total <- sum(df_Cluster_4$Freq)
df_Cluster_5$total <- sum(df_Cluster_5$Freq)
df_Cluster_7$total <- sum(df_Cluster_7$Freq)
df_Cluster_8$total <- sum(df_Cluster_8$Freq)

df_Cluster_3$percent <- round(100*df_Cluster_3$Freq/df_Cluster_3$total,2)
df_Cluster_6$percent <- round(100*df_Cluster_6$Freq/df_Cluster_6$total,2)
df_Cluster_1$percent <- round(100*df_Cluster_1$Freq/df_Cluster_1$total,2)
df_Cluster_2$percent <- round(100*df_Cluster_2$Freq/df_Cluster_2$total,2)
df_Cluster_4$percent <- round(100*df_Cluster_4$Freq/df_Cluster_4$total,2)
df_Cluster_5$percent <- round(100*df_Cluster_5$Freq/df_Cluster_5$total,2)
df_Cluster_7$percent <- round(100*df_Cluster_7$Freq/df_Cluster_7$total,2)
df_Cluster_8$percent <- round(100*df_Cluster_8$Freq/df_Cluster_8$total,2)

df_Cluster_3$Cluster <- "3"
df_Cluster_6$Cluster <- "6"
df_Cluster_1$Cluster <- "1"
df_Cluster_2$Cluster <- "2"
df_Cluster_4$Cluster <- "4"
df_Cluster_5$Cluster <- "5"
df_Cluster_7$Cluster <- "7"
df_Cluster_8$Cluster <- "8"

df <- rbind(df_Cluster_3,df_Cluster_6,df_Cluster_1,df_Cluster_2,df_Cluster_4,df_Cluster_5,df_Cluster_7,df_Cluster_8)
df$Var1 <- factor(df$Var1)
df$Cluster <- factor(df$Cluster, levels = rev(c("1","2","3","4","5","6","7","8")))

ggplot(df, 
       aes(fill=Var1, y=Cluster, x=percent)) + 
  geom_bar(position="fill", stat="identity",colour="black")+theme_classic()+ theme(text = element_text(size = 15),plot.title = element_text(hjust = 0.5),
                                                                                   plot.subtitle = element_text(hjust = 0.5))+
  labs(x = "Patient contribution", y="Monocle Clusters")+scale_x_continuous(labels=scales::percent)+labs(fill = "Patient")+ theme(text = element_text(size = 17))+ ggtitle("Patient contributions to Monocle clusters") # Export 5x8



