# Preprocessing of the MBC single cell RNA sequencing data.

# This script is divided into 7 parts. Parts 1-4 are subdivided into sup-parts (1.1 / 1.2 / 1.3). 
# Parts 1-4 are highly repetitive, since they execute the same preprocessing procedure for all 4 datasets that are present for the MBC project.
# Part 5 contains dataset integration and WNN anlaysis
# Part 6 contains Dim reductions and some plotting of UMAPS
# Part 7 produces the heatmap plots, which show the gene segment usages for Spike RBD binders and non binders.

# In order to run the code without errors, some of the file locations need to be adapted.
# All the input files needed are present in the zenodo repository entry of the dataset, where the names of the individual datasets are 
# the same as in the titles of the parts below (Set 1, Set 2, Set 3, Set4).


# Loading packages
library(rgl)
library(ggforce)
library(Seurat)
library(scRepertoire)
library(airr)
library(tidyverse)
library(patchwork)
library(umap)
library(SeuratData)
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
  
} # Highlight choose cells in the FeatureScatter plot
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
    return <- "\u2217"
  }
  if(prop.test.result$p.value<=0.01){
    return <- "\u2217\u2217"
  }
  if(prop.test.result$p.value<=0.001){
    return <- "\u2217\u2217\u2217"
  }
  if(prop.test.result$p.value<=0.0001){
    return <- "\u2217\u2217\u2217\u2217"
  }
  return(return)
} # Used in Part 7 to plot significance asterics for gene segment usage heatmap.
asterics.function2 <- function(prop.test.result) {
  
  if(prop.test.result>0.05){
    return <- "ns"
  }
  if(prop.test.result<=0.05){
    return <- "\u2217"
  }
  if(prop.test.result<=0.01){
    return <- "\u2217\u2217"
  }
  if(prop.test.result<=0.001){
    return <- "\u2217\u2217\u2217"
  }
  if(prop.test.result<=0.0001){
    return <- "\u2217\u2217\u2217\u2217"
  }
  return(return)
} # Used in Part 7 to plot significance asterics for gene segment usage heatmap.
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
} # Outputs some general information about a seurat object
`%notin%` <- Negate(`%in%`) # Creating a useful operator

######################################################################################################################
# Part 1.1: Second Set - Loading the data set, filtering, Hashing demultiplexing, addition of metadata columns
######################################################################################################################
# Loading the data into R
setwd("~/NAS/Jan/Experiments and Data/210801 - Exp 035 - Profiling of COVID vaccinated patient samples/Memory B cell Second Experiment 211127/multi_exp035_MemoryBcellSecondExperiment211127")  

Bmem_dataset <- "./outs/per_sample_outs/multi_exp035_MemBSecondExp211127/count/sample_feature_bc_matrix"
Bmem_dataset <- Read10X(data.dir = Bmem_dataset)

# Seurat object is created
Bmem <- CreateSeuratObject(counts = Bmem_dataset$`Gene Expression`)
rownames(Bmem_dataset$`Antibody Capture`)

# Assays for Baiting constructs, Protein quantification and hashing are created
Baiting_assay <- CreateAssayObject(counts = Bmem_dataset$`Antibody Capture`[c(11:15),])
Protein_assay <- CreateAssayObject(counts = Bmem_dataset$`Antibody Capture`[c(16:20),])
Hashing_assay <- CreateAssayObject(counts = Bmem_dataset$`Antibody Capture`[c(1:10),])

# Now the assays are added to the previously created Seurat object
Bmem[["Baiting"]] <- Baiting_assay
Bmem[["Protein"]] <- Protein_assay
Bmem[["Hashing"]] <- Hashing_assay

# Removal of assay objects
remove(Baiting_assay,Hashing_assay,Protein_assay,Bmem_dataset)

# QC and selecting cells for further analysis
Bmem[["percent.mt"]] <- PercentageFeatureSet(Bmem, pattern = "^MT-")
VlnPlot(Bmem, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
hist(Bmem@meta.data$percent.mt,breaks = c(0:99))
hist(Bmem@meta.data$nFeature_RNA,breaks = c(0:2600))
quantile(Bmem@meta.data$nFeature_RNA,probs = seq(0, 1, 1/20))
Bmem <- subset(Bmem, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 10)

# Now, normalization of the Gene expression data and of the assay data is done.
Bmem <- NormalizeData(Bmem)
Bmem <- NormalizeData(Bmem, assay = "Hashing", normalization.method = "CLR", margin = 2)
Bmem <- NormalizeData(Bmem, assay= "Protein", normalization.method = "CLR", margin = 2)

# Finding variable features (before any further subsetting is done and before dataset integration happens)
Bmem <- FindVariableFeatures(Bmem, assay = "RNA")
Bmem <- FindVariableFeatures(Bmem, assay = "Protein")

# Demultiplexing the HTO data
Bmem <- HTODemux(Bmem, assay = "Hashing", positive.quantile = 0.99)
RidgePlot(Bmem, assay = "Hashing", features = rownames(Bmem[["Hashing"]])[1:10], ncol = 2)

# Filtering out cells, that are not singlets
Bmem <- subset(Bmem, subset = Hashing_classification.global=="Singlet")

# Adding Timepoint and Patient IDs as columns
Bmem@meta.data$Timepoint <- "unknown"
Bmem@meta.data$Timepoint[Bmem@meta.data$hash.ID =="PFCL1-179308-1988-2"|Bmem@meta.data$hash.ID =="PFCL1-199474-1994-2"|
                           Bmem@meta.data$hash.ID =="PFCL1-196878-1691-2"|Bmem@meta.data$hash.ID =="PFCL1-LIM828246-2"|
                           Bmem@meta.data$hash.ID =="PFCL1-154583-1943-2"] <- "6 months after infection"

Bmem@meta.data$Timepoint[Bmem@meta.data$hash.ID =="PFCL1-199474-1994-3" | Bmem@meta.data$hash.ID =="PFCL1-196878-1691-3"|
                           Bmem@meta.data$hash.ID=="PFCL1-LIM828246-3"| Bmem@meta.data$hash.ID =="PFCL1-154583-1943-3"|
                           Bmem@meta.data$hash.ID =="PFCL1-179308-1988-3"] <- "12 months after infection"

Bmem@meta.data$Patient <-  "unknown" 
Bmem@meta.data$Patient[Bmem@meta.data$hash.ID =="PFCL1-179308-1988-2" | Bmem@meta.data$hash.ID =="PFCL1-179308-1988-3"] <- "PFCL1-179308-1988"
Bmem@meta.data$Patient[Bmem@meta.data$hash.ID =="PFCL1-199474-1994-2" | Bmem@meta.data$hash.ID == "PFCL1-199474-1994-3"] <- "PFCL1-199474-1994"
Bmem@meta.data$Patient[Bmem@meta.data$hash.ID =="PFCL1-196878-1691-2" | Bmem@meta.data$hash.ID == "PFCL1-196878-1691-3"] <- "PFCL1-196878-1691"
Bmem@meta.data$Patient[Bmem@meta.data$hash.ID =="PFCL1-LIM828246-2" | Bmem@meta.data$hash.ID == "PFCL1-LIM828246-3"] <- "PFCL1-LIM828246"
Bmem@meta.data$Patient[Bmem@meta.data$hash.ID =="PFCL1-154583-1943-2" | Bmem@meta.data$hash.ID == "PFCL1-154583-1943-3"] <- "PFCL1-154583-1943"

# Adding the info which samples are vaccinated
Bmem@meta.data$Vaccinated <- "no"
Bmem@meta.data$Vaccinated[Bmem@meta.data$hash.ID=="PFCL1-154583-1943-3"|Bmem@meta.data$hash.ID=="PFCL1-179308-1988-3"|
                            Bmem@meta.data$hash.ID=="PFCL1-LIM828246-3"|Bmem@meta.data$hash.ID=="PFCL1-199474-1994-3"] <- "yes"

# Adding the info when 12 months sample was taken with respect to the vaccination
Bmem@meta.data$time_after_vac <- "not vaccinated"
Bmem@meta.data$time_after_vac[Bmem@meta.data$hash.ID=="PFCL1-154583-1943-3"] <- "23d"
Bmem@meta.data$time_after_vac[Bmem@meta.data$hash.ID=="PFCL1-179308-1988-3"] <- "85d"
Bmem@meta.data$time_after_vac[Bmem@meta.data$hash.ID=="PFCL1-LIM828246-3"] <- "87d"
Bmem@meta.data$time_after_vac[Bmem@meta.data$hash.ID=="PFCL1-199474-1994-3"] <- "108d"

# Adding a column that combines timepoint and vaccination status
Bmem@meta.data$sample_group <- "6 months"
Bmem@meta.data$sample_group[Bmem@meta.data$Timepoint=="12 months after infection" & Bmem@meta.data$Vaccinated=="no"] <- "12 months not vaccinated"
Bmem@meta.data$sample_group[Bmem@meta.data$Timepoint=="12 months after infection" & Bmem@meta.data$Vaccinated=="yes"] <- "12 months vaccinated"

######################################################################################################################
# Part 1.2: Second Set - Adding VDJ data, filtering of cells with BCR, Addition of Isotype
######################################################################################################################
# Loading the data
contigs <- read.csv("./outs/per_sample_outs/multi_exp035_MemBSecondExp211127/vdj_b/filtered_contig_annotations.csv")
contig.list <- createHTOContigList(contigs, Bmem, group.by = "hash.ID")

#Fixing the parseBCR function from the scRepertoire version that I use (1.3.5)
my_parseBCR <- function (Con.df, unique_df, data2)
{
  for (y in seq_along(unique_df)) {
    barcode.i <- Con.df$barcode[y]
    location.i <- which(barcode.i == data2$barcode)
    if (length(location.i) == 2) {
      if (is.na(data2[location.i[1], c("IGHct")])) {
        if (is.na(data2[location.i[2], c("IGLct")]) &
            is.na(data2[location.i[2], c("IGKct")])) {
          Con.df[y, heavy_lines] <- data2[location.i[2],
                                          h_lines]
        }
        if (is.na(data2[location.i[1], c("IGKct")])) {
          Con.df[y, light_lines] <- data2[location.i[1],
                                          l_lines]
        }
        else if (!is.na(data2[location.i[1], c("IGKct")])) {
          Con.df[y, light_lines] <- data2[location.i[1],
                                          k_lines]
        }
      }
      else {
        if (is.na(data2[location.i[2], c("IGLct")]) &
            is.na(data2[location.i[2], c("IGKct")])) {
          Con.df[y, heavy_lines] <- data2[location.i[1],
                                          h_lines]
        }
        if (!is.na(data2[location.i[2], c("IGKct")])) {
          Con.df[y, light_lines] <- data2[location.i[2],
                                          k_lines]
          Con.df[y, heavy_lines] <- data2[location.i[1], ###
                                          h_lines] ###
        }
        else {
          Con.df[y, light_lines] <- data2[location.i[2],
                                          l_lines]
          Con.df[y, heavy_lines] <- data2[location.i[1], ###
                                          h_lines] ###
        }
      }
    }
    else if (length(location.i) == 1) {
      chain.i <- data2$chain[location.i]
      if (chain.i == "IGH") {
        Con.df[y, heavy_lines] <- data2[location.i[1],
                                        h_lines]
      }
      else if (chain.i == "IGL") {
        Con.df[y, light_lines] <- data2[location.i[1],
                                        l_lines]
      }
      else {
        Con.df[y, light_lines] <- data2[location.i[1],
                                        k_lines]
      }
    }
  }
  return(Con.df)
}

# Making my adatped version the active version to be used by the package
environment(my_parseBCR) <- asNamespace('scRepertoire')
assignInNamespace("parseBCR", my_parseBCR, ns = "scRepertoire")

# Executing the scRepertoire function
combined <- combineBCR(contig.list,
                       samples = c("PFCL1-179308-198_2", "PFCL1-179308-1988_3", "PFCL1-199474-1994_2", "PFCL1-199474-1994_3",
                                   "PFCL1-196878-1691_2","PFCL1-196878-1691_3","PFCL1-LIM828246_2","PFCL1-LIM828246_3",
                                   "PFCL1-154583-1943_2","PFCL1-154583-1943_3"))

remove(contigs,contig.list)

# The row names of BCRcombined and Bmem need to match, otherwise "combineExpression" fails
head(rownames(Bmem@meta.data))
head(combined$`PFCL1-179308-198_2`[,1])

for (i in seq_along(combined)) {
  combined[[i]] <- stripBarcode(combined[[i]], 
                                column = 1, connector = "_", num_connects = 3)
}

Bmem <- combineExpression(combined, Bmem, cloneCall="aa", cloneTypes=c(Single=1, Small=5, Medium=20, Large=100, Hyperexpanded=500),proportion = FALSE)

# Add info whether cell has BCR info
Bmem@meta.data$BCR.known <- "yes"
Bmem@meta.data$BCR.known[is.na(Bmem@meta.data$CTaa)] <- "no"

# Add info whether cell has both chains
Bmem@meta.data$full.BCR.known <- "yes"
Bmem@meta.data$full.BCR.known[grep("NA",Bmem@meta.data$CTgene)] <- "no"
Bmem@meta.data$full.BCR.known[is.na(Bmem@meta.data$CTaa)]<- "no"

#Removing cells without any BCR info
Bmem <- subset(Bmem, subset = BCR.known =="yes")

# Adding Isotype Info
heavychains <- unlist(strsplit(Bmem@meta.data$CTgene, "[_]"))[seq(1, length(unlist(strsplit(Bmem@meta.data$CTgene, "[_]"))), 2)]
heavychains <- as.data.frame(heavychains)
heavychains$Isotype <- "unknown"
heavychains$Isotype[grep("IGHD",heavychains$heavychains)] <- "IGHD"
heavychains$Isotype[grep("IGHA",heavychains$heavychains)] <- "IGHA"
heavychains$Isotype[grep("IGHM",heavychains$heavychains)] <- "IGHM"
heavychains$Isotype[grep("IGHG",heavychains$heavychains)] <- "IGHG"
Bmem@meta.data$Isotype <- heavychains$Isotype
remove(heavychains)

######################################################################################################################
# Part 1.3: Second Set - Processing of baiting counts (cutoffs, normalization) and visualization
######################################################################################################################
# Pre-processing baiting counts
Baiting_df <- as.data.frame(Bmem@assays$Baiting@counts)
Baiting_df <- t(Baiting_df)
Baiting_df <- as.data.frame(Baiting_df)
colnames(Baiting_df) <- c("wt_Spike", "B.1.351","B.1.617.2","RBD","NegControl")

# Density plots
melted_Baiting_df <- reshape2::melt(Baiting_df)
melted_Baiting_df <- melted_Baiting_df[melted_Baiting_df$value>0,]
ggplot(melted_Baiting_df, aes(x = value, y = variable,fill =variable)) +
  geom_density_ridges(scale = 4)+
  scale_x_continuous(trans='log2',breaks = c(1,5,10,20,40,100,200,500,700,1000,1500,3000,5000))+
  labs(title = "Density plots after negative control subtraction") +
  theme_ridges(font_size = 13, grid = TRUE) + 
  theme(axis.title.y = element_blank(),axis.title.x = element_blank())
remove(melted_Baiting_df)

# Correcting by using the negative control
Baiting_df$cor_wt_Spike <- Baiting_df$wt_Spike-Baiting_df$NegControl
Baiting_df$cor_B.1.351 <- Baiting_df$B.1.351-Baiting_df$NegControl
Baiting_df$cor_B.1.617.2 <- Baiting_df$B.1.617.2-Baiting_df$NegControl
Baiting_df$cor_RBD <- Baiting_df$RBD-Baiting_df$NegControl

# Density plots after negative control subtraction
melted_Baiting_df <- Baiting_df[6:9]
melted_Baiting_df <- reshape2::melt(melted_Baiting_df)
melted_Baiting_df <- melted_Baiting_df[melted_Baiting_df$value>0,]
ggplot(melted_Baiting_df, aes(x = value, y = variable,fill =variable)) +
  geom_density_ridges(scale = 4)+
  scale_x_continuous(trans='log2',breaks = c(1,5,10,20,40,100,200,500,700,1000,1500,3000,5000))+
  labs(title = "") +
  theme_ridges(font_size = 13, grid = TRUE) + 
  theme(axis.title.y = element_blank(),axis.title.x = element_blank())
remove(melted_Baiting_df)

# Classify binders
Baiting_df$cor_wt_Spike[Baiting_df$cor_wt_Spike<1] <- 0
Baiting_df$cor_B.1.351[Baiting_df$cor_B.1.351<1] <- 0
Baiting_df$cor_B.1.617.2[Baiting_df$cor_B.1.617.2<1] <- 0
Baiting_df$cor_RBD[Baiting_df$cor_RBD<1] <- 0

Baiting_df$wt_Spike.classification <- "Negative"
Baiting_df$B.1.351.classification <- "Negative"
Baiting_df$B.1.617.2.classification <- "Negative"
Baiting_df$RBD.classification <- "Negative"

# These cutoffs are chosen after inspection of the density plots
Baiting_df$wt_Spike.classification[Baiting_df$cor_wt_Spike>10] <- "Positive"
Baiting_df$B.1.351.classification[Baiting_df$cor_B.1.351>5] <- "Positive"
Baiting_df$B.1.617.2.classification[Baiting_df$cor_B.1.617.2>5] <- "Positive"
Baiting_df$RBD.classification[Baiting_df$cor_RBD>100] <- "Positive"

# Adding a column that tells, whether a cell is positive for at least one bait and for RBD + full length Spike
Baiting_df$bait.positive <- "no"
Baiting_df$bait.positive[Baiting_df$wt_Spike.classification=="Positive"|
                           Baiting_df$B.1.351.classification=="Positive"|
                           Baiting_df$B.1.617.2.classification=="Positive"|
                           Baiting_df$RBD.classification=="Positive"] <- "yes"

Baiting_df$wt_Spike_RBD.positive <- "no"
Baiting_df$wt_Spike_RBD.positive[Baiting_df$wt_Spike.classification=="Positive"&
                                   Baiting_df$RBD.classification=="Positive"] <- "yes"

# Setting corrected baiting counts to NA if classified as binding negative
Baiting_df$cor_wt_Spike[Baiting_df$wt_Spike.classification=="Negative"] <- NA
Baiting_df$cor_B.1.351[Baiting_df$B.1.351.classification =="Negative"] <- NA
Baiting_df$cor_B.1.617.2[Baiting_df$B.1.617.2.classification =="Negative"] <- NA
Baiting_df$cor_RBD[Baiting_df$RBD.classification=="Negative"] <- NA

# Density plots after cutoff correction
melted_Baiting_df <- Baiting_df[6:9]
melted_Baiting_df <- reshape2::melt(melted_Baiting_df)
melted_Baiting_df <- na.omit(melted_Baiting_df)
ggplot(melted_Baiting_df, aes(x = value, y = variable,fill =variable)) +
  geom_density_ridges(scale = 4)+
  scale_x_continuous(trans='log2',breaks = c(1,5,10,20,40,100,200,500,700,1000,1500,3000,5000))+
  labs(title = "Density plots after cutoff correction") +
  theme_ridges(font_size = 13, grid = TRUE) + 
  theme(axis.title.y = element_blank(),axis.title.x = element_blank())
remove(melted_Baiting_df)

# In order to use the Seurat Normalization function for corrected baiting counts, we need to have the corrected counts as assay object
# For this, the Baiting_df needs to be transposed
Baiting_df <- as.data.frame(t(Baiting_df[,6:9]))
Cor_Baiting_assay <- CreateAssayObject(counts = Baiting_df)
Bmem[["Cor_Baiting"]] <- Cor_Baiting_assay
remove(Cor_Baiting_assay)
Bmem <- NormalizeData(Bmem, assay = "Cor_Baiting", normalization.method = 'CLR', margin = 1)
normalized_across_features <- as.data.frame(t(Bmem@assays$Cor_Baiting@data))
colnames(normalized_across_features) <- c("cor_wt_Spike","cor_B.1.351","cor_B.1.617.2","cor_RBD")
remove(Baiting_df)
Bmem@meta.data$nCount_Cor_Baiting <- NULL
Bmem@meta.data$nFeature_Cor_Baiting <- NULL

# Density plots after normalization
normalized_across_features.melted <- reshape2::melt(normalized_across_features)
normalized_across_features.melted <- na.omit(normalized_across_features.melted)
ggplot(normalized_across_features.melted, aes(x = value, y = variable,fill =variable)) +
  geom_density_ridges(scale = 4)+
  labs(title = "Density plots after normalization") +
  theme_ridges(font_size = 13, grid = TRUE) + 
  theme(axis.title.y = element_blank(),axis.title.x = element_blank())
remove(normalized_across_features.melted)

# For scaling, I cannot use Seurat Scaling function because it gives 0s as result when NAs are present.
# I also need to transform again, then each column is scaled to 0:1 individually using my own scaling function.
normalized_across_features <- myscaling(normalized_across_features)  
normalized_across_features.melted <- reshape2::melt(normalized_across_features)
normalized_across_features.melted <- na.omit(normalized_across_features.melted)
normalized_across_features.melted <- normalized_across_features.melted[normalized_across_features.melted$value>0,]
ggplot(normalized_across_features.melted, aes(x = value, y = variable,fill =variable)) +
  geom_density_ridges(scale = 4)+
  labs(title = "Density plots after scaling") +
  theme_ridges(font_size = 13, grid = TRUE) + 
  theme(axis.title.y = element_blank(),axis.title.x = element_blank())
remove(normalized_across_features.melted)

# These are the final LIBRA scores
# Adding this to the Seurat object
Bmem@meta.data <- merge(Bmem@meta.data,normalized_across_features, by=0)
rownames(Bmem@meta.data) <- Bmem@meta.data$Row.names
remove(normalized_across_features)

# Feature scatters
# NAs are now set to 0
Bmem@meta.data$cor_wt_Spike[is.na(Bmem@meta.data$cor_wt_Spike)] <- 0
Bmem@meta.data$cor_B.1.351[is.na(Bmem@meta.data$cor_B.1.351)] <- 0
Bmem@meta.data$cor_B.1.617.2[is.na(Bmem@meta.data$cor_B.1.617.2)] <- 0
Bmem@meta.data$cor_RBD[is.na(Bmem@meta.data$cor_RBD)] <- 0

FeatureScatter(Bmem,  feature1 = "cor_wt_Spike", feature2 = "cor_RBD", pt.size = 1, group.by = "Isotype")+
  ggtitle("Wt Spike and RBD scores")+xlab("wt-Spike scores")+ylab("RBD scores")+
  
  FeatureScatter(Bmem,  feature1 = "cor_wt_Spike", feature2 = "cor_B.1.351", pt.size = 1, group.by = "Isotype")+
  ggtitle("Wt Spike and Beta variant scores")+xlab("wt-Spike scores")+ylab("Beta variant scores")+
  
  FeatureScatter(Bmem,  feature1 = "cor_wt_Spike", feature2 = "cor_B.1.617.2", pt.size = 1, group.by = "Isotype")+
  ggtitle("Wt Spike and Delta variant scores")+xlab("wt-Spike scores")+ylab("Delta variant scores")+
  
  FeatureScatter(Bmem,  feature1 = "cor_B.1.351", feature2 = "cor_B.1.617.2", pt.size = 1, group.by = "Isotype")+
  ggtitle("Beta and Delta variant scores")+xlab("Beta scores")+ylab("Delta variant scores")

# Addition of bait positive classification columns in metadata
Bmem@meta.data$bait.positive <- "no"
Bmem@meta.data$bait.positive[Bmem@meta.data$cor_wt_Spike>0 | 
                               Bmem@meta.data$cor_B.1.351>0 | 
                               Bmem@meta.data$cor_B.1.617.2 >0 | 
                               Bmem@meta.data$cor_RBD > 0] <- "yes"
Bmem@meta.data$spike_RBD_positive <- "no"
Bmem@meta.data$spike_RBD_positive[Bmem@meta.data$cor_wt_Spike>0 &
                                    Bmem@meta.data$cor_RBD >0] <- "yes"

Bmem@meta.data$B.1.351_positive <- "no"
Bmem@meta.data$B.1.351_positive[Bmem@meta.data$cor_B.1.351>0] <- "yes"

Bmem@meta.data$B.1.617.2_positive <- "no"
Bmem@meta.data$B.1.617.2_positive[Bmem@meta.data$cor_B.1.617.2>0] <- "yes"

Bmem@meta.data$all_positive <- "no"
Bmem@meta.data$all_positive[Bmem@meta.data$cor_wt_Spike>0 & 
                              Bmem@meta.data$cor_B.1.351>0 & 
                              Bmem@meta.data$cor_B.1.617.2 >0 & 
                              Bmem@meta.data$cor_RBD > 0] <- "yes"

# Plots to show Isotypes of bait positives
barplot(table(Bmem@meta.data$Isotype[Bmem@meta.data$bait.positive=="yes"]), main = "Isotypes of bait positive BCRs")
barplot(table(Bmem@meta.data$Isotype[Bmem@meta.data$spike_RBD_positive=="yes"]), main = "Isotypes of RBD + full length Spike positive BCRs")

# Renaming of cells
Bmem <- RenameCells(object = Bmem, add.cell.id = "Dataset_2")
Bmem$Dataset <- "Second"

# Subsetting of cells so that we only include bait positive ones:
# This decision was made based on the observation that RNAseq quality in this dataset is considerably worse than in the other three datasets.
# Cells belonging to this dataset will not be included in any transctiptional profiling.
Bmem <- subset(Bmem, subset = bait.positive =="yes")
Bmem@meta.data$nCount_Cor_Baiting <- NULL
Bmem@meta.data$nFeature_Cor_Baiting <- NULL

######################################################################################################################
# Part 2.1: First Set - Loading the data set, filtering, Hashing demultiplexing, addition of metadata columns
######################################################################################################################
setwd("~/NAS/Jan/Experiments and Data/210801 - Exp 035 - Profiling of COVID vaccinated patient samples/Memory B cell Pilot Data 211010/multi_SpikeMemoryBCells")  

Bmem_pilot_dataset <- "./outs/per_sample_outs/multi_SpikeMemoryBCells/count/sample_feature_bc_matrix"
Bmem_pilot_dataset <- Read10X(data.dir = Bmem_pilot_dataset)

# Seurat object is created
Bmem_pilot <- CreateSeuratObject(counts = Bmem_pilot_dataset$`Gene Expression`)

# Assays for Baiting constructs, Protein quantification and hashing are created
Baiting_assay <- CreateAssayObject(counts = Bmem_pilot_dataset$`Antibody Capture`[c(1:4),])
Protein_assay <- CreateAssayObject(counts = Bmem_pilot_dataset$`Antibody Capture`[c(5:9),])
Hashing_assay <- CreateAssayObject(counts = Bmem_pilot_dataset$`Antibody Capture`[c(10:17),])

# Now the assays are added to the previously created Seurat object
Bmem_pilot[["Baiting"]] <- Baiting_assay
Bmem_pilot[["Protein"]] <- Protein_assay
Bmem_pilot[["Hashing"]] <- Hashing_assay

# Removal of assay objects
remove(Baiting_assay,Hashing_assay,Protein_assay,Bmem_pilot_dataset)

# QC and selecting cells for further analysis
Bmem_pilot[["percent.mt"]] <- PercentageFeatureSet(Bmem_pilot, pattern = "^MT-")
VlnPlot(Bmem_pilot, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
hist(Bmem_pilot@meta.data$percent.mt,breaks = c(0:99))
Bmem_pilot <- subset(Bmem_pilot, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 10)

# Now, normalization of the Gene expression data and of the assay data is done.
Bmem_pilot <- NormalizeData(Bmem_pilot)
Bmem_pilot <- NormalizeData(Bmem_pilot, assay = "Hashing", normalization.method = "CLR", margin = 2)
Bmem_pilot <- NormalizeData(Bmem_pilot, assay= "Protein", normalization.method = "CLR", margin = 2)

# Finding variable features (before any further subsetting is done and before dataset integration happens)
Bmem_pilot <- FindVariableFeatures(Bmem_pilot, assay = "RNA")
Bmem_pilot <- FindVariableFeatures(Bmem_pilot, assay = "Protein")

# Demultiplexing the HTO data
Bmem_pilot <- HTODemux(Bmem_pilot, assay = "Hashing", positive.quantile = 0.99)
RidgePlot(Bmem_pilot, assay = "Hashing", features = rownames(Bmem_pilot[["Hashing"]])[1:7], ncol = 2)

# Filtering out cells, that are not singlets
Bmem_pilot <- subset(Bmem_pilot, subset = Hashing_classification.global=="Singlet")

# Adding Timepoint and Patient IDs as columns
Bmem_pilot@meta.data$Timepoint <- "unknown"
Bmem_pilot@meta.data$Timepoint[Bmem_pilot@meta.data$hash.ID =="PFCL1-LIM-313000-2"|Bmem_pilot@meta.data$hash.ID =="PFCL1-UST-190762-2"|
                                 Bmem_pilot@meta.data$hash.ID =="PFCL1-LIM-674950-2"|Bmem_pilot@meta.data$hash.ID =="PFCL1-LIM-137402-2"] <- "6 months after infection"

Bmem_pilot@meta.data$Timepoint[Bmem_pilot@meta.data$hash.ID =="PFCL1-LIM-313000-3" | Bmem_pilot@meta.data$hash.ID =="PFCL1-UST-190762-3"|
                           Bmem_pilot@meta.data$hash.ID=="PFCL1-LIM-674950-3"| Bmem_pilot@meta.data$hash.ID =="PFCL1-LIM-137402-3"] <- "12 months after infection"

Bmem_pilot@meta.data$Patient <-  "unknown" 
Bmem_pilot@meta.data$Patient[Bmem_pilot@meta.data$hash.ID =="PFCL1-LIM-313000-2" | Bmem_pilot@meta.data$hash.ID =="PFCL1-LIM-313000-3"] <- "PFCL1-LIM-313000"
Bmem_pilot@meta.data$Patient[Bmem_pilot@meta.data$hash.ID =="PFCL1-UST-190762-2" | Bmem_pilot@meta.data$hash.ID == "PFCL1-UST-190762-3"] <- "PFCL1-UST-190762"
Bmem_pilot@meta.data$Patient[Bmem_pilot@meta.data$hash.ID =="PFCL1-LIM-674950-2" | Bmem_pilot@meta.data$hash.ID == "PFCL1-LIM-674950-3"] <- "PFCL1-LIM-674950"
Bmem_pilot@meta.data$Patient[Bmem_pilot@meta.data$hash.ID =="PFCL1-LIM-137402-2" | Bmem_pilot@meta.data$hash.ID == "PFCL1-LIM-137402-3"] <- "PFCL1-LIM-137402"

# Adding the info which samples are vaccinated
Bmem_pilot@meta.data$Vaccinated <- "no"
Bmem_pilot@meta.data$Vaccinated[Bmem_pilot@meta.data$hash.ID=="PFCL1-LIM-313000-3"|Bmem_pilot@meta.data$hash.ID=="PFCL1-LIM-137402-3"] <- "yes"

# Adding the info when 12 months sample was taken with respect to the vaccination
Bmem_pilot@meta.data$time_after_vac <- "not vaccinated"
Bmem_pilot@meta.data$time_after_vac[Bmem_pilot@meta.data$hash.ID=="PFCL1-LIM-313000-3"] <- "15d"
Bmem_pilot@meta.data$time_after_vac[Bmem_pilot@meta.data$hash.ID=="PFCL1-LIM-137402-3"] <- "15d"

# Adding a column that combines timepoint and vaccination status
Bmem_pilot@meta.data$sample_group <- "6 months"
Bmem_pilot@meta.data$sample_group[Bmem_pilot@meta.data$Timepoint=="12 months after infection" & Bmem_pilot@meta.data$Vaccinated=="no"] <- "12 months not vaccinated"
Bmem_pilot@meta.data$sample_group[Bmem_pilot@meta.data$Timepoint=="12 months after infection" & Bmem_pilot@meta.data$Vaccinated=="yes"] <- "12 months vaccinated"

######################################################################################################################
# Part 2.2: First Set - Adding VDJ data, filtering of cells with BCR, Addition of Isotype
######################################################################################################################
# Addition of the VDJ data
contigs <- read.csv("./outs/per_sample_outs/multi_SpikeMemoryBCells/vdj_b/filtered_contig_annotations.csv")
contig.list <- createHTOContigList(contigs, Bmem_pilot, group.by = "hash.ID")

combined_pilot <- combineBCR(contig.list, 
                       samples = c("PFCL1-LIM-313000_2", "PFCL1-LIM-313000_3", "PFCL1-LIM-674950_2", "PFCL1-UST-190762_2", 
                                   "PFCL1-UST-190762_3","PFCL1-LIM-137402_3","PFCL1-LIM-674950_3","PFCL1-LIM-137402_2"))

remove(contigs,contig.list)

# The row names of BCRcombined and Bmem need to match, otherwise "combineExpression" fails
for (i in seq_along(combined_pilot)) {
  combined_pilot[[i]] <- stripBarcode(combined_pilot[[i]], 
                                column = 1, connector = "_", num_connects = 3)
}

Bmem_pilot <- combineExpression(combined_pilot, Bmem_pilot, cloneCall="aa", cloneTypes=c(Single=1, Small=5, Medium=20, Large=100, Hyperexpanded=500),proportion = FALSE)

# Add info whether cell has BCR info
Bmem_pilot@meta.data$BCR.known <- "yes"
Bmem_pilot@meta.data$BCR.known[is.na(Bmem_pilot@meta.data$CTaa)] <- "no"
table(Bmem_pilot@meta.data$BCR.known)

# Add info whether cell has both chains
Bmem_pilot@meta.data$full.BCR.known <- "yes"
Bmem_pilot@meta.data$full.BCR.known[grep("NA",Bmem_pilot@meta.data$CTgene)] <- "no"
Bmem_pilot@meta.data$full.BCR.known[is.na(Bmem_pilot@meta.data$CTaa)]<- "no"

#Removing cells without BCR --> Could be skipped?
Bmem_pilot <- subset(Bmem_pilot, subset = BCR.known =="yes")

# Adding Isotype Info
# Adding a "delete_me" mock CTaa to all cells with no BCR info
Bmem_pilot@meta.data$CTgene[is.na(Bmem_pilot@meta.data$CTgene)] <- "delete_me"
heavychains <- unlist(strsplit(Bmem_pilot@meta.data$CTgene, "[_]"))[seq(1, length(unlist(strsplit(Bmem_pilot@meta.data$CTgene, "[_]"))), 2)]
heavychains <- as.data.frame(heavychains)
heavychains$Isotype <- "unknown"
heavychains$Isotype[grep("IGHD",heavychains$heavychains)] <- "IGHD"
heavychains$Isotype[grep("IGHA",heavychains$heavychains)] <- "IGHA"
heavychains$Isotype[grep("IGHM",heavychains$heavychains)] <- "IGHM"
heavychains$Isotype[grep("IGHG",heavychains$heavychains)] <- "IGHG"
Bmem_pilot@meta.data$Isotype <- heavychains$Isotype
Bmem_pilot@meta.data$CTgene[Bmem_pilot@meta.data$CTgene=="delete_me"] <- NA
remove(heavychains)


######################################################################################################################
# Part 2.3: First Set - Processing of baiting counts (cutoffs, normalization) and visualization
######################################################################################################################
# Pre-processing baiting counts
Pilot_Baiting_df <- as.data.frame(Bmem_pilot@assays$Baiting@counts)
Pilot_Baiting_df <- t(Pilot_Baiting_df)
Pilot_Baiting_df <- as.data.frame(Pilot_Baiting_df)
colnames(Pilot_Baiting_df) <- c("wt_Spike", "B.1.351","B.1.617.2","RBD")

# Density plots
melted_Baiting_df <- reshape2::melt(Pilot_Baiting_df)
melted_Baiting_df <- melted_Baiting_df[melted_Baiting_df$value>0,]
ggplot(melted_Baiting_df, aes(x = value, y = variable,fill =variable)) +
  geom_density_ridges(scale = 4)+
  scale_x_continuous(trans='log2',breaks = c(1,5,10,20,40,100,200,500,700,1000,1500,3000,5000))+
  labs(title = "Density plots after negative control subtraction") +
  theme_ridges(font_size = 13, grid = TRUE) + 
  theme(axis.title.y = element_blank(),axis.title.x = element_blank())
remove(melted_Baiting_df)

# Classify binders
Pilot_Baiting_df$cor_wt_Spike <- Pilot_Baiting_df$wt_Spike
Pilot_Baiting_df$cor_B.1.351 <- Pilot_Baiting_df$B.1.351
Pilot_Baiting_df$cor_B.1.617.2 <- Pilot_Baiting_df$B.1.617.2
Pilot_Baiting_df$cor_RBD <- Pilot_Baiting_df$RBD

Pilot_Baiting_df$wt_Spike.classification <- "Negative"
Pilot_Baiting_df$B.1.351.classification <- "Negative"
Pilot_Baiting_df$B.1.617.2.classification <- "Negative"
Pilot_Baiting_df$RBD.classification <- "Negative"

Pilot_Baiting_df$wt_Spike.classification[Pilot_Baiting_df$cor_wt_Spike>30] <- "Positive"
Pilot_Baiting_df$B.1.351.classification[Pilot_Baiting_df$cor_B.1.351>10] <- "Positive"
Pilot_Baiting_df$B.1.617.2.classification[Pilot_Baiting_df$cor_B.1.617.2>10] <- "Positive"
Pilot_Baiting_df$RBD.classification[Pilot_Baiting_df$cor_RBD>500] <- "Positive"

# Adding a column that tells, whether a cell is positive for at least one bait and for RBD + full length Spike
Pilot_Baiting_df$bait.positive <- "no"
Pilot_Baiting_df$bait.positive[Pilot_Baiting_df$wt_Spike.classification=="Positive"|
                                 Pilot_Baiting_df$B.1.351.classification=="Positive"|
                                 Pilot_Baiting_df$B.1.617.2.classification=="Positive"|
                                 Pilot_Baiting_df$RBD.classification=="Positive"] <- "yes"

Pilot_Baiting_df$wt_Spike_RBD.positive <- "no"
Pilot_Baiting_df$wt_Spike_RBD.positive[Pilot_Baiting_df$wt_Spike.classification=="Positive"&
                                         Pilot_Baiting_df$RBD.classification=="Positive"] <- "yes"

# Setting corrected baiting counts to NA if classified as binding negative
Pilot_Baiting_df$cor_wt_Spike[Pilot_Baiting_df$wt_Spike.classification=="Negative"] <- NA
Pilot_Baiting_df$cor_B.1.351[Pilot_Baiting_df$B.1.351.classification =="Negative"] <- NA
Pilot_Baiting_df$cor_B.1.617.2[Pilot_Baiting_df$B.1.617.2.classification =="Negative"] <- NA
Pilot_Baiting_df$cor_RBD[Pilot_Baiting_df$RBD.classification=="Negative"] <- NA

# Density plots after cutoff correction
melted_Baiting_df <- Pilot_Baiting_df[5:8]
melted_Baiting_df <- reshape2::melt(melted_Baiting_df)
melted_Baiting_df <- na.omit(melted_Baiting_df)
ggplot(melted_Baiting_df, aes(x = value, y = variable,fill =variable)) +
  geom_density_ridges(scale = 4)+
  scale_x_continuous(trans='log2',breaks = c(1,5,10,20,40,100,200,500,700,1000,1500,3000,5000))+
  labs(title = "Density plots after negative control subtraction") +
  theme_ridges(font_size = 13, grid = TRUE) + 
  theme(axis.title.y = element_blank(),axis.title.x = element_blank())
remove(melted_Baiting_df)

# In order to use the Seurat Normalization function, we need to have the corrected counts as assay object
# For this, the Baiting_df needs to be transposed
Pilot_Baiting_df <- as.data.frame(t(Pilot_Baiting_df[,5:8]))
Cor_Baiting_assay <- CreateAssayObject(counts = Pilot_Baiting_df)
Bmem_pilot[["Cor_Baiting"]] <- Cor_Baiting_assay
remove(Cor_Baiting_assay)
Bmem_pilot <- NormalizeData(Bmem_pilot, assay = "Cor_Baiting", normalization.method = 'CLR', margin = 1)
normalized_across_cells <- as.data.frame(t(Bmem_pilot@assays$Cor_Baiting@data))
colnames(normalized_across_cells) <- c("cor_wt_Spike","cor_B.1.351","cor_B.1.617.2","cor_RBD")
remove(Pilot_Baiting_df)
Bmem_pilot@meta.data$nCount_Cor_Baiting <- NULL
Bmem_pilot@meta.data$nFeature_Cor_Baiting <- NULL

# For scaling, I cannot use Seurat Scaling function because it gives 0s as result when NAs are present.
# I also need to transform again, then each column is scaled to 0:1 individually
normalized_across_cells <- myscaling(normalized_across_cells)  

# These are the final LIBRA scores
# Adding this to the Seurat object
Bmem_pilot@meta.data <- merge(Bmem_pilot@meta.data,normalized_across_cells, by=0)
rownames(Bmem_pilot@meta.data) <- Bmem_pilot@meta.data$Row.names
remove(normalized_across_cells)

# Feature scatters
# NAs are now set to 0
Bmem_pilot@meta.data$cor_wt_Spike[is.na(Bmem_pilot@meta.data$cor_wt_Spike)] <- 0
Bmem_pilot@meta.data$cor_B.1.351[is.na(Bmem_pilot@meta.data$cor_B.1.351)] <- 0
Bmem_pilot@meta.data$cor_B.1.617.2[is.na(Bmem_pilot@meta.data$cor_B.1.617.2)] <- 0
Bmem_pilot@meta.data$cor_RBD[is.na(Bmem_pilot@meta.data$cor_RBD)] <- 0

FeatureScatter(Bmem_pilot,  feature1 = "cor_wt_Spike", feature2 = "cor_RBD", pt.size = 1, group.by = "Isotype")+
  ggtitle("Wt Spike and RBD scores")+xlab("wt-Spike scores")+ylab("RBD scores")+
  
  FeatureScatter(Bmem_pilot,  feature1 = "cor_wt_Spike", feature2 = "cor_B.1.351", pt.size = 1, group.by = "Isotype")+
  ggtitle("Wt Spike and Beta variant scores")+xlab("wt-Spike scores")+ylab("Beta variant scores")+
  
  FeatureScatter(Bmem_pilot,  feature1 = "cor_wt_Spike", feature2 = "cor_B.1.617.2", pt.size = 1, group.by = "Isotype")+
  ggtitle("Wt Spike and Delta variant scores")+xlab("wt-Spike scores")+ylab("Delta variant scores")+
  
  FeatureScatter(Bmem_pilot,  feature1 = "cor_B.1.351", feature2 = "cor_B.1.617.2", pt.size = 1, group.by = "Isotype")+
  ggtitle("Beta and Delta variant scores")+xlab("Beta scores")+ylab("Delta variant scores")


# Addition of bait positive classification columns in metadata
Bmem_pilot@meta.data$bait.positive <- "no"
Bmem_pilot@meta.data$bait.positive[Bmem_pilot@meta.data$cor_wt_Spike>0 | 
                               Bmem_pilot@meta.data$cor_B.1.351>0 | 
                               Bmem_pilot@meta.data$cor_B.1.617.2 >0 | 
                               Bmem_pilot@meta.data$cor_RBD > 0] <- "yes"
Bmem_pilot@meta.data$spike_RBD_positive <- "no"
Bmem_pilot@meta.data$spike_RBD_positive[Bmem_pilot@meta.data$cor_wt_Spike>0 &
                                    Bmem_pilot@meta.data$cor_RBD >0] <- "yes"

Bmem_pilot@meta.data$B.1.351_positive <- "no"
Bmem_pilot@meta.data$B.1.351_positive[Bmem_pilot@meta.data$cor_B.1.351>0] <- "yes"

Bmem_pilot@meta.data$B.1.617.2_positive <- "no"
Bmem_pilot@meta.data$B.1.617.2_positive[Bmem_pilot@meta.data$cor_B.1.617.2>0] <- "yes"

Bmem_pilot@meta.data$all_positive <- "no"
Bmem_pilot@meta.data$all_positive[Bmem_pilot@meta.data$cor_wt_Spike>0 & 
                              Bmem_pilot@meta.data$cor_B.1.351>0 & 
                              Bmem_pilot@meta.data$cor_B.1.617.2 >0 & 
                              Bmem_pilot@meta.data$cor_RBD > 0] <- "yes"

# Plots to show Isotypes of bait positives
barplot(table(Bmem_pilot@meta.data$Isotype[Bmem_pilot@meta.data$bait.positive=="yes"]), main = "Isotypes of bait positive BCRs")
barplot(table(Bmem_pilot@meta.data$Isotype[Bmem_pilot@meta.data$spike_RBD_positive=="yes"]), main = "Isotypes of RBD + full length Spike positive BCRs")

# Renaming of cells
Bmem_pilot <- RenameCells(object = Bmem_pilot, add.cell.id = "Dataset_1")
Bmem_pilot$Dataset <- "First"

######################################################################################################################
# Part 3.1: Third Set - Loading the data set, filtering, Hashing demultiplexing, addition of metadata columns
######################################################################################################################
# Loading the data into R
setwd("~/NAS/Jan/Experiments and Data/210801 - Exp 035 - Profiling of COVID vaccinated patient samples/Memory B cell Third Experiment 220214/multi_exp035_3")  

Bmem_third_dataset <- "./outs/per_sample_outs/multi_exp035_3/count/sample_feature_bc_matrix"
Bmem_third_dataset <- Read10X(data.dir = Bmem_third_dataset)

# Seurat object is created
Bmem_third <- CreateSeuratObject(counts = Bmem_third_dataset$`Gene Expression`)
rownames(Bmem_third_dataset$`Antibody Capture`)

# Assays for Baiting constructs, Protein quantification and hashing are created
Baiting_assay <- CreateAssayObject(counts = Bmem_third_dataset$`Antibody Capture`[c(5:9),])
Protein_assay <- CreateAssayObject(counts = Bmem_third_dataset$`Antibody Capture`[c(10:14),])
Hashing_assay <- CreateAssayObject(counts = Bmem_third_dataset$`Antibody Capture`[c(1:4),])

# Now the assays are added to the previously created Seurat object
Bmem_third[["Baiting"]] <- Baiting_assay
Bmem_third[["Protein"]] <- Protein_assay
Bmem_third[["Hashing"]] <- Hashing_assay

# Removal of assay objects
remove(Baiting_assay,Hashing_assay,Protein_assay,Bmem_third_dataset)

# QC and selecting cells for further analysis
Bmem_third[["percent.mt"]] <- PercentageFeatureSet(Bmem_third, pattern = "^MT-")
VlnPlot(Bmem_third, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
hist(Bmem_third@meta.data$percent.mt,breaks = c(0:99))
hist(Bmem_third@meta.data$nFeature_RNA,breaks = c(0:7000))
quantile(Bmem_third@meta.data$nFeature_RNA,probs = seq(0, 1, 1/20))
Bmem_third <- subset(Bmem_third, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 10) # I am increasing the normal nFeature cutoff here.

# Now, normalization of the Gene expression data and of the assay data is done.
Bmem_third <- NormalizeData(Bmem_third)
Bmem_third <- NormalizeData(Bmem_third, assay = "Hashing", normalization.method = "CLR", margin = 2)
Bmem_third <- NormalizeData(Bmem_third, assay= "Protein", normalization.method = "CLR", margin = 2)

# Finding variable features (before any further subsetting is done and before dataset integration happens)
Bmem_third <- FindVariableFeatures(Bmem_third, assay = "RNA")
Bmem_third <- FindVariableFeatures(Bmem_third, assay = "Protein")

# Demultiplexing the HTO data
Bmem_third <- HTODemux(Bmem_third, assay = "Hashing", positive.quantile = 0.99)
RidgePlot(Bmem_third, assay = "Hashing", features = rownames(Bmem_third[["Hashing"]])[1:4], ncol = 2)

# Filtering out cells, that are not singlets
Bmem_third <- subset(Bmem_third, subset = Hashing_classification.global=="Singlet")

# Adding Timepoint and Patient IDs as columns
Bmem_third@meta.data$Timepoint <- "12 months after infection"

Bmem_third@meta.data$Patient <-  "unknown" 
Bmem_third@meta.data$Patient[Bmem_third@meta.data$hash.ID =="PFCL1-LIM-137402-3"] <- "PFCL1-LIM-137402"
Bmem_third@meta.data$Patient[Bmem_third@meta.data$hash.ID =="PFCL1-196878-1691-3"] <- "PFCL1-196878-1691"
Bmem_third@meta.data$Patient[Bmem_third@meta.data$hash.ID =="PFCL1-179308-1988-3"] <- "PFCL1-179308-1988"
Bmem_third@meta.data$Patient[Bmem_third@meta.data$hash.ID =="PFCL1-199474-1994-3"] <- "PFCL1-199474-1994"

# Adding the info which samples are vaccinated
Bmem_third@meta.data$Vaccinated <- "no"
Bmem_third@meta.data$Vaccinated[Bmem_third@meta.data$hash.ID=="PFCL1-LIM-137402-3"|Bmem_third@meta.data$hash.ID=="PFCL1-179308-1988-3"|
                                  Bmem_third@meta.data$hash.ID=="PFCL1-199474-1994-3"] <- "yes"

# Adding the info when 12 months sample was taken with respect to the vaccination
Bmem_third@meta.data$time_after_vac <- "not vaccinated"
Bmem_third@meta.data$time_after_vac[Bmem_third@meta.data$hash.ID=="PFCL1-LIM-137402-3"] <- "15d"
Bmem_third@meta.data$time_after_vac[Bmem_third@meta.data$hash.ID=="PFCL1-179308-1988-3"] <- "85d"
Bmem_third@meta.data$time_after_vac[Bmem_third@meta.data$hash.ID=="PFCL1-199474-1994-3"] <- "108d"

# Adding a column that combines timepoint and vaccination status
Bmem_third@meta.data$sample_group <- "12 months vaccinated"
Bmem_third@meta.data$sample_group[Bmem_third@meta.data$Timepoint=="12 months after infection" & Bmem_third@meta.data$Vaccinated=="no"] <- "12 months not vaccinated"


######################################################################################################################
# Part 3.2: Third Set - Adding VDJ data, filtering of cells with BCR, Addition of Isotype
######################################################################################################################
# Loading the data
contigs <- read.csv("./outs/per_sample_outs/multi_exp035_3/vdj_b/filtered_contig_annotations.csv")
contig.list <- createHTOContigList(contigs, Bmem_third, group.by = "hash.ID")
summary(contig.list)

combined_third <- combineBCR(contig.list,
                       samples = c("PFCL1-LIM-137402_3", "PFCL1-199474-1994_3", "PFCL1-179308-1988_3", "PFCL1-196878-1691_3"))

remove(contigs,contig.list)

# The row names of BCRcombined and Bmem need to match, otherwise "combineExpression" fails
head(rownames(Bmem_third@meta.data))
head(combined_third$`PFCL1-LIM-137402_3`$barcode)

for (i in seq_along(combined_third)) {
  combined_third[[i]] <- stripBarcode(combined_third[[i]], 
                                column = 1, connector = "_", num_connects = 3)
}

Bmem_third <- combineExpression(combined_third, Bmem_third, cloneCall="aa", cloneTypes=c(Single=1, Small=5, Medium=20, Large=100, Hyperexpanded=500),proportion = FALSE)

# Add info whether cell has BCR info
Bmem_third@meta.data$BCR.known <- "yes"
Bmem_third@meta.data$BCR.known[is.na(Bmem_third@meta.data$CTaa)] <- "no"
table(Bmem_third@meta.data$BCR.known)

# Add info whether cell has both chains
Bmem_third@meta.data$full.BCR.known <- "yes"
Bmem_third@meta.data$full.BCR.known[grep("NA",Bmem_third@meta.data$CTgene)] <- "no"
Bmem_third@meta.data$full.BCR.known[is.na(Bmem_third@meta.data$CTaa)]<- "no"

#Removing cells without BCR --> Could be skipped?
Bmem_third <- subset(Bmem_third, subset = BCR.known =="yes")

# Adding Isotype Info
# Adding a "delete_me" mock CTaa to all cells with no BCR info
Bmem_third@meta.data$CTgene[is.na(Bmem_third@meta.data$CTgene)] <- "delete_me"
heavychains <- unlist(strsplit(Bmem_third@meta.data$CTgene, "[_]"))[seq(1, length(unlist(strsplit(Bmem_third@meta.data$CTgene, "[_]"))), 2)]
heavychains <- as.data.frame(heavychains)
heavychains$Isotype <- "unknown"
heavychains$Isotype[grep("IGHD",heavychains$heavychains)] <- "IGHD"
heavychains$Isotype[grep("IGHA",heavychains$heavychains)] <- "IGHA"
heavychains$Isotype[grep("IGHM",heavychains$heavychains)] <- "IGHM"
heavychains$Isotype[grep("IGHG",heavychains$heavychains)] <- "IGHG"
Bmem_third@meta.data$Isotype <- heavychains$Isotype
Bmem_third@meta.data$CTgene[Bmem_third@meta.data$CTgene=="delete_me"] <- NA
remove(heavychains)


######################################################################################################################
# Part 3.3: Third Set - Processing of baiting counts (cutoffs, normalization) and visualization
######################################################################################################################
# Pre-processing baiting counts
Baiting_df <- as.data.frame(Bmem_third@assays$Baiting@counts)
Baiting_df <- t(Baiting_df)
Baiting_df <- as.data.frame(Baiting_df)
colnames(Baiting_df) <- c("wt_Spike", "B.1.351","B.1.617.2","RBD","NegControl")

# Density plots
melted_Baiting_df <- reshape2::melt(Baiting_df)
melted_Baiting_df <- melted_Baiting_df[melted_Baiting_df$value>0,]
ggplot(melted_Baiting_df, aes(x = value, y = variable,fill =variable)) +
  geom_density_ridges(scale = 4)+
  scale_x_continuous(trans='log2',breaks = c(1,5,10,20,40,100,200,500,700,1000,1500,3000,5000))+
  labs(title = "Density plots after negative control subtraction") +
  theme_ridges(font_size = 13, grid = TRUE) + 
  theme(axis.title.y = element_blank(),axis.title.x = element_blank())
remove(melted_Baiting_df)

# Correcting by using the negative control
Baiting_df$cor_wt_Spike <- Baiting_df$wt_Spike-Baiting_df$NegControl
Baiting_df$cor_B.1.351 <- Baiting_df$B.1.351-Baiting_df$NegControl
Baiting_df$cor_B.1.617.2 <- Baiting_df$B.1.617.2-Baiting_df$NegControl
Baiting_df$cor_RBD <- Baiting_df$RBD-Baiting_df$NegControl

# Density plots after negative control subtraction
melted_Baiting_df <- Baiting_df[6:9]
melted_Baiting_df <- reshape2::melt(melted_Baiting_df)
melted_Baiting_df <- melted_Baiting_df[melted_Baiting_df$value>0,]
ggplot(melted_Baiting_df, aes(x = value, y = variable,fill =variable)) +
  geom_density_ridges(scale = 4)+
  scale_x_continuous(trans='log2',breaks = c(1,5,10,20,40,100,200,500,700,1000,1500,3000,5000))+
  labs(title = "Density plots after negative control subtraction") +
  theme_ridges(font_size = 13, grid = TRUE) + 
  theme(axis.title.y = element_blank(),axis.title.x = element_blank())
remove(melted_Baiting_df)

# Classify binders
Baiting_df$cor_wt_Spike[Baiting_df$cor_wt_Spike<1] <- 0
Baiting_df$cor_B.1.351[Baiting_df$cor_B.1.351<1] <- 0
Baiting_df$cor_B.1.617.2[Baiting_df$cor_B.1.617.2<1] <- 0
Baiting_df$cor_RBD[Baiting_df$cor_RBD<1] <- 0

Baiting_df$wt_Spike.classification <- "Negative"
Baiting_df$B.1.351.classification <- "Negative"
Baiting_df$B.1.617.2.classification <- "Negative"
Baiting_df$RBD.classification <- "Negative"

Baiting_df$wt_Spike.classification[Baiting_df$cor_wt_Spike>40] <- "Positive"
Baiting_df$B.1.351.classification[Baiting_df$cor_B.1.351>10] <- "Positive"
Baiting_df$B.1.617.2.classification[Baiting_df$cor_B.1.617.2>10] <- "Positive"
Baiting_df$RBD.classification[Baiting_df$cor_RBD>200] <- "Positive"

# Adding a column that tells, whether a cell is positive for at least one bait and for RBD + full length Spike
Baiting_df$bait.positive <- "no"
Baiting_df$bait.positive[Baiting_df$wt_Spike.classification=="Positive"|
                           Baiting_df$B.1.351.classification=="Positive"|
                           Baiting_df$B.1.617.2.classification=="Positive"|
                           Baiting_df$RBD.classification=="Positive"] <- "yes"

Baiting_df$wt_Spike_RBD.positive <- "no"
Baiting_df$wt_Spike_RBD.positive[Baiting_df$wt_Spike.classification=="Positive"&
                                   Baiting_df$RBD.classification=="Positive"] <- "yes"

# Setting corrected baiting counts to NA if classified as binding negative
Baiting_df$cor_wt_Spike[Baiting_df$wt_Spike.classification=="Negative"] <- NA
Baiting_df$cor_B.1.351[Baiting_df$B.1.351.classification =="Negative"] <- NA
Baiting_df$cor_B.1.617.2[Baiting_df$B.1.617.2.classification =="Negative"] <- NA
Baiting_df$cor_RBD[Baiting_df$RBD.classification=="Negative"] <- NA

# Density plots after cutoff correction
melted_Baiting_df <- Baiting_df[6:9]
melted_Baiting_df <- reshape2::melt(melted_Baiting_df)
melted_Baiting_df <- na.omit(melted_Baiting_df)
ggplot(melted_Baiting_df, aes(x = value, y = variable,fill =variable)) +
  geom_density_ridges(scale = 4)+
  scale_x_continuous(trans='log2',breaks = c(1,5,10,20,40,100,200,500,700,1000,1500,3000,5000))+
  labs(title = "Density plots after negative control subtraction") +
  theme_ridges(font_size = 13, grid = TRUE) + 
  theme(axis.title.y = element_blank(),axis.title.x = element_blank())
remove(melted_Baiting_df)


# In order to use the Seurat Normalization function, we need to have the corrected counts as assay object
# For this, the Baiting_df needs to be transposed
Baiting_df <- as.data.frame(t(Baiting_df[,6:9]))
Cor_Baiting_assay <- CreateAssayObject(counts = Baiting_df)
Bmem_third[["Cor_Baiting"]] <- Cor_Baiting_assay
remove(Cor_Baiting_assay)
Bmem_third <- NormalizeData(Bmem_third, assay = "Cor_Baiting", normalization.method = 'CLR', margin = 1) # Decide here: Normalize across features or cells.
normalized_across_cells <- as.data.frame(t(Bmem_third@assays$Cor_Baiting@data))
colnames(normalized_across_cells) <- c("cor_wt_Spike","cor_B.1.351","cor_B.1.617.2","cor_RBD")
remove(Baiting_df)
Bmem_third@meta.data$nCount_Cor_Baiting <- NULL
Bmem_third@meta.data$nFeature_Cor_Baiting <- NULL

# For scaling, I cannot use Seurat Scaling function because it gives 0s as result when NAs are present.
# I also need to transform again, then each column is scaled to 0:1 individually
normalized_across_cells <- myscaling(normalized_across_cells)  

# These are the final LIBRA scores
# Adding this to the Seurat object
Bmem_third@meta.data <- merge(Bmem_third@meta.data,normalized_across_cells, by=0)
rownames(Bmem_third@meta.data) <- Bmem_third@meta.data$Row.names
remove(normalized_across_cells)

# Feature scatters
# NAs are now set to 0
Bmem_third@meta.data$cor_wt_Spike[is.na(Bmem_third@meta.data$cor_wt_Spike)] <- 0
Bmem_third@meta.data$cor_B.1.351[is.na(Bmem_third@meta.data$cor_B.1.351)] <- 0
Bmem_third@meta.data$cor_B.1.617.2[is.na(Bmem_third@meta.data$cor_B.1.617.2)] <- 0
Bmem_third@meta.data$cor_RBD[is.na(Bmem_third@meta.data$cor_RBD)] <- 0

FeatureScatter(Bmem_third,  feature1 = "cor_wt_Spike", feature2 = "cor_RBD", pt.size = 1, group.by = "Isotype")+
  ggtitle("Wt Spike and RBD scores")+xlab("wt-Spike scores")+ylab("RBD scores")+
  
  FeatureScatter(Bmem_third,  feature1 = "cor_wt_Spike", feature2 = "cor_B.1.351", pt.size = 1, group.by = "Isotype")+
  ggtitle("Wt Spike and Beta variant scores")+xlab("wt-Spike scores")+ylab("Beta variant scores")+
  
  FeatureScatter(Bmem_third,  feature1 = "cor_wt_Spike", feature2 = "cor_B.1.617.2", pt.size = 1, group.by = "Isotype")+
  ggtitle("Wt Spike and Delta variant scores")+xlab("wt-Spike scores")+ylab("Delta variant scores")+
  
  FeatureScatter(Bmem_third,  feature1 = "cor_B.1.351", feature2 = "cor_B.1.617.2", pt.size = 1, group.by = "Isotype")+
  ggtitle("Beta and Delta variant scores")+xlab("Beta scores")+ylab("Delta variant scores")

# Addition of bait positive classification columns in metadata
Bmem_third@meta.data$bait.positive <- "no"
Bmem_third@meta.data$bait.positive[Bmem_third@meta.data$cor_wt_Spike>0 | 
                                     Bmem_third@meta.data$cor_B.1.351>0 | 
                                     Bmem_third@meta.data$cor_B.1.617.2 >0 | 
                                     Bmem_third@meta.data$cor_RBD > 0] <- "yes"
Bmem_third@meta.data$spike_RBD_positive <- "no"
Bmem_third@meta.data$spike_RBD_positive[Bmem_third@meta.data$cor_wt_Spike>0 &
                                          Bmem_third@meta.data$cor_RBD >0] <- "yes"

Bmem_third@meta.data$B.1.351_positive <- "no"
Bmem_third@meta.data$B.1.351_positive[Bmem_third@meta.data$cor_B.1.351>0] <- "yes"

Bmem_third@meta.data$B.1.617.2_positive <- "no"
Bmem_third@meta.data$B.1.617.2_positive[Bmem_third@meta.data$cor_B.1.617.2>0] <- "yes"

Bmem_third@meta.data$all_positive <- "no"
Bmem_third@meta.data$all_positive[Bmem_third@meta.data$cor_wt_Spike>0 & 
                                    Bmem_third@meta.data$cor_B.1.351>0 & 
                                    Bmem_third@meta.data$cor_B.1.617.2 >0 & 
                                    Bmem_third@meta.data$cor_RBD > 0] <- "yes"

# Plots to show Isotypes of bait positives
barplot(table(Bmem_third@meta.data$Isotype[Bmem_third@meta.data$bait.positive=="yes"]), main = "Isotypes of bait positive BCRs")
barplot(table(Bmem_third@meta.data$Isotype[Bmem_third@meta.data$spike_RBD_positive=="yes"]), main = "Isotypes of RBD + full length Spike positive BCRs")

# Renaming of cells
Bmem_third <- RenameCells(object = Bmem_third, add.cell.id = "Dataset_3")
Bmem_third$Dataset <- "Third"


######################################################################################################################
# Part 4.1: Fourth Set - Loading the data set, filtering, Hashing demultiplexing, addition of metadata columns
######################################################################################################################
# Loading the data into R
setwd("~/NAS/Jan/Experiments and Data/210801 - Exp 035 - Profiling of COVID vaccinated patient samples/Memory B cell Fourth Experiment 220323/multi_exp035_4")  

Bmem_fourth_dataset <- "./outs/per_sample_outs/multi_exp035_4/count/sample_feature_bc_matrix"
Bmem_fourth_dataset <- Read10X(data.dir = Bmem_fourth_dataset)

# Seurat object is created
Bmem_fourth <- CreateSeuratObject(counts = Bmem_fourth_dataset$`Gene Expression`)
rownames(Bmem_fourth_dataset$`Antibody Capture`)

# Assays for Baiting constructs, Protein quantification and hashing are created
Baiting_assay <- CreateAssayObject(counts = Bmem_fourth_dataset$`Antibody Capture`[c(8:12),])
Protein_assay <- CreateAssayObject(counts = Bmem_fourth_dataset$`Antibody Capture`[c(13:18),])
Hashing_assay <- CreateAssayObject(counts = Bmem_fourth_dataset$`Antibody Capture`[c(1:7),])

# Now the assays are added to the previously created Seurat object
Bmem_fourth[["Baiting"]] <- Baiting_assay
Bmem_fourth[["Protein"]] <- Protein_assay
Bmem_fourth[["Hashing"]] <- Hashing_assay

# Removal of assay objects
remove(Baiting_assay,Hashing_assay,Protein_assay,Bmem_fourth_dataset)

# QC and selecting cells for further analysis
Bmem_fourth[["percent.mt"]] <- PercentageFeatureSet(Bmem_fourth, pattern = "^MT-")
VlnPlot(Bmem_fourth, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
hist(Bmem_fourth@meta.data$percent.mt,breaks = c(0:99))
hist(Bmem_fourth@meta.data$nFeature_RNA,breaks = c(0:7000))
quantile(Bmem_fourth@meta.data$nFeature_RNA,probs = seq(0, 1, 1/20))
Bmem_fourth <- subset(Bmem_fourth, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 10) # I am increasing the normal nFeature cutoff here.

# Now, normalization of the Gene expression data and of the assay data is done.
DefaultAssay(Bmem_fourth)
Bmem_fourth <- NormalizeData(Bmem_fourth)
Bmem_fourth <- NormalizeData(Bmem_fourth, assay = "Hashing", normalization.method = "CLR", margin = 2)
Bmem_fourth <- NormalizeData(Bmem_fourth, assay= "Protein", normalization.method = "CLR", margin = 2)

# Finding variable features (before any further subsetting is done and before dataset integration happens)
Bmem_fourth <- FindVariableFeatures(Bmem_fourth, assay = "RNA")
Bmem_fourth <- FindVariableFeatures(Bmem_fourth, assay = "Protein")

# Demultiplexing the HTO data
Bmem_fourth <- HTODemux(Bmem_fourth, assay = "Hashing", positive.quantile = 0.99)
RidgePlot(Bmem_fourth, assay = "Hashing", features = rownames(Bmem_fourth[["Hashing"]])[1:7], ncol = 2)

# Filtering out cells, that are not singlets
Bmem_fourth <- subset(Bmem_fourth, subset = Hashing_classification.global=="Singlet")
nrow(Bmem_fourth@meta.data)

# Adding Timepoint and Patient IDs as columns
Bmem_fourth@meta.data$Timepoint <- "12 months after infection"
Bmem_fourth@meta.data$Timepoint[Bmem_fourth@meta.data$hash.ID=="PFCL1-199474-1994-2"|
                                  Bmem_fourth@meta.data$hash.ID=="PFCL1-179308-1988-2"|
                                  Bmem_fourth@meta.data$hash.ID=="PFCL1-LIM-313000-2"|
                                  Bmem_fourth@meta.data$hash.ID=="PFCL1-LIM828246-2"] <- "6 months after infection"

Bmem_fourth@meta.data$Patient <-  "unknown" 
Bmem_fourth@meta.data$Patient[Bmem_fourth@meta.data$hash.ID =="PFCL1-199474-1994-2"] <- "PFCL1-199474-1994"
Bmem_fourth@meta.data$Patient[Bmem_fourth@meta.data$hash.ID =="PFCL1-199474-1994-3"] <- "PFCL1-199474-1994"
Bmem_fourth@meta.data$Patient[Bmem_fourth@meta.data$hash.ID =="PFCL1-179308-1988-2"] <- "PFCL1-179308-1988"
Bmem_fourth@meta.data$Patient[Bmem_fourth@meta.data$hash.ID =="PFCL1-179308-1988-3"] <- "PFCL1-179308-1988"
Bmem_fourth@meta.data$Patient[Bmem_fourth@meta.data$hash.ID =="PFCL1-LIM-313000-2"] <- "PFCL1-LIM-313000"
Bmem_fourth@meta.data$Patient[Bmem_fourth@meta.data$hash.ID =="PFCL1-LIM828246-2"] <- "PFCL1-LIM828246"
Bmem_fourth@meta.data$Patient[Bmem_fourth@meta.data$hash.ID =="PFCL1-LIM828246-3"] <- "PFCL1-LIM828246"

# Adding the info which samples are vaccinated
Bmem_fourth@meta.data$Vaccinated <- "no"
Bmem_fourth@meta.data$Vaccinated[Bmem_fourth@meta.data$hash.ID=="PFCL1-199474-1994-3"|Bmem_fourth@meta.data$hash.ID=="PFCL1-179308-1988-3"|
                                   Bmem_fourth@meta.data$hash.ID=="PFCL1-LIM828246-3"] <- "yes"

# Adding the info when 12 months sample was taken with respect to the vaccination
Bmem_fourth@meta.data$time_after_vac <- "not vaccinated"
Bmem_fourth@meta.data$time_after_vac[Bmem_fourth@meta.data$hash.ID=="PFCL1-199474-1994-3"] <- "108d"
Bmem_fourth@meta.data$time_after_vac[Bmem_fourth@meta.data$hash.ID=="PFCL1-179308-1988-3"] <- "85d"
Bmem_fourth@meta.data$time_after_vac[Bmem_fourth@meta.data$hash.ID=="PFCL1-LIM828246-3"] <- "87d"

# Adding a column that combines timepoint and vaccination status
Bmem_fourth@meta.data$sample_group <- "6 months"
Bmem_fourth@meta.data$sample_group[Bmem_fourth@meta.data$Timepoint=="12 months after infection" & Bmem_fourth@meta.data$Vaccinated=="yes"] <- "12 months vaccinated"

######################################################################################################################
# Part 4.2: Fourth Set - Adding VDJ data, filtering of cells with BCR, Addition of Isotype
######################################################################################################################
# Loading the data
contigs <- read.csv("./outs/per_sample_outs/multi_exp035_4/vdj_b/filtered_contig_annotations.csv")
contig.list <- createHTOContigList(contigs, Bmem_fourth, group.by = "hash.ID")
summary(contig.list)

#Fixing the parseBCR function
combined_fourth <- combineBCR(contig.list,
                             samples = c("PFCL1-LIM828246_3", "PFCL1-179308-1988_2", "PFCL1-179308-1988_3", "PFCL1-199474-1994_2",
                                         "PFCL1-LIM828246_2","PFCL1-199474-1994_3","PFCL1-LIM-313000_2"))

remove(contigs,contig.list)

# The row names of BCRcombined and Bmem need to match, otherwise "combineExpression" fails
head(rownames(Bmem_fourth@meta.data))
head(combined_fourth$`PFCL1-LIM828246_3`$barcode)

for (i in seq_along(combined_fourth)) {
  combined_fourth[[i]] <- stripBarcode(combined_fourth[[i]], 
                                      column = 1, connector = "_", num_connects = 3)
}

Bmem_fourth <- combineExpression(combined_fourth, Bmem_fourth, cloneCall="aa", cloneTypes=c(Single=1, Small=5, Medium=20, Large=100, Hyperexpanded=500),proportion = FALSE)
table(Bmem_fourth@meta.data$Frequency)

# Add info whether cell has BCR info
Bmem_fourth@meta.data$BCR.known <- "yes"
Bmem_fourth@meta.data$BCR.known[is.na(Bmem_fourth@meta.data$CTaa)] <- "no"
table(Bmem_fourth@meta.data$BCR.known)

# Add info whether cell has both chains
Bmem_fourth@meta.data$full.BCR.known <- "yes"
Bmem_fourth@meta.data$full.BCR.known[grep("NA",Bmem_fourth@meta.data$CTgene)] <- "no"
Bmem_fourth@meta.data$full.BCR.known[is.na(Bmem_fourth@meta.data$CTaa)]<- "no"
table(Bmem_fourth@meta.data$full.BCR.known)

#Removing cells without BCR --> Could be skipped?
Bmem_fourth <- subset(Bmem_fourth, subset = BCR.known =="yes")

# Adding Isotype Info
# Adding a "delete_me" mock CTaa to all cells with no BCR info
Bmem_fourth@meta.data$CTgene[is.na(Bmem_fourth@meta.data$CTgene)] <- "delete_me"
heavychains <- unlist(strsplit(Bmem_fourth@meta.data$CTgene, "[_]"))[seq(1, length(unlist(strsplit(Bmem_fourth@meta.data$CTgene, "[_]"))), 2)]
heavychains <- as.data.frame(heavychains)
heavychains$Isotype <- "unknown"
heavychains$Isotype[grep("IGHD",heavychains$heavychains)] <- "IGHD"
heavychains$Isotype[grep("IGHA",heavychains$heavychains)] <- "IGHA"
heavychains$Isotype[grep("IGHM",heavychains$heavychains)] <- "IGHM"
heavychains$Isotype[grep("IGHG",heavychains$heavychains)] <- "IGHG"
Bmem_fourth@meta.data$Isotype <- heavychains$Isotype
Bmem_fourth@meta.data$CTgene[Bmem_fourth@meta.data$CTgene=="delete_me"] <- NA
remove(heavychains)

######################################################################################################################
# Part 4.3: Fourth Set - Processing of baiting counts (cutoffs, normalization) and visualization
######################################################################################################################
# Pre-processing baiting counts
Baiting_df <- as.data.frame(Bmem_fourth@assays$Baiting@counts)
Baiting_df <- t(Baiting_df)
Baiting_df <- as.data.frame(Baiting_df)
colnames(Baiting_df) <- c("wt_Spike", "B.1.351","B.1.617.2","RBD","NegControl")

# Density plots
melted_Baiting_df <- reshape2::melt(Baiting_df)
melted_Baiting_df <- melted_Baiting_df[melted_Baiting_df$value>0,]
ggplot(melted_Baiting_df, aes(x = value, y = variable,fill =variable)) +
  geom_density_ridges(scale = 4)+
  scale_x_continuous(trans='log2',breaks = c(1,5,10,20,40,100,200,500,700,1000,1500,3000,5000))+
  labs(title = "Density plots before negative control subtraction") +
  theme_ridges(font_size = 13, grid = TRUE) + 
  theme(axis.title.y = element_blank(),axis.title.x = element_blank())
remove(melted_Baiting_df)

# Correcting by using the negative control
Baiting_df$cor_wt_Spike <- Baiting_df$wt_Spike-Baiting_df$NegControl
Baiting_df$cor_B.1.351 <- Baiting_df$B.1.351-Baiting_df$NegControl
Baiting_df$cor_B.1.617.2 <- Baiting_df$B.1.617.2-Baiting_df$NegControl
Baiting_df$cor_RBD <- Baiting_df$RBD-Baiting_df$NegControl

# Density plots after negative control subtraction
melted_Baiting_df <- Baiting_df[6:9]
melted_Baiting_df <- reshape2::melt(melted_Baiting_df)
melted_Baiting_df <- melted_Baiting_df[melted_Baiting_df$value>0,]
ggplot(melted_Baiting_df, aes(x = value, y = variable,fill =variable)) +
  geom_density_ridges(scale = 4)+
  scale_x_continuous(trans='log2',breaks = c(1,5,10,20,40,100,200,500,700,1000,1500,3000,5000))+
  labs(title = "Density plots after negative control subtraction") +
  theme_ridges(font_size = 13, grid = TRUE) + 
  theme(axis.title.y = element_blank(),axis.title.x = element_blank())

# Preparing a density plot for the publication's supplementary figures.
melted_Baiting_df <- melted_Baiting_df[melted_Baiting_df$variable=="cor_wt_Spike",]
melted_Baiting_df$variable <- "Wuhan-Hu-1 full \nlength spike construct"
ggplot(melted_Baiting_df, aes(x = value,fill=variable)) +
  geom_density(color="darkblue", fill="dodgerblue3",alpha=0.4)+
  scale_x_continuous(trans='log10',breaks = c(1,5,10,20,40,100,200,500,700,1000,1500,3000,5000))+
  ggtitle("Density of count frequencies across all cells ") +
  theme_classic()+center.title()+
  theme(axis.text.x = element_text(angle = 45,hjust = 1),text = element_text(size = 15))+ guides(fill="none")+
  xlab("Barcode counts per cell (log10)")+ylab("Density")+black.axis.text()+geom_vline(xintercept = 40,color = "red", size=1,linetype="dashed")+
  annotate(geom="text", x=65, y=0.8, label="Cutoff", color="red",size=5)

remove(melted_Baiting_df)

# Classify binders
Baiting_df$cor_wt_Spike[Baiting_df$cor_wt_Spike<1] <- 0
Baiting_df$cor_B.1.351[Baiting_df$cor_B.1.351<1] <- 0
Baiting_df$cor_B.1.617.2[Baiting_df$cor_B.1.617.2<1] <- 0
Baiting_df$cor_RBD[Baiting_df$cor_RBD<1] <- 0

Baiting_df$wt_Spike.classification <- "Negative"
Baiting_df$B.1.351.classification <- "Negative"
Baiting_df$B.1.617.2.classification <- "Negative"
Baiting_df$RBD.classification <- "Negative"

Baiting_df$wt_Spike.classification[Baiting_df$cor_wt_Spike>40] <- "Positive"
Baiting_df$B.1.351.classification[Baiting_df$cor_B.1.351>10] <- "Positive"
Baiting_df$B.1.617.2.classification[Baiting_df$cor_B.1.617.2>10] <- "Positive"
Baiting_df$RBD.classification[Baiting_df$cor_RBD>200] <- "Positive"

# Adding a column that tells, whether a cell is positive for at least one bait and for RBD + full length Spike
Baiting_df$bait.positive <- "no"
Baiting_df$bait.positive[Baiting_df$wt_Spike.classification=="Positive"|
                           Baiting_df$B.1.351.classification=="Positive"|
                           Baiting_df$B.1.617.2.classification=="Positive"|
                           Baiting_df$RBD.classification=="Positive"] <- "yes"

Baiting_df$wt_Spike_RBD.positive <- "no"
Baiting_df$wt_Spike_RBD.positive[Baiting_df$wt_Spike.classification=="Positive"&
                                   Baiting_df$RBD.classification=="Positive"] <- "yes"

# Setting corrected baiting counts to NA if classified as binding negative
Baiting_df$cor_wt_Spike[Baiting_df$wt_Spike.classification=="Negative"] <- NA
Baiting_df$cor_B.1.351[Baiting_df$B.1.351.classification =="Negative"] <- NA
Baiting_df$cor_B.1.617.2[Baiting_df$B.1.617.2.classification =="Negative"] <- NA
Baiting_df$cor_RBD[Baiting_df$RBD.classification=="Negative"] <- NA

# Density plots after cutoff correction
melted_Baiting_df <- Baiting_df[6:9]
melted_Baiting_df <- reshape2::melt(melted_Baiting_df)
melted_Baiting_df <- na.omit(melted_Baiting_df)
ggplot(melted_Baiting_df, aes(x = value, y = variable,fill =variable)) +
  geom_density_ridges(scale = 4)+
  scale_x_continuous(trans='log2',breaks = c(1,5,10,20,40,100,200,500,700,1000,1500,3000,5000))+
  labs(title = "Density plots after negative control subtraction") +
  theme_ridges(font_size = 13, grid = TRUE) + 
  theme(axis.title.y = element_blank(),axis.title.x = element_blank())
remove(melted_Baiting_df)

# In order to use the Seurat Normalization function, we need to have the corrected counts as assay object
# For this, the Baiting_df needs to be transposed
Baiting_df <- as.data.frame(t(Baiting_df[,6:9]))
Cor_Baiting_assay <- CreateAssayObject(counts = Baiting_df)
Bmem_fourth[["Cor_Baiting"]] <- Cor_Baiting_assay
remove(Cor_Baiting_assay)
Bmem_fourth <- NormalizeData(Bmem_fourth, assay = "Cor_Baiting", normalization.method = 'CLR', margin = 1) # Decide here: Normalize across features or cells.
normalized_across_cells <- as.data.frame(t(Bmem_fourth@assays$Cor_Baiting@data))
colnames(normalized_across_cells) <- c("cor_wt_Spike","cor_B.1.351","cor_B.1.617.2","cor_RBD")
remove(Baiting_df)
Bmem_fourth@meta.data$nCount_Cor_Baiting <- NULL
Bmem_fourth@meta.data$nFeature_Cor_Baiting <- NULL

# For scaling, I cannot use Seurat Scaling function because it gives 0s as result when NAs are present.
# I also need to transform again, then each column is scaled to 0:1 individually
normalized_across_cells <- myscaling(normalized_across_cells)  

# These are the final LIBRA scores
# Adding this to the Seurat object
Bmem_fourth@meta.data <- merge(Bmem_fourth@meta.data,normalized_across_cells, by=0)
rownames(Bmem_fourth@meta.data) <- Bmem_fourth@meta.data$Row.names
remove(normalized_across_cells)

# Feature scatters
# NAs are now set to 0
Bmem_fourth@meta.data$cor_wt_Spike[is.na(Bmem_fourth@meta.data$cor_wt_Spike)] <- 0
Bmem_fourth@meta.data$cor_B.1.351[is.na(Bmem_fourth@meta.data$cor_B.1.351)] <- 0
Bmem_fourth@meta.data$cor_B.1.617.2[is.na(Bmem_fourth@meta.data$cor_B.1.617.2)] <- 0
Bmem_fourth@meta.data$cor_RBD[is.na(Bmem_fourth@meta.data$cor_RBD)] <- 0

FeatureScatter(Bmem_fourth,  feature1 = "cor_wt_Spike", feature2 = "cor_RBD", pt.size = 1, group.by = "Isotype")+
  ggtitle("Wt Spike and RBD scores")+xlab("wt-Spike scores")+ylab("RBD scores")+
  
  FeatureScatter(Bmem_fourth,  feature1 = "cor_wt_Spike", feature2 = "cor_B.1.351", pt.size = 1, group.by = "Isotype")+
  ggtitle("Wt Spike and Beta variant scores")+xlab("wt-Spike scores")+ylab("Beta variant scores")+
  
  FeatureScatter(Bmem_fourth,  feature1 = "cor_wt_Spike", feature2 = "cor_B.1.617.2", pt.size = 1, group.by = "Isotype")+
  ggtitle("Wt Spike and Delta variant scores")+xlab("wt-Spike scores")+ylab("Delta variant scores")+
  
  FeatureScatter(Bmem_fourth,  feature1 = "cor_B.1.351", feature2 = "cor_B.1.617.2", pt.size = 1, group.by = "Isotype")+
  ggtitle("Beta and Delta variant scores")+xlab("Beta scores")+ylab("Delta variant scores")


# Addition of bait positive classification columns in metadata
Bmem_fourth@meta.data$bait.positive <- "no"
Bmem_fourth@meta.data$bait.positive[Bmem_fourth@meta.data$cor_wt_Spike>0 | 
                                      Bmem_fourth@meta.data$cor_B.1.351>0 | 
                                      Bmem_fourth@meta.data$cor_B.1.617.2 >0 | 
                                      Bmem_fourth@meta.data$cor_RBD > 0] <- "yes"
Bmem_fourth@meta.data$spike_RBD_positive <- "no"
Bmem_fourth@meta.data$spike_RBD_positive[Bmem_fourth@meta.data$cor_wt_Spike>0 &
                                           Bmem_fourth@meta.data$cor_RBD >0] <- "yes"

Bmem_fourth@meta.data$B.1.351_positive <- "no"
Bmem_fourth@meta.data$B.1.351_positive[Bmem_fourth@meta.data$cor_B.1.351>0] <- "yes"

Bmem_fourth@meta.data$B.1.617.2_positive <- "no"
Bmem_fourth@meta.data$B.1.617.2_positive[Bmem_fourth@meta.data$cor_B.1.617.2>0] <- "yes"

Bmem_fourth@meta.data$all_positive <- "no"
Bmem_fourth@meta.data$all_positive[Bmem_fourth@meta.data$cor_wt_Spike>0 & 
                                     Bmem_fourth@meta.data$cor_B.1.351>0 & 
                                     Bmem_fourth@meta.data$cor_B.1.617.2 >0 & 
                                     Bmem_fourth@meta.data$cor_RBD > 0] <- "yes"

# Plots to show Isotypes of bait positives
barplot(table(Bmem_fourth@meta.data$Isotype[Bmem_fourth@meta.data$bait.positive=="yes"]), main = "Isotypes of bait positive BCRs")
barplot(table(Bmem_fourth@meta.data$Isotype[Bmem_fourth@meta.data$spike_RBD_positive=="yes"]), main = "Isotypes of RBD + full length Spike positive BCRs")

# Renaming of cells
Bmem_fourth <- RenameCells(object = Bmem_fourth, add.cell.id = "Dataset_4")
Bmem_fourth$Dataset <- "Fourth"


######################################################################################################################
# Part 4.4 Identifying and excluding naive B cells in the Dataset 4
######################################################################################################################
# The following section is used to define naive B cells present in Dataset 4:
# For this, I am following this strategy:
# Make a new Seurat object that can me manipulated
Bmem_fourth.IgD <- Bmem_fourth

# Cluster cells by WNN, using all TotalSeqs
Bmem_fourth.IgD <- FindVariableFeatures(Bmem_fourth.IgD)
Bmem_fourth.IgD <- ScaleData(Bmem_fourth.IgD)
Bmem_fourth.IgD <- RunPCA(Bmem_fourth.IgD)

DefaultAssay(Bmem_fourth.IgD) <- "Protein"
VariableFeatures(Bmem_fourth.IgD) <- rownames(Bmem_fourth.IgD[["Protein"]])
Bmem_fourth.IgD <-ScaleData(Bmem_fourth.IgD) 
Bmem_fourth.IgD <- RunPCA(Bmem_fourth.IgD, reduction.name = 'apca')

Bmem_fourth.IgD <- FindMultiModalNeighbors(
  Bmem_fourth.IgD, reduction.list = list("pca", "apca"), 
  dims.list = list(1:30, 1:5), modality.weight.name = "RNA.weight")

Bmem_fourth.IgD <- RunUMAP(Bmem_fourth.IgD, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
Bmem_fourth.IgD <- FindClusters(Bmem_fourth.IgD, graph.name = "wsnn", algorithm = 3, resolution = 0.5, verbose = FALSE)

# Check UMAP
DimPlot(Bmem_fourth.IgD, reduction = 'wnn.umap', label = TRUE, repel = TRUE, label.size = 2.5) # Export as 4x7
DefaultAssay(Bmem_fourth.IgD) <- "RNA"
#Markers <- FindAllMarkers(Bmem_fourth.IgD)
Markers <- Markers %>% group_by(cluster) %>% slice_max(n = 10, order_by = avg_log2FC)
Markers <- Markers[Markers$p_val_adj<0.05,]
table(Bmem_fourth.IgD@meta.data$seurat_clusters, Bmem_fourth.IgD@meta.data$bait.positive)

# Check possibilty of patient driven clustering
DimPlot(Bmem_fourth.IgD, reduction = 'wnn.umap', label = TRUE, repel = TRUE, label.size = 2.5, group.by = "Patient")

# Check which clusters can be exluded based on TotalSeqs
VlnPlot(Bmem_fourth.IgD, features = c("FcRL5","CD27.1","IgD","CD71","CXCR5.1","CD21"),assay = "Protein",pt.size = 0) # Export 4x7

# I decide to mark cells of clusters 0 and 1 as naive. This decision is made upon inspection of protein expression levels and RNA Marker genes.
Naives <- rownames(Bmem_fourth.IgD@meta.data[Bmem_fourth.IgD@meta.data$seurat_clusters==0 |
                                               Bmem_fourth.IgD@meta.data$seurat_clusters==1,])
writeLines(Naives, "Naives.txt")
#saveRDS(Bmem_fourth.IgD, "Bmem_fourth.IgD.rds")
remove(Bmem_fourth.IgD, Markers)

# Keeping only non-naive cells from Dataset 4
Bmem_fourth@meta.data$naive <- "no"
Bmem_fourth@meta.data$naive[rownames(Bmem_fourth@meta.data) %in% Naives] <- "yes"
Bmem_fourth <- subset(Bmem_fourth, subset= naive =="no")
Bmem_fourth@meta.data$naive <- NULL
Bmem_fourth@meta.data$nCount_Cor_Baiting <- NULL
Bmem_fourth@meta.data$nFeature_Cor_Baiting <- NULL

######################################################################################################################
# Part 5: Data set integration and WNN + Addition of mutational load information for cells
######################################################################################################################
# From here on, I am following the instructions for integration from:
# https://satijalab.org/seurat/articles/integration_introduction.html
# During the IntegrateData() function, dimensional reductions of the data are needed. I have tested two different strategies here:
# In option A, PCA is automatically being calculated during the IntegrateData() function, as no dimensional reduction was calculated before.
# But this gave a warning message.
# In option B, PCA is calculated before Integration is done (only for the protein data). --> No warning and very similar results with better P values.
# In the following, I therefore only kept option B.
# Create a list containing the two datasets
Bmem <- ScaleData(Bmem, assay = "Protein")
Bmem <- RunPCA(Bmem, assay = "Protein" ,approx = F)
Bmem_pilot <- ScaleData(Bmem_pilot, assay = "Protein")
Bmem_pilot <- RunPCA(Bmem_pilot, assay = "Protein", approx=F)
Bmem_third <- ScaleData(Bmem_third, assay = "Protein")
Bmem_third <- RunPCA(Bmem_third, assay = "Protein", approx=F)
Bmem_fourth <- ScaleData(Bmem_fourth, assay = "Protein")
Bmem_fourth <- RunPCA(Bmem_fourth, assay = "Protein", approx=F)

Bmem_dataset_list <- list(Bmem_pilot,Bmem,Bmem_third,Bmem_fourth)

# Select features that are repeatedly variable across datasets for integration
# The integration is done independently for the RNA and the Protein assay.
features_RNA <- SelectIntegrationFeatures(object.list = Bmem_dataset_list,assay = c("RNA","RNA","RNA","RNA"))
features_Protein <- SelectIntegrationFeatures(object.list = Bmem_dataset_list,assay = c("Protein","Protein","Protein","Protein"))

# Perform integration. This command requires long time / a lot of computational resources to run.
# Because of only 5 protein features, dims is set to 4 for the anchor finding of protein data -> https://github.com/satijalab/seurat/issues/5089
anchors_RNA <- FindIntegrationAnchors(object.list = Bmem_dataset_list, anchor.features = features_RNA,assay = c("RNA","RNA","RNA","RNA"))
anchors_Protein <- FindIntegrationAnchors(object.list = Bmem_dataset_list, anchor.features = features_Protein,assay = c("Protein","Protein","Protein","Protein"),dims=1:4)

# This command creates an 'integrated' data assay in a new Seurat object.
# Because of only 5 protein features, dims is set to 4 for the integration of protein data -> https://github.com/satijalab/seurat/issues/5089
# The way I am doing this is inspired by https://github.com/timoast/signac/discussions/438.
# See the respective power point to see the results of some comparisons of different ways to integrate that I did.
Bmem.combined.RNA <- IntegrateData(anchorset = anchors_RNA,new.assay.name = "integratedRNA")
Bmem.combined.RNA <- FindVariableFeatures(Bmem.combined.RNA, assay = "Protein")
Bmem.combined.RNA <- ScaleData(Bmem.combined.RNA, assay = "Protein")
Bmem.combined.RNA <- RunPCA(Bmem.combined.RNA, assay = "Protein", reduction.name="pca_protein_all_cells")
Bmem.combined.Protein <- IntegrateEmbeddings(anchorset = anchors_Protein, reductions = Bmem.combined.RNA@reductions$pca_protein_all_cells)
remove(anchors_Protein,anchors_RNA,Bmem_dataset_list,Bmem_pilot,Bmem,Bmem_third,Bmem_fourth)

# Now WNN
# Following:
# https://github.com/satijalab/seurat/issues/2706
# https://github.com/satijalab/seurat/issues/3645
# https://github.com/satijalab/seurat/issues/3843
# https://github.com/satijalab/seurat/issues/3890
# https://github.com/satijalab/seurat/issues/4922
# Now I want to perform WNN. For this, I need to transfer the integrated protein dim reduction from the Bmem.combined.Protein object to the
# Bmem.combined.RNA object.
# But first, I am preforming data scaling and PCA on the combined RNA project.
DefaultAssay(Bmem.combined.RNA)
DefaultAssay(Bmem.combined.Protein)

Bmem.combined.RNA <- ScaleData(Bmem.combined.RNA, verbose = FALSE)
Bmem.combined.RNA <- RunPCA(Bmem.combined.RNA, npcs = 30, verbose = FALSE)

# Now we transfer the integrated protein dim reduction
Bmem.combined.RNA@reductions$IntegrateEmbeddings.pca <- Bmem.combined.Protein@reductions$integrated_dr

# I double check that the correct assays are associated with the correct PCAs.
Bmem.combined.RNA[['pca']]@assay.used
Bmem.combined.RNA[['IntegrateEmbeddings.pca']]@assay.used

# Finally, the FindMultiModalNeighbors function is called, using the PCAs from above.
Bmem.combined.RNA <- FindMultiModalNeighbors(
  Bmem.combined.RNA, reduction.list = list("pca", "IntegrateEmbeddings.pca"), 
  dims.list = list(1:30, 1:4), modality.weight.name = c("RNA.weight","Protein.weight"))

Bmem <- Bmem.combined.RNA
Bmem@meta.data$Full.Row.names <- rownames(Bmem@meta.data)
DefaultAssay(Bmem) <- "RNA"
Bmem <- ScaleData(Bmem)
remove(features_Protein,features_RNA,Bmem.combined.Protein,Bmem.combined.RNA)

# Read in the mutational load file
# I am loading a file which we have prepared with help of the Immcanation pipeline, which contains the mutational counts of all cells in the dataset.
# This file is included in the Zenodo upload of the dataset.
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

Mutations <- Mutations[,c(63,64,50)]

# Do the merge
Bmem@meta.data <- merge(Bmem@meta.data, Mutations, by = 0, all.x = T)
Bmem@meta.data <- Bmem@meta.data[,2:ncol(Bmem@meta.data)]
rownames(Bmem@meta.data) <- Bmem@meta.data$Full.Row.names

remove(Mutations,i)


######################################################################################################################
# Part 6: Clustering, UMAP visualization and Phenotyping of cells
######################################################################################################################
# First clustering, UMAP and Dimplot on RNA data only
Bmem <- FindNeighbors(Bmem, dims = 1:30,graph.name = "integratedRNA_nn",assay="integratedRNA",reduction = "pca")
Bmem <- FindClusters(Bmem, graph.name = "integratedRNA_nn", algorithm = 1, resolution = 0.5, verbose = FALSE)
Bmem <- RunUMAP(Bmem, reduction = 'pca', dims = 1:30, assay = 'integratedRNA', reduction.name = 'rna.umap', reduction.key = 'rnaUMAP_')

DimPlot(Bmem, reduction = 'rna.umap', label = F, repel = TRUE, label.size = 2.5, pt.size = 0.6,group.by = "seurat_clusters", 
        cols = MetBrewer::met.brewer(name="Egypt", n=length(unique(Bmem@meta.data$seurat_clusters))))+ggtitle("RNA UMAP")+center.title() # Export 4x6

# Then clustering, UMAP and Dimplot on Protein data only
Bmem <- FindNeighbors(Bmem, dims = 1:4,graph.name = "integratedADT_nn",assay="Protein",reduction = "IntegrateEmbeddings.pca")
Bmem <- FindClusters(Bmem, graph.name = "integratedADT_nn", algorithm = 1, resolution = 0.5, verbose = FALSE)
Bmem <- RunUMAP(Bmem, reduction = 'IntegrateEmbeddings.pca', dims = 1:3, assay = 'Protein', reduction.name = 'adt.umap', reduction.key = 'adtUMAP_')

DimPlot(Bmem, reduction = 'adt.umap', label = F, repel = TRUE, label.size = 2.5, pt.size = 0.6, 
        cols = MetBrewer::met.brewer(name="Egypt", n=length(unique(Bmem@meta.data$seurat_clusters))))+ggtitle("Protein UMAP")+center.title()

# Then on the WNN data - find multimodal neighbours done already in Part 5.
Bmem <- FindClusters(Bmem, graph.name = "wsnn", algorithm = 1, resolution = 0.4, verbose = FALSE)
a <- 1.4
b <- 0.75
Bmem <- RunUMAP(Bmem, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_",a=a ,b=b)

# Supp.Fig.5.B.1
DimPlot(Bmem, reduction = 'wnn.umap', label = T, repel = TRUE, label.size = 2.5, pt.size = 0.5, 
        cols = MetBrewer::met.brewer(name="Egypt", n=length(unique(Bmem@meta.data$seurat_clusters))))+
  theme(text = element_text(size = 20)) +ggtitle(paste("WNN UMAP"))+center.title() # Export 4x7
DimPlot(Bmem, reduction = 'wnn.umap', label = F, repel = TRUE, label.size = 2.5, pt.size = 0.6,
        cols = MetBrewer::met.brewer(name="Egypt", n=length(unique(Bmem@meta.data$seurat_clusters))))+NoAxes() # Export 4x7

DefaultAssay(Bmem) <- "Protein"
VlnPlot(Bmem, features = c("FcRL5", "CD21", "CD27.1", "CXCR5.1", "CD71"), pt.size = 0.03, 
        cols = MetBrewer::met.brewer(name="Egypt", n=length(unique(Bmem@meta.data$seurat_clusters)))) # Export 8x14


# Based on the expression of surface markers, Isotypes and CD19 I give the following names to the WNN clusters:
Bmem@meta.data$named.clusters <- ""
Bmem@meta.data$named.clusters[Bmem@meta.data$seurat_clusters==0]<-"Unswitched"
Bmem@meta.data$named.clusters[Bmem@meta.data$seurat_clusters==1]<-"CD27low RM"
Bmem@meta.data$named.clusters[Bmem@meta.data$seurat_clusters==2]<-"CD27high RM"
Bmem@meta.data$named.clusters[Bmem@meta.data$seurat_clusters==3]<-"CD27high RM"
Bmem@meta.data$named.clusters[Bmem@meta.data$seurat_clusters==4]<-"Activated"
Bmem@meta.data$named.clusters[Bmem@meta.data$seurat_clusters==5]<-"CD27high RM"
Bmem@meta.data$named.clusters[Bmem@meta.data$seurat_clusters==6]<-"Unswitched"
Bmem@meta.data$named.clusters[Bmem@meta.data$seurat_clusters==7]<-"Atypical"
Bmem@meta.data$named.clusters[Bmem@meta.data$seurat_clusters==8]<-"CD27high RM"
Bmem@meta.data$named.clusters[Bmem@meta.data$seurat_clusters==9]<-"Activated"

# In addition, I set all class switched cells from clusters 0 and 6 as "CD27low RM"
Bmem@meta.data$named.clusters[Bmem@meta.data$seurat_clusters==0 & Bmem@meta.data$Isotype=="IGHA"] <- "CD27low RM"
Bmem@meta.data$named.clusters[Bmem@meta.data$seurat_clusters==0 & Bmem@meta.data$Isotype=="IGHG"] <- "CD27low RM"
Bmem@meta.data$named.clusters[Bmem@meta.data$seurat_clusters==6 & Bmem@meta.data$Isotype=="IGHA"] <- "CD27low RM"
Bmem@meta.data$named.clusters[Bmem@meta.data$seurat_clusters==6 & Bmem@meta.data$Isotype=="IGHG"] <- "CD27low RM"

# Fig.5.A.1
DimPlot(Bmem,reduction = "wnn.umap",group.by = "named.clusters", pt.size = 0.5, 
        cols = c("#c2c1c0","#228833","#4477AA","#AA3377","#DBC35E"))+
  theme(text = element_text(size = 20)) + NoAxes()+ggtitle("")# Export 4x8

Bmem@meta.data$named.clusters <- factor(Bmem@meta.data$named.clusters, levels = c("Unswitched","CD27low RM","CD27high RM","Atypical","Activated"))
Idents(Bmem) <- "named.clusters"
VlnPlot(Bmem, assay = "Protein",features = c("CD27.1","CD21","FcRL5"),pt.size = 0, 
        cols = rev(MetBrewer::met.brewer(name="Austria", n=length(unique(Bmem@meta.data$seurat_clusters)))[1:5]))# Export 5x6

# The object that is saved at this point, is the "final", preprocessed Seurat object, which is used for all downstream analysis.
saveRDS(Bmem, "Bmem.rds")

######################################################################################################################
# Part 7: Heatmap for segment usages - Fig_1_G and Supp.Fig.2.D
######################################################################################################################
# From here on we start to make the heatmap that compares segment usage between groups of cells.
# This part is still part of the preprocessing script because I need objects which are created during the pre processing (the "combined" objects).
# Rownames of Bmem don't match with "combined" any more after integration. I am adding the required prefixes to the "combined" lists.
i <- 1
for (i in 1:length(combined)) {
  combined[[i]]$barcode <- paste0("Dataset_2_",combined[[i]]$barcode)
  i <- i+1
}
i <- 1
for (i in 1:length(combined_pilot)) {
  combined_pilot[[i]]$barcode <- paste0("Dataset_1_",combined_pilot[[i]]$barcode)
  i <- i+1
}
i <- 1
for (i in 1:length(combined_third)) {
  combined_third[[i]]$barcode <- paste0("Dataset_3_",combined_third[[i]]$barcode)
  i <- i+1
}
i <- 1
for (i in 1:length(combined_fourth)) {
  combined_fourth[[i]]$barcode <- paste0("Dataset_4_",combined_fourth[[i]]$barcode)
  i <- i+1
}

# Gene segment usage with integrated dataset
combined_list <- c(combined_pilot,combined,combined_third,combined_fourth)

# Use of the vizGenes function to compare the gene usages of RBD positive and negative cells.
i <- 1
df <- combined_list[[1]][FALSE,]
for (i in 1:length(combined_list)) {
  df_loop <- combined_list[[i]]
  df <- rbind(df,df_loop)
  i <- i+1
}

df_positives <- df[df$barcode %in% rownames(Bmem@meta.data[Bmem@meta.data$spike_RBD_positive=="yes",]),]
df_negatives <- df[df$barcode %in% rownames(Bmem@meta.data[Bmem@meta.data$bait.positive=="no",]),]

df_positives$sample <- "positive"
df_negatives$sample <- "negative"
df_all <- rbind(df_positives,df_negatives)
df_positives_list <- list(df_positives)
df_negatives_list <- list(df_negatives)

trace(vizGenes, edit = T) # Change output to df instead of plot - needs to be change the first time of every session.

# Decide here which gene segments should be looked at. Specify chain "IGH" or "IGL" for heavy and light chain.
positives <- vizGenes(df_positives_list, gene = "V", chain = "IGH", plot = "bar", order = "variance", scale = T)
negatives <- vizGenes(df_negatives_list, gene = "V", chain = "IGH", plot = "bar", order = "variance", scale = T)

totalforstats_positives <- unique(positives$sum)
totalforstats_negatives <- unique(negatives$sum)

positives <- unique(positives[, c("Var1", "Var2", "sd", "mean")])
negatives <- unique(negatives[, c("Var1", "Var2", "sd", "mean")])
positives$mean <- positives$mean*100
negatives$mean <- negatives$mean*100
positives$Set <- "RBD binder"
negatives$Set <- "Non binder"

adder <- as.data.frame(negatives[negatives$Var1 %notin% positives$Var1,])
adder$mean <- 0
adder$Set <- "RBD binder"
total <- rbind(positives,adder,negatives)

# Plotting the result as bar chart
ggplot(total, aes(x = Var1, y = mean, fill=Set)) + geom_bar(stat = "identity",position=position_dodge()) + 
  geom_errorbar(aes(ymin = mean, ymax = mean + sd), 
                width = 0.2, position = position_dodge(0.9)) + 
  theme_classic() + theme(axis.title.x = element_blank(), 
                          axis.title.y = element_blank(), axis.ticks.x = element_blank(), 
                          axis.text.x = element_text(angle = 90, vjust = 0.5, 
                                                     hjust = 1, size = rel(0.5)))+
  theme(text = element_text(size=30),plot.title = element_text(hjust = 0.5),legend.title=element_blank())+
  scale_fill_manual(values=met.brewer("Austria", 2))+ggtitle("V Light Frequency (%)")+
  scale_x_discrete(expand = c(0, 0))+scale_y_continuous(expand = c(0,0) )

# Show results as heatmap
total <- total[,c(1,4,5)]
total <- as.data.frame(total)
total_Binders <- total[total$Set=="RBD binder",]
total_non_Binders <- total[total$Set=="Non binder",]
rownames(total_Binders) <- total_Binders$Var1
rownames(total_non_Binders) <- total_non_Binders$Var1
total.new <- merge(total_Binders, total_non_Binders, all.x = T, by=0)
rownames(total.new) <- total.new$Row.names
total.new <- total.new[,c(3,6)]
colnames(total.new) <- c("RBD Binder", "Non Binder")
total.new$`RBD non Binder`[is.na(total.new$`RBD non Binder`)] <- 0
total.new.scaled <-myscaling(total.new)

# Subset to only top 30 RBD binding segments
total.new <- total.new[order(total.new$`RBD Binder`, decreasing = T),]
total.new <- total.new[1:30,]

make_bold_names <- function(mat, rc_fun, rc_names) {
  bold_names <- rc_fun(mat)
  ids <- rc_names %>% match(rc_fun(mat))
  ids %>%
    walk(
      function(i)
        bold_names[i] <<-
        bquote(bold(.(rc_fun(mat)[i]))) %>%
        as.expression()
    )
  bold_names
}


# Calculating statistics for all rows of the heatmap and adding them to the plot
total.new.stats <- total.new
total.new.stats$RBD_Binder_absolute <- totalforstats_positives*(total.new.stats$`RBD Binder`/100)
total.new.stats$RBD_non_Binder_absolute <- totalforstats_negatives*(total.new.stats$`Non Binder`/100)
total.new.stats$pvalue <- ""
total.new.stats$pvalue.number <- ""								   

i <- 1
for (i in 1:nrow(total.new.stats)) {
  test <- prop.test(x = c(total.new.stats$RBD_Binder_absolute[i], 
                        total.new.stats$RBD_non_Binder_absolute[i]), 
                        n = c(totalforstats_positives,totalforstats_negatives))
total.new.stats$pvalue[i] <- asterics.function(test)
total.new.stats$pvalue.number[i] <- test$p.value												
i <- i+1
}

# Multiple testing comparison
total.new.stats$pvalue.number.corrected <- p.adjust(total.new.stats$pvalue.number, "bonferroni")

total.new.stats$pvalue.corrected <- ""
i <- 1
for (i in 1:nrow(total.new.stats)) {
  total.new.stats$pvalue.corrected[i] <- asterics.function2(as.numeric(total.new.stats$pvalue.number.corrected[i]))
}

# In order to show the significance asteriks in the heatmap, I am creating the display_numbers matrix
display_numbers <- matrix(data="", ncol = 2, nrow = nrow(total.new.stats))
colnames(display_numbers) <- colnames(total.new.stats[,1:2])

# In the following loop, the display_numbers matrix is filled with the correct entries.
i <- 1
for (i in 1:nrow(display_numbers)) {
  if(total.new.stats$pvalue.corrected[i]=="ns") next
  else{
    if(total.new.stats$`RBD Binder`[i]>total.new.stats$`Non Binder`[i]){
      display_numbers[i,1] <- total.new.stats[i,"pvalue.corrected"]
    }
    else {
      display_numbers[i,2] <- total.new.stats[i,"pvalue.corrected"]
    }
  }
}

# Fig. 1.G and Supp.Fig.2.D
# Plotting the heatmap. The plot is also exported as pdf.
cairo_pdf("Gene_usage_heavy_chain.pdf",width = 4, height = 14)
pheatmap(total.new,color =  inferno(80),cutree_rows = 4,cluster_cols=F,fontsize_col=12,
         labels_row = make_bold_names(total.new, rownames, c("IGHV3-30","IGHV3-53","IGHV3-13","IGHV3-66","IGHV1-69D",
                                                             "IGKV1D-39","IGKV1-33","IGLV3-21","IGLV6-57","IGKV1-5","IGKV1-9")),
         display_numbers = display_numbers, fontsize_number=12, number_color = "black",cellwidth = 50)+theme_classic()
grid.text("% VH Segment Usage",0.81,0.92,rot = 90) # Export 14x4
dev.off()


remove(combined_list, adder, df, df_all,df_loop,df_negatives,df_negatives_list,df_positives,df_positives_list,negatives,positives,total,
       total_Binders, total_non_Binders, total.new, total.new.scaled)
remove(combined,combined_fourth,combined_pilot,combined_third,i)

# This is the end of the preprocessing file.