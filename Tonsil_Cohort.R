# This script is used to analyze the scRNAseq data from the MBC project investigating the convalescent patient cohort.
# This is part of our experiments for the revision of the manuscript.
# Two 10x lanes were run on the 16th November 2022. Both lanes contain the same sample. The samples are matched samples of 
# PBMC and Tonsils from the same patients. Four patients are included.

# Loading necessary packagess
library(ggforce)
library(Seurat)
library(tidyverse)
library(EnhancedVolcano)
library(Biostrings)
library(alakazam)
library(shazam)
library(patchwork)
library(ggallin)
library(yarrr)
library(RColorBrewer)
library(umap)
library(ComplexHeatmap)
library(circlize)
library(VennDiagram)
library(cowplot)
library(UpSetR)
library(dplyr)
library(data.table)
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
library(airr)
library(bruceR)
library(clustree)
library(divo)

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
  
} # Highlight chosen cells in the FeatureScatter plot
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
Make.stacked.percentage.plot <- function(Seurat.Object, groups.on.Y.axis, groups.in.legend,X.Axis.title,Y.Axis.title,Legend.title,Plot.title){
  colnum.y.axis <- which( colnames(Seurat.Object@meta.data)==groups.on.Y.axis )
  colnum.groups <- which( colnames(Seurat.Object@meta.data)==groups.in.legend )
  df <- as.data.frame(table(Seurat.Object@meta.data[,colnum.y.axis],Seurat.Object@meta.data[,colnum.groups]))
  df <- df %>% group_by(Var1) %>% add_tally(Freq, name = "Cells.per.Y.axis.category")
  df$groups.in.legend.per.Y.axis.category <- 100*round(df$Freq/df$`Cells.per.Y.axis.category`,3)
  df$Var2 <- factor(df$Var2)
  
  p <- ggplot(df, aes(fill=Var2, y=Var1, x=groups.in.legend.per.Y.axis.category)) + 
    geom_bar(position="fill", stat="identity",colour="black")+theme_classic()+ theme(text = element_text(size = 15))+
    labs(x = X.Axis.title, y=Y.Axis.title)+scale_x_continuous(labels=scales::percent)+labs(fill = Legend.title)+
    scale_fill_manual(values=met.brewer("Austria", length(unique(df$Var2))))+ theme(text = element_text(size = 17))+ggtitle(Plot.title)+
    center.title()+black.axis.text()
  plot(p)
}
Show.clones.on.UMAP <- function(Seurat.Object, clone_IDs, reduction,DimplotGroup){
  # Getting the cells which will be connected
  Cells.to.connect <- Seurat.Object@meta.data$Row.names[Seurat.Object@meta.data$clone_id %in% clone_IDs& !is.na(Seurat.Object@meta.data$clone_id)]
  
  #Getting the UMAP coordinates
  UMAP.coords <- Seurat.Object[["wnn.umap"]]@cell.embeddings
  UMAP.coords <- UMAP.coords[rownames(UMAP.coords) %in% Cells.to.connect,]
  UMAP.coords <- as.data.frame(UMAP.coords)
  
  # Preparing the information of which cells need to be connected with each other
  df.start <- Seurat.Object@meta.data[Cells.to.connect,c("Tissue","Row.names","clone_id")]
  
  # How many arrows will be needed? I want an arrow going from each early timepoint to each later timepoint
  # Therefore the total number of arrows will be (per clone) :
  # n cells in timeopint 1 * n cells in timepoint 2 + n cells in timepoint * n cells in timepoint 3 + n cells in timepoint 3 * n cells in timepoint 4
  # If several clones are given as input, we have to do this for each clone individually
  # Therefore, I will go through the clones one by one.
  df_final <- data.frame()
  col_vector = c("brown","gray21", "purple",
                 "gray46","#FDBF6F","khaki2","maroon","orchid1","deeppink1","blue1","steelblue4",
                 "darkturquoise", "green1", "yellow4", "yellow3","darkorange4", "brown")
  
  i <- 1
  for (i in 1:length(unique(df.start$clone_id))) {
    this.rounds.clone <- unique(df.start$clone_id)[i]
    df.subset <- df.start[df.start$clone_id==this.rounds.clone,]
    
    # Now create all arrows one by one
    # Separate the df.subset into all timepoints
    df.subset_1 <- df.subset[df.subset$Tissue=="PBMC",]
    df.subset_2 <- df.subset[df.subset$Tissue=="Tonsil",]
    
    matrix = matrix(nrow = nrow(df.subset_1)*nrow(df.subset_2), ncol = 2)
    df.one = data.frame(matrix) 
    
    if(nrow(df.one)>0){
      c <- 1
      for (c in 1:nrow(df.subset_1)) {
        df.one[(c*nrow(df.subset_2)-(nrow(df.subset_2)-1)):((c*nrow(df.subset_2)-(nrow(df.subset_2)-1))+nrow(df.subset_2)-1),1] <- df.subset_1$Row.names[c]
        df.one[(c*nrow(df.subset_2)-(nrow(df.subset_2)-1)):((c*nrow(df.subset_2)-(nrow(df.subset_2)-1))+nrow(df.subset_2)-1),2] <- df.subset_2$Row.names
      }
    }
    
    df <- df.one
    df <- drop_na(df)
    
    # Preparing a dataframe that will be given as inputs to geom_line.
    # The dataframe will have three columns. Column x for all the x coordinates, column y for all the y coordinates and column z for the grouping.
    # Each arrow gets its own group. 
    # Adding the groups at this point:
    df$group <- paste0(i,"_",1:nrow(df))
    df$group <- as.character(df$group)
    
    # Then preparing the final dataframe as input for geomline step by step
    colnames(df) <- c("First_cell","Second_cell","group")
    df_first_cell <- df[,c(1,3)]
    df_second_cell <- df[,c(2,3)]
    colnames(df_first_cell) <- c("Cell","group")
    colnames(df_second_cell) <- c("Cell","group")
    UMAP.coords$Cell <- rownames(UMAP.coords)
    df_first_cell <- merge(df_first_cell,UMAP.coords, all.x=T, by="Cell")
    df_second_cell <- merge(df_second_cell,UMAP.coords, all.x=T, by="Cell")
    
    df_first_cell$clonecolor <- col_vector[i]
    df_second_cell$clonecolor <- col_vector[i]
    
    df_final <- rbind(df_final,df_first_cell,df_second_cell)
    
    i <- i+1
  }
  
  Dimplot.colors <- c(yarrr::transparent("#a40000", trans.val = .7),yarrr::transparent("#16317d", trans.val = .7),
                      yarrr::transparent("#007e2f", trans.val = .7),yarrr::transparent("#ffcd12", trans.val = .7),
                      yarrr::transparent("#b86092", trans.val = .7),yarrr::transparent("#721b3e", trans.val = .7),
                      yarrr::transparent("#00b7a7", trans.val = .7))
  
  p <- DimPlot(Seurat.Object, reduction = reduction,group.by = DimplotGroup, cols = Dimplot.colors[1:length(unique(Seurat.Object@meta.data[,DimplotGroup]))])+
    labs(color = DimplotGroup)
  p+geom_path(data=df_final, aes(x=wnnUMAP_1, y=wnnUMAP_2, group=group),col=df_final$clonecolor,size = 0.8)+ggtitle("Clonal connections")
  
}
Show.clones.on.UMAP.one.color <- function(Seurat.Object, clone_IDs, reduction,DimplotGroup,connection.color){
  # Getting the cells which will be connected
  Cells.to.connect <- Seurat.Object@meta.data$Row.names[Seurat.Object@meta.data$clone_id %in% clone_IDs& !is.na(Seurat.Object@meta.data$clone_id)]
  
  #Getting the UMAP coordinates
  UMAP.coords <- Seurat.Object[["wnn.umap"]]@cell.embeddings
  UMAP.coords <- UMAP.coords[rownames(UMAP.coords) %in% Cells.to.connect,]
  UMAP.coords <- as.data.frame(UMAP.coords)
  
  # Preparing the information of which cells need to be connected with each other
  df.start <- Seurat.Object@meta.data[Cells.to.connect,c("Tissue","Row.names","clone_id")]
  
  # How many arrows will be needed? I want an arrow going from each early timepoint to each later timepoint
  # Therefore the total number of arrows will be (per clone) :
  # n cells in timeopint 1 * n cells in timepoint 2 + n cells in timepoint * n cells in timepoint 3 + n cells in timepoint 3 * n cells in timepoint 4
  # If several clones are given as input, we have to do this for each clone individually
  # Therefore, I will go through the clones one by one.
  df_final <- data.frame()
  col_vector = c("brown","gray21", "purple",
                 "gray46","#FDBF6F","khaki2","maroon","orchid1","deeppink1","blue1","steelblue4",
                 "darkturquoise", "green1", "yellow4", "yellow3","darkorange4", "brown")
  
  i <- 1
  for (i in 1:length(unique(df.start$clone_id))) {
    this.rounds.clone <- unique(df.start$clone_id)[i]
    df.subset <- df.start[df.start$clone_id==this.rounds.clone,]
    
    # Now create all arrows one by one
    # Separate the df.subset into all timepoints
    df.subset_1 <- df.subset[df.subset$Tissue=="PBMC",]
    df.subset_2 <- df.subset[df.subset$Tissue=="Tonsil",]
    
    matrix = matrix(nrow = nrow(df.subset_1)*nrow(df.subset_2), ncol = 2)
    df.one = data.frame(matrix) 
    
    if(nrow(df.one)>0){
      c <- 1
      for (c in 1:nrow(df.subset_1)) {
        df.one[(c*nrow(df.subset_2)-(nrow(df.subset_2)-1)):((c*nrow(df.subset_2)-(nrow(df.subset_2)-1))+nrow(df.subset_2)-1),1] <- df.subset_1$Row.names[c]
        df.one[(c*nrow(df.subset_2)-(nrow(df.subset_2)-1)):((c*nrow(df.subset_2)-(nrow(df.subset_2)-1))+nrow(df.subset_2)-1),2] <- df.subset_2$Row.names
      }
    }
    
    df <- df.one
    df <- drop_na(df)
    
    # Preparing a dataframe that will be given as inputs to geom_line.
    # The dataframe will have three columns. Column x for all the x coordinates, column y for all the y coordinates and column z for the grouping.
    # Each arrow gets its own group. 
    # Adding the groups at this point:
    df$group <- paste0(i,"_",1:nrow(df))
    df$group <- as.character(df$group)
    
    # Then preparing the final dataframe as input for geomline step by step
    colnames(df) <- c("First_cell","Second_cell","group")
    df_first_cell <- df[,c(1,3)]
    df_second_cell <- df[,c(2,3)]
    colnames(df_first_cell) <- c("Cell","group")
    colnames(df_second_cell) <- c("Cell","group")
    UMAP.coords$Cell <- rownames(UMAP.coords)
    df_first_cell <- merge(df_first_cell,UMAP.coords, all.x=T, by="Cell")
    df_second_cell <- merge(df_second_cell,UMAP.coords, all.x=T, by="Cell")
    
    df_first_cell$clonecolor <- col_vector[i]
    df_second_cell$clonecolor <- col_vector[i]
    
    df_final <- rbind(df_final,df_first_cell,df_second_cell)
    
    i <- i+1
  }
  
  Dimplot.colors <- c(yarrr::transparent("#a40000", trans.val = .7),yarrr::transparent("#16317d", trans.val = .7),
                      yarrr::transparent("#007e2f", trans.val = .7),yarrr::transparent("#ffcd12", trans.val = .7),
                      yarrr::transparent("#b86092", trans.val = .7),yarrr::transparent("#721b3e", trans.val = .7),
                      yarrr::transparent("#00b7a7", trans.val = .7))
  
  p <- DimPlot(Seurat.Object, reduction = reduction,group.by = DimplotGroup, cols = Dimplot.colors[1:length(unique(Seurat.Object@meta.data[,DimplotGroup]))])+
    labs(color = DimplotGroup)
  p+geom_path(data=df_final, aes(x=wnnUMAP_1, y=wnnUMAP_2, group=group),col=connection.color,size = 0.8)+ggtitle("Clonal connections")
  
}
Show.clones.on.UMAP.color.by.timepoint <- function(Seurat.Object, clone_IDs, reduction,DimplotGroup){
  # Getting the cells which will be connected
  Cells.to.connect <- Seurat.Object@meta.data$Row.names[Seurat.Object@meta.data$clone_id %in% clone_IDs& !is.na(Seurat.Object@meta.data$clone_id)]
  
  #Getting the UMAP coordinates
  UMAP.coords <- Seurat.Object[["wnn.umap"]]@cell.embeddings
  UMAP.coords <- UMAP.coords[rownames(UMAP.coords) %in% Cells.to.connect,]
  UMAP.coords <- as.data.frame(UMAP.coords)
  
  # Preparing the information of which cells need to be connected with each other
  df.start <- Seurat.Object@meta.data[Cells.to.connect,c("Timepoint","Row.names","clone_id")]
  
  # How many arrows will be needed? I want an arrow going from each early timepoint to each later timepoint
  # Therefore the total number of arrows will be (per clone) :
  # n cells in timeopint 1 * n cells in timepoint 2 + n cells in timepoint * n cells in timepoint 3 + n cells in timepoint 3 * n cells in timepoint 4
  # If several clones are given as input, we have to do this for each clone individually
  # Therefore, I will go through the clones one by one.
  df_final <- data.frame()
  col_vector = c("brown","gray21", "purple",
                 "gray46","#FDBF6F","khaki2","maroon","orchid1","deeppink1","blue1","steelblue4",
                 "darkturquoise", "green1", "yellow4", "yellow3","darkorange4", "brown")
  
  i <- 1
  for (i in 1:length(unique(df.start$clone_id))) {
    this.rounds.clone <- unique(df.start$clone_id)[i]
    df.subset <- df.start[df.start$clone_id==this.rounds.clone,]
    
    # Now create all arrows one by one
    # Separate the df.subset into all timepoints
    df.subset_1 <- df.subset[df.subset$Timepoint=="Post 2nd vacc",]
    df.subset_2 <- df.subset[df.subset$Timepoint=="6M post 2nd vacc",]
    df.subset_3 <- df.subset[df.subset$Timepoint=="Post boost vacc",]
    df.subset_4 <- df.subset[df.subset$Timepoint=="Post reinfection",]
    
    matrix = matrix(nrow = nrow(df.subset_1)*nrow(df.subset_2), ncol = 2)
    df.one = data.frame(matrix) 
    matrix = matrix(nrow = nrow(df.subset_2)*nrow(df.subset_3), ncol = 2)
    df.two = data.frame(matrix)
    matrix = matrix(nrow = nrow(df.subset_3)*nrow(df.subset_4), ncol = 2)
    df.three = data.frame(matrix)
    matrix = matrix(nrow = nrow(df.subset_1)*nrow(df.subset_3), ncol = 2)
    df.four = data.frame(matrix)
    matrix = matrix(nrow = nrow(df.subset_1)*nrow(df.subset_4), ncol = 2)
    df.five = data.frame(matrix)
    matrix = matrix(nrow = nrow(df.subset_2)*nrow(df.subset_4), ncol = 2)
    df.six = data.frame(matrix)
    
    if(nrow(df.one)>0){
      c <- 1
      for (c in 1:nrow(df.subset_1)) {
        df.one[(c*nrow(df.subset_2)-(nrow(df.subset_2)-1)):((c*nrow(df.subset_2)-(nrow(df.subset_2)-1))+nrow(df.subset_2)-1),1] <- df.subset_1$Row.names[c]
        df.one[(c*nrow(df.subset_2)-(nrow(df.subset_2)-1)):((c*nrow(df.subset_2)-(nrow(df.subset_2)-1))+nrow(df.subset_2)-1),2] <- df.subset_2$Row.names
        df.one$color <- col_vector[1]
      }
    }
    if(nrow(df.two)>0){
      c <- 1
      for (c in 1:nrow(df.subset_2)) {
        df.two[(c*nrow(df.subset_3)-(nrow(df.subset_3)-1)):((c*nrow(df.subset_3)-(nrow(df.subset_3)-1))+nrow(df.subset_3)-1),1] <- df.subset_2$Row.names[c]
        df.two[(c*nrow(df.subset_3)-(nrow(df.subset_3)-1)):((c*nrow(df.subset_3)-(nrow(df.subset_3)-1))+nrow(df.subset_3)-1),2] <- df.subset_3$Row.names
        df.two$color <- col_vector[2]
      }
    }
    if(nrow(df.three)>0){
      c <- 1
      for (c in 1:nrow(df.subset_3)) {
        df.three[(c*nrow(df.subset_4)-(nrow(df.subset_4)-1)):((c*nrow(df.subset_4)-(nrow(df.subset_4)-1))+nrow(df.subset_4)-1),1] <- df.subset_3$Row.names[c]
        df.three[(c*nrow(df.subset_4)-(nrow(df.subset_4)-1)):((c*nrow(df.subset_4)-(nrow(df.subset_4)-1))+nrow(df.subset_4)-1),2] <- df.subset_4$Row.names
        df.three$color <- col_vector[3]
      }
    }
    if(nrow(df.four)>0&nrow(df.one)==0){
      c <- 1
      for (c in 1:nrow(df.subset_1)) {
        df.four[(c*nrow(df.subset_3)-(nrow(df.subset_3)-1)):((c*nrow(df.subset_3)-(nrow(df.subset_3)-1))+nrow(df.subset_3)-1),1] <- df.subset_1$Row.names[c]
        df.four[(c*nrow(df.subset_3)-(nrow(df.subset_3)-1)):((c*nrow(df.subset_3)-(nrow(df.subset_3)-1))+nrow(df.subset_3)-1),2] <- df.subset_3$Row.names
        df.four$color <- col_vector[4]
      }
    }
    if(nrow(df.five)>0&nrow(df.one)==0&nrow(df.four)==0&nrow(df.two)==0){
      c <- 1
      for (c in 1:nrow(df.subset_1)) {
        df.five[(c*nrow(df.subset_4)-(nrow(df.subset_4)-1)):((c*nrow(df.subset_4)-(nrow(df.subset_4)-1))+nrow(df.subset_4)-1),1] <- df.subset_1$Row.names[c]
        df.five[(c*nrow(df.subset_4)-(nrow(df.subset_4)-1)):((c*nrow(df.subset_4)-(nrow(df.subset_4)-1))+nrow(df.subset_4)-1),2] <- df.subset_4$Row.names
        df.five$color <- col_vector[5]
      }
    }
    if(nrow(df.six)>0&nrow(df.two)==0){
      c <- 1
      for (c in 1:nrow(df.subset_2)) {
        df.six[(c*nrow(df.subset_4)-(nrow(df.subset_4)-1)):((c*nrow(df.subset_4)-(nrow(df.subset_4)-1))+nrow(df.subset_4)-1),1] <- df.subset_2$Row.names[c]
        df.six[(c*nrow(df.subset_4)-(nrow(df.subset_4)-1)):((c*nrow(df.subset_4)-(nrow(df.subset_4)-1))+nrow(df.subset_4)-1),2] <- df.subset_4$Row.names
        df.six$color <- col_vector[6]
      }
    }
    
    if(nrow(df.one)==0|
       is.na(df.one$X1[1])){
      df.one <- data.frame(X1 = NA,X2=NA,color=NA)
    }
    if(nrow(df.two)==0|
       is.na(df.two$X1[1])){
      df.two <- data.frame(X1 = NA,X2=NA,color=NA)
    }
    if(nrow(df.three)==0|
       is.na(df.three$X1[1])){
      df.three <- data.frame(X1 = NA,X2=NA,color=NA)
    }
    if(nrow(df.four)==0|
       is.na(df.four$X1[1])){
      df.four <- data.frame(X1 = NA,X2=NA,color=NA)
    }
    if(nrow(df.five)==0|
       is.na(df.five$X1[1])){
      df.five <- data.frame(X1 = NA,X2=NA,color=NA)
    }
    if(nrow(df.six)==0 |
       is.na(df.six$X1[1])){
      df.six <- data.frame(X1 = NA,X2=NA,color=NA)
    }
    df <- rbind(df.one,df.two,df.three,df.four,df.five,df.six)
    df <- drop_na(df)
    
    # Preparing a dataframe that will be given as inputs to geom_line.
    # The dataframe will have three columns. Column x for all the x coordinates, column y for all the y coordinates and column z for the grouping.
    # Each arrow gets its own group. 
    # Adding the groups at this point:
    df$group <- paste0(i,"_",1:nrow(df))
    df$group <- as.character(df$group)
    
    # Then preparing the final dataframe as input for geomline step by step
    colnames(df) <- c("First_cell","Second_cell","color","group")
    df_first_cell <- df[,c(1,3,4)]
    df_second_cell <- df[,c(2,3,4)]
    colnames(df_first_cell) <- c("Cell","color","group")
    colnames(df_second_cell) <- c("Cell","color","group")
    UMAP.coords$Cell <- rownames(UMAP.coords)
    df_first_cell <- merge(df_first_cell,UMAP.coords, all.x=T, by="Cell")
    df_second_cell <- merge(df_second_cell,UMAP.coords, all.x=T, by="Cell")
    
    df_final <- rbind(df_final,df_first_cell,df_second_cell)
    
    i <- i+1
  }
  
  Dimplot.colors <- c(yarrr::transparent("#a40000", trans.val = .8),yarrr::transparent("#16317d", trans.val = .8),
                      yarrr::transparent("#007e2f", trans.val = .8),yarrr::transparent("#ffcd12", trans.val = .8),
                      yarrr::transparent("#b86092", trans.val = .8),yarrr::transparent("#721b3e", trans.val = .8),
                      yarrr::transparent("#00b7a7", trans.val = .8))
  
  p <- DimPlot(Seurat.Object, reduction = reduction,group.by = DimplotGroup,cols = Dimplot.colors[1:length(unique(Seurat.Object@meta.data[,DimplotGroup]))])
  p+geom_path(data=df_final, aes(x=wnnUMAP_1, y=wnnUMAP_2, group=group),col=df_final$color,size = 0.8,arrow=arrow(type = "open", angle = 30, length = unit(0.1, "inches")))+ggtitle("Clonal paths")
  
}
Cutoff_function <- function(Seurat_object, bait, left_border, right_border,melted_df){
  
  melted_df <- melted_df[melted_df$variable==bait,]
  p <- ggplot(melted_df, aes(x = value, y = variable,fill =variable)) +
    geom_density_ridges(scale = 4)+
    scale_x_continuous(trans='log10',breaks = c(1,5,10,20,40,100,200,500,700,1000,1500,3000,5000))+
    labs(title = "Density plot") +
    theme_ridges(font_size = 13, grid = TRUE) + 
    theme(axis.title.y = element_blank(),axis.title.x = element_blank())
  
  plot.data <-ggplot_build(p)
  plot.data <- as.data.frame(plot.data$data)
  plot.data <- plot.data[plot.data$x > log10(left_border)&
                           plot.data$x < log10(right_border),]
  
  intersect <- round(10^plot.data$x[plot.data$density==min(plot.data$density)],0)
  
  p <- ggplot(melted_df, aes(x = value, y = variable,fill =variable)) +
    geom_density_ridges(scale = 4)+
    scale_x_continuous(trans='log10',breaks = c(1,5,10,20,40,100,200,500,700,1000,1500,3000,5000))+
    labs(title = paste("Density plot",bait)) +
    theme_ridges(font_size = 13, grid = TRUE) + 
    theme(axis.title.y = element_blank(),axis.title.x = element_blank(),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    geom_vline(xintercept=intersect)
  
  print(p)
  print(paste("The cutoff for",bait, "is",intersect))
  positive_cells <- names(Seurat_object@assays$Baiting@counts[bait,Seurat_object@assays$Baiting@counts[bait,]> intersect])
  metadata_column <- paste0(bait,"_specific")
  Seurat_object@meta.data[rownames(Seurat_object@meta.data) %in% positive_cells,metadata_column] <- "yes"
  
  return(Seurat_object)
}
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
} # Used to plot significance asterics for waffle plot.
`%notin%` <- Negate(`%in%`) # Creating a useful operator

######################################################################################################################
# Part 1: Loading the data sets, filtering, Hashing demultiplexing, addition of metadata columns
######################################################################################################################
# Loading the first dataset into R
setwd("~/NAS/Jan/Experiments and Data/210801 - Exp 035 - Profiling of COVID vaccinated patient samples/Memory B cell Fifth Experiment 221115/multi_exp035_revision_2_S1")  
Tonsil_dataset <- "./outs/per_sample_outs/multi_exp035_revision_2_S1/count/sample_feature_bc_matrix"
Tonsil_dataset <- Read10X(data.dir = Tonsil_dataset)

# Loading the second dataset into R
setwd("~/NAS/Jan/Experiments and Data/210801 - Exp 035 - Profiling of COVID vaccinated patient samples/Memory B cell Fifth Experiment 221115/multi_exp035_revision_2_S2")  
Tonsil_dataset_2 <- "./outs/per_sample_outs/multi_exp035_revision_2_S2/count/sample_feature_bc_matrix"
Tonsil_dataset_2 <- Read10X(data.dir = Tonsil_dataset_2)

# Seurat object is created
Tonsil_1 <- CreateSeuratObject(counts = Tonsil_dataset$`Gene Expression`)
Tonsil_2 <- CreateSeuratObject(counts = Tonsil_dataset_2$`Gene Expression`)

# Assays for Baiting constructs, Protein quantification and hashing are created
Baiting_assay <- CreateAssayObject(counts = Tonsil_dataset$`Antibody Capture`[c(9:13),])
Baiting_assay_2 <- CreateAssayObject(counts = Tonsil_dataset_2$`Antibody Capture`[c(9:13),])
Protein_assay <- CreateAssayObject(counts = Tonsil_dataset$`Antibody Capture`[c(14:25),])
Protein_assay_2 <- CreateAssayObject(counts = Tonsil_dataset_2$`Antibody Capture`[c(14:25),])
Hashing_assay <- CreateAssayObject(counts = Tonsil_dataset$`Antibody Capture`[c(1:8),])
Hashing_assay_2 <- CreateAssayObject(counts = Tonsil_dataset_2$`Antibody Capture`[c(1:8),])

# Now the assays are added to the previously created Seurat object
Tonsil_1[["Baiting"]] <- Baiting_assay
Tonsil_1[["Protein"]] <- Protein_assay
Tonsil_1[["Hashing"]] <- Hashing_assay
Tonsil_2[["Baiting"]] <- Baiting_assay_2
Tonsil_2[["Protein"]] <- Protein_assay_2
Tonsil_2[["Hashing"]] <- Hashing_assay_2

# Removal of assay objects
remove(Baiting_assay,Hashing_assay,Protein_assay,Tonsil_dataset,Baiting_assay_2,Hashing_assay_2,Protein_assay_2,Tonsil_dataset_2)

# QC and selecting cells for further analysis - first dataset
Tonsil_1[["percent.mt"]] <- PercentageFeatureSet(Tonsil_1, pattern = "^MT-")
VlnPlot(Tonsil_1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
quantile(Tonsil_1@meta.data$nFeature_RNA,probs = seq(0, 1, 1/20))
quantile(Tonsil_1@meta.data$percent.mt,probs = seq(0, 1, 1/20))
nrow(Tonsil_1@meta.data[Tonsil_1@meta.data$percent.mt < 7.5 &
                          Tonsil_1@meta.data$nFeature_RNA > 200 &
                          Tonsil_1@meta.data$nFeature_RNA < 4000,])/nrow(Tonsil_1@meta.data)

Tonsil_1 <- subset(Tonsil_1, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 7.5)

Tonsil_1 <- NormalizeData(Tonsil_1)

# QC and selecting cells for further analysis - second dataset
Tonsil_2[["percent.mt"]] <- PercentageFeatureSet(Tonsil_2, pattern = "^MT-")
VlnPlot(Tonsil_2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
quantile(Tonsil_2@meta.data$nFeature_RNA,probs = seq(0, 1, 1/20))
quantile(Tonsil_2@meta.data$percent.mt,probs = seq(0, 1, 1/20))
nrow(Tonsil_2@meta.data[Tonsil_2@meta.data$percent.mt < 7.5 &
                          Tonsil_2@meta.data$nFeature_RNA > 200 &
                          Tonsil_2@meta.data$nFeature_RNA < 4000,])/nrow(Tonsil_2@meta.data)

Tonsil_2 <- subset(Tonsil_2, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 7.5)

Tonsil_2 <- NormalizeData(Tonsil_2)

# Demultiplexing the HTO data
Tonsil_1 <- NormalizeData(Tonsil_1, assay = "Hashing", normalization.method = "CLR", margin = 2)
Tonsil_2 <- NormalizeData(Tonsil_2, assay = "Hashing", normalization.method = "CLR", margin = 2)
Tonsil_1 <- HTODemux(Tonsil_1, assay = "Hashing", positive.quantile = 0.99)
Tonsil_2 <- HTODemux(Tonsil_2, assay = "Hashing", positive.quantile = 0.99)
table(Tonsil_1@meta.data$Hashing_classification.global)
table(Tonsil_2@meta.data$Hashing_classification.global)
Tonsil_1 <- subset(Tonsil_1, subset = Hashing_classification.global=="Singlet")
Tonsil_2 <- subset(Tonsil_2, subset = Hashing_classification.global=="Singlet")
RidgePlot(Tonsil_1, assay = "Hashing", features = rownames(Tonsil_1[["Hashing"]])[1:10], ncol = 2)
RidgePlot(Tonsil_2, assay = "Hashing", features = rownames(Tonsil_2[["Hashing"]])[1:10], ncol = 2)
table(Tonsil_1@meta.data$Hashing_classification, useNA="always")
table(Tonsil_2@meta.data$Hashing_classification, useNA="always")

# Introducing metadata columns for Patient, Tissue, Vaccinated, Reconvalescent
Tonsil_1@meta.data$Patient <- "None"
Tonsil_2@meta.data$Patient <- "None"
Tonsil_1@meta.data$Patient[grep("243651", Tonsil_1@meta.data$Hashing_classification)] <- "243651-2001"
Tonsil_2@meta.data$Patient[grep("243651", Tonsil_2@meta.data$Hashing_classification)] <- "243651-2001"
Tonsil_1@meta.data$Patient[grep("253340", Tonsil_1@meta.data$Hashing_classification)] <- "253340-1997"
Tonsil_2@meta.data$Patient[grep("253340", Tonsil_2@meta.data$Hashing_classification)] <- "253340-1997"
Tonsil_1@meta.data$Patient[grep("254869", Tonsil_1@meta.data$Hashing_classification)] <- "254869-1996"
Tonsil_2@meta.data$Patient[grep("254869", Tonsil_2@meta.data$Hashing_classification)] <- "254869-1996"
Tonsil_1@meta.data$Patient[grep("277926", Tonsil_1@meta.data$Hashing_classification)] <- "277926-1975"
Tonsil_2@meta.data$Patient[grep("277926", Tonsil_2@meta.data$Hashing_classification)] <- "277926-1975"

Tonsil_1@meta.data$Tissue <- "None"
Tonsil_2@meta.data$Tissue <- "None"
Tonsil_1@meta.data$Tissue[grep("Tonsil",Tonsil_1@meta.data$Hashing_classification)] <- "Tonsil"
Tonsil_2@meta.data$Tissue[grep("Tonsil",Tonsil_2@meta.data$Hashing_classification)] <- "Tonsil"
Tonsil_1@meta.data$Tissue[grep("PBMC",Tonsil_1@meta.data$Hashing_classification)] <- "PBMC"
Tonsil_2@meta.data$Tissue[grep("PBMC",Tonsil_2@meta.data$Hashing_classification)] <- "PBMC"

Tonsil_1@meta.data$Vaccinated <- "no"
Tonsil_2@meta.data$Vaccinated <- "no"
Tonsil_1@meta.data$Vaccinated[grep("Vac",Tonsil_1@meta.data$Hashing_classification)] <- "yes"
Tonsil_2@meta.data$Vaccinated[grep("Vac",Tonsil_2@meta.data$Hashing_classification)] <- "yes"

Tonsil_1@meta.data$Reconvalescent <- "no"
Tonsil_2@meta.data$Reconvalescent <- "no"
Tonsil_1@meta.data$Reconvalescent[grep("Rec",Tonsil_1@meta.data$Hashing_classification)] <- "yes"
Tonsil_2@meta.data$Reconvalescent[grep("Rec",Tonsil_2@meta.data$Hashing_classification)] <- "yes"

######################################################################################################################
# Part 2: Processing of baiting counts (cutoffs, normalization) and visualization
######################################################################################################################
# Adding metadata columns to the Seurat object which will indicate if a cell is binding an antigen
Tonsil_1@meta.data$wtSpike_specific <- "no"
Tonsil_1@meta.data$RBD_specific <- "no"
Tonsil_1@meta.data$B.1.529_specific <- "no"
Tonsil_1@meta.data$HE_specific <- "no"
Tonsil_2@meta.data$wtSpike_specific <- "no"
Tonsil_2@meta.data$RBD_specific <- "no"
Tonsil_2@meta.data$B.1.529_specific <- "no"
Tonsil_2@meta.data$HE_specific <- "no"

# Pre-processing baiting counts
Baiting_df <- as.data.frame(Tonsil_1@assays$Baiting@counts)
Baiting_df <- t(Baiting_df)
Baiting_df <- as.data.frame(Baiting_df)
Baiting_df_2 <- as.data.frame(Tonsil_2@assays$Baiting@counts)
Baiting_df_2 <- t(Baiting_df_2)
Baiting_df_2 <- as.data.frame(Baiting_df_2)

# Subtracting the negative control counts
Baiting_df$wtSpike <- Baiting_df$wtSpike-Baiting_df$NegControl
Baiting_df$B.1.529 <- Baiting_df$B.1.529-Baiting_df$NegControl
Baiting_df$RBD <- Baiting_df$RBD-Baiting_df$NegControl
Baiting_df$HE <- Baiting_df$HE-Baiting_df$NegControl
Baiting_df_2$wtSpike <- Baiting_df_2$wtSpike-Baiting_df_2$NegControl
Baiting_df_2$B.1.529 <- Baiting_df_2$B.1.529-Baiting_df_2$NegControl
Baiting_df_2$RBD <- Baiting_df_2$RBD-Baiting_df_2$NegControl
Baiting_df_2$HE <- Baiting_df_2$HE-Baiting_df_2$NegControl

# Setting negative bait counts to 0
Baiting_df$wtSpike[Baiting_df$wtSpike<0] <- 0
Baiting_df$B.1.529[Baiting_df$B.1.529<0] <- 0
Baiting_df$HE[Baiting_df$HE<0] <- 0
Baiting_df$RBD[Baiting_df$RBD<0] <- 0
Baiting_df_2$wtSpike[Baiting_df_2$wtSpike<0] <- 0
Baiting_df_2$B.1.529[Baiting_df_2$B.1.529<0] <- 0
Baiting_df_2$HE[Baiting_df_2$HE<0] <- 0
Baiting_df_2$RBD[Baiting_df_2$RBD<0] <- 0

melted_Baiting_df <- reshape2::melt(Baiting_df)
melted_Baiting_df <- melted_Baiting_df[melted_Baiting_df$value>0,]
melted_Baiting_df_2 <- reshape2::melt(Baiting_df_2)
melted_Baiting_df_2 <- melted_Baiting_df_2[melted_Baiting_df_2$value>0,]

# Plot the data and select an x axis window, in which the cutoff should be determined.
ggplot(melted_Baiting_df, aes(x = value, y = variable,fill =variable)) +
  geom_density_ridges(scale = 4)+
  scale_x_continuous(trans='log10',breaks = c(1,5,10,20,40,100,200,500,700,1000,1500,3000,5000))+
  labs(title = "Density plot") +
  theme_ridges(font_size = 13, grid = TRUE) + 
  theme(axis.title.y = element_blank(),axis.title.x = element_blank(),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggplot(melted_Baiting_df_2, aes(x = value, y = variable,fill =variable)) +
  geom_density_ridges(scale = 4)+
  scale_x_continuous(trans='log10',breaks = c(1,5,10,20,40,100,200,500,700,1000,1500,3000,5000))+
  labs(title = "Density plot") +
  theme_ridges(font_size = 13, grid = TRUE) + 
  theme(axis.title.y = element_blank(),axis.title.x = element_blank(),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# A function that determines the cutoff and plots density chart including the cutoff line
Tonsil_1 <- Cutoff_function(Seurat_object = Tonsil_1, bait="wtSpike",left_border = 1,right_border = 7,melted_df = melted_Baiting_df)
Tonsil_1 <- Cutoff_function(Seurat_object = Tonsil_1, bait="RBD",left_border = 5,right_border = 100,melted_df = melted_Baiting_df)
Tonsil_2 <- Cutoff_function(Seurat_object = Tonsil_2, bait="wtSpike",left_border = 1,right_border = 7,melted_df = melted_Baiting_df)
Tonsil_2 <- Cutoff_function(Seurat_object = Tonsil_2, bait="RBD",left_border = 5,right_border = 100,melted_df = melted_Baiting_df)

remove(melted_Baiting_df,melted_Baiting_df_2)

# Setting baiting counts below determined thresholds to NA
Baiting_df$wtSpike[Baiting_df$wtSpike<5] <- NA
Baiting_df$RBD[Baiting_df$RBD<52] <- NA

Baiting_df_2$wtSpike[Baiting_df_2$wtSpike<5] <- NA
Baiting_df_2$RBD[Baiting_df_2$RBD<52] <- NA

# For now, subset Baiting_dfs to only wt full length and RBD
Baiting_df <- Baiting_df[,c("wtSpike","RBD")]
Baiting_df_2 <- Baiting_df_2[,c("wtSpike","RBD")]

# In order to use the Seurat Normalization function for corrected baiting counts, we need to have the corrected counts as assay object
# For this, the Baiting_df needs to be transposed
Baiting_df <- as.data.frame(t(Baiting_df))
Baiting_df_2 <- as.data.frame(t(Baiting_df_2))
Cor_Baiting_assay <- CreateAssayObject(counts = Baiting_df)
Cor_Baiting_assay_2 <- CreateAssayObject(counts = Baiting_df_2)
Tonsil_1[["Cor_Baiting"]] <- Cor_Baiting_assay
Tonsil_2[["Cor_Baiting"]] <- Cor_Baiting_assay_2
remove(Cor_Baiting_assay,Cor_Baiting_assay_2)
Tonsil_1 <- NormalizeData(Tonsil_1, assay = "Cor_Baiting", normalization.method = 'CLR', margin = 1)
Tonsil_2 <- NormalizeData(Tonsil_2, assay = "Cor_Baiting", normalization.method = 'CLR', margin = 1)
normalized_across_features <- as.data.frame(t(Tonsil_1@assays$Cor_Baiting@data))
normalized_across_features_2 <- as.data.frame(t(Tonsil_2@assays$Cor_Baiting@data))
remove(Baiting_df,Baiting_df_2)

# For scaling, I need to transform again, then each column is scaled to 0:1 individually using my own scaling function.
normalized_across_features <- myscaling(normalized_across_features)  
normalized_across_features_2 <- myscaling(normalized_across_features_2)  

# These are the final LIBRA scores
# Adding this to the Seurat object
Tonsil_1@meta.data <- merge(Tonsil_1@meta.data,normalized_across_features, by=0)
Tonsil_2@meta.data <- merge(Tonsil_2@meta.data,normalized_across_features_2, by=0)
rownames(Tonsil_1@meta.data) <- Tonsil_1@meta.data$Row.names
rownames(Tonsil_2@meta.data) <- Tonsil_2@meta.data$Row.names
remove(normalized_across_features,normalized_across_features_2)

# NAs are now set to 0
Tonsil_1@meta.data$wtSpike[is.na(Tonsil_1@meta.data$wtSpike)] <- 0
Tonsil_1@meta.data$RBD[is.na(Tonsil_1@meta.data$RBD)] <- 0
Tonsil_2@meta.data$wtSpike[is.na(Tonsil_2@meta.data$wtSpike)] <- 0
Tonsil_2@meta.data$RBD[is.na(Tonsil_2@meta.data$RBD)] <- 0

# FeatureScatter plots for first dataset
FeatureScatter(Tonsil_1,  feature1 = "wtSpike", feature2 = "RBD", pt.size = 1, group.by = "Patient")+
  ggtitle("Wt Spike and RBD scores")+xlab("wt-Spike scores")+ylab("RBD scores")

# FeatureScatter plots for first dataset
FeatureScatter(Tonsil_2,  feature1 = "wtSpike", feature2 = "RBD", pt.size = 1, group.by = "Patient")+
  ggtitle("Wt Spike and RBD scores")+xlab("wt-Spike scores")+ylab("RBD scores")


######################################################################################################################
# Part 3: Dataset integration
######################################################################################################################
# Renaming of cells
Tonsil_1 <- RenameCells(object = Tonsil_1, add.cell.id = "Dataset_1")
Tonsil_2 <- RenameCells(object = Tonsil_2, add.cell.id = "Dataset_2")
Tonsil_1@meta.data$Dataset <- "First"
Tonsil_2@meta.data$Dataset <- "Second"

# During the IntegrateData() function, dimensional reductions of the data are needed. PCA is calculated before Integration is done (only for the protein data)
DefaultAssay(Tonsil_1) <- "Protein"
DefaultAssay(Tonsil_2) <- "Protein"
VariableFeatures(Tonsil_1) <- rownames(Tonsil_1[["Protein"]])
VariableFeatures(Tonsil_2) <- rownames(Tonsil_2[["Protein"]])
Tonsil_1 <- NormalizeData(Tonsil_1, normalization.method = 'CLR', margin = 2) %>% ScaleData() %>% RunPCA(reduction.name = 'apca',approx=FALSE) # approx=F is used to remove the warning that comes otherwise
Tonsil_2 <- NormalizeData(Tonsil_2, normalization.method = 'CLR', margin = 2) %>% ScaleData() %>% RunPCA(reduction.name = 'apca',approx=FALSE)

# Create a list containing the two datasets
Tonsil_dataset_list <- list(Tonsil_1,Tonsil_2)

# Select features that are repeatedly variable across datasets for integration
# The integration is done independently for the RNA and the Protein assay.
features_RNA <- SelectIntegrationFeatures(object.list = Tonsil_dataset_list,assay = c("RNA","RNA")) # IG variable regions pop up!
features_Protein <- SelectIntegrationFeatures(object.list = Tonsil_dataset_list,assay = c("Protein","Protein"))

# Remove HLA, Immunoglobulin, RNA, MT, and RP genes based on HUGO gene names
var_regex = '^HLA-|^IG[HJKL]|^RNA|^MT|^RP'
features_RNA <- grep(var_regex, features_RNA, invert=T, value=T)
remove(var_regex)

# Perform integration. This command requires long time / a lot of computational resources to run.
# Because of only 6 protein features, dims is set to 5 for the anchor finding of protein data
anchors_RNA <- FindIntegrationAnchors(object.list = Tonsil_dataset_list, anchor.features = features_RNA,assay = c("RNA","RNA"))
anchors_Protein <- FindIntegrationAnchors(object.list = Tonsil_dataset_list, anchor.features = features_Protein,assay = c("Protein","Protein"))

# This command creates an 'integrated' data assay in a new Seurat object.
Tonsil.combined.RNA <- IntegrateData(anchorset = anchors_RNA,new.assay.name = "integratedRNA")
Tonsil.combined.RNA <- FindVariableFeatures(Tonsil.combined.RNA, assay = "Protein")
Tonsil.combined.RNA <- ScaleData(Tonsil.combined.RNA, assay = "Protein")
Tonsil.combined.RNA <- RunPCA(Tonsil.combined.RNA, assay = "Protein", reduction.name="pca_protein_all_cells",approx=FALSE)
Tonsil.combined.Protein <- IntegrateEmbeddings(anchorset = anchors_Protein, reductions = Tonsil.combined.RNA@reductions$pca_protein_all_cells)
remove(anchors_Protein,anchors_RNA,Tonsil_dataset_list,Tonsil_1,Tonsil_2)

Tonsil <- Tonsil.combined.RNA
remove(Tonsil.combined.RNA)

######################################################################################################################
# Part 4: Creating Immcantation Input files (all cells)
######################################################################################################################
all.cells <- rownames(Tonsil@meta.data)

# Changing the barcodes
# AA -> dataset1
# GG -> dataset2
Immcantation.Input.Barcode.Function <- function(Barcode.Vector){
  i <- 1
  for (i in 1:length(Barcode.Vector)) {
    if(substr(Barcode.Vector[i],9,9)=="1"){
      Barcode.Vector[i] <- paste0("AA",substr(Barcode.Vector[i],11,28))
    }
    if(substr(Barcode.Vector[i],9,9)=="2"){
      Barcode.Vector[i] <- paste0("GG",substr(Barcode.Vector[i],11,28))
    }
    i <- i+1
  } 
  return(Barcode.Vector)
}

all.cells <- Immcantation.Input.Barcode.Function(all.cells)

write(all.cells, file = "~/NAS/Jan/Experiments and Data/210801 - Exp 035 - Profiling of COVID vaccinated patient samples/Memory B cell Fifth Experiment 221115/Immcantation/Input_data/all.cells.txt")

# Prepare the contig.csv file per patient
contigs_Dataset_1 <- read.csv("~/NAS/Jan/Experiments and Data/210801 - Exp 035 - Profiling of COVID vaccinated patient samples/Memory B cell Fifth Experiment 221115/multi_exp035_revision_2_S1/outs/per_sample_outs/multi_exp035_revision_2_S1/vdj_b/filtered_contig_annotations.csv")
contigs_Dataset_2 <- read.csv("~/NAS/Jan/Experiments and Data/210801 - Exp 035 - Profiling of COVID vaccinated patient samples/Memory B cell Fifth Experiment 221115/multi_exp035_revision_2_S2/outs/per_sample_outs/multi_exp035_revision_2_S2/vdj_b/filtered_contig_annotations.csv")

# Changing barcodes
contigs_Dataset_1$barcode <- paste0("AA",contigs_Dataset_1$barcode)
contigs_Dataset_2$barcode <- paste0("GG",contigs_Dataset_2$barcode)

# Changing contig_id
contigs_Dataset_1$contig_id <- paste0("AA",contigs_Dataset_1$contig_id)
contigs_Dataset_2$contig_id <- paste0("GG",contigs_Dataset_2$contig_id)

# Merging
all.contigs <- rbind(contigs_Dataset_1,contigs_Dataset_2)
remove(contigs_Dataset_1,contigs_Dataset_2)

# Creating contig objects
contigs_all_cells <- all.contigs[all.contigs$barcode %in% all.cells,]

write_csv(contigs_all_cells, file = "~/NAS/Jan/Experiments and Data/210801 - Exp 035 - Profiling of COVID vaccinated patient samples/Memory B cell Fifth Experiment 221115/Immcantation/Input_data/all.cells.csv")


######################################################################################################################
# Part 5: WNN
######################################################################################################################
# Now WNN
# For this, I need to transfer the integrated protein dim reduction from the Bmem.combined.Protein object to the Bmem.combined.RNA object.
# But first, I am preforming data scaling and PCA on the combined RNA project.
DefaultAssay(Tonsil)
DefaultAssay(Tonsil.combined.Protein)

Tonsil <- ScaleData(Tonsil, verbose = FALSE)
Tonsil <- RunPCA(Tonsil, npcs = 30, verbose = FALSE)

# Now we transfer the integrated protein dim reduction
Tonsil@reductions$IntegrateEmbeddings.pca <- Tonsil.combined.Protein@reductions$integrated_dr

# I double check that the correct assays are associated with the correct PCAs.
Tonsil[['pca']]@assay.used
Tonsil[['IntegrateEmbeddings.pca']]@assay.used

# Finally, the FindMultiModalNeighbors function is called, using the PCAs from above.
ElbowPlot(Tonsil, reduction="pca")
Tonsil <- FindMultiModalNeighbors(
  Tonsil, reduction.list = list("pca", "IntegrateEmbeddings.pca"), 
  dims.list = list(1:15, 1:12), modality.weight.name = c("RNA.weight","Protein.weight"))


Tonsil@meta.data$Full.Row.names <- rownames(Tonsil@meta.data)
DefaultAssay(Tonsil) <- "RNA"
Tonsil <- ScaleData(Tonsil)
remove(features_Protein,features_RNA,Tonsil.combined.Protein)


######################################################################################################################
# Part 6.1: Clustering and UMAP visualization before subsetting
######################################################################################################################
# First clustering, UMAP and Dimplot on RNA data only
Tonsil <- FindNeighbors(Tonsil, dims = 1:15,graph.name = "integratedRNA_nn",assay="integratedRNA",reduction = "pca")
Tonsil <- FindClusters(Tonsil, graph.name = "integratedRNA_nn", algorithm = 1, verbose = FALSE,resolution = 0.5)
Tonsil <- RunUMAP(Tonsil, reduction = 'pca', dims = 1:15, assay = 'integratedRNA', reduction.name = 'rna.umap', reduction.key = 'rnaUMAP_')

DimPlot(Tonsil, reduction = 'rna.umap', label = F, repel = TRUE, label.size = 2.5, pt.size = 0.2,group.by = "integratedRNA_nn_res.0.5", 
        cols = MetBrewer::met.brewer(name="Egypt", n=length(unique(Tonsil@meta.data$seurat_clusters))))+ggtitle("RNA UMAP")+center.title() + NoAxes() # Export 4x6


# Then clustering, UMAP and Dimplot on Protein data only
Tonsil <- FindNeighbors(Tonsil, dims = 1:12,graph.name = "integratedADT_nn",assay="Protein",reduction = "IntegrateEmbeddings.pca")
Tonsil <- FindClusters(Tonsil, graph.name = "integratedADT_nn", algorithm = 1, verbose = FALSE)
Tonsil <- RunUMAP(Tonsil, reduction = 'IntegrateEmbeddings.pca', dims = 1:12, assay = 'Protein', reduction.name = 'adt.umap', reduction.key = 'adtUMAP_')

DimPlot(Tonsil, reduction = 'adt.umap', label = F, repel = TRUE, label.size = 2.5, pt.size = 0.6, 
        cols = MetBrewer::met.brewer(name="Egypt", n=length(unique(Tonsil@meta.data$seurat_clusters))))+ggtitle("Protein UMAP")+center.title()


# Then on the WNN data - find multimodal neighbours done already in Part 3.
DefaultAssay(Tonsil)
Tonsil <- FindClusters(Tonsil, graph.name = "wsnn", algorithm = 1, verbose = FALSE,resolution = 0.4)
Tonsil <- RunUMAP(Tonsil, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")

DimPlot(Tonsil, reduction = 'wnn.umap', label = T, repel = TRUE, label.size = 4, pt.size = 0.5, 
        cols = MetBrewer::met.brewer(name="Egypt", n=length(unique(Tonsil@meta.data$seurat_clusters))))+
  theme(text = element_text(size = 20)) +ggtitle(paste("WNN UMAP"))+center.title()+NoAxes() #5x7


p <- DimPlot(Tonsil, reduction = 'wnn.umap', label = T, repel = TRUE, label.size = 2.5, pt.size = 0.5, 
             cols = MetBrewer::met.brewer(name="Egypt", n=length(unique(Tonsil@meta.data$seurat_clusters))))+
  theme(text = element_text(size = 20)) +ggtitle(paste("WNN UMAP"))+center.title()+
  DimPlot(Tonsil,reduction = "wnn.umap",group.by = "wtSpike_specific",
          cols = MetBrewer::met.brewer(name="Austria", n=2))+ 
  DimPlot(Tonsil,reduction = "wnn.umap",group.by = "Dataset", 
          cols = MetBrewer::met.brewer(name="Austria", n=2))+ 
  DimPlot(Tonsil, reduction = "wnn.umap", group.by = "Vaccinated", 
          cols = MetBrewer::met.brewer(name="Austria", n=2))+ 
  DimPlot(Tonsil, reduction = "wnn.umap", group.by = "Patient", 
          cols = MetBrewer::met.brewer(name="Austria", n=4))+ 
  DimPlot(Tonsil, reduction = "wnn.umap", group.by = "Reconvalescent", 
          cols = MetBrewer::met.brewer(name="Austria", n=2))+ 
  DimPlot(Tonsil, reduction = "wnn.umap", group.by = "Tissue", 
          cols = MetBrewer::met.brewer(name="Austria", n=2))

p & NoAxes()
remove(p)

Markers <- rev(c("IGHD","IGHM","IGHA1","IGHA2","IGHG1","CD19","MS4A1","CR2","CD24","CD27","CD38","BACH2","FCER2","TCL1A","IL4R",
                 "AICDA","BCL6","MKI67","PRDM1","XBP1","IRF4","TBX21","FCRL5","ITGAX","FCRL4","ENTPD1","CD69"))

Idents(Tonsil) <- "seurat_clusters"
DotPlot(Tonsil,assay = "RNA",features = Markers)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ylab("Cluster") # 4x9

DefaultAssay(Tonsil) <- "Protein"
FeaturePlot(Tonsil,features = rownames(Tonsil@assays$Protein),reduction = "wnn.umap")

VlnPlot(Tonsil, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


# Annotation Part 1
Tonsil@meta.data$Subset <- "Memory"
Tonsil@meta.data$Subset[Tonsil@meta.data$seurat_clusters==2|
                          Tonsil@meta.data$seurat_clusters==4] <-"Naive"
Tonsil@meta.data$Subset[Tonsil@meta.data$seurat_clusters==6] <- "Activated and Unswitched Naive/Memory"
Tonsil@meta.data$Subset[Tonsil@meta.data$seurat_clusters==8|
                          Tonsil@meta.data$seurat_clusters==10] <- "GC cells"
Tonsil@meta.data$Subset[Tonsil@meta.data$seurat_clusters==7] <- "PB/PC"
Tonsil@meta.data$Subset[Tonsil@meta.data$seurat_clusters==11] <- "Non B cells"

# For the following plots, I don't want to include the PB/PC cluster because it contains a lot of cells with bad quality.
Tonsil.subset <- subset(Tonsil, subset=Subset == "Naive"|
                          Subset=="Activated and Unswitched Naive/Memory"|
                          Subset=="GC cells"|
                          Subset=="Memory"|
                          Subset=="Non B cells")
Tonsil.subset@meta.data$seurat_clusters <- as.numeric(as.character(Tonsil.subset@meta.data$seurat_clusters))
Tonsil.subset@meta.data$seurat_clusters[Tonsil.subset@meta.data$seurat_clusters<8] <- Tonsil.subset@meta.data$seurat_clusters[Tonsil.subset@meta.data$seurat_clusters<8]+1
Tonsil.subset@meta.data$seurat_clusters <- factor(Tonsil.subset@meta.data$seurat_clusters)
Idents(Tonsil.subset) <- "seurat_clusters"
DimPlot(Tonsil.subset, reduction = 'wnn.umap', label = F, repel = TRUE, label.size = 2.5, pt.size = 0.5,group.by = "Subset",
        cols = MetBrewer::met.brewer(name="Egypt", n=length(unique(Tonsil@meta.data$Subset)))) + NoAxes() #5x9
DimPlot(Tonsil.subset, reduction = 'wnn.umap', label = F, repel = TRUE, label.size = 2.5, pt.size = 0.5,group.by = "seurat_clusters",
        cols = MetBrewer::met.brewer(name="Egypt", n=length(unique(Tonsil@meta.data$seurat_clusters)))) + NoAxes()+ggtitle("WNN Clusters") #5x9

Tonsil.subset@meta.data$Subset <- factor(Tonsil.subset@meta.data$Subset, levels = rev(c("Naive","Activated and Unswitched Naive/Memory","GC cells","Memory","Non B cells")))
Idents(Tonsil.subset)  <- "Subset"
DotPlot(Tonsil.subset,assay = "RNA",features = Markers)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ylab("Celltype") # 4x10.5

remove(Tonsil.subset,Markers)

# Getting the cell names of cells to exclude from the further analysis:
Kick <- rownames(Tonsil@meta.data[Tonsil@meta.data$Subset=="PB/PC"|
                            Tonsil@meta.data$Subset=="Naive"|
                            Tonsil@meta.data$Subset=="GC cells"|
                            Tonsil@meta.data$Subset=="Activated and Unswitched Naive/Memory"|
                              Tonsil@meta.data$Subset=="Myeloid",])


Tonsil.before.subset <- Tonsil

######################################################################################################################
# Part 6.2: Loading the raw data sets again, filtering, Hashing demultiplexing, addition of metadata columns
######################################################################################################################
# Loading the first dataset into R
setwd("~/NAS/Jan/Experiments and Data/210801 - Exp 035 - Profiling of COVID vaccinated patient samples/Memory B cell Fifth Experiment 221115/multi_exp035_revision_2_S1")  
Tonsil_dataset <- "./outs/per_sample_outs/multi_exp035_revision_2_S1/count/sample_feature_bc_matrix"
Tonsil_dataset <- Read10X(data.dir = Tonsil_dataset)

# Loading the second dataset into R
setwd("~/NAS/Jan/Experiments and Data/210801 - Exp 035 - Profiling of COVID vaccinated patient samples/Memory B cell Fifth Experiment 221115/multi_exp035_revision_2_S2")  
Tonsil_dataset_2 <- "./outs/per_sample_outs/multi_exp035_revision_2_S2/count/sample_feature_bc_matrix"
Tonsil_dataset_2 <- Read10X(data.dir = Tonsil_dataset_2)

# Seurat object is created
Tonsil_1 <- CreateSeuratObject(counts = Tonsil_dataset$`Gene Expression`)
Tonsil_2 <- CreateSeuratObject(counts = Tonsil_dataset_2$`Gene Expression`)

# Assays for Baiting constructs, Protein quantification and hashing are created
Baiting_assay <- CreateAssayObject(counts = Tonsil_dataset$`Antibody Capture`[c(9:13),])
Baiting_assay_2 <- CreateAssayObject(counts = Tonsil_dataset_2$`Antibody Capture`[c(9:13),])
Protein_assay <- CreateAssayObject(counts = Tonsil_dataset$`Antibody Capture`[c(14:25),])
Protein_assay_2 <- CreateAssayObject(counts = Tonsil_dataset_2$`Antibody Capture`[c(14:25),])
Hashing_assay <- CreateAssayObject(counts = Tonsil_dataset$`Antibody Capture`[c(1:8),])
Hashing_assay_2 <- CreateAssayObject(counts = Tonsil_dataset_2$`Antibody Capture`[c(1:8),])

# Now the assays are added to the previously created Seurat object
Tonsil_1[["Baiting"]] <- Baiting_assay
Tonsil_1[["Protein"]] <- Protein_assay
Tonsil_1[["Hashing"]] <- Hashing_assay
Tonsil_2[["Baiting"]] <- Baiting_assay_2
Tonsil_2[["Protein"]] <- Protein_assay_2
Tonsil_2[["Hashing"]] <- Hashing_assay_2

# Removal of assay objects
remove(Baiting_assay,Hashing_assay,Protein_assay,Tonsil_dataset,Baiting_assay_2,Hashing_assay_2,Protein_assay_2,Tonsil_dataset_2)


# QC and selecting cells for further analysis - first dataset
Tonsil_1[["percent.mt"]] <- PercentageFeatureSet(Tonsil_1, pattern = "^MT-")
VlnPlot(Tonsil_1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
quantile(Tonsil_1@meta.data$nFeature_RNA,probs = seq(0, 1, 1/20))
quantile(Tonsil_1@meta.data$percent.mt,probs = seq(0, 1, 1/20))
nrow(Tonsil_1@meta.data[Tonsil_1@meta.data$percent.mt < 7.5 &
                          Tonsil_1@meta.data$nFeature_RNA > 200 &
                          Tonsil_1@meta.data$nFeature_RNA < 4000,])/nrow(Tonsil_1@meta.data)

Tonsil_1 <- subset(Tonsil_1, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 7.5)

Tonsil_1 <- NormalizeData(Tonsil_1)

# QC and selecting cells for further analysis - second dataset
Tonsil_2[["percent.mt"]] <- PercentageFeatureSet(Tonsil_2, pattern = "^MT-")
VlnPlot(Tonsil_2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
quantile(Tonsil_2@meta.data$nFeature_RNA,probs = seq(0, 1, 1/20))
quantile(Tonsil_2@meta.data$percent.mt,probs = seq(0, 1, 1/20))
nrow(Tonsil_2@meta.data[Tonsil_2@meta.data$percent.mt < 7.5 &
                          Tonsil_2@meta.data$nFeature_RNA > 200 &
                          Tonsil_2@meta.data$nFeature_RNA < 4000,])/nrow(Tonsil_2@meta.data)

Tonsil_2 <- subset(Tonsil_2, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 7.5)

Tonsil_2 <- NormalizeData(Tonsil_2)

# Demultiplexing the HTO data
Tonsil_1 <- NormalizeData(Tonsil_1, assay = "Hashing", normalization.method = "CLR", margin = 2)
Tonsil_2 <- NormalizeData(Tonsil_2, assay = "Hashing", normalization.method = "CLR", margin = 2)
Tonsil_1 <- HTODemux(Tonsil_1, assay = "Hashing", positive.quantile = 0.99)
Tonsil_2 <- HTODemux(Tonsil_2, assay = "Hashing", positive.quantile = 0.99)
table(Tonsil_1@meta.data$Hashing_classification.global)
table(Tonsil_2@meta.data$Hashing_classification.global)
Tonsil_1 <- subset(Tonsil_1, subset = Hashing_classification.global=="Singlet")
Tonsil_2 <- subset(Tonsil_2, subset = Hashing_classification.global=="Singlet")
RidgePlot(Tonsil_1, assay = "Hashing", features = rownames(Tonsil_1[["Hashing"]])[1:10], ncol = 2)
RidgePlot(Tonsil_2, assay = "Hashing", features = rownames(Tonsil_2[["Hashing"]])[1:10], ncol = 2)
table(Tonsil_1@meta.data$Hashing_classification, useNA="always")
table(Tonsil_2@meta.data$Hashing_classification, useNA="always")

# Introducing metadata columns for Patient, Tissue, Vaccinated, Reconvalescent
Tonsil_1@meta.data$Patient <- "None"
Tonsil_2@meta.data$Patient <- "None"
Tonsil_1@meta.data$Patient[grep("243651", Tonsil_1@meta.data$Hashing_classification)] <- "243651-2001"
Tonsil_2@meta.data$Patient[grep("243651", Tonsil_2@meta.data$Hashing_classification)] <- "243651-2001"
Tonsil_1@meta.data$Patient[grep("253340", Tonsil_1@meta.data$Hashing_classification)] <- "253340-1997"
Tonsil_2@meta.data$Patient[grep("253340", Tonsil_2@meta.data$Hashing_classification)] <- "253340-1997"
Tonsil_1@meta.data$Patient[grep("254869", Tonsil_1@meta.data$Hashing_classification)] <- "254869-1996"
Tonsil_2@meta.data$Patient[grep("254869", Tonsil_2@meta.data$Hashing_classification)] <- "254869-1996"
Tonsil_1@meta.data$Patient[grep("277926", Tonsil_1@meta.data$Hashing_classification)] <- "277926-1975"
Tonsil_2@meta.data$Patient[grep("277926", Tonsil_2@meta.data$Hashing_classification)] <- "277926-1975"

Tonsil_1@meta.data$Tissue <- "None"
Tonsil_2@meta.data$Tissue <- "None"
Tonsil_1@meta.data$Tissue[grep("Tonsil",Tonsil_1@meta.data$Hashing_classification)] <- "Tonsil"
Tonsil_2@meta.data$Tissue[grep("Tonsil",Tonsil_2@meta.data$Hashing_classification)] <- "Tonsil"
Tonsil_1@meta.data$Tissue[grep("PBMC",Tonsil_1@meta.data$Hashing_classification)] <- "PBMC"
Tonsil_2@meta.data$Tissue[grep("PBMC",Tonsil_2@meta.data$Hashing_classification)] <- "PBMC"

Tonsil_1@meta.data$Vaccinated <- "no"
Tonsil_2@meta.data$Vaccinated <- "no"
Tonsil_1@meta.data$Vaccinated[grep("Vac",Tonsil_1@meta.data$Hashing_classification)] <- "yes"
Tonsil_2@meta.data$Vaccinated[grep("Vac",Tonsil_2@meta.data$Hashing_classification)] <- "yes"

Tonsil_1@meta.data$Reconvalescent <- "no"
Tonsil_2@meta.data$Reconvalescent <- "no"
Tonsil_1@meta.data$Reconvalescent[grep("Rec",Tonsil_1@meta.data$Hashing_classification)] <- "yes"
Tonsil_2@meta.data$Reconvalescent[grep("Rec",Tonsil_2@meta.data$Hashing_classification)] <- "yes"

######################################################################################################################
# Part 6.3: Subsetting of data and dataset integration
######################################################################################################################
# Renaming of cells
Tonsil_1 <- RenameCells(object = Tonsil_1, add.cell.id = "Dataset_1")
Tonsil_2 <- RenameCells(object = Tonsil_2, add.cell.id = "Dataset_2")
Tonsil_1@meta.data$Dataset <- "First"
Tonsil_2@meta.data$Dataset <- "Second"
Tonsil_1@meta.data$Row.names <- rownames(Tonsil_1@meta.data)
Tonsil_2@meta.data$Row.names <- rownames(Tonsil_2@meta.data)

# Removing of cells selected in Part 6.1
Tonsil_1 <- subset(Tonsil_1, subset = Row.names %notin% Kick)
Tonsil_2 <- subset(Tonsil_2, subset = Row.names %notin% Kick)

remove(Kick)

# During the IntegrateData() function, dimensional reductions of the data are needed. PCA is calculated before Integration is done (only for the protein data)
DefaultAssay(Tonsil_1) <- "Protein"
DefaultAssay(Tonsil_2) <- "Protein"
VariableFeatures(Tonsil_1) <- rownames(Tonsil_1[["Protein"]])
VariableFeatures(Tonsil_2) <- rownames(Tonsil_2[["Protein"]])
Tonsil_1 <- NormalizeData(Tonsil_1, normalization.method = 'CLR', margin = 2) %>% ScaleData() %>% RunPCA(reduction.name = 'apca',approx=FALSE) # approx=F is used to remove the warning that comes otherwise
Tonsil_2 <- NormalizeData(Tonsil_2, normalization.method = 'CLR', margin = 2) %>% ScaleData() %>% RunPCA(reduction.name = 'apca',approx=FALSE)

# Create a list containing the two datasets
Tonsil_dataset_list <- list(Tonsil_1,Tonsil_2)

# Select features that are repeatedly variable across datasets for integration
# The integration is done independently for the RNA and the Protein assay.
features_RNA <- SelectIntegrationFeatures(object.list = Tonsil_dataset_list,assay = c("RNA","RNA")) # IG variable regions pop up!
features_Protein <- SelectIntegrationFeatures(object.list = Tonsil_dataset_list,assay = c("Protein","Protein"))

# Remove HLA, Immunoglobulin, RNA, MT, and RP genes based on HUGO gene names
var_regex = '^HLA-|^IG[HJKL]|^RNA|^MT|^RP'
features_RNA <- grep(var_regex, features_RNA, invert=T, value=T)
remove(var_regex)

# Perform integration. This command requires long time / a lot of computational resources to run.
# Because of only 6 protein features, dims is set to 5 for the anchor finding of protein data
anchors_RNA <- FindIntegrationAnchors(object.list = Tonsil_dataset_list, anchor.features = features_RNA,assay = c("RNA","RNA"))
anchors_Protein <- FindIntegrationAnchors(object.list = Tonsil_dataset_list, anchor.features = features_Protein,assay = c("Protein","Protein"))

# This command creates an 'integrated' data assay in a new Seurat object.
Tonsil.combined.RNA <- IntegrateData(anchorset = anchors_RNA,new.assay.name = "integratedRNA")
Tonsil.combined.RNA <- FindVariableFeatures(Tonsil.combined.RNA, assay = "Protein")
Tonsil.combined.RNA <- ScaleData(Tonsil.combined.RNA, assay = "Protein")
Tonsil.combined.RNA <- RunPCA(Tonsil.combined.RNA, assay = "Protein", reduction.name="pca_protein_all_cells",approx=FALSE)
Tonsil.combined.Protein <- IntegrateEmbeddings(anchorset = anchors_Protein, reductions = Tonsil.combined.RNA@reductions$pca_protein_all_cells)
Tonsil <- Tonsil.combined.RNA
remove(anchors_Protein,anchors_RNA,Tonsil_dataset_list,Tonsil_1,Tonsil_2,Tonsil.combined.RNA)



######################################################################################################################
# Part 6.4: WNN
######################################################################################################################
# Now WNN
# For this, I need to transfer the integrated protein dim reduction from the Bmem.combined.Protein object to the Bmem.combined.RNA object.
# But first, I am preforming data scaling and PCA on the combined RNA project.
DefaultAssay(Tonsil)
DefaultAssay(Tonsil.combined.Protein)

Tonsil <- ScaleData(Tonsil, verbose = FALSE)
Tonsil <- RunPCA(Tonsil, npcs = 30, verbose = FALSE)

# Now we transfer the integrated protein dim reduction
Tonsil@reductions$IntegrateEmbeddings.pca <- Tonsil.combined.Protein@reductions$integrated_dr

# I double check that the correct assays are associated with the correct PCAs.
Tonsil[['pca']]@assay.used
Tonsil[['IntegrateEmbeddings.pca']]@assay.used

# Finally, the FindMultiModalNeighbors function is called, using the PCAs from above.
ElbowPlot(Tonsil, reduction="pca")
Tonsil <- FindMultiModalNeighbors(
  Tonsil, reduction.list = list("pca", "IntegrateEmbeddings.pca"), 
  dims.list = list(1:15, 1:12), modality.weight.name = c("RNA.weight","Protein.weight"))

Tonsil@meta.data$Full.Row.names <- rownames(Tonsil@meta.data)
DefaultAssay(Tonsil) <- "RNA"
Tonsil <- ScaleData(Tonsil)
remove(features_Protein,features_RNA,Tonsil.combined.Protein)

# Finally, I am transfering some columns from the Tonsil.before.subset object to the Tonsil object
Tonsil@meta.data$Row.names <- NULL
Tonsil@meta.data <- merge(Tonsil@meta.data,Tonsil.before.subset@meta.data[,c("wtSpike_specific","RBD_specific","B.1.529_specific",
                                                                             "HE_specific","wtSpike","RBD")],by=0,all.x=T)
rownames(Tonsil@meta.data) <- Tonsil@meta.data$Row.names

######################################################################################################################
# Part 7: Creating Immcantation Input files (Per Patient)
######################################################################################################################
# Prepare the .txt file per patient
Patients <- unique(Tonsil@meta.data$Patient)
Tonsil@meta.data$Row.names <- rownames(Tonsil@meta.data)
positives_243651_2001 <- Tonsil@meta.data$Row.names[Tonsil@meta.data$Patient == "243651-2001" & Tonsil@meta.data$wtSpike_specific=="yes"] 
positives_253340_1997 <- Tonsil@meta.data$Row.names[Tonsil@meta.data$Patient == "253340-1997" & Tonsil@meta.data$wtSpike_specific=="yes"]
positives_277926_1975 <- Tonsil@meta.data$Row.names[Tonsil@meta.data$Patient == "277926-1975" & Tonsil@meta.data$wtSpike_specific=="yes"]
positives_254869_1996 <- Tonsil@meta.data$Row.names[Tonsil@meta.data$Patient == "254869-1996" & Tonsil@meta.data$wtSpike_specific=="yes"]

all_243651_2001 <- Tonsil@meta.data$Row.names[Tonsil@meta.data$Patient == "243651-2001"] 
all_253340_1997 <- Tonsil@meta.data$Row.names[Tonsil@meta.data$Patient == "253340-1997"]
all_277926_1975 <- Tonsil@meta.data$Row.names[Tonsil@meta.data$Patient == "277926-1975"]
all_254869_1996 <- Tonsil@meta.data$Row.names[Tonsil@meta.data$Patient == "254869-1996"]

# Changing the barcodes
# AA -> dataset1
# GG -> dataset2
Immcantation.Input.Barcode.Function <- function(Barcode.Vector){
  i <- 1
  for (i in 1:length(Barcode.Vector)) {
    if(substr(Barcode.Vector[i],9,9)=="1"){
      Barcode.Vector[i] <- paste0("AA",substr(Barcode.Vector[i],11,28))
    }
    if(substr(Barcode.Vector[i],9,9)=="2"){
      Barcode.Vector[i] <- paste0("GG",substr(Barcode.Vector[i],11,28))
    }
    i <- i+1
  } 
  return(Barcode.Vector)
}
positives_243651_2001 <- Immcantation.Input.Barcode.Function(positives_243651_2001)
positives_253340_1997 <- Immcantation.Input.Barcode.Function(positives_253340_1997)
positives_277926_1975 <- Immcantation.Input.Barcode.Function(positives_277926_1975)
positives_254869_1996 <- Immcantation.Input.Barcode.Function(positives_254869_1996)

all_243651_2001 <- Immcantation.Input.Barcode.Function(all_243651_2001)
all_253340_1997 <- Immcantation.Input.Barcode.Function(all_253340_1997)
all_277926_1975 <- Immcantation.Input.Barcode.Function(all_277926_1975)
all_254869_1996 <- Immcantation.Input.Barcode.Function(all_254869_1996)

write(positives_243651_2001, file = "~/NAS/Jan/Experiments and Data/210801 - Exp 035 - Profiling of COVID vaccinated patient samples/Memory B cell Fifth Experiment 221115/Immcantation/Input_data/positives_243651_2001.txt")
write(positives_253340_1997, file = "~/NAS/Jan/Experiments and Data/210801 - Exp 035 - Profiling of COVID vaccinated patient samples/Memory B cell Fifth Experiment 221115/Immcantation/Input_data/positives_253340_1997.txt")
write(positives_277926_1975, file = "~/NAS/Jan/Experiments and Data/210801 - Exp 035 - Profiling of COVID vaccinated patient samples/Memory B cell Fifth Experiment 221115/Immcantation/Input_data/positives_277926_1975.txt")
write(positives_254869_1996, file = "~/NAS/Jan/Experiments and Data/210801 - Exp 035 - Profiling of COVID vaccinated patient samples/Memory B cell Fifth Experiment 221115/Immcantation/Input_data/positives_254869_1996.txt")

write(all_243651_2001, file = "~/NAS/Jan/Experiments and Data/210801 - Exp 035 - Profiling of COVID vaccinated patient samples/Memory B cell Fifth Experiment 221115/Immcantation/Input_data/all_243651_2001.txt")
write(all_253340_1997, file = "~/NAS/Jan/Experiments and Data/210801 - Exp 035 - Profiling of COVID vaccinated patient samples/Memory B cell Fifth Experiment 221115/Immcantation/Input_data/all_253340_1997.txt")
write(all_277926_1975, file = "~/NAS/Jan/Experiments and Data/210801 - Exp 035 - Profiling of COVID vaccinated patient samples/Memory B cell Fifth Experiment 221115/Immcantation/Input_data/all_277926_1975.txt")
write(all_254869_1996, file = "~/NAS/Jan/Experiments and Data/210801 - Exp 035 - Profiling of COVID vaccinated patient samples/Memory B cell Fifth Experiment 221115/Immcantation/Input_data/all_254869_1996.txt")

# Prepare the contig.csv file per patient
contigs_Dataset_1 <- read.csv("~/NAS/Jan/Experiments and Data/210801 - Exp 035 - Profiling of COVID vaccinated patient samples/Memory B cell Fifth Experiment 221115/multi_exp035_revision_2_S1/outs/per_sample_outs/multi_exp035_revision_2_S1/vdj_b/filtered_contig_annotations.csv")
contigs_Dataset_2 <- read.csv("~/NAS/Jan/Experiments and Data/210801 - Exp 035 - Profiling of COVID vaccinated patient samples/Memory B cell Fifth Experiment 221115/multi_exp035_revision_2_S2/outs/per_sample_outs/multi_exp035_revision_2_S2/vdj_b/filtered_contig_annotations.csv")

# Changing barcodes
contigs_Dataset_1$barcode <- paste0("AA",contigs_Dataset_1$barcode)
contigs_Dataset_2$barcode <- paste0("GG",contigs_Dataset_2$barcode)

# Changing contig_id
contigs_Dataset_1$contig_id <- paste0("AA",contigs_Dataset_1$contig_id)
contigs_Dataset_2$contig_id <- paste0("GG",contigs_Dataset_2$contig_id)

# Merging
all.contigs <- rbind(contigs_Dataset_1,contigs_Dataset_2)
remove(contigs_Dataset_1,contigs_Dataset_2)

# Creating patient specific contig objects
contigs_positives_243651_2001 <- all.contigs[all.contigs$barcode %in% positives_243651_2001,]
contigs_positives_253340_1997 <- all.contigs[all.contigs$barcode %in% positives_253340_1997,]
contigs_positives_277926_1975 <- all.contigs[all.contigs$barcode %in% positives_277926_1975,]
contigs_positives_254869_1996 <- all.contigs[all.contigs$barcode %in% positives_254869_1996,]

contigs_all_243651_2001 <- all.contigs[all.contigs$barcode %in% all_243651_2001,]
contigs_all_253340_1997 <- all.contigs[all.contigs$barcode %in% all_253340_1997,]
contigs_all_277926_1975 <- all.contigs[all.contigs$barcode %in% all_277926_1975,]
contigs_all_254869_1996 <- all.contigs[all.contigs$barcode %in% all_254869_1996,]

write_csv(contigs_positives_243651_2001, file = "~/NAS/Jan/Experiments and Data/210801 - Exp 035 - Profiling of COVID vaccinated patient samples/Memory B cell Fifth Experiment 221115/Immcantation/Input_data/contigs_positives_243651_2001.csv")
write_csv(contigs_positives_253340_1997, file = "~/NAS/Jan/Experiments and Data/210801 - Exp 035 - Profiling of COVID vaccinated patient samples/Memory B cell Fifth Experiment 221115/Immcantation/Input_data/contigs_positives_253340_1997.csv")
write_csv(contigs_positives_277926_1975, file = "~/NAS/Jan/Experiments and Data/210801 - Exp 035 - Profiling of COVID vaccinated patient samples/Memory B cell Fifth Experiment 221115/Immcantation/Input_data/contigs_positives_277926_1975.csv")
write_csv(contigs_positives_254869_1996, file = "~/NAS/Jan/Experiments and Data/210801 - Exp 035 - Profiling of COVID vaccinated patient samples/Memory B cell Fifth Experiment 221115/Immcantation/Input_data/contigs_positives_254869_1996.csv")

write_csv(contigs_all_243651_2001, file = "~/NAS/Jan/Experiments and Data/210801 - Exp 035 - Profiling of COVID vaccinated patient samples/Memory B cell Fifth Experiment 221115/Immcantation/Input_data/contigs_all_243651_2001.csv")
write_csv(contigs_all_253340_1997, file = "~/NAS/Jan/Experiments and Data/210801 - Exp 035 - Profiling of COVID vaccinated patient samples/Memory B cell Fifth Experiment 221115/Immcantation/Input_data/contigs_all_253340_1997.csv")
write_csv(contigs_all_277926_1975, file = "~/NAS/Jan/Experiments and Data/210801 - Exp 035 - Profiling of COVID vaccinated patient samples/Memory B cell Fifth Experiment 221115/Immcantation/Input_data/contigs_all_277926_1975.csv")
write_csv(contigs_all_254869_1996, file = "~/NAS/Jan/Experiments and Data/210801 - Exp 035 - Profiling of COVID vaccinated patient samples/Memory B cell Fifth Experiment 221115/Immcantation/Input_data/contigs_all_254869_1996.csv")

remove(all.contigs, contigs_all_cells,contigs_positives_243651_2001,contigs_positives_253340_1997,contigs_positives_254869_1996,contigs_positives_277926_1975,
       all.cells,Patients,positives_243651_2001,positives_253340_1997,positives_254869_1996,positives_277926_1975,
       all_243651_2001,all_253340_1997,all_254869_1996,all_277926_1975,contigs_all_254869_1996,contigs_all_277926_1975,
       contigs_all_253340_1997,contigs_all_243651_2001)


######################################################################################################################
# Part 5: Immcantation Output analysis
######################################################################################################################
Tonsil.backup <- Tonsil
Tonsil <- Tonsil.backup

# Loading the output files
Patient_243651_2001.Immcantation_output <- as.data.frame(read_airr("~/NAS/Jan/Experiments and Data/210801 - Exp 035 - Profiling of COVID vaccinated patient samples/Memory B cell Fifth Experiment 221115/Immcantation/Output_data/positives_243651_2001/outputs/positives_243651_2001_0.05/positives_243651_2001_0.05_heavy_germ-pass.tsv"))
Patient_253340_1997.Immcantation_output <- as.data.frame(read_airr("~/NAS/Jan/Experiments and Data/210801 - Exp 035 - Profiling of COVID vaccinated patient samples/Memory B cell Fifth Experiment 221115/Immcantation/Output_data/positives_253340_1997/outputs/positives_253340_1997_0.05/positives_253340_1997_0.05_heavy_germ-pass.tsv"))
Patient_277926_1975.Immcantation_output <- as.data.frame(read_airr("~/NAS/Jan/Experiments and Data/210801 - Exp 035 - Profiling of COVID vaccinated patient samples/Memory B cell Fifth Experiment 221115/Immcantation/Output_data/positives_277926_1975/outputs/positives_277926_1975_0.05/positives_277926_1975_0.05_heavy_germ-pass.tsv"))
Patient_254869_1996.Immcantation_output <- as.data.frame(read_airr("~/NAS/Jan/Experiments and Data/210801 - Exp 035 - Profiling of COVID vaccinated patient samples/Memory B cell Fifth Experiment 221115/Immcantation/Output_data/positives_254869_1996/outputs/positives_254869_1996_0.05/positives_254869_1996_0.05_heavy_germ-pass.tsv"))

Patient_243651_2001.Immcantation_output_all <- as.data.frame(read_airr("~/NAS/Jan/Experiments and Data/210801 - Exp 035 - Profiling of COVID vaccinated patient samples/Memory B cell Fifth Experiment 221115/Immcantation/Output_data/all_243651_2001/outputs/all_243651_2001_0.05/all_243651_2001_0.05_heavy_germ-pass.tsv"))
Patient_253340_1997.Immcantation_output_all <- as.data.frame(read_airr("~/NAS/Jan/Experiments and Data/210801 - Exp 035 - Profiling of COVID vaccinated patient samples/Memory B cell Fifth Experiment 221115/Immcantation/Output_data/all_253340_1997/outputs/all_253340_1997_0.05/all_253340_1997_0.05_heavy_germ-pass.tsv"))
Patient_277926_1975.Immcantation_output_all <- as.data.frame(read_airr("~/NAS/Jan/Experiments and Data/210801 - Exp 035 - Profiling of COVID vaccinated patient samples/Memory B cell Fifth Experiment 221115/Immcantation/Output_data/all_277926_1975/outputs/all_277926_1975_0.05/all_277926_1975_0.05_heavy_germ-pass.tsv"))
Patient_254869_1996.Immcantation_output_all <- as.data.frame(read_airr("~/NAS/Jan/Experiments and Data/210801 - Exp 035 - Profiling of COVID vaccinated patient samples/Memory B cell Fifth Experiment 221115/Immcantation/Output_data/all_254869_1996/outputs/all_254869_1996_0.05/all_254869_1996_0.05_heavy_germ-pass.tsv"))

All.cells <- read.table("~/NAS/Jan/Experiments and Data/210801 - Exp 035 - Profiling of COVID vaccinated patient samples/Memory B cell Fifth Experiment 221115/Immcantation/Output_data/all.cells/outputs/all.cells_0.05/db_d-masked_germ-pass.tsv", sep = "\t", header = T)

# Taking care of duplicate row names in the All.cells file
All.cells <- do.call("rbind", by(All.cells, All.cells$cell_id, function(x) x[which.max(x$umi_count), ]))

Immcantation.Output.Barcode.Function <- function(Barcode.Vector){
  i <- 1
  for (i in 1:length(Barcode.Vector)) {
    if(substr(Barcode.Vector[i],1,2)=="AA"){
      Barcode.Vector[i] <- paste0("Dataset_1_",substr(Barcode.Vector[i],3,28))
    }
    if(substr(Barcode.Vector[i],1,2)=="GG"){
      Barcode.Vector[i] <- paste0("Dataset_2_",substr(Barcode.Vector[i],3,28))
    }
    i <- i+1
  }
  return(Barcode.Vector)
}
rownames(All.cells) <- Immcantation.Output.Barcode.Function(rownames(All.cells))
All.cells$cell_id <- Immcantation.Output.Barcode.Function(All.cells$cell_id)
Patient_243651_2001.Immcantation_output$cell_id <- Immcantation.Output.Barcode.Function(Patient_243651_2001.Immcantation_output$cell_id)
Patient_253340_1997.Immcantation_output$cell_id <- Immcantation.Output.Barcode.Function(Patient_253340_1997.Immcantation_output$cell_id)
Patient_277926_1975.Immcantation_output$cell_id <- Immcantation.Output.Barcode.Function(Patient_277926_1975.Immcantation_output$cell_id)
Patient_254869_1996.Immcantation_output$cell_id <- Immcantation.Output.Barcode.Function(Patient_254869_1996.Immcantation_output$cell_id)

Patient_243651_2001.Immcantation_output_all$cell_id <- Immcantation.Output.Barcode.Function(Patient_243651_2001.Immcantation_output_all$cell_id)
Patient_253340_1997.Immcantation_output_all$cell_id <- Immcantation.Output.Barcode.Function(Patient_253340_1997.Immcantation_output_all$cell_id)
Patient_277926_1975.Immcantation_output_all$cell_id <- Immcantation.Output.Barcode.Function(Patient_277926_1975.Immcantation_output_all$cell_id)
Patient_254869_1996.Immcantation_output_all$cell_id <- Immcantation.Output.Barcode.Function(Patient_254869_1996.Immcantation_output_all$cell_id)

# Changing the rownames to cell_ids (not needed for All.cells file)
rownames(Patient_243651_2001.Immcantation_output) <- Patient_243651_2001.Immcantation_output$cell_id
rownames(Patient_253340_1997.Immcantation_output) <- Patient_253340_1997.Immcantation_output$cell_id
rownames(Patient_277926_1975.Immcantation_output) <- Patient_277926_1975.Immcantation_output$cell_id
rownames(Patient_254869_1996.Immcantation_output) <- Patient_254869_1996.Immcantation_output$cell_id

rownames(Patient_243651_2001.Immcantation_output_all) <- Patient_243651_2001.Immcantation_output_all$cell_id
rownames(Patient_253340_1997.Immcantation_output_all) <- Patient_253340_1997.Immcantation_output_all$cell_id
rownames(Patient_277926_1975.Immcantation_output_all) <- Patient_277926_1975.Immcantation_output_all$cell_id
rownames(Patient_254869_1996.Immcantation_output_all) <- Patient_254869_1996.Immcantation_output_all$cell_id

# Making the clone IDs unique to patients
Patient_243651_2001.Immcantation_output$clone_id <- paste0(Patient_243651_2001.Immcantation_output$clone_id,"_1")
Patient_253340_1997.Immcantation_output$clone_id <- paste0(Patient_253340_1997.Immcantation_output$clone_id,"_2")
Patient_277926_1975.Immcantation_output$clone_id <- paste0(Patient_277926_1975.Immcantation_output$clone_id,"_3")
Patient_254869_1996.Immcantation_output$clone_id <- paste0(Patient_254869_1996.Immcantation_output$clone_id,"_3")

Patient_243651_2001.Immcantation_output_all$clone_id <- paste0(Patient_243651_2001.Immcantation_output_all$clone_id,"_1_all")
Patient_253340_1997.Immcantation_output_all$clone_id <- paste0(Patient_253340_1997.Immcantation_output_all$clone_id,"_2_all")
Patient_277926_1975.Immcantation_output_all$clone_id <- paste0(Patient_277926_1975.Immcantation_output_all$clone_id,"_3_all")
Patient_254869_1996.Immcantation_output_all$clone_id <- paste0(Patient_254869_1996.Immcantation_output_all$clone_id,"_4_all")

# Merging the Immcantation output files and adding them to the Seurat object
Immcantation_output <- rbind(Patient_243651_2001.Immcantation_output,Patient_253340_1997.Immcantation_output,
                             Patient_277926_1975.Immcantation_output,Patient_254869_1996.Immcantation_output)
Immcantation_output_all <- rbind(Patient_243651_2001.Immcantation_output_all,Patient_253340_1997.Immcantation_output_all,
                                 Patient_254869_1996.Immcantation_output_all,Patient_277926_1975.Immcantation_output_all)
Tonsil@meta.data$Row.names <- NULL
Tonsil.all.Immcantation.output <- Tonsil
Tonsil@meta.data <- merge(Tonsil@meta.data,Immcantation_output[,c("cdr3","clone_id")],by=0, all.x=T)
Tonsil.all.Immcantation.output@meta.data <- merge(Tonsil.all.Immcantation.output@meta.data,Immcantation_output_all[,c("cdr3","clone_id")],by=0,all.x=T)
rownames(Tonsil@meta.data) <- Tonsil@meta.data$Row.names
rownames(Tonsil.all.Immcantation.output@meta.data) <- Tonsil.all.Immcantation.output@meta.data$Row.names

# Doing the same for the object from before subsetting:
Tonsil.before.subset@meta.data$Row.names <- NULL
Tonsil.before.subset.all.Immcantation.output <- Tonsil.before.subset
Tonsil.before.subset@meta.data <- merge(Tonsil.before.subset@meta.data,Immcantation_output[,c("cdr3","clone_id")],by=0, all.x=T)
Tonsil.before.subset.all.Immcantation.output@meta.data <- merge(Tonsil.before.subset.all.Immcantation.output@meta.data,Immcantation_output_all[,c("cdr3","clone_id")],by=0,all.x=T)
rownames(Tonsil.before.subset@meta.data) <- Tonsil.before.subset@meta.data$Row.names
rownames(Tonsil.before.subset.all.Immcantation.output@meta.data) <- Tonsil.before.subset.all.Immcantation.output@meta.data$Row.names

remove(Immcantation_output,Patient_243651_2001.Immcantation_output,Patient_253340_1997.Immcantation_output,
       Patient_254869_1996.Immcantation_output,Patient_277926_1975.Immcantation_output)
remove(Immcantation_output_all,Patient_243651_2001.Immcantation_output_all,Patient_253340_1997.Immcantation_output_all,
       Patient_254869_1996.Immcantation_output_all,Patient_277926_1975.Immcantation_output_all)

# This gives the counts of mutations for the All.cells file
obs <- observedMutations(All.cells, sequenceColumn="sequence_alignment",
                         germlineColumn="germline_alignment_d_mask",
                         regionDefinition=NULL,
                         frequency=FALSE,
                         combine=TRUE,
                         nproc=1)
obs <- obs[,c("cell_id","mu_count")]
rownames(obs) <- obs$cell_id

# This gives the frequency of mutations for the All.cells file
freq <- observedMutations(All.cells, sequenceColumn="sequence_alignment",
                          germlineColumn="germline_alignment_d_mask",
                          regionDefinition=NULL,
                          frequency=TRUE, 
                          combine=TRUE,
                          nproc=1)
freq <- freq[,c("cell_id","mu_freq")]
rownames(freq) <- freq$cell_id

obs_freq <- merge(obs,freq,by=0)
rownames(obs_freq) <- obs_freq$Row.names
obs_freq <- obs_freq[,c("mu_count","mu_freq")]

All.cells <- merge(All.cells,obs_freq,by=0)
rownames(All.cells) <- All.cells$Row.names


# Now I am adding Isotypes and heavy chain info to the cells of the Seurat object. Also for the Seurat Object Tonsil.before.subset
Tonsil@meta.data$Row.names <- NULL
Tonsil.all.Immcantation.output@meta.data$Row.names <- NULL
Tonsil.before.subset@meta.data$Row.names <- NULL
Tonsil.before.subset.all.Immcantation.output@meta.data$Row.names <- NULL

Tonsil@meta.data <- merge(Tonsil@meta.data,All.cells[,c("v_call","d_call","j_call","c_call","mu_count","mu_freq")],by=0,all.x=T)
Tonsil.all.Immcantation.output@meta.data <- merge(Tonsil.all.Immcantation.output@meta.data,All.cells[,c("v_call","d_call","j_call","c_call","mu_count","mu_freq")],by=0,all.x=T)
Tonsil.before.subset@meta.data <- merge(Tonsil.before.subset@meta.data,All.cells[,c("v_call","d_call","j_call","c_call","mu_count","mu_freq")],by=0,all.x=T)
Tonsil.before.subset.all.Immcantation.output@meta.data <- merge(Tonsil.before.subset.all.Immcantation.output@meta.data,All.cells[,c("v_call","d_call","j_call","c_call","mu_count","mu_freq")],by=0,all.x=T)

rownames(Tonsil@meta.data) <- Tonsil@meta.data$Row.names
rownames(Tonsil.all.Immcantation.output@meta.data) <- Tonsil.all.Immcantation.output@meta.data$Row.names
rownames(Tonsil.before.subset@meta.data) <- Tonsil.before.subset@meta.data$Row.names
rownames(Tonsil.before.subset.all.Immcantation.output@meta.data) <- Tonsil.before.subset.all.Immcantation.output@meta.data$Row.names

Tonsil@meta.data$c_call[Tonsil@meta.data$c_call!="IGHA1" &
                          Tonsil@meta.data$c_call!="IGHA2"&
                          Tonsil@meta.data$c_call!="IGHD"&
                          Tonsil@meta.data$c_call!="IGHG1"&
                          Tonsil@meta.data$c_call!="IGHG2"&
                          Tonsil@meta.data$c_call!="IGHG3"&
                          Tonsil@meta.data$c_call!="IGHG4"&
                          Tonsil@meta.data$c_call!="IGHM"&
                          !is.na(Tonsil@meta.data$c_call)] <- "Unknown"
Tonsil.all.Immcantation.output@meta.data$c_call[Tonsil.all.Immcantation.output@meta.data$c_call!="IGHA1" &
                                                  Tonsil.all.Immcantation.output@meta.data$c_call!="IGHA2"&
                                                  Tonsil.all.Immcantation.output@meta.data$c_call!="IGHD"&
                                                  Tonsil.all.Immcantation.output@meta.data$c_call!="IGHG1"&
                                                  Tonsil.all.Immcantation.output@meta.data$c_call!="IGHG2"&
                                                  Tonsil.all.Immcantation.output@meta.data$c_call!="IGHG3"&
                                                  Tonsil.all.Immcantation.output@meta.data$c_call!="IGHG4"&
                                                  Tonsil.all.Immcantation.output@meta.data$c_call!="IGHM"&
                                                  !is.na(Tonsil.all.Immcantation.output@meta.data$c_call)] <- "Unknown"
Tonsil.before.subset@meta.data$c_call[Tonsil.before.subset@meta.data$c_call!="IGHA1" &
                                        Tonsil.before.subset@meta.data$c_call!="IGHA2"&
                                        Tonsil.before.subset@meta.data$c_call!="IGHD"&
                                        Tonsil.before.subset@meta.data$c_call!="IGHG1"&
                                        Tonsil.before.subset@meta.data$c_call!="IGHG2"&
                                        Tonsil.before.subset@meta.data$c_call!="IGHG3"&
                                        Tonsil.before.subset@meta.data$c_call!="IGHG4"&
                                        Tonsil.before.subset@meta.data$c_call!="IGHM"&
                          !is.na(Tonsil.before.subset@meta.data$c_call)] <- "Unknown"
Tonsil.before.subset.all.Immcantation.output@meta.data$c_call[Tonsil.before.subset.all.Immcantation.output@meta.data$c_call!="IGHA1" &
                                                                Tonsil.before.subset.all.Immcantation.output@meta.data$c_call!="IGHA2"&
                                                                Tonsil.before.subset.all.Immcantation.output@meta.data$c_call!="IGHD"&
                                                                Tonsil.before.subset.all.Immcantation.output@meta.data$c_call!="IGHG1"&
                                                                Tonsil.before.subset.all.Immcantation.output@meta.data$c_call!="IGHG2"&
                                                                Tonsil.before.subset.all.Immcantation.output@meta.data$c_call!="IGHG3"&
                                                                Tonsil.before.subset.all.Immcantation.output@meta.data$c_call!="IGHG4"&
                                                                Tonsil.before.subset.all.Immcantation.output@meta.data$c_call!="IGHM"&
                                                  !is.na(Tonsil.before.subset.all.Immcantation.output@meta.data$c_call)] <- "Unknown"

Tonsil@meta.data$c_call[is.na(Tonsil@meta.data$c_call)] <- "Unknown"
Tonsil.all.Immcantation.output@meta.data$c_call[is.na(Tonsil.all.Immcantation.output@meta.data$c_call)] <- "Unknown"
Tonsil.before.subset@meta.data$c_call[is.na(Tonsil.before.subset@meta.data$c_call)] <- "Unknown"
Tonsil.before.subset.all.Immcantation.output@meta.data$c_call[is.na(Tonsil.before.subset.all.Immcantation.output@meta.data$c_call)] <- "Unknown"

remove(All.cells, obs, freq, obs_freq)


# Testing if a patient has cross-tissue clones - clones that we find in the PBMCs and the Tonsil
# The for loop also adds the "cross.tissue.clone" info to the Seurat object meta data.
Tonsil.subset <- subset(Tonsil, subset= clone_id != is.na(clone_id))
Tonsil.all.Immcantation.output.subset <- subset(Tonsil.all.Immcantation.output,subset=clone_id != is.na(clone_id))
Tonsil.subset.list <- Tonsil.subset@meta.data %>% dplyr::group_by(clone_id) %>% group_split()
Tonsil.all.Immcantation.output.subset.list <- Tonsil.all.Immcantation.output.subset@meta.data %>% dplyr::group_by(clone_id) %>% group_split()

i <- 1
Tonsil@meta.data$cross.tissue.clone <- "no"
Tonsil.all.Immcantation.output@meta.data$cross.tissue.clone <- "no"
for (i in 1:length(Tonsil.subset.list)) {
  
  this.round.clone <- as.data.frame(Tonsil.subset.list[[i]])
  if(length(unique(this.round.clone$Tissue))>1){
    Tonsil@meta.data$cross.tissue.clone[Tonsil@meta.data$Row.names %in% Tonsil.subset.list[[i]]$Row.names] <- "yes"
  } else {next}
  
}
i <- 1
for (i in 1:length(Tonsil.all.Immcantation.output.subset.list)) {
  
  this.round.clone <- as.data.frame(Tonsil.all.Immcantation.output.subset.list[[i]])
  if(length(unique(this.round.clone$Tissue))>1){
    Tonsil.all.Immcantation.output@meta.data$cross.tissue.clone[Tonsil.all.Immcantation.output@meta.data$Row.names %in% Tonsil.all.Immcantation.output.subset.list[[i]]$Row.names] <- "yes"
  } else {next}
  
}
remove(Tonsil.subset,Tonsil.subset.list,this.round.clone,i,Tonsil.all.Immcantation.output.subset,Tonsil.all.Immcantation.output.subset.list)

# Check within the cross tissue clones which clone_id they are and how many cells they have in which tissue
Cross.tissues <- subset(Tonsil, subset=cross.tissue.clone=="yes")
Cross.tissues.all <- subset(Tonsil.all.Immcantation.output,subset=cross.tissue.clone=="yes")
ncol(table(Cross.tissues@meta.data$Tissue,Cross.tissues@meta.data$clone_id))
ncol(table(Cross.tissues.all@meta.data$Tissue,Cross.tissues.all@meta.data$clone_id))

remove(Cross.tissues,Cross.tissues.all)

# Adding a clonal frequency column
Tonsil@meta.data <- Tonsil@meta.data %>% group_by(clone_id)%>% add_tally(name = "clonal.frequency") %>% ungroup() %>% as.data.frame()
Tonsil.all.Immcantation.output@meta.data <- Tonsil.all.Immcantation.output@meta.data %>% group_by(clone_id)%>% add_tally(name = "clonal.frequency") %>% ungroup() %>% as.data.frame()
rownames(Tonsil@meta.data) <- Tonsil@meta.data$Row.names
rownames(Tonsil.all.Immcantation.output@meta.data) <- Tonsil.all.Immcantation.output@meta.data$Row.names
Tonsil@meta.data$clonal.frequency[is.na(Tonsil@meta.data$clone_id)] <- NA
Tonsil.all.Immcantation.output@meta.data$clonal.frequency[is.na(Tonsil.all.Immcantation.output@meta.data$clone_id)] <- NA

table(Tonsil@meta.data$clonal.frequency)
table(Tonsil.all.Immcantation.output@meta.data$clonal.frequency)

remove(Tonsil.backup)


######################################################################################################################
# Part 8: Clustering and UMAP visualization after subsetting
######################################################################################################################
# First clustering, UMAP and Dimplot on RNA data only
Tonsil <- FindNeighbors(Tonsil, dims = 1:15,graph.name = "integratedRNA_nn",assay="integratedRNA",reduction = "pca")
Tonsil <- FindClusters(Tonsil, graph.name = "integratedRNA_nn", algorithm = 1, verbose = FALSE,resolution = 0.5)
Tonsil <- RunUMAP(Tonsil, reduction = 'pca', dims = 1:15, assay = 'integratedRNA', reduction.name = 'rna.umap', reduction.key = 'rnaUMAP_')

DimPlot(Tonsil, reduction = 'rna.umap', label = F, repel = TRUE, label.size = 2.5, pt.size = 0.2,group.by = "integratedRNA_nn_res.0.5", 
        cols = MetBrewer::met.brewer(name="Egypt", n=length(unique(Tonsil@meta.data$seurat_clusters))))+ggtitle("RNA UMAP")+center.title() + NoAxes() # Export 4x6


# Then clustering, UMAP and Dimplot on Protein data only
Tonsil <- FindNeighbors(Tonsil, dims = 1:12,graph.name = "integratedADT_nn",assay="Protein",reduction = "IntegrateEmbeddings.pca")
Tonsil <- FindClusters(Tonsil, graph.name = "integratedADT_nn", algorithm = 1, verbose = FALSE)
Tonsil <- RunUMAP(Tonsil, reduction = 'IntegrateEmbeddings.pca', dims = 1:12, assay = 'Protein', reduction.name = 'adt.umap', reduction.key = 'adtUMAP_')

DimPlot(Tonsil, reduction = 'adt.umap', label = F, repel = TRUE, label.size = 2.5, pt.size = 0.6, 
        cols = MetBrewer::met.brewer(name="Egypt", n=length(unique(Tonsil@meta.data$seurat_clusters))))+ggtitle("Protein UMAP")+center.title()


# Then on the WNN data - find multimodal neighbours done already in Part 3.
DefaultAssay(Tonsil)
Tonsil <- FindClusters(Tonsil, graph.name = "wsnn", algorithm = 1, verbose = FALSE,resolution = 0.4)
Tonsil <- RunUMAP(Tonsil, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")

Tonsil@meta.data$seurat_clusters <- as.numeric(as.character(Tonsil@meta.data$seurat_clusters))
Tonsil@meta.data$seurat_clusters<- Tonsil@meta.data$seurat_clusters+1
Tonsil@meta.data$seurat_clusters <- factor(Tonsil@meta.data$seurat_clusters)
Idents(Tonsil) <- "seurat_clusters"

DimPlot(Tonsil, reduction = 'wnn.umap', label = F, repel = TRUE, label.size = 2.5, pt.size = 0.5, 
        cols = MetBrewer::met.brewer(name="Egypt", n=length(unique(Tonsil@meta.data$seurat_clusters))))+
  theme(text = element_text(size = 20)) +ggtitle(paste("WNN UMAP"))+center.title()+NoAxes() #4x6

DimPlot(Tonsil, reduction = 'wnn.umap', label = F, repel = TRUE, label.size = 2.5, group.by = "Tissue", pt.size = 0.5, 
        cols = MetBrewer::met.brewer(name="Egypt", n=length(unique(Tonsil@meta.data$Tissue))))+
  theme(text = element_text(size = 20)) +ggtitle(paste("WNN UMAP"))+center.title()+NoAxes() #4x6.5

Make.stacked.percentage.plot(Seurat.Object = Tonsil, groups.on.Y.axis = "seurat_clusters",groups.in.legend="c_call",
                             X.Axis.title = "Isotypes",Y.Axis.title = "Clusters",Legend.title = "Isotypes",Plot.title = "Isotypes across Clusters") # 5x7

Tonsil.subset <- subset(Tonsil, subset=c_call != "Unknown")
Make.stacked.percentage.plot(Seurat.Object = Tonsil.subset, groups.on.Y.axis = "seurat_clusters",groups.in.legend="c_call",
                             X.Axis.title = "Isotypes",Y.Axis.title = "Clusters",Legend.title = "Isotypes",Plot.title = "Isotypes across Clusters") # 5x7

Make.stacked.percentage.plot(Seurat.Object = Tonsil, groups.on.Y.axis = "seurat_clusters",groups.in.legend="Patient",
                             X.Axis.title = "Patients",Y.Axis.title = "Clusters",Legend.title = "Patients",Plot.title = "Patients across Clusters") # 5x7
Tonsil.subset <- subset(Tonsil, subset=wtSpike_specific=="yes")
Make.stacked.percentage.plot(Seurat.Object = Tonsil.subset, groups.on.Y.axis = "seurat_clusters",groups.in.legend="Patient",
                             X.Axis.title = "Patients",Y.Axis.title = "Clusters",Legend.title = "Patients",Plot.title = "Patients across Clusters \n Only wtSpike binding cells") # 5x7
remove(Tonsil.subset)

p <- DimPlot(Tonsil, reduction = 'wnn.umap', label = T, repel = TRUE, label.size = 2.5, pt.size = 0.5, 
             cols = MetBrewer::met.brewer(name="Egypt", n=length(unique(Tonsil@meta.data$seurat_clusters))))+
  theme(text = element_text(size = 20)) +ggtitle(paste("WNN UMAP"))+center.title()+
  DimPlot(Tonsil,reduction = "wnn.umap",group.by = "wtSpike_specific",
          cols = MetBrewer::met.brewer(name="Austria", n=2))+ 
  DimPlot(Tonsil,reduction = "wnn.umap",group.by = "Dataset", 
          cols = MetBrewer::met.brewer(name="Austria", n=2))+ 
  DimPlot(Tonsil, reduction = "wnn.umap", group.by = "Vaccinated", 
          cols = MetBrewer::met.brewer(name="Austria", n=2))+ 
  DimPlot(Tonsil, reduction = "wnn.umap", group.by = "Patient", 
          cols = MetBrewer::met.brewer(name="Austria", n=4))+ 
  DimPlot(Tonsil, reduction = "wnn.umap", group.by = "Reconvalescent", 
          cols = MetBrewer::met.brewer(name="Austria", n=2))+ 
  DimPlot(Tonsil, reduction = "wnn.umap", group.by = "Tissue", 
          cols = MetBrewer::met.brewer(name="Austria", n=2))+ 
  DimPlot(Tonsil, reduction = "wnn.umap", group.by = "c_call", 
          cols = MetBrewer::met.brewer(name="Austria", n=9))


p & NoAxes()
remove(p)

DefaultAssay(Tonsil) <- "Protein"
FeaturePlot(Tonsil, features = rownames(Tonsil@assays$Protein),reduction = "wnn.umap") & NoAxes() #7x12

Genelist <- c("MS4A1","CR2","SELL","CD69","CXCR4","CCR7","KLF2","TNFRSF13B","ENTPD1","FCRL4","IGHD","IGHM","IGHA1","IGHA2","IGHG1")

DotPlot(Tonsil, assay = "RNA",features = Genelist)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ylab("Cluster")+xlab("Gene Markers")+
  DotPlot(Tonsil, assay = "Protein",features = rownames(Tonsil@assays$Protein@counts),cols = c("lightgrey","darkgreen"))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ylab("Cluster")+xlab("TotalSeq Markers") #5x10.5

remove(Genelist)

# Highlighting Spike binders on the UMAP
highlights <- rownames(Tonsil@meta.data[Tonsil@meta.data$wtSpike_specific=="yes",])
DimPlot(Tonsil, reduction = 'wnn.umap', label = F, repel = TRUE, sizes.highlight = 0.5, pt.size = 0.5,cells.highlight = highlights)+
  theme(text = element_text(size = 20)) +ggtitle(paste("wtSpike binding cells\n highlighted on UMAP"))+center.title()+NoAxes() #4x7
remove(highlights)

#Marker <- FindAllMarkers(Tonsil, assay = "RNA")
#Marker <- Marker[Marker$p_val_adj< 0.0005&
#                  Marker$avg_log2FC >2,]



######################################################################################################################
# Part 9: Dataset Exploration
######################################################################################################################
# Comparing SHMs accross tissues - spike specific cells including naive cells
Naives <- subset(Tonsil.before.subset,subset=Subset=="Naive")
df2 <- Naives@meta.data[,c("Subset","mu_count")]
colnames(df2) <- c("Tissue","mu_count")
Tonsil.subset <- subset(Tonsil, subset=wtSpike_specific=="yes")
df <- Tonsil.subset@meta.data[,c("Tissue","mu_count")]
df$Tissue[df$Tissue=="PBMC"] <- "Circulating"
df <- rbind(df,df2)
df$Tissue <- factor(df$Tissue, levels = c("Naive","Circulating","Tonsil"))
stat.test <- df %>% rstatix::wilcox_test(mu_count~Tissue,detailed = T)

ggplot(df, aes(x=Tissue,y=mu_count))+ 
  geom_boxplot(df, mapping= aes(x = Tissue, y = mu_count,fill=Tissue))+
  theme_classic()+stat_pvalue_manual(stat.test, label = "p.adj", 
                                     y.position = c(42,48,54))+
  labs(x="",y="Mutational count",fill = "Subset")+ black.axis.text()+theme(text = element_text(size = 14))+          
  scale_fill_manual(values=met.brewer(name="Egypt", n=length(unique(df$Tissue))))+ggtitle("Spike Wuhan-Hu-1 binding cells")+
  center.title()# Export 5x7

remove(df,df2,Naives,Tonsil.subset,stat.test)


# Looking at the connection of longitudinal clones on the UMAP
clone_IDs <- unique(Tonsil@meta.data$clone_id[Tonsil@meta.data$cross.tissue.clone=="yes"])
DimPlot(Tonsil,reduction = "wnn.umap",group.by = "wtSpike_specific",
        cols = MetBrewer::met.brewer(name="Austria", n=2))+labs(color = "wtSpike specific")+
  Show.clones.on.UMAP.one.color(Seurat.Object = Tonsil,clone_IDs = clone_IDs, reduction = "wnn.umap",
                                DimplotGroup = "Tissue",connection.color = "gray34")+NoAxes()+
  DimPlot(Tonsil,reduction = "wnn.umap",group.by = "seurat_clusters",
          cols = MetBrewer::met.brewer(name="Egypt", n=8))+labs(color="Tissue origin") # Export 4x14 or 3x4 clonal connections only 4.5x7


DimPlot(Tonsil,reduction = "wnn.umap",group.by = "wtSpike_specific",
        cols = MetBrewer::met.brewer(name="Austria", n=2))+labs(color = "wtSpike specific")+
  Show.clones.on.UMAP(Seurat.Object = Tonsil,clone_IDs = clone_IDs, reduction = "wnn.umap",
                      DimplotGroup = "Tissue")+
  DimPlot(Tonsil,reduction = "wnn.umap",group.by = "Tissue",
          cols = MetBrewer::met.brewer(name="Egypt", n=2))+labs(color="Tissue origin") # Export 4x14

remove(clone_IDs)

# Plotting the clonal connection per patient
clone_IDs_243651_2001 <- unique(Tonsil@meta.data$clone_id[Tonsil@meta.data$cross.tissue.clone=="yes" & Tonsil@meta.data$Patient=="243651-2001"])
clone_IDs_253340_1997 <- unique(Tonsil@meta.data$clone_id[Tonsil@meta.data$cross.tissue.clone=="yes" & Tonsil@meta.data$Patient=="253340-1997"])
clone_IDs_254869_1996 <- unique(Tonsil@meta.data$clone_id[Tonsil@meta.data$cross.tissue.clone=="yes" & Tonsil@meta.data$Patient=="254869-1996"])
clone_IDs_277926_1975 <- unique(Tonsil@meta.data$clone_id[Tonsil@meta.data$cross.tissue.clone=="yes" & Tonsil@meta.data$Patient=="277926-1975"])

Show.clones.on.UMAP.one.color(Seurat.Object = Tonsil,clone_IDs = clone_IDs_253340_1997, reduction = "wnn.umap",
                              DimplotGroup = "Tissue",connection.color = "gray34")
Show.clones.on.UMAP.one.color(Seurat.Object = Tonsil,clone_IDs = clone_IDs_254869_1996, reduction = "wnn.umap",
                              DimplotGroup = "Tissue",connection.color = "gray34")
Show.clones.on.UMAP.one.color(Seurat.Object = Tonsil,clone_IDs = clone_IDs_277926_1975, reduction = "wnn.umap",
                              DimplotGroup = "Tissue",connection.color = "gray34")

# Checking the hashing reads for individual cross tissue clones to exclude the possibility of doublets / wrong sample calls
Tonsil@assays$Hashing@counts[,rownames(Tonsil@meta.data[Tonsil@meta.data$clone_id=="201_8_2"&
                                                          !is.na(Tonsil@meta.data$clone_id),])]


# Making a venn diagram showing the clonal overlap of tissues - over all patients
Patient_243651_2001 <- subset(Tonsil, subset=Patient=="243651-2001")
Patient_253340_1997 <- subset(Tonsil, subset=Patient=="253340-1997")
Patient_254869_1996 <- subset(Tonsil, subset=Patient=="254869-1996")
Patient_277926_1975 <- subset(Tonsil, subset=Patient=="277926-1975")

Patient_243651_2001 <- Patient_243651_2001@meta.data[!is.na(Patient_243651_2001@meta.data$clone_id),c("clone_id","Tissue","cross.tissue.clone")]
Patient_253340_1997 <- Patient_253340_1997@meta.data[!is.na(Patient_253340_1997@meta.data$clone_id),c("clone_id","Tissue","cross.tissue.clone")]
Patient_254869_1996 <- Patient_254869_1996@meta.data[!is.na(Patient_254869_1996@meta.data$clone_id),c("clone_id","Tissue","cross.tissue.clone")]
Patient_277926_1975 <- Patient_277926_1975@meta.data[!is.na(Patient_277926_1975@meta.data$clone_id),c("clone_id","Tissue","cross.tissue.clone")]

Patientlist <- list("Patient_243651_2001"=Patient_243651_2001,"Patient_253340_1997"=Patient_253340_1997,
                    "Patient_254869_1996"=Patient_254869_1996,"Patient_277926_1975"=Patient_277926_1975)

i <- 1
df <- data.frame("Patient"=c("Test","Test","Test","Test"),"n.clones.tonsil"=c(0,0,0,0),"n.clones.PBMC"=c(0,0,0,0),"n.clones.cross.tissue"=c(0,0,0,0))
for (i in 1:length(Patientlist)) {
  df[i,1] <- names(Patientlist[i])
  df[i,2] <- length(unique(Patientlist[[i]][Patientlist[[i]]$Tissue=="Tonsil","clone_id"]))
  df[i,3] <- length(unique(Patientlist[[i]][Patientlist[[i]]$Tissue=="PBMC","clone_id"]))
  df[i,4] <- length(unique(Patientlist[[i]][Patientlist[[i]]$cross.tissue.clone=="yes","clone_id"]))
}

# Creating the Venn diagramm using eulerr
# See documentation https://cran.r-project.org/web/packages/eulerr/vignettes/gallery.html
library(eulerr)
vd <- euler(c("Tonsil" = sum(df$n.clones.tonsil), "PBMC" = sum(df$n.clones.PBMC), "Tonsil&PBMC" = sum(df$n.clones.cross.tissue)))
plot(vd, counts = T,lwd = 2,
     fill=c("#004488","#997700"),
     opacity = .7,
     legend = list( space= "right", columns=1),
     quantities = list(type="counts"),
     main="Clonal overlap between Tissues") # Export 4x6


# What percentage of the expanded clones belongs to cross tissue clones?
Percentage.cross <- 100*(length(unique(Tonsil@meta.data$clone_id[Tonsil@meta.data$cross.tissue.clone=="yes"&
                                               Tonsil@meta.data$clonal.frequency>1&
                                               !is.na(Tonsil@meta.data$clonal.frequency)])))/
  (length(unique(Tonsil@meta.data$clone_id[Tonsil@meta.data$clonal.frequency>1&
                                             !is.na(Tonsil@meta.data$clonal.frequency)])))

# Calculating this for the object that contains the Immcantation output using all cells from each patient.
Percentage.cross2 <- 100*(length(unique(Tonsil.all.Immcantation.output@meta.data$clone_id[Tonsil.all.Immcantation.output@meta.data$cross.tissue.clone=="yes"&
                                                                                           Tonsil.all.Immcantation.output@meta.data$clonal.frequency>1&
                                                                   !is.na(Tonsil.all.Immcantation.output@meta.data$clonal.frequency)&
                                                                     Tonsil.all.Immcantation.output@meta.data$wtSpike_specific=="yes"])))/
  (length(unique(Tonsil.all.Immcantation.output@meta.data$clone_id[Tonsil.all.Immcantation.output@meta.data$clonal.frequency>1&
                                             !is.na(Tonsil.all.Immcantation.output@meta.data$clonal.frequency)&
                                               Tonsil.all.Immcantation.output@meta.data$wtSpike_specific=="yes"])))

Percentage.cross2 <- round(Percentage.cross2,2)
Percentage.not.cross <- 100-Percentage.cross2

data <- data.frame(
  Cross=c("Cross tissue clone","Tissue unique clone"),
  value=c(Percentage.cross2,Percentage.not.cross))
data$Cross <- factor(data$Cross, levels = c("Tissue unique clone","Cross tissue clone"))
colnames(data) <- c("Clone.distribution","value")

ggplot(data, aes(x="", y=value, fill=Clone.distribution)) +
  geom_bar(stat="identity", width=1,colour="black",size=1) +
  coord_polar("y", start=0)+theme_void()+
  scale_fill_manual(values=met.brewer("Lakota", 2))+ggtitle("Expanded and Wuhan-Hu-1 binding clones\n found in Tonsil and PBMC")+
  geom_text(x=1.1, y=15, label=paste0(Percentage.cross2,"%"), colour="black",size=6)+center.title()+labs(fill="Clone distribution") # Export 5.5x5.5

remove(Percentage.cross,Percentage.not.cross,data,Percentage.cross2)


# Before making the more complex Venn diagram that includes non spike specific cells, I want to know:
# How many clones contain cells with mixed "wtSpike_specific" call?
# Number of total clones:
length(unique(Tonsil.all.Immcantation.output@meta.data$clone_id))

# Number of clones with at least one wtSpike_specific cell:
length(unique(Tonsil.all.Immcantation.output@meta.data$clone_id[Tonsil.all.Immcantation.output@meta.data$wtSpike_specific=="yes"]))

# Number of clones with mixed "wtSpuke_specific" call:
Tonsil.all.Immcantation.output@meta.data$mixed.call <- "not.mixed"
Tonsil.all.Immcantation.output@meta.data$mixed.call[is.na(Tonsil.all.Immcantation.output@meta.data$clone_id)] <- NA
All.clones <- unique(Tonsil.all.Immcantation.output@meta.data$clone_id)
i <- 1
for (i in 1:length(All.clones)) {
  if(length(unique(Tonsil.all.Immcantation.output@meta.data$wtSpike_specific[Tonsil.all.Immcantation.output@meta.data$clone_id %in% All.clones[i]]))>1){
    Tonsil.all.Immcantation.output@meta.data$mixed.call[Tonsil.all.Immcantation.output@meta.data$clone_id==All.clones[i]] <- "mixed"
  } else{
    next
  }
}
length(unique(Tonsil.all.Immcantation.output@meta.data$clone_id[Tonsil.all.Immcantation.output@meta.data$mixed.call=="mixed"&
                                                                  !is.na(Tonsil.all.Immcantation.output@meta.data$mixed.call)]))

# Make 4 vectors, one for each of the 4 categories that will be represented in the venn diagram
# 1: Spike specific clones that are found in the Tonsil
# 2: Spike specific clones that are found in the Circulation
# 3: Non-spike specific clones that are found in the Tonsil
# 4: Non-spike specific clones that are found in the Circulation
Tonsil.Binding <- unique(Tonsil.all.Immcantation.output@meta.data$clone_id[Tonsil.all.Immcantation.output@meta.data$Tissue=="Tonsil"&
                                                                           Tonsil.all.Immcantation.output@meta.data$wtSpike_specific=="yes"&
                                                                           !is.na(Tonsil.all.Immcantation.output@meta.data$clone_id)])
Circulation.Binding <- unique(Tonsil.all.Immcantation.output@meta.data$clone_id[Tonsil.all.Immcantation.output@meta.data$Tissue=="PBMC"&
                                                                            Tonsil.all.Immcantation.output@meta.data$wtSpike_specific=="yes"&
                                                                            !is.na(Tonsil.all.Immcantation.output@meta.data$clone_id)])
Tonsil.Non.Binding <- unique(Tonsil.all.Immcantation.output@meta.data$clone_id[Tonsil.all.Immcantation.output@meta.data$Tissue=="Tonsil"&
                                                                           Tonsil.all.Immcantation.output@meta.data$wtSpike_specific=="no"&
                                                                           !is.na(Tonsil.all.Immcantation.output@meta.data$clone_id)])
Circulation.Non.Binding <- unique(Tonsil.all.Immcantation.output@meta.data$clone_id[Tonsil.all.Immcantation.output@meta.data$Tissue=="PBMC"&
                                                                            Tonsil.all.Immcantation.output@meta.data$wtSpike_specific=="no"&
                                                                            !is.na(Tonsil.all.Immcantation.output@meta.data$clone_id)])
# Plotting -> https://r-graph-gallery.com/14-venn-diagramm.html
setwd("~/NAS/Jan/Experiments and Data/210801 - Exp 035 - Profiling of COVID vaccinated patient samples/Review nature immunology/Plots")
venn.diagram(
  x = list(First.vector, Second.vector, Third.vector,Fourth.vector),
  category.names = c("Tonsil\n spike specific" , "Circulating\n spike specific" , "Tonsil\n non spike specific","Circulating\n non spike specific"), filename = 'Clonal.overlap.Tonsils.png', output=T,
  imagetype = "png", height = 3000, width = 3000,fill=c("#DBC35E","#228833","#4477AA","#EE6677"),fontface = "bold",cat.fontface = "bold")


# UpSet Plot 
lt = list(Tonsil.Binding = Tonsil.Binding,
          Circulation.Binding = Circulation.Binding,
          Tonsil.Non.Binding=Tonsil.Non.Binding,
          Circulation.Non.Binding=Circulation.Non.Binding)
m <- make_comb_mat(lt)
set_name(m) <- c("Tonsil Binding","Circulation Binding", "Tonsil Non Binding", "Circulation Non Binding")
UpSet(m)
UpSet(m[comb_size(m) <= 50],top_annotation = upset_top_annotation(m[comb_size(m) <= 50], add_numbers = TRUE),
      right_annotation = upset_right_annotation(m[comb_size(m) <= 50], add_numbers = TRUE))


# Doing exactly the same for the 4 patients individually
# Make 4 vectors, one for each of the 4 categories that will be represented in the venn diagram
# 1: Spike specific clones that are found in the Tonsil
# 2: Spike specific clones that are found in the Circulation
# 3: Non-spike specific clones that are found in the Tonsil
# 4: Non-spike specific clones that are found in the Circulation
unique(Tonsil@meta.data$Patient)
Tonsil.all.Immcantation.output.Patientsubset <- subset(Tonsil.all.Immcantation.output, subset=Patient=="254869-1996")
First.vector <- unique(Tonsil.all.Immcantation.output.Patientsubset@meta.data$clone_id[Tonsil.all.Immcantation.output.Patientsubset@meta.data$Tissue=="Tonsil"&
                                                                           Tonsil.all.Immcantation.output.Patientsubset@meta.data$wtSpike_specific=="yes"&
                                                                           !is.na(Tonsil.all.Immcantation.output.Patientsubset@meta.data$clone_id)])
Second.vector <- unique(Tonsil.all.Immcantation.output.Patientsubset@meta.data$clone_id[Tonsil.all.Immcantation.output.Patientsubset@meta.data$Tissue=="PBMC"&
                                                                            Tonsil.all.Immcantation.output.Patientsubset@meta.data$wtSpike_specific=="yes"&
                                                                            !is.na(Tonsil.all.Immcantation.output.Patientsubset@meta.data$clone_id)])
Third.vector <- unique(Tonsil.all.Immcantation.output.Patientsubset@meta.data$clone_id[Tonsil.all.Immcantation.output.Patientsubset@meta.data$Tissue=="Tonsil"&
                                                                           Tonsil.all.Immcantation.output.Patientsubset@meta.data$wtSpike_specific=="no"&
                                                                           !is.na(Tonsil.all.Immcantation.output.Patientsubset@meta.data$clone_id)])
Fourth.vector <- unique(Tonsil.all.Immcantation.output.Patientsubset@meta.data$clone_id[Tonsil.all.Immcantation.output.Patientsubset@meta.data$Tissue=="PBMC"&
                                                                            Tonsil.all.Immcantation.output.Patientsubset@meta.data$wtSpike_specific=="no"&
                                                                            !is.na(Tonsil.all.Immcantation.output.Patientsubset@meta.data$clone_id)])
# Plotting -> https://r-graph-gallery.com/14-venn-diagramm.html
setwd("~/NAS/Jan/Experiments and Data/210801 - Exp 035 - Profiling of COVID vaccinated patient samples/Review nature immunology/Plots")
venn.diagram(
  x = list(First.vector, Second.vector, Third.vector,Fourth.vector),
  category.names = c("Tonsil\n spike specific" , "Circulating\n spike specific" , "Tonsil\n non spike specific","Circulating\n non spike specific"), filename = 'Clonal.overlap.Tonsils.254869-1996.png', output=T,
  imagetype = "png", height = 3000, width = 3000,fill=c("#DBC35E","#228833","#4477AA","#EE6677"),fontface = "bold",cat.fontface = "bold")



# Stacked Barplot - Cells belonging to each Tissue per Cluster
Idents(Tonsil)
Make.stacked.percentage.plot(Seurat.Object = Tonsil,groups.on.Y.axis = "seurat_clusters",groups.in.legend = "Tissue",
                             Legend.title = "Tissue",X.Axis.title = "Tissue",Y.Axis.title = "Clusters",Plot.title = "Tissues across Clusters") # 5x7

# Stacked Barplot - WtSpike specific cells per Cluster
Make.stacked.percentage.plot(Seurat.Object = Tonsil,groups.on.Y.axis = "seurat_clusters",groups.in.legend = "wtSpike_specific",
                             Legend.title = "Spike Binding",X.Axis.title = "Spike Binding",Y.Axis.title = "Clusters",Plot.title = "WtSpike specific cells per Cluster") # 5x7


# Isotypes of wtSpike specific cells per Tissue
Tonsil.subset <- subset(Tonsil, subset=wtSpike_specific=="yes")
Kick <- rownames(Tonsil.subset@meta.data[Tonsil.subset@meta.data$c_call=="Unknown",])
Tonsil.subset <- subset(Tonsil.subset, subset=Row.names %notin% Kick)
Make.stacked.percentage.plot(Seurat.Object = Tonsil.subset,groups.on.Y.axis = "Tissue",groups.in.legend = "c_call",
                             Legend.title = "Isotype",X.Axis.title = "Isotype",Y.Axis.title = "Tissues",
                             Plot.title = "Isotypes of wtSpike binding cells per Tissue origin") # 4.5x6
remove(Kick, Tonsil.subset)


# Showing in which clusters the spike binding cells of each tissue are.
Tonsil.subset <- subset(Tonsil, subset=wtSpike_specific=="yes")
Make.stacked.percentage.plot(Seurat.Object = Tonsil.subset, groups.on.Y.axis = "Tissue",groups.in.legend = "seurat_clusters",
                             X.Axis.title = "Cluster",Y.Axis.title = "Tissues",Legend.title = "Cluster",Plot.title = "Cluster distribution of \n Spike binding cells in Tonsil and Circulation")
remove(Tonsil.subset)



# How do the cross tissue clones distribute among the patients?
Tonsil.subset <- subset(Tonsil, subset=cross.tissue.clone=="yes")
df <- Tonsil.subset@meta.data[,c("Patient","clone_id")]
df <- unique(df)
table(df$Patient)




######################################################################################################################
# Part 10: DGEA
######################################################################################################################
# Differences between spike binders in circulation vs. Tonsil
Tonsil.subset <- subset(Tonsil, subset=wtSpike_specific=="yes")
Idents(Tonsil.subset) <- "Tissue"
Markers <- FindAllMarkers(Tonsil.subset, assay = "RNA")
Markers <- Markers[Markers$cluster=="Tonsil",]
Markers <- Markers[Markers$p_val_adj<1,]

write.xlsx(Markers, "Tonsil.vs.Circulating.DGEA.xlsx")
#Markers_2 <- Markers[Markers$avg_log2FC>2,]
#write.xlsx(Markers_2, "Circulating.vs.Tonsil.DGEA.avg_log2FC.bigger.than.two.xlsx")

# Volcano Plot
EnhancedVolcano(Markers, lab = rownames(Markers),x = 'avg_log2FC',y = 'p_val_adj', title = 'Tonsil vs Circulating',
                pCutoff = 0.05, FCcutoff = 0.3, selectLab = c("MS4A1","CD19","CR2","CCR6","CD44",
                              "KLF2","CD69","FGR","NOTCH2","CD52","LTB"),
                boxedLabels = F, drawConnectors = TRUE, legendPosition = 'none', subtitle="", labSize=4,
                arrowheads = F, col = c("black", "black", "black", "red2"), caption = "", widthConnectors = 0.8)+
  black.axis.text()+center.title() +xlim(-2,2)# Export 6x6


remove(Markers,Markers_2)

