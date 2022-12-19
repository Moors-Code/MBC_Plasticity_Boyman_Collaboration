# This script is used to analyze the scRNAseq data from the MBC project investigating the vaccinated patient cohort.

# Loading necessary packagess
library(ggforce)
library(Seurat)
library(tidyverse)
library(Biostrings)
library(alakazam)
library(shazam)
library(patchwork)
library(yarrr)
library(RColorBrewer)
library(umap)
library(ComplexHeatmap)
library(circlize)
library(cowplot)
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
      }
    }
    if(nrow(df.two)>0){
      c <- 1
      for (c in 1:nrow(df.subset_2)) {
        df.two[(c*nrow(df.subset_3)-(nrow(df.subset_3)-1)):((c*nrow(df.subset_3)-(nrow(df.subset_3)-1))+nrow(df.subset_3)-1),1] <- df.subset_2$Row.names[c]
        df.two[(c*nrow(df.subset_3)-(nrow(df.subset_3)-1)):((c*nrow(df.subset_3)-(nrow(df.subset_3)-1))+nrow(df.subset_3)-1),2] <- df.subset_3$Row.names
      }
    }
    if(nrow(df.three)>0){
      c <- 1
      for (c in 1:nrow(df.subset_3)) {
        df.three[(c*nrow(df.subset_4)-(nrow(df.subset_4)-1)):((c*nrow(df.subset_4)-(nrow(df.subset_4)-1))+nrow(df.subset_4)-1),1] <- df.subset_3$Row.names[c]
        df.three[(c*nrow(df.subset_4)-(nrow(df.subset_4)-1)):((c*nrow(df.subset_4)-(nrow(df.subset_4)-1))+nrow(df.subset_4)-1),2] <- df.subset_4$Row.names
      }
    }
    if(nrow(df.four)>0&nrow(df.one)==0){
      c <- 1
      for (c in 1:nrow(df.subset_1)) {
        df.four[(c*nrow(df.subset_3)-(nrow(df.subset_3)-1)):((c*nrow(df.subset_3)-(nrow(df.subset_3)-1))+nrow(df.subset_3)-1),1] <- df.subset_1$Row.names[c]
        df.four[(c*nrow(df.subset_3)-(nrow(df.subset_3)-1)):((c*nrow(df.subset_3)-(nrow(df.subset_3)-1))+nrow(df.subset_3)-1),2] <- df.subset_3$Row.names
      }
    }
    if(nrow(df.five)>0&nrow(df.one)==0&nrow(df.four)==0&nrow(df.two)==0){
      c <- 1
      for (c in 1:nrow(df.subset_1)) {
        df.five[(c*nrow(df.subset_4)-(nrow(df.subset_4)-1)):((c*nrow(df.subset_4)-(nrow(df.subset_4)-1))+nrow(df.subset_4)-1),1] <- df.subset_1$Row.names[c]
        df.five[(c*nrow(df.subset_4)-(nrow(df.subset_4)-1)):((c*nrow(df.subset_4)-(nrow(df.subset_4)-1))+nrow(df.subset_4)-1),2] <- df.subset_4$Row.names
      }
    }
    if(nrow(df.six)>0&nrow(df.two)==0){
      c <- 1
      for (c in 1:nrow(df.subset_2)) {
        df.six[(c*nrow(df.subset_4)-(nrow(df.subset_4)-1)):((c*nrow(df.subset_4)-(nrow(df.subset_4)-1))+nrow(df.subset_4)-1),1] <- df.subset_2$Row.names[c]
        df.six[(c*nrow(df.subset_4)-(nrow(df.subset_4)-1)):((c*nrow(df.subset_4)-(nrow(df.subset_4)-1))+nrow(df.subset_4)-1),2] <- df.subset_4$Row.names
      }
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
  
  Dimplot.colors <- c(yarrr::transparent("#a40000", trans.val = .8),yarrr::transparent("#16317d", trans.val = .8),
                      yarrr::transparent("#007e2f", trans.val = .8),yarrr::transparent("#ffcd12", trans.val = .8),
                      yarrr::transparent("#b86092", trans.val = .8),yarrr::transparent("#721b3e", trans.val = .8),
                      yarrr::transparent("#00b7a7", trans.val = .8))
  
  p <- DimPlot(Seurat.Object, reduction = reduction,group.by = DimplotGroup, cols = Dimplot.colors[1:length(unique(Seurat.Object@meta.data[,DimplotGroup]))])
  p+geom_path(data=df_final, aes(x=wnnUMAP_1, y=wnnUMAP_2, group=group),col=df_final$clonecolor,size = 0.8,arrow=arrow(type = "open", angle = 30, length = unit(0.1, "inches")))+ggtitle("Clonal paths")
  
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
# Part 1: Loading the data set, filtering, Hashing demultiplexing, addition of metadata columns
######################################################################################################################
# Loading the first dataset into R
setwd("~/NAS/Jan/Experiments and Data/220920 - Exp 040 - MBC Vacc Cohort/Pilot 072022/multi_exp035_5_vacc_prj")  
MBC_dataset <- "./outs/per_sample_outs/multi_exp035_5_vacc_prj/count/sample_feature_bc_matrix"
MBC_dataset <- Read10X(data.dir = MBC_dataset)

# Loading the second dataset into R
setwd("~/NAS/Jan/Experiments and Data/220920 - Exp 040 - MBC Vacc Cohort/Run 05102022/multi_exp040_MBC_vacc_cohort")  
MBC_dataset_2 <- "./outs/per_sample_outs/multi_exp040_MBC_vacc_cohort/count/sample_feature_bc_matrix"
MBC_dataset_2 <- Read10X(data.dir = MBC_dataset_2)
# The following line is needed because I did a mistake in the preprocessing of the data where in the "feature_reference" sheet I did not specify the sample id correctly.
rownames(MBC_dataset_2$`Antibody Capture`)[5:8] <- c("Patient_3_Post_2nd_Vacc","Patient_3_6M_post_vaccination","Patient_3_Post_Boost","Patient_3_Post_Reinfection")

# Seurat object is created
MBC <- CreateSeuratObject(counts = MBC_dataset$`Gene Expression`)
MBC_2 <- CreateSeuratObject(counts = MBC_dataset_2$`Gene Expression`)

# Assays for Baiting constructs, Protein quantification and hashing are created
Baiting_assay <- CreateAssayObject(counts = MBC_dataset$`Antibody Capture`[c(9:15),])
Baiting_assay_2 <- CreateAssayObject(counts = MBC_dataset_2$`Antibody Capture`[c(9:15),])
Protein_assay <- CreateAssayObject(counts = MBC_dataset$`Antibody Capture`[c(16:21),])
Protein_assay_2 <- CreateAssayObject(counts = MBC_dataset_2$`Antibody Capture`[c(16:21),])
Hashing_assay <- CreateAssayObject(counts = MBC_dataset$`Antibody Capture`[c(1:8),])
Hashing_assay_2 <- CreateAssayObject(counts = MBC_dataset_2$`Antibody Capture`[c(1:8),])

# Now the assays are added to the previously created Seurat object
MBC[["Baiting"]] <- Baiting_assay
MBC[["Protein"]] <- Protein_assay
MBC[["Hashing"]] <- Hashing_assay
MBC_2[["Baiting"]] <- Baiting_assay_2
MBC_2[["Protein"]] <- Protein_assay_2
MBC_2[["Hashing"]] <- Hashing_assay_2

# Removal of assay objects
remove(Baiting_assay,Hashing_assay,Protein_assay,MBC_dataset,Baiting_assay_2,Hashing_assay_2,Protein_assay_2,MBC_dataset_2)

# QC and selecting cells for further analysis - first dataset
MBC[["percent.mt"]] <- PercentageFeatureSet(MBC, pattern = "^MT-")
VlnPlot(MBC, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
quantile(MBC@meta.data$nFeature_RNA,probs = seq(0, 1, 1/20))
quantile(MBC@meta.data$percent.mt,probs = seq(0, 1, 1/20))
nrow(MBC@meta.data[MBC@meta.data$percent.mt < 10 &
                     MBC@meta.data$nFeature_RNA > 200 &
                     MBC@meta.data$nFeature_RNA < 4000,])/nrow(MBC@meta.data)

MBC <- subset(MBC, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 10)

# QC and selecting cells for further analysis - second dataset
MBC_2[["percent.mt"]] <- PercentageFeatureSet(MBC_2, pattern = "^MT-")
VlnPlot(MBC_2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
quantile(MBC_2@meta.data$nFeature_RNA,probs = seq(0, 1, 1/20))
quantile(MBC_2@meta.data$percent.mt,probs = seq(0, 1, 1/20))
nrow(MBC_2@meta.data[MBC_2@meta.data$percent.mt < 10 &
                     MBC_2@meta.data$nFeature_RNA > 200 &
                     MBC_2@meta.data$nFeature_RNA < 4000,])/nrow(MBC_2@meta.data)

MBC_2 <- subset(MBC_2, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 10)

# Demultiplexing the HTO data
MBC <- NormalizeData(MBC, assay = "Hashing", normalization.method = "CLR", margin = 2)
MBC_2 <- NormalizeData(MBC_2, assay = "Hashing", normalization.method = "CLR", margin = 2)
MBC <- HTODemux(MBC, assay = "Hashing", positive.quantile = 0.99)
MBC_2 <- HTODemux(MBC_2, assay = "Hashing", positive.quantile = 0.99)
table(MBC@meta.data$Hashing_classification.global)
table(MBC_2@meta.data$Hashing_classification.global)
MBC <- subset(MBC, subset = Hashing_classification.global=="Singlet")
MBC_2 <- subset(MBC_2, subset = Hashing_classification.global=="Singlet")
RidgePlot(MBC, assay = "Hashing", features = rownames(MBC[["Hashing"]])[1:10], ncol = 2)
RidgePlot(MBC_2, assay = "Hashing", features = rownames(MBC_2[["Hashing"]])[1:10], ncol = 2)
table(MBC@meta.data$Hashing_classification)
table(MBC_2@meta.data$Hashing_classification)

# Introducing metadata columns
MBC@meta.data$Patient <- "Patient 1"
MBC@meta.data$Patient[grep("Patient-2",MBC@meta.data$hash.ID)] <- "Patient 2"
MBC@meta.data$Timepoint <- "6M post 2nd vacc"
MBC@meta.data$Timepoint[grep("Post-2nd",MBC@meta.data$hash.ID)] <- "Post 2nd vacc"
MBC@meta.data$Timepoint[grep("Post-Boost",MBC@meta.data$hash.ID)] <- "Post boost vacc"
MBC@meta.data$Timepoint[grep("Post-Reinfection",MBC@meta.data$hash.ID)] <- "Post reinfection"

MBC_2@meta.data$Patient <- "Patient 1"
MBC_2@meta.data$Patient[grep("Patient-3",MBC_2@meta.data$hash.ID)] <- "Patient 3"
MBC_2@meta.data$Timepoint <- "6M post 2nd vacc"
MBC_2@meta.data$Timepoint[grep("Post-2nd",MBC_2@meta.data$hash.ID)] <- "Post 2nd vacc"
MBC_2@meta.data$Timepoint[grep("Post-Boost",MBC_2@meta.data$hash.ID)] <- "Post boost vacc"
MBC_2@meta.data$Timepoint[grep("Post-Reinfection",MBC_2@meta.data$hash.ID)] <- "Post reinfection"


######################################################################################################################
# Part 2: Processing of baiting counts (cutoffs, normalization) and visualization
######################################################################################################################
# Adding metadata columns to the Seurat object which will indicate if a cell is binding an antigen
MBC@meta.data$Nucleocapsid_specific <- "no"
MBC@meta.data$wtSpike_specific <- "no"
MBC@meta.data$RBD_specific <- "no"
MBC@meta.data$B.1.351_specific <- "no"
MBC@meta.data$B.1.617.2_specific <- "no"
MBC@meta.data$B.1.529_specific <- "no"
MBC_2@meta.data$Nucleocapsid_specific <- "no"
MBC_2@meta.data$wtSpike_specific <- "no"
MBC_2@meta.data$RBD_specific <- "no"
MBC_2@meta.data$B.1.351_specific <- "no"
MBC_2@meta.data$B.1.617.2_specific <- "no"
MBC_2@meta.data$B.1.529_specific <- "no"

# Pre-processing baiting counts
Baiting_df <- as.data.frame(MBC@assays$Baiting@counts)
Baiting_df <- t(Baiting_df)
Baiting_df <- as.data.frame(Baiting_df)
Baiting_df_2 <- as.data.frame(MBC_2@assays$Baiting@counts)
Baiting_df_2 <- t(Baiting_df_2)
Baiting_df_2 <- as.data.frame(Baiting_df_2)

# Subtracting the negative control counts
Baiting_df$wtSpike <- Baiting_df$wtSpike-Baiting_df$NegControl
Baiting_df$B.1.529 <- Baiting_df$B.1.529-Baiting_df$NegControl
Baiting_df$Nucleocapsid <- Baiting_df$Nucleocapsid-Baiting_df$NegControl
Baiting_df$RBD <- Baiting_df$RBD-Baiting_df$NegControl
Baiting_df$B.1.351 <- Baiting_df$B.1.351-Baiting_df$NegControl
Baiting_df$B.1.617.2 <- Baiting_df$B.1.617.2-Baiting_df$NegControl
Baiting_df_2$wtSpike <- Baiting_df_2$wtSpike-Baiting_df_2$NegControl
Baiting_df_2$B.1.529 <- Baiting_df_2$B.1.529-Baiting_df_2$NegControl
Baiting_df_2$Nucleocapsid <- Baiting_df_2$Nucleocapsid-Baiting_df_2$NegControl
Baiting_df_2$RBD <- Baiting_df_2$RBD-Baiting_df_2$NegControl
Baiting_df_2$B.1.351 <- Baiting_df_2$B.1.351-Baiting_df_2$NegControl
Baiting_df_2$B.1.617.2 <- Baiting_df_2$B.1.617.2-Baiting_df_2$NegControl

# Setting negative bait counts to 0
Baiting_df$wtSpike[Baiting_df$wtSpike<0] <- 0
Baiting_df$B.1.529[Baiting_df$B.1.529<0] <- 0
Baiting_df$Nucleocapsid[Baiting_df$Nucleocapsid<0] <- 0
Baiting_df$RBD[Baiting_df$RBD<0] <- 0
Baiting_df$B.1.351[Baiting_df$B.1.351<0] <- 0
Baiting_df$B.1.617.2[Baiting_df$B.1.617.2<0] <- 0
Baiting_df_2$wtSpike[Baiting_df_2$wtSpike<0] <- 0
Baiting_df_2$B.1.529[Baiting_df_2$B.1.529<0] <- 0
Baiting_df_2$Nucleocapsid[Baiting_df_2$Nucleocapsid<0] <- 0
Baiting_df_2$RBD[Baiting_df_2$RBD<0] <- 0
Baiting_df_2$B.1.351[Baiting_df_2$B.1.351<0] <- 0
Baiting_df_2$B.1.617.2[Baiting_df_2$B.1.617.2<0] <- 0

# Density plot of raw data
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
MBC <- Cutoff_function(Seurat_object = MBC, bait="wtSpike",left_border = 5,right_border = 40,melted_df = melted_Baiting_df)
MBC <- Cutoff_function(Seurat_object = MBC, bait="RBD",left_border = 10,right_border = 500,melted_df = melted_Baiting_df)
MBC <- Cutoff_function(Seurat_object = MBC, bait="B.1.529",left_border = 5,right_border = 40,melted_df = melted_Baiting_df)
MBC <- Cutoff_function(Seurat_object = MBC, bait="B.1.351",left_border = 5,right_border = 40,melted_df = melted_Baiting_df)
MBC <- Cutoff_function(Seurat_object = MBC, bait="B.1.617.2",left_border = 5,right_border = 100,melted_df = melted_Baiting_df)
MBC <- Cutoff_function(Seurat_object = MBC, bait="Nucleocapsid",left_border = 1,right_border = 10,melted_df = melted_Baiting_df)

MBC_2 <- Cutoff_function(Seurat_object = MBC_2, bait="wtSpike",left_border = 2,right_border = 10,melted_df = melted_Baiting_df_2)
MBC_2 <- Cutoff_function(Seurat_object = MBC_2, bait="RBD",left_border = 10,right_border = 200,melted_df = melted_Baiting_df_2)
MBC_2 <- Cutoff_function(Seurat_object = MBC_2, bait="B.1.529",left_border = 5,right_border = 20,melted_df = melted_Baiting_df_2)
MBC_2 <- Cutoff_function(Seurat_object = MBC_2, bait="B.1.351",left_border = 5,right_border = 40,melted_df = melted_Baiting_df_2)
MBC_2 <- Cutoff_function(Seurat_object = MBC_2, bait="B.1.617.2",left_border = 5,right_border = 100,melted_df = melted_Baiting_df_2)
MBC_2 <- Cutoff_function(Seurat_object = MBC_2, bait="Nucleocapsid",left_border = 1,right_border = 10,melted_df = melted_Baiting_df_2)

remove(melted_Baiting_df,melted_Baiting_df_2)

# Setting baiting counts below determined thresholds to NA
Baiting_df$wtSpike[Baiting_df$wtSpike<18] <- NA
Baiting_df$RBD[Baiting_df$RBD<64] <- NA
Baiting_df$B.1.529[Baiting_df$B.1.529<12] <- NA
Baiting_df$B.1.351[Baiting_df$B.1.351<14] <- NA
Baiting_df$B.1.617.2[Baiting_df$B.1.617.2<16] <- NA
Baiting_df$Nucleocapsid[Baiting_df$Nucleocapsid<4] <- NA
Baiting_df_2$wtSpike[Baiting_df_2$wtSpike<4] <- NA
Baiting_df_2$RBD[Baiting_df_2$RBD<79] <- NA
Baiting_df_2$B.1.529[Baiting_df_2$B.1.529<9] <- NA
Baiting_df_2$B.1.351[Baiting_df_2$B.1.351<10] <- NA
Baiting_df_2$B.1.617.2[Baiting_df_2$B.1.617.2<14] <- NA
Baiting_df_2$Nucleocapsid[Baiting_df_2$Nucleocapsid<3] <- NA

# In order to use the Seurat Normalization function for corrected baiting counts, we need to have the corrected counts as assay object
# For this, the Baiting_df needs to be transposed
Baiting_df <- as.data.frame(t(Baiting_df[,c(1,2,3,4,6,7)]))
Baiting_df_2 <- as.data.frame(t(Baiting_df_2[,c(1,2,3,4,6,7)]))
Cor_Baiting_assay <- CreateAssayObject(counts = Baiting_df)
Cor_Baiting_assay_2 <- CreateAssayObject(counts = Baiting_df_2)
MBC[["Cor_Baiting"]] <- Cor_Baiting_assay
MBC_2[["Cor_Baiting"]] <- Cor_Baiting_assay_2
remove(Cor_Baiting_assay,Cor_Baiting_assay_2)
MBC <- NormalizeData(MBC, assay = "Cor_Baiting", normalization.method = 'CLR', margin = 1)
MBC_2 <- NormalizeData(MBC_2, assay = "Cor_Baiting", normalization.method = 'CLR', margin = 1)
normalized_across_features <- as.data.frame(t(MBC@assays$Cor_Baiting@data))
normalized_across_features_2 <- as.data.frame(t(MBC_2@assays$Cor_Baiting@data))
remove(Baiting_df,Baiting_df_2)

# For scaling, I need to transform again, then each column is scaled to 0:1 individually using my own scaling function.
normalized_across_features <- myscaling(normalized_across_features)  
normalized_across_features_2 <- myscaling(normalized_across_features_2)  

# These are the final LIBRA scores
# Adding this to the Seurat object
MBC@meta.data <- merge(MBC@meta.data,normalized_across_features, by=0)
MBC_2@meta.data <- merge(MBC_2@meta.data,normalized_across_features_2, by=0)
rownames(MBC@meta.data) <- MBC@meta.data$Row.names
rownames(MBC_2@meta.data) <- MBC_2@meta.data$Row.names
remove(normalized_across_features,normalized_across_features_2)

# NAs are now set to 0
MBC@meta.data$wtSpike[is.na(MBC@meta.data$wtSpike)] <- 0
MBC@meta.data$B.1.529[is.na(MBC@meta.data$B.1.529)] <- 0
MBC@meta.data$Nucleocapsid[is.na(MBC@meta.data$Nucleocapsid)] <- 0
MBC@meta.data$RBD[is.na(MBC@meta.data$RBD)] <- 0
MBC@meta.data$B.1.351[is.na(MBC@meta.data$B.1.351)] <- 0
MBC@meta.data$B.1.617.2[is.na(MBC@meta.data$B.1.617.2)] <- 0
MBC_2@meta.data$wtSpike[is.na(MBC_2@meta.data$wtSpike)] <- 0
MBC_2@meta.data$B.1.529[is.na(MBC_2@meta.data$B.1.529)] <- 0
MBC_2@meta.data$Nucleocapsid[is.na(MBC_2@meta.data$Nucleocapsid)] <- 0
MBC_2@meta.data$RBD[is.na(MBC_2@meta.data$RBD)] <- 0
MBC_2@meta.data$B.1.351[is.na(MBC_2@meta.data$B.1.351)] <- 0
MBC_2@meta.data$B.1.617.2[is.na(MBC_2@meta.data$B.1.617.2)] <- 0

# FeatureScatter plots for first dataset
FeatureScatter(MBC,  feature1 = "wtSpike", feature2 = "RBD", pt.size = 1, group.by = "Patient")+
  ggtitle("Wt Spike and RBD scores")+xlab("wt-Spike scores")+ylab("RBD scores")+
  
  FeatureScatter(MBC,  feature1 = "wtSpike", feature2 = "B.1.351", pt.size = 1, group.by = "Patient")+
  ggtitle("Wt Spike and Beta variant scores")+xlab("wt-Spike scores")+ylab("Beta variant scores")+
  
  FeatureScatter(MBC,  feature1 = "wtSpike", feature2 = "Nucleocapsid", pt.size = 1, group.by = "Patient")+
  ggtitle("Wt Spike and Nucleocapsid scores")+xlab("wt-Spike scores")+ylab("Nucleocapsid scores")+
  
  FeatureScatter(MBC,  feature1 = "wtSpike", feature2 = "B.1.529", pt.size = 1, group.by = "Patient")+
  ggtitle("Beta and Omciron variant scores")+xlab("wt-Spike scores")+ylab("Omicron variant scores")

FeatureScatter(MBC,  feature1 = "wtSpike", feature2 = "RBD", pt.size = 1, group.by = "Patient")+
  ggtitle("Wt Spike and RBD scores")+xlab("wt-Spike scores")+ylab("RBD scores")+
  
  FeatureScatter(MBC,  feature1 = "wtSpike", feature2 = "B.1.351", pt.size = 1, group.by = "Patient")+
  ggtitle("Wt Spike and Beta variant scores")+xlab("wt-Spike scores")+ylab("Beta variant scores")+
  
  FeatureScatter(MBC,  feature1 = "wtSpike", feature2 = "B.1.617.2", pt.size = 1, group.by = "Patient")+
  ggtitle("Wt Spike and Delta variant scores")+xlab("wt-Spike scores")+ylab("Delta variant scores")+
  
  FeatureScatter(MBC,  feature1 = "wtSpike", feature2 = "B.1.529", pt.size = 1, group.by = "Patient")+
  ggtitle("Beta and Omciron variant scores")+xlab("wt-Spike scores")+ylab("Omicron variant scores")

# FeatureScatter plots for second ataset
FeatureScatter(MBC_2,  feature1 = "wtSpike", feature2 = "RBD", pt.size = 1, group.by = "Patient")+
  ggtitle("Wt Spike and RBD scores")+xlab("wt-Spike scores")+ylab("RBD scores")+
  
  FeatureScatter(MBC_2,  feature1 = "wtSpike", feature2 = "B.1.351", pt.size = 1, group.by = "Patient")+
  ggtitle("Wt Spike and Beta variant scores")+xlab("wt-Spike scores")+ylab("Beta variant scores")+
  
  FeatureScatter(MBC_2,  feature1 = "wtSpike", feature2 = "Nucleocapsid", pt.size = 1, group.by = "Patient")+
  ggtitle("Wt Spike and Nucleocapsid scores")+xlab("wt-Spike scores")+ylab("Nucleocapsid scores")+
  
  FeatureScatter(MBC_2,  feature1 = "wtSpike", feature2 = "B.1.529", pt.size = 1, group.by = "Patient")+
  ggtitle("Beta and Omciron variant scores")+xlab("wt-Spike scores")+ylab("Omicron variant scores")

FeatureScatter(MBC_2,  feature1 = "wtSpike", feature2 = "RBD", pt.size = 1, group.by = "Patient")+
  ggtitle("Wt Spike and RBD scores")+xlab("wt-Spike scores")+ylab("RBD scores")+
  
  FeatureScatter(MBC_2,  feature1 = "wtSpike", feature2 = "B.1.351", pt.size = 1, group.by = "Patient")+
  ggtitle("Wt Spike and Beta variant scores")+xlab("wt-Spike scores")+ylab("Beta variant scores")+
  
  FeatureScatter(MBC_2,  feature1 = "wtSpike", feature2 = "B.1.617.2", pt.size = 1, group.by = "Patient")+
  ggtitle("Wt Spike and Delta variant scores")+xlab("wt-Spike scores")+ylab("Delta variant scores")+
  
  FeatureScatter(MBC_2,  feature1 = "wtSpike", feature2 = "B.1.529", pt.size = 1, group.by = "Patient")+
  ggtitle("Beta and Omciron variant scores")+xlab("wt-Spike scores")+ylab("Omicron variant scores")



######################################################################################################################
# Part 3: Dataset integration and WNN
######################################################################################################################
# Renaming of cells
MBC <- RenameCells(object = MBC, add.cell.id = "Dataset_1")
MBC_2 <- RenameCells(object = MBC_2, add.cell.id = "Dataset_2")
MBC@meta.data$Dataset <- "First"
MBC_2@meta.data$Dataset <- "Second"

# During the IntegrateData() function, dimensional reductions of the data are needed. PCA is calculated before Integration is done (only for the protein data)
DefaultAssay(MBC) <- "Protein"
DefaultAssay(MBC_2) <- "Protein"
VariableFeatures(MBC) <- rownames(MBC[["Protein"]])
VariableFeatures(MBC_2) <- rownames(MBC_2[["Protein"]])
MBC <- NormalizeData(MBC, normalization.method = 'CLR', margin = 2) %>% ScaleData() %>% RunPCA(reduction.name = 'apca',approx=FALSE) # approx=F is used to remove the warning that comes otherwise
MBC_2 <- NormalizeData(MBC_2, normalization.method = 'CLR', margin = 2) %>% ScaleData() %>% RunPCA(reduction.name = 'apca',approx=FALSE)

# Create a list containing the two datasets
MBC_dataset_list <- list(MBC,MBC_2)

# Select features that are repeatedly variable across datasets for integration
# The integration is done independently for the RNA and the Protein assay.
features_RNA <- SelectIntegrationFeatures(object.list = MBC_dataset_list,assay = c("RNA","RNA")) # IG variable regions pop up!
features_Protein <- SelectIntegrationFeatures(object.list = MBC_dataset_list,assay = c("Protein","Protein"))

# Remove HLA, Immunoglobulin, RNA, MT, and RP genes based on HUGO gene names
var_regex = '^HLA-|^IG[HJKL]|^RNA|^MT|^RP'
features_RNA <- grep(var_regex, features_RNA, invert=T, value=T)
remove(var_regex)

# Perform integration. This command requires long time / a lot of computational resources to run.
# Because of only 6 protein features, dims is set to 5 for the anchor finding of protein data
anchors_RNA <- FindIntegrationAnchors(object.list = MBC_dataset_list, anchor.features = features_RNA,assay = c("RNA","RNA"))
anchors_Protein <- FindIntegrationAnchors(object.list = MBC_dataset_list, anchor.features = features_Protein,assay = c("Protein","Protein"),dims=1:5)

# This command creates an 'integrated' data assay in a new Seurat object.
MBC.combined.RNA <- IntegrateData(anchorset = anchors_RNA,new.assay.name = "integratedRNA")
MBC.combined.RNA <- FindVariableFeatures(MBC.combined.RNA, assay = "Protein")
MBC.combined.RNA <- ScaleData(MBC.combined.RNA, assay = "Protein")
MBC.combined.RNA <- RunPCA(MBC.combined.RNA, assay = "Protein", reduction.name="pca_protein_all_cells",approx=FALSE)
MBC.combined.Protein <- IntegrateEmbeddings(anchorset = anchors_Protein, reductions = MBC.combined.RNA@reductions$pca_protein_all_cells)
remove(anchors_Protein,anchors_RNA,MBC_dataset_list,MBC,MBC_2)

# Now WNN
# For this, I need to transfer the integrated protein dim reduction from the Bmem.combined.Protein object to the Bmem.combined.RNA object.
# But first, I am preforming data scaling and PCA on the combined RNA project.
DefaultAssay(MBC.combined.RNA)
DefaultAssay(MBC.combined.Protein)

MBC.combined.RNA <- ScaleData(MBC.combined.RNA, verbose = FALSE)
MBC.combined.RNA <- RunPCA(MBC.combined.RNA, npcs = 30, verbose = FALSE)

# Now we transfer the integrated protein dim reduction
MBC.combined.RNA@reductions$IntegrateEmbeddings.pca <- MBC.combined.Protein@reductions$integrated_dr

# I double check that the correct assays are associated with the correct PCAs.
MBC.combined.RNA[['pca']]@assay.used
MBC.combined.RNA[['IntegrateEmbeddings.pca']]@assay.used


# Finally, the FindMultiModalNeighbors function is called, using the PCAs from above.
ElbowPlot(MBC.combined.RNA, reduction="pca")
MBC.combined.RNA <- FindMultiModalNeighbors(
  MBC.combined.RNA, reduction.list = list("pca", "IntegrateEmbeddings.pca"), 
  dims.list = list(1:15, 1:5), modality.weight.name = c("RNA.weight","Protein.weight"))

MBC <- MBC.combined.RNA
MBC@meta.data$Full.Row.names <- rownames(MBC@meta.data)
DefaultAssay(MBC) <- "RNA"
MBC <- ScaleData(MBC)
remove(features_Protein,features_RNA,MBC.combined.Protein,MBC.combined.RNA)

MBC@meta.data$Timepoint <- factor(MBC@meta.data$Timepoint, levels = c("Post 2nd vacc", "6M post 2nd vacc", "Post boost vacc","Post reinfection"))

######################################################################################################################
# Part 4: Creating Immcantation Input files (Per Patient)
######################################################################################################################
# Prepare the .txt file per patient
Patients <- unique(MBC@meta.data$Patient)
MBC@meta.data$Row.names <- rownames(MBC@meta.data)
Patient_1_positives <- MBC@meta.data$Row.names[MBC@meta.data$Patient == Patients[2] & MBC@meta.data$wtSpike_specific=="yes"] 
Patient_2_positives <- MBC@meta.data$Row.names[MBC@meta.data$Patient == Patients[1] & MBC@meta.data$wtSpike_specific=="yes"]
Patient_3_positives <- MBC@meta.data$Row.names[MBC@meta.data$Patient == Patients[3] & MBC@meta.data$wtSpike_specific=="yes"]
all.cells <- rownames(MBC@meta.data)

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

Patient_1_positives <- Immcantation.Input.Barcode.Function(Patient_1_positives)
Patient_2_positives <- Immcantation.Input.Barcode.Function(Patient_2_positives)
Patient_3_positives <- Immcantation.Input.Barcode.Function(Patient_3_positives)
all.cells <- Immcantation.Input.Barcode.Function(all.cells)

write(Patient_1_positives, file = "~/NAS/Jan/Experiments and Data/220920 - Exp 040 - MBC Vacc Cohort/Immcantation/Input_data/Patient_1_positives.txt")
write(Patient_2_positives, file = "~/NAS/Jan/Experiments and Data/220920 - Exp 040 - MBC Vacc Cohort/Immcantation/Input_data/Patient_2_positives.txt")
write(Patient_3_positives, file = "~/NAS/Jan/Experiments and Data/220920 - Exp 040 - MBC Vacc Cohort/Immcantation/Input_data/Patient_3_positives.txt")
write(all.cells, file = "~/NAS/Jan/Experiments and Data/220920 - Exp 040 - MBC Vacc Cohort/Immcantation/Input_data/all.cells.txt")

# Prepare the contig.csv file per patient
contigs_Dataset_1 <- read.csv("~/NAS/Jan/Experiments and Data/220920 - Exp 040 - MBC Vacc Cohort/Pilot 072022/multi_exp035_5_vacc_prj/outs/per_sample_outs/multi_exp035_5_vacc_prj/vdj_b/filtered_contig_annotations.csv")
contigs_Dataset_2 <- read.csv("~/NAS/Jan/Experiments and Data/220920 - Exp 040 - MBC Vacc Cohort/Run 05102022/multi_exp040_MBC_vacc_cohort/outs/per_sample_outs/multi_exp040_MBC_vacc_cohort/vdj_b/filtered_contig_annotations.csv")

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
contigs_Patient_1_positives <- all.contigs[all.contigs$barcode %in% Patient_1_positives,]
contigs_Patient_2_positives <- all.contigs[all.contigs$barcode %in% Patient_2_positives,]
contigs_Patient_3_positives <- all.contigs[all.contigs$barcode %in% Patient_3_positives,]
contigs_all_cells <- all.contigs[all.contigs$barcode %in% all.cells,]

write_csv(contigs_Patient_1_positives, file = "~/NAS/Jan/Experiments and Data/220920 - Exp 040 - MBC Vacc Cohort/Immcantation/Input_data/Patient_1_positives.csv")
write_csv(contigs_Patient_2_positives, file = "~/NAS/Jan/Experiments and Data/220920 - Exp 040 - MBC Vacc Cohort/Immcantation/Input_data/Patient_2_positives.csv")
write_csv(contigs_Patient_3_positives, file = "~/NAS/Jan/Experiments and Data/220920 - Exp 040 - MBC Vacc Cohort/Immcantation/Input_data/Patient_3_positives.csv")
write_csv(contigs_all_cells, file = "~/NAS/Jan/Experiments and Data/220920 - Exp 040 - MBC Vacc Cohort/Immcantation/Input_data/all.cells.csv")

remove(Patient_1_positives,Patient_2_positives,Patient_3_positives,all.cells,contigs_all_cells,Patients,
       all.contigs,contigs_Patient_1_positives,contigs_Patient_2_positives,contigs_Patient_3_positives)

######################################################################################################################
# Part 5: Immcantation Output analysis including Donut charts (not in Manuskript)
######################################################################################################################
# Loading the output files
Patient_1.Immcantation_output <- as.data.frame(read_airr("~/NAS/Jan/Experiments and Data/220920 - Exp 040 - MBC Vacc Cohort/Immcantation/Output_data/Patient_1_positives/outputs/Patient_1_positives_0.20/Patient_1_positives_0.20_heavy_germ-pass.tsv"))
Patient_2.Immcantation_output <- as.data.frame(read_airr("~/NAS/Jan/Experiments and Data/220920 - Exp 040 - MBC Vacc Cohort/Immcantation/Output_data/Patient_2_positives/outputs/Patient_2_positives_0.20/Patient_2_positives_0.20_heavy_germ-pass.tsv"))
Patient_3.Immcantation_output <- as.data.frame(read_airr("~/NAS/Jan/Experiments and Data/220920 - Exp 040 - MBC Vacc Cohort/Immcantation/Output_data/Patient_3_positives/outputs/Patient_3_positives_0.20/Patient_3_positives_0.20_heavy_germ-pass.tsv"))
All.cells <- read.table("~/NAS/Jan/Experiments and Data/220920 - Exp 040 - MBC Vacc Cohort/Immcantation/Output_data/all.cells/outputs/all.cells_0.20/db_d-masked_germ-pass.tsv", sep = "\t", header = T)

# Taking care of duplicate row names in the All.cells file
All.cells <- do.call("rbind", by(All.cells, All.cells$cell_id, function(x) x[which.max(x$umi_count), ]))

# Translating the cell barcodes back
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
Patient_1.Immcantation_output$cell_id <- Immcantation.Output.Barcode.Function(Patient_1.Immcantation_output$cell_id)
Patient_2.Immcantation_output$cell_id <- Immcantation.Output.Barcode.Function(Patient_2.Immcantation_output$cell_id)
Patient_3.Immcantation_output$cell_id <- Immcantation.Output.Barcode.Function(Patient_3.Immcantation_output$cell_id)

# Changing the rownames to cell_ids (not needed for All.cells file)
rownames(Patient_1.Immcantation_output) <- Patient_1.Immcantation_output$cell_id
rownames(Patient_2.Immcantation_output) <- Patient_2.Immcantation_output$cell_id
rownames(Patient_3.Immcantation_output) <- Patient_3.Immcantation_output$cell_id

# Making the clone IDs unique to patients
Patient_1.Immcantation_output$clone_id <- paste0(Patient_1.Immcantation_output$clone_id,"_1")
Patient_2.Immcantation_output$clone_id <- paste0(Patient_2.Immcantation_output$clone_id,"_2")
Patient_3.Immcantation_output$clone_id <- paste0(Patient_3.Immcantation_output$clone_id,"_3")

# Merging the Immcantation output files and adding them to the Seurat object
Immcantation_output <- rbind(Patient_1.Immcantation_output,Patient_2.Immcantation_output,Patient_3.Immcantation_output)
MBC@meta.data$Row.names <- NULL
MBC@meta.data <- merge(MBC@meta.data,Immcantation_output[,c("cdr3","clone_id")],by=0, all.x=T)
rownames(MBC@meta.data) <- MBC@meta.data$Row.names
remove(Immcantation_output)

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


# Now I am adding Isotypes and heavy chain info to the cells of the Seurat object.
MBC@meta.data$Row.names <- NULL
MBC@meta.data <- merge(MBC@meta.data,All.cells[,c("v_call","d_call","j_call","c_call","mu_count","mu_freq")],by=0,all.x=T)
rownames(MBC@meta.data) <- MBC@meta.data$Row.names
MBC@meta.data$c_call[MBC@meta.data$c_call!="IGHA1" &
                       MBC@meta.data$c_call!="IGHA2"&
                       MBC@meta.data$c_call!="IGHD"&
                       MBC@meta.data$c_call!="IGHG1"&
                       MBC@meta.data$c_call!="IGHG2"&
                       MBC@meta.data$c_call!="IGHG3"&
                       MBC@meta.data$c_call!="IGHG4"&
                       MBC@meta.data$c_call!="IGHM"&
                       !is.na(MBC@meta.data$c_call)] <- "Unknown"
MBC@meta.data$c_call[is.na(MBC@meta.data$c_call)] <- "Unknown"
remove(All.cells, obs, freq, obs_freq)

# Testing if a patient has longitudinal clones - clones that we find in more than one time point
# The for loop also adds the lonngitudinal clone info to the Seurat object meta data.
MBC.subset <- subset(MBC, subset= clone_id != is.na(clone_id))
MBC.subset.list <- MBC.subset@meta.data %>% dplyr::group_by(clone_id) %>% group_split()

i <- 1
MBC@meta.data$longitudinal.clone <- "no"
for (i in 1:length(MBC.subset.list)) {
  
  this.round.clone <- as.data.frame(MBC.subset.list[[i]])
  if(length(unique(this.round.clone$Timepoint))>1){
    print(paste0(this.round.clone$clone_id," has cells in timepoints ",unique(this.round.clone$Timepoint)))
    print(paste0(this.round.clone$clone_id," has cells in ",length(unique(this.round.clone$Timepoint))," timepoints"))
    MBC@meta.data$longitudinal.clone[MBC@meta.data$Row.names %in% MBC.subset.list[[i]]$Row.names] <- "yes"
  } else {next}
  
}

remove(MBC.subset,MBC.subset.list,this.round.clone,i)

# Check within the longitudinal clones which timepoints they cover
Longitudinals <- subset(MBC, subset=longitudinal.clone=="yes")
table(Longitudinals@meta.data$Timepoint,Longitudinals@meta.data$clone_id)

# Adding a clonal frequency column
MBC@meta.data <- MBC@meta.data %>% group_by(clone_id)%>% add_tally(name = "clonal.frequency") %>% ungroup() %>% as.data.frame()
rownames(MBC@meta.data) <- MBC@meta.data$Row.names
MBC@meta.data$clonal.frequency[is.na(MBC@meta.data$clone_id)] <- NA

# Making donut charts - one for each timepoint and patient (four per patient).
# First define the Patient, for which you want to plot the donut charts and accordingly, the title of the donut chart.
Patient <- Patient_1.Immcantation_output
Donuttitle <- "Patient 1"

# Then execute everything below here
find.longitudinal.clones <- function(Immcantation.Output){
  IDs <- unique(Immcantation.Output$clone_id)
  output.vector <- c("")
  i <- 1
  for(i in 1:length(IDs)){
    candidates <- Immcantation.Output$cell_id[Immcantation.Output$clone_id ==IDs[i]]
    timepoints <-MBC@meta.data$Timepoint[rownames(MBC@meta.data)%in%candidates]
    if(length(unique(timepoints))>1){
      output.vector <- c(output.vector,IDs[i])
    }
  }
  
  output.vector <- output.vector[2:length(output.vector)]
  return(output.vector)
  
}
Longitudinals <- find.longitudinal.clones(Patient)

annotate.clones <- function(Longitudinals,Immcantation.Output,Seuratobject){
  IDs <- Longitudinals
  output.df <- data.frame(matrix(ncol = 2, nrow = 0))
  x <- c("Clone", "TruePersistent")
  colnames(output.df) <- x
  i <- 1
  for(i in 1:length(IDs)){
    candidates <- Immcantation.Output$cell_id[Immcantation.Output$clone_id ==IDs[i]]
    Bmem.Inc <- Seuratobject
    candidates.metadata <-Bmem.Inc@meta.data[rownames(Bmem.Inc@meta.data)%in%candidates,c("Row.names","Timepoint","Patient")]
    if(nrow(candidates.metadata)>2){ # in the MBC code (Zurbuchen, Michler et al) this will be "yes" if the clone has freq >2 in one of the two timepoints.
      output.df[i,1] <- IDs[i]
      output.df[i,2] <- "yes"
    } else {
      output.df[i,1] <- IDs[i]
      output.df[i,2] <- "no"
    }
  }
  output.df <- output.df[output.df$TruePersistent=="yes",]
  Immcantation.Output <- Immcantation.Output %>% group_by(clone_id) %>% add_count() %>% ungroup()
  Immcantation.Output$persistence <- "Singlet"
  `%notin%` <- Negate(`%in%`)
  Immcantation.Output$persistence[Immcantation.Output$clone_id %notin% Longitudinals&
                                    Immcantation.Output$n >1] <- "Unique clone"
  Immcantation.Output$persistence[Immcantation.Output$clone_id %in% Longitudinals &
                                    Immcantation.Output$clone_id %notin% output.df$Clone] <- "Persistent clone" # This line together with the line below has the effect that the if else part from above (commented section) has no effect.
  Immcantation.Output$persistence[Immcantation.Output$clone_id %in% Longitudinals &
                                    Immcantation.Output$clone_id %in% output.df$Clone] <- "Persistent clone"
  return(Immcantation.Output)
  
}
Patient_1.Immcantation_output_annotated <- annotate.clones(Longitudinals, Patient, MBC)

# We will mainly use the Immcantation output file here but the file needs to be updated with a column that contains the time point
Patient.cells <- Patient$cell_id
timepoints <- MBC@meta.data[rownames(MBC@meta.data)%in%Patient.cells,c("Timepoint","Row.names")]

Patient_1.Immcantation_output_annotated <- as.data.frame(Patient_1.Immcantation_output_annotated)
rownames(Patient_1.Immcantation_output_annotated) <- Patient_1.Immcantation_output_annotated$cell_id
Patient_1.Immcantation_output_annotated <- merge(Patient_1.Immcantation_output_annotated, timepoints, by=0, all.x = T)

# Adding the colors that the pie chart segments will have later. Options are Singlet, Unique clone and persistent clone
Patient_1.Immcantation_output_annotated$color <- NA
Patient_1.Immcantation_output_annotated$color[Patient_1.Immcantation_output_annotated$persistence=="Singlet"] <- "white"
Patient_1.Immcantation_output_annotated$color[Patient_1.Immcantation_output_annotated$persistence=="Unique clone"] <- "grey"

# For the category Persistent clone we need to give colors.
Persistents <- unique(Patient_1.Immcantation_output_annotated$clone_id[Patient_1.Immcantation_output_annotated$persistence=="Persistent clone"])

qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col_vector <- unique(col_vector)
col_vector <- col_vector[1:length(Persistents)]
df <- data.frame("clone_id"= Persistents, "color"=col_vector)

Patient_1.Immcantation_output_annotated <- merge(Patient_1.Immcantation_output_annotated, df, by="clone_id", all.x = T)
Patient_1.Immcantation_output_annotated$color <- coalesce(Patient_1.Immcantation_output_annotated$color.x,Patient_1.Immcantation_output_annotated$color.y)
Patient_1.Immcantation_output_annotated$color.x <- NULL
Patient_1.Immcantation_output_annotated$color.y <- NULL

Patient_Post.2nd.vacc <- Patient_1.Immcantation_output_annotated[Patient_1.Immcantation_output_annotated$Timepoint=="Post 2nd vacc",]
Patient_6M.post.2nd.vacc <- Patient_1.Immcantation_output_annotated[Patient_1.Immcantation_output_annotated$Timepoint=="6M post 2nd vacc",]
Patient_Post.boost.vacc <- Patient_1.Immcantation_output_annotated[Patient_1.Immcantation_output_annotated$Timepoint=="Post boost vacc",]
Patient_Post.reinfection <- Patient_1.Immcantation_output_annotated[Patient_1.Immcantation_output_annotated$Timepoint=="Post reinfection",]

# Now the four dataframes are prepared and the donut charts can be created.
Patient_Post.2nd.vacc.donut <- as.data.frame(table(Patient_Post.2nd.vacc$clone_id))
Patient_6M.post.2nd.vacc.donut <- as.data.frame(table(Patient_6M.post.2nd.vacc$clone_id))
Patient_Post.boost.vacc.donut <- as.data.frame(table(Patient_Post.boost.vacc$clone_id))
Patient_Post.reinfection.donut <- as.data.frame(table(Patient_Post.reinfection$clone_id))

Patient_Post.2nd.vacc.donut$Var1 <- as.character(Patient_Post.2nd.vacc.donut$Var1)
Patient_6M.post.2nd.vacc.donut$Var1 <- as.character(Patient_6M.post.2nd.vacc.donut$Var1)
Patient_Post.boost.vacc.donut$Var1 <- as.character(Patient_Post.boost.vacc.donut$Var1)
Patient_Post.reinfection.donut$Var1 <- as.character(Patient_Post.reinfection.donut$Var1)

colnames(Patient_Post.2nd.vacc.donut) <- c("clone_id","Freq")
colnames(Patient_6M.post.2nd.vacc.donut) <- c("clone_id","Freq")
colnames(Patient_Post.boost.vacc.donut) <- c("clone_id","Freq")
colnames(Patient_Post.reinfection.donut) <- c("clone_id","Freq")

# Now I am adding the colors
Patient_Post.2nd.vacc.donut <- merge(Patient_Post.2nd.vacc.donut,unique(Patient_1.Immcantation_output_annotated[Patient_1.Immcantation_output_annotated$Timepoint=="Post 2nd vacc",c("clone_id","color","persistence")]), by="clone_id", all.x=T)
Patient_6M.post.2nd.vacc.donut <- merge(Patient_6M.post.2nd.vacc.donut,unique(Patient_1.Immcantation_output_annotated[Patient_1.Immcantation_output_annotated$Timepoint=="6M post 2nd vacc",c("clone_id","color","persistence")]), by="clone_id", all.x=T)
Patient_Post.boost.vacc.donut <- merge(Patient_Post.boost.vacc.donut,unique(Patient_1.Immcantation_output_annotated[Patient_1.Immcantation_output_annotated$Timepoint=="Post boost vacc",c("clone_id","color","persistence")]), by="clone_id", all.x=T)
Patient_Post.reinfection.donut <- merge(Patient_Post.reinfection.donut,unique(Patient_1.Immcantation_output_annotated[Patient_1.Immcantation_output_annotated$Timepoint=="Post reinfection",c("clone_id","color","persistence")]), by="clone_id", all.x=T)

# Order by clone size
Patient_Post.2nd.vacc.donut <- Patient_Post.2nd.vacc.donut[order(Patient_Post.2nd.vacc.donut$Freq, decreasing = TRUE),]
Patient_6M.post.2nd.vacc.donut <- Patient_6M.post.2nd.vacc.donut[order(Patient_6M.post.2nd.vacc.donut$Freq, decreasing = TRUE),]
Patient_Post.boost.vacc.donut <- Patient_Post.boost.vacc.donut[order(Patient_Post.boost.vacc.donut$Freq, decreasing = TRUE),]
Patient_Post.reinfection.donut <- Patient_Post.reinfection.donut[order(Patient_Post.reinfection.donut$Freq, decreasing = TRUE),]

# Now I need to add the info, which clones should be put together in one large segment in the pie chart
# These are the clones that have the persistence "Singlet" assigned
n_singlets_2nd.vacc <- nrow(Patient_Post.2nd.vacc.donut[Patient_Post.2nd.vacc.donut$persistence=="Singlet",])
n_singlets_6M.post.2nd.vacc <- nrow(Patient_6M.post.2nd.vacc.donut[Patient_6M.post.2nd.vacc.donut$persistence=="Singlet",])
n_singlets_Post.boost <- nrow(Patient_Post.boost.vacc.donut[Patient_Post.boost.vacc.donut$persistence=="Singlet",])
n_singlets_Post.reinfection <- nrow(Patient_Post.reinfection.donut[Patient_Post.reinfection.donut$persistence=="Singlet",])

Patient_Post.2nd.vacc.donut <- Patient_Post.2nd.vacc.donut[Patient_Post.2nd.vacc.donut$persistence != "Singlet",]
Patient_6M.post.2nd.vacc.donut <- Patient_6M.post.2nd.vacc.donut[Patient_6M.post.2nd.vacc.donut$persistence != "Singlet",]
Patient_Post.boost.vacc.donut <- Patient_Post.boost.vacc.donut[Patient_Post.boost.vacc.donut$persistence != "Singlet",]
Patient_Post.reinfection.donut <- Patient_Post.reinfection.donut[Patient_Post.reinfection.donut$persistence != "Singlet",]

Patient_Post.2nd.vacc.donut[nrow(Patient_Post.2nd.vacc.donut)+1,] <- c("Singlets",n_singlets_2nd.vacc,"white","Singlet")
Patient_6M.post.2nd.vacc.donut[nrow(Patient_6M.post.2nd.vacc.donut)+1,] <- c("Singlets",n_singlets_6M.post.2nd.vacc,"white","Singlet")
Patient_Post.boost.vacc.donut[nrow(Patient_Post.boost.vacc.donut)+1,] <- c("Singlets",n_singlets_Post.boost,"white","Singlet")
Patient_Post.reinfection.donut[nrow(Patient_Post.reinfection.donut)+1,] <- c("Singlets",n_singlets_Post.reinfection,"white","Singlet")

Patient_Post.2nd.vacc.donut$Freq <- as.integer(Patient_Post.2nd.vacc.donut$Freq)
Patient_6M.post.2nd.vacc.donut$Freq <- as.integer(Patient_6M.post.2nd.vacc.donut$Freq)
Patient_Post.boost.vacc.donut$Freq <- as.integer(Patient_Post.boost.vacc.donut$Freq)
Patient_Post.reinfection.donut$Freq <- as.integer(Patient_Post.reinfection.donut$Freq)

# Order by clone size > persistence status
# First factor the persistence column with desired level order
# Then order and put the first row to the last position
Patient_Post.2nd.vacc.donut$persistence <- factor(Patient_Post.2nd.vacc.donut$persistence, levels = c("Singlet","Unique clone","Persistent clone"))
Patient_6M.post.2nd.vacc.donut$persistence <- factor(Patient_6M.post.2nd.vacc.donut$persistence, levels = c("Singlet","Unique clone","Persistent clone"))
Patient_Post.boost.vacc.donut$persistence <- factor(Patient_Post.boost.vacc.donut$persistence, levels = c("Singlet","Unique clone","Persistent clone"))
Patient_Post.reinfection.donut$persistence <- factor(Patient_Post.reinfection.donut$persistence, levels = c("Singlet","Unique clone","Persistent clone"))

Patient_Post.2nd.vacc.donut <- Patient_Post.2nd.vacc.donut[order(Patient_Post.2nd.vacc.donut$Freq,Patient_Post.2nd.vacc.donut$persistence,decreasing = T),]
Patient_6M.post.2nd.vacc.donut <- Patient_6M.post.2nd.vacc.donut[order(Patient_6M.post.2nd.vacc.donut$Freq,Patient_6M.post.2nd.vacc.donut$persistence,decreasing = T),]
Patient_Post.boost.vacc.donut <- Patient_Post.boost.vacc.donut[order(Patient_Post.boost.vacc.donut$Freq,Patient_Post.boost.vacc.donut$persistence,decreasing = T),]
Patient_Post.reinfection.donut <- Patient_Post.reinfection.donut[order(Patient_Post.reinfection.donut$Freq,Patient_Post.reinfection.donut$persistence,decreasing = T),]

Patient_Post.2nd.vacc.donut[nrow(Patient_Post.2nd.vacc.donut)+1,] <- Patient_Post.2nd.vacc.donut[1,]
Patient_6M.post.2nd.vacc.donut[nrow(Patient_6M.post.2nd.vacc.donut)+1,] <- Patient_6M.post.2nd.vacc.donut[1,]
Patient_Post.boost.vacc.donut[nrow(Patient_Post.boost.vacc.donut)+1,] <- Patient_Post.boost.vacc.donut[1,]
Patient_Post.reinfection.donut[nrow(Patient_Post.reinfection.donut)+1,] <- Patient_Post.reinfection.donut[1,]

Patient_Post.2nd.vacc.donut <- Patient_Post.2nd.vacc.donut[2:nrow(Patient_Post.2nd.vacc.donut),]
Patient_6M.post.2nd.vacc.donut <- Patient_6M.post.2nd.vacc.donut[2:nrow(Patient_6M.post.2nd.vacc.donut),]
Patient_Post.boost.vacc.donut <- Patient_Post.boost.vacc.donut[2:nrow(Patient_Post.boost.vacc.donut),]
Patient_Post.reinfection.donut <- Patient_Post.reinfection.donut[2:nrow(Patient_Post.reinfection.donut),]

# Now create the donut
# Compute percentages
Patient_Post.2nd.vacc.donut$fraction <- Patient_Post.2nd.vacc.donut$Freq / sum(Patient_Post.2nd.vacc.donut$Freq)
Patient_6M.post.2nd.vacc.donut$fraction <- Patient_6M.post.2nd.vacc.donut$Freq / sum(Patient_6M.post.2nd.vacc.donut$Freq)
Patient_Post.boost.vacc.donut$fraction <- Patient_Post.boost.vacc.donut$Freq / sum(Patient_Post.boost.vacc.donut$Freq)
Patient_Post.reinfection.donut$fraction <- Patient_Post.reinfection.donut$Freq / sum(Patient_Post.reinfection.donut$Freq)

# Compute the cumulative percentages (top of each rectangle)
Patient_Post.2nd.vacc.donut$ymax = cumsum(Patient_Post.2nd.vacc.donut$fraction)
Patient_6M.post.2nd.vacc.donut$ymax = cumsum(Patient_6M.post.2nd.vacc.donut$fraction)
Patient_Post.boost.vacc.donut$ymax = cumsum(Patient_Post.boost.vacc.donut$fraction)
Patient_Post.reinfection.donut$ymax = cumsum(Patient_Post.reinfection.donut$fraction)

# Compute the bottom of each rectangle
Patient_Post.2nd.vacc.donut$ymin = c(0, head(Patient_Post.2nd.vacc.donut$ymax, n=-1))
Patient_6M.post.2nd.vacc.donut$ymin = c(0, head(Patient_6M.post.2nd.vacc.donut$ymax, n=-1))
Patient_Post.boost.vacc.donut$ymin = c(0, head(Patient_Post.boost.vacc.donut$ymax, n=-1))
Patient_Post.reinfection.donut$ymin = c(0, head(Patient_Post.reinfection.donut$ymax, n=-1))

Patient_Post.2nd.vacc.donut$percentage <- Patient_Post.2nd.vacc.donut$fraction*100
Patient_6M.post.2nd.vacc.donut$percentage <- Patient_6M.post.2nd.vacc.donut$fraction*100
Patient_Post.boost.vacc.donut$percentage <- Patient_Post.boost.vacc.donut$fraction*100
Patient_Post.reinfection.donut$percentage <- Patient_Post.reinfection.donut$fraction*100

Patient_Post.2nd.vacc.donut$rounded <- round(Patient_Post.2nd.vacc.donut$percentage,digits = 1)
Patient_6M.post.2nd.vacc.donut$rounded <- round(Patient_6M.post.2nd.vacc.donut$percentage,digits = 1)
Patient_Post.boost.vacc.donut$rounded <- round(Patient_Post.boost.vacc.donut$percentage,digits = 1)
Patient_Post.reinfection.donut$rounded <- round(Patient_Post.reinfection.donut$percentage,digits = 1)

# Compute label position
Patient_Post.2nd.vacc.donut$labelPosition <- (Patient_Post.2nd.vacc.donut$ymax + Patient_Post.2nd.vacc.donut$ymin) / 2
Patient_6M.post.2nd.vacc.donut$labelPosition <- (Patient_6M.post.2nd.vacc.donut$ymax + Patient_6M.post.2nd.vacc.donut$ymin) / 2
Patient_Post.boost.vacc.donut$labelPosition <- (Patient_Post.boost.vacc.donut$ymax + Patient_Post.boost.vacc.donut$ymin) / 2
Patient_Post.reinfection.donut$labelPosition <- (Patient_Post.reinfection.donut$ymax + Patient_Post.reinfection.donut$ymin) / 2

Patient_Post.2nd.vacc.donut$Legend <- paste0(Patient_Post.2nd.vacc.donut$Clonotype," (",Patient_Post.2nd.vacc.donut$rounded,"%)")
Patient_6M.post.2nd.vacc.donut$Legend <- paste0(Patient_6M.post.2nd.vacc.donut$Clonotype," (",Patient_6M.post.2nd.vacc.donut$rounded,"%)")
Patient_Post.boost.vacc.donut$Legend <- paste0(Patient_Post.boost.vacc.donut$Clonotype," (",Patient_Post.boost.vacc.donut$rounded,"%)")
Patient_Post.reinfection.donut$Legend <- paste0(Patient_Post.reinfection.donut$Clonotype," (",Patient_Post.reinfection.donut$rounded,"%)")

# Calculating the outer border
outer_border_Post.2nd.vacc <- Patient_Post.2nd.vacc.donut[Patient_Post.2nd.vacc.donut$Freq>1,]
outer_border_6M.post.2nd.vacc <- Patient_6M.post.2nd.vacc.donut[Patient_6M.post.2nd.vacc.donut$Freq>1,]
outer_border_Post.boost.vacc <- Patient_Post.boost.vacc.donut[Patient_Post.boost.vacc.donut$Freq>1,]
outer_border_Post.reinfection <- Patient_Post.reinfection.donut[Patient_Post.reinfection.donut$Freq>1,]

outer_border_Post.2nd.vacc <- outer_border_Post.2nd.vacc[1:nrow(outer_border_Post.2nd.vacc)-1,]
outer_border_6M.post.2nd.vacc <- outer_border_6M.post.2nd.vacc[1:nrow(outer_border_6M.post.2nd.vacc)-1,]
outer_border_Post.boost.vacc <- outer_border_Post.boost.vacc[1:nrow(outer_border_Post.boost.vacc)-1,]
outer_border_Post.reinfection <- outer_border_Post.reinfection[1:nrow(outer_border_Post.reinfection)-1,]

# Make the plot
p <- ggplot(Patient_Post.2nd.vacc.donut, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3))+
  geom_rect(color='black', fill=Patient_Post.2nd.vacc.donut$color)+NoLegend()
if(nrow(outer_border_Post.2nd.vacc)>0){
  p <- p + geom_rect(data = outer_border_Post.2nd.vacc,  aes(xmin=4.0,xmax=4.05), fill = 'black')+
    coord_polar(theta="y")+
    xlim(c(2.3, 4.4)) +
    theme_void()+
    geom_text(x = 2.3, y = 0.5, label = sum(Patient_Post.2nd.vacc.donut$Freq), size = 10)+
    geom_text(x = 4.4, y = 0.1, label = paste(round(100*(sum(outer_border_Post.2nd.vacc$Freq)/sum(Patient_Post.2nd.vacc.donut$Freq)),1),"%"), size = 5)+
    guides(
      fill = guide_legend(
        title = "Clonotypes",
        override.aes = aes(label = "")))+ggtitle(paste(Donuttitle,"Post 2nd vacc"))+
    theme(plot.title = element_text(hjust = 0.5,vjust=-3.5,size=18),legend.position="none")
}
if(nrow(outer_border_Post.2nd.vacc)==0){
  p <-p + geom_rect(data = Patient_Post.2nd.vacc.donut,  aes(xmin=4.0,xmax=4.01), fill = 'black')+
    coord_polar(theta="y")+
    xlim(c(2.3, 4.4)) +
    theme_void()+
    geom_text(x = 2.3, y = 0.5, label = sum(Patient_Post.2nd.vacc.donut$Freq), size = 10)+
    geom_text(x = 4.4, y = 0.1, label = paste(round(100*(sum(outer_border_Post.2nd.vacc$Freq)/sum(Patient_Post.2nd.vacc.donut$Freq)),1),"%"), size = 5)+
    guides(
      fill = guide_legend(
        title = "Clonotypes",
        override.aes = aes(label = "")))+ggtitle(paste(Donuttitle,"Post 2nd vacc"))+
    theme(plot.title = element_text(hjust = 0.5,vjust=-3.5,size=18),legend.position="none")
}
p

q <- ggplot(Patient_6M.post.2nd.vacc.donut, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3))+
  geom_rect(color='black', fill=Patient_6M.post.2nd.vacc.donut$color)+NoLegend()
if(nrow(outer_border_6M.post.2nd.vacc)>0){
  q <- q + geom_rect(data = outer_border_6M.post.2nd.vacc,  aes(xmin=4.0,xmax=4.05), fill = 'black')+
    coord_polar(theta="y")+
    xlim(c(2.3, 4.4)) +
    theme_void()+
    geom_text(x = 2.3, y = 0.5, label = sum(Patient_6M.post.2nd.vacc.donut$Freq), size = 10)+
    geom_text(x = 4.4, y = 0.1, label = paste(round(100*(sum(outer_border_6M.post.2nd.vacc$Freq)/sum(Patient_6M.post.2nd.vacc.donut$Freq)),1),"%"), size = 5)+
    guides(
      fill = guide_legend(
        title = "Clonotypes",
        override.aes = aes(label = "")))+ggtitle(paste(Donuttitle,"6M post 2nd vacc"))+
    theme(plot.title = element_text(hjust = 0.5,vjust=-3.5,size=18),legend.position="none")
}
if(nrow(outer_border_6M.post.2nd.vacc)==0){
  q <-q + geom_rect(data = Patient_6M.post.2nd.vacc.donut,  aes(xmin=4.0,xmax=4.01), fill = 'black')+
    coord_polar(theta="y")+
    xlim(c(2.3, 4.4)) +
    theme_void()+
    geom_text(x = 2.3, y = 0.5, label = sum(Patient_6M.post.2nd.vacc.donut$Freq), size = 10)+
    geom_text(x = 4.4, y = 0.1, label = paste(round(100*(sum(outer_border_6M.post.2nd.vacc$Freq)/sum(Patient_6M.post.2nd.vacc.donut$Freq)),1),"%"), size = 5)+
    guides(
      fill = guide_legend(
        title = "Clonotypes",
        override.aes = aes(label = "")))+ggtitle(paste(Donuttitle,"6M post 2nd vacc"))+
    theme(plot.title = element_text(hjust = 0.5,vjust=-3.5,size=18),legend.position="none")
}
q

z <- ggplot(Patient_Post.boost.vacc.donut, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3))+
  geom_rect(color='black', fill=Patient_Post.boost.vacc.donut$color)+NoLegend()
if(nrow(outer_border_Post.boost.vacc)>0){
  z <- z + geom_rect(data = outer_border_Post.boost.vacc,  aes(xmin=4.0,xmax=4.05), fill = 'black')+
    coord_polar(theta="y")+
    xlim(c(2.3, 4.4)) +
    theme_void()+
    geom_text(x = 2.3, y = 0.5, label = sum(Patient_Post.boost.vacc.donut$Freq), size = 10)+
    geom_text(x = 4.4, y = 0.1, label = paste(round(100*(sum(outer_border_Post.boost.vacc$Freq)/sum(Patient_Post.boost.vacc.donut$Freq)),1),"%"), size = 5)+
    guides(
      fill = guide_legend(
        title = "Clonotypes",
        override.aes = aes(label = "")))+ggtitle(paste(Donuttitle,"Post boost vacc"))+
    theme(plot.title = element_text(hjust = 0.5,vjust=-3.5,size=18),legend.position="none")
}
if(nrow(outer_border_Post.boost.vacc)==0){
  z <-z + geom_rect(data = Patient_Post.boost.vacc.donut,  aes(xmin=4.0,xmax=4.01), fill = 'black')+
    coord_polar(theta="y")+
    xlim(c(2.3, 4.4)) +
    theme_void()+
    geom_text(x = 2.3, y = 0.5, label = sum(Patient_Post.boost.vacc.donut$Freq), size = 10)+
    geom_text(x = 4.4, y = 0.1, label = paste(round(100*(sum(outer_border_Post.boost.vacc$Freq)/sum(Patient_Post.boost.vacc.donut$Freq)),1),"%"), size = 5)+
    guides(
      fill = guide_legend(
        title = "Clonotypes",
        override.aes = aes(label = "")))+ggtitle(paste(Donuttitle,"Post boost vacc"))+
    theme(plot.title = element_text(hjust = 0.5,vjust=-3.5,size=18),legend.position="none")
}
z

r <- ggplot(Patient_Post.reinfection.donut, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3))+
  geom_rect(color='black', fill=Patient_Post.reinfection.donut$color)+NoLegend()
if(nrow(outer_border_Post.reinfection)>0){
  r <- r + geom_rect(data = outer_border_Post.reinfection,  aes(xmin=4.0,xmax=4.05), fill = 'black')+
    coord_polar(theta="y")+
    xlim(c(2.3, 4.4)) +
    theme_void()+
    geom_text(x = 2.3, y = 0.5, label = sum(Patient_Post.reinfection.donut$Freq), size = 10)+
    geom_text(x = 4.4, y = 0.1, label = paste(round(100*(sum(outer_border_Post.reinfection$Freq)/sum(Patient_Post.reinfection.donut$Freq)),1),"%"), size = 5)+
    guides(
      fill = guide_legend(
        title = "Clonotypes",
        override.aes = aes(label = "")))+ggtitle(paste(Donuttitle,"Post reinfection"))+
    theme(plot.title = element_text(hjust = 0.5,vjust=-3.5,size=18),legend.position="none")
}
if(nrow(outer_border_Post.reinfection)==0){
  r <-r + geom_rect(data = Patient_Post.reinfection.donut,  aes(xmin=4.0,xmax=4.01), fill = 'black')+
    coord_polar(theta="y")+
    xlim(c(2.3, 4.4)) +
    theme_void()+
    geom_text(x = 2.3, y = 0.5, label = sum(Patient_Post.reinfection.donut$Freq), size = 10)+
    geom_text(x = 4.4, y = 0.1, label = paste(round(100*(sum(outer_border_Post.reinfection$Freq)/sum(Patient_Post.reinfection.donut$Freq)),1),"%"), size = 5)+
    guides(
      fill = guide_legend(
        title = "Clonotypes",
        override.aes = aes(label = "")))+ggtitle(paste(Donuttitle,"Post reinfection"))+
    theme(plot.title = element_text(hjust = 0.5,vjust=-3.5,size=18),legend.position="none")
}
r

p+q+z+r

remove(df,outer_border_6M.post.2nd.vacc,outer_border_Post.2nd.vacc,outer_border_Post.boost.vacc,outer_border_Post.reinfection,p,q,r,z,Patient,
       Patient_1.Immcantation_output_annotated,Patient_6M.post.2nd.vacc,Patient_6M.post.2nd.vacc.donut,Patient_Post.2nd.vacc,Patient_Post.2nd.vacc.donut,
       Patient_Post.boost.vacc,Patient_Post.boost.vacc.donut,Patient_Post.reinfection,Patient_Post.reinfection.donut,qual_col_pals,timepoints,
       col_vector,Donuttitle,Longitudinals,n_singlets_2nd.vacc,n_singlets_6M.post.2nd.vacc,n_singlets_Post.boost,n_singlets_Post.reinfection,Patient.cells,
       Persistents)

remove(Patient_1.Immcantation_output,Patient_2.Immcantation_output,Patient_3.Immcantation_output)

######################################################################################################################
# Part 6: Dataset Exploration
######################################################################################################################
# Boxplot SHM across timepoints (spike wt specific cells)
MBC.subset <- subset(MBC, subset = wtSpike_specific=="yes" &
                       Timepoint !="Post reinfection")
MBC.subset <- MBC.subset@meta.data[,c("Timepoint","mu_count")]
MBC.subset <- droplevels(MBC.subset)
stat.test <- MBC.subset %>% rstatix::wilcox_test(mu_count~Timepoint)

ggplot(MBC.subset, aes(x=Timepoint,y=mu_count))+ 
  geom_boxplot(MBC.subset, mapping= aes(x = Timepoint, y = mu_count,fill=Timepoint))+
  theme_classic()+stat_pvalue_manual(stat.test, label = "p.adj", 
                                     y.position = c(42,48,54))+
  labs(x="",y="Mutational count",fill = "Subsets")+ black.axis.text()+theme(text = element_text(size = 14))+          
  scale_fill_manual(values=met.brewer(name="Egypt", n=length(unique(MBC@meta.data$Subset))))+ggtitle("Spike Wuhan-Hu-1 binding cells")+
  center.title()# Export 5x6.5

remove(MBC.subset, stat.test)







