# This script is used for the downstream clonal analysis of the dataset.
# It comes after dataset preprocessing and after Immcantation input files preparation (using the Immcantation_Input.R script) as well as processing within Immcantation.
# Input files are the preprocessed Bmem.rds object and the Immcantation output files.

# Loading packages
library(airr)
library(alakazam)
library(ggpubr)
library(ggalluvial)
library(dplyr)
library(MetBrewer)
library(ggpubr)		   
library(openxlsx)
library(VennDiagram)
library(ggplot2)				
library(circlize)
library(Seurat)
library(RColorBrewer)
rm(list = ls())
black.axis.text <- function(){theme(axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"))} # For plotting
center.title <- function(){theme(plot.title = element_text(hjust = 0.5))} # For plotting

###################################################################################################################################
# Part 1 - Automated analysis of Immcantation output
#          This part loads an output file (from the Immcantation changeo-10x pipeline) that has to be specified and then automatically
#             - determines which are longitudinal clones in that output file
#             - annotates the clones as being persistent, unique clone, persistent singlet or unique singlet
#             - uses the annotation and creates donut charts inspired by the Nussenzweig et al papers.
#             - for Fig.6.A
###################################################################################################################################
# Write a function that takes an Immcantation output and returns a list of clone_IDs which contain longitudinal clones
Bmem.Inc <- readRDS("~/NAS/Jan/Experiments and Data/210801 - Exp 035 - Profiling of COVID vaccinated patient samples/Memory B cell Fourth Experiment 220323/multi_exp035_4/Bmem.rds")
setwd("~/NAS/Jan/Experiments and Data/210801 - Exp 035 - Profiling of COVID vaccinated patient samples/Immcantation Analysis/Immcantation Outputs/update_020522")

# The following section is divided into two parts.
# Part A: For all patients which do not have longitudinal clones.
# Part B: For patients which have longitudinal clones.

# Part A:
# If a patient does not have longitudinal clones, we just take the frequencies of all expanded clones directly.
annotate.clones.with.timepoint <- function(Immcantaion.Output, Seuratobject){
  
  b <- 1
  for (b in 1:nrow(Immcantaion.Output)) {
    if(substr(Immcantaion.Output$cell_id[b],1,2)=="CC"){
      Immcantaion.Output$cell_id[b] <- paste0("Dataset_4_",substr(Immcantaion.Output$cell_id[b],3,nchar(Immcantaion.Output$cell_id[b])))
    }
    if(substr(Immcantaion.Output$cell_id[b],1,2)=="AA"){
      Immcantaion.Output$cell_id[b] <- paste0("Dataset_1_",substr(Immcantaion.Output$cell_id[b],3,nchar(Immcantaion.Output$cell_id[b])))
    }
    if(substr(Immcantaion.Output$cell_id[b],1,2)=="GG"){
      Immcantaion.Output$cell_id[b] <- paste0("Dataset_2_",substr(Immcantaion.Output$cell_id[b],3,nchar(Immcantaion.Output$cell_id[b])))
    }
    if(substr(Immcantaion.Output$cell_id[b],1,2)=="TT"){
      Immcantaion.Output$cell_id[b] <- paste0("Dataset_3_",substr(Immcantaion.Output$cell_id[b],3,nchar(Immcantaion.Output$cell_id[b])))
    }
  }
  
  Seurat.subset <- Seuratobject@meta.data[,c("Full.Row.names","Timepoint")]
  colnames(Seurat.subset) <- c("cell_id","Timepoint")
  Immcantaion.Output <- merge(Immcantaion.Output, Seurat.subset, by="cell_id", all.x=T)

}
Patient <- "PFCL1_LIM_313000" # -> Specify patient that should be analyzed here!
Patient.file <- read_airr(paste0("./",Patient,"_positives/outputs/",Patient,"_positives_0.20/",Patient,"_positives_0.20_heavy_germ-pass.tsv"))
Patient.file <- annotate.clones.with.timepoint(Patient.file, Bmem.Inc)
Patient.file <- Patient.file %>% group_by(clone_id) %>% add_count() %>% ungroup()
Patient.file$color <- "white"
Patient.file$color[Patient.file$n>1] <- "grey"
Patient.file$persistence <- "Singlet"
Patient.file$persistence[Patient.file$n>1] <- "Unique clone"
LIM828246.annotated.6mo <- Patient.file[Patient.file$Timepoint=="6 months after infection",]  # Change in object name to "LIM828246.annotated" not meaningful, just done to be able to re-use older code.
LIM828246.annotated.12mo <- Patient.file[Patient.file$Timepoint =="12 months after infection",] # Change in object name to "LIM828246.annotated" not meaningful, just done to be able to re-use older code.
LIM828246.annotated <- Patient.file
# Now continue here *** for patients without a longitudinal clone.


# Part B: 
# If a patient has longitudinal (persistent) clones, calculations of clonal expansions are a bit complicated.
# This is handled in the code further below.
# Annotate clones in the Immcantation output file.
# We want to make donut charts with a similar style as Nussenzweig uses them.
# Individual pie chart pieces for all expanded clones that are persistent -> in color
# Here, make sure that all cells belonging to persistent clones get the same color (the same clone should have the same color in both timepoints)
# Individual pie chart pieces for all expanded clones that are unique -> grey
# Individual pie chart pieces for all singlets that are persistent (repeating sequences isolated only once per time point)
Patient <- "PFCL1_LIM_313000" # -> Specify patient that should be analyzed here!
LIM828246 <- read_airr(paste0("./",Patient,"_positives/outputs/",Patient,"_positives_0.20/",Patient,"_positives_0.20_heavy_germ-pass.tsv"))
find.longitudinal.clones <- function(Immcantation.Output){
  IDs <- unique(Immcantation.Output$clone_id)
  output.vector <- c("")
  i <- 1
  for(i in 1:length(IDs)){
    
    candidates <- Immcantation.Output$cell_id[Immcantation.Output$clone_id ==IDs[i]]
    b <- 1
    for (b in 1:length(candidates)) {
      if(substr(candidates[b],1,2)=="CC"){
        candidates[b] <- paste0("Dataset_4_",substr(candidates[b],3,nchar(candidates[b])))
      }
      if(substr(candidates[b],1,2)=="AA"){
        candidates[b] <- paste0("Dataset_1_",substr(candidates[b],3,nchar(candidates[b])))
      }
      if(substr(candidates[b],1,2)=="GG"){
        candidates[b] <- paste0("Dataset_2_",substr(candidates[b],3,nchar(candidates[b])))
      }
      if(substr(candidates[b],1,2)=="TT"){
        candidates[b] <- paste0("Dataset_3_",substr(candidates[b],3,nchar(candidates[b])))
      }
    }
    
    timepoints <-Bmem.Inc@meta.data$Timepoint[rownames(Bmem.Inc@meta.data)%in%candidates]
    if(length(unique(timepoints))>1){
      output.vector <- c(output.vector,IDs[i])
    }
  }
  
  output.vector <- output.vector[2:length(output.vector)]
  return(output.vector)
  
}
LIM828246.longitudinals <- find.longitudinal.clones(LIM828246)

annotate.clones <- function(Longitudinals,Immcantation.Output,Seuratobject){
  IDs <- Longitudinals
  output.df <- data.frame(matrix(ncol = 2, nrow = 0))
  x <- c("Clone", "TruePersistent")
  colnames(output.df) <- x
  i <- 1
  for(i in 1:length(IDs)){
    
    candidates <- Immcantation.Output$cell_id[Immcantation.Output$clone_id ==IDs[i]]
    b <- 1
    for (b in 1:length(candidates)) {
      if(substr(candidates[b],1,2)=="CC"){
        candidates[b] <- paste0("Dataset_4_",substr(candidates[b],3,nchar(candidates[b])))
      }
      if(substr(candidates[b],1,2)=="AA"){
        candidates[b] <- paste0("Dataset_1_",substr(candidates[b],3,nchar(candidates[b])))
      }
      if(substr(candidates[b],1,2)=="GG"){
        candidates[b] <- paste0("Dataset_2_",substr(candidates[b],3,nchar(candidates[b])))
      }
      if(substr(candidates[b],1,2)=="TT"){
        candidates[b] <- paste0("Dataset_3_",substr(candidates[b],3,nchar(candidates[b])))
      }
    }
    
    Bmem.Inc <- Seuratobject
    candidates.metadata <-Bmem.Inc@meta.data[rownames(Bmem.Inc@meta.data)%in%candidates,c(1,18,19)]
    if(nrow(candidates.metadata)>2){
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
                                    Immcantation.Output$clone_id %notin% output.df$Clone] <- "Persistent clone"
  Immcantation.Output$persistence[Immcantation.Output$clone_id %in% Longitudinals &
                                    Immcantation.Output$clone_id %in% output.df$Clone] <- "Persistent clone"
  return(Immcantation.Output)
  
}
LIM828246.annotated <- annotate.clones(LIM828246.longitudinals, LIM828246, Bmem.Inc)

# Now we make the donut charts, one for each timepoint and patient (two per patient).
# We will mainly use the Immcantation output file here but the file needs to be updated with a column that contains the time point
LIM828246.cells <- LIM828246$cell_id
i <- 1
for (i in 1:length(LIM828246.cells)) {
  if(substr(LIM828246.cells[i],1,2)=="CC"){
    LIM828246.cells[i] <- paste0("Dataset_4_",substr(LIM828246.cells[i],3,nchar(LIM828246.cells[i])))
  }
  if(substr(LIM828246.cells[i],1,2)=="AA"){
    LIM828246.cells[i] <- paste0("Dataset_1_",substr(LIM828246.cells[i],3,nchar(LIM828246.cells[i])))
  }
  if(substr(LIM828246.cells[i],1,2)=="GG"){
    LIM828246.cells[i] <- paste0("Dataset_2_",substr(LIM828246.cells[i],3,nchar(LIM828246.cells[i])))
  }
  if(substr(LIM828246.cells[i],1,2)=="TT"){
    LIM828246.cells[i] <- paste0("Dataset_3_",substr(LIM828246.cells[i],3,nchar(LIM828246.cells[i])))
  }
}
timepoints <- Bmem.Inc@meta.data[rownames(Bmem.Inc@meta.data)%in%LIM828246.cells,c("Timepoint","Full.Row.names")]

i <- 1
for (i in 1:nrow(timepoints)) {
  if(substr(rownames(timepoints)[i],1,10)=="Dataset_4_"){
    rownames(timepoints)[i] <- paste0("CC",substr(rownames(timepoints)[i],11,nchar(rownames(timepoints)[i])))
  }
  if(substr(rownames(timepoints)[i],1,10)=="Dataset_1_"){
    rownames(timepoints)[i] <- paste0("AA",substr(rownames(timepoints)[i],11,nchar(rownames(timepoints)[i])))
  }
  if(substr(rownames(timepoints)[i],1,10)=="Dataset_2_"){
    rownames(timepoints)[i] <- paste0("GG",substr(rownames(timepoints)[i],11,nchar(rownames(timepoints)[i])))
  }
  if(substr(rownames(timepoints)[i],1,10)=="Dataset_3_"){
    rownames(timepoints)[i] <- paste0("TT",substr(rownames(timepoints)[i],11,nchar(rownames(timepoints)[i])))
  }
}

LIM828246.annotated <- as.data.frame(LIM828246.annotated)
rownames(LIM828246.annotated) <- LIM828246.annotated$cell_id
LIM828246.annotated <- merge(LIM828246.annotated, timepoints, by=0, all.x = T)


# Adding the colors that the pie chart segments will have later. Options are Singlet, Unique clone, Persistent clone and Persistent singlet
LIM828246.annotated$color <- NA
LIM828246.annotated$color[LIM828246.annotated$persistence=="Singlet"] <- "white"
LIM828246.annotated$color[LIM828246.annotated$persistence=="Unique clone"] <- "grey"
LIM828246.annotated$color[LIM828246.annotated$persistence=="Persistent singlet"] <- "white"

# For the category Persistent clone we need to give colors.
Persistents <- unique(LIM828246.annotated$clone_id[LIM828246.annotated$persistence=="Persistent clone"])

qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col_vector <- col_vector[1:length(Persistents)]
df <- data.frame("clone_id"= Persistents, "color"=col_vector)

LIM828246.annotated <- merge(LIM828246.annotated, df, by="clone_id", all.x = T)
LIM828246.annotated$color <- coalesce(LIM828246.annotated$color.x,LIM828246.annotated$color.y)
LIM828246.annotated$color.x <- NULL
LIM828246.annotated$color.y <- NULL

LIM828246.annotated.6mo <- LIM828246.annotated[LIM828246.annotated$Timepoint=="6 months after infection",]
LIM828246.annotated.12mo <- LIM828246.annotated[LIM828246.annotated$Timepoint =="12 months after infection",]


# End of the separation between patients with or without longitudinal clones. Both parts (A and B) continue here.
# Here continuing no matter if patient has persistent clones or not! (Continue after ***)
# Now the two dataframes are prepared and the donut charts can be created.
LIM828246.annotated.6mo.donut <- as.data.frame(table(LIM828246.annotated.6mo$clone_id))
LIM828246.annotated.12mo.donut <- as.data.frame(table(LIM828246.annotated.12mo$clone_id))
LIM828246.annotated.6mo.donut$Var1 <- as.character(LIM828246.annotated.6mo.donut$Var1)
LIM828246.annotated.12mo.donut$Var1 <- as.character(LIM828246.annotated.12mo.donut$Var1)
colnames(LIM828246.annotated.6mo.donut) <- c("clone_id","Freq")
colnames(LIM828246.annotated.12mo.donut) <- c("clone_id","Freq")

# Now I am adding the colors
LIM828246.annotated.6mo.donut <- merge(LIM828246.annotated.6mo.donut,unique(LIM828246.annotated[LIM828246.annotated$Timepoint=="6 months after infection",c("clone_id","color","persistence")]), by="clone_id", all.x=T)
LIM828246.annotated.12mo.donut <- merge(LIM828246.annotated.12mo.donut,unique(LIM828246.annotated[LIM828246.annotated$Timepoint=="12 months after infection",c("clone_id","color","persistence")]), by="clone_id", all.x=T)

# Order by clone size
LIM828246.annotated.6mo.donut <- LIM828246.annotated.6mo.donut[order(LIM828246.annotated.6mo.donut$Freq, decreasing = TRUE),]
LIM828246.annotated.12mo.donut <- LIM828246.annotated.12mo.donut[order(LIM828246.annotated.12mo.donut$Freq, decreasing = TRUE),]

# Now I need to add the info, which clones should be put together in one large segment in the pie chart
# These are the clones that have the persistence "Singlet" assigned
n_singlets_6mo <- nrow(LIM828246.annotated.6mo.donut[LIM828246.annotated.6mo.donut$persistence=="Singlet",])
n_singlets_12mo <- nrow(LIM828246.annotated.12mo.donut[LIM828246.annotated.12mo.donut$persistence=="Singlet",])

LIM828246.annotated.6mo.donut <- LIM828246.annotated.6mo.donut[LIM828246.annotated.6mo.donut$persistence != "Singlet",]
LIM828246.annotated.12mo.donut <- LIM828246.annotated.12mo.donut[LIM828246.annotated.12mo.donut$persistence != "Singlet",]
LIM828246.annotated.6mo.donut[nrow(LIM828246.annotated.6mo.donut)+1,] <- c("Singlets",n_singlets_6mo,"white","Singlet")
LIM828246.annotated.12mo.donut[nrow(LIM828246.annotated.12mo.donut)+1,] <- c("Singlets",n_singlets_12mo,"white","Singlet")

LIM828246.annotated.6mo.donut$Freq <- as.integer(LIM828246.annotated.6mo.donut$Freq)
LIM828246.annotated.12mo.donut$Freq <- as.integer(LIM828246.annotated.12mo.donut$Freq)

# Order by clone size > persistence status
# First factor the persistence column with desired level order
LIM828246.annotated.6mo.donut$persistence <- factor(LIM828246.annotated.6mo.donut$persistence, 
                                                    levels = c("Singlet","Unique clone","Persistent clone"))
LIM828246.annotated.12mo.donut$persistence <- factor(LIM828246.annotated.12mo.donut$persistence, 
                                                     levels = c("Singlet","Unique clone","Persistent clone"))

# Then order and put the first row to the last position
LIM828246.annotated.6mo.donut <- LIM828246.annotated.6mo.donut[order(LIM828246.annotated.6mo.donut$Freq,
                                                                     LIM828246.annotated.6mo.donut$persistence,decreasing = T),]
LIM828246.annotated.12mo.donut <- LIM828246.annotated.12mo.donut[order(LIM828246.annotated.12mo.donut$Freq,
                                                                       LIM828246.annotated.12mo.donut$persistence,decreasing = T),]
LIM828246.annotated.6mo.donut[nrow(LIM828246.annotated.6mo.donut)+1,] <- LIM828246.annotated.6mo.donut[1,]
LIM828246.annotated.12mo.donut[nrow(LIM828246.annotated.12mo.donut)+1,] <- LIM828246.annotated.12mo.donut[1,]
LIM828246.annotated.6mo.donut <- LIM828246.annotated.6mo.donut[2:nrow(LIM828246.annotated.6mo.donut),]
LIM828246.annotated.12mo.donut <- LIM828246.annotated.12mo.donut[2:nrow(LIM828246.annotated.12mo.donut),]


# Now create the donut
# Compute percentages
LIM828246.annotated.6mo.donut$fraction <- LIM828246.annotated.6mo.donut$Freq / sum(LIM828246.annotated.6mo.donut$Freq)
LIM828246.annotated.12mo.donut$fraction <- LIM828246.annotated.12mo.donut$Freq / sum(LIM828246.annotated.12mo.donut$Freq)

# Compute the cumulative percentages (top of each rectangle)
LIM828246.annotated.6mo.donut$ymax = cumsum(LIM828246.annotated.6mo.donut$fraction)
LIM828246.annotated.12mo.donut$ymax = cumsum(LIM828246.annotated.12mo.donut$fraction)

# Compute the bottom of each rectangle
LIM828246.annotated.6mo.donut$ymin = c(0, head(LIM828246.annotated.6mo.donut$ymax, n=-1))
LIM828246.annotated.12mo.donut$ymin = c(0, head(LIM828246.annotated.12mo.donut$ymax, n=-1))

LIM828246.annotated.6mo.donut$percentage <- LIM828246.annotated.6mo.donut$fraction*100
LIM828246.annotated.12mo.donut$percentage <- LIM828246.annotated.12mo.donut$fraction*100

LIM828246.annotated.6mo.donut$rounded <- round(LIM828246.annotated.6mo.donut$percentage,digits = 1)
LIM828246.annotated.12mo.donut$rounded <- round(LIM828246.annotated.12mo.donut$percentage,digits = 1)

# Compute label position
LIM828246.annotated.6mo.donut$labelPosition <- (LIM828246.annotated.6mo.donut$ymax + LIM828246.annotated.6mo.donut$ymin) / 2
LIM828246.annotated.12mo.donut$labelPosition <- (LIM828246.annotated.12mo.donut$ymax + LIM828246.annotated.12mo.donut$ymin) / 2

LIM828246.annotated.6mo.donut$Legend <- paste0(LIM828246.annotated.6mo.donut$Clonotype," (",LIM828246.annotated.6mo.donut$rounded,"%)")
LIM828246.annotated.12mo.donut$Legend <- paste0(LIM828246.annotated.12mo.donut$Clonotype," (",LIM828246.annotated.12mo.donut$rounded,"%)")

outer_border_6mo.test <- LIM828246.annotated.6mo.donut[1:(nrow(LIM828246.annotated.6mo.donut)-1),] # Wrong for true outer border, because persistent singlets are then counted as clonal
outer_border_12mo.test <- LIM828246.annotated.12mo.donut[1:(nrow(LIM828246.annotated.12mo.donut)-1),] # Wrong for true outer border, because persistent singlets are then counted as clonal
outer_border_6mo <- LIM828246.annotated.6mo.donut[LIM828246.annotated.6mo.donut$Freq>1,]
outer_border_12mo <- LIM828246.annotated.12mo.donut[LIM828246.annotated.12mo.donut$Freq>1,]
outer_border_6mo <- outer_border_6mo[1:nrow(outer_border_6mo)-1,]
outer_border_12mo <- outer_border_12mo[1:nrow(outer_border_12mo)-1,]


# Make the plot
# This happens in several steps. The first step is the same for all cases
p <- ggplot(LIM828246.annotated.6mo.donut, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3))+
  geom_rect(color='black', fill=LIM828246.annotated.6mo.donut$color)+NoLegend()

# Then, it matters whether there are expanded clones or not. If there are not expanded clones, then there will be no outer border.
if(nrow(outer_border_6mo)>0){
p <- p + geom_rect(data = outer_border_6mo,  aes(xmin=4.0,xmax=4.05), fill = 'black')+
  coord_polar(theta="y")+
  xlim(c(2.3, 4.4)) +
  theme_void()+
  geom_text(x = 2.3, y = 0.5, label = sum(LIM828246.annotated.6mo.donut$Freq), size = 10)+
  geom_text(x = 4.4, y = 0.1, label = paste(round(100*(sum(outer_border_6mo$Freq)/sum(LIM828246.annotated.6mo.donut$Freq)),1),"%"), size = 5)+
  guides(
    fill = guide_legend(
      title = "Clonotypes",
      override.aes = aes(label = "")))+ggtitle(paste(Patient,"6 months"))+
  theme(plot.title = element_text(hjust = 0.5,vjust=-3.5,size=18),legend.position="none")
}
if(nrow(outer_border_6mo)==0){
  p <-p + geom_rect(data = LIM828246.annotated.6mo.donut,  aes(xmin=4.0,xmax=4.01), fill = 'black')+
    coord_polar(theta="y")+
    xlim(c(2.3, 4.4)) +
    theme_void()+
    geom_text(x = 2.3, y = 0.5, label = sum(LIM828246.annotated.6mo.donut$Freq), size = 10)+
    geom_text(x = 4.4, y = 0.1, label = paste(round(100*(sum(outer_border_6mo$Freq)/sum(LIM828246.annotated.6mo.donut$Freq)),1),"%"), size = 5)+
    guides(
      fill = guide_legend(
        title = "Clonotypes",
        override.aes = aes(label = "")))+ggtitle(paste(Patient,"6 months"))+
    theme(plot.title = element_text(hjust = 0.5,vjust=-3.5,size=18),legend.position="none")
}
  
p
# The same that was done for the 6 month timepoint is now done for the 12 month timepoint.
q <- ggplot(LIM828246.annotated.12mo.donut, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3)) +
  geom_rect(color='black', fill=LIM828246.annotated.12mo.donut$color)+NoLegend()

if(nrow(outer_border_12mo)>0){
q <- q + geom_rect(data = outer_border_12mo,  aes(xmin=4.0,xmax=4.05), fill = 'black')+
  coord_polar(theta="y")+
  xlim(c(2.3, 4.4)) +
  theme_void()+
  geom_text(x = 2.3, y = 0.5, label = sum(LIM828246.annotated.12mo.donut$Freq), size = 10)+
  geom_text(x = 4.4, y = 0.1, label = paste(round(100*(sum(outer_border_12mo$Freq)/sum(LIM828246.annotated.12mo.donut$Freq)),1),"%"), size = 5)+
  guides(
    fill = guide_legend(
      title = "Clonotypes",
      override.aes = aes(label = "")))+ggtitle(paste(Patient,"12 months"))+
  theme(plot.title = element_text(hjust = 0.5,vjust=-3.5,size=18),legend.position="none")
}
if(nrow(outer_border_6mo)==0){
  q <-q + geom_rect(data = LIM828246.annotated.12mo.donut,  aes(xmin=4.0,xmax=4.01), fill = 'black')+
    coord_polar(theta="y")+
    xlim(c(2.3, 4.4)) +
    theme_void()+
    geom_text(x = 2.3, y = 0.5, label = sum(LIM828246.annotated.12mo.donut$Freq), size = 10)+
    geom_text(x = 4.4, y = 0.1, label = paste(round(100*(sum(outer_border_12mo$Freq)/sum(LIM828246.annotated.12mo.donut$Freq)),1),"%"), size = 5)+
    guides(
      fill = guide_legend(
        title = "Clonotypes",
        override.aes = aes(label = "")))+ggtitle(paste(Patient,"12 months"))+
    theme(plot.title = element_text(hjust = 0.5,vjust=-3.5,size=18),legend.position="none")
}
q

# Finally, the plots are shown together, above each other.
ggarrange(p,q, ncol = 1) # Export 6.2x4.5

remove(df, LIM828246,LIM828246.annotated,LIM828246.annotated.12mo,LIM828246.annotated.6mo,LIM828246.annotated.12mo.donut,LIM828246.annotated.6mo.donut,
        LIM828246.cells,LIM828246.longitudinals,outer_border_12mo,outer_border_6mo,p,q,qual_col_pals,timepoints,col_vector,n_singlets_12mo,n_singlets_6mo,
       Persistents,Patient, Bmem.Inc, outer_border_12mo.test,outer_border_6mo.test,six_mo_border,twlve_mo_border,i)


###################################################################################################################################
# Part 2 - Automated analysis to decide a cutoff value for the distance
#          This part is using a part of the code from Part 2 but puts it into a bigger loop that allows for automated cutoff
#          testing across all patients.
###################################################################################################################################
Patientvector <- c("PFCL1_154583_1943","PFCL1_179308_1988","PFCL1_196878_1691","PFCL1_199474_1994","PFCL1_LIM_137402",
                   "PFCL1_LIM_313000","PFCL1_LIM_674950","PFCL1_LIM828246","PFCL1_UST_190762")
Bmem.Inc <- readRDS("~/NAS/Jan/Experiments and Data/210801 - Exp 035 - Profiling of COVID vaccinated patient samples/Memory B cell Fourth Experiment 220323/multi_exp035_4/Bmem.rds")
outputfile <- data.frame(Patient=as.character(),Cutoff=as.numeric(),n_longitudinals=as.numeric(),n_persistents=as.numeric(),
                         n_unique_clones=as.numeric(),n_persistent_clones=as.numeric(),n_singlets = as.numeric())
i <- 1
for (i in 1:length(Patientvector)) {
  thisroundspatient <- Patientvector[i]
  thisroundsoutputfile <- data.frame(matrix(ncol = 7, nrow = 9))
  colnames(thisroundsoutputfile) <- c("Patient","Cutoff","n_longitudinals","n_persistens",
                                     "n_unique_clones","n_persistent_singlets","n_singlets")
  setwd(paste0("~/NAS/Jan/Experiments and Data/210801 - Exp 035 - Profiling of COVID vaccinated patient samples/Immcantation Analysis/Immcantation Outputs/update_020522/",thisroundspatient,"_positives"))
  
  cutoffs <- c("0.15","0.16","0.17","0.18","0.20","0.22","0.23","0.25","0.27")
  
  o <- 1
  for (o in 1:length(cutoffs)) {
    thisroundscutoff <- cutoffs[o]
    
    thisroundspatientfile <- read_airr(paste0("./outputs/",thisroundspatient,"_positives_",thisroundscutoff,"/",thisroundspatient,"_positives_",thisroundscutoff,"_heavy_germ-pass.tsv"))

  find.longitudinal.clones <- function(Immcantation.Output){
  IDs <- unique(Immcantation.Output$clone_id)
  output.vector <- c("")
  i <- 1
  for(i in 1:length(IDs)){
    
    candidates <- Immcantation.Output$cell_id[Immcantation.Output$clone_id ==IDs[i]]
    b <- 1
    for (b in 1:length(candidates)) {
      if(substr(candidates[b],1,2)=="CC"){
        candidates[b] <- paste0("Dataset_4_",substr(candidates[b],3,nchar(candidates[b])))
      }
      if(substr(candidates[b],1,2)=="AA"){
        candidates[b] <- paste0("Dataset_1_",substr(candidates[b],3,nchar(candidates[b])))
      }
      if(substr(candidates[b],1,2)=="GG"){
        candidates[b] <- paste0("Dataset_2_",substr(candidates[b],3,nchar(candidates[b])))
      }
      if(substr(candidates[b],1,2)=="TT"){
        candidates[b] <- paste0("Dataset_3_",substr(candidates[b],3,nchar(candidates[b])))
      }
    }
    
    timepoints <-Bmem.Inc@meta.data$Timepoint[rownames(Bmem.Inc@meta.data)%in%candidates]
    if(length(unique(timepoints))>1){
      output.vector <- c(output.vector,IDs[i])
    }
  }
  
  output.vector <- output.vector[2:length(output.vector)]
  return(output.vector)
  
}
    
    Longitudinals <- find.longitudinal.clones(thisroundspatientfile)
    
    
    
      
    annotate.clones <- function(Longitudinals,Immcantation.Output){
        IDs <- Longitudinals
        if(!NA %in% IDs) {
          output.df <- data.frame(matrix(ncol = 2, nrow = 0))
        x <- c("Clone", "TruePersistent")
        colnames(output.df) <- x
        i <- 1
        for(i in 1:length(IDs)){
          
          candidates <- Immcantation.Output$cell_id[Immcantation.Output$clone_id ==IDs[i]]
          b <- 1
          for (b in 1:length(candidates)) {
            if(substr(candidates[b],1,2)=="CC"){
              candidates[b] <- paste0("Dataset_4_",substr(candidates[b],3,nchar(candidates[b])))
            }
            if(substr(candidates[b],1,2)=="AA"){
              candidates[b] <- paste0("Dataset_1_",substr(candidates[b],3,nchar(candidates[b])))
            }
            if(substr(candidates[b],1,2)=="GG"){
              candidates[b] <- paste0("Dataset_2_",substr(candidates[b],3,nchar(candidates[b])))
            }
            if(substr(candidates[b],1,2)=="TT"){
              candidates[b] <- paste0("Dataset_3_",substr(candidates[b],3,nchar(candidates[b])))
            }
          }
          
          
          candidates.metadata <-Bmem.Inc@meta.data[rownames(Bmem.Inc@meta.data)%in%candidates,c(1,18,19)]
          if(nrow(candidates.metadata)>2){
            output.df[i,1] <- IDs[i]
            output.df[i,2] <- "yes"
          } else {
            output.df[i,1] <- IDs[i]
            output.df[i,2] <- "no"
          }
        }
        output.df <- output.df[output.df$TruePersistent=="yes",]
        Immcantation.Output <- Immcantation.Output %>% group_by(clone_id) %>% add_count() %>% ungroup()
        Immcantation.Output$persistence <- "Unique singlet"
        `%notin%` <- Negate(`%in%`)
        Immcantation.Output$persistence[Immcantation.Output$clone_id %notin% Longitudinals&
                                          Immcantation.Output$n >1] <- "Unique clone"
        Immcantation.Output$persistence[Immcantation.Output$clone_id %in% Longitudinals &
                                          Immcantation.Output$clone_id %notin% output.df$Clone] <- "Persistent singlet"
        Immcantation.Output$persistence[Immcantation.Output$clone_id %in% Longitudinals &
                                          Immcantation.Output$clone_id %in% output.df$Clone] <- "Persistent clone"
        return(Immcantation.Output)
        
      }
      
        else{
          Immcantation.Output <- Immcantation.Output %>% group_by(clone_id) %>% add_count() %>% ungroup()
          Immcantation.Output$persistence <- "Unique singlet"
          Immcantation.Output$persistence[Immcantation.Output$n >1] <- "Unique clone"
          return(Immcantation.Output)
        }
      }
      
        
        
      Annotated <- annotate.clones(Longitudinals, thisroundspatientfile)
      
      # subsetting Annotated file to one line per clone
      Annotated <- Annotated[,c("clone_id","n","persistence")]
      Annotated <- unique(Annotated)
      
      thisroundsoutputfile$Patient[o] <- thisroundspatient
      thisroundsoutputfile$Cutoff[o] <- thisroundscutoff
      thisroundsoutputfile$n_longitudinals[o] <- length(Longitudinals[Longitudinals!=""& !is.na(Longitudinals)])
      thisroundsoutputfile$n_persistens[o] <- nrow(Annotated[Annotated$persistence=="Persistent clone",])
      thisroundsoutputfile$n_unique_clones[o] <- nrow(Annotated[Annotated$persistence=="Unique clone",])
      thisroundsoutputfile$n_persistent_singlets[o] <- nrow(Annotated[Annotated$persistence=="Persistent singlet",])
      thisroundsoutputfile$n_singlets[o] <- nrow(Annotated[Annotated$persistence=="Unique singlet",])
      
      
      }
        
      
      outputfile <- rbind(outputfile,thisroundsoutputfile)
    
    
    
    
  }
  
outputfile <- outputfile %>% group_by(Cutoff) %>% add_tally(n_longitudinals, name = "Longitudinals per Cutoff")

# I decide for now to go with a distance cutoff of 0.2 
#write.csv2(outputfile, "~/NAS/Jan/Experiments and Data/210801 - Exp 035 - Profiling of COVID vaccinated patient samples/Memory B cell Fourth Experiment 220323/multi_exp035_4/Distance cutoff effects.csv")

remove(Annotated, Bmem.Inc,outputfile,thisroundsoutputfile,thisroundspatientfile,cutoffs,i,Longitudinals,o,Patientvector,thisroundscutoff,thisroundspatient)

###################################################################################################################################
# Part 3 - Getting an overview over the persistent clones
#          After Part 3 I have decided to pick the distance cutoff 0.2 for the downstream analysis
#          To which phenotypes do longitudinal clones belong? - One file for all patients
###################################################################################################################################
Bmem.Inc <- readRDS("~/NAS/Jan/Experiments and Data/210801 - Exp 035 - Profiling of COVID vaccinated patient samples/Memory B cell Fourth Experiment 220323/multi_exp035_4/Bmem.rds")
setwd("~/NAS/Jan/Experiments and Data/210801 - Exp 035 - Profiling of COVID vaccinated patient samples/Immcantation Analysis/Immcantation Outputs/update_020522/")

Patientvector <- c("PFCL1_154583_1943","PFCL1_179308_1988","PFCL1_196878_1691","PFCL1_199474_1994","PFCL1_LIM_137402",
                   "PFCL1_LIM_313000","PFCL1_LIM_674950","PFCL1_LIM828246","PFCL1_UST_190762")
outputfile <- data.frame(cell_id=as.character(),clone_id=as.character(),Tmiepoint=as.character(),Patient=as.character(),
                         Vacinated=as.character(),Isotype=as.character(),Phenotype = as.character(),mu_count = as.numeric())
i <- 1
for (i in 1:length(Patientvector)) {
  
  thisroundspatient <- Patientvector[i]
  patientfile <- read_airr(paste0("./",thisroundspatient,"_positives/outputs/",thisroundspatient,"_positives_0.20/",thisroundspatient,"_positives_0.20_heavy_germ-pass.tsv"))
  find.longitudinal.clones <- function(Immcantation.Output){
  IDs <- unique(Immcantation.Output$clone_id)
  output.vector <- c("")
  i <- 1
  for(i in 1:length(IDs)){
    
    candidates <- Immcantation.Output$cell_id[Immcantation.Output$clone_id ==IDs[i]]
    b <- 1
    for (b in 1:length(candidates)) {
      if(substr(candidates[b],1,2)=="CC"){
        candidates[b] <- paste0("Dataset_4_",substr(candidates[b],3,nchar(candidates[b])))
      }
      if(substr(candidates[b],1,2)=="AA"){
        candidates[b] <- paste0("Dataset_1_",substr(candidates[b],3,nchar(candidates[b])))
      }
      if(substr(candidates[b],1,2)=="GG"){
        candidates[b] <- paste0("Dataset_2_",substr(candidates[b],3,nchar(candidates[b])))
      }
      if(substr(candidates[b],1,2)=="TT"){
        candidates[b] <- paste0("Dataset_3_",substr(candidates[b],3,nchar(candidates[b])))
      }
    }
    
    timepoints <-Bmem.Inc@meta.data$Timepoint[rownames(Bmem.Inc@meta.data)%in%candidates]
    if(length(unique(timepoints))>1){
      output.vector <- c(output.vector,IDs[i])
    }
  }
  
  output.vector <- output.vector[2:length(output.vector)]
  return(output.vector)
  
} # This function can handle all patients, even those that don't have longitudinal clones.
  rounds.longitudinals <- find.longitudinal.clones(patientfile)
  
  if(NA %in% rounds.longitudinals){
    next
  }
  
  else{
    
    Longitudinal.cells.imc <- patientfile[patientfile$clone_id %in% rounds.longitudinals,c(49,52)]
    Longitudinal.cells.imc <- as.data.frame(Longitudinal.cells.imc)
    
    b <- 1
    for (b in 1:nrow(Longitudinal.cells.imc)) {
    if(substr(Longitudinal.cells.imc[b,1],1,2)=="CC"){
      Longitudinal.cells.imc[b,1] <- paste0("Dataset_4_",substr(Longitudinal.cells.imc[b,1],3,nchar(Longitudinal.cells.imc[b,1])))
    }
    if(substr(Longitudinal.cells.imc[b,1],1,2)=="AA"){
      Longitudinal.cells.imc[b,1] <- paste0("Dataset_1_",substr(Longitudinal.cells.imc[b,1],3,nchar(Longitudinal.cells.imc[b,1])))
    }
    if(substr(Longitudinal.cells.imc[b,1],1,2)=="GG"){
      Longitudinal.cells.imc[b,1] <- paste0("Dataset_2_",substr(Longitudinal.cells.imc[b,1],3,nchar(Longitudinal.cells.imc[b,1])))
    }
    if(substr(Longitudinal.cells.imc[b,1],1,2)=="TT"){
      Longitudinal.cells.imc[b,1] <- paste0("Dataset_3_",substr(Longitudinal.cells.imc[b,1],3,nchar(Longitudinal.cells.imc[b,1])))
    }
  }
    
    Longitudinal.cells.seurat <- Bmem.Inc@meta.data[Bmem.Inc@meta.data$Full.Row.names %in% Longitudinal.cells.imc$cell_id,]
    Longitudinal.cells.seurat <- Longitudinal.cells.seurat[,c(18,19,20,32,47,53)]
    rownames(Longitudinal.cells.imc) <- Longitudinal.cells.imc$cell_id
    Longitudinals <- merge(Longitudinal.cells.imc,Longitudinal.cells.seurat, by=0)
    Longitudinals <- Longitudinals[2:9]
    
  }
  
  outputfile <- rbind(outputfile,Longitudinals)
  
}
write.csv2(outputfile, "~/NAS/Jan/Experiments and Data/210801 - Exp 035 - Profiling of COVID vaccinated patient samples/Memory B cell Fourth Experiment 220323/multi_exp035_4/Phenotypes of cells belonging to longitudinal clones.csv")# Created in part 3 of this script

remove(Bmem.Inc, Longitudinal.cells.imc, Longitudinal.cells.seurat, Longitudinals,patientfile,Patientvector,b,i,
       rounds.longitudinals,thisroundspatient,outputfile)

###################################################################################################################################
# Part 4 - Checking the IGG sub-Isotypes from bait specific cells - also depending on the cell subset. Supp.Fig.5.E
###################################################################################################################################
setwd("~/NAS/Jan/Experiments and Data/210801 - Exp 035 - Profiling of COVID vaccinated patient samples/Memory B cell Fourth Experiment 220323/multi_exp035_4")  
Bmem <- readRDS("Bmem.rds")

setwd("~/NAS/Jan/Experiments and Data/210801 - Exp 035 - Profiling of COVID vaccinated patient samples/Immcantation Analysis/Immcantation Outputs/update_020522/all_cells")
all.cells <- as.data.frame(read_airr("all_cells_heavy_clone-pass_no-duplicate.tsv"))

Bmem.Inc <- Bmem
b <- 1
for (b in 1:nrow(all.cells)) {
  if(substr(all.cells[b,49],1,2)=="CC"){
    all.cells[b,49] <- paste0("Dataset_4_",substr(all.cells[b,49],3,nchar(all.cells[b,49])))
  }
  if(substr(all.cells[b,49],1,2)=="AA"){
    all.cells[b,49] <- paste0("Dataset_1_",substr(all.cells[b,49],3,nchar(all.cells[b,49])))
  }
  if(substr(all.cells[b,49],1,2)=="GG"){
    all.cells[b,49] <- paste0("Dataset_2_",substr(all.cells[b,49],3,nchar(all.cells[b,49])))
  }
  if(substr(all.cells[b,49],1,2)=="TT"){
    all.cells[b,49] <- paste0("Dataset_3_",substr(all.cells[b,49],3,nchar(all.cells[b,49])))
  }
}

ignore.positivity <- all.cells
ignore.positivity <- ignore.positivity[,c(49,50)]

ignore.positivity <- ignore.positivity %>% group_by(c_call) %>% add_tally()%>% ungroup()
ignore.positivity$percent <- round(100*(ignore.positivity$n/nrow(ignore.positivity)),2)
ignore.positivity <- ignore.positivity[,2:4]
colnames(ignore.positivity) <- c("Isotype", "n","percent")
ignore.positivity <- unique(ignore.positivity)
ignore.positivity$Isotype <- factor(ignore.positivity$Isotype, levels = c("IGHG4","IGHG3","IGHG2","IGHG1","IGHA1","IGHM","IGHD"))


a <- ggplot(ignore.positivity, aes(fill=Isotype, y=percent, x="All cells")) + 
  geom_bar(position="fill", stat="identity")+ theme_classic()+labs(x="" ,y="Fraction")+
  scale_y_continuous(expand = c(0,0))+theme(legend.position = "none",axis.text.x = element_text(color="black"))+scale_x_discrete(labels = "All cells")+scale_fill_manual(values=met.brewer("Austria", 8))+
  ggtitle("All B cells")+ center.title()

remove(ignore.positivity)

positive.cells <- rownames(Bmem.Inc@meta.data)[Bmem.Inc@meta.data$Full.Row.names %in% all.cells$cell_id & Bmem.Inc@meta.data$bait.positive=="yes"]

positive.cells.seurat <- Bmem.Inc@meta.data[rownames(Bmem@meta.data)%in% positive.cells,c("named.clusters","Timepoint")]
positive.cells.imc <- all.cells[all.cells$cell_id %in% positive.cells,c("cell_id","c_call")] 
rownames(positive.cells.imc) <- positive.cells.imc$cell_id
positive.cells <- merge(positive.cells.seurat, positive.cells.imc,by=0)
positive.cells <- positive.cells %>% group_by(c_call) %>% add_tally()%>% ungroup()
positive.cells$percent <- round(100*(positive.cells$n/nrow(positive.cells)),2)


ignore.phenotypes <- positive.cells[,c(5,7)]
ignore.phenotypes <- as.data.frame(ignore.phenotypes)
ignore.phenotypes <- unique(ignore.phenotypes)
colnames(ignore.phenotypes) <- c("Isotype","percent")
ignore.phenotypes$Isotype <- factor(ignore.phenotypes$Isotype, levels = c("IGHG4","IGHG3","IGHG2","IGHG1","IGHA1","IGHM","IGHD"))

b <- ggplot(ignore.phenotypes, aes(fill=Isotype, y=percent, x="Bait positive cells")) + 
  geom_bar(position="fill", stat="identity")+ theme_classic()+labs(x="",y="Fraction")+
  scale_y_continuous(expand = c(0,0))+theme(legend.position = "none",axis.text.x = element_text(color="black"))+scale_x_discrete(labels = "Bait positive cells")+scale_fill_manual(values=met.brewer("Austria", 8))+
  ggtitle("Bait positive cells")+ center.title()

remove(ignore.phenotypes)



positive.cells <- rownames(Bmem.Inc@meta.data)[Bmem.Inc@meta.data$Full.Row.names %in% all.cells$cell_id & Bmem.Inc@meta.data$bait.positive=="yes"]
positive.cells.seurat <- Bmem.Inc@meta.data[rownames(Bmem@meta.data)%in% positive.cells,c("named.clusters","Timepoint")]
positive.cells.imc <- all.cells[all.cells$cell_id %in% positive.cells,c("cell_id","c_call")] 
rownames(positive.cells.imc) <- positive.cells.imc$cell_id
positive.cells <- merge(positive.cells.seurat, positive.cells.imc,by=0)
positive.cells <- positive.cells[,c(2,5)]
colnames(positive.cells) <- c("named_cluster","Isotype")
positive.cells <- positive.cells %>% group_by(named_cluster,Isotype) %>% add_tally() %>% ungroup
positive.cells <- as.data.frame(positive.cells)
positive.cells <- positive.cells %>% group_by(named_cluster) %>% add_tally() %>% ungroup()
positive.cells$percent <- round(100*(positive.cells$n/positive.cells$nn),2)
positive.cells <- unique(positive.cells)
positive.cells$Isotype <- factor(positive.cells$Isotype, levels = c("IGHG4","IGHG3","IGHG2","IGHG1","IGHA1","IGHM","IGHD"))
positive.cells$named_cluster <- factor(positive.cells$named_cluster, levels=c("Unswitched","CD27low RM",
                                                                      "CD27high RM","Activated","Atypical"))

c <- ggplot(positive.cells, aes(fill=Isotype, y=percent, x=named_cluster)) + 
  geom_bar(position="fill", stat="identity")+ theme_classic()+labs(x="",y="Fraction")+
  scale_y_continuous(expand = c(0,0))+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, color="black"))+scale_x_discrete()+scale_fill_manual(values=met.brewer("Austria", 8))+
  ggtitle("Bait positive cells\n(per Named Cluster)")+ center.title() # Export 5x10

a+b+c
remove(a,all.cells,b,Bmem.Inc,c,positive.cells,positive.cells.imc,positive.cells.seurat)


###################################################################################################################################
# Part 5 - Making a plot that shows the moving of longitudinal clones between phenotypes - Fig.6.E
###################################################################################################################################
Persistents <- read.csv2("~/NAS/Jan/Experiments and Data/210801 - Exp 035 - Profiling of COVID vaccinated patient samples/Memory B cell Fourth Experiment 220323/multi_exp035_4/Phenotypes of cells belonging to longitudinal clones.csv") # Created in part 3 of this script
Persistents <- Persistents[,c(3,4,9)]
Persistents <- Persistents[with(Persistents, order(clone_id, Timepoint)),]

# Taking out clone 38_2 because it is not from a vaccinated patient
Persistents <- Persistents[Persistents$clone_id != "38_2",]

sectors <- c("Activated","Atypical","Unswitched","CD27+ RM","CD27- RM")
Matrix <- matrix(1:10, nrow = 5, ncol = 2)
rownames(Matrix) <- sectors
Matrix[1,] <- c(0,22)  # Activated
Matrix[2,] <- c(23,38) # Atypical
Matrix[3,] <- c(38,41) # Unswitched
Matrix[4,] <- c(41,75) # CD27+ RM
Matrix[5,] <- c(77,86) # CD27- RM

circos.initialize(sectors = sectors, xlim = Matrix)

circos.track(ylim=c(0,1), bg.col=c("#EE6677","#AA3377","#DBC35E","#4477AA","cornsilk"))
circos.text(11.5, 1.5, "Activated", sector.index = "Activated", track.index = 1,cex=1.2, facing = "bending")
circos.text(30.5, 1.5, "Atypical", sector.index = "Atypical", track.index = 1,cex=1.2, facing = "bending")
circos.text(39.5, 1.5, "Unswitched", sector.index = "Unswitched", track.index = 1,cex=1.2, facing = "bending")
circos.text(60.5, 1.5, "CD27+ RM", sector.index = "CD27+ RM", track.index = 1,cex=1.2, facing = "bending")
circos.text(81.5, 1.5, "CD27- RM", sector.index = "CD27- RM", track.index = 1,cex=1.2, facing = "bending")

# In the following, the circos plot is created quite manually in a clone by clone way.
# 117_16
circos.link(sector.index1 = "CD27+ RM",point1 = c(41,42),sector.index2 = "Atypical",point2 = c(23,24),directional = 1,arr.type="big.arrow",col=add_transparency("#AA3377", 0.3),border = "black")
# 143_294
circos.link(sector.index1 = "Activated",point1 = c(0,1),sector.index2 = "Atypical",point2 = c(24,25),directional = 1,arr.type="big.arrow",col=add_transparency("#AA3377", 0.3),border = "black")
# 147_8 splitter
circos.link(sector.index1 = "CD27+ RM",point1 = c(42,43),sector.index2 = "CD27+ RM",point2 = c(60,61),directional = 1,arr.type="big.arrow",col=add_transparency("#4477AA", 0.3),border = "black")
circos.link(sector.index1 = "CD27+ RM",point1 = c(42,43),sector.index2 = "Activated",point2 = c(1,2),directional = 1,arr.type="big.arrow",col=add_transparency("#EE6677", 0.3),border = "black")
# 149_10
circos.link(sector.index1 = "CD27+ RM",point1 = c(43,44),sector.index2 = "Atypical",point2 = c(25,27),directional = 1,arr.type="big.arrow",col=add_transparency("#AA3377", 0.3),border = "black")
# 159_414 splitter
circos.link(sector.index1 = "CD27+ RM",point1 = c(44,45),sector.index2 = "Atypical",point2 = c(27,28),directional = 1,arr.type="big.arrow",col=add_transparency("#AA3377", 0.3),border = "black")
circos.link(sector.index1 = "CD27+ RM",point1 = c(44,45),sector.index2 = "Activated",point2 = c(2,3),directional = 1,arr.type="big.arrow",col=add_transparency("#EE6677", 0.3),border = "black")
# 160_313
circos.link(sector.index1 = "CD27+ RM",point1 = c(45,46),sector.index2 = "Activated",point2 = c(3,4),directional = 1,arr.type="big.arrow",col=add_transparency("#EE6677", 0.3),border = "black")
# 205_139
circos.link(sector.index1 = "CD27- RM",point1 = c(77,78),sector.index2 = "Activated",point2 = c(4,5),directional = 1,arr.type="big.arrow",col=add_transparency("#EE6677", 0.3),border = "black")
# 233_4
circos.link(sector.index1 = "CD27+ RM",point1 = c(46,47),sector.index2 = "Atypical",point2 = c(28,29),directional = 1,arr.type="big.arrow",col=add_transparency("#AA3377", 0.3),border = "black")
# 238_43
circos.link(sector.index1 = "CD27+ RM",point1 = c(47,48),sector.index2 = "Activated",point2 = c(5,6),directional = 1,arr.type="big.arrow",col=add_transparency("#EE6677", 0.3),border = "black")
# 251_325
circos.link(sector.index1 = "Activated",point1 = c(6,7),sector.index2 = "Activated",point2 = c(18,19),directional = 1,arr.type="big.arrow",col=add_transparency("#EE6677", 0.3),border = "black")
# 267_293
circos.link(sector.index1 = "CD27+ RM",point1 = c(48,49),sector.index2 = "CD27+ RM",point2 = c(65,66),directional = 1,arr.type="big.arrow",col=add_transparency("#4477AA", 0.3),border = "black")
# 285_88
circos.link(sector.index1 = "Activated",point1 = c(7,8),sector.index2 = "Atypical",point2 = c(29,30),directional = 1,arr.type="big.arrow",col=add_transparency("#AA3377", 0.3),border = "black")
# 29_0 splitter
circos.link(sector.index1 = "CD27+ RM",point1 = c(49,50),sector.index2 = "CD27+ RM",point2 = c(66,67),directional = 1,arr.type="big.arrow",col=add_transparency("#4477AA", 0.3),border = "black")
circos.link(sector.index1 = "CD27+ RM",point1 = c(49,50),sector.index2 = "Activated",point2 = c(8,9),directional = 1,arr.type="big.arrow",col=add_transparency("#EE6677", 0.3),border = "black")
# 29_7
circos.link(sector.index1 = "CD27+ RM",point1 = c(50,51),sector.index2 = "CD27+ RM",point2 = c(71,72),directional = 1,arr.type="big.arrow",col=add_transparency("#4477AA", 0.3),border = "black")
# 313_89
circos.link(sector.index1 = "Activated",point1 = c(9,10),sector.index2 = "CD27+ RM",point2 = c(51,52),directional = 1,arr.type="big.arrow",col=add_transparency("#4477AA", 0.3),border = "black")
# 314_256
circos.link(sector.index1 = "CD27- RM",point1 = c(78,79),sector.index2 = "Atypical",point2 = c(30,31),directional = 1,arr.type="big.arrow",col=add_transparency("#AA3377", 0.3),border = "black")
# 322_198
circos.link(sector.index1 = "CD27- RM",point1 = c(79,80),sector.index2 = "Atypical",point2 = c(31,32),directional = 1,arr.type="big.arrow",col=add_transparency("#AA3377", 0.3),border = "black")
# 351_492
circos.link(sector.index1 = "CD27- RM",point1 = c(80,81),sector.index2 = "Activated",point2 = c(10,11),directional = 1,arr.type="big.arrow",col=add_transparency("#EE6677", 0.3),border = "black")
# 357_26
circos.link(sector.index1 = "Activated",point1 = c(11,12),sector.index2 = "Atypical",point2 = c(32,33),directional = 1,arr.type="big.arrow",col=add_transparency("#AA3377", 0.3),border = "black")
# 358_69 splitter
circos.link(sector.index1 = "Activated",point1 = c(12,13),sector.index2 = "Unswitched",point2 = c(38,39),directional = 1,arr.type="big.arrow",col=add_transparency("#DBC35E", 0.3),border = "black")
circos.link(sector.index1 = "Activated",point1 = c(12,13),sector.index2 = "CD27+ RM",point2 = c(52,54),directional = 1,arr.type="big.arrow",col=add_transparency("#4477AA", 0.3),border = "black")
circos.link(sector.index1 = "Activated",point1 = c(12,13),sector.index2 = "Activated",point2 = c(20,22),directional = 1,arr.type="big.arrow",col=add_transparency("#EE6677", 0.3),border = "black")
# 368_22 splitter
circos.link(sector.index1 = "Unswitched",point1 = c(39,40),sector.index2 = "Unswitched",point2 = c(40,41),directional = 1,arr.type="big.arrow",col=add_transparency("#DBC35E", 0.3),border = "black")
circos.link(sector.index1 = "Unswitched",point1 = c(39,40),sector.index2 = "CD27+ RM",point2 = c(54,56),directional = 1,arr.type="big.arrow",col=add_transparency("#4477AA", 0.3),border = "black")
circos.link(sector.index1 = "Unswitched",point1 = c(39,40),sector.index2 = "Activated",point2 = c(13,15),directional = 1,arr.type="big.arrow",col=add_transparency("#EE6677", 0.3),border = "black")
circos.link(sector.index1 = "Unswitched",point1 = c(39,40),sector.index2 = "Atypical",point2 = c(33,34),directional = 1,arr.type="big.arrow",col=add_transparency("#AA3377", 0.3),border = "black")
# 402_82 splitter
circos.link(sector.index1 = "CD27- RM",point1 = c(81,82),sector.index2 = "CD27+ RM",point2 = c(56,57),directional = 1,arr.type="big.arrow",col=add_transparency("#4477AA", 0.3),border = "black")
circos.link(sector.index1 = "CD27- RM",point1 = c(81,82),sector.index2 = "Activated",point2 = c(15,16),directional = 1,arr.type="big.arrow",col=add_transparency("#EE6677", 0.3),border = "black")
circos.link(sector.index1 = "CD27- RM",point1 = c(81,82),sector.index2 = "Atypical",point2 = c(34,35),directional = 1,arr.type="big.arrow",col=add_transparency("#AA3377", 0.3),border = "black")
# 41_363
circos.link(sector.index1 = "CD27+ RM",point1 = c(57,58),sector.index2 = "Atypical",point2 = c(35,36),directional = 1,arr.type="big.arrow",col=add_transparency("#AA3377", 0.3),border = "black")
# 427_160
circos.link(sector.index1 = "CD27- RM",point1 = c(82,83),sector.index2 = "CD27+ RM",point2 = c(58,60),directional = 1,arr.type="big.arrow",col=add_transparency("#4477AA", 0.3),border = "black")
# 428_306
circos.link(sector.index1 = "CD27- RM",point1 = c(83,84),sector.index2 = "CD27+ RM",point2 = c(61,62),directional = 1,arr.type="big.arrow",col=add_transparency("#4477AA", 0.3),border = "black")
# 47_128
circos.link(sector.index1 = "CD27+ RM",point1 = c(62,63),sector.index2 = "Activated",point2 = c(16,18),directional = 1,arr.type="big.arrow",col=add_transparency("#EE6677", 0.3),border = "black")
# 48_229
circos.link(sector.index1 = "CD27+ RM",point1 = c(63,64),sector.index2 = "Atypical",point2 = c(36,37),directional = 1,arr.type="big.arrow",col=add_transparency("#AA3377", 0.3),border = "black")
circos.link(sector.index1 = "CD27- RM",point1 = c(84,85),sector.index2 = "Atypical",point2 = c(36,37),directional = 1,arr.type="big.arrow",col=add_transparency("#AA3377", 0.3),border = "black")
# 507_397
circos.link(sector.index1 = "CD27+ RM",point1 = c(67,70),sector.index2 = "CD27+ RM",point2 = c(72,75),directional = 1,arr.type="big.arrow",col=add_transparency("#4477AA", 0.3),border = "black")
# 53_63
circos.link(sector.index1 = "CD27+ RM",point1 = c(64,65),sector.index2 = "Atypical",point2 = c(37,38),directional = 1,arr.type="big.arrow",col=add_transparency("#AA3377", 0.3),border = "black")
# 85_142 splitter
circos.link(sector.index1 = "CD27- RM",point1 = c(85,86),sector.index2 = "CD27+ RM",point2 = c(70,71),directional = 1,arr.type="big.arrow",col=add_transparency("#4477AA", 0.3),border = "black")
circos.link(sector.index1 = "CD27- RM",point1 = c(85,86),sector.index2 = "Activated",point2 = c(19,20),directional = 1,arr.type="big.arrow",col=add_transparency("#EE6677", 0.3),border = "black")



# Second version of the plot, making only clones that split into different subsets colorful. The other clones are white.
circos.initialize(sectors = sectors, xlim = Matrix)

circos.track(ylim=c(0,1), bg.col=c("#EE6677","#AA3377","#DBC35E","#4477AA","cornsilk"))
circos.text(11.5, 1.5, "Activated", sector.index = "Activated", track.index = 1,cex=1.2, facing = "bending")
circos.text(30.5, 1.5, "Atypical", sector.index = "Atypical", track.index = 1,cex=1.2, facing = "bending")
circos.text(39.5, 1.5, "Unswitched", sector.index = "Unswitched", track.index = 1,cex=1.2, facing = "bending")
circos.text(60.5, 1.5, "CD27+ RM", sector.index = "CD27+ RM", track.index = 1,cex=1.2, facing = "bending")
circos.text(85, 1.5, "CD27- RM", sector.index = "CD27- RM", track.index = 1,cex=1.2, facing = "bending")


# 117_16
circos.link(sector.index1 = "CD27+ RM",point1 = c(41,42),sector.index2 = "Atypical",point2 = c(23,24),directional = 1,arr.type="big.arrow",col=add_transparency("gray97", 0.3),border = "black")
# 143_294
circos.link(sector.index1 = "Activated",point1 = c(0,1),sector.index2 = "Atypical",point2 = c(24,25),directional = 1,arr.type="big.arrow",col=add_transparency("gray97", 0.3),border = "black")
# 147_8 splitter
circos.link(sector.index1 = "CD27+ RM",point1 = c(42,43),sector.index2 = "CD27+ RM",point2 = c(61,62),directional = 1,arr.type="big.arrow",col=add_transparency("#4477AA", 0.3),border = "black")
circos.link(sector.index1 = "CD27+ RM",point1 = c(42,43),sector.index2 = "Activated",point2 = c(1,2),directional = 1,arr.type="big.arrow",col=add_transparency("#EE6677", 0.3),border = "black")
# 149_10
circos.link(sector.index1 = "CD27+ RM",point1 = c(43,44),sector.index2 = "Atypical",point2 = c(25,27),directional = 1,arr.type="big.arrow",col=add_transparency("gray97", 0.3),border = "black")
# 159_414 splitter
circos.link(sector.index1 = "CD27+ RM",point1 = c(44,45),sector.index2 = "Atypical",point2 = c(27,28),directional = 1,arr.type="big.arrow",col=add_transparency("#AA3377", 0.3),border = "black")
circos.link(sector.index1 = "CD27+ RM",point1 = c(44,45),sector.index2 = "Activated",point2 = c(2,3),directional = 1,arr.type="big.arrow",col=add_transparency("#EE6677", 0.3),border = "black")
# 160_313
circos.link(sector.index1 = "CD27+ RM",point1 = c(45,46),sector.index2 = "Activated",point2 = c(3,4),directional = 1,arr.type="big.arrow",col=add_transparency("gray97", 0.3),border = "black")
# 205_139
circos.link(sector.index1 = "CD27- RM",point1 = c(80,81),sector.index2 = "Activated",point2 = c(4,5),directional = 1,arr.type="big.arrow",col=add_transparency("gray97", 0.3),border = "black")
# 233_4
circos.link(sector.index1 = "CD27+ RM",point1 = c(46,47),sector.index2 = "Atypical",point2 = c(28,29),directional = 1,arr.type="big.arrow",col=add_transparency("gray97", 0.3),border = "black")
# 238_43
circos.link(sector.index1 = "CD27+ RM",point1 = c(47,48),sector.index2 = "Activated",point2 = c(5,6),directional = 1,arr.type="big.arrow",col=add_transparency("gray97", 0.3),border = "black")
# 251_325
circos.link(sector.index1 = "Activated",point1 = c(6,7),sector.index2 = "Activated",point2 = c(18,19),directional = 1,arr.type="big.arrow",col=add_transparency("gray97", 0.3),border = "black")
# 267_293
circos.link(sector.index1 = "CD27+ RM",point1 = c(48,49),sector.index2 = "CD27+ RM",point2 = c(66,67),directional = 1,arr.type="big.arrow",col=add_transparency("gray97", 0.3),border = "black")
# 285_88
circos.link(sector.index1 = "Activated",point1 = c(7,8),sector.index2 = "Atypical",point2 = c(29,30),directional = 1,arr.type="big.arrow",col=add_transparency("gray97", 0.3),border = "black")
# 29_0 splitter
circos.link(sector.index1 = "CD27+ RM",point1 = c(49,50),sector.index2 = "CD27+ RM",point2 = c(67,68),directional = 1,arr.type="big.arrow",col=add_transparency("#4477AA", 0.3),border = "black")
circos.link(sector.index1 = "CD27+ RM",point1 = c(49,50),sector.index2 = "Activated",point2 = c(8,9),directional = 1,arr.type="big.arrow",col=add_transparency("#EE6677", 0.3),border = "black")
# 29_7
circos.link(sector.index1 = "CD27+ RM",point1 = c(50,51),sector.index2 = "CD27+ RM",point2 = c(72,73),directional = 1,arr.type="big.arrow",col=add_transparency("gray97", 0.3),border = "black")
# 313_89
circos.link(sector.index1 = "Activated",point1 = c(9,10),sector.index2 = "CD27+ RM",point2 = c(51,52),directional = 1,arr.type="big.arrow",col=add_transparency("gray97", 0.3),border = "black")
# 314_256
circos.link(sector.index1 = "CD27- RM",point1 = c(81,82),sector.index2 = "Atypical",point2 = c(30,31),directional = 1,arr.type="big.arrow",col=add_transparency("gray97", 0.3),border = "black")
# 322_198
circos.link(sector.index1 = "CD27- RM",point1 = c(82,83),sector.index2 = "Atypical",point2 = c(31,32),directional = 1,arr.type="big.arrow",col=add_transparency("gray97", 0.3),border = "black")
# 351_492
circos.link(sector.index1 = "CD27- RM",point1 = c(83,84),sector.index2 = "Activated",point2 = c(10,11),directional = 1,arr.type="big.arrow",col=add_transparency("gray97", 0.3),border = "black")
# 357_26
circos.link(sector.index1 = "Activated",point1 = c(11,12),sector.index2 = "Atypical",point2 = c(32,33),directional = 1,arr.type="big.arrow",col=add_transparency("gray97", 0.3),border = "black")
# 358_69 splitter
circos.link(sector.index1 = "Activated",point1 = c(12,13),sector.index2 = "Unswitched",point2 = c(38,39),directional = 1,arr.type="big.arrow",col=add_transparency("#DBC35E", 0.3),border = "black")
circos.link(sector.index1 = "Activated",point1 = c(12,13),sector.index2 = "CD27+ RM",point2 = c(52,54),directional = 1,arr.type="big.arrow",col=add_transparency("#4477AA", 0.3),border = "black")
circos.link(sector.index1 = "Activated",point1 = c(12,13),sector.index2 = "Activated",point2 = c(20,22),directional = 1,arr.type="big.arrow",col=add_transparency("#EE6677", 0.3),border = "black")
# 368_22 splitter
circos.link(sector.index1 = "Unswitched",point1 = c(39,40),sector.index2 = "Unswitched",point2 = c(40,41),directional = 1,arr.type="big.arrow",col=add_transparency("#DBC35E", 0.3),border = "black")
circos.link(sector.index1 = "Unswitched",point1 = c(39,40),sector.index2 = "CD27+ RM",point2 = c(54,56),directional = 1,arr.type="big.arrow",col=add_transparency("#4477AA", 0.3),border = "black")
circos.link(sector.index1 = "Unswitched",point1 = c(39,40),sector.index2 = "Activated",point2 = c(13,15),directional = 1,arr.type="big.arrow",col=add_transparency("#EE6677", 0.3),border = "black")
circos.link(sector.index1 = "Unswitched",point1 = c(39,40),sector.index2 = "Atypical",point2 = c(33,34),directional = 1,arr.type="big.arrow",col=add_transparency("#AA3377", 0.3),border = "black")
# 38_2
circos.link(sector.index1 = "CD27+ RM",point1 = c(56,57),sector.index2 = "CD27+ RM",point2 = c(73,75),directional = 1,arr.type="big.arrow",col=add_transparency("gray97", 0.3),border = "black")
# 402_82 splitter
circos.link(sector.index1 = "CD27- RM",point1 = c(84,85),sector.index2 = "CD27+ RM",point2 = c(57,58),directional = 1,arr.type="big.arrow",col=add_transparency("#4477AA", 0.3),border = "black")
circos.link(sector.index1 = "CD27- RM",point1 = c(84,85),sector.index2 = "Activated",point2 = c(15,16),directional = 1,arr.type="big.arrow",col=add_transparency("#EE6677", 0.3),border = "black")
circos.link(sector.index1 = "CD27- RM",point1 = c(84,85),sector.index2 = "Atypical",point2 = c(34,35),directional = 1,arr.type="big.arrow",col=add_transparency("#AA3377", 0.3),border = "black")
# 41_363
circos.link(sector.index1 = "CD27+ RM",point1 = c(58,59),sector.index2 = "Atypical",point2 = c(35,36),directional = 1,arr.type="big.arrow",col=add_transparency("gray97", 0.3),border = "black")
# 427_160
circos.link(sector.index1 = "CD27- RM",point1 = c(85,86),sector.index2 = "CD27+ RM",point2 = c(59,61),directional = 1,arr.type="big.arrow",col=add_transparency("gray97", 0.3),border = "black")
# 428_306
circos.link(sector.index1 = "CD27- RM",point1 = c(86,87),sector.index2 = "CD27+ RM",point2 = c(62,63),directional = 1,arr.type="big.arrow",col=add_transparency("gray97", 0.3),border = "black")
# 47_128
circos.link(sector.index1 = "CD27+ RM",point1 = c(63,64),sector.index2 = "Activated",point2 = c(16,18),directional = 1,arr.type="big.arrow",col=add_transparency("gray97", 0.3),border = "black")
# 48_229
circos.link(sector.index1 = "CD27+ RM",point1 = c(64,65),sector.index2 = "Atypical",point2 = c(36,37),directional = 1,arr.type="big.arrow",col=add_transparency("gray97", 0.3),border = "black")
circos.link(sector.index1 = "CD27- RM",point1 = c(87,88),sector.index2 = "Atypical",point2 = c(36,37),directional = 1,arr.type="big.arrow",col=add_transparency("gray97", 0.3),border = "black")
# 507_397
circos.link(sector.index1 = "CD27+ RM",point1 = c(68,71),sector.index2 = "CD27+ RM",point2 = c(75,78),directional = 1,arr.type="big.arrow",col=add_transparency("gray97", 0.3),border = "black")
# 53_63
circos.link(sector.index1 = "CD27+ RM",point1 = c(65,66),sector.index2 = "Atypical",point2 = c(37,38),directional = 1,arr.type="big.arrow",col=add_transparency("gray97", 0.3),border = "black")
# 85_142 splitter
circos.link(sector.index1 = "CD27- RM",point1 = c(88,89),sector.index2 = "CD27+ RM",point2 = c(71,72),directional = 1,arr.type="big.arrow",col=add_transparency("#4477AA", 0.3),border = "black")
circos.link(sector.index1 = "CD27- RM",point1 = c(88,89),sector.index2 = "Activated",point2 = c(19,20),directional = 1,arr.type="big.arrow",col=add_transparency("#EE6677", 0.3),border = "black")



# Third version of the Circos plot, only splitting clones are shown, the other clones are invisible.
circos.initialize(sectors = sectors, xlim = Matrix)

circos.track(ylim=c(0,1), bg.col=c("#EE6677","#AA3377","#DBC35E","#4477AA","cornsilk"))
circos.text(11.5, 1.5, "Activated", sector.index = "Activated", track.index = 1,cex=1.2, facing = "bending")
circos.text(30.5, 1.5, "Atypical", sector.index = "Atypical", track.index = 1,cex=1.2, facing = "bending")
circos.text(39.5, 1.5, "Unswitched", sector.index = "Unswitched", track.index = 1,cex=1.2, facing = "bending")
circos.text(60.5, 1.5, "CD27+ RM", sector.index = "CD27+ RM", track.index = 1,cex=1.2, facing = "bending")
circos.text(81.5, 1.5, "CD27- RM", sector.index = "CD27- RM", track.index = 1,cex=1.2, facing = "bending")


# 147_8 splitter
circos.link(sector.index1 = "CD27+ RM",point1 = c(42,43),sector.index2 = "CD27+ RM",point2 = c(60,61),directional = 1,arr.type="big.arrow",col=add_transparency("#4477AA", 0.3),border = "black")
circos.link(sector.index1 = "CD27+ RM",point1 = c(42,43),sector.index2 = "Activated",point2 = c(1,2),directional = 1,arr.type="big.arrow",col=add_transparency("#EE6677", 0.3),border = "black")
# 159_414 splitter
circos.link(sector.index1 = "CD27+ RM",point1 = c(44,45),sector.index2 = "Atypical",point2 = c(27,28),directional = 1,arr.type="big.arrow",col=add_transparency("#AA3377", 0.3),border = "black")
circos.link(sector.index1 = "CD27+ RM",point1 = c(44,45),sector.index2 = "Activated",point2 = c(2,3),directional = 1,arr.type="big.arrow",col=add_transparency("#EE6677", 0.3),border = "black")
# 29_0 splitter
circos.link(sector.index1 = "CD27+ RM",point1 = c(49,50),sector.index2 = "CD27+ RM",point2 = c(66,67),directional = 1,arr.type="big.arrow",col=add_transparency("#4477AA", 0.3),border = "black")
circos.link(sector.index1 = "CD27+ RM",point1 = c(49,50),sector.index2 = "Activated",point2 = c(8,9),directional = 1,arr.type="big.arrow",col=add_transparency("#EE6677", 0.3),border = "black")
# 358_69 splitter
circos.link(sector.index1 = "Activated",point1 = c(12,13),sector.index2 = "Unswitched",point2 = c(38,39),directional = 1,arr.type="big.arrow",col=add_transparency("#DBC35E", 0.3),border = "black")
circos.link(sector.index1 = "Activated",point1 = c(12,13),sector.index2 = "CD27+ RM",point2 = c(52,54),directional = 1,arr.type="big.arrow",col=add_transparency("#4477AA", 0.3),border = "black")
circos.link(sector.index1 = "Activated",point1 = c(12,13),sector.index2 = "Activated",point2 = c(20,22),directional = 1,arr.type="big.arrow",col=add_transparency("#EE6677", 0.3),border = "black")
# 368_22 splitter
circos.link(sector.index1 = "Unswitched",point1 = c(39,40),sector.index2 = "Unswitched",point2 = c(40,41),directional = 1,arr.type="big.arrow",col=add_transparency("#DBC35E", 0.3),border = "black")
circos.link(sector.index1 = "Unswitched",point1 = c(39,40),sector.index2 = "CD27+ RM",point2 = c(54,56),directional = 1,arr.type="big.arrow",col=add_transparency("#4477AA", 0.3),border = "black")
circos.link(sector.index1 = "Unswitched",point1 = c(39,40),sector.index2 = "Activated",point2 = c(13,15),directional = 1,arr.type="big.arrow",col=add_transparency("#EE6677", 0.3),border = "black")
circos.link(sector.index1 = "Unswitched",point1 = c(39,40),sector.index2 = "Atypical",point2 = c(33,34),directional = 1,arr.type="big.arrow",col=add_transparency("#AA3377", 0.3),border = "black")
# 402_82 splitter
circos.link(sector.index1 = "CD27- RM",point1 = c(81,82),sector.index2 = "CD27+ RM",point2 = c(56,57),directional = 1,arr.type="big.arrow",col=add_transparency("#4477AA", 0.3),border = "black")
circos.link(sector.index1 = "CD27- RM",point1 = c(81,82),sector.index2 = "Activated",point2 = c(15,16),directional = 1,arr.type="big.arrow",col=add_transparency("#EE6677", 0.3),border = "black")
circos.link(sector.index1 = "CD27- RM",point1 = c(81,82),sector.index2 = "Atypical",point2 = c(34,35),directional = 1,arr.type="big.arrow",col=add_transparency("#AA3377", 0.3),border = "black")
# 85_142 splitter
circos.link(sector.index1 = "CD27- RM",point1 = c(85,86),sector.index2 = "CD27+ RM",point2 = c(70,71),directional = 1,arr.type="big.arrow",col=add_transparency("#4477AA", 0.3),border = "black")
circos.link(sector.index1 = "CD27- RM",point1 = c(85,86),sector.index2 = "Activated",point2 = c(19,20),directional = 1,arr.type="big.arrow",col=add_transparency("#EE6677", 0.3),border = "black")

remove(Matrix,Persistents,sectors)

###################################################################################################################################
# Part 6 - Paired Boxplots / Violinplots - Fig.6.C
###################################################################################################################################
# Paired plot looking at Gene expression levels of selected marker genes.
setwd("~/NAS/Jan/Experiments and Data/210801 - Exp 035 - Profiling of COVID vaccinated patient samples/Memory B cell Fourth Experiment 220323/multi_exp035_4")  
Bmem <- readRDS("Bmem.rds")
Persistents <- read.csv2("~/NAS/Jan/Experiments and Data/210801 - Exp 035 - Profiling of COVID vaccinated patient samples/Memory B cell Fourth Experiment 220323/multi_exp035_4/Phenotypes of cells belonging to longitudinal clones.csv")# Created in part 3 of this script
Persistents <- Persistents[Persistents$Patient != "PFCL1-LIM-674950",] # Kicking out the non vaccinated patient.
Persistents <- Persistents[,c(2:4,9)]

# Select the gene to analyze.
gene <- "FCRL5"

# 1) Select the groups of cells you want to compare
#    In the first comparison we want to show how cells of persistent clones acquire expression of atypical markers after vaccination
#    So in the plot there will be two violins, one for each timepoint. In the right violin there will be one dot for each clone which 
#    has atypical cells in the 12 month timepoint. The position of the dot will represent the average expression level of a feature for all cells
#    of that clone in this timepoint. In the right violin there will be the same clones represented, again one dot for each. The position of 
#    the dot will again represent the average expression level of a feature for all cells of that clone in this timepoint.

Clones.involved <- unique(Persistents$clone_id[Persistents$named.clusters=="Atypical" & Persistents$Timepoint=="12 months after infection"])
six_mo_cells <- Persistents$cell_id[Persistents$clone_id %in% Clones.involved & Persistents$Timepoint=="6 months after infection"]
twelve_mo_cells <- Persistents$cell_id[Persistents$clone_id %in% Clones.involved & Persistents$Timepoint=="12 months after infection"]

six_mo_cells_df <- Persistents[Persistents$cell_id %in% six_mo_cells,1:2]
twelve_mo_cells_df <- Persistents[Persistents$cell_id %in% twelve_mo_cells,1:2]

six_mo_cells_df$expressionlevel <- NA
twelve_mo_cells_df$expressionlevel <- NA

rownames(six_mo_cells_df) <- c(1:nrow(six_mo_cells_df))
rownames(twelve_mo_cells_df) <- c(1:nrow(twelve_mo_cells_df))

# 2) Get expression levels of a gene for the cells selected
i <- 1
for (i in i:nrow(six_mo_cells_df)) {
  six_mo_cells_df[i,3] <- Bmem@assays$RNA@data[gene,six_mo_cells_df[i,1]]
}
i <- 1
for (i in i:nrow(twelve_mo_cells_df)) {
  twelve_mo_cells_df[i,3] <- Bmem@assays$RNA@data[gene,twelve_mo_cells_df[i,1]]
}

# 3) Calculate the means per clone_id and rbind frames
six_mo_cells_df <- aggregate(six_mo_cells_df[, 3], list(six_mo_cells_df$clone_id), mean)
twelve_mo_cells_df <- aggregate(twelve_mo_cells_df[, 3], list(twelve_mo_cells_df$clone_id), mean)
six_mo_cells_df$Timepoint <- "6 months"
twelve_mo_cells_df$Timepoint <- "12 months"
cells_df <- rbind(six_mo_cells_df,twelve_mo_cells_df)
cells_df$Timepoint <- factor(cells_df$Timepoint, levels=c("6 months","12 months"))

# 4) Plot # Export 5x4
my_comparisons <- list( c("6 months", "12 months"))
ggplot(cells_df, aes(x = Timepoint, y = x,group=Group.1)) + geom_violin(cells_df, mapping= aes(x = Timepoint, y = x,group=Timepoint))+
  geom_line(alpha = 0.8) + geom_point(alpha = 0.7, size = 3)+ 
  theme_classic()+labs(x="",y="Average Expression Level")+ stat_compare_means(comparisons = my_comparisons, method = "t.test",paired = T)+
  black.axis.text()+ theme(plot.title = element_text(face="bold"),text = element_text(size = 20))+ggtitle(gene)+ center.title()

write.xlsx(cells_df, paste0("~/NAS/Jan/Experiments and Data/210801 - Exp 035 - Profiling of COVID vaccinated patient samples/Memory B cell Fourth Experiment 220323/multi_exp035_4/",gene,"_expression_change_clonal_level.xlsx"))

remove(six_mo_cells_df,twelve_mo_cells_df,six_mo_cells,twelve_mo_cells,my_comparisons,Persistents,Clones.involved,cells_df,gene,i)


###################################################################################################################################
# Part 7 - Clonal diversity and expansion analysis - clonotyping done with heavy and light chain - Fig.4.D
###################################################################################################################################
# I work with the same expansion numbers as they are calculated in the donut charts. I was thinking about whether it would be better to treat all 
# the timepoints individually here (as in an older version of this script - see backups) or to follow the same rational as for the donut charts 
# --> one clonotyping for all cells of a patient together.
# I came to the conclusion that this second way is also a correct way of doing it because: All cells from a patient which are clonally
# related within a timepoint, should also be clonally related when both timepoints are mixed - no clonal expansion should be lost through the
# merging of timepoints.
donut_clonality  <- data.frame("Patient"=c("A","B","C","D","E","F"),"Time" = c("6 months","6 months","6 months","6 months","6 months","6 months",
                                                                               "12 months + vaccinated","12 months + vaccinated","12 months + vaccinated",
                                                                               "12 months + vaccinated","12 months + vaccinated","12 months + vaccinated"),"Percent.clonal"=c(0,2.1,19.1,2.9,0,0,13.1,6.1,39.6,18.4,13.8,0))
donut_clonality$Time <- factor(donut_clonality$Time, levels = c("6 months","12 months + vaccinated"))
my_comparisons <- list(c("6 months","12 months + vaccinated"))
ggplot(donut_clonality, aes(x=Time, y=Percent.clonal,group=Patient)) + geom_line(alpha = 0.8) + 
  geom_point(alpha = 0.7,size = 1.5) + theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+black.axis.text()+labs(x="",y="% Clonal")+
  ggtitle("Clonality of Spike specific cells \nVaccinated patients")+ center.title()+
  stat_compare_means(comparisons = my_comparisons, method = "t.test",paired = T) # Export 5x3

write.xlsx(donut_clonality, "~/NAS/Jan/Experiments and Data/210801 - Exp 035 - Profiling of COVID vaccinated patient samples/Memory B cell Fourth Experiment 220323/multi_exp035_4/Clonal.expansion.vaccinated.patients.xlsx")

###################################################################################################################################
# Part 8 - Analysis of cells belonging to persistent clones and the 12 months time point - Fig.6.B and part of Fig.6.D
###################################################################################################################################
# What are the frequencies of different cell states pre vac vs post vac in persistent clones INCLUDING NON PERSISTERS AFTER 12 MONTHS VACCINATED
# Way too much code than would be needed probably, but I am reusing code from other sections to save coding time.
Persistents <- read.csv2("~/NAS/Jan/Experiments and Data/210801 - Exp 035 - Profiling of COVID vaccinated patient samples/Memory B cell Fourth Experiment 220323/multi_exp035_4/Phenotypes of cells belonging to longitudinal clones.csv")# Created in part 3 of this script
# Kicking out the patient that is not vaccinated!
Persistents <- Persistents[Persistents$Patient != "PFCL1-LIM-674950",]
Persistents <- Persistents %>% group_by(Timepoint, named.clusters) %>% add_count(name = 'Cells per state and time') %>% ungroup()
Persistents <- Persistents %>% group_by(Timepoint) %>% add_count(name = 'Cells per time') %>% ungroup()
Persistents$Percent <- 100*round(Persistents$`Cells per state and time`/Persistents$`Cells per time`,2)
Persistents <- Persistents[,c("Timepoint","named.clusters","Percent")]
Persistents <- unique(Persistents)
True_Persistents <- Persistents

setwd("~/NAS/Jan/Experiments and Data/210801 - Exp 035 - Profiling of COVID vaccinated patient samples/Memory B cell Fourth Experiment 220323/multi_exp035_4")  
Bmem <- readRDS("Bmem.rds")
Persistents <- read.csv2("~/NAS/Jan/Experiments and Data/210801 - Exp 035 - Profiling of COVID vaccinated patient samples/Memory B cell Fourth Experiment 220323/multi_exp035_4/Phenotypes of cells belonging to longitudinal clones.csv")# Created in part 3 of this script
Patients <- c("PFCL1_LIM_313000","PFCL1_LIM828246","PFCL1_199474_1994","PFCL1_179308_1988")

setwd("~/NAS/Jan/Experiments and Data/210801 - Exp 035 - Profiling of COVID vaccinated patient samples/Immcantation Analysis/Immcantation Outputs/update_020522")

z <- 1
totalfile <- as.data.frame(matrix(nrow = 1,ncol = 63))
for (z in 1:length(Patients)) {
  
  Patient <- Patients[z]
  Thisroundfile<- as.data.frame(read_airr(paste0("./",Patient,"_positives/outputs/",Patient,"_positives_0.20/",Patient,"_positives_0.20_heavy_germ-pass.tsv")))
  Thisroundfile$Patient <- Patient
  colnames(totalfile) <- colnames(Thisroundfile)
  totalfile <- rbind(totalfile,Thisroundfile)
}
totalfile <- totalfile[2:nrow(totalfile),]

# Translating the cell_id column
i <- 1
for(i in 1:nrow(totalfile)){
  
  
  if(substr(totalfile[i,"cell_id"],1,2)=="CC"){
    totalfile[i,"cell_id"] <- paste0("Dataset_4_",substr(totalfile[i,"cell_id"],3,nchar(totalfile[i,"cell_id"])))
  }
  if(substr(totalfile[i,"cell_id"],1,2)=="AA"){
    totalfile[i,"cell_id"] <- paste0("Dataset_1_",substr(totalfile[i,"cell_id"],3,nchar(totalfile[i,"cell_id"])))
  }
  if(substr(totalfile[i,"cell_id"],1,2)=="GG"){
    totalfile[i,"cell_id"] <- paste0("Dataset_2_",substr(totalfile[i,"cell_id"],3,nchar(totalfile[i,"cell_id"])))
  }
  if(substr(totalfile[i,"cell_id"],1,2)=="TT"){
    totalfile[i,"cell_id"] <- paste0("Dataset_3_",substr(totalfile[i,"cell_id"],3,nchar(totalfile[i,"cell_id"])))
  }
}

# Making clear which cell belongs to a persistent clone
Persistents$Persistence <- "Persistent clones"
Persistents <- Persistents[,c("cell_id","Persistence")]
totalfile <- merge(totalfile,Persistents, by="cell_id",all.x = T)
totalfile$Persistence[is.na(totalfile$Persistence)] <- "New clones"

# Getting the timepoints and the cell state and filtering cells from 6 months
totalfile$Full.Row.names <- totalfile$cell_id
totalfile <- merge(totalfile, Bmem@meta.data, by="Full.Row.names",all.x = T)
totalfile$Timepoint
totalfile$cell_id
totalfile$named.clusters

# Simplifying the columns, only cell name and cell state are needed
totalfile <- totalfile[,c("cell_id","named.clusters","Timepoint","Persistence")]

# Kicking the ones from 6 months
totalfile <- totalfile[totalfile$Timepoint=="12 months after infection",]

# Grouping, calculating percentages
totalfile <- totalfile %>% group_by(Persistence) %>% add_count(name = "Cells_per_Persistence") %>% ungroup()
totalfile <- totalfile %>% group_by(Persistence,named.clusters) %>% add_count(name="Cells_per_Persistence_and_state") %>% ungroup()
totalfile <- totalfile[,2:ncol(totalfile)]
totalfile <- unique(totalfile)
totalfile$Percent <- round(100*totalfile$Cells_per_Persistence_and_state/totalfile$Cells_per_Persistence,2)

totalfile <- totalfile[totalfile$Persistence=="New clones",]
totalfile <- totalfile[,c(2,1,6)]
totalfile$Timepoint <- "New clones 12 mo. after infect. and vacc."

True_Persistents <- rbind(True_Persistents,totalfile)

True_Persistents$Timepoint <- factor(True_Persistents$Timepoint, levels = c("New clones 12 mo. after infect. and vacc.","12 months after infection","6 months after infection"))

ggplot(True_Persistents, aes(fill=factor(named.clusters, levels=c("Atypical","Activated",
                                                             "CD27high RM","CD27low RM","Unswitched")), y=Timepoint, x=Percent)) + 
  geom_bar(position="fill", stat="identity",colour="black")+theme_classic()+ 
  theme(text = element_text(size = 15),plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))+
  labs(x = "Subset percentages", y="")+scale_x_continuous(labels=scales::percent)+labs(fill = "Subset")+
  scale_fill_manual(values=rev(c("#DBC35E","#228833","#4477AA","#EE6677","#AA3377")))+ theme(text = element_text(size = 17))+ 
  ggtitle("Subsets of cells belonging to persistent clones \n before and after vaccination")+black.axis.text() # Export 5x10


# Little addition for the dendrograms: Getting the chains of the persistent clones
setwd("~/NAS/Jan/Experiments and Data/210801 - Exp 035 - Profiling of COVID vaccinated patient samples/Memory B cell Fourth Experiment 220323/multi_exp035_4")  
Bmem <- readRDS("Bmem.rds")
Persistents <- read.csv2("~/NAS/Jan/Experiments and Data/210801 - Exp 035 - Profiling of COVID vaccinated patient samples/Memory B cell Fourth Experiment 220323/multi_exp035_4/Phenotypes of cells belonging to longitudinal clones.csv")# Created in part 3 of this script
clone_402_82 <- Persistents[Persistents$clone_id=="402_82",]
clone_368_22 <- Persistents[Persistents$clone_id=="368_22",]
clone_29_0 <- Persistents[Persistents$clone_id=="29_0",]
clone_48_229 <- Persistents[Persistents$clone_id=="48_229",]
clone_159_414 <- Persistents[Persistents$clone_id=="159_414",]
clone_427_160 <- Persistents[Persistents$clone_id=="427_160",]

Segments_402_82 <- Bmem@meta.data[Bmem@meta.data$Full.Row.names %in% clone_402_82$cell_id,"CTgene"]
Segments_368_22 <- Bmem@meta.data[Bmem@meta.data$Full.Row.names %in% clone_368_22$cell_id,"CTgene"]
Segments_29_0 <- Bmem@meta.data[Bmem@meta.data$Full.Row.names %in% clone_29_0$cell_id,"CTgene"]
Segments_48_229 <- Bmem@meta.data[Bmem@meta.data$Full.Row.names %in% clone_48_229$cell_id,"CTgene"]
Segments_159_414 <- Bmem@meta.data[Bmem@meta.data$Full.Row.names %in% clone_159_414$cell_id,"CTgene"]

Bmem@meta.data[Bmem@meta.data$CTgene == "IGHV3-53..IGHJ6.IGHG1_IGKV1-9.IGKJ4.IGKC" &
                 Bmem@meta.data$Full.Row.names %in% Persistents$cell_id,c("CTgene","Full.Row.names")]
Bmem@meta.data[Bmem@meta.data$CTgene == "IGHV3-53..IGHJ4.IGHG1_IGKV1D-39.IGKJ2.IGKC" &
                 Bmem@meta.data$Full.Row.names %in% Persistents$cell_id,c("CTgene","Full.Row.names")]
Bmem@meta.data[Bmem@meta.data$CTgene == "IGHV3-30..IGHJ3.IGHG1_IGLV3-21.IGLJ2.IGLC2" &
                 Bmem@meta.data$Full.Row.names %in% Persistents$cell_id,c("CTgene","Full.Row.names")]
Bmem@meta.data[Bmem@meta.data$CTgene == "IGHV3-30.IGHD3-22.IGHJ4.IGHG1_IGKV3-15.IGKJ2.IGKC" &
                 Bmem@meta.data$Full.Row.names %in% Persistents$cell_id,c("CTgene","Full.Row.names")]
Bmem@meta.data[Bmem@meta.data$CTgene == "IGHV3-30..IGHJ6.IGHA1_IGKV1D-39.IGKJ3.IGKC" &
                 Bmem@meta.data$Full.Row.names %in% Persistents$cell_id,c("CTgene","Full.Row.names")]
Bmem@meta.data[Bmem@meta.data$CTgene == "IGHV3-53..IGHJ4.IGHG1_IGKV1D-39.IGKJ2.IGKC" &
                 Bmem@meta.data$Full.Row.names %in% Persistents$cell_id,c("CTgene","Full.Row.names")]



# What are the antigen specificites of the persistent clones?
Persistents <- read.csv2("~/NAS/Jan/Experiments and Data/210801 - Exp 035 - Profiling of COVID vaccinated patient samples/Memory B cell Fourth Experiment 220323/multi_exp035_4/Phenotypes of cells belonging to longitudinal clones.csv")# Created in part 3 of this script
Persistents_specificities <- Bmem@meta.data[Bmem@meta.data$Full.Row.names %in% Persistents$cell_id,c("cor_wt_Spike","cor_RBD","bait.positive","spike_RBD_positive",
                                                                                              "B.1.351_positive","B.1.617.2_positive","Timepoint")]
Persistents_specificities$Wt_spike_pos <- "no"
Persistents_specificities$Wt_spike_pos[Persistents_specificities$cor_wt_Spike>0] <- "yes"
Persistents_specificities$RBD_pos <- "no"
Persistents_specificities$RBD_pos[Persistents_specificities$cor_RBD>0] <- "yes"
rownames(Persistents) <- Persistents$cell_id
Persistents_specificities <- merge(Persistents_specificities,Persistents[,c("clone_id","cell_id")],by=0, all.x=T)

Persistents_specificities <- Persistents_specificities[,c("bait.positive","Wt_spike_pos","RBD_pos","spike_RBD_positive","B.1.351_positive","B.1.617.2_positive","Timepoint","clone_id")]
write.xlsx(Persistents_specificities, "Persistents_specificities.xlsx")


###################################################################################################################################
# Part 9 - Creating a Venn Diagrams - Supp.Fig.6.A and Supp.Fig.6.B
###################################################################################################################################
# Venn diagram showing the clonal overlap between time points in all patients
# For this I can reuse code from part 2:
Bmem.Inc <- readRDS("~/NAS/Jan/Experiments and Data/210801 - Exp 035 - Profiling of COVID vaccinated patient samples/Memory B cell Fourth Experiment 220323/multi_exp035_4/Bmem.rds")
setwd("~/NAS/Jan/Experiments and Data/210801 - Exp 035 - Profiling of COVID vaccinated patient samples/Immcantation Analysis/Immcantation Outputs/update_020522")

# Putting everything into a loop
Patients <- c("PFCL1_LIM_313000","PFCL1_LIM_137402","PFCL1_LIM828246","PFCL1_154583_1943","PFCL1_199474_1994","PFCL1_179308_1988",
              "PFCL1_196878_1691","PFCL1_LIM_674950","PFCL1_UST_190762")

# Creating output df
output_df <- data.frame(matrix(ncol = 4, nrow = 6))

#provide column names
colnames(output_df) <- c('Patient', 'n_non_persistents_6_months', 'n_non_persistents_12_months','n_persistents')

z <- 1
for (z in 1:length(Patients)) {

Patient <- Patients[z]
LIM828246 <- read_airr(paste0("./",Patient,"_positives/outputs/",Patient,"_positives_0.20/",Patient,"_positives_0.20_heavy_germ-pass.tsv"))
find.longitudinal.clones <- function(Immcantation.Output){
  IDs <- unique(Immcantation.Output$clone_id)
  output.vector <- c("")
  i <- 1
  for(i in 1:length(IDs)){
    
    candidates <- Immcantation.Output$cell_id[Immcantation.Output$clone_id ==IDs[i]]
    b <- 1
    for (b in 1:length(candidates)) {
      if(substr(candidates[b],1,2)=="CC"){
        candidates[b] <- paste0("Dataset_4_",substr(candidates[b],3,nchar(candidates[b])))
      }
      if(substr(candidates[b],1,2)=="AA"){
        candidates[b] <- paste0("Dataset_1_",substr(candidates[b],3,nchar(candidates[b])))
      }
      if(substr(candidates[b],1,2)=="GG"){
        candidates[b] <- paste0("Dataset_2_",substr(candidates[b],3,nchar(candidates[b])))
      }
      if(substr(candidates[b],1,2)=="TT"){
        candidates[b] <- paste0("Dataset_3_",substr(candidates[b],3,nchar(candidates[b])))
      }
    }
    
    timepoints <-Bmem.Inc@meta.data$Timepoint[rownames(Bmem.Inc@meta.data)%in%candidates]
    if(length(unique(timepoints))>1){
      output.vector <- c(output.vector,IDs[i])
    }
  }
  
  output.vector <- output.vector[2:length(output.vector)]
  return(output.vector)
  
}
LIM828246.longitudinals <- find.longitudinal.clones(LIM828246)


#
# This is a little addition to the copy pasted code from above. If there is no longitudinal clone found here, then all what comes can be ignored 
#
if(!is.na(LIM828246.longitudinals[1])){
  


# Annotate clones in the Immcantation output file.
# We want to make donut charts with the same style as Nussenzweig uses them.
# Individual pie chart pieces for all expanded clones that are persistent -> in color
# Here, make sure that all cells belonging to persistent clones get the same color (the same clone should have the same color in both timepoints)
# Individual pie chart pieces for all expanded clones that are unique -> grey
# Individual pie chart pieces for all singlets that are persistent (repeating sequences isolated only once per time point)
annotate.clones <- function(Longitudinals,Immcantation.Output,Seuratobject){
  IDs <- Longitudinals
  output.df <- data.frame(matrix(ncol = 2, nrow = 0))
  x <- c("Clone", "TruePersistent")
  colnames(output.df) <- x
  i <- 1
  for(i in 1:length(IDs)){
    
    candidates <- Immcantation.Output$cell_id[Immcantation.Output$clone_id ==IDs[i]]
    b <- 1
    for (b in 1:length(candidates)) {
      if(substr(candidates[b],1,2)=="CC"){
        candidates[b] <- paste0("Dataset_4_",substr(candidates[b],3,nchar(candidates[b])))
      }
      if(substr(candidates[b],1,2)=="AA"){
        candidates[b] <- paste0("Dataset_1_",substr(candidates[b],3,nchar(candidates[b])))
      }
      if(substr(candidates[b],1,2)=="GG"){
        candidates[b] <- paste0("Dataset_2_",substr(candidates[b],3,nchar(candidates[b])))
      }
      if(substr(candidates[b],1,2)=="TT"){
        candidates[b] <- paste0("Dataset_3_",substr(candidates[b],3,nchar(candidates[b])))
      }
    }
    
    Bmem.Inc <- Seuratobject
    candidates.metadata <-Bmem.Inc@meta.data[rownames(Bmem.Inc@meta.data)%in%candidates,c(1,18,19)]
    if(nrow(candidates.metadata)>2){
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
                                    Immcantation.Output$clone_id %notin% output.df$Clone] <- "Persistent clone"
  Immcantation.Output$persistence[Immcantation.Output$clone_id %in% Longitudinals &
                                    Immcantation.Output$clone_id %in% output.df$Clone] <- "Persistent clone"
  return(Immcantation.Output)
  
}

LIM828246.annotated <- annotate.clones(LIM828246.longitudinals, LIM828246, Bmem.Inc)


# Now we make the donut charts, one for each timepoint and patient (two per patient).
# We will mainly use the Immcantation output file here but the file needs to be updated with a column that contains the time point
LIM828246.cells <- LIM828246$cell_id
i <- 1
for (i in 1:length(LIM828246.cells)) {
  if(substr(LIM828246.cells[i],1,2)=="CC"){
    LIM828246.cells[i] <- paste0("Dataset_4_",substr(LIM828246.cells[i],3,nchar(LIM828246.cells[i])))
  }
  if(substr(LIM828246.cells[i],1,2)=="AA"){
    LIM828246.cells[i] <- paste0("Dataset_1_",substr(LIM828246.cells[i],3,nchar(LIM828246.cells[i])))
  }
  if(substr(LIM828246.cells[i],1,2)=="GG"){
    LIM828246.cells[i] <- paste0("Dataset_2_",substr(LIM828246.cells[i],3,nchar(LIM828246.cells[i])))
  }
  if(substr(LIM828246.cells[i],1,2)=="TT"){
    LIM828246.cells[i] <- paste0("Dataset_3_",substr(LIM828246.cells[i],3,nchar(LIM828246.cells[i])))
  }
}
timepoints <- Bmem.Inc@meta.data[rownames(Bmem.Inc@meta.data)%in%LIM828246.cells,c("Timepoint","Full.Row.names")]

i <- 1
for (i in 1:nrow(timepoints)) {
  if(substr(rownames(timepoints)[i],1,10)=="Dataset_4_"){
    rownames(timepoints)[i] <- paste0("CC",substr(rownames(timepoints)[i],11,nchar(rownames(timepoints)[i])))
  }
  if(substr(rownames(timepoints)[i],1,10)=="Dataset_1_"){
    rownames(timepoints)[i] <- paste0("AA",substr(rownames(timepoints)[i],11,nchar(rownames(timepoints)[i])))
  }
  if(substr(rownames(timepoints)[i],1,10)=="Dataset_2_"){
    rownames(timepoints)[i] <- paste0("GG",substr(rownames(timepoints)[i],11,nchar(rownames(timepoints)[i])))
  }
  if(substr(rownames(timepoints)[i],1,10)=="Dataset_3_"){
    rownames(timepoints)[i] <- paste0("TT",substr(rownames(timepoints)[i],11,nchar(rownames(timepoints)[i])))
  }
}

LIM828246.annotated <- as.data.frame(LIM828246.annotated)
rownames(LIM828246.annotated) <- LIM828246.annotated$cell_id
LIM828246.annotated <- merge(LIM828246.annotated, timepoints, by=0, all.x = T)


# Adding the colors that the pie chart segments will have later. Options are Singlet, Unique clone, Persistent clone and Persistent singlet
LIM828246.annotated$color <- NA
LIM828246.annotated$color[LIM828246.annotated$persistence=="Singlet"] <- "white"
LIM828246.annotated$color[LIM828246.annotated$persistence=="Unique clone"] <- "grey"
LIM828246.annotated$color[LIM828246.annotated$persistence=="Persistent singlet"] <- "white"

# For the category Persistent clone we need to give colors.
Persistents <- unique(LIM828246.annotated$clone_id[LIM828246.annotated$persistence=="Persistent clone"])
library(RColorBrewer)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col_vector <- col_vector[1:length(Persistents)]
df <- data.frame("clone_id"= Persistents, "color"=col_vector)

LIM828246.annotated <- merge(LIM828246.annotated, df, by="clone_id", all.x = T)
LIM828246.annotated$color <- coalesce(LIM828246.annotated$color.x,LIM828246.annotated$color.y)
LIM828246.annotated$color.x <- NULL
LIM828246.annotated$color.y <- NULL

LIM828246.annotated.6mo <- LIM828246.annotated[LIM828246.annotated$Timepoint=="6 months after infection",]
LIM828246.annotated.12mo <- LIM828246.annotated[LIM828246.annotated$Timepoint =="12 months after infection",]


# At this point we change to new code
# The persistence column in the dataframes at this point can have different entries: Persistent clone, Singlet, Unique clone
# The venn diagram should have overlapping regions showing persistent clones, singlets and unique clones just count as non-persistent

LIM828246.annotated.6mo$Venn_class <- "Non_persistent"
LIM828246.annotated.12mo$Venn_class <- "Non_persistent"
LIM828246.annotated.6mo$Venn_class[LIM828246.annotated.6mo$persistence=="Persistent clone"] <- "Persistent"
LIM828246.annotated.12mo$Venn_class[LIM828246.annotated.12mo$persistence=="Persistent clone"] <- "Persistent"

# We only need three numbers per patient: The number of non persistent clones per time point and the number of persistent clones.
LIM828246.annotated.6_nonpersistents <- length(unique(LIM828246.annotated.6mo$clone_id[LIM828246.annotated.6mo$Venn_class=="Non_persistent"]))
LIM828246.annotated.12_nonpersistents <- length(unique(LIM828246.annotated.12mo$clone_id[LIM828246.annotated.12mo$Venn_class=="Non_persistent"]))
LIM828246.annotated_persistents <- length(unique(LIM828246.annotated.6mo$clone_id[LIM828246.annotated.6mo$Venn_class=="Persistent"]))

# What follows here is for the patients where we don't have persistent clones
} else {
  
# Finding the timepoint for all cells of this patient's Immcantation output file

i <- 1
for(i in 1:nrow(LIM828246)){
    
    
      if(substr(LIM828246[i,"cell_id"],1,2)=="CC"){
        LIM828246[i,"cell_id"] <- paste0("Dataset_4_",substr(LIM828246[i,"cell_id"],3,nchar(LIM828246[i,"cell_id"])))
      }
  if(substr(LIM828246[i,"cell_id"],1,2)=="AA"){
    LIM828246[i,"cell_id"] <- paste0("Dataset_1_",substr(LIM828246[i,"cell_id"],3,nchar(LIM828246[i,"cell_id"])))
      }
  if(substr(LIM828246[i,"cell_id"],1,2)=="GG"){
    LIM828246[i,"cell_id"] <- paste0("Dataset_2_",substr(LIM828246[i,"cell_id"],3,nchar(LIM828246[i,"cell_id"])))
      }
  if(substr(LIM828246[i,"cell_id"],1,2)=="TT"){
    LIM828246[i,"cell_id"] <- paste0("Dataset_3_",substr(LIM828246[i,"cell_id"],3,nchar(LIM828246[i,"cell_id"])))
      }
    }
    

# Getting the timepoints
LIM828246$Full.Row.names <- LIM828246$cell_id
LIM828246 <- merge(LIM828246, Bmem.Inc@meta.data, by="Full.Row.names", all.x = T)
table(LIM828246$Timepoint)
 
LIM828246.annotated.6_nonpersistents <- nrow(LIM828246[LIM828246$Timepoint=="6 months after infection",])
LIM828246.annotated.12_nonpersistents<- nrow(LIM828246[LIM828246$Timepoint=="12 months after infection",])
LIM828246.annotated_persistents <- 0
}

output_df[z,1] <- Patient
output_df[z,2] <- LIM828246.annotated.6_nonpersistents
output_df[z,3] <- LIM828246.annotated.12_nonpersistents
output_df[z,4] <- LIM828246.annotated_persistents

}


# Creating the Venn diagramm using eulerr
# See documentation https://cran.r-project.org/web/packages/eulerr/vignettes/gallery.html
library(eulerr)

vd <- euler(c("6 months" = sum(output_df$n_non_persistents_6_months), "12 months" = sum(output_df$n_non_persistents_12_months), "6 months&12 months" = sum(output_df$n_persistents)))
plot(vd, counts = T,lwd = 2,
     fill=c("#004488","#997700"),
     opacity = .7,
     legend = list( space= "right", columns=1),
     quantities = list(type="counts"),
     main="Clonal overlap between timepoints") # Export 4x6

remove(Bmem.Inc,df,LIM828246,LIM828246.annotated,LIM828246.annotated.12mo,LIM828246.annotated.6mo,output_df,qual_col_pals,timepoints,vd,
       col_vector,i,LIM828246.annotated_persistents,LIM828246.annotated.12_nonpersistents,LIM828246.annotated.6_nonpersistents,LIM828246.cells,
       LIM828246.longitudinals,Patient,Patients,Persistents,z)


# Venn diagram showing the clonal overlap between cell subsets at the 12 months time point (all patients)
# Get all the clones and cells from the 12 month timepoint - use the Immcantation output where both timepoints have been merged for each patient
# In the end there has to be a data frame, which has all cells in rows, with columns for clone (make sure no overlapping clones for patients!) and cell subset. 
setwd("~/NAS/Jan/Experiments and Data/210801 - Exp 035 - Profiling of COVID vaccinated patient samples/Memory B cell Fourth Experiment 220323/multi_exp035_4")  
Bmem <- readRDS("Bmem.rds")
Patients <- c("PFCL1_LIM_313000","PFCL1_LIM828246","PFCL1_199474_1994","PFCL1_179308_1988","PFCL1_LIM_674950","PFCL1_UST_190762","PFCL1_LIM_137402",
              "PFCL1_154583_1943","PFCL1_196878_1691")

setwd("~/NAS/Jan/Experiments and Data/210801 - Exp 035 - Profiling of COVID vaccinated patient samples/Immcantation Analysis/Immcantation Outputs/update_020522")

z <- 1
totalfile <- as.data.frame(matrix(nrow = 1,ncol = 63))
for (z in 1:length(Patients)) {
  
  Patient <- Patients[z]
  Thisroundfile<- as.data.frame(read_airr(paste0("./",Patient,"_positives/outputs/",Patient,"_positives_0.20/",Patient,"_positives_0.20_heavy_germ-pass.tsv")))
  Thisroundfile$Patient <- Patient
  colnames(totalfile) <- colnames(Thisroundfile)
  totalfile <- rbind(totalfile,Thisroundfile)
}
totalfile <- totalfile[2:nrow(totalfile),]

# Translating the cell_id column
i <- 1
for(i in 1:nrow(totalfile)){
  
  
  if(substr(totalfile[i,"cell_id"],1,2)=="CC"){
    totalfile[i,"cell_id"] <- paste0("Dataset_4_",substr(totalfile[i,"cell_id"],3,nchar(totalfile[i,"cell_id"])))
  }
  if(substr(totalfile[i,"cell_id"],1,2)=="AA"){
    totalfile[i,"cell_id"] <- paste0("Dataset_1_",substr(totalfile[i,"cell_id"],3,nchar(totalfile[i,"cell_id"])))
  }
  if(substr(totalfile[i,"cell_id"],1,2)=="GG"){
    totalfile[i,"cell_id"] <- paste0("Dataset_2_",substr(totalfile[i,"cell_id"],3,nchar(totalfile[i,"cell_id"])))
  }
  if(substr(totalfile[i,"cell_id"],1,2)=="TT"){
    totalfile[i,"cell_id"] <- paste0("Dataset_3_",substr(totalfile[i,"cell_id"],3,nchar(totalfile[i,"cell_id"])))
  }
}

# Making the clone_ID unique so that no overlap between patient is possible
totalfile$clone_id <- paste0(totalfile$clone_id,totalfile$Patient)


# Get the cell subsets
totalfile$Full.Row.names <- totalfile$cell_id
totalfile <- merge(totalfile, Bmem@meta.data, by="Full.Row.names", all.x = T)
totalfile <- totalfile[totalfile$Timepoint=="12 months after infection",]
totalfile <- totalfile[,c("named.clusters","clone_id")]

# We should not have the same clone twice with the same cell state.
totalfile <- unique(totalfile)

# Make the Venn diagram 
# First a vector for each cell subset, containing the clones
Unswitched.vector <- totalfile[totalfile$named.clusters=="Unswitched",2]
CD27low_RM.vector <- totalfile[totalfile$named.clusters=="CD27low RM",2]
CD27high_RM.vector <- totalfile[totalfile$named.clusters=="CD27high RM",2]
Activated.vector <-totalfile[totalfile$named.clusters=="Activated",2]
Atypical.vector <- totalfile[totalfile$named.clusters=="Atypical",2]

# Plotting -> https://r-graph-gallery.com/14-venn-diagramm.html
setwd("~/NAS/Jan/Experiments and Data/210801 - Exp 035 - Profiling of COVID vaccinated patient samples/Memory B cell Fourth Experiment 220323/multi_exp035_4")
venn.diagram(
  x = list(Unswitched.vector, CD27low_RM.vector, CD27high_RM.vector,Activated.vector,Atypical.vector),
  category.names = c("Unswitched" , "CD27low RM" , "CD27high RM","Activated","Atypical"), filename = 'Clonal.overlap.between.subsets.12months.png', output=T,
  imagetype = "png", height = 3000, width = 3000, cat.pos = c(0,0,-115,115,0),cat.dist=c(0.2,0.2,0.25,0.2,0.2),fill=c("#DBC35E","#228833","#4477AA","#EE6677","#AA3377"))

# Only CD27hig RM vs Activated
venn.diagram(
  x = list(CD27high_RM.vector,Activated.vector),
  category.names = c("CD27high RM","Activated"), filename = 'Clonal.overlap.between.RM.and.Activated.12months.png', output=T,
  imagetype = "png", height = 3000, width = 3000,fill=c("#4477AA","#EE6677"),cat.pos=c(0,0))

# Only CD27high vs Atypical
venn.diagram(
  x = list(CD27high_RM.vector,Atypical.vector),
  category.names = c("CD27high RM","Atypical"), filename = 'Clonal.overlap.between.RM.and.Atypical.12months.png', output=T,
  imagetype = "png", height = 3000, width = 3000,fill=c("#4477AA","#AA3377"),cat.pos=c(0,0))


remove(Bmem, Thisroundfile,totalfile,Activated.vector,Atypical.vector,CD27high_RM.vector,CD27low_RM.vector,i,Patient,Patients,Unswitched.vector,z)

