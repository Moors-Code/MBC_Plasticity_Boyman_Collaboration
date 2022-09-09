# This script is used to prepare files which are used for clonal analysis in Immcantation.


###################################################################################################################################
# Part 1 - Preparation of .txt files (containing cell barcodes) for all patients individually and for all patients together
###################################################################################################################################
# There will be one vector for each patient. 
# The vectors will contain barcodes from both timepoints, and will represent bait positive cells.

# The Seurat objects after preprocessing is loaded.
# Also, the pre processed dataset containing Set 4 alone, and the list of the naive cells of Set 4 are loaded.
setwd("~/NAS/Jan/Experiments and Data/210801 - Exp 035 - Profiling of COVID vaccinated patient samples/Memory B cell Fourth Experiment 220323/multi_exp035_4")  
Bmem.Inc <- readRDS("Bmem.rds")
Bmem_fourth <- readRDS("Bmem_fourth.IgD.rds")
Naives <- (read.table("Naives.txt"))
Naives <- Naives$V1


# This part removes the Dataset_1 part (int the Row.names column!) and instead adds the following code to the cell barcodes:
# added to the start of the barcode:
# AA -> dataset1
# GG -> dataset2
# TT -> dataset3
# CC -> dataset4
i <-1
for (i in 1:nrow(Bmem.Inc@meta.data)) {
  if(substr(rownames(Bmem.Inc@meta.data)[i],1,9)=="Dataset_1"){
    Bmem.Inc@meta.data$Row.names[i] <- paste0("AA",Bmem.Inc@meta.data$Row.names[i])
  }
  if(substr(rownames(Bmem.Inc@meta.data)[i],1,9)=="Dataset_2"){
    Bmem.Inc@meta.data$Row.names[i] <- paste0("GG",Bmem.Inc@meta.data$Row.names[i])
  }
  if(substr(rownames(Bmem.Inc@meta.data)[i],1,9)=="Dataset_3"){
    Bmem.Inc@meta.data$Row.names[i] <- paste0("TT",Bmem.Inc@meta.data$Row.names[i])
  }
  if(substr(rownames(Bmem.Inc@meta.data)[i],1,9)=="Dataset_4"){
    Bmem.Inc@meta.data$Row.names[i] <- paste0("CC",Bmem.Inc@meta.data$Row.names[i])
  }
}
i <- 1
for (i in 1:nrow(Bmem_fourth@meta.data)) {
  
  if(substr(rownames(Bmem_fourth@meta.data)[i],1,9)=="Dataset_4"){
    Bmem_fourth@meta.data$Row.names[i] <- paste0("CC",Bmem_fourth@meta.data$Row.names[i])
  }
}
i <- 1
for (i in 1:length(Naives)) {
  Naives[i] <- paste0("CC",substr(Naives[i],11,nchar(Naives[i])))
}

# Make a vector with all patient names
Patients <- unique(Bmem.Inc@meta.data$Patient)

# Make the vectors
PFCL1_LIM_313000_positives <- Bmem.Inc@meta.data$Row.names[Bmem.Inc@meta.data$Patient == Patients[1] &
                                                             Bmem.Inc@meta.data$bait.positive=="yes"] 
PFCL1_UST_190762_positives <- Bmem.Inc@meta.data$Row.names[Bmem.Inc@meta.data$Patient == Patients[2] &
                                                             Bmem.Inc@meta.data$bait.positive=="yes"]
PFCL1_LIM_137402_positives <- Bmem.Inc@meta.data$Row.names[Bmem.Inc@meta.data$Patient == Patients[3] &
                                                             Bmem.Inc@meta.data$bait.positive=="yes"]
PFCL1_LIM_674950_positives <- Bmem.Inc@meta.data$Row.names[Bmem.Inc@meta.data$Patient == Patients[4] &
                                                             Bmem.Inc@meta.data$bait.positive=="yes"]
PFCL1_179308_1988_positives <- Bmem.Inc@meta.data$Row.names[Bmem.Inc@meta.data$Patient == Patients[8] &
                                                              Bmem.Inc@meta.data$bait.positive=="yes"]
PFCL1_199474_1994_positives <- Bmem.Inc@meta.data$Row.names[Bmem.Inc@meta.data$Patient == Patients[5] &
                                                              Bmem.Inc@meta.data$bait.positive=="yes"]
PFCL1_196878_1691_positives <- Bmem.Inc@meta.data$Row.names[Bmem.Inc@meta.data$Patient == Patients[9] &
                                                              Bmem.Inc@meta.data$bait.positive=="yes"]
PFCL1_LIM828246_positives <- Bmem.Inc@meta.data$Row.names[Bmem.Inc@meta.data$Patient == Patients[6] &
                                                            Bmem.Inc@meta.data$bait.positive=="yes"]
PFCL1_154583_1943_positives <- Bmem.Inc@meta.data$Row.names[Bmem.Inc@meta.data$Patient == Patients[7] &
                                                              Bmem.Inc@meta.data$bait.positive=="yes"]
All.cells <- Bmem.Inc@meta.data$Row.names
Bmem_fourth <- Bmem_fourth@meta.data$Row.names
Bmem_fourth <- Bmem_fourth[Bmem_fourth%in%Naives]
All.cells <- c(All.cells, Bmem_fourth) # This all cells object then contains all the preprocessed cells - including the naive cells. Used for mutational load analysis.

write(PFCL1_LIM_313000_positives, file = "~/NAS/Jan/Experiments and Data/210801 - Exp 035 - Profiling of COVID vaccinated patient samples/Immcantation Analysis/Immcantation Inputs/PFCL1_LIM_313000_positives.txt")
write(PFCL1_UST_190762_positives, file = "~/NAS/Jan/Experiments and Data/210801 - Exp 035 - Profiling of COVID vaccinated patient samples/Immcantation Analysis/Immcantation Inputs/PFCL1_UST_190762_positives.txt")
write(PFCL1_LIM_137402_positives, file = "~/NAS/Jan/Experiments and Data/210801 - Exp 035 - Profiling of COVID vaccinated patient samples/Immcantation Analysis/Immcantation Inputs/PFCL1_LIM_137402_positives.txt")
write(PFCL1_LIM_674950_positives, file = "~/NAS/Jan/Experiments and Data/210801 - Exp 035 - Profiling of COVID vaccinated patient samples/Immcantation Analysis/Immcantation Inputs/PFCL1_LIM_674950_positives.txt")
write(PFCL1_179308_1988_positives, file = "~/NAS/Jan/Experiments and Data/210801 - Exp 035 - Profiling of COVID vaccinated patient samples/Immcantation Analysis/Immcantation Inputs/PFCL1_179308_1988_positives.txt")
write(PFCL1_199474_1994_positives, file = "~/NAS/Jan/Experiments and Data/210801 - Exp 035 - Profiling of COVID vaccinated patient samples/Immcantation Analysis/Immcantation Inputs/PFCL1_199474_1994_positives.txt")
write(PFCL1_196878_1691_positives, file = "~/NAS/Jan/Experiments and Data/210801 - Exp 035 - Profiling of COVID vaccinated patient samples/Immcantation Analysis/Immcantation Inputs/PFCL1_196878_1691_positives.txt")
write(PFCL1_LIM828246_positives, file = "~/NAS/Jan/Experiments and Data/210801 - Exp 035 - Profiling of COVID vaccinated patient samples/Immcantation Analysis/Immcantation Inputs/PFCL1_LIM828246_positives.txt")
write(PFCL1_154583_1943_positives, file = "~/NAS/Jan/Experiments and Data/210801 - Exp 035 - Profiling of COVID vaccinated patient samples/Immcantation Analysis/Immcantation Inputs/PFCL1_154583_1943_positives.txt")
write(All.cells, file = "~/NAS/Jan/Experiments and Data/210801 - Exp 035 - Profiling of COVID vaccinated patient samples/Immcantation Analysis/Immcantation Inputs/All.cells.txt")


###################################################################################################################################
# Part 2 - Preparation of .csv files (containing BCR contigs) for all patients individually and for all patients together
###################################################################################################################################
# Load the four filtered_contig_annotations, adapt the barcodes, merge the four objects, subset the objects to generate one per patient
# Loading data
contigs_Dataset1 <- read.csv("~/NAS/Jan/Experiments and Data/210801 - Exp 035 - Profiling of COVID vaccinated patient samples/Memory B cell Pilot Data 211010/multi_SpikeMemoryBCells/outs/per_sample_outs/multi_SpikeMemoryBCells/vdj_b/filtered_contig_annotations.csv")
contigs_Dataset2 <- read.csv("~/NAS/Jan/Experiments and Data/210801 - Exp 035 - Profiling of COVID vaccinated patient samples/Memory B cell Second Experiment 211127/multi_exp035_MemoryBcellSecondExperiment211127/outs/per_sample_outs/multi_exp035_MemoryBcellSecondExperiment211127/vdj_b/filtered_contig_annotations.csv")
contigs_Dataset3 <- read.csv("~/NAS/Jan/Experiments and Data/210801 - Exp 035 - Profiling of COVID vaccinated patient samples/Memory B cell Third Experiment 220214/multi_exp035_3/outs/per_sample_outs/multi_exp035_3/vdj_b/filtered_contig_annotations.csv")
contigs_Dataset4 <- read.csv("~/NAS/Jan/Experiments and Data/210801 - Exp 035 - Profiling of COVID vaccinated patient samples/Memory B cell Fourth Experiment 220323/multi_exp035_4/outs/per_sample_outs/multi_exp035_4/vdj_b/filtered_contig_annotations.csv")

# Changing barcodes
contigs_Dataset1$barcode <- paste0("AA",contigs_Dataset1$barcode)
contigs_Dataset2$barcode <- paste0("GG",contigs_Dataset2$barcode)
contigs_Dataset3$barcode <- paste0("TT",contigs_Dataset3$barcode)
contigs_Dataset4$barcode <- paste0("CC",contigs_Dataset4$barcode)

# Changing contig_id
contigs_Dataset1$contig_id <- paste0("AA",contigs_Dataset1$contig_id)
contigs_Dataset2$contig_id <- paste0("GG",contigs_Dataset2$contig_id)
contigs_Dataset3$contig_id <- paste0("TT",contigs_Dataset3$contig_id)
contigs_Dataset4$contig_id <- paste0("CC",contigs_Dataset4$contig_id)


# Merging
all.contigs <- rbind(contigs_Dataset1,contigs_Dataset2,contigs_Dataset3,contigs_Dataset4)
remove(contigs_Dataset1,contigs_Dataset2,contigs_Dataset3,contigs_Dataset4)

# Creating patient specific contig objects
contigs_PFCL1_LIM_313000_positives <- all.contigs[all.contigs$barcode %in% PFCL1_LIM_313000_positives,]
contigs_PFCL1_UST_190762_positives <- all.contigs[all.contigs$barcode %in% PFCL1_UST_190762_positives,]
contigs_PFCL1_LIM_137402_positives <- all.contigs[all.contigs$barcode %in% PFCL1_LIM_137402_positives,]
contigs_PFCL1_LIM_674950_positives <- all.contigs[all.contigs$barcode %in% PFCL1_LIM_674950_positives,]
contigs_PFCL1_179308_1988_positives <- all.contigs[all.contigs$barcode %in% PFCL1_179308_1988_positives,]
contigs_PFCL1_199474_1994_positives <- all.contigs[all.contigs$barcode %in% PFCL1_199474_1994_positives,]
contigs_PFCL1_196878_1691_positives <- all.contigs[all.contigs$barcode %in% PFCL1_196878_1691_positives,]
contigs_PFCL1_LIM828246_positives <- all.contigs[all.contigs$barcode %in% PFCL1_LIM828246_positives,]
contigs_PFCL1_154583_1943_positives <- all.contigs[all.contigs$barcode %in% PFCL1_154583_1943_positives,]
contigs_all.cells <- all.contigs[all.contigs$barcode %in% All.cells,]

# Write the files
write_csv(contigs_PFCL1_LIM_313000_positives, file = "~/NAS/Jan/Experiments and Data/210801 - Exp 035 - Profiling of COVID vaccinated patient samples/Immcantation Analysis/Immcantation Inputs/filtered_contig_annotations_PFCL1_LIM_313000_positives.csv")
write_csv(contigs_PFCL1_UST_190762_positives, file = "~/NAS/Jan/Experiments and Data/210801 - Exp 035 - Profiling of COVID vaccinated patient samples/Immcantation Analysis/Immcantation Inputs/filtered_contig_annotations_PFCL1_UST_190762_positives.csv")
write_csv(contigs_PFCL1_LIM_137402_positives, file = "~/NAS/Jan/Experiments and Data/210801 - Exp 035 - Profiling of COVID vaccinated patient samples/Immcantation Analysis/Immcantation Inputs/filtered_contig_annotations_PFCL1_LIM_137402_positives.csv")
write_csv(contigs_PFCL1_LIM_674950_positives, file = "~/NAS/Jan/Experiments and Data/210801 - Exp 035 - Profiling of COVID vaccinated patient samples/Immcantation Analysis/Immcantation Inputs/filtered_contig_annotations_PFCL1_LIM_674950_positives.csv")
write_csv(contigs_PFCL1_179308_1988_positives, file = "~/NAS/Jan/Experiments and Data/210801 - Exp 035 - Profiling of COVID vaccinated patient samples/Immcantation Analysis/Immcantation Inputs/filtered_contig_annotations_PFCL1_179308_1988_positives.csv")
write_csv(contigs_PFCL1_199474_1994_positives, file = "~/NAS/Jan/Experiments and Data/210801 - Exp 035 - Profiling of COVID vaccinated patient samples/Immcantation Analysis/Immcantation Inputs/filtered_contig_annotations_PFCL1_199474_1994_positives.csv")
write_csv(contigs_PFCL1_196878_1691_positives, file = "~/NAS/Jan/Experiments and Data/210801 - Exp 035 - Profiling of COVID vaccinated patient samples/Immcantation Analysis/Immcantation Inputs/filtered_contig_annotations_PFCL1_196878_1691_positives.csv")
write_csv(contigs_PFCL1_LIM828246_positives, file = "~/NAS/Jan/Experiments and Data/210801 - Exp 035 - Profiling of COVID vaccinated patient samples/Immcantation Analysis/Immcantation Inputs/filtered_contig_annotations_PFCL1_LIM828246_positives.csv")
write_csv(contigs_PFCL1_154583_1943_positives, file = "~/NAS/Jan/Experiments and Data/210801 - Exp 035 - Profiling of COVID vaccinated patient samples/Immcantation Analysis/Immcantation Inputs/filtered_contig_annotations_PFCL1_154583_1943_positives.csv")
write_csv(contigs_all.cells, file = "~/NAS/Jan/Experiments and Data/210801 - Exp 035 - Profiling of COVID vaccinated patient samples/Immcantation Analysis/Immcantation Inputs/contigs_all.cells.csv")

remove(PFCL1_LIM_313000_positives,PFCL1_UST_190762_positives, PFCL1_LIM_137402_positives,PFCL1_LIM_674950_positives,
       PFCL1_179308_1988_positives,PFCL1_199474_1994_positives,PFCL1_196878_1691_positives,PFCL1_154583_1943_positives,
       PFCL1_LIM828246_positives, All.cells, Patients, i, all.contigs, Bmem.Inc, Bmem_fourth,Naives,
       contigs_all.cells,contigs_PFCL1_LIM_313000_positives,contigs_PFCL1_UST_190762_positives,contigs_PFCL1_LIM_137402_positives,
       contigs_PFCL1_LIM_674950_positives,contigs_PFCL1_179308_1988_positives,contigs_PFCL1_199474_1994_positives,contigs_PFCL1_196878_1691_positives,
       contigs_PFCL1_LIM828246_positives,contigs_PFCL1_154583_1943_positives)


###################################################################################################################################
# Part 3 - More automated preparation of files for selected groups of cells
###################################################################################################################################
# Only create a vector of cell barcodes, give this as input to a function which gives the .txt and the .csv as outputs with a specified name
setwd("~/NAS/Jan/Experiments and Data/210801 - Exp 035 - Profiling of COVID vaccinated patient samples/Memory B cell Fourth Experiment 220323/multi_exp035_4")  
Bmem.Inc <- readRDS("Bmem.rds")

# Create the vector of cell barcodes - this vector needs to be updated for the respective group of cells to analyze.
PFCL1_179308_1988_vacc_late <- rownames(Bmem.Inc@meta.data[Bmem.Inc@meta.data$bait.positive=="yes" &
                                           Bmem.Inc@meta.data$Patient=="PFCL1-179308-1988"&
                                             Bmem.Inc@meta.data$time_after_vac=="85d",])
PFCL1_LIM828246_vacc_late <- rownames(Bmem.Inc@meta.data[Bmem.Inc@meta.data$bait.positive=="yes" &
                                                             Bmem.Inc@meta.data$Patient=="PFCL1-LIM828246"&
                                                             Bmem.Inc@meta.data$time_after_vac=="87d",])
PFCL1_199474_1994_vacc_late <- rownames(Bmem.Inc@meta.data[Bmem.Inc@meta.data$bait.positive=="yes" &
                                                           Bmem.Inc@meta.data$Patient=="PFCL1-199474-1994"&
                                                           Bmem.Inc@meta.data$time_after_vac=="108d",])
PFCL1_LIM_137402_vacc_early <- rownames(Bmem.Inc@meta.data[Bmem.Inc@meta.data$bait.positive=="yes" &
                                                             Bmem.Inc@meta.data$Patient=="PFCL1-LIM-137402"&
                                                             Bmem.Inc@meta.data$time_after_vac=="15d",])
PFCL1_LIM_313000_vacc_early <- rownames(Bmem.Inc@meta.data[Bmem.Inc@meta.data$bait.positive=="yes" &
                                                             Bmem.Inc@meta.data$Patient=="PFCL1-LIM-313000"&
                                                             Bmem.Inc@meta.data$time_after_vac=="15d",])
PFCL1_154583_1943_vacc_early <- rownames(Bmem.Inc@meta.data[Bmem.Inc@meta.data$bait.positive=="yes" &
                                                             Bmem.Inc@meta.data$Patient=="PFCL1-154583-1943"&
                                                             Bmem.Inc@meta.data$time_after_vac=="23d",])
PFCL1_179308_1988_6mo <- rownames(Bmem.Inc@meta.data[Bmem.Inc@meta.data$bait.positive=="yes" &
                                                             Bmem.Inc@meta.data$Patient=="PFCL1-179308-1988"&
                                                             Bmem.Inc@meta.data$Timepoint=="6 months after infection",])
PFCL1_LIM828246_6mo <- rownames(Bmem.Inc@meta.data[Bmem.Inc@meta.data$bait.positive=="yes" &
                                                           Bmem.Inc@meta.data$Patient=="PFCL1-LIM828246"&
                                                           Bmem.Inc@meta.data$Timepoint=="6 months after infection",])
PFCL1_199474_1994_6mo <- rownames(Bmem.Inc@meta.data[Bmem.Inc@meta.data$bait.positive=="yes" &
                                                             Bmem.Inc@meta.data$Patient=="PFCL1-199474-1994"&
                                                             Bmem.Inc@meta.data$Timepoint=="6 months after infection",])
PFCL1_LIM_137402_6mo <- rownames(Bmem.Inc@meta.data[Bmem.Inc@meta.data$bait.positive=="yes" &
                                                             Bmem.Inc@meta.data$Patient=="PFCL1-LIM-137402"&
                                                             Bmem.Inc@meta.data$Timepoint=="6 months after infection",])
PFCL1_LIM_313000_6mo <- rownames(Bmem.Inc@meta.data[Bmem.Inc@meta.data$bait.positive=="yes" &
                                                             Bmem.Inc@meta.data$Patient=="PFCL1-LIM-313000"&
                                                             Bmem.Inc@meta.data$Timepoint=="6 months after infection",])
PFCL1_154583_1943_6mo <- rownames(Bmem.Inc@meta.data[Bmem.Inc@meta.data$bait.positive=="yes" &
                                                              Bmem.Inc@meta.data$Patient=="PFCL1-154583-1943"&
                                                              Bmem.Inc@meta.data$Timepoint=="6 months after infection",])
# The non vaccinated ones
PFCL1_196878_1691_6mo <- rownames(Bmem.Inc@meta.data[Bmem.Inc@meta.data$bait.positive=="yes" &
                                                      Bmem.Inc@meta.data$Patient=="PFCL1-196878-1691"&
                                                      Bmem.Inc@meta.data$Timepoint=="6 months after infection",])

PFCL1_LIM_674950_6mo <- rownames(Bmem.Inc@meta.data[Bmem.Inc@meta.data$bait.positive=="yes" &
                                                     Bmem.Inc@meta.data$Patient=="PFCL1-LIM-674950"&
                                                     Bmem.Inc@meta.data$Timepoint=="6 months after infection",])

PFCL1_UST_190762_6mo <- rownames(Bmem.Inc@meta.data[Bmem.Inc@meta.data$bait.positive=="yes" &
                                                     Bmem.Inc@meta.data$Patient=="PFCL1-UST-190762"&
                                                     Bmem.Inc@meta.data$Timepoint=="6 months after infection",])


PFCL1_196878_1691_12mo <- rownames(Bmem.Inc@meta.data[Bmem.Inc@meta.data$bait.positive=="yes" &
                                                       Bmem.Inc@meta.data$Patient=="PFCL1-196878-1691"&
                                                       Bmem.Inc@meta.data$Timepoint=="12 months after infection",])

PFCL1_LIM_674950_12mo <- rownames(Bmem.Inc@meta.data[Bmem.Inc@meta.data$bait.positive=="yes" &
                                                      Bmem.Inc@meta.data$Patient=="PFCL1-LIM-674950"&
                                                      Bmem.Inc@meta.data$Timepoint=="12 months after infection",])

PFCL1_UST_190762_12mo <- rownames(Bmem.Inc@meta.data[Bmem.Inc@meta.data$bait.positive=="yes" &
                                                      Bmem.Inc@meta.data$Patient=="PFCL1-UST-190762"&
                                                      Bmem.Inc@meta.data$Timepoint=="12 months after infection",])



# This function can be used for a more automated preparation of Immcantation input files.
# The function needs as input a vector of row names from the Seurat Object Bmem.Inc. Also, the file name of the exported files has to be defined.
Immcantation_Input_function <- function(x,name_of_file){
  Cell.list <- x
  i <-1
  for (i in 1:length(Cell.list)) {
    if(substr(Cell.list[i],1,9)=="Dataset_1"){
      Cell.list[i] <- paste0("AA",substr(Cell.list[i],11,nchar(Cell.list[i])))
    }
    if(substr(Cell.list[i],1,9)=="Dataset_2"){
      Cell.list[i] <- paste0("GG",substr(Cell.list[i],11,nchar(Cell.list[i])))
    }
    if(substr(Cell.list[i],1,9)=="Dataset_3"){
      Cell.list[i] <- paste0("TT",substr(Cell.list[i],11,nchar(Cell.list[i])))
    }
    if(substr(Cell.list[i],1,9)=="Dataset_4"){
      Cell.list[i] <- paste0("CC",substr(Cell.list[i],11,nchar(Cell.list[i])))
    }
  }
  
  x <- Cell.list
  write(x, file = paste0("~/NAS/Jan/Experiments and Data/210801 - Exp 035 - Profiling of COVID vaccinated patient samples/Immcantation Analysis/Immcantation Inputs/",name_of_file,".txt"))
  
  # Loading data
  contigs_Dataset1 <- read.csv("~/NAS/Jan/Experiments and Data/210801 - Exp 035 - Profiling of COVID vaccinated patient samples/Memory B cell Pilot Data 211010/multi_SpikeMemoryBCells/outs/per_sample_outs/multi_SpikeMemoryBCells/vdj_b/filtered_contig_annotations.csv")
  contigs_Dataset2 <- read.csv("~/NAS/Jan/Experiments and Data/210801 - Exp 035 - Profiling of COVID vaccinated patient samples/Memory B cell Second Experiment 211127/multi_exp035_MemoryBcellSecondExperiment211127/outs/per_sample_outs/multi_exp035_MemoryBcellSecondExperiment211127/vdj_b/filtered_contig_annotations.csv")
  contigs_Dataset3 <- read.csv("~/NAS/Jan/Experiments and Data/210801 - Exp 035 - Profiling of COVID vaccinated patient samples/Memory B cell Third Experiment 220214/multi_exp035_3/outs/per_sample_outs/multi_exp035_3/vdj_b/filtered_contig_annotations.csv")
  contigs_Dataset4 <- read.csv("~/NAS/Jan/Experiments and Data/210801 - Exp 035 - Profiling of COVID vaccinated patient samples/Memory B cell Fourth Experiment 220323/multi_exp035_4/outs/per_sample_outs/multi_exp035_4/vdj_b/filtered_contig_annotations.csv")
  
  # Changing barcodes
  contigs_Dataset1$barcode <- paste0("AA",contigs_Dataset1$barcode)
  contigs_Dataset2$barcode <- paste0("GG",contigs_Dataset2$barcode)
  contigs_Dataset3$barcode <- paste0("TT",contigs_Dataset3$barcode)
  contigs_Dataset4$barcode <- paste0("CC",contigs_Dataset4$barcode)
  
  # Changing contig_id
  contigs_Dataset1$contig_id <- paste0("AA",contigs_Dataset1$contig_id)
  contigs_Dataset2$contig_id <- paste0("GG",contigs_Dataset2$contig_id)
  contigs_Dataset3$contig_id <- paste0("TT",contigs_Dataset3$contig_id)
  contigs_Dataset4$contig_id <- paste0("CC",contigs_Dataset4$contig_id)
  
  # Merging
  all.contigs <- rbind(contigs_Dataset1,contigs_Dataset2,contigs_Dataset3,contigs_Dataset4)
  remove(contigs_Dataset1,contigs_Dataset2,contigs_Dataset3,contigs_Dataset4)
  
  contigs_x <- all.contigs[all.contigs$barcode %in% x,]
  write_csv(contigs_x, file = paste0("~/NAS/Jan/Experiments and Data/210801 - Exp 035 - Profiling of COVID vaccinated patient samples/Immcantation Analysis/Immcantation Inputs/",name_of_file,".csv"))
  
} 

Immcantation_Input_function(PFCL1_UST_190762_12mo,"PFCL1_UST_190762_12mo")

remove(Cell.list)
