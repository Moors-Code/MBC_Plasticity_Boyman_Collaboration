# This is the GSVA analysis of the MBC project. (Fig.5.E)
# A subsetted Seurat object("Bmem.subset.rds") was prepared in part 3 of the main analysis script of the project.
# We are running the GSVA in a separate script only because we could not install all the packages needed for it in the R project which we were using
# for the other analyses.

# Loading packages needed
library("GSVA")
library("GSVAdata")
library("openxlsx")
library("limma")
library("SingleCellExperiment")
library("Seurat")
library("scater")
library("Biobase")
library("viridis")
library("ComplexHeatmap")
library('circlize')
library('GOexpress')
library('MKmisc')

# Clearing the environment
rm(list = ls())
setwd("~/NAS/Jan/Experiments and Data/210801 - Exp 035 - Profiling of COVID vaccinated patient samples/GSVA")

# Loading in the dataset I want to use for GSVA
Bmem.subset <- readRDS("./Bmem.subset.rds")
Bmem.subset <- subset(Bmem.subset, Dataset != "Second")

# Preparing the dataset for GSVA -> Defining which groups to compare and calculation of pseudobulk expression values.
Bmem.subset$named.clusters[Bmem.subset@meta.data$named.clusters=="CD27high RM"|
                             Bmem.subset@meta.data$named.clusters=="CD27low RM"] <- "RM"
Idents(Bmem.subset) <- "named.clusters"

# Using single cells as individual samples for gsva analysis
# The first step is to create an Expression set object, which contains the data to be analyzed.
# I followed the vignette from Michael Love: https://biodatascience.github.io/compbio/bioc/objects.html
input_matrix <- Bmem.subset@assays$RNA@data[,colnames(Bmem.subset@assays$RNA@data) %in% Bmem.subset@meta.data$Full.Row.names]

# Changing column order so that later, when the heatmap is made, the cells of the same subset are next to each other
col.order <- c(rownames(Bmem.subset@meta.data[Bmem.subset@meta.data$named.clusters=="RM",]),
               rownames(Bmem.subset@meta.data[Bmem.subset@meta.data$named.clusters=="Activated",]),
               rownames(Bmem.subset@meta.data[Bmem.subset@meta.data$named.clusters=="Atypical",]))
input_matrix <- input_matrix[,col.order]

my_exprs <- matrix(input_matrix,ncol(input_matrix),nrow=nrow(input_matrix))

my_phenoData <- data.frame(sample=colnames(Bmem.subset@assays$RNA@data[,colnames(Bmem.subset@assays$RNA@data) %in% Bmem.subset@meta.data$Full.Row.names]),
                           condition="RM",treated="no")

Bmem.subset@meta.data$Rownumber <- 1:nrow(Bmem.subset@meta.data)
Activated.positions <- Bmem.subset@meta.data$Rownumber[Bmem.subset@meta.data$named.clusters=="Activated"]
Atypical.positions <- Bmem.subset@meta.data$Rownumber[Bmem.subset@meta.data$named.clusters=="Atypical"]

my_phenoData$condition[Activated.positions] <- "Activated"
my_phenoData$condition[Atypical.positions] <- "Atypical"

# Also here, correct the order of the cells, this time the rows.
my_phenoData$condition <- factor(my_phenoData$condition, levels = c("RM","Activated","Atypical"))
my_phenoData <- my_phenoData[order(my_phenoData$condition),] 
rownames(my_phenoData) <- c(1:nrow(my_phenoData))

my_feature_Data <- data.frame(geneID=rownames(Bmem.subset@assays$RNA@data),geneSymbol=rownames(Bmem.subset@assays$RNA@data))

my_eset <- ExpressionSet(my_exprs,
                         AnnotatedDataFrame(my_phenoData),
                         AnnotatedDataFrame(my_feature_Data))
rownames(my_eset) <- rownames(Bmem.subset@assays$RNA@data)
colnames(my_eset) <- col.order

my_eset
exprs(my_eset)
pData(my_eset)
fData(my_eset)

# After the construction of the ExpressionSet Object, I load the gene sets that we are interested in as .xml files.
gs1 <- getBroadSets(c("./Gene Sets/GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_PEPTIDE_ANTIGEN.v7.5.1.xml",
                      "./Gene Sets/HALLMARK_INTERFERON_GAMMA_RESPONSE.v7.5.1.xml",
                      "./Gene Sets/GOBP_REGULATION_OF_LYMPHOCYTE_ACTIVATION.v7.5.1.xml",
                      "./Gene Sets/GOBP_INTEGRIN_MEDIATED_SIGNALING_PATHWAY.v7.5.1.xml",
                      "./Gene Sets/GOBP_B_CELL_PROLIFERATION.v7.5.1.xml",
                      "./Gene Sets/GOBP_B_CELL_RECEPTOR_SIGNALING_PATHWAY.v7.5.1.xml"))

# Adding the atypical gene set
gs1@.Data[7]<-gs1@.Data[1]
object <- gs1@.Data[7]
atypical_genes <- read.xlsx("./Gene Sets/elife-41641-supp2-v2.xlsx")
object[[1]]@geneIds <- atypical_genes$Sanz_DN2_Up
object[[1]]@setName <- "Atypical_Gene_Set"
object[[1]]@setIdentifier <- "JM1"
gs1@.Data[7] <- object

# Alternative for getting all the gene sets from the hallmark collection
#gs2 <- getGmt("./Gene Sets/h.all.v7.5.1.symbols.gmt")

# Then, the gsva is done. Here, I am following the vignette and explanations from:
# https://github.com/rcastelo/GSVA/blob/master/vignettes/GSVA.Rmd
# https://www.youtube.com/watch?v=ZRet1oeGiUU&t=665s
# https://www.youtube.com/watch?v=Hg1abiNlPE4&t=4s
my_eset <- gsva(my_eset,gs1,kcdf="Gaussian")
head(featureNames(my_eset))


# From here on we are creating pseudo bulk enrichment scores instead of continuing with the single cell resolution.
# The reason for this is, that it is not really fair to use single cells as replicates.
# Due to factors in our data, like for example patients (cells from same patient have high tendency to be similar to each other), 
# the variance across the single cells is much smaller than what would be expected if they were real samples.
# This small variance leads to the phenomenon that small differences between groups are called to be significant - which makes no sense.
# By reducing the resolution to patient pseudo bulks, I am reducing the number of samples greatly.
# The input column into limma will then have only 8 (patients) x 3 (cell subsets) = 24 columns.
# See all the commented lines - this is done because for pseudobulking, we only include patients with n>20 cells in each cell subset.

# For this I will just create a new ExpressionSet Object because this is easier for me than changing the old one.
# Getting the cell names of the 24 groups
PFCL1_LIM_313000_RM <- Bmem.subset@meta.data$Full.Row.names[Bmem.subset@meta.data$Patient=="PFCL1-LIM-313000"& Bmem.subset@meta.data$named.clusters=="RM"]
#PFCL1_UST_190762_RM <- Bmem.subset@meta.data$Full.Row.names[Bmem.subset@meta.data$Patient=="PFCL1-UST-190762"& Bmem.subset@meta.data$named.clusters=="RM"]
PFCL1_LIM_137402_RM <- Bmem.subset@meta.data$Full.Row.names[Bmem.subset@meta.data$Patient=="PFCL1-LIM-137402"& Bmem.subset@meta.data$named.clusters=="RM"]
#PFCL1_LIM_674950_RM <- Bmem.subset@meta.data$Full.Row.names[Bmem.subset@meta.data$Patient=="PFCL1-LIM-674950"& Bmem.subset@meta.data$named.clusters=="RM"]
PFCL1_179308_1988_RM <- Bmem.subset@meta.data$Full.Row.names[Bmem.subset@meta.data$Patient=="PFCL1-179308-1988"& Bmem.subset@meta.data$named.clusters=="RM"]
#PFCL1_196878_1691_RM <- Bmem.subset@meta.data$Full.Row.names[Bmem.subset@meta.data$Patient=="PFCL1-196878-1691"& Bmem.subset@meta.data$named.clusters=="RM"]
PFCL1_199474_1994_RM <- Bmem.subset@meta.data$Full.Row.names[Bmem.subset@meta.data$Patient=="PFCL1-199474-1994"& Bmem.subset@meta.data$named.clusters=="RM"]
PFCL1_LIM828246_RM <- Bmem.subset@meta.data$Full.Row.names[Bmem.subset@meta.data$Patient=="PFCL1-LIM828246"& Bmem.subset@meta.data$named.clusters=="RM"]

PFCL1_LIM_313000_ACT <- Bmem.subset@meta.data$Full.Row.names[Bmem.subset@meta.data$Patient=="PFCL1-LIM-313000"& Bmem.subset@meta.data$named.clusters=="Activated"]
#PFCL1_UST_190762_ACT <- Bmem.subset@meta.data$Full.Row.names[Bmem.subset@meta.data$Patient=="PFCL1-UST-190762"& Bmem.subset@meta.data$named.clusters=="Activated"]
PFCL1_LIM_137402_ACT <- Bmem.subset@meta.data$Full.Row.names[Bmem.subset@meta.data$Patient=="PFCL1-LIM-137402"& Bmem.subset@meta.data$named.clusters=="Activated"]
#PFCL1_LIM_674950_ACT <- Bmem.subset@meta.data$Full.Row.names[Bmem.subset@meta.data$Patient=="PFCL1-LIM-674950"& Bmem.subset@meta.data$named.clusters=="Activated"]
PFCL1_179308_1988_ACT <- Bmem.subset@meta.data$Full.Row.names[Bmem.subset@meta.data$Patient=="PFCL1-179308-1988"& Bmem.subset@meta.data$named.clusters=="Activated"]
#PFCL1_196878_1691_ACT <- Bmem.subset@meta.data$Full.Row.names[Bmem.subset@meta.data$Patient=="PFCL1-196878-1691"& Bmem.subset@meta.data$named.clusters=="Activated"]
PFCL1_199474_1994_ACT <- Bmem.subset@meta.data$Full.Row.names[Bmem.subset@meta.data$Patient=="PFCL1-199474-1994"& Bmem.subset@meta.data$named.clusters=="Activated"]
PFCL1_LIM828246_ACT <- Bmem.subset@meta.data$Full.Row.names[Bmem.subset@meta.data$Patient=="PFCL1-LIM828246"& Bmem.subset@meta.data$named.clusters=="Activated"]

PFCL1_LIM_313000_AT <- Bmem.subset@meta.data$Full.Row.names[Bmem.subset@meta.data$Patient=="PFCL1-LIM-313000"& Bmem.subset@meta.data$named.clusters=="Atypical"]
#PFCL1_UST_190762_AT <- Bmem.subset@meta.data$Full.Row.names[Bmem.subset@meta.data$Patient=="PFCL1-UST-190762"& Bmem.subset@meta.data$named.clusters=="Atypical"]
PFCL1_LIM_137402_AT <- Bmem.subset@meta.data$Full.Row.names[Bmem.subset@meta.data$Patient=="PFCL1-LIM-137402"& Bmem.subset@meta.data$named.clusters=="Atypical"]
#PFCL1_LIM_674950_AT <- Bmem.subset@meta.data$Full.Row.names[Bmem.subset@meta.data$Patient=="PFCL1-LIM-674950"& Bmem.subset@meta.data$named.clusters=="Atypical"]
PFCL1_179308_1988_AT <- Bmem.subset@meta.data$Full.Row.names[Bmem.subset@meta.data$Patient=="PFCL1-179308-1988"& Bmem.subset@meta.data$named.clusters=="Atypical"]
#PFCL1_196878_1691_AT <- Bmem.subset@meta.data$Full.Row.names[Bmem.subset@meta.data$Patient=="PFCL1-196878-1691"& Bmem.subset@meta.data$named.clusters=="Atypical"]
PFCL1_199474_1994_AT <- Bmem.subset@meta.data$Full.Row.names[Bmem.subset@meta.data$Patient=="PFCL1-199474-1994"& Bmem.subset@meta.data$named.clusters=="Atypical"]
PFCL1_LIM828246_AT <- Bmem.subset@meta.data$Full.Row.names[Bmem.subset@meta.data$Patient=="PFCL1-LIM828246"& Bmem.subset@meta.data$named.clusters=="Atypical"]


# Calculating mean Gene set enrichment values for each of these 24 groups by first subsetting the expression set
# See all the commented lines - this is done because for pseudobulking, we only include patients with n>20 cells in each cell subset.
PFCL1_LIM_313000_RM_matrix <- my_eset[,PFCL1_LIM_313000_RM]
#PFCL1_UST_190762_RM_matrix <- my_eset[,PFCL1_UST_190762_RM]
PFCL1_LIM_137402_RM_matrix <- my_eset[,PFCL1_LIM_137402_RM]
#PFCL1_LIM_674950_RM_matrix <- my_eset[,PFCL1_LIM_674950_RM]
PFCL1_179308_1988_RM_matrix <- my_eset[,PFCL1_179308_1988_RM]
#PFCL1_196878_1691_RM_matrix <- my_eset[,PFCL1_196878_1691_RM]
PFCL1_199474_1994_RM_matrix <- my_eset[,PFCL1_199474_1994_RM]
PFCL1_LIM828246_RM_matrix <- my_eset[,PFCL1_LIM828246_RM]

PFCL1_LIM_313000_ACT_matrix <- my_eset[,PFCL1_LIM_313000_ACT]
#PFCL1_UST_190762_ACT_matrix <- my_eset[,PFCL1_UST_190762_ACT]
PFCL1_LIM_137402_ACT_matrix <- my_eset[,PFCL1_LIM_137402_ACT]
#PFCL1_LIM_674950_ACT_matrix <- my_eset[,PFCL1_LIM_674950_ACT]
PFCL1_179308_1988_ACT_matrix <- my_eset[,PFCL1_179308_1988_ACT]
#PFCL1_196878_1691_ACT_matrix <- my_eset[,PFCL1_196878_1691_ACT]
PFCL1_199474_1994_ACT_matrix <- my_eset[,PFCL1_199474_1994_ACT]
PFCL1_LIM828246_ACT_matrix <- my_eset[,PFCL1_LIM828246_ACT]

PFCL1_LIM_313000_AT_matrix <- my_eset[,PFCL1_LIM_313000_AT]
#PFCL1_UST_190762_AT_matrix <- my_eset[,PFCL1_UST_190762_AT]
PFCL1_LIM_137402_AT_matrix <- my_eset[,PFCL1_LIM_137402_AT]
#PFCL1_LIM_674950_AT_matrix <- my_eset[,PFCL1_LIM_674950_AT]
PFCL1_179308_1988_AT_matrix <- my_eset[,PFCL1_179308_1988_AT]
#PFCL1_196878_1691_AT_matrix <- my_eset[,PFCL1_196878_1691_AT]
PFCL1_199474_1994_AT_matrix <- my_eset[,PFCL1_199474_1994_AT]
PFCL1_LIM828246_AT_matrix <- my_eset[,PFCL1_LIM828246_AT]

# Calculating gene set enrichment means
PFCL1_LIM_313000_RM_matrix_means <- rowMeans(exprs(PFCL1_LIM_313000_RM_matrix))
#PFCL1_UST_190762_RM_matrix_means <- rowMeans(exprs(PFCL1_UST_190762_RM_matrix))
PFCL1_LIM_137402_RM_matrix_means <- rowMeans(exprs(PFCL1_LIM_137402_RM_matrix))
#PFCL1_LIM_674950_RM_matrix_means <- rowMeans(exprs(PFCL1_LIM_674950_RM_matrix))
PFCL1_179308_1988_RM_matrix_means <- rowMeans(exprs(PFCL1_179308_1988_RM_matrix))
#PFCL1_196878_1691_RM_matrix_means <- rowMeans(exprs(PFCL1_196878_1691_RM_matrix))
PFCL1_199474_1994_RM_matrix_means <- rowMeans(exprs(PFCL1_199474_1994_RM_matrix))
PFCL1_LIM828246_RM_matrix_means <- rowMeans(exprs(PFCL1_LIM828246_RM_matrix))

PFCL1_LIM_313000_ACT_matrix_means <- rowMeans(exprs(PFCL1_LIM_313000_ACT_matrix))
#PFCL1_UST_190762_ACT_matrix_means <- rowMeans(exprs(PFCL1_UST_190762_ACT_matrix))
PFCL1_LIM_137402_ACT_matrix_means <- rowMeans(exprs(PFCL1_LIM_137402_ACT_matrix))
#PFCL1_LIM_674950_ACT_matrix_means <- rowMeans(exprs(PFCL1_LIM_674950_ACT_matrix))
PFCL1_179308_1988_ACT_matrix_means <-rowMeans(exprs(PFCL1_179308_1988_ACT_matrix))
#PFCL1_196878_1691_ACT_matrix_means <-rowMeans(exprs(PFCL1_196878_1691_ACT_matrix))
PFCL1_199474_1994_ACT_matrix_means <-rowMeans(exprs(PFCL1_199474_1994_ACT_matrix))
PFCL1_LIM828246_ACT_matrix_means <- rowMeans(exprs(PFCL1_LIM828246_ACT_matrix))

PFCL1_LIM_313000_AT_matrix_means <- rowMeans(exprs(PFCL1_LIM_313000_AT_matrix))
#PFCL1_UST_190762_AT_matrix_means <-rowMeans(exprs(PFCL1_UST_190762_AT_matrix))
PFCL1_LIM_137402_AT_matrix_means <-rowMeans(exprs(PFCL1_LIM_137402_AT_matrix))
#PFCL1_LIM_674950_AT_matrix_means <-rowMeans(exprs(PFCL1_LIM_674950_AT_matrix))
PFCL1_179308_1988_AT_matrix_means <- rowMeans(exprs(PFCL1_179308_1988_AT_matrix))
#PFCL1_196878_1691_AT_matrix_means <- rowMeans(exprs(PFCL1_196878_1691_AT_matrix))
PFCL1_199474_1994_AT_matrix_means <- rowMeans(exprs(PFCL1_199474_1994_AT_matrix))
PFCL1_LIM828246_AT_matrix_means <- rowMeans(exprs(PFCL1_LIM828246_AT_matrix))


# Creating the object
input_matrix <- cbind(PFCL1_LIM_313000_RM_matrix_means,PFCL1_LIM_137402_RM_matrix_means,
                      PFCL1_179308_1988_RM_matrix_means,PFCL1_199474_1994_RM_matrix_means,PFCL1_LIM828246_RM_matrix_means,
                      PFCL1_LIM_313000_ACT_matrix_means,PFCL1_LIM_137402_ACT_matrix_means,
                      PFCL1_179308_1988_ACT_matrix_means,PFCL1_199474_1994_ACT_matrix_means,PFCL1_LIM828246_ACT_matrix_means,
                      PFCL1_LIM_313000_AT_matrix_means,PFCL1_LIM_137402_AT_matrix_means,
                      PFCL1_179308_1988_AT_matrix_means,PFCL1_199474_1994_AT_matrix_means,PFCL1_LIM828246_AT_matrix_means)

my_exprs <- matrix(input_matrix,ncol(input_matrix),nrow=nrow(input_matrix))

my_phenoData <- data.frame(sample=colnames(input_matrix),
                           condition="RM",treated="no")

my_phenoData$condition[6:10] <- "Activated"
my_phenoData$condition[11:15] <- "Atypical"
my_feature_Data <- data.frame(geneID=rownames(input_matrix),geneSymbol=rownames(input_matrix))

my_eset <- ExpressionSet(my_exprs,
                         AnnotatedDataFrame(my_phenoData),
                         AnnotatedDataFrame(my_feature_Data))
rownames(my_eset) <- rownames(input_matrix)
colnames(my_eset) <- colnames(input_matrix)

my_eset$condition <- factor(my_eset$condition,levels = c("RM","Activated","Atypical"))

# First for fitting limma (linear models and differential expression for microarray data) model: Creating a design matrix.
# Number of rows is eaual to number of samples. Columns describe parameters in the model.
# I am using the +0 to have only expression levels of parameters in the columns of the model. NOT the expression levels of the second and the third parameter IN RESPECT to the 
# expression level of the parameter in the first column.
mod2 <- model.matrix(~ 0 + my_eset$condition,data = pData(my_eset))
colnames(mod2) <- c("RM", "Activated","Atypical")

# Fitting a linear model to the data.
fit2 <- lmFit(my_eset, mod2)

# Now we have to make a contrast. A contrast is a hypothesis. Basically we define the comparisons that we want to do.
# The second part in all the X-Y expressions is the reference level for the comparisons.
contrast.matrix <- makeContrasts("Activated-RM", "Atypical-RM", "Atypical-Activated", levels=mod2)
fit2c <- contrasts.fit(fit2, contrast.matrix)

# empirical Bayes moderation (eBayes) for calculation of moderated t-statistics, moderated F-statistic, and log-odds of differential expression.
fit2c <- eBayes(fit2c)

# Looking at the results for the indicidual comparisons defined before. 
topTable(fit2c,coef = "Activated-RM")[,7:8]
topTable(fit2c,coef = "Atypical-RM")
topTable(fit2c,coef = "Atypical-Activated")

# Creating a heatmap showing the results (In the example here, we are making a pre-selection to plot only the genesets which are differentially expressed between specific categories):
tt2 <- topTable(fit2c, coef="Atypical-RM", n=Inf)
DEpwys2 <- rownames(tt2)[tt2$adj.P.Val <= 0.05]
DEpwys_es2 <- exprs(my_eset[DEpwys2, ])
colorLegend2 <- c("darkgreen", "darkblue","darkred")
names(colorLegend2) <- c("RM", "Activated","Atypical")
sample.color.map2 <- colorLegend2[pData(my_eset)[, "condition"]]
names(sample.color.map2) <- colnames(DEpwys_es2)

scaled_mat = t(scale(t(DEpwys_es2)))

column_split = rep("RM", ncol(DEpwys_es2))
column_split[6:10] = "Activated"
column_split[11:15] = "Atypical"
column_split <- factor(column_split, levels = c("RM","Activated","Atypical"))

Seurat::PurpleAndYellow()
col_fun = circlize::colorRamp2(c(-2, 0, 2), c("#FF00FF", "#000000", "#FFFF00"))

heatmap <- Heatmap(scaled_mat, cluster_rows = FALSE, cluster_columns = FALSE,show_column_names = FALSE, 
        column_split = column_split,border = TRUE,column_gap = unit(1.5, "mm"),
        row_labels = gsub("_", " ", gsub("^GOBP_|^GOMF_", "",rownames(scaled_mat))),
        name = "Enrichment Score",row_names_side ="left",use_raster=F,col = col_fun)

draw(heatmap, padding = unit(c(2, 79, 2, 2), "mm")) # Export 4x13 or 3x11

