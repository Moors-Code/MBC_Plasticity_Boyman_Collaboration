# MBC_Plasticity_Moor_Boyman_Collaboration

## Citation
This code was used to analyze the single cell sequencing data for the project 

"Fate and plasticity of B cell memory and recall response against SARS-CoV-2 in humans" 

by Y. Zurbuchen and J. Michler et al. 

## Data availability

The data is deposited at zenodo.org under access number: DOI: 10.5281/zenodo.7463895

## Code structure

The code is separated into eight scripts:

**1) MBC_Preprocessing.R:** This part is used for preprocessing of the raw count matrices, dataset integration and dimensional reduction. In the last section of this code, the V(D)J segment usage between Spike RBD Binders and Non-Binders is compared.

**2) MBC_Preprocessing_2.R** This part is very similar to part 1 (MBC_Preprocessing.R). It differs only in the last section, where the code is used to generate V(D)J segment usage comparisons between MBC subsets.

**3) MBC_Analysis.R:** This part contains a big part of the downstream analysis including DGEA, GSEA and Monocle 3 Analysis.

**4) GSVA.Server.R:** This part contains the Gene Set Variation Analysis.

**5) MBC_Immcantation_Input.R:** This part is used to create input files for the Immcantation changeo-10x pipeline.

**6) MBC_Immcantation_Analysis.R:** This part uses Immcantation changeo-10x pipeline outputs and is used for downstream clonal analysis.

**7) Tonsil_Cohort.R**: This code was used to analyze the Tonsil-PBMC paired analysis samples. The code contains preprocessing, Immcantation Input and Output analysis and downstream processing of the dataset.

**8) Vaccination_Cohort.R**: This code was used to analyze the sample of the Vaccination Cohort. The code contains preprocessing, Immcantation Input and Output analysis, including the comparison of SHMs between timepoints in this cohort.

