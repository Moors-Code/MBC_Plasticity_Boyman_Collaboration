# MBC_Plasticity_Moor_Boyman_Collaboration

## Citation
This code was used to analyze the single cell sequencing data for the project 

"Fate and plasticity of SARS-CoV-2-specific B cells during memory and recall in humans" 

by Zurbuchen and Michler et al. 

## Data availability

The data is deposited at zenodo.org under access number:

## Code structure

The code is separated into five scripts:

**1) Memory_B_cell_COVID_Preprocessing.R:** This part is used for preprocessing of the raw count matrices, dataset integration and dimensional reduction.

**2) Memory_B_cell_COVID_Analysis.R:** This part containsa big part of the downstream analysis including DGEA, GSEA and Monocle 3 Analysis.

**3) GSVA.Server.R:** This part contains the Gene Set Variation Analysis.

**4) Immcantation_Input.R:** This part is used to create input files for the Immcantation changeo-10x pipeline.

**5) Outs.analysis_Persistent_singlet_update.R:** This part uses Immcantation changeo-10x pipeline outputs and is used for downstream clonal analysis.

