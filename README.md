# Computer codes for Francois and Jay Nat. Commun. (2020)

This repository contains R codes for simulation and data analyses from (Francois and Jay, Nat. Commun. 2020). Figures refer to the main text of the manuscript.

#### Codes for Figure 1
These codes are in the folders **Simulation_divergence_models** and
**Simulation_single_population_models**

#### Codes for Figure 1 and for Figure 2
These codes are in the folder **Simulation_admixture_models**. These analyses require having ADMIXTOOLS installed on the computer system (Linux is preferable for ADMIXTOOLS analyses). We provided a notebook that explains how any analysis in the study could be reproduced with the packages ``tfa` (this paper), `LEA` (bioconductor), and `admixr` (CRAN).

#### Codes for Figure 3 and related supplementary figures
These codes are in the folder **Data_analysis_704_samples**, reproducing FA and PC projections for all ancient samples in the study. The data must be loaded before analysis. Open the R script and copy the url for downloading them from an external repository (Github does not accept large data sets).

#### Codes for Figure 4 and related supplementary figures. Figures S10 - S12 
These codes are in the folder **Two_way_admixture_analyses**. Two-way admixture analyses of European samples. These analyses require to have ADMIXTOOLS installed, and to place some data files in the "snp" folder (see README file).

#### Codes for Figure 5 and 6 and related supplementary figures
These codes are in the folder **Example_European_samples** and reproduce the analyses in Figure 5 and Figure 6


**Important Note**: Some scripts have been improved by the latest release of the R package 'tfa', and their use may be obsolete. Prefer using the latest version of the package.


