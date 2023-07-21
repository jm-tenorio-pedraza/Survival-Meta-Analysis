# ICB-Meta-Analysis
 Supplementary files for manuscript

This is the repository for the publication "Meta-analysis of preclinical measures of efficacy in immune checkpoint blockade therapies and comparison to clinical efficacy estimates" to be published in Translational Medicine Communications.

This repository is composed of the following folders:
1. Data
2. Functions
3. Scripts
4. Ouput

The data folder contains the .csv files related to the preclinical studies and their respective meta-data (i.e. Cell lines and bias preventing measures), as well as the clinical studies in immuno-oncology.

The functions folder contains an .R script where all the data-processing functions are stored and are required for the main analysis workflows found in the scripts folder.

The scripts folder contains four workflows related to the different analysis presented in the manuscript. The Preclinical_Workflow.R presents the analyses performed using the preclinical data. The Clinical_Workflow.R contains the analyses pertaining to the clinical data. The Translational_Workflow.R details the analyses performed to relate the preclinical estimates to the clinical estimates. The 8_ExpDesign.R contains the simulation setup to optimize the design of preclinical experiments.

The output folder contains all the results from the scripts folder.


