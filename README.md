# Kinetic modules in biochemical networks

The repository contains code, data, (intermediate) results that allow identification of 

  - balanced complexes
  - concordant complexes
  - kinetic modules

The identification of kinetic modules allows to find metabolites of absolute concentration robustness and pairs with absolute concentration ratio robustness in large-scale metabolic networks.

## Dependencies

- Matlab (tested with 2023b / 2024a)
- R (tested with R-4.3.0), packages igraph and R.matlab
- COBRA toolbox (https://opencobra.github.io/cobratoolbox/stable/index.html)

  to compare result with full coupling based on stoichiometry
  - F2C2 tool (https://pubmed.ncbi.nlm.nih.gov/22524245/) (used with glpk)
 
## Main functions and scripts

### To run the whole pipeline

To run the whole pipeline execute run_kinetic_module_workflow.m \
The code executes the workflow for the Arabidopsis Core Model located in Models/original \
To run the code for an other model modify variables files and f as described in the script

### Models
The original model files downloaded from BiGG database or original publication can be found under Models/original/ \
An over view of the original publications is provided in Overview_model_original_publications.xlsx

Model reactions were split assuming a fixed order of substrate binding or random binding \
The resulting models after preprocessing and splitting can be found in Models/models_with_elementary_steps/

### Code

Functions to run individual steps e.g. identification of balanced and concordant complexes can be found in folder Code/

### Results

The folder  Results/ includes the following subfolder:

- concordant/
  the files include the model structure after splitting, the set of balanced complexes (B), concordant complexes (CC)
  for details on model variables check comments in the example workflow or individual functions

- Figures/
  code to generate result figures

- MetDouble/ and MetSingle
  metabolite pairs with absolute concentration ratio robustness and metabolites with absolute concentration robustness, respectively \
  MetDouble_* and MetSingle_* also considers enzymes and enzyme-substrate complexes \
  M_MetDouble_* and M_MetSingle_* includes free metabolites only \
  (obtained from  running script Results/write_M_MetSingle_M_MetDouble.m)

- Overview/
  The files combined resulted in Supplementary Table 1 of the related manuscipt

- Reactions_Giant/
  List of reaction indices that belong to the giant kinetic module

- stoichiometric coupling
  matrices obtained form F2C2

- files *.RData in Results/
  results of running the Upstream algorithm in Code/kinetic_modules/code_kineticModule_analysis.R

_Note: intermediate results of few models are missing as they were to large to be uploaded_


  

