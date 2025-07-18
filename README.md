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

The file run_kinetic_module_workflow.m will execute the functions given under Code. Hereby run_preprocessing.m is the first and will execute itself functions provided in the folder Code/preprocessing/. The folder Code/balanced_concordant_complexes/ includes fucntions related to the identifcation of balanced and concordant complexes. The folder Code/kinetic/modules/ includes the functions used to identify and summarize the kinetic modules in metabolic networks. 

Functions inside the subfoldes of Code/ are not dependent on each other and can be used seperately with one exception which is Code/preprocessing/get_AY_matrix.r that depends on Code/preprocessing/deficiency_sparse. 

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
  (obtained from running script Results/write_M_MetSingle_M_MetDouble.m) \
  If the generated .csv file in folder MetSingle or MetDouble is empty no metabolite with absolute concentration robustness or metabolite pair with absolute concentration ratio robustness could be found for the respective network.

- Overview/
  The files combined resulted in Supplementary Table 1 of the related manuscipt (Table_ED-1.xlsx) 

  Table_ED-1.xlsx includes the following information
  - general information (column A-H, no color)
  - overview of results related to kinetic modules and their size (column I-P, green color)
  - overview of results related to absolute concentration (ratio) robustness (column Q-X, orange color)
  - overview of results from stoichiometric coupling (column Y-AD, blue color)

  column  | information \
  ---------------------------------- \
  A       | Model number for counting \
  B       | Binding type considered (ordered or random) \
  C       | Model name  \
  D       | Organism\
  E       | number of complexes in the model\
  F       | number of reactions in the model\
  G       | number of free metabolites in the model\
  H       | number of free metabolites, enzymes and enzyme-metabolite   complexes in the model \
  I       | number of kinetic modules \
  J       | number of kinetic modules with one complex \
  K       | maximum size (number of complexes) of kinetic module \
  L       | average size of kinetic modules \
  M       | standard deviation of size of kinetic modules \
  N       | % of complexes from all model complexes found in largest kinetic module \
  O       | number of reactions whose substrate is in the largest kinetic module \
  P       | % of reactions with substrate in largest kinetic module from all model reactions \
  Q       | number of components with absolute concentration robustness \
  R       | % of components with absolute concentration robustness from all components (see column H) \
  S       | number of component pairs with absolute concentration ratio robustness \
  T       | % of component pairs with absolute concentration ratio robustness from all unique combinations of two components \
  U       | number of free metabolites with absolute concentration robustness \
  V       | % of free metabolites with absolute concentration robustness from all free metabolites (see column G) \
  W       | number of free metabolite pairs with absolute concentration ratio robustness \
  X       | % of free metabolite pairs with absolute concentration ratio robustness \
  Y       | number of modules identified from stoichiometric coupling \
  Z       | number of modules with one complex identified from stoichiometric coupling \
  AA      | maximum size of module identified from stoichiometric coupling \
  AB      | average size of modules identified from stoichiometric coupling \
  AC      | standard deviation size of modules identified from stoichiometric coupling \
  AD      | fraction size of largest kinetic module to size of largest module identified from stoichiometric coupling

  The script CreateResultTable.m combinds the Overview_*.csv results for the individual networks to one joint table.

- Reactions_Giant/
  List of reaction indices that belong to the giant kinetic module

- stoichiometric coupling
  matrices obtained form F2C2

- files *.RData in Results/
  results of running the Upstream algorithm in Code/kinetic_modules/code_kineticModule_analysis.R

- the scipt write_M_MetSingle_M_MetDouble.m in the Results/ folder filters the free metabolites from the set of all identified components (enzymes, enzyme-substrate complexes, free metabolites) with absolute concentration robustness or pairs with absolute concentration ratio robustness

_Note: intermediate results of few models are missing as they were to large to be uploaded_


  

