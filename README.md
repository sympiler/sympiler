# sympiler
Sympiler: Transforming Sparse Matrix Codes by Decoupling Symbolic Analysis

This repository is only created for Artifact Evaluation of SC17 conference 
and it is not to share with other people. The repository will be shared for 
public after its paperwork finished.

## Installation
### requirements
CHOLMOD, Eigen and OpenBLAS libraries.

### Instruction


## Source Tree Description
### libTest
This folder is inteded to evaluate the two library competitors i.e., CHOLMOD 
and Eigen. 

### matrixDB
The scripts to download the matrices used in the paper.

### symGen
The Sympiler-generated code is stored here by default and is used for testing 
Sympiler in symTest folder. 

### sympiler
This folder contains the source of Sympiler. Running this code generates code
for a specific sparsity.

### symTest
This folder tests the Sympiler-generated code for a given matrix.

