# Sympiler
Kazem Cheshmi, Shoaib Kamil, Michelle Mills Strout, and Maryam Mehri Dehnavi. 2017. Sympiler: Transforming Sparse Matrix Codes by Decoupling
Symbolic Analysis. In Proceedings of SC17, Denver, CO, USA, November 12–17, 2017,  13 pages.  DOI: 10.1145/3126908.3126936


## Installation
### requirements
CHOLMOD, Eigen and OpenBLAS libraries.

### Building the project
You should use cmake to build the LLVM:
```bash
%cd <where you cloned Sympiler>
%mkdir build
%cd build
% cmake ..
% make -j8
```
This will build all parts of the project including Sympiler, Sympiler tests, 
library tests, and required utility functions. 

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

