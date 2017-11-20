# Sympiler
Sympiler is a code generator for transforming sparse matrix methods. 
The current code supports the features explained in the following paper:

Kazem Cheshmi, Shoaib Kamil, Michelle Mills Strout, and Maryam Mehri Dehnavi. 2017. [Sympiler: Transforming Sparse Matrix Codes by Decoupling
Symbolic Analysis](https://dl.acm.org/citation.cfm?id=3126936&CFID=825768759&CFTOKEN=28284703). In Proceedings of SC17, Denver, CO, USA, November 12–17, 2017,  13 pages.  DOI: 10.1145/3126908.3126936


## Installation
### Library requirements
CHOLMOD, Eigen and OpenBLAS libraries need to be installed and their corresponding variables 
need to be set in the CMakeLists.txt file in the root.


### Building the project
You should use cmake to build the LLVM:
```bash
cd where/you/cloned/Sympiler
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make 
```
This will build all parts of the project including Sympiler, Sympiler tests, and 
library tests. 

### Matrix dataset
In order to download and generate the matrices, the following commands need to be ran.
```bash
cd ..
mkdir matrixDB
./scripts/dlMat.sh matrixDB/ build/libTest/cholesky/TriangGen
```
After downloading and generating the required matrices you can evaluate Sympiler and 
other libraries.

###Evaluating Sympiler
After build is done successfully, the following commands can be used 
to evaluate Sympiler generated code:

####Cholesky
```bash
./build/symTest/cholesky/symChol matrixDB/<MATRIX NAME>.mtx matrixDB/ccache/rajat21.mtx

```

####Triangular solve
```bash
./build/symTest/triangular/SymTriang matrixDB/triangular/<MATRIX NAME>_trns.mtx
```

###Evaluating the libraries
The two mentioned libraries i.e., Eigen and CHOLMOD can be evaluated as following:
for Cholesky in CHOLMOD and Eigen:
```bash
./build/libTest/cholesky/CholeskyTestCHOLMOD matrixDB/<MATRIX NAME>.mtx matrixDB/ccache/rajat21.mtx
./build/libTest/cholesky/cholEigen matrixDB/<MATRIX NAME>.mtx matrixDB/ccache/rajat21.mtx
```
For Triangular Solve in Eigen:
```bash
./build/libTest/triangular/trnsEigen matrixDB/triangular/<MATRIX NAME>_trns.mtx
```
## Source Tree Description
### libTest
This folder is inteded to evaluate the two library competitors i.e., CHOLMOD 
and Eigen. 

### scripts
The scripts to download the matrices used in the paper.

### symGen
The Sympiler-generated code is stored here by default and is used for testing 
Sympiler in symTest folder. 

### sympiler
This folder contains the source of Sympiler. Running this code generates code
for a specific sparsity.

### symTest
This folder tests the Sympiler-generated code for a given matrix.

