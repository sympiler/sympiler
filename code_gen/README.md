
## Installation

The code generator is under major changes. We are working on it. 


### Building the Sympiler tests
The first step is to set the environmental variables corresponding
to each library. The following shows how the variables are set in bash.
```bash
export MKLROOT <path to MKL>
export SUITEROOT <path to Suitesparse>
```
After setting the library paths:

```bash
cd where/you/cloned/Sympiler/build/
make
```
This will build the two remaining parts of the project including
Sympiler tests for both Cholesky and Triangular solve.

### Matrix dataset
In order to download and generate the matrices, the following commands need to be ran.
```bash
cd ..
mkdir matrixDB
./scripts/dlMat.sh matrixDB/ build/libTest/cholesky/TriangGen
```
After downloading and generating the required matrices you can evaluate Sympiler.

### Evaluating Sympiler
After build is done successfully, the following commands can be used 
to evaluate Sympiler generated code:

#### Cholesky
```bash
./build/symTest/cholesky/symChol matrixDB/<MATRIX NAME>.mtx matrixDB/ccache/rajat21.mtx

```

#### Triangular solve
```bash
./build/symTest/triangular/SymTriang matrixDB/triangular/<MATRIX NAME>_trns.mtx
```

## Source Tree Description

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

