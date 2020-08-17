# Sympiler
Sympiler is a code generator for transforming sparse matrix methods.
To access the list of publication and resources please visit: http://www.sympiler.com/

**We are actively working on Sympiler code and a major release will come soon.**

**ParSy** is parallel version of Sympiler. The evaluation benchmark for ParSy is
available from ParSy_bench repository: https://github.com/cheshmi/parsy_bench

**NASOQ and LBL** solvers use Sympiler and ParSy code internally. For more information visit [NASOQ Webpage](https://nasoq.github.io/).

## Installation

### Library requirements
Sympiler does not need any external library however for testing the
Sympiler-generated code, Suitesparse and Intel MKL libraries are required.

### Building Sympiler
```bash
cd where/you/cloned/Sympiler
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make sympiler_proj
```


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

