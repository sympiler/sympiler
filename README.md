# Sympiler
Sympiler is a code generator for transforming sparse matrix methods.
To access the list of publication and resources please visit: http://www.sympiler.com/


**ParSy** is parallel version of Sympiler. The evaluation benchmark for ParSy is
available from ParSy_bench repository: https://github.com/cheshmi/parsy_bench


**NASOQ and LBL** solvers use Sympiler and ParSy code internally. For more information visit [NASOQ Webpage](https://nasoq.github.io/).

## Install

### Dependencies

#### METIS
Set the following labels to enable using METIS:
```bash
export METISINC /path/to/metis.h
export METISLIB /path/to/libmetis.so/.a
```
If you don't set these variables, the build will continue without using METIS.

#### Intel MKL
Intel MKL should automatically set `MKLROOT` or `MKL_ROOT`, if not, please set:
```bash
export MKLROOT /path/to/mkl
```
If you don't set these variables, the build will continue without using Intel MKL for benchmark

#### PAPI
PAPI is an open-source profiler that can be obtained from https://icl.utk.edu/papi/
After installation, please set `PAPIROOT` to enable profiling 
```bash
export PAPIROOT /path/to/papi
```
PAPI works only on Linux and is optional for profiling only.

### Build

```bash
cd where/sympiler/is/cloned/
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make 
```


