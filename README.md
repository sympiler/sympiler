![APM](https://badgen.net/github/license/micromatch/micromatch)
![example workflow](https://github.com/sympiler/sympiler/actions/workflows/cmakeUbuntu.yml/badge.svg)
![example workflow](https://github.com/sympiler/sympiler/actions/workflows/cmakeMac.yml/badge.svg)

# Sympiler

Sympiler is a code generator for transforming sparse matrix methods.
To access the list of publication and resources please visit: http://www.sympiler.com/


## Quick Build Guide for Impatient Users

If you have CMake 3.16 or higher and a C++11 compiler, then:

```bash
git  clone --recursive https://github.com/sympiler/sympiler.git
cd sympiler
cmake -DCMAKE_BUILD_TYPE=Release  -S . -B build
cmake --build build --config Release -j 6 
```

For details, please see the table below.

# Table of Contents:

* [Building Sympiler](https://www.sympiler.com/docs/getting-started-sympiler/#building-sympiler)
* [Sympiler Overview](https://www.sympiler.com/docs/)
    * [Iteration-space Pruning](https://www.sympiler.com/docs/prune/)
    * [Loop Tiling](https://www.sympiler.com/docs/tiling/)
    * [Loop Fusion](https://www.sympiler.com/docs/fusion/)
    * [Vectorization](https://www.sympiler.com/docs/vect/)
* [Using Sympiler/Aggregation in C++](https://www.sympiler.com/docs/sympiler-lib/)
* [Benchmarks](https://www.sympiler.com/docs/benchmark/)
* [Publications](https://www.sympiler.com/#publications)
* [Citation](https://www.sympiler.com/docs/citation/)
* [Sympiler Homepage](https://nasoq.github.io/)
* [Sympiler Documentation](https://nasoq.github.io/docs/)
* [GitHub](https://github.com/sympiler/sympiler)
* [Twitter](https://twitter.com/sympiler)
