//
// Created by kazem on 4/23/17.
//

// Always use assert, even if llvm-config defines NDEBUG
#ifdef NDEBUG
#undef NDEBUG
#include <assert.h>
#define NDEBUG
#else
#include <assert.h>
#endif

#ifndef SYMPILER_PROJ_UTIL_H
#define SYMPILER_PROJ_UTIL_H

/** \file
 * Various utility functions used internally Halide. */

#include <cstdint>
#include <utility>
#include <vector>
#include <string>
#include <cstring>

#if defined(_WIN32) && defined(Halide_SHARED)
#ifdef Halide_EXPORTS
#define EXPORT __declspec(dllexport)
#else
#define EXPORT __declspec(dllimport)
#endif
#else
#define EXPORT
#endif

namespace Sympiler {
namespace Internal {

/** Generate a unique name starting with the given prefix. It's unique
* relative to all other strings returned by unique_name in this
* process.
*
* The single-character version always appends a numeric suffix to the
* character.
*
* The string version will either return the input as-is (with high
* probability on the first time it is called with that input), or
* replace any existing '$' characters with underscores, then add a
* '$' sign and a numeric suffix to it.
*
* Note that unique_name('f') therefore differs from
* unique_name("f"). The former returns something like f123, and the
* latter returns either f or f$123.
*/
// @{
std::string unique_name(char prefix);

std::string unique_name(const std::string &prefix);

void reset_counter(char prefix);

    //Matrix related operations
bool readCSCMatrixPattern(std::string fName, int &n, int& m, int &NNZ, int* &col, int* &row);

}
}
#endif //SYMPILER_PROJ_UTIL_H
