//
// Created by george on 2019-11-13.
//

#include <def.h>
#include <test_utils.h>
#include <sparse_utilities.h>
#include <sparse_blas_lib.h>
#include <sparse_fusion.h>
#include "catch.hpp"

using namespace sym_lib;
using namespace std;

TEST_CASE("Check_Jacobi", "[JacobiChecks]") {
 SECTION("unfused") {

 }

 SECTION("fused_spmv") {

 }
}


TEST_CASE("Check GS", "[GaussSeidelChecks]") {
 SECTION("unfused") {

 }

 SECTION("spmv_sptrsv_fuse") {

 }

 SECTION("sptrsv_spmv_fused") {

 }
}
