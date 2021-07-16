//
// Created by george on 2020-02-15.
//

#ifndef FUSION_BCSCMATRIX_H
#define FUSION_BCSCMATRIX_H

#include "def.h"
namespace sym_lib {
 class BCSCMatrix {
  BCSC *M;
 public:
  /// finds all the different supernodes / blocks in the matrix
  int supernodes(CSC *A, int limit, bool isLower);

  /// calculate the total number of nonzero in the matrix including zero padding
  int calcSize(CSC *A);

  /// Fill in the row and numerical values for blocked matrix
  void createFormat(CSC *A);

  // Extract a supernodal CSC from a BCSC
  CSC *compressed_BCSC_to_CSC();

  /// Generates the BCSC arrays from A
  /// \param nA
  /// \param A
  void generateBCSC(CSC *A) {
   M->nnz = calcSize(A);
   M->i = new int[M->nnz]();
   M->x = new double[M->nnz]();
   M->nrows = new int[M->nodes]();
   createFormat(A);
  }

 public:
  BCSCMatrix(CSC *A) {
   M = new BCSC(A);
   M->nodes = supernodes(A, A->n, true);
   generateBCSC(A);
  }

  BCSCMatrix(CSC *A, int nodes, int *supernodes) {
   M = new BCSC(A);
   M->nodes = nodes;
   M->supernodes = new int[nodes + 1]();
   std::memcpy(M->supernodes, supernodes, sizeof(int) * (nodes + 1));
   generateBCSC(A);
  }

  BCSC *getBCSC() { return M; }

  ~BCSCMatrix() {
   delete (M);
  }
 };
}

#endif //FUSION_BCSCMATRIX_H
