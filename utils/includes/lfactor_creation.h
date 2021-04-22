//
// Created by kazem on 6/11/20.
//

#ifndef FUSION_LFACTOR_CREATION_H
#define FUSION_LFACTOR_CREATION_H

#include "def.h"

namespace sym_lib{
 struct LFactorSymbolicInfo{
  int *col_counts, *row_counts;
  int *perm, *iperm;
  int *post, *ipost;
  int *parent;
  double *accessed_nnz_per_col;
  CSC *L;
  size_t flops;
  LFactorSymbolicInfo(int n){
   col_counts = new int[n]();
   row_counts = new int[n]();
   perm  = new int[n]();
   iperm = new int[n]();
   post = new int[n]();
   ipost = new int[n]();
   parent = new int[n]();
   accessed_nnz_per_col = new double[n]();
   L=NULLPNTR; flops=0;
  }

  ~LFactorSymbolicInfo(){
   delete[] col_counts;
   delete[] row_counts;
   delete[] perm;
   delete[] iperm;
   delete[] post;
   delete[] ipost;
   delete[] parent;
   delete[] accessed_nnz_per_col;
   delete L;
  }

 };

 /// Compute the post order of the input tree
 /// \param Parent input size n tree
 /// \param n input size of tree
 /// \param Weight optional input size n specifies the weight of each tree node
 /// \param Post output size n
 /// \return n if it worked well
 int tree_post_order (int *Parent, size_t n, int *Weight, int *Post);


 /// Compute the row/col counts of the L-factor of matrix A
 /// \param A input matrix
 /// \param Parent input etree, Parent [i] = p if p is the parent of i
 /// \param Post input size A->m.  Post [k] = i if i is the kth node in
 //			      the postordered etree.
 /// \param RowCount output size A->m. RowCount [i] = # entries in the ith row of
 //			       L, including the diagonal.
 /// \param ColCount output size A->m. ColCount [i] = # entries in the ith
 //			      column of L, including the diagonal.
 /// \param First output size A->m.  First [i] = k is the least postordering
 //			       of any descendant of i.
 /// \param Level output  size A->m.  Level [i] is the length of the path from
 //			  i to the root, with Level [root] = 0.
 /// \param fl output number of FLOPs
 /// \return
 bool row_col_counts_symmetric(CSC *A, int *Parent, int *Post, int *RowCount,
   int *ColCount, int *First, int *Level, size_t &fl);

 /// Building the CSC matrix (L-factor) from column count
 /// \param n
 /// \param A in original matrix
 /// \param col_count in colcount of L
 /// \param parent in etree
 /// \param accessed_nnz out the number of accessed nonzero values in
 /// left-looking factorization per col
 /// \return the pattern of L
 CSC *build_L_pattern_from_col_counts(int n, CSC *A,
                                      int *col_count, int *parent,
                                      double *accessed_nnz);

 LFactorSymbolicInfo *build_symbolic_simplicial_lfactor(CSC *A);



}
#endif //FUSION_LFACTOR_CREATION_H
