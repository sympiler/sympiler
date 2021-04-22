//
// Created by kazem on 10/9/19.
//

#ifndef FUSION_SPARSE_UTILITIES_H
#define FUSION_SPARSE_UTILITIES_H


#include <string>
#include <iostream>
#include <iomanip>
#include <def.h>

namespace sym_lib {

 ///
 /// \param fname
 /// \param vec
 /// \param length
 void save_vector(char *fname, const double *vec, int length);

 ///
 /// \param nrow
 /// \param ncol
 /// \param Ap
 /// \param Ai
 /// \param Ax
 /// \param rowptr
 /// \param colind
 /// \param values
 /// \return
 int
 csc_to_csr(int nrow, int ncol, int *Ap, int *Ai, double *Ax, int *&rowptr,
            int *&colind, double *&values);

 ///
 /// \param nrow
 /// \param ncol
 /// \param Ai
 /// \param Ap
 /// \param Ax
 /// \param colptr
 /// \param rowind
 /// \param values
 /// \return
 int csr_to_csc(int nrow, int ncol, int *Ai, int *Ap, double *Ax, int *&colptr,
                int *&rowind, double *&values);


 ///
 /// \param A input matrix in CSC format
 /// \return converted to CSR format
 CSR* csc_to_csr(CSC* A);

 ///
 /// \param A input matrix in CSR format
 /// \return converted to CSC format
 CSC* csr_to_csc(CSR* A);

 ///
 /// Genreate a tridiagonal matrix based on the given diagonal values
 ///
 /// \param n rank of matrix to generate
 /// \param a0 above diagonal element
 /// \param a1 diagonal element
 /// \param a2 below diagonal element
 /// \return the tridiagonal matrix
 CSC* tridiag(int n, double a0, double a1, double a2);

/// Converts tree to CSC
/// \param n
/// \param tree
/// \return
 CSC* tree_to_csc(int n, int *tree);

 ///
 /// \param A
 /// \return
 CSC* transpose_general(CSC *A);

 ///
 /// \param An
 /// \param Anz
 /// \param Ap
 /// \param Ai
 /// \param Ax
 /// \param lower if true symmetric with only lower part otherwise the upper.
 /// \return
 CSC* make_half(size_t An, int *Ap, int *Ai, double *Ax, bool lower=true);

 /// Makes a symmetric matrix full
 /// \param A
 /// \return
 CSC* make_full(CSC *A);


 /// transpose matrix A after applying the given permutation, if give.
 /// \param A Input symmetric matrix, either lower or upper part should be stored
 /// \param Perm
 /// \return
 CSC *transpose_symmetric(CSC *A, int *Perm);

 /// Compute the norm of the difference of two dense vectors
 /// \param beg_idx
 /// \param end_idx
 /// \param x0
 /// \param x1
 /// \return
 double residual(int beg_idx, int end_idx, const double *x0, const double *x1);

 /// Computer the norm of a dense vector
 /// \param size
 /// \param vec
 /// \return
 double norm(int size, double *vec);

/// Computes the inverse of a given permutation
/// \param n
/// \param perm
/// \param pinv
/// \return
  int *compute_inv_perm(int n, int *perm, int *pinv);


  /// Makes a copy of matrix A and returns it.
  /// \param A
  /// \return
  CSC *copy_sparse(CSC *A);


  ///
  /// \param beg
  /// \param end
  /// \param src
  /// \param dst
  void copy_vector_dense(size_t beg, size_t end, const double *src, double *dst);

  /// Copy CSC to CSC
  /// \param src
  /// \param dst
 void copy_from_to(CSR *src, CSR *dst);
 void copy_from_to(CSC *src, CSC *dst);

/// Finds the number of empty columns.
/// \param A 
/// \return
  int number_empty_col(CSC *A);


 /// Mixed BFS and topological sort for BCSC
 /// \param n
 /// \param Lp
 /// \param Li_ptr
 /// \param Li
 /// \param sup2col
 /// \param col2sup
 /// \param inDegree
 /// \param visited
 /// \param node2partition
 /// \param levelPtr
 /// \param levelSet
 /// \param bfsLevel
 /// \param newLeveledParList
 /// \return
 int modified_BFS_BCSC(int n, size_t *Lp, size_t *Li_ptr, int* Li,
                      const int* sup2col, const int* col2sup, int *inDegree,
                      bool *visited, int *node2partition, int* &levelPtr,
                      size_t* levelSet, int bfsLevel,
                      std::vector<std::vector<int>> &newLeveledParList);


 /// Mixed BFS and topological sort for CSC format.
 /// \param n
 /// \param Lp
 /// \param Li
 /// \param inDegree
 /// \param visited
 /// \param node2partition
 /// \param levelPtr
 /// \param levelSet
 /// \param bfsLevel
 /// \param newLeveledParList
 /// \return
 int modified_BFS_CSC(int n, int *Lp, int *Li, int *inDegree, bool *visited,
                      int *node2partition, int *&levelPtr, int *levelSet,
                      int bfsLevel,
                      std::vector<std::vector<int>> &newLeveledParList);


 /// Stores the input tree in a list to facilitate accessing to children
 /// \param n
 /// \param eTree
 /// \param childPtr
 /// \param childNo
 /// \param nChild
 void populate_children(int n, const int *eTree, int *childPtr,
                        int *childNo, int *nChild);


 /// Computes the depth of node in tree
 /// \param node
 /// \param n
 /// \param tree
 /// \param weight
 /// \return
 int get_node_depth(int node, int n, const int *tree, int *weight = NULLPNTR);


 /// Returns the height of tree by up-traversing from all nodes.
 /// \param n
 /// \param tree
 /// \param weight
 /// \return
 int get_tree_height_bruteforce(int n, const int *tree, int *weight = NULLPNTR);

 ///  /// Returns the height of tree by up-traversing from all children nodes.
 /// \param n
 /// \param tree
 /// \param nChild1
 /// \param weight
 /// \return
 int get_tree_height(int n, const int *tree, int *nChild1,
                     int *weight = NULLPNTR);
 int get_tree_height_efficient(int n, const int *tree, const int *nChild1,
                               const int *weight=NULLPNTR);

 /// Computes the cost of given tree using the weight vector.
 /// \param n
 /// \param tree
 /// \param nChild
 /// \param weight
 /// \return
 double *compute_subtree_cost(int n, const int *tree, double *weight);


 /// converting sparse csc to dense matrix
 /// \param A input CSC matrix
 /// \param output D array of size m*n
 void sparse2dense(CSC *A, double *D);


 /// Generates a diagonal matrix with val on diagonal
 /// \param n
 /// \param val
 /// \return
 CSC *diagonal(int n, double val);


 /// Extracting diagonal locations of a compressed matrix
 /// \param n
 /// \param Ap
 /// \param Ai
 /// \return
 int *extract_diagonals(const int n, const int *Ap, const int *Ai);

/// Merge multiple sequential graphs into a single graph
/// \param ngraphs number of graphs to fuse
/// \param n dimension of each graph assuming 1-to-1 correspondence
/// \param Gps list of graph pointers
/// \param Gis list of graph indices
/// \param nGp the new graph pointer
/// \param nGi the new graph indices
 // TODO create function to accept a DAG of graph pointers for better fusion
 void merge_graph(int ngraphs, int n, int **Gps, int **Gis, int *&nGp, int *&nGi);


/// Merges multiple graph and returns a CSC matrix
/// \param ngraphs
/// \param n
/// \param Gps
/// \param Gis
/// \return
 CSC* merge_graph(int ngraphs, int n, int **Gps, int **Gis);

 /// Takes ngraphs and n-1 dependence graph and merge them
 /// G1 -> DG1 -> G2 -> DG2 -> G3
 /// \param ngraphs
 /// \param n
 /// \param Gps
 /// \param Gis
 /// \param nd
 /// \param DGps
 /// \param DGis
 /// \return
 CSC* merge_DAGs_with_partial_order(int ngraphs, int n, int **Gps,
   int **Gis, int nd, int **DGps, int **DGis);
 ///
 /// \param inSize
 /// \param inTree
 /// \param inCost
 /// \param inChildPtr
 /// \param inChildNo
 /// \param nChild
 /// \param n
 /// \param partitionNum
 /// \param outSize
 /// \param outCost
 /// \param outNode2Par
 /// \param parList
 /// \return
 int post_order_spliting(int inSize, int *inTree, double *inCost,
                         int *inChildPtr, int *inChildNo,//Children list
                         int *nChild,int n, int partitionNum,
   /*Outputs*/
                         int &outSize, double* outCost,
                         int* outNode2Par,
                         std::vector<std::vector<int>> &parList);


 /// x = x + beta * A(:,j), where x is a dense vector and A(:,j) is sparse
 /// \param A
 /// \param j
 /// \param beta
 /// \param w
 /// \param x
 /// \param mark
 /// \param C
 /// \param nz
 /// \return
 int scatter (const CSC *A, int j, double beta, int *w,
              double *x, int mark, CSC *C, int nz);


 /// Make A symmetric by (A+At)/2
 /// \param A
 /// \param lower if true returns lower part of the symmetrized matrix
 /// \return 
 CSC * make_symmetric(CSC *A, bool lower=true);


 /// Computes nnz per column assuming matrix is general for now, FIXME
 /// \param A
 /// \param nnz_cnt, output
 void compute_nnz_per_col(CSC *A, double *nnz_cnt);


 /// Reorder array arr with a the permutation perm
 /// \param n size of array
 /// \param arr in/out
 /// \param perm in
 /// \param ws workspace temporary array of size n
 void reorder_array(int n, int *arr, int *perm, int *ws);


 /// computing the inverse of permutation
 /// \param n
 /// \param perm in
 /// \param iperm out
 void inv_perm(int n, int *perm, int *iperm);

/// Computes A(p,q)
/// \param A  is CSC
/// \param pinv inverse of p
/// \param q row perm
/// \return
 CSC *permute_general (const CSC *A, const int *pinv,
                       const int *q);
/// Coarsening every k rows/cols in a CSR/CSC matrix.
/// \param n
/// \param nnz
/// \param Ap
/// \param Ai
/// \param stype
/// \param k
/// \return
 CSC *coarsen_k_times(int n, int nnz, int *Ap, int *Ai, int stype, int k);
}
#endif //FUSION_SPARSE_UTILITIES_H
