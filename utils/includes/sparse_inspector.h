//
// Created by Kazem on 10/19/19.
//

#ifndef PROJECT_SPARSE_INSPECTOR_H
#define PROJECT_SPARSE_INSPECTOR_H

#include "def.h"

namespace sym_lib {
 /// ========================== Inspection ===========================

 /// Creates levelset from a CSC matrix,
 /// the matrix should form DAG i.e., no cycle. It returns -1 if
 /// the input graph has a cycle.
 /// \param n
 /// \param nnz
 /// \param Lp
 /// \param Li
 /// \param levelPtr
 /// \param levelSet
 /// \return
 int build_levelSet_CSC(size_t n, int *Lp, int *Li,
                        int *&levelPtr, int *&levelSet);

 ///
 /// \param n
 /// \param nnz
 /// \param Lp
 /// \param Li_ptr
 /// \param Li
 /// \param blockNo
 /// \param sup2col
 /// \param col2sup
 /// \param levelPtr
 /// \param levelSet
 /// \return
 int build_levelSet_BCSC(int n, int nnz, size_t *Lp, size_t *Li_ptr, int *Li,
                         size_t blockNo, const int *sup2col, const int *col2sup,
                         int *&levelPtr, size_t *&levelSet);

 ///
 /// \param n
 /// \param Lp
 /// \param Li
 /// \param levelptr
 /// \param levels
 /// \param n_kern
 /// \return
 int
 level_set_multi_graphs(int n, const int *Lp, const int *Li, int *&levelptr,
                        int *&levels, int n_kern);



 /// ========================= Block Format ============================



 /**
* Generates the supernodes / blocks and matrix L in
* new blocked structure based on the DAG of dependency
*
* @param n column rank of matrix L
* @param nz number non-zero of matrix L
* @param Lp column pointers of matrix L
* @param Li row indices of matrix L
* @param Lx values of matrix L
* @param supernodes_ supernodes to be created
* @param num_nodes_ number of supernodes
* @param newLp_ column pointers of blocked L
* @param newLi_ row indices of blocked L
* @param nrows_ number of rows in each block
* @param newLx_ new values of blocked L
*/
 void block_gen_lsolve(int n, int *Lp, int *Li, double *Lx,
                       int *&supernodes, int &num_nodes, int *&newLp,
                       int *&newLi, int *&nrows, double *&newLx);

 ///
 /// \param n
 /// \param Lp
 /// \param Li
 /// \param levelptr
 /// \param levels
 /// \param supernodes
 /// \param num_nodes
 /// \param n_kern
 /// \param isLower
 /// \return
 int
 level_set_bn(int n, const int *Lp, const int *Li, int *&levelptr, int *&levels,
              int *supernodes, int num_nodes, int n_kern, bool isLower);


 /// Creates level-set from a tree
 /// \param n
 /// \param inTree : array of parents
 /// \param levelPtr
 /// \param levelSet
 /// \return
 int build_level_set_tree(size_t n, const int *inTree, int *levelPtr,
                          int *levelSet);
 int build_level_set_tree_efficient(size_t n, const int *tree,
                                    const int *nChild1,
                                    int *levelPtr, int *levelSet,
                                    int *node2level);
/**
 * Finds a node in list of supernodes
 *
 * @param supernodes list of supernodes
 * @param start the starting index to search
 * @param end the ending index to search
 * @param target target node to search
 * @return the index of the supernode in the list
 */
 int find(const int *supernodes, int start, int end, int target);

 void sup2node_gen(int n, int size,
                   const int *levels, int num_nodes, int *supernodes,
                   int *&sup2node);


 /// Merges BCSC L matrix with CSC A matrix and BCSC L matrix into level set
 /// \param A first matrix in CSC
 /// \param B second matrix in CSC
 /// \param C third matrix in CSC
 /// \param levelptr resulting level pointer
 /// \param levelset resulting level set
 /// \param nodes the number of supernodes in A and C
 /// \param supernodes the supernodes in A and C
 /// \return
 int bcsc_csc_bcsc_levelset(CSC *A, CSC *B, CSC *C, int *&levelptr, int *&levelset,
                         int nodes, int *supernodes);

 /// Extract supernodal CSC graph from (BCSC,CSC,BCSC) fusion graph
 /// \param A first matrix in BCSC
 /// \param B second matrix in CSC
 /// \param C third matrix in BCSC
 /// \return the extracted supernodal CSC
 CSC *merge_graphs(BCSC *A, CSC *B, BCSC *C);

 /// Extract a supernodal represented CSC from BCSC format matrix
 /// \param A BCSC matrix to get compressed supernodal CSC
 /// \return
 CSC *compressed_BCSC_to_CSC(BCSC *A);

 /// Computes the depth of each node in A by doing BFS in reversed DAG
 /// \param A
 /// \param degrees
 void
 compute_depth(CSC *A, int *degrees);


 /// =========================== ETREE COMPUTATION ===============================

///
/// \param k
/// \param i
/// \param Parent
/// \param Ancestor
 void update_etree(
   // inputs, not modified
   int k,  // process the edge (k,i) in the input graph
   int i,
   // inputs, modified on output
   int Parent[], // Parent [t] = p if p is the parent of t
   int Ancestor[] // Ancestor [t] is the ancestor of node t in the
   //partially-constructed etree
 );

/// Find the elimination tree of A or A'*A
/// \param A
/// \param Parent
/// \return
 int compute_etree(
   // input
   CSC *A,
   // output
   int *Parent // size ncol.  Parent [j] = p if p is the parent of j
 );


 ///
 /// \param p
 /// \param k
 /// \param Post
 /// \param Head
 /// \param Next
 /// \param Pstack
 /// \return
 int dfs_tree( int p, int k, int *Post, int *Head, int *Next,int *Pstack);

/// Takes the pattern of a CSC matrix and finds the pattern of row k
/// find nonzero pattern of Cholesky L(k,1:k-1) using etree and triu(A(:,k))
/// \param n in
/// \param Ap in
/// \param Ai in
/// \param k in row number
/// \param parent in etree
/// \param s output size of n
/// \param w output size of n
/// \return
 int ereach(int n, int *Ap, int *Ai, int k, const int *parent,
            int *s, int *w);


}
#endif //PROJECT_SPARSE_INSPECTOR_H
