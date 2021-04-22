//
// Created by kazem on 2/7/20.
//

#include "includes/sparse_inspector.h"
#include <utils.h>

namespace sym_lib{
 void update_etree(
     // inputs, not modified
     int k,		// process the edge (k,i) in the input graph
     int i,
     // inputs, modified on output
     int Parent [ ],	// Parent [t] = p if p is the parent of t
     int Ancestor [ ]	// Ancestor [t] is the ancestor of node t in the
			   //partially-constructed etree
   ){
  int a ;
  // traverse the path from k to the root of the tree
  for ( ; ; ){
   a = Ancestor [k] ;
   if (a == i){
    // final ancestor reached; no change to tree
    return ;
   }
   // perform path compression
   Ancestor [k] = i ;
   if (a == EMPTY){
    // final ancestor undefined; this is a new edge in the tree
    Parent [k] = i ;
    return ;
   }
   // traverse up to the ancestor of k
   k = a ;
  }
 }


 int compute_etree(
   // input
   CSC *A,
   // output
   int *Parent	// size ncol.  Parent [j] = p if p is the parent of j
 ) {
  int *Ap, *Ai, *Anz, *Ancestor, *Prev, *Iwork ;
  int i, j, jprev, p, pend, nrow, ncol, packed, stype ;
  size_t s ;
  int ok = true;
  stype = A->stype ;
  s = add_size_t (A->m, (stype ? 0 : A->n), &ok) ;
  if (!ok){
   //ERROR (CHOLMOD_TOO_LARGE, "problem too large") ;
   return false;
  }
  Iwork = new int[s]();

  // get inputs
  ncol = A->n ;	/* the number of columns of A */
  nrow = A->m ;	/* the number of rows of A */
  Ap = A->p ;		/* size ncol+1, column pointers for A */
  Ai = A->i ;		/* the row indices of A */
  // TODO: following two variables should be removed
  Anz = NULLPNTR ;	/* number of nonzeros in each column of A */
  packed = true; //A->packed ;
  Ancestor = Iwork ;	/* size ncol (i/i/l) */
#pragma omp parallel for
  for (j = 0 ; j < ncol ; j++){
   Parent [j] = EMPTY ;
   Ancestor [j] = EMPTY ;
  }
  /// compute the etree
  if (stype > 0){// symmetric (upper) case: compute etree (A)
   for (j = 0 ; j < ncol ; j++){
    // for each row i in column j of triu(A), excluding the diagonal
    p = Ap [j] ;
    pend = (packed) ? (Ap [j+1]) : (p + Anz [j]) ;
    for ( ; p < pend ; p++){
     i = Ai [p] ;
     if (i < j){
      update_etree (i, j, Parent, Ancestor) ;
     }
    }
   }
  }
  else if (stype == 0){
   /// unsymmetric case: compute etree (A'*A)
   Prev = Iwork + ncol ;	// size nrow (i/i/l)
   for (i = 0 ; i < nrow ; i++){
    Prev [i] = EMPTY ;
   }
   for (j = 0 ; j < ncol ; j++){
    // for each row i in column j of A
    p = Ap [j] ;
    pend = (packed) ? (Ap [j+1]) : (p + Anz [j]) ;
    for ( ; p < pend ; p++){
     /* a graph is constructed dynamically with one path per row
      * of A.  If the ith row of A contains column indices
      * (j1,j2,j3,j4) then the new graph has edges (j1,j2), (j2,j3),
      * and (j3,j4).  When at node i of this path-graph, all edges
      * (jprev,j) are considered, where jprev<j */
     i = Ai [p] ;
     jprev = Prev [i] ;
     if (jprev != EMPTY){
      update_etree (jprev, j, Parent, Ancestor) ;
     }
     Prev [i] = j ;
    }
   }
  }
  else{// symmetric case with lower triangular part not supported
   return false ;
  }
  delete []Iwork;
  return true ;
 }
}
