//
// Created by kazem on 6/11/20.
//

#include <sparse_utilities.h>
#include <metis_interface.h>
#include <sys/param.h>
#include <utils.h>
#include <sparse_io.h>
#include "lfactor_creation.h"
#include "def.h"
#include "sparse_inspector.h"


namespace sym_lib{

 int tree_post_order (int *Parent, size_t n, int *Weight, int *Post) {
  int *Head, *Next, *Pstack, *Iwork;
  int j, p, k, w, nextj;
  size_t s;
  s = 3*n;
  assert(s>0);
  Iwork = new int[s + 1]();
  Head = Iwork;
  Next = Iwork + n + 1;
  Pstack = Iwork + 2 * n + 1;
  for (int i = 0; i < n + 1; ++i) {
   Head[i] = EMPTY;
  }
  /* construct a link list of children for each node */
  if (Weight == NULL) {
   /* in reverse order so children are in ascending order in each list */
   for (j = n - 1; j >= 0; j--) {
    p = Parent[j];
    if (p >= 0 && p < ((int) n)) {
     /* add j to the list of children for node p */
     Next[j] = Head[p];
     Head[p] = j;
    }
   }
   /* Head [p] = j if j is the youngest (least-numbered) child of p */
   /* Next [j1] = j2 if j2 is the next-oldest sibling of j1 */
  } else {
   /* First, construct a set of link lists according to Weight.
    * Whead [w] = j if node j is the first node in bucket w.
    * Next [j1] = j2 if node j2 follows j1 in a link list.
    */
   int *Whead = Pstack;     /* use Pstack as workspace for Whead [ */
   for (w = 0; w < ((int) n); w++) {
    Whead[w] = EMPTY;
   }
   /* do in forward order, so nodes that ties are ordered by node index */
   for (j = 0; j < ((int) n); j++) {
    p = Parent[j];
    if (p >= 0 && p < ((int) n)) {
     w = Weight[j];
     w = MAX (0, w);
     w = MIN (w, ((int) n) - 1);
     /* place node j at the head of link list for weight w */
     Next[j] = Whead[w];
     Whead[w] = j;
    }
   }
   /* traverse weight buckets, placing each node in its parent's list */
   for (w = n - 1; w >= 0; w--) {
    for (j = Whead[w]; j != EMPTY; j = nextj) {
     nextj = Next[j];
     /* put node j in the link list of its parent */
     p = Parent[j];
     assert (p >= 0 && p < ((int) n));
     Next[j] = Head[p];
     Head[p] = j;
    }
   }
   /* Whead no longer needed ] */
   /* Head [p] = j if j is the lightest child of p */
   /* Next [j1] = j2 if j2 is the next-heaviest sibling of j1 */
  }
  /* start a DFS at each root node of the etree */
  k = 0;
  for (j = 0; j < ((int) n); j++) {
   if (Parent[j] == EMPTY) {
    /* j is the root of a tree; start a DFS here */
    k = dfs_tree(j, k, Post, Head, Next, Pstack);
   }
  }
  /* this would normally be EMPTY already, unless Parent is invalid */
  for (j = 0; j < ((int) n); j++) {
   Head[j] = EMPTY;
  }
  delete[]Iwork;
  return k;
 }


/* === initialize_node ====================================================== */
 int initialize_node  /* initial work for kth node in postordered etree */
   (
     int k,  /* at the kth step of the algorithm (and kth node) */
     int Post[], /* Post [k] = i, the kth node in postordered etree */
     int Parent[], /* Parent [i] is the parent of i in the etree */
     int ColCount[], /* ColCount [c] is the current weight of node c */
     int PrevNbr[] /* PrevNbr [u] = k if u was last considered at step k */
   ) {
  int p, parent;
  /* determine p, the kth node in the postordered etree */
  p = Post[k];
  /* adjust the weight if p is not a root of the etree */
  parent = Parent[p];
  if (parent != EMPTY) {
   ColCount[parent]--;
  }
  /* flag node p to exclude self edges (p,p) */
  PrevNbr[p] = k;
  return (p);
 }


/* === process_edge ========================================================= */
/* edge (p,u) is being processed.  p < u is a descendant of its ancestor u in
 * the etree.  node p is the kth node in the postordered etree.  */
 void process_edge
   (
     int p,  /* process edge (p,u) of the matrix */
     int u,
     int k,  /* we are at the kth node in the postordered etree */
     int First[], /* First [i] = k if the postordering of first
			 * descendent of node i is k */
     int PrevNbr[], /* u was last considered at step k = PrevNbr [u] */
     int ColCount[], /* ColCount [c] is the current weight of node c */
     int PrevLeaf[], /* s = PrevLeaf [u] means that s was the last leaf
			 * seen in the subtree rooted at u.  */
     int RowCount[], /* RowCount [i] is # of nonzeros in row i of L,
			 * including the diagonal.  Not computed if NULL. */
     int SetParent[], /* the FIND/UNION data structure, which forms a set
			 * of trees.  A root i has i = SetParent [i].  Following
			 * a path from i to the root q of the subtree containing
			 * i means that q is the SetParent representative of i.
			 * All nodes in the tree could have their SetParent
			 * equal to the root q; the tree representation is used
			 * to save time.  When a path is traced from i to its
			 * root q, the path is re-traversed to set the SetParent
			 * of the whole path to be the root q. */
     int Level[]  /* Level [i] = length of path from node i to root */
   ) {
  int prevleaf, q, s, sparent;
  if (First[p] > PrevNbr[u]) {
   /* p is a leaf of the subtree of u */
   ColCount[p]++;
   prevleaf = PrevLeaf[u];
   if (prevleaf == EMPTY) {
    /* p is the first leaf of subtree of u; RowCount will be incremented
     * by the length of the path in the etree from p up to u. */
    q = u;
   } else {
    /* q = FIND (prevleaf): find the root q of the
     * SetParent tree containing prevleaf */
    for (q = prevleaf; q != SetParent[q]; q = SetParent[q]) { ;
    }
    /* the root q has been found; re-traverse the path and
     * perform path compression */
    s = prevleaf;
    for (s = prevleaf; s != q; s = sparent) {
     sparent = SetParent[s];
     SetParent[s] = q;
    }
    /* adjust the RowCount and ColCount; RowCount will be incremented by
     * the length of the path from p to the SetParent root q, and
     * decrement the ColCount of q by one. */
    ColCount[q]--;
   }
   if (RowCount != NULL) {
    /* if RowCount is being computed, increment it by the length of
     * the path from p to q */
    RowCount[u] += (Level[p] - Level[q]);
   }
   /* p is a leaf of the subtree of u, so mark PrevLeaf [u] to be p */
   PrevLeaf[u] = p;
  }
  /* flag u has having been processed at step k */
  PrevNbr[u] = k;
 }


/* === finalize_node ======================================================== */
 void finalize_node    /* compute UNION (p, Parent [p]) */
   (
     int p,
     int Parent[], /* Parent [p] is the parent of p in the etree */
     int SetParent[] /* see process_edge, above */
   ) {
  /* all nodes in the SetParent tree rooted at p now have as their final
   * root the node Parent [p].  This computes UNION (p, Parent [p]) */
  if (Parent[p] != EMPTY) {
   SetParent[p] = Parent[p];
  }
 }

  bool row_col_counts_symmetric(CSC *A, int *Parent, int *Post, int *RowCount,
                                int *ColCount, int *First, int *Level,
                                size_t &fl){
  double ff;
  int *Ap, *Ai, *PrevNbr, *SetParent, *Head, *PrevLeaf,
    *Iwork;
  int i, j, r, k, len, s, p, pend, stype, anz, parent,
    nrow, ncol;
  size_t w;
  int ok = true;
  stype = A->stype;
  if (stype >= 0)
   return false;
  nrow = A->m; /* the number of rows of A */
  ncol = A->n; /* the number of columns of A */
  /* w = 2*nrow + (stype ? 0 : ncol) */
  w = mult_size_t(nrow, 2, &ok);
  w = add_size_t(w, (stype ? 0 : ncol), &ok);
  if (!ok)
   return false;
  Ap = A->p; /* size ncol+1, column pointers for A */
  Ai = A->i; /* the row indices of A, of size nz=Ap[ncol+1] */
  Iwork = new int[w]();
  SetParent = Iwork;      /* size nrow (i/i/l) */
  PrevNbr = Iwork + nrow;     /* size nrow (i/i/l) */
  PrevLeaf = new int[nrow];
  Head = new int[nrow + 1];
  /* find the first descendant and level of each node in the tree */
  /* First [i] = k if the postordering of first descendent of node i is k */
  /* Level [i] = length of path from node i to the root (Level [root] = 0) */
  for (i = 0; i < nrow; i++) {
   First[i] = EMPTY;
  }
  /* postorder traversal of the etree */
  for (k = 0; k < nrow; k++) {
   /* node i of the etree is the kth node in the postordered etree */
   i = Post[k];
   /* i is a leaf if First [i] is still EMPTY */
   /* ColCount [i] starts at 1 if i is a leaf, zero otherwise */
   ColCount[i] = (First[i] == EMPTY) ? 1 : 0;
   /* traverse the path from node i to the root, stopping if we find a
    * node r whose First [r] is already defined. */
   len = 0;
   for (r = i; (r != EMPTY) && (First[r] == EMPTY); r = Parent[r]) {
    First[r] = k;
    len++;
   }
   if (r == EMPTY) {
    /* we hit a root node, the level of which is zero */
    len--;
   } else {
    /* we stopped at node r, where Level [r] is already defined */
    len += Level[r];
   }
   /* re-traverse the path from node i to r; set the level of each node */
   for (s = i; s != r; s = Parent[s]) {
    Level[s] = len--;
   }
  }
  fl = 0.0;
  /* compute the row counts and node weights */
  if (RowCount != NULL) {
   for (i = 0; i < nrow; i++) {
    RowCount[i] = 1;
   }
  }
  for (i = 0; i < nrow; i++) {
   PrevLeaf[i] = EMPTY;
   PrevNbr[i] = EMPTY;
   SetParent[i] = i; /* every node is in its own set, by itself */
  }
   /* also determine the number of entries in triu(A) */
   anz = nrow;
   for (k = 0; k < nrow; k++) {
    /* j is the kth node in the postordered etree */
    j = initialize_node(k, Post, Parent, ColCount, PrevNbr);
    /* for all nonzeros A(i,j) below the diagonal, in column j of A */
    p = Ap[j];
    pend = Ap[j + 1];
    for (; p < pend; p++) {
     i = Ai[p];
     if (i > j) {
      /* j is a descendant of i in etree(A) */
      anz++;
      process_edge(j, i, k, First, PrevNbr, ColCount,
                   PrevLeaf, RowCount, SetParent, Level);
     }
    }
    /* update SetParent: UNION (j, Parent [j]) */
    finalize_node(j, Parent, SetParent);
   }
  for (j = 0; j < nrow; j++) {
   parent = Parent[j];
   if (parent != EMPTY) {
    /* add the ColCount of j to its parent */
    ColCount[parent] += ColCount[j];
   }
  }
  /* flop count and nnz(L) for subsequent LL' numerical factorization */
  /* use double to avoid integer overflow.  lnz cannot be NaN. */
  fl = 0;
  for (j = 0; j < nrow; j++) {
   ff = (double) (ColCount[j]);
   fl += ff * ff;
  }
  delete[]Iwork;
  delete[]PrevLeaf;
  delete[]Head;
  return true;
 }


 double compute_cost_per_col(int n, int colNo, int *eTree,
                          int *cT, int *rT, int *xi, int &top){
  //int *xi = new int[2*n]();
  double total=0.0;
  top = ereach(n, cT, rT, colNo, eTree, xi, xi+n);
  for (int i = top; i < n; ++i) {
   int spCol = xi[i];
   bool sw=false;
   int ub = 0;
   total += cT[spCol + 1] - cT[spCol];
  }
  //std::cout<<total<<",";
  total += cT[colNo + 1] - cT[colNo];
  return total;
 }

 CSC *build_L_pattern_from_col_counts(int n, CSC *A,
   int *col_count, int *parent, double *accessed_nnz){
  int *Lp = new int[n + 1]();
  Lp[0] = 0;
  for (int l = 1; l < n + 1; ++l) {
   Lp[l] = Lp[l - 1] + col_count[l - 1];
  }
  int nnz = Lp[n];
  int *Li = new int[nnz];
  int top = 0;
  int *xi = new int[2 * n]();
  int *finger = new int[n]();
  //copying original row idx to L
  for (int s = 0; s < n; ++s) {
   //Populating row indices of L
    accessed_nnz[s]=compute_cost_per_col(n, s, parent,
                                  A->p, A->i, xi, top);
   Li[Lp[s]] = s;
   finger[s]++;
   for (int i = top; i < n; ++i) {
    int col_beg = xi[i];
    assert(col_beg < n);
    assert(Lp[col_beg] + finger[col_beg] < nnz);
    Li[Lp[col_beg] + finger[col_beg]] = s;
    finger[col_beg]++;
   }
  }
  CSC *L = new CSC(n,n,nnz,Lp,Li,NULLPNTR);
  L->x = new double[nnz]();
  L->pre_alloc = false;
  delete[]xi;
  delete[]finger;
  return L;
 }


 void reorder_tree(int n, int *parent, int *Post,
                   int *ipost, int *Wi){
  int oldchild, oldparent, newparent;
  for (int newchild = 0; newchild < n; newchild++) {
   oldchild = Post[newchild];
   oldparent = parent[oldchild];
   newparent = (oldparent == EMPTY) ? EMPTY : ipost[oldparent];
   Wi[newchild] = newparent;
  }
  for (int k = 0; k < n; k++) {
   parent[k] = Wi[k];
  }
 }

 LFactorSymbolicInfo *build_symbolic_simplicial_lfactor(CSC *A){
  if(!A)
   return NULLPNTR;
  int n = A->n;
  LFactorSymbolicInfo *lfsi = new LFactorSymbolicInfo(n);
  int *perm, *parent, *post, *col_count, *iperm, *ipost;
  perm = lfsi->perm; iperm = lfsi->iperm; post = lfsi->post;
  ipost = lfsi->ipost; parent = lfsi->parent;
  col_count=lfsi->col_counts;
  int *first = new int[n](), *level = new int[n]();
  int *ws = new int[n];
  int fl, aatfl, lnz;
  /// compute fill-in reducing permutation
#ifdef METIS
  int *perm_metis;
  CSC *A_full = make_full(A);
  metis_perm_general(A_full, perm_metis);
  CSC *At_ord = transpose_symmetric(A, perm_metis);
  CSC *A_ord = transpose_symmetric(At_ord, NULLPNTR);
  std::copy(perm_metis,perm_metis+n, perm);
  delete perm_metis;
  delete A_full;
#else
  for(int i=0; i<n; i++) perm[i] = i; // no permutation
#endif
  /// Building etree
  compute_etree(At_ord, parent);
  int ret = tree_post_order(parent, n, NULLPNTR, post);
  assert(ret == n);
  /// computing column count
  row_col_counts_symmetric(A_ord, parent, post, NULLPNTR, col_count,
    first, level,lfsi->flops);
  ret = tree_post_order(parent, n, col_count, post);
  assert(ret == n);
  /// Combining post with perm
  reorder_array(n, perm, post, ws);
  reorder_array(n, col_count, post, ws);
  inv_perm(n, post, ipost);
  inv_perm(n, perm, iperm);
  reorder_tree(n, parent, post, ipost, ws);
  delete At_ord;
  At_ord = transpose_symmetric(A, perm);
  /// Build the L-factor pattern
  lfsi->L = build_L_pattern_from_col_counts(n, At_ord, col_count, parent,
    lfsi->accessed_nnz_per_col);
  delete A_ord;
  delete At_ord;
  delete []ws;
  delete []first;
  delete []level;
  return lfsi;
 }

}
