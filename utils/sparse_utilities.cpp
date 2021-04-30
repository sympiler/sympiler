//
// Created by kazem on 10/9/19.
//

#include <omp.h>
#include <cstring>
#include <def.h>
#include <cassert>
#include <cmath>
#include "sparse_utilities.h"

namespace sym_lib {
 using namespace std;

 void save_vector(char *fname, const double *vec, int length) {
  FILE *f = fopen(fname, "w");
  for(int i = 0; i < length; i++)
   fprintf(f, "%f\n", vec[i]);
  fclose(f);
 }

 // TODO eventually use CSC / CSR class
 // TODO test for these two format change
 int
 csc_to_csr(int nrow, int ncol, int *Ap, int *Ai, double *Ax, int *&rowptr,
            int *&colind, double *&values) {
  // count row entries to generate row ptr
  int nnz = Ap[ncol];
  int *rowCnt = new int[nrow]();
  for (int i = 0; i < nnz; i++)
   rowCnt[Ai[i]]++;

  rowptr = new int[nrow + 1]();
  int counter = 0;
  for (int i = 0; i < nrow; i++) {
   rowptr[i] = counter;
   counter += rowCnt[i];
  }
  rowptr[nrow] = nnz;

  colind = new int[nnz]();
  values = new double[nnz]();

  memset(rowCnt, 0, sizeof(int) * nrow);
  for (int i = 0; i < ncol; i++) {
   for (int j = Ap[i]; j < Ap[i + 1]; j++) {
    int row = Ai[j];
    int index = rowptr[row] + rowCnt[row];
    colind[index] = i;
    values[index] = Ax[j];
    rowCnt[row]++;
   }
  }
  delete[]rowCnt;

  return 0;
 }

 CSR* csc_to_csr(CSC* A) {
  // count row entries to generate row ptr
  int nnz = A->p[A->n];
  int *rowCnt = new int[A->n]();
  for (int i = 0; i < nnz; i++)
   rowCnt[A->i[i]]++;

  CSR *B = new CSR(A->n,A->m,A->nnz,A->is_pattern);
  int *rowptr = B->p; //new int[nrow + 1]();
  size_t ncol = B->n;
  size_t nrow = B->m;
  int counter = 0;
  for (int i = 0; i < (int)nrow; i++) {
   rowptr[i] = counter;
   counter += rowCnt[i];
  }
  rowptr[nrow] = nnz;

  int *colind = B->i;
  double *values = B->x;

  memset(rowCnt, 0, sizeof(int) * nrow);
  for (int i = 0; i < (int)ncol; i++) {
   for (int j = A->p[i]; j < A->p[i + 1]; j++) {
    int row = A->i[j];
    int index = rowptr[row] + rowCnt[row];
    colind[index] = i;
    if(!B->is_pattern)
     values[index] = A->x[j];
    rowCnt[row]++;
   }
  }
  delete[]rowCnt;
  return B;
 }

 int csr_to_csc(int nrow, int ncol, int *Ai, int *Ap, double *Ax, int *&colptr,
                int *&rowind, double *&values) {
  return csc_to_csr(ncol, nrow, Ap, Ai, Ax, colptr, rowind, values);
 }

 CSC* csr_to_csc(CSR *A) {
  CSC *B = new CSC(A->m,A->n,A->nnz,NULLPNTR,NULLPNTR,A->stype);
  B->pre_alloc = false;
  B->is_pattern = false;
  int st = csc_to_csr(A->m, A->n, A->p, A->i, A->x, B->p, B->i, B->x);
  B->nnz = B->p[B->n];
  return B;
 }

 CSC *tridiag(int n, double a0, double a1, double a2) {
  CSC *A = new CSC(n, n, 3 * n - 2);
  A->is_pattern = true;
  A->stype = -1;

  int ind = 0;
  A->p[0] = 0;
  for(int i = 0; i < n; i++) {
   if(i == 0 || i == n-1) {
    A->p[i+1] = A->p[i] + 2;
    if(i == 0) {
     A->i[ind] = i;
     A->x[ind] = a1;
     A->i[ind+1] = i+1;
     A->x[ind+1] = a2;
    } else {
     A->i[ind] = i-1;
     A->x[ind] = a0;
     A->i[ind+1] = i;
     A->x[ind+1] = a1;
    }
    ind += 2;
   } else {
    A->p[i+1] = A->p[i] + 3;
    A->i[ind + 0] = i-1;
    A->i[ind + 1] = i;
    A->i[ind + 2] = i+1;
    A->x[ind + 0] = a0;
    A->x[ind + 1] = a1;
    A->x[ind + 2] = a2;
    ind += 3;
   }
  }

  return A;
 }


 CSC* tree_to_csc(int n, int *tree){
  int *Ai = new int[n](), *Ap = new int[n+1](), *cnt = new int[n]();
  populate_children(n, tree, Ap, Ai, cnt);
  CSC *T = new CSC(n,n,n-1,Ap,Ai,0);
  delete []cnt;
  return T;
 }

 CSC* transpose_general(CSC *A){
  size_t row = A->m;
  size_t col = A->n;
  size_t  colT = row;
  CSC *AT=NULLPNTR;
  if(row == 0){
   return AT;
  }
  AT = new CSC(col,row,A->nnz,A->is_pattern);
  int *col_cnt = new int[row]();//cols will be rows
  for (int i = 0; i < A->p[col]; ++i) {
   col_cnt[A->i[i]]++;
  }
  AT->p[0]=0;
  //Determining column pointer of AT
  for (int j = 1; j < (int)colT+1; ++j) {
   AT->p[j] = AT->p[j-1] + col_cnt[j-1];
  }
  std::fill_n(col_cnt,colT,0);
  //Determining row pointer of AT
  for (int k = 0; k < (int)col; ++k) {
   for (int i = A->p[k]; i < A->p[k+1]; ++i) {
    int beg = A->i[i];
    assert(AT->p[beg]+col_cnt[beg] < A->nnz);
    AT->i[AT->p[beg]+col_cnt[beg]] = k;
    if(!AT->is_pattern)
     AT->x[AT->p[beg]+col_cnt[beg]] = A->x[i];
    col_cnt[beg]++;
   }
  }
  delete []col_cnt;
  return AT;
 }


 CSC* make_half(size_t An, int *Ap, int *Ai, double *Ax, bool lower){
  size_t Bn = An;
  size_t Bnz = 0;
  for (int i = 0; i < (int)An; ++i) {
   for (int j = Ap[i]; j < Ap[i + 1]; ++j) {
    bool l_u = lower ? Ai[j]>=i : Ai[j]<=i;
    if(l_u){
     Bnz ++;
    }
   }
  }
  CSC *B = new CSC(Bn,Bn,Bnz,Ax==NULLPNTR);
  int *Bp = B->p;
  int *Bi = B->i;
  double *Bx = B->x;
  int nnz_cnt=0;
  for (int i = 0; i < An; ++i) {
   for (int j = Ap[i]; j < Ap[i + 1]; ++j) {
    bool l_u = lower ? Ai[j]>=i : Ai[j]<=i;
    if(l_u){
     Bi[nnz_cnt] = Ai[j];
     if(Ax!=NULLPNTR) Bx[nnz_cnt] = Ax[j];
     nnz_cnt++;
    }
   }
   Bp[i+1] = nnz_cnt;
  }
  B->stype = -1;
  return B;
 }

 CSC *make_full(CSC *A) {
  if(A->stype == 0) {
   std::cerr << "Not symmetric\n";
   return nullptr;
  }

  CSC *Afull = new CSC(A->m, A->n, A->nnz * 2 - A->n);
  auto ind = new int[A->n]();

  for(size_t i = 0; i < A->n; i++) {
   for(size_t p = A->p[i]; p < A->p[i+1]; p++) {
    int row = A->i[p];
    ind[i]++;
    if(row != i)
     ind[row]++;
   }
  }
  Afull->p[0] = 0;
  for(size_t i = 0; i < A->n; i++)
   Afull->p[i+1] = Afull->p[i] + ind[i];

  for(size_t i = 0; i < A->n; i++)
   ind[i] = 0;
  for(size_t i = 0; i < A->n; i++) {
   for(size_t p = A->p[i]; p < A->p[i+1]; p++) {
    int row = A->i[p];
    int index = Afull->p[i] + ind[i];
    Afull->i[index] = row;
    Afull->x[index] = A->x[p];
    ind[i]++;
    if(row != i) {
     index = Afull->p[row] + ind[row];
     Afull->i[index] = i;
     Afull->x[index] = A->x[p];
     ind[row]++;
    }
   }
  }
  delete[]ind;

  return Afull;
 }


 CSC *transpose_symmetric(CSC *A, int *Perm) {
  int *Ap, *Ai, *Fi, *Fp;
  double *Ax, *Fx;
  int  *Wi, *Pinv, *Iwork;
  int p, pend, upper, permute, jold, i, j, k, iold, fp;
  size_t n;
  size_t s;
  CSC *F;
  int stype;
  stype = A->stype;
  int is_np = !A->is_pattern;
  if (stype == 0 || A->m != A->n) {
   return NULLPNTR;
  }
  Ap = A->p;
  Ai = A->i;
  Ax = A->x;
  n = A->m;
  permute = (Perm != NULLPNTR);
  upper = (A->stype > 0);
  F = new CSC(A->n,A->m,A->nnz);
  Fp = F->p;        /* size A->nrow+1, row pointers of F */
  Fi = F->i;
  Fx = F->x;
  s = n +  ((Perm != NULLPNTR) ? n : 0) ;
  Iwork = new int[s]();
  Wi = Iwork;
  Pinv = Iwork + n; // unused if Perm NULL

  //check Perm and construct inverse permutation
  if (permute) {
   for (i = 0; i < n; i++) {
    Pinv[i] = EMPTY;
   }
   for (k = 0; k < n; k++) {
    i = Perm[k];
    if (i < 0 || i > n || Pinv[i] != EMPTY) {
     return NULLPNTR;
    }
    Pinv[i] = k;
   }
  }

  // count the entries in each row of F
  for (i = 0; i < n; i++) {
   Wi[i] = 0;
  }
  if (permute) {
   if (upper) {
    /* packed, permuted, upper */
    for (j = 0; j < n; j++) {
     jold = Perm[j];
     pend = Ap[jold + 1];
     for (p = Ap[jold]; p < pend; p++) {
      iold = Ai[p];
      if (iold <= jold) {
       i = Pinv[iold];
       Wi[std::min(i, j)]++;
      }
     }
    }
   } else {
    /* packed, permuted, lower */
    for (j = 0; j < n; j++) {
     jold = Perm[j];
     pend = Ap[jold + 1];
     for (p = Ap[jold]; p < pend; p++) {
      iold = Ai[p];
      if (iold >= jold) {
       i = Pinv[iold];
       Wi[std::max(i, j)]++;
      }
     }
    }
   }
  } else {
   if (upper) {
    /* packed, unpermuted, upper */
    for (j = 0; j < n; j++) {
     pend = Ap[j + 1];
     for (p = Ap[j]; p < pend; p++) {
      i = Ai[p];
      if (i <= j) {
       Wi[i]++;
      }
     }
    }
   } else {
    /* packed, unpermuted, lower */
    for (j = 0; j < n; j++) {
     pend = Ap[j + 1];
     for (p = Ap[j]; p < pend; p++) {
      i = Ai[p];
      if (i >= j) {
       Wi[i]++;
      }
     }
    }
   }
  }
  // compute the row pointers
  p = 0;
  for (i = 0; i < n; i++) {
   Fp[i] = p;
   p += Wi[i];
  }
  Fp[n] = p;
  for (i = 0; i < n; i++) {
   Wi[i] = Fp[i];
  }
  if (p > F->nnz) {
   return NULLPNTR;
  }

/// Transpose the matrix
  if (permute) {
   if (upper) {
    /* permuted, upper */
    for (j = 0; j < n; j++) {
     jold = Perm[j];
     p = Ap[jold];
     pend = Ap[jold + 1];
     for (; p < pend; p++) {
      iold = Ai[p];
      if (iold <= jold) {
       i = Pinv[iold];
       if (i < j) {
        fp = Wi[i]++;
        Fi[fp] = j;
        if(is_np)
         Fx[fp] = Ax[p];
       } else {
        fp = Wi[j]++;
        Fi[fp] = i;
        if(is_np)
         Fx[fp] = Ax[p];
       }
      }
     }
    }
   } else {
    /* permuted, lower */
    for (j = 0; j < n; j++) {
     jold = Perm[j];
     p = Ap[jold];
     pend = Ap[jold + 1];
     for (; p < pend; p++) {
      iold = Ai[p];
      if (iold >= jold) {
       i = Pinv[iold];
       if (i > j) {
        fp = Wi[i]++;
        Fi[fp] = j;
        if(is_np)
         Fx[fp] = Ax[p];
       } else {
        fp = Wi[j]++;
        Fi[fp] = i;
        if(is_np)
         Fx[fp] = Ax[p];
       }
      }
     }
    }
   }
  } else {
   if (upper) {
    /* unpermuted, upper */
    for (j = 0; j < n; j++) {
     p = Ap[j];
     pend = Ap[j + 1];
     for (; p < pend; p++) {
      i = Ai[p];
      if (i <= j) {
       fp = Wi[i]++;
       Fi[fp] = j;
       if(is_np)
        Fx[fp] = Ax[p];
      }
     }
    }
   } else {
    /* unpermuted, lower */
    for (j = 0; j < n; j++) {
     p = Ap[j];
     pend = Ap[j + 1];
     for (; p < pend; p++) {
      i = Ai[p];
      if (i >= j) {
       fp = Wi[i]++;
       Fi[fp] = j;
       if(is_np)
        Fx[fp] = Ax[p];
      }
     }
    }
   }
  }
  F->stype = -A->stype; // flip the stype
  delete []Iwork;
  return F;
 }

 double residual(int beg_idx, int end_idx, const double *x0, const double *x1) {
  double max = 0.0;
  for (int i = beg_idx; i < end_idx; i++)
   if (std::fabs(x0[i] - x1[i]) > max)
    max = std::fabs(x0[i] - x1[i]);
  return max;
 }

 // FIXME: different types of norms
 double norm(int size, double *vec) {
  double max = 0.0;
  for(int i = 0; i < size; i++) {
   if(std::fabs(vec[i]) > max)
    max = std::fabs(vec[i]);
  }
  return max;
 }


 int *compute_inv_perm(int n, int *perm, int *pinv){
  if (n<=0)
   return NULL;
  for (int i = 0; i < n; ++i) {
   assert(perm[i]>=0 && perm[i]<n);
   pinv[perm[i]] = i;
  }
  return pinv;
 }


 CSC *copy_sparse(CSC *A){
  CSC *clone = new CSC(A->m,A->n,A->nnz);
  for (int i = 0; i < A->n+1; ++i) {
   clone->p[i] = A->p[i];
  }
  for (int j = 0; j < A->nnz; ++j) {
   clone->i[j] = A->i[j];
   clone->x[j] = A->x[j];
  }
  clone->stype = A->stype;
  return clone;
 }

 void copy_from_to(CSR *src, CSR *dst){
  for (int i = 0; i < src->n; ++i) {
   dst->p[i] = dst->p[i];
   for (int j = src->p[i]; j < src->p[i + 1]; ++j) {
    dst->i[j] = src->i[j];
    dst->x[j] = src->x[j];
   }
  }
 }


 void copy_from_to(CSC *src, CSC *dst){
  for (int i = 0; i < src->n; ++i) {
   dst->p[i] = dst->p[i];
   for (int j = src->p[i]; j < src->p[i + 1]; ++j) {
    dst->i[j] = src->i[j];
    dst->x[j] = src->x[j];
   }
  }
 }


 void copy_vector_dense(size_t beg, size_t end, const double *src, double *dst) {
  for(size_t i = beg; i < end; i++)
   dst[i] = src[i];
 }


 int number_empty_col(CSC *A){
  int cnt=0;
  for (int i = 0; i < A->n; ++i) {
   if(A->p[i+1]-A->p[i] == 0){
    cnt++;
   }
  }
  return cnt;
 }


 int modified_BFS_BCSC(int n, size_t *Lp, size_t *Li_ptr, int* Li,
                      const int* sup2col, const int* col2sup, int *inDegree,
                      bool *visited, int *node2partition, int* &levelPtr,
                      size_t* levelSet, int bfsLevel,
                      std::vector<std::vector<int>> &newLeveledParList){

  std::vector<int> queue;
  //Let's do BFS for every leaf node of the
  for (int ii = levelPtr[bfsLevel]; ii < levelPtr[bfsLevel+1]; ++ii) {
   int curNode=levelSet[ii];
   assert(node2partition[curNode]>=0);
   queue.push_back(curNode);
   while (!queue.empty()){
    int popedNode=queue[0];
    queue.erase(queue.begin());
    visited[popedNode]=true;
    newLeveledParList[node2partition[popedNode]].push_back(popedNode);
    //Find the adjacent nodes_
    int curCol = sup2col[popedNode];
    int nxtCol = sup2col[popedNode+1];
    int supWdt = nxtCol-curCol;
    for (int r = Li_ptr[curCol]+supWdt; r < Li_ptr[nxtCol]; ++r) {
     int cn=col2sup[Li[r]];
     inDegree[cn]--;
     if(inDegree[cn]==1 && !visited[cn]){
      queue.push_back(cn);
     }
    }
   }
  }

  return 1;// TODO:
 }


 int modified_BFS_CSC(int n, int *Lp, int *Li, int *inDegree, bool *visited,
                      int *node2partition, int *&levelPtr, int *levelSet,
                      int bfsLevel,
                      std::vector<std::vector<int>> &newLeveledParList){

  std::vector<int> queue;
  //Let's do BFS for every leaf node of the
  for (int ii = levelPtr[bfsLevel]; ii < levelPtr[bfsLevel+1]; ++ii) {
   int curNode=levelSet[ii];
   assert(!visited[curNode]);
   assert(node2partition[curNode]>=0);
   queue.push_back(curNode);
   while (!queue.empty()){
    int popedNode=queue[0];
    assert(popedNode < n);
    queue.erase(queue.begin());
    visited[popedNode]=true;
    assert(node2partition[popedNode] < newLeveledParList.size());
    newLeveledParList[node2partition[popedNode]].push_back(popedNode);
#if 0
    if(node2partition[popedNode] == 805){
     std::cout<<" "<<popedNode<<"\n";
    }
#endif
    //Find the adjacent nodes_
    for (int r = Lp[popedNode]; r < Lp[popedNode+1]; ++r) {
     int cn=Li[r];
     inDegree[cn]--;
     if(inDegree[cn]==1 && !visited[cn]){
      queue.push_back(cn);
     }
    }
   }
  }
  return 1; //TODO
 }


 void populate_children(int n, const int *eTree, int *childPtr,
                        int *childNo, int *nChild) {
  int *childCnt = new int[n]();
  for (int k = 0; k < n; ++k) {
   if (eTree[k] >= 0)
    nChild[eTree[k]]++;
  }
  childPtr[0] = 0;
  for (int k = 1; k < n + 1; ++k) {
   childPtr[k] = childPtr[k - 1] + nChild[k - 1];
  }
  for (int l = 0; l < n; ++l) {
   int p = eTree[l];
   if (p >= 0) {
    childNo[childPtr[p] + childCnt[p]] = l;
    childCnt[p]++;
   }
  }
  delete[]childCnt;
 }


 int get_node_depth(int node, int n, const int *tree, int *weight) {
  int level = 0;
  if (weight != NULLPNTR)
   level += weight[node];
  while (tree[node] >= 0) {
   node = tree[node];
   if (weight == NULLPNTR)
    level++;
   else
    level += weight[node];
  }
  return level;
 }


 int get_tree_height_bruteforce(int n, const int *tree, int *weight) {
  int maxLen = 0;
  for (int i = 0; i < n; ++i) {
   int ltmp = get_node_depth(i, n, tree, weight);
   if (ltmp > maxLen)
    maxLen = ltmp;
  }
  return maxLen;
 }


 int get_tree_height(int n, const int *tree, int *nChild1, int *weight) {
  int maxLen = 0;
  for (int i = 0; i < n; ++i) {
   if (nChild1[i] == 0) {
    int ltmp = get_node_depth(i, n, tree, weight);
    if (ltmp > maxLen)
     maxLen = ltmp;
   }
  }
  return maxLen;
 }


 int get_tree_height_efficient(int n, const int *tree, const int *nChild1,
   const int *weight) {
  int maxLen = 0, level = 0;
  auto *touch = new int[n]();
  for (int i = 0; i < n; ++i) {
   if (nChild1[i] == 0) {
    int node = i; level = 0;
    while (tree[node] >= 0) {
     node = tree[node];
     if (weight == NULLPNTR)
      level++;
     else
      level += weight[node];
     if(level < touch[node])
      break; // a longer path visted this node before
     else
      touch[node] = level;// touched this node with level
    }
    if (level > maxLen)
     maxLen = level;
   }
  }
  delete []touch;
  return maxLen;
 }

 double *compute_subtree_cost(int n, const int *tree, double *weight) {
  double *subTreeCost = new double[n];
  for (int i = 0; i < n; ++i) {//each subtree has a node
   subTreeCost[i] = weight[i];
  }
  for (int i = 0; i < n; ++i) {
   int par = tree[i];
   if (par > 0) {
    subTreeCost[par] += subTreeCost[i];
   }
  }
  return subTreeCost;
 }


 void sparse2dense(CSC *A, double *D){
  std::fill_n(D,A->n*A->m,0);
  for (int i = 0; i < A->n; ++i) {
   for (int j = A->p[i]; j < A->p[i+1]; ++j) {
    int r = A->i[j];
    D[i*A->m+r] = A->x[j];
    if(A->stype==-1 || A->stype==1){
     D[r*A->n+i] = A->x[j];
    }
   }
  }
 }


 CSC *diagonal(int n, double val){
  CSC *A = new CSC(n,n,n);
  A->p[0] = 0;
  for (int i = 0; i < n; ++i) {
   A->i[i] = i;
   A->p[i+1] = A->p[i] + 1;
   A->x[i] = val;
  }
  assert(A->p[n] == n);
  return A;
 }


 int *extract_diagonals(const int n, const int *Ap, const int *Ai){
  auto *ret_ptr = new int[n];
  for (int i = 0; i < n; ++i) {
   for (int j = Ap[i]; j < Ap[i+1]; ++j) {
    if(Ai[j] == i){
     ret_ptr[i] = j;
     break;
    }
   }
  }
  return ret_ptr;
 }


 void merge_graph(int ngraphs, int n, int **Gps, int **Gis,
   int *&nGp, int *&nGi) {
  const int *Gp, *Gi;

  int nnz = 0;
  for(int i = 0; i < ngraphs; i++)
   nnz += Gps[i][n];
  nnz += (ngraphs-1) * n;
  /** allocate new graph space **/
  nGp = new int[ngraphs * n + 1]();
  nGi = new int[nnz]();

  int p_counter = 1;
  int i_counter = 0;
  nGp[0] = 0;
  /** first <ngraphs>-1 graphs**/
  for(int i = 0; i < ngraphs-1; i++) {
   Gp = Gps[i];
   Gi = Gis[i];

   for(int j = 0; j < n; j++) {
    int diff = Gp[j+1] - Gp[j] + 1;
    nGp[p_counter] = nGp[p_counter-1] + diff;
    p_counter++;
    for(int p = Gp[j]; p < Gp[j+1]; p++) {
     int row = Gi[p] + (i * n);
     nGi[i_counter] = row;
     i_counter++;
    }
    nGi[i_counter] = j + (i+1) * n;
    i_counter++;
   }
  }
  /** last graph **/
  Gp = Gps[ngraphs-1];
  Gi = Gis[ngraphs-1];
  for(int j = 0; j < n; j++) {
   int diff = Gp[j+1] - Gp[j];
   nGp[p_counter] = nGp[p_counter-1] + diff;
   p_counter++;
   for(int p = Gp[j]; p < Gp[j+1]; p++) {
    int row = Gi[p] + (ngraphs-1) * n;
    nGi[i_counter] = row;
    i_counter++;
   }
  }
 }

 CSC* merge_graph(int ngraphs, int n, int **Gps, int **Gis) {
  const int *Gp, *Gi;
  int nnz = 0;
  for(int i = 0; i < ngraphs; i++)
   nnz += Gps[i][n];
  nnz += (ngraphs-1) * n;
  CSC *merged_graph = new CSC(n*ngraphs,n*ngraphs,nnz, true);
  /** allocate new graph space **/
  int *nGp = merged_graph->p;
  int *nGi = merged_graph->i;

  int p_counter = 1;
  int i_counter = 0;
  nGp[0] = 0;
  /** first <ngraphs>-1 graphs**/
  for(int i = 0; i < ngraphs-1; i++) {
   Gp = Gps[i];
   Gi = Gis[i];

   for(int j = 0; j < n; j++) {
    int diff = Gp[j+1] - Gp[j] + 1;
    nGp[p_counter] = nGp[p_counter-1] + diff;
    p_counter++;
    for(int p = Gp[j]; p < Gp[j+1]; p++) {
     int row = Gi[p] + (i * n);
     nGi[i_counter] = row;
     i_counter++;
    }
    nGi[i_counter] = j + (i+1) * n;
    i_counter++;
   }
  }
  /** last graph **/
  Gp = Gps[ngraphs-1];
  Gi = Gis[ngraphs-1];
  for(int j = 0; j < n; j++) {
   int diff = Gp[j+1] - Gp[j];
   nGp[p_counter] = nGp[p_counter-1] + diff;
   p_counter++;
   for(int p = Gp[j]; p < Gp[j+1]; p++) {
    int row = Gi[p] + (ngraphs-1) * n;
    nGi[i_counter] = row;
    i_counter++;
   }
  }
  return merged_graph;
 }

/*
 * Takes ngraphs and n-1 dependence graph and merge them
 * G1 -> DG1 -> G2 -> DG2 -> G3
 */
 CSC* merge_DAGs_with_partial_order(int ngraphs, int n, int **Gps, int **Gis,
   int nd, int **DGps, int **DGis){
  const int *Gp, *Gi, *DGp, *DGi;
  int nnz = 0;
  for(int i = 0; i < ngraphs; i++){
   nnz += Gps[i][n];
  }
  for(int i = 0; i < ngraphs-1; i++){
   nnz += DGps[i][n];
  }
  nnz += (ngraphs-1) * n;
  CSC *merged_graph = new CSC(n*ngraphs,n*ngraphs,nnz, true);

/** allocate new graph space **/

  int *nGp = merged_graph->p;
  int *nGi = merged_graph->i;

  int p_counter = 1;
  int i_counter = 0;
  nGp[0] = 0;

/** first <ngraphs>-1 graphs**/

  for(int i = 0; i < ngraphs-1; i++) {
   Gp = Gps[i];
   Gi = Gis[i];
   DGp = DGps[i];
   DGi = DGis[i];

   for(int j = 0; j < n; j++) {
    int diff = Gp[j+1] - Gp[j];
    diff += (DGp[j+1] - DGp[j]);
    nGp[p_counter] = nGp[p_counter-1] + diff;
    p_counter++;
    for(int p = Gp[j]; p < Gp[j+1]; p++) {
     int row = Gi[p] + (i * n);
     nGi[i_counter] = row;
     i_counter++;
    }
    for(int p = DGp[j]; p < DGp[j+1]; p++) {
     int row = DGi[p] + (i+1) * n;
     nGi[i_counter] = row;
     i_counter++;
    }
   }
  }

/** last graph **/

  Gp = Gps[ngraphs-1];
  Gi = Gis[ngraphs-1];
  for(int j = 0; j < n; j++) {
   int diff = Gp[j+1] - Gp[j];
   nGp[p_counter] = nGp[p_counter-1] + diff;
   p_counter++;
   for(int p = Gp[j]; p < Gp[j+1]; p++) {
    int row = Gi[p] + (ngraphs-1) * n;
    nGi[i_counter] = row;
    i_counter++;
   }
  }
  return merged_graph;
 }



 int post_order_spliting(int inSize, int *inTree, double *inCost,
                       int *inChildPtr, int *inChildNo,//Children list
                       int *nChild,int n, int partitionNum,
   /*Outputs*/
                       int &outSize, double* outCost,
                       int* outNode2Par, std::vector<std::vector<int>> &parList){
  auto *visited = new bool[inSize]();
  auto *mark = new bool[inSize]();
  int Threshold, k=0;
  Threshold=n/partitionNum;//Almost the same numjber of columns in each partition
  //Threshold=1;
  std::vector<int> stack;
  //Initial partitioning, each partition one node.
  //parList.resize(inSize);
  std::vector<int> extraDim;
  parList.push_back(extraDim);
#if 0
  for (int l = 0; l < inSize; ++l) {
        std::cout<<l<<": "<<inTree[l]<<",";
    }
    std::cout<<"\n";
#endif
  int curPart=0;
  outSize=0;
  for (int curNode = 0; curNode < inSize; ++curNode) {//k is in partitioned node
   if(inTree[curNode]==-2)
    continue;
   if(!visited[curNode]){
    stack.push_back(curNode);
    while (stack.size()!=0 ) {
     k=stack[0];
     if (nChild[k] == 0) {
      //Add it to current part
      stack.erase(stack.begin());
      parList[curPart].push_back(k);
      outCost[curPart]+=inCost[k];
      visited[k]= true;//mark as visited
      //nChild[inTree[k]] =k>=0? nChild[inTree[k]]-1:nChild[inTree[k]];
      if(inTree[k]>=0) {//The node is a single node
       nChild[inTree[k]]--;
       if(!visited[inTree[k]]){
        if(!mark[inTree[k]])//if it is not in the stack
         stack.insert(stack.begin(),inTree[k]);//Add its parent
        mark[inTree[k]]=true;//There is at least a children of this node in this partition
       }
      }
      else//There is no other node to check, check next one
       break;
     } else{
      //adding its children for further evaluation
//                    std::cout<<"The children of: "<<k<<"\n";
      for (int i = inChildPtr[k]; i < inChildPtr[k+1]; ++i) {
       int tmpChild=inChildNo[i];
       if(!visited[tmpChild])
        if(nChild[tmpChild]==0){//First add leaves
         //    stack.insert(stack.begin(),tmpChild);
         //    tmpFront++;
         if(!visited[k]){
          int tmpk=k;
          if(stack.size()>0){
           while (tmpk!=stack[stack.size()-1]){//Mark all parents upwards, they should be in the partiotion
            mark[tmpk]=true;//There is at least a children of this node in this partition
            tmpk=inTree[tmpk];
           }
          }
          mark[tmpk]=true;
         }
         parList[curPart].push_back(tmpChild);
         outCost[curPart]+=inCost[tmpChild];
         visited[tmpChild]= true;//mark as visited
         nChild[inTree[tmpChild]]--;//we know inTree[tmpChild] is k
        }else{
         stack.insert(stack.begin(),tmpChild);
        }
      }
     }//End else
/*    if(outCost[curPart]>Threshold){
     for (int s = 0; s < stack.size(); ++s) {
      int leftover=stack[s];
      if(mark[leftover]){
       parList[curPart].push_back(leftover);
       outCost[curPart]+=inCost[leftover];
       visited[leftover]= true;//mark as visited
       mark[leftover]=false;//for further iterations
       if(inTree[leftover]>=0)
        nChild[inTree[leftover]]--;
      }
     }
     stack.erase(stack.begin(),stack.end());
     break;
    }*/
    }
    //It is either reaches the defined size
    // or there is no other node to add (stack is empty)
    curPart++;
    std::vector<int> extraDim;
    parList.push_back(extraDim);
   }
  }
  //remove the last empty element
  parList.erase(parList.begin()+parList.size()-1);
  outSize=curPart;
  for (int i = 0; i < curPart; ++i) {
   for (int j = 0; j < parList[i].size(); ++j) {
    assert(parList[i][j] < inSize && parList[i][j] >=0);
    outNode2Par[parList[i][j]]=i;
   }
  }

#if 0
  for (int i = 0; i < curPart; ++i) {
        std::cout<<"Partition #"<<i<<" :";
        for (int j = 0; j < parList[i].size(); ++j) {
            std::cout<<parList[i][j]<<",";
        }
        std::cout<<"\n";
    }
#endif
#if 0
  auto *check = new bool[inSize]();
    for (int i = 0; i < curPart; ++i) {
        for (int j = 0; j < parList[i].size(); ++j) {
                check[parList[i][j]]=true;
        }
    }
    for (int i = 0; i < inSize; ++i) {
        assert(check[i]);
    }
    delete []check;
#endif
  delete []mark;
  delete []visited;
  return outSize;
 }

 /* x = x + beta * A(:,j), where x is a dense vector and A(:,j) is sparse */
 int scatter (const CSC *A, int j, double beta, int *w,
              double *x, int mark, CSC *C, int nz){
  int i, p, *Ap, *Ai, *Ci ;
  double *Ax ;
  Ap = A->p ; Ai = A->i ; Ax = A->x ; Ci = C->i ;
  for (p = Ap [j] ; p < Ap [j+1] ; p++){
   i = Ai [p] ;                            /* A(i,j) is nonzero */
   if (w [i] < mark){
    w [i] = mark ;                      /* i is new entry in column j */
    Ci [nz++] = i ;                     /* add i to pattern of C(:,j) */
    if (x) x [i] = beta * Ax [p] ;      /* x(i) = beta*A(i,j) */
   }
   else if (x) x [i] += beta * Ax [p] ;    /* i exists in C(:,j) already */
  }
  return (nz) ;
 }

 // TODO: replace the following with a more efficient way of making matrix symmetric
 // TODO:  it is copied here to remove sparseblas dependency
 CSC* add_tmp(CSC *A, CSC *B, double alpha, double beta, bool sort=true){
  if(A->stype != B->stype)
   return NULLPNTR;
  int p, j, nz = 0, anz, *Cp, *Ci, *Bp, m, n, bnz, *w ;
  bool values = A->is_pattern || B->is_pattern; //true if either is pattern
  double *x, *Bx, *Cx ;
  m = A->m ; anz = A->p [A->n] ;
  n = B->n ; Bp = B->p ; Bx = B->x ; bnz = Bp [n] ;
  w = new int[m]();
  x = !values ? new double[m]() : NULLPNTR;
  CSC *C = new CSC(m,n,anz+bnz,values);
  C->stype = A->stype;
  C->is_pattern = A->is_pattern;
  Cp = C->p ; Ci = C->i ; Cx = C->x ;
  for (j = 0 ; j < n ; j++){
   Cp [j] = nz ;                   /* column j of C starts here */
   nz = scatter (A, j, alpha, w, x, j+1, C, nz) ;   /* alpha*A(:,j)*/
   nz = scatter (B, j, beta, w, x, j+1, C, nz) ;    /* beta*B(:,j) */
   if(sort) std::sort(&C->i[Cp[j]], &C->i[nz]);
   if (!values) for (p = Cp [j] ; p < nz ; p++) Cx [p] = x [Ci [p]] ;
  }
  Cp [n] = nz ;
  delete []x;
  delete []w;
  return C;
 }

 CSC * make_symmetric(CSC *A, bool lower){
  CSC *ApAt, *At;
  if((A->stype>0 && !lower) || (A->stype<0 && lower))
   return A;
  if((A->stype<0 && !lower) || (A->stype>0 && lower) )
   return transpose_symmetric(A, NULLPNTR);

  At = transpose_general(A);
  ApAt = add_tmp(A,At,0.5,0.5);//(A+At)/2
  CSC* As = make_half(ApAt->n, ApAt->p, ApAt->i, ApAt->x, lower);
  As->stype = lower ? -1 : 1;
  delete At;
  delete ApAt;
  return As;
 }


 void compute_nnz_per_col(CSC *A, double *nnz_cnt){
  for (int i = 0; i < A->n; ++i) {
   nnz_cnt[i] = A->p[i+1] - A->p[i];
  }
 }


 void reorder_array(int n, int *arr, int *perm, int *ws){
  for (int i = 0; i < n; ++i) {
   ws[i] = arr[perm[i]];
  }
  for (int j = 0; j < n; ++j) {
   arr[j] = ws[j];
  }
 }

 void inv_perm(int n, int *perm, int *iperm){
  for (int i = 0; i < n; ++i) {
   iperm[perm[i]] = i;
  }
 }


 CSC *permute_general (const CSC *A, const int *pinv,
   const int *q)
 {
  int nz = 0, n, *Ap, *Ai, *Cp, *Ci ;
  double *Cx, *Ax ;
  int m = A->m; n = A->n ; Ap = A->p ; Ai = A->i ; Ax = A->x ;
  CSC *C = new CSC(m, n, A->nnz);
  Cp = C->p ; Ci = C->i ; Cx = C->x ;
  for (int k = 0 ; k < n ; k++)
  {
   Cp [k] = nz ; //* column k of C is column q[k] of A */
   int j = q ? (q [k]) : k ;
   for (int t = Ap [j] ; t < Ap [j+1] ; t++)
   {
    if (Cx) Cx [nz] = Ax [t] ;  /* row i of A is row pinv[i] of C */
    Ci [nz++] = pinv ? (pinv [Ai [t]]) : Ai [t] ;
   }
  }
  Cp [n] = nz ; //* finalize the last column of C */
  return C;
 }


 CSC *coarsen_k_times(int n, int nnz, int *Ap, int *Ai, int stype, int k){
  int new_n = (n%k == 0) ? n/k : (n/k) + 1;
  CSC *c_mat = new CSC(new_n, new_n, nnz, false);
  c_mat->stype = stype;
  int *Bp = c_mat->p; int *Bi = c_mat->i; int cur_col=1, cur_nnz=0;
  Bp[0] = 0;
  std::vector<int> faled_vals;
  auto is_exist = new bool[n]();
  for (int j = 0; j < n; j+=k) {
   int bnd = std::min<int>(j + k,n);
   for (int l = j; l < bnd; ++l) {
    for (int m = Ap[l]; m < Ap[l + 1]; ++m) {
     auto tmp = Ai[m];
     int new_idx = tmp / k;
     if(!is_exist[new_idx]){
      is_exist[new_idx] = true;
      assert(new_idx < new_n);
      Bi[cur_nnz] = new_idx;
      cur_nnz++;
      faled_vals.push_back(new_idx);
     }
    }
   }
   Bp[cur_col] = cur_nnz;
   cur_col++;
   for (auto i : faled_vals ) {
    is_exist[i] = false;
   }
   faled_vals.clear();
   //std::fill_n(is_exist,n,false);
  }
  assert(new_n == cur_col-1);
  c_mat->nnz = cur_nnz;
  delete is_exist;
  return c_mat;
 }

}

