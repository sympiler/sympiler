//
// Created by george on 2019-10-14.
//

#include <cstring>
#include <queue>
#include <cassert>
#include "sparse_utilities.h"
#include "includes/def.h"
#include "includes/BCSCMatrix.h"
#include "includes/sparse_inspector.h"

namespace sym_lib {

 int build_levelSet_CSC(size_t n, int *Lp, int *Li,
                        int *&levelPtr, int *&levelSet) {
  int begin = 0, end = n - 1;
  int cur_level = 0, cur_levelCol = 0;
  levelPtr = new int[n + 1]();
  levelSet = new int[n]();
  int *inDegree = new int[n]();
  bool *visited = new bool[n]();
  for (int i = 0; i < Lp[n]; ++i) {//O(nnz)
   inDegree[Li[i]]++;
  }
  //print_vec("dd\n",0,n,inDegree);
  while (begin <= end) {
   for (int i = begin; i <= end; ++i) {//For level cur_level
    if (inDegree[i] == 1 && !visited[i]) {//if no incoming edge
     visited[i] = true;
     levelSet[cur_levelCol] = i; //add it to current level
     cur_levelCol++;//Adding to level-set
    }
   }
   cur_level++;//all nodes_ with zero indegree are processed.
   //assert(cur_level < n);
   if (cur_level >= n)
    return -1; // The input graph has a cycle
   levelPtr[cur_level] = cur_levelCol;
   while (inDegree[begin] == 1) {
    begin++;
    if (begin >= n)
     break;
   }
   while (inDegree[end] == 1 && begin <= end)
    end--;
   //Updating degrees after removing the nodes_
   for (int l = levelPtr[cur_level - 1]; l < levelPtr[cur_level]; ++l) {
    int cc = levelSet[l];
    for (int j = Lp[cc]; j < Lp[cc + 1]; ++j) {
     if (Li[j] != cc) //skip diagonals
      inDegree[Li[j]]--;//removing corresponding edges
    }
   }
   //print_vec("dd\n",0,n,inDegree);
  }
  delete[]inDegree;
  delete[]visited;
  return cur_level;//return number of levels
 }

 int build_levelSet_BCSC(int n, int nnz, size_t *Lp, size_t *Li_ptr, int *Li,
                         size_t blockNo, const int *sup2col, const int *col2sup,
                         int *&levelPtr, size_t *&levelSet) {
  int begin = 0, end = blockNo - 1;
  int cur_level = 0, cur_levelCol = 0, curCol, nxtCol, supWdt;

  levelPtr = new int[blockNo + 1]();
  levelSet = new size_t[blockNo]();
  int *inDegree = new int[blockNo]();
  bool *visited = new bool[blockNo]();
  int *node2Level = new int[blockNo]();
  for (size_t i = 0; i < blockNo; ++i) {
   inDegree[i] = 0;
   visited[0] = false;
  }
  for (size_t s = 1;
       s <= blockNo; ++s) { //In degree computation in assembly tree
   curCol = s != 0 ? sup2col[s - 1] : 0;
   nxtCol = sup2col[s];
   supWdt = nxtCol - curCol;
   for (size_t r = Li_ptr[curCol] + supWdt - 1; r < Li_ptr[nxtCol]; ++r) {
    inDegree[col2sup[Li[r]]]++;
   }
  }
  while (begin <= end) {
   for (int i = begin; i <= end; ++i) {//For level cur_level
    if (inDegree[i] == 1 && !visited[i]) {//if no incoming edge
     visited[i] = true;
     assert(i >= 0);
     levelSet[cur_levelCol] = i; //add it to current level
     node2Level[i] = cur_level;
     cur_levelCol++;//Adding to level-set
    }
   }
   cur_level++;//all nodes_ with 1 indegree are processed.
   assert((size_t) cur_levelCol >= 0);
   assert((size_t) cur_level <= blockNo);
   levelPtr[cur_level] = cur_levelCol;
   while (inDegree[begin] == 1)
    begin++;
   while (inDegree[end] == 1 && begin <= end)
    end--;
   //Updating degrees after removing the nodes_
   for (int l = levelPtr[cur_level - 1]; l < levelPtr[cur_level]; ++l) {
    size_t cc = levelSet[l] + 1;
    curCol = cc != 0 ? sup2col[cc - 1] : 0;
    nxtCol = sup2col[cc];
    supWdt = nxtCol - curCol;
    for (size_t j = Li_ptr[curCol] + supWdt; j < Li_ptr[nxtCol]; ++j) {
     inDegree[col2sup[Li[j]]]--;//removing corresponding edges
    }
   }
  }
  delete[]inDegree;
  delete[]visited;
  return cur_level;//return number of levels
 }


 /// other version of levelset CSC
 int level_set_multi_graphs(int n, const int *Lp, const int *Li, int *&levelptr,
                            int *&levels, int n_kern) {
  int *levelPtr = new int[n_kern * n + 1]();
  int *degrees = new int[n_kern * n]();
  int *visited = new int[n_kern * n]();
  bool *inqueue = new bool[n_kern * n]();
  levels = new int[n_kern * n]();
  memset(degrees + n, 1, sizeof(int) * n * (n_kern - 1));

  /* find source nodes_ and append in initial queue */
  std::queue<int> queue;
  for (int i = 0; i < n; i++) {
   for (int p = Lp[i] + 1; p < Lp[i + 1]; p++) {
    for (int q = 0; q < n_kern; q++) {
     degrees[Li[p] + n * q]++;
    }
   }
  }
  for (int i = 0; i < n; i++) {
   if (degrees[i] == 0) {
    queue.push(i);
    inqueue[i] = true;
   }
  }

  /* dequeue each level while enqueue each level */
  int ptr_num = 0, lev_num = 0;
  int cur_level = 0;
  int cur_n = queue.size();
  int next_n = 0;

  levelPtr[0] = 0;
  while (!queue.empty()) {
   levelPtr[cur_level + 1] = levelPtr[cur_level] + cur_n;

   for (int i = 0; i < cur_n; i++) {
    int index = queue.front();
    queue.pop();
    visited[index]++;
    inqueue[index] = false;

    levels[lev_num] = index;
    lev_num++;

    // iterate through all out-going edges
    int kern_n = index / n;
    index -= kern_n * n;

    for (int p = Lp[index] + 1; p < Lp[index + 1]; p++) {
     // if not visited, add to current level
     int offset = Li[p] + n * kern_n;
     if (visited[offset] == degrees[offset] - 1 && !inqueue[offset]) {
      queue.push(offset);
      inqueue[offset] = true;
      next_n++;
     } else {
      visited[offset]++;
     }
    }

    // set up for next level
    if (kern_n < n_kern - 1) {
     queue.push(index + n * (kern_n + 1));
     inqueue[index + n * (kern_n + 1)] = true;
     next_n++;
    }
   }

   cur_n = next_n;
   ptr_num += cur_n;
   next_n = 0;

   cur_level++;
  }

  // copy information
  levelptr = new int[cur_level + 1];
  std::memcpy(levelptr, levelPtr, sizeof(int) * (cur_level + 1));

  delete[]levelPtr;
  delete[]visited;
  delete[]degrees;
  delete[]inqueue;
  return cur_level;
 }


 int
 level_set_bn(int n, const int *Lp, const int *Li, int *&levelptr, int *&levels,
              int *supernodes, int num_nodes, int n_kern, bool isLower) {
  int *levelPtr = new int[n_kern * num_nodes + 1]();
  levels = new int[n_kern * num_nodes]();
  int *degrees = new int[n_kern * n]();
  int *visited = new int[n_kern * n]();
  int *supers = new int[n]();
  bool *inqueue = new bool[n_kern * n]();
  memset(degrees + n, 1, sizeof(int) * n * (n_kern - 1));

  int cur_level = 0;
  if (isLower) {
   // initialize the super nodes_ of all nodes_
   int counter = 0;
   for (int i = 0; i < num_nodes; i++) {
    while (counter >= supernodes[i] && counter < supernodes[i + 1]) {
     supers[counter] = supernodes[i];
     counter++;
    }
   }

   // initialize degrees of all supernodes
   for (int i = 0; i < num_nodes; i++) {
    int index = supernodes[i];
    for (int p = Lp[index] + 1; p < Lp[index + 1]; p++) {
     int super = supers[Li[p]];
     if (index != super) {
      for (int q = 0; q < n_kern; q++) {
       degrees[super + q * n]++;
      }
     }
    }
   }

   /* find source nodes_ and append in initial queue */
   std::queue<int> queue;
   for (int i = 0; i < num_nodes; i++) {
    if (degrees[supernodes[i]] == 0) {
     queue.push(supernodes[i]);
     inqueue[supernodes[i]] = true;
    }
   }

   /* dequeue each level while enqueue each level */
   int ptr_num = 0, lev_num = 0;
   int cur_n = queue.size();
   int next_n = 0;

   levelPtr[0] = 0;
   while (!queue.empty()) {
    levelPtr[cur_level + 1] = levelPtr[cur_level] + cur_n;

    for (int i = 0; i < cur_n; i++) {
     int index = queue.front();
     queue.pop();
     visited[index]++;
     inqueue[index] = false;

     // make index be in [0, n)
     levels[lev_num] = index;
     int kern_num = index / n;
     index -= kern_num * n;
     lev_num++;

     // iterate through all out-going edges
     for (int p = Lp[index] + 1; p < Lp[index + 1]; p++) {
      int super = supers[Li[p]] + n * kern_num;
      visited[super]++;

      // if not visited, add to current level
      if (visited[super] == degrees[super] && !inqueue[super]) {
       queue.push(super);
       inqueue[super] = true;
       next_n++;
      }
     }

     // setup for next level
     if (kern_num < n_kern - 1) {
      queue.push(index + n * (kern_num + 1));
      inqueue[index + n * (kern_num + 1)] = true;
      next_n++;
     }
    }

    cur_n = next_n;
    ptr_num += cur_n;
    next_n = 0;

    cur_level++;
   }
  } else {
   // initialize the super nodes_ of all nodes_
   int counter = n - 1;
   for (int i = 0; i < num_nodes; i++) {
    while (counter <= supernodes[i] && counter > supernodes[i + 1]) {
     supers[counter] = supernodes[i];
     counter--;
    }
   }

   // initialize degrees of all supernodes
   for (int i = 0; i < num_nodes; i++) {
    int index = supernodes[i];
    for (int p = Lp[index + 1] - 1; p >= Lp[index]; p--) {
     int super = supers[Li[p]];
     if (index != super) {
      for (int q = 0; q < n_kern; q++) {
       degrees[super + q * n]++;
      }
     }
    }
   }

   /* find source nodes_ and append in initial queue */
   std::queue<int> queue;
   for (int i = 0; i < num_nodes; i++) {
    if (degrees[supernodes[i]] == 0) {
     queue.push(supernodes[i]);
     inqueue[supernodes[i]] = true;
    }
   }

   /* dequeue each level while enqueue each level */
   int ptr_num = 0, lev_num = 0;
   int cur_n = queue.size();
   int next_n = 0;

   levelPtr[0] = 0;
   while (!queue.empty()) {
    levelPtr[cur_level + 1] = levelPtr[cur_level] + cur_n;

    for (int i = 0; i < cur_n; i++) {
     int index = queue.front();
     queue.pop();
     visited[index]++;
     inqueue[index] = false;

     // make index be in [0, n)
     levels[lev_num] = index;
     int kern_num = index / n;
     index -= kern_num * n;
     lev_num++;

     // iterate through all out-going edges
     for (int p = Lp[index + 1] - 1; p >= Lp[index]; p--) {
      int super = supers[Li[p]] + n * kern_num;
      visited[super]++;

      // if not visited, add to current level
      if (visited[super] == degrees[super] && !inqueue[super]) {
       queue.push(super);
       inqueue[super] = true;
       next_n++;
      }
     }

     // setup for next level
     if (kern_num < n_kern - 1) {
      queue.push(index + n * (kern_num + 1));
      inqueue[index + n * (kern_num + 1)] = true;
      next_n++;
     }
    }

    cur_n = next_n;
    ptr_num += cur_n;
    next_n = 0;

    cur_level++;
   }
  }

  // copy information
  levelptr = new int[cur_level + 1];
  std::memcpy(levelptr, levelPtr, sizeof(int) * (cur_level + 1));

  delete[]levelPtr;
  delete[]visited;
  delete[]degrees;
  delete[]inqueue;
  delete []supers;
  return cur_level;
 }


 int build_level_set_tree_efficient(size_t n, const int *tree,
   const int *nChild1,
   int *levelPtr, int *levelSet, int *node2level) {
  int maxLen = 0, level = 0, num_levels=0;
  auto *level_cnt = new int[n]();
  for (int i = 0; i < n; ++i) {
   if (nChild1[i] == 0) {
    levelSet[levelPtr[1]] = i;
    levelPtr[1] ++;
    int node = i; level = 0;
    while (tree[node] >= 0) {
     node = tree[node];
     level++;
     if(level < node2level[node])
      break; // a longer path visted this node before
     else
      node2level[node] = level;// touched this node with level
    }
    if (level > maxLen)
     maxLen = level;
   }
  }
  num_levels = maxLen + 1;
  for (int j = 0; j < n; ++j) {
   level_cnt[node2level[j]]++;
  }
  assert(level_cnt[0] == levelPtr[1]);
  for (int k = 2; k < num_levels+1; ++k) {
   levelPtr[k] = levelPtr[k-1] + level_cnt[k-1];
   level_cnt[k-1] = 0;
  }
  level_cnt[1]=0;
  for (int l = 0; l < n; ++l) {
   auto lev = node2level[l];
   if(!lev) continue;
   levelSet[levelPtr[lev]+level_cnt[lev]] = l;
   level_cnt[lev]++;
  }
  delete []level_cnt;
  return num_levels;
 }


 int build_level_set_tree(size_t n, const int *inTree, int *levelPtr,
                          int *levelSet) {
//Naive code for generating level levelSet from ETree, it is not part of
// the inspector. In a real code, we return a node or supernode at
// a time so, it will in O(n)
  int begin = 0, end = (int) n - 1;
  auto *nChild = new int[n]();
  auto *visited = new bool[n]();
  //Counting the number of children
  for (int k = 0; k < n; ++k) {
   if (inTree[k] >= 0)
    nChild[inTree[k]]++;
  }
  int curLevel = 0, curLevelCnt = 0;
  levelPtr[0] = 0;

  while (begin <= end) {
   for (int i = begin; i <= end; ++i) {//For level curLevel
    if (nChild[i] == 0 && !visited[i]) {//if no incoming edge
     visited[i] = true;
     levelSet[curLevelCnt] = i; //add it to current level
     curLevelCnt++;//Adding to level-set
    }
   }
   curLevel++;//all nodes with zero indegree are processed.
   levelPtr[curLevel] = curLevelCnt;
   if (curLevelCnt == n)
    break;
   while (nChild[begin] == 0)
    begin++;
   while (nChild[end] == 0 && begin <= end)
    end--;
   //Updating degrees after removing the nodes
   for (int l = levelPtr[curLevel - 1]; l < levelPtr[curLevel]; ++l) {
    int cc = levelSet[l];
    if (inTree[cc] >= 0)
     nChild[inTree[cc]]--;
   }
  }
#if DEBUG >= 2
  std::cout<<"FinalSet\n";
 for (int l = 0; l < curLevel; ++l) {
  for (int lp = levelPtr[l]; lp < levelPtr[l+1]; ++lp) {
   std::cout<<levelSet[lp]<<",";
  }
  std::cout<<"\n\n";
 }
 std::cout<<"\n";
#endif
  delete[]nChild;
  delete[]visited;
  return curLevel;

 }

 int find(const int *supernodes, int start, int end, int target) {
  if (start == end && supernodes[start] == target)
   return start;
  if (end - start == 1) {
   if (supernodes[start] == target) return start;
   else if (supernodes[end] == target) return end;
  }

  int m = (end + start) / 2;
  if (supernodes[m] == target)
   return m;
  else if (target < supernodes[m]) {
   return find(supernodes, start, m, target);
  } else {
   return find(supernodes, m, end, target);
  }
 }

 void sup2node_gen(int n, int size, const int *levels, int num_nodes,
                   int *supernodes, int *&sup2node) {
  sup2node = new int[size]();

  for (int i = 0; i < size; i++) {
   int super = levels[i];
   int real = super - (super / n) * n;
   sup2node[i] = find(supernodes, 0, num_nodes, real);
  }
 }


 // assumes that A, C have the same supernodes
 // FIXME: revert back?
 int bcsc_csc_bcsc_levelset(CSC *A, CSC *B, CSC *C, int *&levelptr, int *&levels,
                         int nodes, int *supernodes) {
  // check that the supernodes are the same (using size right now)
  // check dimensions of inputs
  assert(A->n == B->n);
  assert(B->n == C->n);
  assert(B->n == B->m);

  int n = A->n;

  auto degrees = new int[2 * n]();
  auto visited = new int[2 * n]();
  auto inqueue = new bool[2 * n]();
  auto levelptr_ = new int[2 * n + 1]();
  levels = new int[2 * nodes]();


  // initialize the super nodes_ of all nodes_
  int counter = 0;
  auto supers = new int[n](); // A, B, C have same supernodal structure
  for (int i = 0; i < nodes; i++) {
   while (counter >= supernodes[i] && counter < supernodes[i + 1]) {
    supers[counter] = supernodes[i];
    counter++;
   }
  }

  // this is the same if B is symmetric, only works for SpMV
//  CSR* row_dep = csc_to_csr(B);

  // initialize in-degree structure
  // First BCSC and Second BCSC internal dependency
  for (int i = 0; i < nodes; i++) {
   int index = supernodes[i];
   for (int p = A->p[index] + 1; p < A->p[index + 1]; p++) {
    int super = supers[A->i[p]];
    if (index != super) {
     degrees[super]++; // for A
    }
   }
   for (int p = C->p[index] + 1; p < C->p[index + 1]; p++) {
    int super = supers[C->i[p]];
    if (index != super)
     degrees[super + n]++;
   }
  }

  for (int i = 0; i < nodes; i++) {
   int index = supernodes[i];
   int width = supernodes[i + 1] - supernodes[i];

   // each SpMV-CSC only depends on previous BCSC node => <width> number of nodes
   for (int p = B->p[index]; p < B->p[index + width]; p++) {
    int super = supers[B->i[p]];
    degrees[super + n]++;
   }
  }

  /* find source nodes_ and append in initial queue */
  std::queue<int> queue;
  for (int i = 0; i < nodes; i++) {
   if (degrees[supernodes[i]] == 0) {
    queue.push(supernodes[i]);
    inqueue[supernodes[i]] = true;
   }
  }

  int cur_level = 0;
  int ptr_num = 0, lev_num = 0;
  int cur_n = queue.size();
  int next_n = 0;

  levelptr_[0] = 0;
  while (!queue.empty()) {
   levelptr_[cur_level + 1] = levelptr_[cur_level] + cur_n;

   for (int i = 0; i < cur_n; i++) {
    int index = queue.front();
    queue.pop();
    visited[index]++;
    inqueue[index] = false;

    // make index be in [0, n)
    levels[lev_num] = index;
    int kern_num = index / n;
    index -= kern_num * n;
    lev_num++;

    if (kern_num == 0) { // BFS for A
     // subsequent SpTRSV nodes
     for (int p = A->p[index] + 1; p < A->p[index + 1]; p++) {
      int super = supers[A->i[p]];
      visited[super]++;

      if (visited[super] == degrees[super] && !inqueue[super]) {
       queue.push(super);
       inqueue[super] = true;
       next_n++;
      }
     }

     // next SpTRSV node (fuse SpMV into current one)
     int start = find(supernodes, 0, nodes, index);
     int width = supernodes[start + 1] - supernodes[start];

     for (int p = B->p[index]; p < B->p[index + width]; p++) {
      int super = supers[B->i[p]] + n;
      visited[super]++;

      if (visited[super] == degrees[super] && !inqueue[super]) {
       queue.push(super);
       inqueue[super] = true;
       next_n++;
      }
     }
    } else {
     for (int p = C->p[index] + 1; p < C->p[index + 1]; p++) {
      int super = supers[C->i[p]] + n;
      visited[super]++;

      if (visited[super] == degrees[super] && !inqueue[super]) {
       queue.push(super);
       inqueue[super] = true;
       next_n++;
      }
     }
    }
   }
   cur_n = next_n;
   ptr_num += cur_n;
   next_n = 0;

   cur_level++;
  }

  delete[]visited;
  delete[]degrees;
  delete[]inqueue;

  levelptr = new int[cur_level + 1];
  std::memcpy(levelptr, levelptr_, sizeof(int) * (cur_level + 1));

  delete[]levelptr_;
  return cur_level;
 }


 CSC *merge_graphs(BCSC *A, CSC *B, BCSC *C) {
  int nodes = A->nodes;
  int *supernodes = A->supernodes;

  auto nB_M = new BCSCMatrix(B, nodes, supernodes);
  auto nB = nB_M->getBCSC();

  int counter = 0;
  auto supers = new int[A->n](); // A, B, C have same supernodal structure
  for (int i = 0; i < nodes; i++) {
   while (counter >= supernodes[i] && counter < supernodes[i + 1]) {
    supers[counter] = i;
    counter++;
   }
  }

  int nnz = 0;

  // Count # nnz in A
  bool *inplace = new bool[nodes]();
  int count = 0;
  int *temp = new int[nodes]();

  // Count # nnz in A
  for (int i = 0; i < nodes; i++) {
   for (int r = 0; r < A->nrows[i]; r++) {
    int row = supers[A->i[A->p[i] + r]];
    if (!inplace[row]) {
     inplace[row] = true;
     temp[count] = row;
     count++;
     nnz++;
    }
   }
   for (int j = 0; j < count; j++)
    inplace[temp[j]] = false;
   count = 0;

   nnz++; // from A to B
  }
  // Count # nnz in B
  for (int i = 0; i < nodes; i++) {
   for (int r = 0; r < nB->nrows[i]; r++) {
    int row = supers[nB->i[nB->p[i]]];
    if (!inplace[row]) {
     inplace[row] = true;
     temp[count] = row;
     count++;
     nnz++;
    }
   }
   for (int j = 0; j < count; j++)
    inplace[temp[j]] = false;
   count = 0;
   nnz++;
  }
  // Count # nnz in C
  for (int i = 0; i < nodes; i++) {
   for (int r = 0; r < C->nrows[i]; r++) {
    int row = supers[C->i[C->p[i] + r]];
    if (!inplace[row]) {
     inplace[row] = true;
     temp[count] = row;
     count++;
     nnz++;
    }
   }
   for (int j = 0; j < count; j++)
    inplace[temp[j]] = false;
   count = 0;
  }

  CSC *M = new CSC(3 * nodes, 3 * nodes, nnz);
  auto Mp = M->p;
  auto Mi = M->i;

  int p = 1, i = 0;
  Mp[0] = 0;
  for (int j = 0; j < nodes; j++) {
   Mp[p] = Mp[p - 1];
   for (int r = 0; r < A->nrows[j]; r++) {
    int row = supers[A->i[A->p[j] + r]];
    if (!inplace[row]) {
     inplace[row] = true;
     temp[count] = row;
     count++;
     Mp[p]++;
     Mi[i] = row;
     i++;
    }
   }
   for (int k = 0; k < count; k++)
    inplace[temp[k]] = false;
   count = 0;

   Mi[i] = j + nodes;
   Mp[p]++;
   i++;
   p++;
  }

  for (int j = 0; j < nodes; j++) {
   Mp[p] = Mp[p - 1];
   Mp[p]++;
   Mi[i] = j + nodes;
   i++;
   for (int r = 0; r < nB->nrows[j]; r++) {
    int row = supers[nB->i[nB->p[j]]];
    if (!inplace[row]) {
     inplace[row] = true;
     temp[count] = row;
     count++;

     Mp[p]++;
     Mi[i] = row + 2 * nodes;
     i++;
    }
   }
   for (int k = 0; k < count; k++)
    inplace[temp[k]] = false;
   count = 0;
   p++;
  }

  for (int j = 0; j < nodes; j++) {
   Mp[p] = Mp[p - 1];
   for (int r = 0; r < C->nrows[j]; r++) {
    int row = supers[C->i[C->p[j] + r]];
    if (!inplace[row]) {
     inplace[row] = true;
     temp[count] = row;
     count++;

     Mp[p]++;
     Mi[i] = row + 2 * nodes;
     i++;
    }
   }
   for (int k = 0; k < count; k++)
    inplace[temp[k]] = false;
   count = 0;
   p++;
  }

  delete(nB_M);
  delete[]supers;
  delete[]inplace;
  delete[]temp;

  return M;
 }


 void
 compute_depth(CSC *A, int *degrees) {
  int *level_ptr, *level_set, level_num;
  CSC *AT = transpose_general(A); //level-set accepts lower triangular only
  //CSR *tmp =  csc_to_csr(A);
  assert(A->p[A->n] == AT->p[AT->n]);
  //level_num = level_set_multi_graphs(A->n,AT->p,AT->i,level_ptr,level_set,1);
  //print_csc("AT\n",AT->n,AT->p,AT->i,AT->x);
  level_num = build_levelSet_CSC(A->n, AT->p, AT->i, level_ptr, level_set);
  for (int i = 0; i < level_num; ++i) {
   for (int j = level_ptr[i]; j < level_ptr[i + 1]; ++j) {
    int cn = level_set[j];
    degrees[cn] = i;
   }
  }
  delete AT;
  delete[]level_ptr;
  delete[]level_set;
 }


 int dfs_tree( int p, int k, int *Post, int *Head, int *Next,int *Pstack){
  int j, phead;
  /* put the root node on the stack */
  Pstack[0] = p;
  phead = 0;

  /* while the stack is not empty, do: */
  while (phead >= 0) {
   /* grab the node p from top of the stack and get its youngest child j */
   p = Pstack[phead];
   j = Head[p];
   if (j == EMPTY) {
    /* all children of p ordered.  remove p from stack and order it */
    phead--;
    Post[k++] = p; /* order node p as the kth node */
   } else {
    /* leave p on the stack.  Start a DFS at child node j by putting
     * j on the stack and removing j from the list of children of p. */
    Head[p] = Next[j];
    Pstack[++phead] = j;
   }
  }
  return (k); /* the next node will be numbered k */
 }

/* find nonzero pattern of Cholesky L(k,1:k-1) using etree and triu(A(:,k)) */
 int ereach(int n, int *Ap, int *Ai, int k, const int *parent,
            int *s, int *w) {
  int i, p, len, top;
  if (!Ap || !Ai || !parent || !s || !w) return (-1);   /* check inputs */
  top = n;
  CS_MARK (w, k);                /* mark node k as visited */
  for (p = Ap[k]; p < Ap[k + 1]; p++) {
   i = Ai[p];                /* A(i,k) is nonzero */
   if (i > k) continue;       /* only use upper triangular part of A */
   for (len = 0; !CS_MARKED (w, i); i = parent[i]) /* traverse up etree*/
   {
    s[len++] = i;         /* L(k,i) is nonzero */
    CS_MARK (w, i);        /* mark i as visited */
   }
   while (len > 0) s[--top] = s[--len]; /* push path onto stack */
  }
  for (p = top; p < n; p++) CS_MARK (w, s[p]);    /* unmark all nodes */
  CS_MARK (w, k);                /* unmark node k */
  return (top);                  /* s [top..n-1] contains pattern of L(k,:)*/
 }

}