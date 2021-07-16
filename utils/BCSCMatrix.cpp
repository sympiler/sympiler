//
// Created by george on 2020-02-15.
//

#include "BCSCMatrix.h"
namespace sym_lib {
 int BCSCMatrix::supernodes(CSC *A, int limit, bool isLower) {
  int n = A->n;

  int num_nodes = 0;
  int *temp_nodes = new int[n]();

  if (limit == 0)
   limit = n;

  int max = 0;
  int pre, cur;

  if (isLower) {
   for (int i = 0; i < n; i++) {
    // start of new block
    temp_nodes[num_nodes] = i;
    num_nodes++;

    pre = i;
    cur = i + 1;

    int sim = 0;
    while (cur < n && sim < limit) {
     int diff = cur - pre;

     // compare number of off diagonals
     int off_pre = A->p[pre + 1] - A->p[pre] - diff;
     int off_cur = A->p[cur + 1] - A->p[cur];

     if (off_pre != off_cur)
      break;

     // now examine row pattern
     bool col_sim = true;
     for (int j = 0; j < off_cur; j++) {
      int pre_index = A->p[pre] + diff + j;
      int cur_index = A->p[cur] + j;

      if (A->i[pre_index] != A->i[cur_index]) {
       col_sim = false;
       break;
      }
     }

     if (col_sim) {
      sim++;
      cur++;
     } else
      break;
    }
    i = pre + sim;

    if (sim > max)
     max = sim;
   }
  } else {
   for (int i = n - 1; i >= 0; i--) {
    // start of new block
    temp_nodes[num_nodes] = i;
    num_nodes++;

    pre = i;
    cur = i - 1;

    int sim = 0;
    while (cur >= 0 && sim < limit) {
     int diff = pre - cur;

     // compare number of off diagonals
     int off_pre = A->p[pre + 1] - A->p[pre] - diff;
     int off_cur = A->p[cur + 1] - A->p[cur];

     if (off_pre != off_cur)
      break;

     // now examine row pattern
     bool col_sim = true;
     for (int j = 0; j < off_cur; j++) {
      int pre_index = A->p[pre + 1] - diff - j - 1;
      int cur_index = A->p[cur + 1] - j - 1;

      if (A->i[pre_index] != A->i[cur_index]) {
       col_sim = false;
       break;
      }
     }

     if (col_sim) {
      sim++;
      cur--;
     } else
      break;
    }
    i = pre - sim;

    if (sim > max)
     max = sim;
   }
  }

  M->supernodes = new int[num_nodes + 1];
  std::memcpy(M->supernodes, temp_nodes, sizeof(int) * num_nodes);
  if (isLower)
   M->supernodes[num_nodes] = n;
  else
   M->supernodes[num_nodes] = -1;

  delete[]temp_nodes;
  return num_nodes;
 }


 int BCSCMatrix::calcSize(CSC *A) {
  int i, j, p;
  int n = A->n;
  int *Ap = A->p;

  int nnz = 0;
  bool *in_block = new bool[n]();
  int *temprows = new int[n]();
  M->p = new int[n + 1]();

  for (i = 0; i < M->nodes; i++) {
   M->p[i] = nnz;

   // find the entire block size
   int nrow = 0;
   int width = M->supernodes[i + 1] - M->supernodes[i];
   for (j = M->supernodes[i]; j < M->supernodes[i + 1]; j++) {
    for (p = Ap[j]; p < Ap[j + 1]; p++) {
     int row_i = A->i[p];

     if (!in_block[row_i]) {
      in_block[row_i] = true;
      nnz += width;
      temprows[nrow] = row_i;
      nrow++;
     }
    }
   }
   for (j = 0; j < nrow; j++)
    in_block[temprows[j]] = false;
  }
  M->p[M->nodes] = nnz;

  delete[]temprows;
  delete[]in_block;
  return nnz;
 }


 void BCSCMatrix::createFormat(CSC *A) {
  int i, j, p;
  int counter = 0;

  int n = A->n;
  int *Ap = A->p;
  int *Ai = A->i;
  double *Ax = A->x;

  // temporary dynamic array to store Ai and Ax
  int *tempvec = new int[n]();
  bool *in_block = new bool[n]();
  std::priority_queue<int, std::vector<int>, std::greater<int>> queue;

  for (i = 0; i < M->nodes; i++) {
   // find the entire block size
   int nrow = 0;
   for (j = M->supernodes[i]; j < M->supernodes[i + 1]; j++) {
    for (p = Ap[j]; p < Ap[j + 1]; p++) {
     int row_i = Ai[p];

     if (!in_block[row_i]) {
      in_block[row_i] = true;
      queue.push(row_i);
      nrow++;
     }
    }
   }
   M->nrows[i] = nrow;

   // obtain order of the rows
   int temp_cnt = 0;
   while (!queue.empty()) {
    int row_i = queue.top();
    queue.pop();
    in_block[row_i] = false;
    tempvec[temp_cnt] = row_i;
    temp_cnt++;
   }

   // initialize temp_Ai and temp_Ax
   for (j = M->supernodes[i]; j < M->supernodes[i + 1]; j++) {
    int index = Ap[j];
    for (p = 0; p < temp_cnt; p++) {
     int row_i = tempvec[p];
     M->i[counter] = row_i;
     if (row_i == Ai[index] && index < Ap[j + 1]) {
      M->x[counter] = Ax[index];
      index++;
     } else {
      M->x[counter] = 0.0;
     }
     counter++;
    }
   }
  }

  delete[]tempvec;
  delete[]in_block;
 }


 CSC *BCSCMatrix::compressed_BCSC_to_CSC() {
  int nodes = M->nodes;
  const int *supernodes = M->supernodes;

  int counter = 0;
  auto supers = new int[M->n](); // A, B, C have same supernodal structure
  for (int i = 0; i < nodes; i++) {
   while (counter >= supernodes[i] && counter < supernodes[i + 1]) {
    supers[counter] = i;
    counter++;
   }
  }

  int nnz = 0;
  auto inplace = new bool[nodes]();

  int count = 0;
  auto temp = new int[nodes]();

  for (int i = 0; i < nodes; i++) {
   for (int r = 0; r < M->nrows[i]; r++) {
    int row = supers[M->i[M->p[i] + r]];
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

  CSC *Acsc = new CSC(nodes, nodes, nnz);
  auto Ap = Acsc->p;
  auto Ai = Acsc->i;

  int p = 1, i = 0;
  Ap[0] = 0;
  for (int j = 0; j < nodes; j++) {
   Ap[p] = Ap[p - 1];
   for (int r = 0; r < M->nrows[j]; r++) {
    int row = supers[M->i[M->p[j] + r]];
    if (!inplace[row]) {
     inplace[row] = true;
     temp[count] = row;
     count++;

     Ap[p]++;
     Ai[i] = row;
     i++;
    }
   }
   for (int k = 0; k < count; k++)
    inplace[temp[k]] = false;
   count = 0;

   Ai[i] = j + nodes;
   Ap[p]++;
   i++;
   p++;
  }

  delete[]inplace;
  delete[]supers;
  delete[]temp;
  return Acsc;
 }
}