//
// Created by kazem on 6/13/20.
//
namespace sym_lib{



/*
 * L^T x = b
 */
 int ltsolve(int n, int *Lp, int *Li, double *Lx, double *x) {
  int p, j;
  if (!Lp || !Li || !x) return (0);                      /* check inputs */
  for (j = n - 1; j >= 0; j--) {
   for (p = Lp[j] + 1; p < Lp[j + 1]; p++) {
    x[j] -= Lx[p] * x[Li[p]];
   }
   x[j] /= Lx[Lp[j]];
  }
  return (1);
 }
}
