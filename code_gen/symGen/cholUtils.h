//
// Created by kazem on 1/3/17.
//

#ifndef LEFT_LU_COL_CHOLUTILS_H
#define LEFT_LU_COL_CHOLUTILS_H
#include <math.h>
#include <mkl.h>

#ifdef MYBLAS
#define trsm_blas dlsolve_blas_nonUnit
#define matvec dmatvec_blas
#endif

#ifdef OPENBLAS
#define syrk dsyrk_
#define gemm dgemm_
#define potrf dpotrf_
#define trsm dtrsm_
#endif


#ifdef MKL
#define syrk(X1,X2,X3,X4,X5,X6,X7,X8) dsyrk("L","N",X1,X2,X3,X4,X5,X6,X7,X8)
#define gemm(X1,X2,X3,X4,X5,X6,X7,X8,X9,X10,X11) dgemm("N","C",X1,X2,X3,X4,X5,X6,X7,X8,X9,X10,X11)
#define potrf(X1,X2,X3,X4) dpotrf("L",X1,X2,X3,X4)
#define trsm(X1,X2,X3,X4,X5,X6,X7) dtrsm("R", "L", "C", "N",X1,X2,X3,X4,X5,X6,X7)
#endif


/*
 * Symbolic Calls
 */

int lSolve_dense(int colSize,int col, double *M, double *rhs){
    for (int i = 0; i < col; ++i) {
        rhs[i]/=M[i*colSize+i];
        for (int j = i+1; j < col; ++j) {
            rhs[j]-=M[i*colSize+j]*rhs[i];
        }
    }
    return 1;
}

int lSolve_dense_col(int colSize,int col, double *M, double *rhs){
    for (int i = 0; i < col; ++i) {
        rhs[i*colSize]/=M[i*colSize+i];
        for (int j = i+1; j < col; ++j) {
            rhs[j*colSize]-=M[i*colSize+j]*rhs[i*colSize];
        }
    }
    return 1;
}

void dlsolve_blas ( int ldm, int ncol, double *M, double *rhs )//unit triangular solver
{
    int k;
    double x0, x1, x2, x3, x4, x5, x6, x7;
    double *M0;
    register double *Mki0, *Mki1, *Mki2, *Mki3, *Mki4, *Mki5, *Mki6, *Mki7;
    register int firstcol = 0;

    M0 = &M[0];

    while ( firstcol < ncol - 7 ) { /* Do 8 columns */
        Mki0 = M0 + 1;
        Mki1 = Mki0 + ldm + 1;
        Mki2 = Mki1 + ldm + 1;
        Mki3 = Mki2 + ldm + 1;
        Mki4 = Mki3 + ldm + 1;
        Mki5 = Mki4 + ldm + 1;
        Mki6 = Mki5 + ldm + 1;
        Mki7 = Mki6 + ldm + 1;

        x0 = rhs[firstcol];
        x1 = rhs[firstcol+1] - x0 * *Mki0++;
        x2 = rhs[firstcol+2] - x0 * *Mki0++ - x1 * *Mki1++;
        x3 = rhs[firstcol+3] - x0 * *Mki0++ - x1 * *Mki1++ - x2 * *Mki2++;
        x4 = rhs[firstcol+4] - x0 * *Mki0++ - x1 * *Mki1++ - x2 * *Mki2++
             - x3 * *Mki3++;
        x5 = rhs[firstcol+5] - x0 * *Mki0++ - x1 * *Mki1++ - x2 * *Mki2++
             - x3 * *Mki3++ - x4 * *Mki4++;
        x6 = rhs[firstcol+6] - x0 * *Mki0++ - x1 * *Mki1++ - x2 * *Mki2++
             - x3 * *Mki3++ - x4 * *Mki4++ - x5 * *Mki5++;
        x7 = rhs[firstcol+7] - x0 * *Mki0++ - x1 * *Mki1++ - x2 * *Mki2++
             - x3 * *Mki3++ - x4 * *Mki4++ - x5 * *Mki5++
             - x6 * *Mki6++;

        rhs[++firstcol] = x1;
        rhs[++firstcol] = x2;
        rhs[++firstcol] = x3;
        rhs[++firstcol] = x4;
        rhs[++firstcol] = x5;
        rhs[++firstcol] = x6;
        rhs[++firstcol] = x7;
        ++firstcol;

        for (k = firstcol; k < ncol; k++)
            rhs[k] = rhs[k] - x0 * *Mki0++ - x1 * *Mki1++
                     - x2 * *Mki2++ - x3 * *Mki3++
                     - x4 * *Mki4++ - x5 * *Mki5++
                     - x6 * *Mki6++ - x7 * *Mki7++;

        M0 += 8 * ldm + 8;
    }

    while ( firstcol < ncol - 3 ) { /* Do 4 columns */
        Mki0 = M0 + 1;
        Mki1 = Mki0 + ldm + 1;
        Mki2 = Mki1 + ldm + 1;
        Mki3 = Mki2 + ldm + 1;

        x0 = rhs[firstcol];
        x1 = rhs[firstcol+1] - x0 * *Mki0++;
        x2 = rhs[firstcol+2] - x0 * *Mki0++ - x1 * *Mki1++;
        x3 = rhs[firstcol+3] - x0 * *Mki0++ - x1 * *Mki1++ - x2 * *Mki2++;

        rhs[++firstcol] = x1;
        rhs[++firstcol] = x2;
        rhs[++firstcol] = x3;
        ++firstcol;

        for (k = firstcol; k < ncol; k++)
            rhs[k] = rhs[k] - x0 * *Mki0++ - x1 * *Mki1++
                     - x2 * *Mki2++ - x3 * *Mki3++;

        M0 += 4 * ldm + 4;
    }

    if ( firstcol < ncol - 1 ) { /* Do 2 columns */
        Mki0 = M0 + 1;
        Mki1 = Mki0 + ldm + 1;

        x0 = rhs[firstcol];
        x1 = rhs[firstcol+1] - x0 * *Mki0++;

        rhs[++firstcol] = x1;
        ++firstcol;

        for (k = firstcol; k < ncol; k++)
            rhs[k] = rhs[k] - x0 * *Mki0++ - x1 * *Mki1++;
    }
}


void dlsolve_blas_nonUnit ( int ldm, int ncol, double *M, double *rhs )//general triangular solver
{
    int k;
    double x0, x1, x2, x3, x4, x5, x6, x7;
    double *M0;
    register double *Mki0, *Mki1, *Mki2, *Mki3, *Mki4, *Mki5, *Mki6, *Mki7;
    register int firstcol = 0;

    M0 = &M[0];

    while ( firstcol < ncol - 7 ) { /* Do 8 columns */
        Mki0 = M0 ;
        Mki1 = Mki0 + ldm + 1;
        Mki2 = Mki1 + ldm + 1;
        Mki3 = Mki2 + ldm + 1;
        Mki4 = Mki3 + ldm + 1;
        Mki5 = Mki4 + ldm + 1;
        Mki6 = Mki5 + ldm + 1;
        Mki7 = Mki6 + ldm + 1;

        x0 = rhs[firstcol]/ *Mki0++ ;
        x1 = (rhs[firstcol+1] - x0 * *Mki0++)/ *Mki1++;
        x2 = (rhs[firstcol+2] - x0 * *Mki0++ - x1 * *Mki1++)/ *Mki2++;
        x3 = (rhs[firstcol+3] - x0 * *Mki0++ - x1 * *Mki1++ - x2 * *Mki2++)/ *Mki3++;
        x4 = (rhs[firstcol+4] - x0 * *Mki0++ - x1 * *Mki1++ - x2 * *Mki2++
             - x3 * *Mki3++)/ *Mki4++;
        x5 = (rhs[firstcol+5] - x0 * *Mki0++ - x1 * *Mki1++ - x2 * *Mki2++
             - x3 * *Mki3++ - x4 * *Mki4++)/ *Mki5++;
        x6 = (rhs[firstcol+6] - x0 * *Mki0++ - x1 * *Mki1++ - x2 * *Mki2++
             - x3 * *Mki3++ - x4 * *Mki4++ - x5 * *Mki5++)/ *Mki6++;
        x7 = (rhs[firstcol+7] - x0 * *Mki0++ - x1 * *Mki1++ - x2 * *Mki2++
             - x3 * *Mki3++ - x4 * *Mki4++ - x5 * *Mki5++
             - x6 * *Mki6++)/ *Mki7++;

        rhs[firstcol++] = x0;
        rhs[firstcol++] = x1;
        rhs[firstcol++] = x2;
        rhs[firstcol++] = x3;
        rhs[firstcol++] = x4;
        rhs[firstcol++] = x5;
        rhs[firstcol++] = x6;
        rhs[firstcol++] = x7;

        for (k = firstcol; k < ncol; k++)
            rhs[k] = rhs[k] - x0 * *Mki0++ - x1 * *Mki1++
                     - x2 * *Mki2++ - x3 * *Mki3++
                     - x4 * *Mki4++ - x5 * *Mki5++
                     - x6 * *Mki6++ - x7 * *Mki7++;

        M0 += 8 * ldm + 8;
    }

    while ( firstcol < ncol - 3 ) { /* Do 4 columns */
        Mki0 = M0 ;
        Mki1 = Mki0 + ldm + 1 ;
        Mki2 = Mki1 + ldm + 1;
        Mki3 = Mki2 + ldm + 1;

        x0 = rhs[firstcol]/ *Mki0++;
        x1 = (rhs[firstcol+1] - x0 * *Mki0++) / *Mki1++;
        x2 = (rhs[firstcol+2] - x0 * *Mki0++ - x1 * *Mki1++) / *Mki2++;
        x3 = (rhs[firstcol+3] - x0 * *Mki0++ - x1 * *Mki1++ - x2 * *Mki2++) / *Mki3++;

        rhs[firstcol++] = x0;
        rhs[firstcol++] = x1;
        rhs[firstcol++] = x2;
        rhs[firstcol++] = x3;

        for (k = firstcol; k < ncol; k++)
            rhs[k] = rhs[k] - x0 * *Mki0++ - x1 * *Mki1++
                     - x2 * *Mki2++ - x3 * *Mki3++;

        M0 += 4 * ldm + 4;
    }

    if ( firstcol < ncol - 1 ) { /* Do 2 columns */
        Mki0 = M0 ;
        Mki1 = Mki0 + ldm + 1;

        x0 = rhs[firstcol]/ *Mki0++;
        x1 = (rhs[firstcol+1] - x0 * *Mki0++)/ *Mki1++;

        rhs[firstcol++] = x0;
        rhs[firstcol++] = x1;

        for (k = firstcol; k < ncol; k++)
            rhs[k] = rhs[k] - x0 * *Mki0++ - x1 * *Mki1++;
        M0 += 2 * ldm + 2;
    }

    if ( firstcol == ncol - 1 ) { /* Do 1 columns */
        Mki0 = M0 ;
        x0 = rhs[firstcol]/ *Mki0;
        rhs[firstcol] = x0;
    }
}

int matVector(int colSize, int row, int col, double *M, double *vec, double *resVec){
    for (int i = 0; i < col; ++i) {
        for (int j = 0; j < row; ++j) {
            resVec[j] += M[i*colSize+j]*vec[i];
        }
    }
    return 1;
}

int matVectorT(int colSize, int row, int col, double *M, double *vec, double *resVec){
    for (int i = 0; i < col; ++i) {
        for (int j = 0; j < row; ++j) {
            resVec[j] += M[i*colSize+j]*vec[col-i];
        }
    }
    return 1;
}

void dmatvec_blas (
        int ldm,	/* in -- leading dimension of M */
        int nrow,	/* in */
        int ncol,	/* in */
        double *M,	/* in */
        double *vec,	/* in */
        double *Mxvec	/* in/out */
)
{
    double vi0, vi1, vi2, vi3, vi4, vi5, vi6, vi7;
    double *M0;
    register double *Mki0, *Mki1, *Mki2, *Mki3, *Mki4, *Mki5, *Mki6, *Mki7;
    register int firstcol = 0;
    int k;

    M0 = &M[0];
    while ( firstcol < ncol - 7 ) {	/* Do 8 columns */

        Mki0 = M0;
        Mki1 = Mki0 + ldm;
        Mki2 = Mki1 + ldm;
        Mki3 = Mki2 + ldm;
        Mki4 = Mki3 + ldm;
        Mki5 = Mki4 + ldm;
        Mki6 = Mki5 + ldm;
        Mki7 = Mki6 + ldm;

        vi0 = vec[firstcol++];
        vi1 = vec[firstcol++];
        vi2 = vec[firstcol++];
        vi3 = vec[firstcol++];
        vi4 = vec[firstcol++];
        vi5 = vec[firstcol++];
        vi6 = vec[firstcol++];
        vi7 = vec[firstcol++];

        for (k = 0; k < nrow; k++)
            Mxvec[k] += vi0 * *Mki0++ + vi1 * *Mki1++
                        + vi2 * *Mki2++ + vi3 * *Mki3++
                        + vi4 * *Mki4++ + vi5 * *Mki5++
                        + vi6 * *Mki6++ + vi7 * *Mki7++;

        M0 += 8 * ldm;
    }

    while ( firstcol < ncol - 3 ) {	/* Do 4 columns */

        Mki0 = M0;
        Mki1 = Mki0 + ldm;
        Mki2 = Mki1 + ldm;
        Mki3 = Mki2 + ldm;

        vi0 = vec[firstcol++];
        vi1 = vec[firstcol++];
        vi2 = vec[firstcol++];
        vi3 = vec[firstcol++];
        for (k = 0; k < nrow; k++)
            Mxvec[k] += vi0 * *Mki0++ + vi1 * *Mki1++
                        + vi2 * *Mki2++ + vi3 * *Mki3++ ;

        M0 += 4 * ldm;
    }

    while ( firstcol < ncol ) {		/* Do 1 column */

        Mki0 = M0;
        vi0 = vec[firstcol++];
        for (k = 0; k < nrow; k++)
            Mxvec[k] += vi0 * *Mki0++;

        M0 += ldm;
    }

}

void cholesky(long N, double *A, double *diag)
{

/* Subroutine cholesky:

   Arguments:
              N  dimension of A

              A
                 on entry: the N by N matrix to be decomposed
                 on exit: upper triangle is still A
                          lower sub-triangle is the sub-trangle
                          of the Cholesky factor L

              diag
                 on entry: an arbitrary vector of length N
                 on exit: the diagonal of the Cholesky factor L

*/
    long i,j,k;
    for(j=0;j<N;j++)
        diag[j] = A[N*j+j];
    for(j=0;j<N;j++)
    {
        for(k=0;k<j;k++)
            diag[j] -= A[N*k+j]*A[N*k+j];
        diag[j] = sqrt(diag[j]);
        for(i=j+1;i<N;i++)
        {
            for(k=0;k<j;k++)
                A[N*j+i] -= A[N*k+i]*A[N*k+j];
            A[N*j+i]/=diag[j];
        }
    }
}
void Cholesky_col(int n, int dim, double* a){
    double tmp = 0;
    for (int j = 0; j < dim; ++j) {
        for (int k = 0; k < j; ++k) {
            tmp = a[k*n+j];
            for (int i = j; i < dim; ++i) {
                a[j*n+i] = a[j*n+i] - a[k*n+i]*tmp;
            }
        }
        tmp = sqrt(a[j*n+j]);
        for (int k = j+1; k < dim; ++k) {
            a[j*n+k] = a[j*n+k] / tmp;
        }
        a[j*n+j] = tmp;
    }
}

// Cholesky requires the matrix to be symmetric positive-definite
void Cholesky(int d,double*S,double*D){
    for(int k=0;k<d;++k){
        double sum=0.;
        for(int p=0;p<k;++p)sum+=D[k*d+p]*D[k*d+p];
        D[k*d+k]=sqrt(S[k*d+k]-sum);
        for(int i=k+1;i<d;++i){
            double sum=0.;
            for(int p=0;p<k;++p)sum+=D[i*d+p]*D[k*d+p];
            D[i*d+k]=(S[i*d+k]-sum)/D[k*d+k];
        }
    }
}
// This version could be more efficient on some architectures
// Use solveCholesky for both Cholesky decompositions
void CholeskyRow(int d,double*S,double*D){
    for(int k=0;k<d;++k){
        for(int j=0;j<d;++j){
            double sum=0.;
            for(int p=0;p<j;++p)sum+=D[k*d+p]*D[j*d+p];
            D[k*d+j]=(S[k*d+j]-sum)/D[j*d+j];
        }
        double sum=0.;
        for(int p=0;p<k;++p)sum+=D[k*d+p]*D[k*d+p];
        D[k*d+k]=sqrt(S[k*d+k]-sum);
    }
}



/*
 ****** Serial implementation
 */
int lsolve (int n, int* Lp, int* Li, double* Lx, double *x){
    int p, j;
    if (!Lp || !Li || !x) return (0) ;                     /* check inputs */
    for (j = 0 ; j < n ; j++)
    {
        x [j] /= Lx [Lp [j]] ;
        for (p = Lp [j]+1 ; p < Lp [j+1] ; p++)
        {
            x [Li [p]] -= Lx [p] * x [j] ;
        }
    }
    return (1) ;
}
//void creatBlockedFormat()
#endif //LEFT_LU_COL_CHOLUTILS_H
