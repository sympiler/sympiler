//
// Created by kazem on 1/3/17.
//

#ifndef LEFT_LU_COL_MYBLAS_H
#define LEFT_LU_COL_MYBLAS_H
#ifdef MYBLAS
#define trsm_blas dlsolve_blas_nonUnit
#define matvec dmatvec_blas
#endif

#ifdef OPENBLAS
#endif


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

void superNodeDetection(int n, int *col, int *row,int *col2sup,int &supNo){
    int prev, cur;
    bool sim;
    supNo=0;
    col2sup[0]=0;
    for (int i = 1; i < n; ++i) {
        sim=true;
        for (prev = col[i-1], cur = col[i]; 1; ) {
            if((row[prev] == i-1 && prev < col[i]) || (row[prev] == i && prev < col[i])){
                //skip diagonal block of prev col
                ++prev;
                continue;
            }
            if(row[cur]==i){//skip diagonal block of cur col
                ++cur;
                continue;
            }
            if(prev - col[i] != cur - col[i+1]){ // the off-diagonals length
                sim=false;
                break;
            }else if(prev - col[i]==0)
                break;
            if(row[prev] != row[cur]){//now off-diagonals
                sim=false;
                break;
            }else{
                ++prev;++cur;
                if (prev >= col[i] || cur >= col[i+1])
                    break;
            }
        }
#if DEBUG >0
        if(i<0 || i>=n)
            printf("TTT \n");
#endif
        if(sim ){//col cur and nxt are similar
            col2sup[i]=supNo;
        }else{
            supNo++;
            col2sup[i]=supNo;
        }
    }
    supNo+= sim ;
}

void calcSize(int n, int *col, int *newCol, int *col2sup, int *sup2col, int supNo,
              int& newRowSize, int& newNNZ){
    newNNZ=0;newRowSize=0;
    newCol[0]=0;
    int curCol = 0, cnt=0, firstCol=0, tmpSize=0;
    for (int j = 0; j < supNo; ++j) {
        for (cnt=0; col2sup[curCol]==j && curCol<n; ++curCol, ++cnt);
        sup2col[j]=curCol;
        firstCol= j!=0 ? sup2col[j-1] : 0;
        cnt--;//To find the last col offset of supernode
        tmpSize=col[firstCol+cnt+1]-col[firstCol+cnt]+cnt;
        newRowSize+=tmpSize;//diagonal dense part + off diagonal sparse one
        for (int i = firstCol; i < curCol; ++i) {
            newNNZ+=tmpSize;
            newCol[i+1]=newNNZ;
        }
    }
#if DEBUG>0
    for (int i = 0; i < supNo; ++i) {//printing the result
        std::cout<<"Snode: "<<i<<" col: "<<sup2col[i]<<"\n";
    }
#endif
}
void createFormat(int n, int *col, int *row, double *val, int NNZ, int* &newRow, int newRowSize, double* &newVal,
                  int* &row_ptr, int* &newCol, int* &col2sup, int* &sup2col, int &supNo){
    int l=0, curOrig, wdth=0, c, origRow, firstCol;
    int firstRow;
    int offDiagofFirst;
    //Creating the new blocked format
    for (int s = 0; s < supNo; ++s) {
        //specifying the compressed row
        firstCol=s!=0 ? sup2col[s-1] : 0;
        firstRow = row[col[firstCol]];
        wdth=sup2col[s]-firstCol;
        row_ptr[firstCol]=l;
        offDiagofFirst=0;
        for (int r = 0; r < wdth; ++r, ++l, ++firstRow) {//row indices for the dense part
            newRow[l]=firstRow;
            if(row[col[firstCol]+offDiagofFirst]<=firstRow)//zero in dense part of first col
                offDiagofFirst++;
        }
        for (int r = col[firstCol]+offDiagofFirst; r < col[firstCol+1]; ++r, ++l) {
            newRow[l]=row[r];
        }
        for (int i = firstCol; i < sup2col[s]; ++i) {
            row_ptr[i] = row_ptr[firstCol];//copying the row number

            firstRow = row[col[firstCol]];
            for (c = newCol[i], curOrig=col[i]; c < newCol[i]+wdth; ++c, ++firstRow) {//Diagonal zero padding
                if(row[curOrig] == firstRow){
                    newVal[c] = val[curOrig];
                    curOrig++;
                }else{
                    newVal[c]=0;
                }
            }
            for (c = newCol[i]+wdth; c < newCol[i+1]; ++curOrig, ++c) {
                newVal[c]=val[curOrig];
            }
        }
    }
    row_ptr[n]=newRowSize;
}


#define CS_FLIP(i) (-(i)-2)
#define CS_UNFLIP(i) (((i) < 0) ? CS_FLIP(i) : (i))
#define CS_MARKED(w,j) (w [j] < 0)
#define CS_MARK(w,j) { w [j] = CS_FLIP (w [j]) ; }

/* depth-first-search of the graph of a matrix, starting at node j */
int dfs (int j, int *Gp, int *Gi, int top, int *xi, int *pstack, const int *pinv)
{
    int i, p, p2, done, jnew, head = 0 ;
    //if (!CS_CSC (G) || !xi || !pstack) return (-1) ;    /* check inputs */
    //Gp = G->p ; Gi = G->i ;
    xi [0] = j ;                /* initialize the recursion stack */

    while (head >= 0)
    {
        j = xi [head] ;         /* get j from the top of the recursion stack */
        jnew = pinv ? (pinv [j]) : j ;
        if (!CS_MARKED (Gp, j))
        {
            CS_MARK (Gp, j) ;       /* mark node j as visited */
            pstack [head] = (jnew < 0) ? 0 : CS_UNFLIP (Gp [jnew]) ;
        }
        done = 1 ;                  /* node j done if no unvisited neighbors */
        p2 = (jnew < 0) ? 0 : CS_UNFLIP (Gp [jnew+1]) ;
        for (p = pstack [head] ; p < p2 ; p++)  /* examine all neighbors of j */
        {
            i = Gi [p] ;            /* consider neighbor node i */
            if (CS_MARKED (Gp, i)) continue ;   /* skip visited node i */
            pstack [head] = p ;     /* pause depth-first search of node j */
            xi [++head] = i ;       /* start dfs at node i */
            done = 0 ;              /* node j is not done */
            break ;                 /* break, to start dfs (i) */
        }
        if (done)               /* depth-first search at node j is done */
        {
            head-- ;            /* remove j from the recursion stack */
            xi [--top] = j ;    /* and place in the output stack */
        }
    }
    return (top) ;
}

/* xi [top...n-1] = nodes reachable from graph of G*P' via nodes in B(:,k).
 * xi [n...2n-1] used as workspace */
int reach (int n, int* Gp, int *Gi, int* Bp, int* Bi, int k, int *pruneVal, const int *pinv)
{
    int p, top, setSize=0 ;
    int *xi = new int[2*n]();
    if (!Gp || !Gi || !Bp || !Bi || !xi) return (-1) ;    /* check inputs */
    //n = G->n ; Bp = B->p ; Bi = B->i ; Gp = G->p ;
    top = n ;
    for (p = Bp [k] ; p < Bp [k+1] ; p++)
    {
        if (!CS_MARKED (Gp, Bi [p]))    /* start a dfs at unmarked node i */
        {
            top = dfs (Bi [p], Gp, Gi, top, xi, xi+n, pinv) ;
        }
    }
    for (p = top ; p < n ; p++) CS_MARK (Gp, xi [p]) ;  /* restore G */
    for (int i = top; i < n; ++i) {
        pruneVal[setSize++] = xi[i];
    }
    delete []xi;
    return (setSize) ;
}

int reach_sn (int n, int* Gp, int *Gi, int* Bp, int* Bi,
              int k, int *PBset, const int *pinv, int sn,
              int *col2sup ){
    int p, top, PBsize=0 ;
    double analysis;
    bool *checked = new bool[sn]();
    int *xi = new int[2*n]();
    std::chrono::time_point<std::chrono::system_clock> start, end;
    std::chrono::duration<double> elapsed_seconds;
    start = std::chrono::system_clock::now();
    if (!Gp || !Gi || !Bp || !Bi || !xi) return (-1) ;    /* check inputs */
    //n = G->n ; Bp = B->p ; Bi = B->i ; Gp = G->p ;
    top = n ;
    for (p = Bp [k] ; p < Bp [k+1] ; p++)
    {
        if (!CS_MARKED (Gp, Bi [p]))    /* start a dfs at unmarked node i */
        {
            top = dfs (Bi [p], Gp, Gi, top, xi, xi+n, pinv) ;
        }
    }
    for (p = top ; p < n ; p++) CS_MARK (Gp, xi [p]) ;  /* restore G */
    end = std::chrono::system_clock::now();
    elapsed_seconds = end-start;
    analysis=elapsed_seconds.count();
    for (int i = top; i < n; ++i) {
        if(!checked[col2sup[xi[i]]]){
            checked[col2sup[xi[i]]] = true;
            PBset[PBsize++] = col2sup[xi[i]];
        }
    }
    delete []checked;
    delete []xi;
    return (PBsize) ;
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
#endif //LEFT_LU_COL_MYBLAS_H
