//
// Created by kazem on 7/18/17.
//

#ifndef TRIANGOPENMP_UTIL_H
#define TRIANGOPENMP_UTIL_H

#include <fstream>

/*
 * reading a CSC matrix from a coordinate file, stored col-ordered
 */
bool readMatrix(std::string fName, int &n, int &NNZ, int* &col,
                int* &row, double* &val){
    /*This function reads the input matrix from "fName" file and
     * allocate memory for matrix A, L and U.
     * - The input file is a coordinate version and e
     * ach row of the file shows (col, row, nnz)
     * - The matrices are zero-indexed
     */
    std::ifstream inFile;
    inFile.open(fName);
    inFile >> n;
    inFile >> n;
    inFile>>NNZ;
    int factorSize= (n * n) / 2;//Worst case assumption
    if(n <= 0 || NNZ <= 0)
        return false;
    col = new int[n + 1]();
    // colL = new int[n + 1]; colU = new int[n + 1];
    row = new int[NNZ];
    // rowL = new int[factorSize]; rowU = new int[factorSize];
    val = new double[NNZ];
    // valL = new double[factorSize]; valU = new double[factorSize];
    if(!val || !col || !row)
        return false;
    //Initializing the result vector
    int y, x, colCnt=0, nnzCnt=0;
    double value;

    col[0]=0;
    for (int i = 0; nnzCnt<NNZ; ) {//Reading from file row by row
        inFile>>x;x--;
        inFile>>y;y--;//zero indexing
        inFile>>value;
        if(y > n)
            return false;
        if(y==i){
            val[nnzCnt]=value;
            row[nnzCnt]=x;
            colCnt++; nnzCnt++;
        }
        else{//New col
            col[i+1]=col[i]+colCnt;
            i++;//next iteration
            colCnt=1;
            val[nnzCnt]=value;
            row[nnzCnt]=x;
            nnzCnt++;
        }

    }
    col[n]= col[n - 1] + colCnt;//last col

    return true;
}

/*
 *
 */

void rhsInit(int n, int *Ap, int *Ai, double *Ax, double *b){
    /*generating a rhs that produces a result of all 1 vector*/
    for (int j = 0; j < n; ++j) {
        b[j]=0;
    }
    for (int c = 0; c < n ; ++c) {
        for (int cc = Ap[c]; cc < Ap[c + 1]; ++cc) {
            b[Ai[cc]]+=Ax[cc];
        }
    }
}

/*
 * Avergae number of nnz in levels of a triangular solve, not needed!
 */

/*int avgNNZperLevel(int n, int *Lp, int *Li, int *Li_ptr, int *sup2col,
                   int nLevels, int *levelPtr, int *levelSet){
    int averageNNZperLevel=0;
    for (int i = 0; i < nLevels; ++i) {
        for (int l = levelPtr[i]; l < levelPtr[i + 1]; ++l) {
            int cc = levelSet[l];
            int curCol = cc != 0 ? sup2col[cc - 1] : 0;
            int nxtCol = sup2col[cc];
            int supWdt = nxtCol - curCol;
            int nSupR = Li_ptr[nxtCol] - Li_ptr[curCol];//row size of supernode
            averageNNZperLevel+=0;
        }
    }
}*/
#endif //TRIANGOPENMP_UTIL_H
