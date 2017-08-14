//
// Created by kazem on 7/25/17.
//

#ifndef CHOLOPENMP_UTIL_H
#define CHOLOPENMP_UTIL_H

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


bool enableColdCache(int n, std::ifstream &f){
    /*
     * n specifies the size of data for double computation. It depends
     * on the cache size
     */
    //TODO check file during read
    assert(!f.fail());
    double curVal;
    double **waste=new double*[n];
    for (int i = 0; i < n; ++i) {
        waste[i] = new double[n];
    }
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            f>>curVal;
            waste[i][j]=curVal;
        }
    }
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            for (int k = 0; k < n; ++k) {
                waste[i][j] += waste[i][k]*waste[k][j];
            }
        }
    }
    for (int i = 0; i < n; ++i) {
        delete waste[i];
    }
    delete waste;
    return true;
}

#endif //CHOLOPENMP_UTIL_H
