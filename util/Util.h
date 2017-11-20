//
// Created by kazem on 7/25/17.
//

#ifndef CHOLOPENMP_UTIL_H
#define CHOLOPENMP_UTIL_H

#include <fstream>
#include <cassert>
#include <sstream>

/*
 * Useful conversion functions int <-->size_t
 */
/*int size_t2int(size_t val) {
    return (val <= INT_MAX) ? (int)((ssize_t)val) : -1;
}

size_t int2size_t(int val) {
    return (val < 0) ? __SIZE_MAX__ : (size_t)((unsigned)val);
}*/

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
    std::string line,banner, mtx, crd, arith, sym;
    /*  File format:
     *    %%MatrixMarket matrix coordinate real general/symmetric/...
     *    % ...
     *    % (optional comments)
     *    % ...
     *    #rows    #non-zero
     *    Triplet in the rest of lines: row    col    value
     */
    std::getline(inFile,line);
    for (unsigned i=0; i<line.length(); line[i]=tolower(line[i]),i++);
    std::istringstream iss(line);
    if (!(iss >> banner >> mtx >> crd >> arith >> sym)){
        std::cout<<"Invalid header (first line does not contain 5 tokens)\n";
        return false;
    }

    if(banner.compare("%%matrixmarket")) {
        std::cout<<"Invalid header (first token is not \"%%%%MatrixMarket\")\n";
        return false;
    }
    if(mtx.compare("matrix")) {
        std::cout<<"Not a matrix; this driver cannot handle that.\"\n";
        return false;
    }
    if(crd.compare("coordinate")) {
        std::cout<<"Not in coordinate format; this driver cannot handle that.\"\n";
        return false;
    }
    if(arith.compare("real")) {
        if(!arith.compare("complex")) {
            std::cout<<"Complex matrix; use zreadMM instead!\n";
            return false;
        }
        else if(!arith.compare("pattern")) {
            std::cout<<"Pattern matrix; values are needed!\n";
            return false;
        }
        else {
            std::cout<<"Unknown arithmetic\n";
            return false;
        }
    }
    while (!line.compare(0,1,"%"))
    {
        std::getline(inFile, line);
    }
    std::istringstream issDim(line);
    if (!(issDim >> n >> n >> NNZ)){
        std::cout<<"The matrix dimension is missing\n";
        return false;
    }
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

#endif //CHOLOPENMP_UTIL_H
