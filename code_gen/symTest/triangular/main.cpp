#include <omp.h>
#include <cstdio>
#include <fstream>
#include <iostream>
#include "cholUtils.h"
#include <chrono>
#include "cbuckle_trns_gen.h"
#include "../../util/Util.h"
#include "../../sympiler/NumericalUtils.h"


#define CSC
#undef RUNALL


int main(int argc, char *argv[])  {
    if(argc<2){
        printf("Please enter a path for the input matrix");
        return -1;
    }

    std::string fName = argv[1];

    int32_t *col, *row;
    double  *y, *val;
    int32_t n, nnz;
    std::chrono::time_point<std::chrono::system_clock> start, end;
    std::chrono::duration<double> elapsed_seconds;
    double duration4 = 0 ,duration3 = 0, duration2=0, duration1=0;
    if (!readMatrix(fName,n,nnz,col,row,val))
        return -1;
    double *x1=new double[n]();
    double spFactor = 0.05;
    int rhsPercent = spFactor*n;

    //***************CSC baseline code
    rhsInit(n,col,row,val,x1);
    for (int i = 0; i < n-rhsPercent-1; ++i) {
        x1[i]=0;
    }
    start = std::chrono::system_clock::now();
//    lsolve(n,col,row,val,x1);
    end = std::chrono::system_clock::now();
    elapsed_seconds = end-start;
    duration1=elapsed_seconds.count();
//    std::cout<<duration1<<",";

    //*****************setting up the RHS
    int *Bp = new int[2]; Bp[0]=0; Bp[1]=rhsPercent+1;
    int *Bi = new int[rhsPercent+1];
    int *pruneSet = new int[2*n]();
    double *xOut = new double[n]();

    for (int i = n-rhsPercent-1,cnt=0; i < n; ++i) {
        Bi[cnt++]=i;
    }

#ifdef BLOCKED
    //***************Sympiler-generated code
    int *col2sup = new int[n];
    int supNo=0, newNNZ=0, newRowSize=0;
    superNodeDetection(n,col,row,col2sup,supNo);
    int *sup2col = new int[supNo];
    int *newCol = new int[n+1];
    calcSize(n,col,newCol,col2sup,sup2col,supNo,newRowSize,newNNZ);
    //int average = averageSupNode(sup2col,supNo);
    int *newRow = new int[newRowSize+1];
    double *newVal = new double[newNNZ];
    int *rowP = new int[n+1];
    createFormat(n,col,row,val,nnz,newRow,newRowSize,newVal,rowP,newCol,
                 col2sup,sup2col,supNo);

    int top = reach_sn(n, col, row, Bp, Bi, 0, pruneSet, 0,supNo,col2sup);
    //creatBlockFormat(n, col, row,);
    start = std::chrono::system_clock::now();
    //lsolve_sup(n,newCol,newRow,newVal,NNZ,rowP,col2sup,sup2col,supNo,x1);
    trns(n,newCol,newRow,newVal,rowP,
         0, Bp,Bi,&x2[n-rhsPercent-1],NULL,
         xOut,
         NULL,pruneSet,top,sup2col,supNo);
    end = std::chrono::system_clock::now();
    elapsed_seconds = end-start;
    duration1=elapsed_seconds.count();
    std::cout<<duration1<<",";

    delete []col2sup;
    delete []sup2col;
    delete []newCol;
#elif PRUNE
    //***************Sympiler-generated code

    int top = reach(n, col, row, Bp, Bi, 0, pruneSet, 0);
    start = std::chrono::system_clock::now();

    trns(n,col,row,val,NULL,
         0, Bp,Bi,&x2[n-rhsPercent-1],NULL,
         xOut,
         NULL,pruneSet,top,NULL,NULL);
#else
    int top = 0;
    start = std::chrono::system_clock::now();
    trns(n,col,row,val,NULL,
         0, Bp,Bi,&x1[n-rhsPercent-1],NULL,
         xOut);
    end = std::chrono::system_clock::now();

    double max = 0;
    double *temp = new double[n]();
    for (int i = 0; i < n; i++) {
     for (int j = col[i]; j < col[i + 1]; j++) {
      temp[row[j]] += val[j] * xOut[i];
     }
    }
    for(int i = 0; i < n; i++) {
     if(std::abs(temp[i] - x1[i]) > max)
      max = std::abs(temp[i] - x1[i]);
    }
    std::cout << max << "\n";

    elapsed_seconds = end-start;
    duration1=elapsed_seconds.count();
//    std::cout<<duration1<<",";



#endif


#if 0
    //Testing
    int test=0;
    for (int i = 0; i < n; ++i) {
        if(x1[i]-xOut[i]>0.001)
            test++;
    }
    if(test>1)
        std::cout<<"Error margin is high:"<<test<<",";
    /*else
        std::cout<<"WELL DONE!\n";*/
#endif
    delete []pruneSet;
    delete []Bp;
    delete []Bi;
    delete []xOut;
    delete []x1;
//    delete []x2;

}