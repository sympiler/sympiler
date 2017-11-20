//
// Created by kazem on 11/20/17.
//



#include <cstdlib>
#include <chrono>
#include <string>
#include <fstream>
#include <iostream>
#include "Eigen/Sparse"
#include "../../util/Util.h"

using namespace std;


template <typename valType>
void load( int row, int *colP, int *rowIdx, double *val,
          std::vector<Eigen::Triplet<valType>>& triplets){
 for (int i = 0; i < row; ++i) {
  for (int j = colP[i]; j < colP[i+1]; ++j) {
   int curRow, curCol;
   valType curVal;
   curRow=rowIdx[j]; curCol=i; curVal=val[j];
   Eigen::Triplet<valType> tmp(curRow,curCol,curVal);
   triplets.push_back(tmp);
  }
 }
}


int main(int argc, char *argv[]){
    //std::string fName="/home/kazem/UFDB/Triangular/cbuckle.mtx";
 if(argc<2){
  printf("Please enter a path for the input matrix");
  return -1;
 }

 std::string fName = argv[1];

 int *col, *row;
 double  *val;
 int n, nnz;
 std::chrono::time_point<std::chrono::system_clock> start, end;
 std::chrono::duration<double> elapsed_seconds;
 if (!readMatrix(fName,n,nnz,col,row,val))
  return -1;
 double spFactor = 0.05;
 int rhsPercent = spFactor*n;

 int rown=n;
 double durationEigen=0;
 std::vector<Eigen::Triplet<double>> triplets;
 double *x=new double[n];
 //Eigen
 load(rown,col,row,val,triplets);
 Eigen::VectorXd b(rown),x1(rown);
 for (int i = 0; i < n; ++i) {
     b[i] = 0;
 }
 //b[1]=1;b[3]=20;b[4]=1;
 for (int i = n-rhsPercent; i < n; ++i) {
     b[i]=1;
 }
 Eigen::SparseMatrix<double> sm1(rown,rown);
 sm1.setFromTriplets(triplets.begin(),triplets.end());//col major
 start = std::chrono::system_clock::now();
 x1=sm1.triangularView<Eigen::Lower>().solve(b);
 end = std::chrono::system_clock::now();
 elapsed_seconds = end-start;
 durationEigen=elapsed_seconds.count();
 cout<<durationEigen<<",";

#if 0
    for(int i=0; i<20; ++i){
        cout<<x[i]<<",";
    }
    std::cout<<"\n";
#endif
    delete []x;
    return 1;
}
